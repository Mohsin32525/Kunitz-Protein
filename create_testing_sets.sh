#Before to start with the analysis we have to select the ids of the proteins that are maintained by Structural alignment to remove them from the test set!
grep ">" pdb_kunitz_rp_strali.sto | sed 's/>//' | sed 's/:/_/' > ids_to_keep_forseq.txt

#Use the ids to extract the 23 kunitz sequences aligned by structural alignment
awk 'BEGIN {
  while ((getline < "ids_to_keep_forseq.txt") > 0) ids[$1] = 1
}
{
  if ($0 ~ /^>/) {
    id = substr($0, 2)
    keep = ids[id]
  }
  if (keep) print
}' pdb_kunitz_rp.fasta > pdb_kunitz_rp_clean.fasta

rm -f ids_to_keep_forseq.txt #remove the temporary file that contain the ids

# Create a BLAST database from all Kunitz proteins
makeblastdb -in all_kunitz.fasta -input_type fasta -dbtype prot -out all_kunitz.fasta

# Perform BLAST search of the 23 representative sequences against the full Kunitz dataset
blastp -query pdb_kunitz_rp_clean.fasta -db all_kunitz.fasta -out tmp_pdb_kunitz_rp_clean.blast -outfmt 7

# Extract Uniprot IDs of sequences with high identity (≥95%) and alignment length ≥50 to remove them from the training/testing pool
grep -v "^#" tmp_pdb_kunitz_rp_clean.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u | cut -d "|" -f 2 > tmp_pdb_kunitz_rp_clean_ids.txt

# Extract all IDs from the full Kunitz dataset
grep ">" all_kunitz.fasta | cut -d "|" -f 2 > tmp_all_kunitz.id

echo 'Creating positive test files (set1 and set2)...'

# Remove the overlapping IDs and keep only the remaining Kunitz proteins
comm -23 <(sort tmp_all_kunitz.id) <(sort tmp_pdb_kunitz_rp_clean_ids.txt) > tmp_to_keep.ids
sort -R tmp_to_keep.ids > tmp_random_ok_kunitz.id
head -n 184 tmp_random_ok_kunitz.id > tmp_pos_1.id
tail -n 184 tmp_random_ok_kunitz.id > tmp_pos_2.id
python3 get_seq.py tmp_pos_1.id uniprot_sprot.fasta > pos_1.fasta
python3 get_seq.py tmp_pos_2.id uniprot_sprot.fasta > pos_2.fasta

echo 'Creating negative test files (set1 and set2)...'

# Extract all Swiss-Prot IDs
grep ">" uniprot_sprot.fasta | cut -d "|" -f 2 > tmp_sp.id

# Remove the Kunitz proteins to get the negative candidates
comm -23 <(sort tmp_sp.id) <(sort tmp_all_kunitz.id) > tmp_sp_negs.ids
sort -R tmp_sp_negs.ids > tmp_random_sp_negs.id
head -n 286286 tmp_random_sp_negs.id > tmp_neg_1.id
tail -n 286286 tmp_random_sp_negs.id > tmp_neg_2.id
python3 get_seq.py tmp_neg_1.id uniprot_sprot.fasta > neg_1.fasta
python3 get_seq.py tmp_neg_2.id uniprot_sprot.fasta > neg_2.fasta

echo 'FASTA files created.'


#STRUCTURAL ALIGNMENT HMM RESULTS
# Run hmmsearch on all positive and negative FASTA files and create tabular output. Generates .out files that contains the E-values computed on the HMM of the sequence alignment.
hmmsearch -Z 1000 --max --tblout pos_1_strali.out pdb_kunitz_rp_strali.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2_strali.out pdb_kunitz_rp_strali.hmm pos_2.fasta
hmmsearch -Z 1000 --max --tblout neg_1_strali.out pdb_kunitz_rp_strali.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2_strali.out pdb_kunitz_rp_strali.hmm neg_2.fasta
# Extract E-values and build .class files with: ID, label (1/0), full sequence E-value ($5), and best domain E-value ($8)
grep -v "^#" pos_1_strali.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}' > pos_1_strali.class
grep -v "^#" pos_2_strali.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}' > pos_2_strali.class
grep -v "^#" neg_1_strali.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}' > neg_1_strali.class
grep -v "^#" neg_2_strali.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}' > neg_2_strali.class

# Add true negatives not detected by HMM (missing from the .out), with a fake high E-value (10.0)
comm -23 <(sort tmp_neg_1.id) <(cut -f 1 neg_1_strali.class | sort) | awk '{print $1"\t"0"\t"10.0"\t"10.0}' >> neg_1_strali.class
comm -23 <(sort tmp_neg_2.id) <(cut -f 1 neg_2_strali.class | sort) | awk '{print $1"\t"0"\t"10.0"\t"10.0}' >> neg_2_strali.class

# Combine positive and negative sets into final training/testing files
cat pos_1_strali.class neg_1_strali.class > set_1_strali.class
cat pos_2_strali.class neg_2_strali.class > set_2_strali.class
cat set_1_strali.class set_2_strali.class > temp_overall_strali.class

# Automatically determine the best thresholds (based on highest MCC)
# SET 1
set1_best_evalue_full_seq=$(
    for i in $(seq 1 9); do
        python3 performance.py set_1_strali.class 1e-"$i"
    done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)
set1_best_evalue_one_domain=$(
    for i in $(seq 1 9); do
        python3 performance.py set_1_strali.class 1e-"$i"
    done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)

# SET 2
set2_best_evalue_full_seq=$(
    for i in $(seq 1 9); do
        python3 performance.py set_2_strali.class 1e-"$i"
    done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)
set2_best_evalue_one_domain=$(
    for i in $(seq 1 9); do
        python3 performance.py set_2_strali.class 1e-"$i"
    done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)




# Run final performance evaluations using the best thresholds
# Identify false positives: non-Kunitz proteins with E-value below threshold (misclassified as positive)
# Identify false negatives: Kunitz proteins with E-value above threshold (misclassified as negative)
#SET_1 BEST THRESHOLD USING FULL SEQUENCE E-VALUE
echo -e "BEST THRESHOLD OBTAINED FROM SET_1 FULL SEQUENCE: $set1_best_evalue_full_seq" > hmm_results_strali.txt
    #SET_2 TEST
echo -e "\nPERFORMANCES SET_2 USING E-VALUE THRESHOLD OF SET_1 - FULL SEQUENCES" >> hmm_results_strali.txt
python3 performance.py set_2_strali.class "$set1_best_evalue_full_seq" 1 >> hmm_results_strali.txt
        #false positives set_2
echo -e "False positives for set 2 considering full sequence set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set1_best_evalue_full_seq" '$3 < num {print $1, $2, $3}' neg_2_strali.class | sort -grk 3 >> hmm_results_strali.txt
        #false negatives set_2
echo -e "False negatives for set 2 considering full sequence set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set1_best_evalue_full_seq" '$3 > num {print $1, $2, $3}' pos_2_strali.class | sort -grk 3 >> hmm_results_strali.txt
    #OVERALL TEST
echo -e "\nOVERALL PERFORMANCES USING E-VALUE THRESHOLD OF SET_1 - FULL SEQUENCES" >> hmm_results_strali.txt
python3 performance.py temp_overall_strali.class "$set1_best_evalue_full_seq" 1 >> hmm_results_strali.txt
        #false positives overall
echo -e "False positives for overall considering full sequence set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set1_best_evalue_full_seq" '$3 < num  {print $1, $2, $3}' neg_2_strali.class neg_1_strali.class | sort -grk 3 >> hmm_results_strali.txt
        #false negatives overall
echo -e "False negatives for overall considering full sequence set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set1_best_evalue_full_seq" '$3 > num {print $1, $2, $3}' pos_1_strali.class pos_2_strali.class | sort -grk 3 >> hmm_results_strali.txt 

#SET_2 BEST THRESHOLD USING FULL SEQUENCE E-VALUE
echo -e "\n\n\nBEST THRESHOLD OBTAINED FROM SET_2 FULL SEQUENCE: $set2_best_evalue_full_seq" >> hmm_results_strali.txt
    #SET_1 TEST
echo -e "\nPERFORMANCES SET_1 USING E-VALUE THRESHOLD OF SET_2 - FULL SEQUENCES" >> hmm_results_strali.txt
python3 performance.py set_1_strali.class "$set2_best_evalue_full_seq" 1 >> hmm_results_strali.txt
        #false positives set_1
echo -e "False positives for set 1 considering full sequence set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set2_best_evalue_full_seq" '$3 < num {print $1, $2, $3}' neg_1_strali.class | sort -grk 3 >> hmm_results_strali.txt
        #false negatives set_1
echo -e "False negatives for set 1 considering full sequence set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set2_best_evalue_full_seq" '$3 > num {print $1, $2, $3}' pos_1_strali.class | sort -grk 3 >> hmm_results_strali.txt
    #OVERALL TEST
echo -e "\nOVERALL PERFORMANCES USING E-VALUE THRESHOLD OF SET_2 - FULL SEQUENCES" >> hmm_results_strali.txt
python3 performance.py temp_overall_strali.class "$set2_best_evalue_full_seq" 1 >> hmm_results_strali.txt
        #false positives overall
echo -e "False positives for overall considering full sequence set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set2_best_evalue_full_seq" '$3 < num  {print $1, $2, $3}' neg_2_strali.class neg_1_strali.class | sort -grk 3 >> hmm_results_strali.txt
        #false negatives overall
echo -e "False negatives for overall considering full sequence set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set2_best_evalue_full_seq" '$3 > num {print $1, $2, $3}' pos_1_strali.class pos_2_strali.class | sort -grk 3 >> hmm_results_strali.txt

#SET_1 BEST THRESHOLD USING SINGLE DOMAIN E-VALUE
echo -e "\n\n\n BEST THRESHOLD OBTAINED FROM SET_1 SINGLE DOMAIN: $set1_best_evalue_one_domain" >> hmm_results_strali.txt
    #SET_2 TEST
echo -e "\nPERFORMANCES SET_2 USING E-VALUE THRESHOLD OF SET_1 - SINGLE DOMAIN" >> hmm_results_strali.txt
python3 performance.py set_2_strali.class "$set1_best_evalue_one_domain" 2 >> hmm_results_strali.txt
        #false positives set_2
echo -e "False positives for set 2 considering single domain set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set1_best_evalue_one_domain" '$4 < num {print $1, $2, $4}' neg_2_strali.class | sort -grk 3 >> hmm_results_strali.txt
        #false negatives set_2
echo -e "False negatives for set 2 considering single domain set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set1_best_evalue_one_domain" '$4 > num {print $1, $2, $4}' pos_2_strali.class | sort -grk 3 >> hmm_results_strali.txt
    #OVERALL
echo -e "\nOVERALL PERFORMANCES USING E-VALUE THRESHOLD OF SET_1 - SINGLE DOMAIN" >> hmm_results_strali.txt
python3 performance.py temp_overall_strali.class "$set1_best_evalue_one_domain" 2 >> hmm_results_strali.txt
        #false positives overall
echo -e "False positives for overall considering single domain set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set1_best_evalue_one_domain" '$4 < num  {print $1, $2, $4}' neg_2_strali.class neg_1_strali.class | sort -grk 3 >> hmm_results_strali.txt   
        #false negatives overall     
echo -e "False negatives for overall considering single domain set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set1_best_evalue_one_domain" '$4 > num {print $1, $2, $4}' pos_1_strali.class pos_2_strali.class | sort -grk 3 >> hmm_results_strali.txt

#SET_2 BEST THRESHOLD USING SINGLE DOMAIN E-VALUE
echo -e "\n\n\n BEST THRESHOLD OBTAINED FROM SET_2 SINGLE DOMAIN: $set2_best_evalue_one_domain" >> hmm_results_strali.txt
    #SET_1 TEST
echo -e "\nPERFORMANCES SET_1 USING E-VALUE THRESHOLD OF SET_2 - SINGLE DOMAIN" >> hmm_results_strali.txt
python3 performance.py set_1_strali.class "$set2_best_evalue_one_domain" 2 >> hmm_results_strali.txt
        #false positives set_1
echo -e "False positives for set 1 considering single domain set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set2_best_evalue_one_domain" '$4 < num {print $1, $2, $4}' neg_1_strali.class | sort -grk 3 >> hmm_results_strali.txt
        #false negatives
echo -e "False negatives for set 1 considering single domain set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set2_best_evalue_one_domain" '$4 > num {print $1, $2, $4}' pos_1_strali.class | sort -grk 3 >> hmm_results_strali.txt
    #OVERALL
echo -e "\nOVERALL PERFORMANCES USING E-VALUE THRESHOLD OF SET_2 - SINGLE DOMAIN" >> hmm_results_strali.txt
python3 performance.py temp_overall_strali.class "$set2_best_evalue_one_domain" 2 >> hmm_results_strali.txt
        #false positives overall
echo -e "False positives for overall considering single domain set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set2_best_evalue_one_domain" '$4 < num  {print $1, $2, $4}' neg_2_strali.class neg_1_strali.class | sort -grk 3 >> hmm_results_strali.txt
        #false negatives overall
echo -e "False negatives for overall considering single domain set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_strali.txt
awk -v num="$set2_best_evalue_one_domain" '$4 > num {print $1, $2, $4}' pos_1_strali.class pos_2_strali.class | sort -grk 3 >> hmm_results_strali.txt





#CLEANING PROCESS
read -p "Do you want to delete temporary files? (y/n): " choice
if [[ "$choice" == "y" || "$choice" == "Y" ]]; then
    #CLEANING PROCEDURE
    #NOTE: This allows the user to delete temporary useless files. If the user has to do many trials it is suggested do not delete temporary files.
    echo "Cleaning temporary and intermediate files..."

    # Cluster-related temporary files
    rm -f pdb_kunitz_customreported.fasta
    rm -f pdb_kunitz_customreported.clstr
    rm -f pdb_kunitz_customreported.clstr.clstr
    rm -f pdb_kunitz_rp.ids
    rm -f pdb_kunitz_rp.fasta
    
    # All intermediate temporary files
    rm -f tmp*.*


    # BLAST database files created by makeblastdb (share prefix 'all_kunitz.pdb')
    rm -f all_kunitz.fasta.p* all_kunitz.fasta.*


    # .out found with hmmsearch
    rm -f *.out

    # Processed .fasta files (will be regenerated if needed)
    rm -f pos*.class neg*.class

    echo "Cleanup completed"
else
    echo "Temporary files kept."
fi
