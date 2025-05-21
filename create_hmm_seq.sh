#WARNING! This script uses some files that are created by the creating_testing_sets.sh so you have to run before all the other scripts and this one if you want a performance for an HMM model using sequence alignment

# Make the multiple sequence alignment using the rapresentative pdb sequences cleaned from structural alignment (to balance the evaluation)
muscle -align pdb_kunitz_rp_clean.fasta -output pdb_kunitz_rp_seqali.fasta
#is the file aligned sequencially  with the most representatives proteins. Contains the training set to generate the HMM model.

# Build the HMM model from the aligned sequences
hmmbuild pdb_kunitz_rp_seqali.hmm pdb_kunitz_rp_seqali.fasta




#SEQUENCE ALIGNMENT HMM RESULTS
# Run hmmsearch on all positive and negative FASTA files and create tabular output. Generates .out files that contains the E-values computed on the HMM of the sequence alignment.
hmmsearch -Z 1000 --max --tblout pos_1_seqali.out pdb_kunitz_rp_seqali.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2_seqali.out pdb_kunitz_rp_seqali.hmm pos_2.fasta
hmmsearch -Z 1000 --max --tblout neg_1_seqali.out pdb_kunitz_rp_seqali.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2_seqali.out pdb_kunitz_rp_seqali.hmm neg_2.fasta
# Extract E-values and build .class files with: ID, label (1/0), full sequence E-value ($5), and best domain E-value ($8)
grep -v "^#" pos_1_seqali.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}' > pos_1_seqali.class
grep -v "^#" pos_2_seqali.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}' > pos_2_seqali.class
grep -v "^#" neg_1_seqali.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}' > neg_1_seqali.class
grep -v "^#" neg_2_seqali.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}' > neg_2_seqali.class

# Add true negatives not detected by HMM (missing from the .out), with a fake high E-value (10.0)
comm -23 <(awk '$2 == 0 {print $1}' set_1_strali.class | sort) <(cut -f 1 neg_1_seqali.class | sort) | awk '{print $1"\t"0"\t"10.0"\t"10.0}' >> neg_1_seqali.class
comm -23 <(awk '$2 == 0 {print $1}' set_2_strali.class | sort) <(cut -f 1 neg_2_seqali.class | sort) | awk '{print $1"\t"0"\t"10.0"\t"10.0}' >> neg_2_seqali.class

# Combine positive and negative sets into final training/testing files
cat pos_1_seqali.class neg_1_seqali.class > set_1_seqali.class
cat pos_2_seqali.class neg_2_seqali.class > set_2_seqali.class
cat set_1_seqali.class set_2_seqali.class > temp_overall_seqali.class

# Automatically determine the best thresholds (based on highest MCC)
# SET 1
set1_best_evalue_full_seq=$(
    for i in $(seq 1 9); do
        python3 performance.py set_1_seqali.class 1e-"$i"
    done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)
set1_best_evalue_one_domain=$(
    for i in $(seq 1 9); do
        python3 performance.py set_1_seqali.class 1e-"$i"
    done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)

# SET 2
set2_best_evalue_full_seq=$(
    for i in $(seq 1 9); do
        python3 performance.py set_2_seqali.class 1e-"$i"
    done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)
set2_best_evalue_one_domain=$(
    for i in $(seq 1 9); do
        python3 performance.py set_2_seqali.class 1e-"$i"
    done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)


# Run final performance evaluations using the best thresholds
# Identify false positives: non-Kunitz proteins with E-value below threshold (misclassified as positive)
# Identify false negatives: Kunitz proteins with E-value above threshold (misclassified as negative)
#SET_1 BEST THRESHOLD USING FULL SEQUENCE E-VALUE
echo -e "BEST THRESHOLD OBTAINED FROM SET_1 FULL SEQUENCE: $set1_best_evalue_full_seq" > hmm_results_seqali.txt
    #SET_2 TEST
echo -e "\nPERFORMANCES SET_2 USING E-VALUE THRESHOLD OF SET_1 - FULL SEQUENCES" >> hmm_results_seqali.txt
python3 performance.py set_2_seqali.class "$set1_best_evalue_full_seq" 1 >> hmm_results_seqali.txt
        #false positives set_2
echo -e "False positives for set 2 considering full sequence set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set1_best_evalue_full_seq" '$3 < num {print $1, $2, $3}' neg_2_seqali.class | sort -grk 3 >> hmm_results_seqali.txt
        #false negatives set_2
echo -e "False negatives for set 2 considering full sequence set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set1_best_evalue_full_seq" '$3 > num {print $1, $2, $3}' pos_2_seqali.class | sort -grk 3 >> hmm_results_seqali.txt
    #OVERALL TEST
echo -e "\nOVERALL PERFORMANCES USING E-VALUE THRESHOLD OF SET_1 - FULL SEQUENCES" >> hmm_results_seqali.txt
python3 performance.py temp_overall_seqali.class "$set1_best_evalue_full_seq" 1 >> hmm_results_seqali.txt
        #false positives overall
echo -e "False positives for overall considering full sequence set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set1_best_evalue_full_seq" '$3 < num  {print $1, $2, $3}' neg_2_seqali.class neg_1_seqali.class | sort -grk 3 >> hmm_results_seqali.txt
        #false negatives overall
echo -e "False negatives for overall considering full sequence set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set1_best_evalue_full_seq" '$3 > num {print $1, $2, $3}' pos_1_seqali.class pos_2_seqali.class | sort -grk 3 >> hmm_results_seqali.txt 

#SET_2 BEST THRESHOLD USING FULL SEQUENCE E-VALUE
echo -e "\n\n\nBEST THRESHOLD OBTAINED FROM SET_2 FULL SEQUENCE: $set2_best_evalue_full_seq" >> hmm_results_seqali.txt
    #SET_1 TEST
echo -e "\nPERFORMANCES SET_1 USING E-VALUE THRESHOLD OF SET_2 - FULL SEQUENCES" >> hmm_results_seqali.txt
python3 performance.py set_1_seqali.class "$set2_best_evalue_full_seq" 1 >> hmm_results_seqali.txt
        #false positives set_1
echo -e "False positives for set 1 considering full sequence set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set2_best_evalue_full_seq" '$3 < num {print $1, $2, $3}' neg_1_seqali.class | sort -grk 3 >> hmm_results_seqali.txt
        #false negatives set_1
echo -e "False negatives for set 1 considering full sequence set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set2_best_evalue_full_seq" '$3 > num {print $1, $2, $3}' pos_1_seqali.class | sort -grk 3 >> hmm_results_seqali.txt
    #OVERALL TEST
echo -e "\nOVERALL PERFORMANCES USING E-VALUE THRESHOLD OF SET_2 - FULL SEQUENCES" >> hmm_results_seqali.txt
python3 performance.py temp_overall_seqali.class "$set2_best_evalue_full_seq" 1 >> hmm_results_seqali.txt
        #false positives overall
echo -e "False positives for overall considering full sequence set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set2_best_evalue_full_seq" '$3 < num  {print $1, $2, $3}' neg_2_seqali.class neg_1_seqali.class | sort -grk 3 >> hmm_results_seqali.txt
        #false negatives overall
echo -e "False negatives for overall considering full sequence set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set2_best_evalue_full_seq" '$3 > num {print $1, $2, $3}' pos_1_seqali.class pos_2_seqali.class | sort -grk 3 >> hmm_results_seqali.txt

#SET_1 BEST THRESHOLD USING SINGLE DOMAIN E-VALUE
echo -e "\n\n\n BEST THRESHOLD OBTAINED FROM SET_1 SINGLE DOMAIN: $set1_best_evalue_one_domain" >> hmm_results_seqali.txt
    #SET_2 TEST
echo -e "\nPERFORMANCES SET_2 USING E-VALUE THRESHOLD OF SET_1 - SINGLE DOMAIN" >> hmm_results_seqali.txt
python3 performance.py set_2_seqali.class "$set1_best_evalue_one_domain" 2 >> hmm_results_seqali.txt
        #false positives set_2
echo -e "False positives for set 2 considering single domain set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set1_best_evalue_one_domain" '$4 < num {print $1, $2, $4}' neg_2_seqali.class | sort -grk 3 >> hmm_results_seqali.txt
        #false negatives set_2
echo -e "False negatives for set 2 considering single domain set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set1_best_evalue_one_domain" '$4 > num {print $1, $2, $4}' pos_2_seqali.class | sort -grk 3 >> hmm_results_seqali.txt
    #OVERALL
echo -e "\nOVERALL PERFORMANCES USING E-VALUE THRESHOLD OF SET_1 - SINGLE DOMAIN" >> hmm_results_seqali.txt
python3 performance.py temp_overall_seqali.class "$set1_best_evalue_one_domain" 2 >> hmm_results_seqali.txt
        #false positives overall
echo -e "False positives for overall considering single domain set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set1_best_evalue_one_domain" '$4 < num  {print $1, $2, $4}' neg_2_seqali.class neg_1_seqali.class | sort -grk 3 >> hmm_results_seqali.txt   
        #false negatives overall     
echo -e "False negatives for overall considering single domain set 1 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set1_best_evalue_one_domain" '$4 > num {print $1, $2, $4}' pos_1_seqali.class pos_2_seqali.class | sort -grk 3 >> hmm_results_seqali.txt

#SET_2 BEST THRESHOLD USING SINGLE DOMAIN E-VALUE
echo -e "\n\n\n BEST THRESHOLD OBTAINED FROM SET_2 SINGLE DOMAIN: $set2_best_evalue_one_domain" >> hmm_results_seqali.txt
    #SET_1 TEST
echo -e "\nPERFORMANCES SET_1 USING E-VALUE THRESHOLD OF SET_2 - SINGLE DOMAIN" >> hmm_results_seqali.txt
python3 performance.py set_1_seqali.class "$set2_best_evalue_one_domain" 2 >> hmm_results_seqali.txt
        #false positives set_1
echo -e "False positives for set 1 considering single domain set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set2_best_evalue_one_domain" '$4 < num {print $1, $2, $4}' neg_1_seqali.class | sort -grk 3 >> hmm_results_seqali.txt
        #false negatives
echo -e "False negatives for set 1 considering single domain set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set2_best_evalue_one_domain" '$4 > num {print $1, $2, $4}' pos_1_seqali.class | sort -grk 3 >> hmm_results_seqali.txt
    #OVERALL
echo -e "\nOVERALL PERFORMANCES USING E-VALUE THRESHOLD OF SET_2 - SINGLE DOMAIN" >> hmm_results_seqali.txt
python3 performance.py temp_overall_seqali.class "$set2_best_evalue_one_domain" 2 >> hmm_results_seqali.txt
        #false positives overall
echo -e "False positives for overall considering single domain set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set2_best_evalue_one_domain" '$4 < num  {print $1, $2, $4}' neg_2_seqali.class neg_1_seqali.class | sort -grk 3 >> hmm_results_seqali.txt
        #false negatives overall
echo -e "False negatives for overall considering single domain set 2 e-value threshold:\nUniprotId|True Class|E-value" >> hmm_results_seqali.txt
awk -v num="$set2_best_evalue_one_domain" '$4 > num {print $1, $2, $4}' pos_1_seqali.class pos_2_seqali.class | sort -grk 3 >> hmm_results_seqali.txt


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


