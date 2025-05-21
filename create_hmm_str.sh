# Convert the PDBefold alignment to FASTA format (uppercase, one-line), removing 'PDB:' prefix
awk '{if (substr($1,1,1)==">") {print "\n"toupper($1)} else {printf "%s",toupper($1)}}' pdb_kunitz_rp.ali | sed s/PDB:// | tail -n +2 > pdb_kunitz_rp_strali.fasta 
#is the file aligned structurally and with the most representatives. Contains the training set to generate the HMM model.

# Build the HMM model from the aligned sequences
hmmbuild pdb_kunitz_rp_strali.hmm pdb_kunitz_rp_strali.sto


