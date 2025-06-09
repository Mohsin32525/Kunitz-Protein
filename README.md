# Project Overview
This project is about building a pipeline (a step-by-step method) to detect the Kunitz-type protease inhibitor domain (Pfam ID: PF00014) in proteins using Profile Hidden Markov Models (HMMs).
Developed by Mohsin Nazir Bhat, and combines two main approaches:
Structural bioinformatics – using 3D structure information of proteins.
Sequence modeling – using the protein sequences directly.
The main goal is to:
Create an HMM based on 3D structure alignments (structure-guided HMM),
Compare it with a regular sequence-based HMM (made from a multiple sequence alignment),
And see which one works better at correctly detecting proteins that have the Kunitz domain and avoiding those that don't.
##  HMMs for Protein Domains
Protein domains like the Kunitz domain often keep their core structure and important functional parts, even in proteins that aren't very similar in sequence.
Profile Hidden Markov Models (HMMs) help us detect these domains by learning:
Which amino acids tend to stay the same (called conserved residues),
Where insertions or deletions (gaps) usually happen,
And the probabilities of each amino acid appearing at different positions.
By training an HMM using multiple aligned examples of the Kunitz domain (based on their 3D structures), we can build a model that can spot other similar proteins, even if their sequences look quite different.

## Required Tools
Install these via conda:

conda install -c bioconda cd-hit hmmer blast muscle
conda install -c conda-forge biopython

## Tool	 Use
 1. CD-HIT	Remove redundancy in training sequences
 2. HMMER 	Build and use HMMs (hmmbuild, hmmsearch)
 3. BLAST+	Identify and filter sequence similarity to prevent overlap
4.  Biopython	FASTA manipulation via get_seq.py
## Pipeline Steps
0. Extract Representative Kunitz Sequences from PDB

Use advanced search: Data Collection Resolution <= 3.5 AND ( Identifier = "PF00014" AND Annotation Type = "Pfam" ) AND Polymer Entity Sequence Length <= 80 AND Polymer 

Entity Sequence Length >= 45

Press the custom report with the following flags:

Entry ID

PDB ID
Entity ID
Auth Asym ID

Sequence

Annotation Identifier

Data collection resolution

This will output a .CSV file. Then execute the following bash script:
```bash
bash script_recover_representative_kunitz.sh
```
## This script:

Extracts PF00014 sequences from PDB

Clusters sequences at 90% identity using CD-HIT

Extracts representative sequences

Generates tmp_pdb_efold_ids.txt in PDB:CHAIN format

Before proceeding: Manually inspect tmp_pdb_efold_ids.txt and remove sequences that are too long.


1. Download Swiss-Prot Dataset
Go to UniProt and download the full Swiss-Prot FASTA file:
uniprot_sprot.fasta


 Download All Kunitz Proteins
From UniProt, search and download all Kunitz domain proteins:
all_kunitz.fasta

3. Structural Alignment
Use PDBeFold with tmp_pdb_efold_ids.txt
Output: pdb_kunitz_rp.ali

4. Build Structural HMM
bash
Copy
Edit
bash create_hmm_str.sh
bash create_testing_sets.sh
Builds structural HMM from PDBeFold alignment

Removes training sequences from test set

Generates balanced positive and negative test sets

Finds best E-value thresholds using MCC (2-fold CV)

Outputs detailed performance to hmm_results_strali.txt

5. Evaluate Performance
   Run:

```bash
python performance.py
```
Metrics computed:

MCC (Matthews Correlation Coefficient)

Q2 (Accuracy)

TPR (True Positive Rate)

PPV (Precision)

Plot Confusion Matrix
```bash
python confusion_matrix.py
```

## Output Files Summary
hmm_results_strali.txt and hmm_results_seqali.txt contain:
The best E-value thresholds selected by maximizing the Matthews Correlation Coefficient (MCC) -Performance metrics for each test set and overall, calculated using the E-value that yielded the highest MCC, based on either full sequence or best single domain evaluations -Lists of false positives and false negatives
neg_1.fasta and neg_2.fasta: FASTA files of non-Kunitz sequences used respectively as the negative set 1 and set 2.
pos_1.fasta and pos_2.fasta: FASTA files of Kunitz (positive) sequences used in set 1 and set 2.
pdb_kunitz_rp_clean.fasta: Cleaned representative Kunitz sequences used for both structural and sequence alignments (after filtering for length).
pdb_kunitz_rp_seqali.fasta and pdb_kunitz_rp_seqali.hmm: The multiple sequence alignment and resulting HMM model produced using MUSCLE.
pdb_kunitz_rp_strali.fasta and pdb_kunitz_rp_strali.hmm: The multiple structure alignment and resulting HMM model built from PDBeFold.
set_1_strali.class and set_2_strali.class: Classification results from hmmsearch for set 1 and set 2 (structure-based model).
temp_overall.class: Classification results for the combined dataset (positive + negative), used to estimate overall performance for each threshold.
 Author

## Mohsin Nazir Bhat 
 MSc Bioinformatics – University of Bologna
Laboratory of Bioinformatics 1 (LB1)

You're welcome to suggest improvements or fork this project to extend its functionality!
