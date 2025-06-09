# Building an HMM for the BPTI/Kunitz Protease Inhibitor Domain

This repository contains a pipeline to detect the BPTI/Kunitz-type protease inhibitor domain (Pfam: PF00014) using profile Hidden Markov Models (HMMs). Inspired by the work of Andrea Lenti, this project combines structural bioinformatics and sequence-based modeling to distinguish proteins containing the Kunitz domain from those that do not.

## 📌 Goal
To build and validate a **structure-guided HMM** and compare its performance with a **sequence-based HMM** constructed from multiple sequence alignment.

---

## 🔍 Why Use HMMs for Protein Domains?

Protein domains like **Kunitz** maintain structural and functional features even when their sequences diverge. Profile HMMs are ideal because they model:

- **Conserved residues (emissions)**
- **Insertions/deletions (gaps)**
- **Position-specific probabilities**

Training an HMM on structurally aligned domains helps detect homologs even with low sequence identity.

---

## 🧰 Required Tools

Install via Conda:

```bash
conda install -c bioconda cd-hit hmmer blast muscle
conda install -c conda-forge biopython
Tool	Purpose
CD-HIT	Redundancy removal in training sequences
HMMER	Build/search HMMs (hmmbuild, hmmsearch)
BLAST+	Similarity filtering to avoid sequence overlap
Biopython	FASTA file manipulation (get_seq.py)
Tool	Purpose
CD-HIT	Redundancy removal in training sequences
HMMER	Build/search HMMs (hmmbuild, hmmsearch)
BLAST+	Similarity filtering to avoid sequence overlap
Biopython	FASTA file manipulation (get_seq.py)
Pipeline Steps
0. Extract Representative Kunitz Sequences from PDB
Advanced PDB Search:
Resolution ≤ 3.5Å, Pfam ID = PF00014, Sequence length 45–80

Custom Report columns:

Entry ID

PDB ID

Entity ID

Auth Asym ID

Sequence

Annotation Identifier

Resolution

Run the following:

bash
Copy
Edit
bash script_recover_representative_kunitz.sh
This script:

Extracts PF00014 sequences from PDB

Clusters sequences at 90% identity using CD-HIT

Extracts representative sequences

Generates tmp_pdb_efold_ids.txt in PDB:CHAIN format

🔎 Before proceeding: Manually inspect tmp_pdb_efold_ids.txt and remove sequences that are too long.

1. Download Swiss-Prot Dataset
Go to UniProt and download the full Swiss-Prot FASTA file:
uniprot_sprot.fasta

2. Download All Kunitz Proteins
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

bash
Copy
Edit
python performance.py
Metrics computed:

MCC (Matthews Correlation Coefficient)

Q2 (Accuracy)

TPR (True Positive Rate)

PPV (Precision)
6. Plot Confusion Matrix
Use:

bash
Copy
Edit
python confusion_matrix.py
Edit the script to insert values for:

python
Copy
Edit
TP = ...
FP = ...
TN = ...
FN = ...
The script will generate .png confusion matrix plots.

📂 Output Summary
File	Description
hmm_results_strali.txt	Structure-based HMM evaluation report
hmm_results_seqali.txt	Sequence-based HMM evaluation report
neg_1.fasta, neg_2.fasta	Negative datasets
pos_1.fasta, pos_2.fasta	Positive datasets
pdb_kunitz_rp_clean.fasta	Cleaned domain sequences
pdb_kunitz_rp_seqali.hmm	HMM from MUSCLE alignment
pdb_kunitz_rp_strali.hmm	HMM from PDBeFold alignment
set_1_strali.class, set_2_strali.class	Classification results from HMMs
temp_overall.class	Combined classification results


Mohsin Nazir Bhat
MSc Bioinformatics – University of Bologna
Course: Laboratory of Bioinformatics 1 (LB)
