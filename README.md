# Project Overview
This project is about building a pipeline (a step-by-step method) to detect the Kunitz-type protease inhibitor domain (Pfam ID: PF00014) in proteins using Profile Hidden Markov Models (HMMs).
Developed by Andrea Lenti, and combines two main approaches:
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

(bash script_recover_representative_kunitz.sh)

