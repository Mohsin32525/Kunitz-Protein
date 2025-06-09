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

bash
conda install -c bioconda cd-hit hmmer blast muscle
conda install -c conda-forge biopython
Tool	Purpose
CD-HIT  	  Redundancy removal in training sequences
HMMER      	Build/search HMMs (hmmbuild, hmmsearch)
BLAST+	    Similarity filtering to avoid sequence overlap
Biopython 	FASTA file manipulation (get_seq.py)
