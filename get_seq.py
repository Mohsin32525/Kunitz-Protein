from sys import argv,stdout
from Bio import SeqIO

"""
Extracts sequences from a UniProt-formatted FASTA file based on a list of UniProt accession IDs,
and prints them to standard output in FASTA format.

Usage:
    python script.py id_list.txt input_fasta.fasta > output.fasta

Arguments:
    argv[1] = Text file containing one UniProt accession ID per line (e.g., A0A1D0BND9)
    argv[2] = FASTA file with UniProt headers (e.g., >sp|A0A1D0BND9|Protein_Name ...)

Notes:
- The ID list must contain only plain accession IDs, without any ">" symbol.
- The FASTA file must follow UniProt header format:
    >sp|ACCESSION|ENTRY_NAME ...
  The script extracts the second field (ACCESSION) for matching.
- Matching sequences are printed to stdout that are redirected to a file using `>` in the shell.
"""




def get_ids(file):
    # Carica gli ID da cercare
    with open(to_get, 'r') as file:
        ids_to_get = {line.strip() for line in file if line.strip()}
    return ids_to_get


def get_seq(ids_to_get, datasetfile):
    with open(datasetfile, 'r') as file1:
        for prot in SeqIO.parse(file1, 'fasta'):
            # Estrai l’ID Swiss-Prot, cioè la seconda parte tra i "|"
            try:
                uniprot_id = prot.id.split('|')[1] if '|' in prot.id else prot.id
            except IndexError:
                continue

            if uniprot_id in ids_to_get:
                SeqIO.write(prot, stdout, 'fasta')
                
if __name__ == '__main__':
    to_get = argv[1]  # File containing the list of IDs to keep
    dataset = argv[2]  # FASTA file from which to extract the matching sequences
    pidlist = get_ids(to_get)
    get_seq(pidlist, dataset)