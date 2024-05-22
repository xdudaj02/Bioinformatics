import itertools
import pandas as pd
from Bio import SeqIO


OHE_DICT = {
    'a': [1, 0, 0, 0],
    'c': [0, 1, 0, 0],
    'g': [0, 0, 1, 0],
    't': [0, 0, 0, 1]
}
PROTEIN_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
KMER_SIZE = 2

def dna_ohe(seq):
    return [OHE_DICT[base] for base in seq]

def word_to_kmer(word, k = KMER_SIZE):
    kmers =[''.join(p) for p in itertools.product(PROTEIN_ALPHABET, repeat=k)]
    kmer_dict = {x: 0 for x in kmers} if kmers else {}
    for i in range(len(word) - k + 1):
        kmer = word[i:i + k]
        try:
            kmer_dict[kmer] += 1
        except KeyError:
            # ignore invalid k-mers (e.g. containing X, B, Z, J, O or U)
            continue
    return kmer_dict

def file_to_kmer_table(file_name):
    records = []
    with open(file_name, 'r', encoding='utf-8') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            records.append(record)
    return pd.DataFrame([word_to_kmer(str(record.seq), KMER_SIZE) for record in records])


if __name__ == '__main__':
    print(file_to_kmer_table('../proj_3/globin.fasta'))
