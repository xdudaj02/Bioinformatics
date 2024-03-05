# Bioinformatics - Group Assignment 1
Python script for the first group assignment for the Bioinformatics cours at FCUP 23/24. Solves all the tasks listed in the [assignment](group_assignment_1.pdf). 

## Authors
- Jakub Duda (202311235)

## Usage
Script is written in *Python 3.11.0* and requires the *Biopython* library to be installed.

To run the script, use the following command:
```bash
python yeast_orfs_chr1.py file_1.fasta file_2.gtf
```
, where 
- `file_1.fasta` is a file in the FASTA format containing the sequence of the yeast genome and
- `file_2.gtf` is a file in the GTF format containing annotations of coordinates of ORFs in the yeast genome.
Both arguments are required.

## Output
The script produces output to the standard output and to two files: `all_potential_proteins.txt` and `orf_coordinates.txt`.

#### Standard output
First the script prints to the standard output basic statistics about the input sequence and its codons. Then it prints all the ORFs from the annotations file and the percentage of maximum overlap with any ORF found in the sequence. 

Example output:
```txt
1. Length of the sequence: 30000
2. Frequency (in %) of A: 31.18%,  G: 19.35%,  C: 17.60%,  T: 31.87%
3. GC content: 36.95%
4. Number of Start (AUG) Codons: 228
5. Number of Stop Codons (UAA, UAG, UGA): 560
6. Most frequent codon(s): TTT  Least frequent codon(s): CGC

9. Overlaps with the annotated ORFs:
NM_001180043.1  37.85%
NM_001184582.1  100.00%
NM_001178208.1  16.34%
NM_001179897.1  28.76%
NM_001180042.1  100.00%
NM_001180041.1  44.74%
NM_001178206.2  100.00%
NM_001184642.1  100.00%
```

#### Output files
The file `all_potential_proteins.txt` contains all the potential proteins found in the sequence printed on separate lines. Each line contains one protein sequence. Stop codons are represented by `*`.

The file `orf_coordinates.txt` contains all the ORFs found in the sequence printed on separate lines. Each line contains one ORF and has the following format:
`Start_X, End_X, ORFX`, where **X** represents the index of the ORF and **Start_X** and **End_X** are the start and end coordinates of this **X**-th ORF.
