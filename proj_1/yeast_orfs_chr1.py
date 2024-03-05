# File: yeast_orfs_chr1.py
# Author: Jakub Duda (202311235)
# Date: 5. 3. 2024
# Description: This script reads a genomic sequence from a FASTA file and a GTF file with annotations. It then calculates various statistics about the sequence, including its length, frequency of each base, GC content, number of start and stop codons, and the most and least frequent codons. Then it finds all open reading frames (ORFs) in the sequence and writes them to a file, along with their coordinates. It also finds all potential proteins and writes them to a file. Finally, it finds the percentage of the biggest overlap between each annotated ORF and any of the ORFs found in the sequence and prints the results.
# Usage: python yeast_orfs_chr1.py <file_name_1.fasta> <file_name_2.gtf>

import sys
from dataclasses import dataclass

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]

@dataclass
class Orf:
    '''
    A class representing an open reading frame (ORF).
    
    Attributes:
    - start: an integer representing the start position of the ORF
    - end: an integer representing the end position of the ORF
    - protein_sequence: a string representing the protein sequence translated from the ORF
    '''
    start: int
    end: int
    protein_sequence: str

@dataclass
class Annotation:
    '''
    A class representing an annotation of an ORF from a GTF file.
    
    Attributes:
    - start: an integer representing the start position of the annotated ORF
    - end: an integer representing the end position of the annotated ORF
    - strand: a string representing the strand of the annotated ORF
    - gene_id: a string representing the gene ID of the annotated ORF
    '''
    start: int
    end: int
    strand: str
    gene_id: str


def find_orfs(sequence: Seq, start_codon: str = START_CODON, stop_codons: list[str] = STOP_CODONS) -> list[Orf]:
    '''
    Find all open reading frames (ORFs) in a sequence and return them as a list of tuples (start, end, protein_sequence).
    
    Args:
    - sequence: a Bio.Seq object representing the sequence
    - start_codon: a string representing the start codon (default: "ATG")
    - stop_codons: a list of strings representing the stop codons (default: ["TAA", "TAG", "TGA"])

    Returns:
    - a list of Orf objects representing the ORFs found in the sequence
    '''
    orfs = []
    for frame in range(3):
        frame_sequence = sequence[frame:]
        for start in nt_search(str(frame_sequence), start_codon)[1:]:
            end = None
            for i in range(start + 3, len(frame_sequence), 3):
                codon = frame_sequence[i:i + 3]
                if codon in stop_codons:
                    end = i + 3
                    break
            if end:
                orf_sequence = frame_sequence[start : end]
                orfs.append(Orf(start + frame, end + frame, str(orf_sequence.translate())))
    return orfs

def get_genomic_sequence(fasta_filename: str, limit: int = 30_000) -> Seq:
    '''
    Read a genomic sequence from a FASTA file and return it as a Bio.Seq object.
    
    Args:
    - fasta_filename: a string representing the name of the FASTA file
    - limit: an integer representing the maximum length of the sequence to read (default: 30,000)

    Returns:
    - a Bio.Seq object representing the genomic sequence
    '''
    record = SeqIO.read(fasta_filename, "fasta")
    return record.seq[:limit]

def sequence_length(sequence: Seq):
    '''
    Get the length of a sequence.
    
    Args:
    - sequence: a Bio.Seq object representing the sequence
    '''
    return len(sequence)

def frequency(sequence: Seq, base: str | list[str]) -> float:
    '''
    Get the frequency of a base in a sequence as a percentage.
    
    Args:
    - sequence: a Bio.Seq object representing the sequence
    - base: a string representing of the base or a list of strings representing multiple bases

    Returns:
    - a float representing the frequency of the base in the sequence as a percentage
    '''
    base_count = sequence.count(base) if isinstance(base, str) else sum(sequence.count(b) for b in base)
    return (base_count / len(sequence)) * 100

def gc_content(sequence: Seq) -> float:
    '''
    Get the GC content of a sequence as a percentage.

    Args:
    - sequence: a Bio.Seq object representing the sequence

    Returns:
    - a float representing the GC content of the sequence as a percentage
    '''
    return (frequency(sequence, ["G", "C"]))

def codon_counts(sequence: Seq) -> dict[str, int]:
    '''
    Get the count of all codons in a sequence.

    Args:
    - sequence: a Bio.Seq object representing the sequence

    Returns:
    - a dictionary with keys representing codons and values representing their counts
    '''
    codon_counts_dict = {}
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i : i + 3]
        codon_counts_dict[codon] = codon_counts_dict.get(codon, 0) + 1
    return codon_counts_dict

def extremal_codons(codon_counts: dict[str, int], extremal_func: callable) -> list[str]:
    '''
    Get the least or most frequent codons in a sequence.

    Args:
    - codon_counts: a dictionary with keys representing codons and values representing their counts
    - extremal_func: a callable representing the function to find the least or most frequent codons (e.g. min or max)

    Returns:
    - a list of strings representing the least or most frequent codons
    '''
    extremal_count = extremal_func(codon_counts.values())
    return [str(codon) for codon, count in codon_counts.items() if count == extremal_count]

def print_statistics(sequence: Seq, codon_counts: dict[str, int]) -> None:
    '''
    Print statistics about a sequence.
    
    Args:
    - sequence: a Bio.Seq object representing the sequence
    - codon_counts: a dictionary with keys representing codons and values representing their counts
    
    Returns:
    - None
    '''
    print(f"1. Length of the sequence: {sequence_length(sequence)}")
    print(f"2. Frequency (in %) of A: {frequency(sequence, 'A'):.2f}%, ",
                                 f"G: {frequency(sequence, 'G'):.2f}%, ",
                                 f"C: {frequency(sequence, 'C'):.2f}%, ",
                                 f"T: {frequency(sequence, 'T'):.2f}%")
    print(f"3. GC content: {gc_content(sequence):.2f}%")
    print(f"4. Number of Start (AUG) Codons: {codon_counts.get(START_CODON, 0)}")
    print(f"5. Number of Stop Codons (UAA, UAG, UGA): {sum(codon_counts.get(codon, 0) for codon in STOP_CODONS)}")
    print(f"6. Most frequent codon(s): {', '.join(extremal_codons(codon_counts, max))} ",
             f"Least frequent codon(s): {', '.join(extremal_codons(codon_counts, min))}")

def get_orfs(sequence: Seq) -> list[Orf]:
    '''
    Find all open reading frames (ORFs) in a sequence and its reverse complement and return them as a list of tuples (start, end, protein_sequence).
    
    Args:
    - sequence: a Bio.Seq object representing the sequence
    
    Returns:
    - a list of Orf objects representing the ORFs found in the sequence
    '''
    return find_orfs(sequence) + find_orfs(sequence.reverse_complement())

def output_to_files(orfs: list[Orf]) -> None:
    '''
    Write ORF information to files.
    
    Args:
    - orfs: a list of tuples (start, end, protein_sequence) representing the ORFs found in the sequence
    '''
    with open("all_potential_proteins.txt", "w") as protein_file, \
            open("orf_coordinates.txt", "w") as coordinates_file:
        for i, orf in enumerate(orfs, start=1):
            protein_file.write(f"{orf.protein_sequence}\n")
            coordinates_file.write(f"{orf.start}, {orf.end}, ORF{i}\n")

def get_annotations(gtf_filename: str) -> list[Annotation]:
    '''
    Get annotations from a GTF file and return them as a list of tuples (start, end, strand, gene_id).
    
    Args:
    - gtf_filename: a string representing the name of the GTF file
    
    
    Returns:
    - a list of Annotation objects representing the annotations found in the GTF file
    '''
    annotations = []
    with open(gtf_filename) as gtf_file:
        for line in gtf_file:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                if fields[2] == "exon":
                    annotations.append(Annotation(int(fields[3]), int(fields[4]), fields[6], fields[8].split("\"")[1]))
    return annotations

def get_overlap(orf: Orf, annotation: Annotation) -> float:
    '''
    Calculate the percentage of overlap between an annotated ORF and an ORF found in the sequence.
    
    Args:
    - orf: an Orf object representing the ORF found in the sequence
    - annotation: an Annotation object representing the annotated ORF
    
    Returns:
    - a float representing the percentage of overlap between the two ORFs
    '''
    if orf.start > annotation.end or orf.end < annotation.start:
        return 0
    overlap_start = max(orf.start, annotation.start)
    overlap_end = min(orf.end, annotation.end)
    return (overlap_end - overlap_start) / (annotation.end - annotation.start) * 100

def get_best_overlap(annotation: Annotation, orfs: list[Orf]) -> tuple[Orf, float]:
    '''
    Get the best overlap for a given annotation.
    
    Args:
    - annotation: an Annotation object representing the annotated ORF
    - orfs: a list of Orf objects representing the ORFs found in the sequence
    
    Returns:
    - a tuple of an Orf object representing the best ORF overlap for the annotation and a float representing the percentage of overlap
    '''
    best_orf = None
    best_overlap = 0
    for orf in orfs:
        overlap = get_overlap(orf, annotation)
        if overlap > best_overlap:
            best_overlap = overlap
            best_orf = orf
    return (best_orf, best_overlap)

def find_overlaps(annotations: list[Annotation], orfs: list[Orf]) -> dict[str, tuple[Orf, float]]:
    '''
    Find the best ORF overlap for each annotation.
    
    Args:
    - annotations: a list of Annotation objects representing the annotations found in the GTF file
    - orfs: a list of Orf objects representing the ORFs found in the sequence
    
    Returns:
    - a dictionary with keys representing annotations by ids and values representing the best ORF overlap for each annotation
    '''
    return {annotation.gene_id: get_best_overlap(annotation, orfs) for annotation in annotations}

def print_overlaps(overlaps: dict[str, tuple[Orf, float]]) -> None:
    '''
    Print overlaps between the annotated ORFs and the ORFs found in the sequence.
    
    Args:
    - overlaps: a dictionary with keys representing annotations by ids and values representing the best ORF overlap for each annotation

    Returns:
    - None
    '''
    print("9. Overlaps with the annotated ORFs:")
    for annotation_id, (orf, overlap) in overlaps.items():
        if orf:
            print(f"{annotation_id}\t{overlap:.2f}%")
        else:
            print(f"{annotation_id}\tNo overlap")

def parse_args() -> tuple[str, str]:
    '''
    Parse command-line arguments.

    Returns:
    - a tuple of strings representing the names of the FASTA and GTF files

    Raises:
    - SystemExit: if the number of arguments is not 3 or if the first file is not a FASTA file or if the second file is not a GTF file
    '''
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python yeast_orfs_chr1.py <file_name_1.fasta> <file_name_2.gtf>")
    if not sys.argv[1].endswith(".fasta"):
        raise SystemExit("The first file must be a FASTA file")
    if not sys.argv[2].endswith(".gtf"):
        raise SystemExit("The second file must be a GTF file")
    
    return sys.argv[1], sys.argv[2]


# Run the script
if __name__ == "__main__":
    fasta_filename, gtf_filename = parse_args()

    gen_seq = get_genomic_sequence(fasta_filename)
    cod_cnts = codon_counts(gen_seq)

    print_statistics(gen_seq, cod_cnts)

    orf_list = get_orfs(gen_seq)

    print()
    output_to_files(orf_list)
    print_overlaps(find_overlaps(get_annotations(gtf_filename), orf_list))
