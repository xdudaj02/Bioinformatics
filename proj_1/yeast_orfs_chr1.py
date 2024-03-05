# File: yeast_orfs_chr1.py
# Author: Jakub Duda (202311235)
# Date: 5. 3. 2024
# Usage: python yeast_orfs_chr1.py <file_name_1.fasta> <file_name_2.gtf>

import sys

SEQUENCE_SIZE = 30_000
MIN_ORF_LENGTH = 150
START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]
CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '',  'TAG': '',
    'TGC': 'C', 'TGT': 'C', 'TGA': '',  'TGG': 'W',
}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


class Orf:
    '''
    A class representing an open reading frame (ORF).
    
    Attributes:
    - start: an integer representing the start position of the ORF
    - end: an integer representing the end position of the ORF
    - strand: a string representing the strand of the ORF
    - frame: an integer representing the frame of the ORF
    - sequence: a string representing the sequence of the ORF
    '''
    start: int
    end: int
    strand: str
    frame: int
    sequence: str

    def __init__(self, start: int, end: int, strand: str, frame: int, sequence: str) -> None:
        '''
        Initialize the ORF object.
        
        Args:
        - start: an integer representing the start position of the ORF
        - end: an integer representing the end position of the ORF
        - strand: a string representing the strand of the ORF
        - frame: an integer representing the frame of the ORF
        - sequence: a string representing the sequence of the ORF
        '''
        self.start = start
        self.end = end
        self.strand = strand
        self.frame = frame
        self.sequence = sequence

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

    def __init__(self, start: int, end: int, strand: str, gene_id: str) -> None:
        '''
        Initialize the Annotation object.
        
        Args:
        - start: an integer representing the start position of the annotated ORF
        - end: an integer representing the end position of the annotated ORF
        - strand: a string representing the strand of the annotated ORF
        - gene_id: a string representing the gene ID of the annotated ORF
        '''
        self.start = start
        self.end = end
        self.strand = strand
        self.gene_id = gene_id


def reverse_complement(sequence: str) -> str:
    '''
    Get the reverse complement of a sequence.
    
    Args:
    - sequence: a string representing the sequence
    
    Returns:
    - a string representing the reverse complement of the sequence
    '''
    return ''.join(COMPLEMENT[base] for base in sequence[::-1])

def translate_sequence(sequence: str) -> str:
    '''
    Translate a sequence to a protein sequence using the standard genetic code.
    
    Args:
    - sequence: a string representing the sequence
    
    Returns:
    - a string representing the protein sequence translated from the sequence
    '''
    return ''.join(CODON_TABLE[sequence[i : i + 3]] for i in range(0, len(sequence), 3))

def get_orfs(sequence: str, start_codon: str = START_CODON, stop_codons: list[str] = STOP_CODONS, min_length: int = MIN_ORF_LENGTH) -> list[Orf]:
    '''
    Find all ORFs in a sequence with a length greater than or equal to min_length. The ORFs are found on both strands and in all frames.
    
    Args:
    - sequence: a string representing the sequence
    - start_codon: a string representing the start codon (default: "ATG")
    - stop_codons: a list of strings representing the stop codons (default: ["TAA", "TAG", "TGA"])
    - min_length: an integer representing the minimum length of the ORF (default: 150)

    Returns:
    - a list of Orf objects representing the ORFs found in the sequence
    '''
    sequence_size = len(sequence)
    orfs = []
    for strand in ["+", "-"]:
        end_indices = set()
        sequence = sequence if strand == "+" else reverse_complement(sequence)
        for frame in range(3):
            frame_sequence = sequence[frame:]
            start_codon_occurrences = [i for i in range(len(frame_sequence) - 2) if frame_sequence[i : i + 3] == start_codon]
            for start in start_codon_occurrences:
                end = None
                for i in range(start + 3, len(frame_sequence), 3):
                    codon = frame_sequence[i:i + 3]
                    if codon in stop_codons:
                        end = i + 3
                        break
                if end:
                    orf_sequence = frame_sequence[start : end]
                    if len(orf_sequence) >= min_length and (end + frame) not in end_indices:
                        if strand == "+":
                            orfs.append(Orf(start + frame + 1, end + frame, strand, frame, orf_sequence))
                        else:
                            orfs.append(Orf(sequence_size - end + 1, sequence_size - start, strand, frame, orf_sequence))
                        end_indices.add(end + frame)

    return orfs

def get_genomic_sequence(fasta_filename: str, limit: int = SEQUENCE_SIZE) -> str:
    '''
    Read a genomic sequence from a FASTA file and return it as a string.
    
    Args:
    - fasta_filename: a string representing the name of the FASTA file
    - limit: an integer representing the maximum length of the sequence to read (default: 30,000)

    Returns:
    - a string representing the genomic sequence

    Raises:
    - SystemExit: if the file does not exist
    '''
    try:
        with open(fasta_filename) as f:
            sequence = ''.join(line.strip() for line in f.read().split('\n')[1:])
            return sequence[:limit]
    except FileNotFoundError:
        raise SystemExit(f"The file {fasta_filename} does not exist.")

def sequence_length(sequence: str) -> int:
    '''
    Get the length of a sequence.
    
    Args:
    - sequence: a string representing the sequence
    '''
    return len(sequence)

def frequency(sequence: str, base: str | list[str]) -> float:
    '''
    Get the frequency of a base in a sequence as a percentage.
    
    Args:
    - sequence: a string representing the sequence
    - base: a string representing of the base or a list of strings representing multiple bases

    Returns:
    - a float representing the frequency of the base in the sequence as a percentage
    '''
    base_count = sequence.count(base) if isinstance(base, str) else sum(sequence.count(b) for b in base)
    return (base_count / len(sequence)) * 100

def gc_content(sequence: str) -> float:
    '''
    Get the GC content of a sequence as a percentage.

    Args:
    - sequence: a string representing the sequence

    Returns:
    - a float representing the GC content of the sequence as a percentage
    '''
    return (frequency(sequence, ["G", "C"]))

def codon_counts(sequence: str) -> dict[str, int]:
    '''
    Get the count of all codons in a sequence.

    Args:
    - sequence: a string representing the sequence

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

def print_statistics(sequence: str, codon_counts: dict[str, int]) -> None:
    '''
    Print statistics about a sequence.
    
    Args:
    - sequence: a string representing the sequence
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

def output_to_files(orfs: list[Orf]) -> None:
    '''
    Write ORF information to files.
    
    Args:
    - orfs: a list of tuples (start, end, protein_sequence) representing the ORFs found in the sequence
    '''
    with open("all_potential_proteins.txt", "w") as protein_file, \
            open("orf_coordinates.txt", "w") as coordinates_file:
        for i, orf in enumerate(orfs, start=1):
            protein_file.write(f"{translate_sequence(orf.sequence)}\n")
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
            print(f"{annotation_id}\t{overlap:.0f}%")
        else:
            print(f"{annotation_id}\t0%")

def parse_args() -> tuple[str, str]:
    '''
    Parse command-line arguments.

    Returns:
    - a tuple of strings representing the names of the FASTA and GTF files

    Raises:
    - SystemExit: if the number of arguments is not 3 or if the first file is not a FASTA file or if the second file is not a GTF file
    '''
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python yeast_orfs_chr1.py <file_name_1.fasta> <file_name_2.gtf>.")
    if not sys.argv[1].endswith(".fasta"):
        raise SystemExit("The first file must be a FASTA file.")
    if not sys.argv[2].endswith(".gtf"):
        raise SystemExit("The second file must be a GTF file.")
    
    return sys.argv[1], sys.argv[2]


# Run the script
if __name__ == "__main__":
    fasta_filename, gtf_filename = parse_args()

    import time
    t0 = time.time()
    gen_seq = get_genomic_sequence(fasta_filename)
    cod_cnts = codon_counts(gen_seq)

    # tasks 1-6
    print_statistics(gen_seq, cod_cnts)

    orf_list = get_orfs(gen_seq)

    print()
    # tasks 7-8
    output_to_files(orf_list)
    # task 9
    print_overlaps(find_overlaps(get_annotations(gtf_filename), orf_list))
