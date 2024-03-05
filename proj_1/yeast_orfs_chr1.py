# TODO: add description, authors, date, etc.
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]

# Function to find ORFs in a sequence
def find_orfs(sequence: Seq, start_codon: str = START_CODON, stop_codons: list[str] = STOP_CODONS) -> list[tuple[int, int, str]]:
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
                orfs.append((start + frame, end + frame, str(orf_sequence.translate())))
    return orfs

# Read a genomic sequence from a FASTA file
def get_genomic_sequence(fasta_filename: str, limit: int = 30_000) -> Seq:
    record = SeqIO.read(fasta_filename, "fasta")
    return record.seq[:limit]

# Get length of a sequence
def sequence_length(sequence: Seq):
    return len(sequence)

# Get frequency (in %) of a base
def frequency(sequence: Seq, base: str | list[str]) -> float:
    base_count = sequence.count(base) if isinstance(base, str) else sum(sequence.count(b) for b in base)
    return (base_count / len(sequence)) * 100

# Get GC content
def gc_content(sequence: Seq) -> float:
    return (frequency(sequence, ["G", "C"]))

# Get count of all codons
def codon_counts(sequence: Seq) -> dict[str, int]:
    codon_counts_dict = {}
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i : i + 3]
        codon_counts_dict[codon] = codon_counts_dict.get(codon, 0) + 1
    return codon_counts_dict

# Get least and most frequent codons (may be multiple)
def extremal_codons(codon_counts: dict[str, int], extremal_func: callable) -> list[str]:
    extremal_count = extremal_func(codon_counts.values())
    return [str(codon) for codon, count in codon_counts.items() if count == extremal_count]

# Print the statistics
def print_statistics(sequence: Seq, codon_counts: dict[str, int]) -> None:
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

# Find ORFs in both strands
def get_orfs(sequence: Seq) -> list[tuple[int, int, str]]:
    return find_orfs(sequence) + find_orfs(sequence.reverse_complement())

# Write ORF information to files
def output_to_files(orfs: list[tuple[int, int, str]]):
    with open("all_potential_proteins.txt", "w") as protein_file, \
            open("orf_coordinates.txt", "w") as coordinates_file:
        for i, (start, end, orf) in enumerate(orfs, start=1):
            protein_file.write(f"{orf}\n")
            coordinates_file.write(f"{start}, {end}, ORF{i}\n")

# Get annotations from a GTF file
def get_annotations(gtf_filename: str) -> list[tuple[int, int, str]]:
    annotations = []
    with open(gtf_filename) as gtf_file:
        for line in gtf_file:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                if fields[2] == "exon":
                    annotations.append((int(fields[3]), int(fields[4]), fields[6], fields[8].split("\"")[1]))
    return annotations

# Calculate overlap between an ORF and an annotation
def get_overlap(orf: tuple[int, int, str], annotation: tuple[int, int, str, str]) -> float:
    orf_start, orf_end, orf_seq = orf
    annot_start, annot_end, *_ = annotation
    if orf_start > annot_end or orf_end < annot_start:
        return 0
    overlap_start = max(orf_start, annot_start)
    overlap_end = min(orf_end, annot_end)
    return (overlap_end - overlap_start) / (annot_end - annot_start) * 100

# Get the best overlap for a given annotation
def get_best_overlap(annotation: tuple[int, int, str, str], orfs: list[tuple[int, int, str]]) -> tuple[tuple[int, int, str, str], tuple[int, int, str], float]:
    best_orf = None
    best_overlap = 0
    for orf in orfs:
        overlap = get_overlap(orf, annotation)
        if overlap > best_overlap:
            best_overlap = overlap
            best_orf = orf
    return (annotation, best_orf, best_overlap)

# Find the best ORF overlap for each annotation
def find_overlaps(annotations: list[tuple[int, int, str, str]], orfs: list[tuple[int, int, str]]) -> list[tuple[tuple[int, int, str, str], tuple[int, int, str], float]]:
    return [get_best_overlap(annot, orfs) for annot in annotations]

# Print overlaps
def print_overlaps(overlaps: list[tuple[tuple[int, int, str, str], tuple[int, int, str], float]]):
    print("9. Overlaps with the annotated ORFs:")
    for annotation, orf, overlap in overlaps:
        if orf:
            print(f"{annotation[3]}\t{overlap:.2f}%")
        else:
            print(f"{annotation[3]}\tNo overlap")

# Parse command-line arguments
def parse_args() -> tuple[str, str]:
    if len(sys.argv) != 3:
        print("Usage: python yeast_orfs_chr1.py <file_name_1.fasta> <file_name_2.gtf>", file=sys.stderr)
        sys.exit(1)
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
