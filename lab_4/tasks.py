import sys
import re

from statistics import multimode


### CLASSES ###
class Sequence:
    full_id: str
    db: str
    uid: str
    entry_name: str
    protein_name: str
    organism_name: str # OS
    organism_id: str # OX
    gene_name: str | None # GN
    protein_existence: str # PE
    sequence_version: str # SV
    sequence: str

    def __init__(self, full_id: str, db: str, uid: str, entry_name: str, protein_name: str, organism_name: str, organism_id: str, gene_name: str | None, protein_existence: str, sequence_version: str, sequence: str):
        self.full_id = full_id
        self.db = db
        self.uid = uid
        self.entry_name = entry_name
        self.protein_name = protein_name
        self.organism_name = organism_name
        self.organism_id = organism_id
        self.gene_name = gene_name
        self.protein_existence = protein_existence
        self.sequence_version = sequence_version
        self.sequence = sequence

    @classmethod
    def from_fasta(cls, header: str, sequence: str) -> 'Sequence':
        header_pattern = r'>((.{2})\|([^|]+)\|(.*?)) (.*?) OS=(.*?) OX=(.*?) (?:GN=(.*?) )?PE=(.*?) SV=(\d+)'
        match = re.match(header_pattern, header)
        if not match:
            raise ValueError('Invalid FASTA header')
        return cls(*match.groups(), sequence)


### UTILITY FUNCTIONS ###
def parse_args() -> tuple[int, str]:
    if len(sys.argv) != 3:
        raise SystemExit('Usage: python3 task_3.py <task_no> <fasta_filename>')
    if not sys.argv[1] in ['2', '3', '4', '5', '6']:
        raise SystemExit('Task number must be 2, 3, 4, 5 or 6.')
    return (int(sys.argv[1]), sys.argv[2])

def read_fasta_file(filename: str) -> list[Sequence]:
    sequences = []
    header = ''
    sequence = ''
    with open(filename, 'r') as file:
        while (line := file.readline()):
            if line.startswith('>') and sequence:
                header = line
                sequences.append(Sequence.from_fasta(header, sequence))
                sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(Sequence.from_fasta(header, sequence))
    return sequences

def prosite_to_regex(prosite_pattern: str) -> str:
    prosite_pattern = prosite_pattern.replace('-', '')
    prosite_pattern = prosite_pattern.replace('x', '.')
    prosite_pattern = prosite_pattern.replace('(', '{')
    prosite_pattern = prosite_pattern.replace(')', '}')
    prosite_pattern = f'({prosite_pattern})'
    return prosite_pattern


### TASKS ###
def task_2(filename: str):
    # filename = 'Q8RXD4.fasta'
    sequence = next(iter(read_fasta_file(filename).values()))
    pattern = prosite_to_regex('C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]')
    sequence = re.sub(pattern.lower(), lambda x: x.group(1).upper(), sequence.lower())
    print(sequence)

def task_3(filename: str):
    # filename = 'PS00727.fasta'
    sequences = read_fasta_file(filename)
    pattern = prosite_to_regex('N-x-G-x-R-[LIVM]-D-[LIVMFYH]-x-[LV]-x-S')
    for seq in sequences:
        not_str = '' if bool(re.search(pattern, seq.sequence)) else 'NOT '
        print(f'{seq.full_id} {not_str}MATCH')

def task_4(filename: str):
    # filename = 'PS00727.fasta'
    sequences = read_fasta_file(filename)
    sequence_id_list = [seq.full_id for seq in sequences]
    most_frequent_species = multimode([seq.organism_name for seq in sequences])
    sequence_dict = {seq.uid: (seq.organism_name, seq.organism_id, seq.gene_name, seq.protein_existence, seq.sequence_version) for seq in sequences}
    print(f'1) Sequence IDs: {sequence_id_list}')
    print(f'2) Most frequent species: {most_frequent_species}')
    print(f'3) Sequence info:')
    for uid, info in sequence_dict.items():
        print(f'{uid}: {info}')

def task_5(filename: str):
    # filename = 'sequence.gb'
    with open(filename, 'r') as file:
        content = file.read()
    pubmed_ids = re.findall(r'PUBMED\s+(\d+)', content)
    gene_data = re.findall(r'CDS\s+\d+..\d+\n\s+/gene=\"([^\"]+)\"\n\s+(?:/[a-z_]+=[^\n]+\n\s+)*?/translation=\"([^\"]+)\"', content)
    gene_data = [(gene[0], gene[1].replace('\n', '')) for gene in gene_data]
    gene_data = [(gene[0], gene[1].replace(' ', '')) for gene in gene_data]
    genes = [gene[0] for gene in gene_data]
    protein_ids = re.findall(r'/protein_id="(.*?)"', content)
    translated_sequences = {gene[0]: gene[1] for gene in gene_data}
    print(f'1) PubMed IDs: {pubmed_ids}')
    print(f'2) Genes: {genes}')
    print(f'3) Protein IDs: {protein_ids}')
    print(f'4) Translated sequences:')
    for gene, sequence in translated_sequences.items():
        print(f'{gene}: {sequence}')

def task_6(filename: str):
    # filename = 'genome_unformatted.fas'
    with open(filename, 'r') as file:
        content = file.read()
    content = re.sub(r'[\s\d]', '', content)
    print(content)


if __name__ == '__main__':
    task_no, filename = parse_args()

    match task_no:
        case 2:
            task_2(filename)
        case 3:
            task_3(filename)
        case 4:
            task_4(filename)
        case 5:
            task_5(filename)
        case 6:
            task_6(filename)
