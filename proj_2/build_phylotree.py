'''Build a phylogenetic tree for a given protein sequence identifier'''

import sys
import time
from io import StringIO
from typing import Optional

import requests
from Bio.Blast import NCBIWWW
from Bio import AlignIO, SeqIO, SearchIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


UNIPROT_OUT_FN = "sequence.fasta"
BLAST_OUT_FN = "sequences_to_analyse.fasta"
MSA_OUT_FN = "alignment.txt"


def get_args():
    '''
    Get command line arguments.

    Returns:
        dict: Dictionary containing command line arguments

    Raises:
        SystemExit: If the number of command line arguments is not equal to 2
    '''
    if len(sys.argv) != 2:
        print("Usage: python build_phylotree.py <prot_seq_id>")
        sys.exit(1)
    # todo: validate prot_seq_id (or not maybe - too complicated)
    return {'prot_seq_id': sys.argv[1]}

def get_protein_seq(prot_seq_id: str, filename: str = UNIPROT_OUT_FN) -> SeqIO.SeqRecord:
    '''
    Fetch the protein sequence and save it to a file.

    Args:
        protein_seq_id (str): Identifier of the protein sequence
        filename (str): Name of the file to save the sequence to

    Returns:
        SeqRecord: Protein sequence

    Raises:
        ValueError: If the data for the protein sequence fails to be fetched
    '''
    url = f'https://rest.uniprot.org/uniprotkb/{prot_seq_id}.fasta'
    response = requests.get(url, timeout=100)
    if not response.ok:
        raise ValueError(f'Failed to fetch data for {prot_seq_id}')
    prot_seq = SeqIO.read(StringIO(response.text), "fasta")
    SeqIO.write(prot_seq, filename, "fasta")
    return prot_seq

def run_blast_search(
        seq: Optional[SeqIO.SeqRecord] = None,
        filename: str = BLAST_OUT_FN
    ) -> list[SeqIO.SeqRecord]:
    '''
    Run a BLAST search for the protein sequence.

    Args:
        sequence (SeqRecord): Protein sequence
        filename (str): Name of the file to save the BLAST search results to

    Returns:
        list[SeqRecord]: List of protein sequences representing the BLAST search results
    '''
    if not seq:
        with open(UNIPROT_OUT_FN, encoding="utf-8") as h:
            seq = next(SeqIO.parse(h, "fasta"))
    eq = 'NOT "Homo sapiens"[Organism] NOT "synthetic construct"[Organism]'
    with NCBIWWW.qblast("blastp", "nr", seq.seq, entrez_query=eq, hitlist_size=10) as result_handle:
        blast_results = SearchIO.parse(result_handle, "blast-xml")
        blast_records = []
        for query_result in blast_results:
            for hit in query_result.hits:
                for fragment in hit.fragments:
                    # todo: (maybe) change id to non-human for mulit-id ones
                    # todo: make sure to keep only one match for each organism
                    #           (need to query for more than 10)
                    blast_records.append(fragment.hit)
        SeqIO.write(blast_records, filename, "fasta")
    return blast_records

def perform_msa(
        sequences: Optional[list[SeqIO.SeqRecord]] = None,
        filename: str = MSA_OUT_FN
    ) -> str:
    '''
    Perform Multiple Sequence Alignment (MSA) using the Clustal Omega algorithm.

    Args:
        sequences (list[SeqRecord]): List of protein sequences to be aligned
        filename (str): Name of the file to save the MSA results to

    Returns:
        str: Result of the MSA

    Raises:
        ValueError: If the MSA fails to be performed or an error occurs when retrieving the results
    '''
    # todo: citation (https://europepmc.org/article/MED/35412617)
    if not sequences:
        with open(UNIPROT_OUT_FN, encoding="utf-8") as f:
            sequences = [next(SeqIO.parse(f, "fasta"))]
        with open(BLAST_OUT_FN, encoding="utf-8") as f:
            sequences.extend(SeqIO.parse(f, "fasta"))
    url = 'https://www.ebi.ac.uk/Tools/services/rest/clustalo'
    form_data = {
        'email': 'up202311235@edu.fc.up.pt',
        'outfmt': 'nexus',
        'sequence': '\n'.join([seq.format('fasta') for seq in sequences]),
        'stype': 'protein',
    }
    response = requests.post(f'{url}/run', data=form_data, timeout=100)
    if not response.ok:
        print(response.text)
        raise ValueError('Failed to perform MSA')
    job_id = response.text
    job_url = f'{url}/status/{job_id}'
    while True:
        response = requests.get(job_url, timeout=100)
        if not response.ok:
            raise ValueError('Failed to get MSA status')
        if response.text == 'FINISHED':
            break
        time.sleep(1)
    response = requests.get(f'{url}/result/{job_id}/aln-nexus', timeout=100)
    if not response.ok:
        raise ValueError('Failed to get MSA results')
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(response.text)

def build_phylotree(alignment_data: Optional[str] = None):
    '''
    Build a phylogenetic tree using the MSA results.

    Args:
        alignment_data (str): MSA results
    '''
    if not alignment_data:
        with open(MSA_OUT_FN, encoding='utf-8') as f:
            alignment_data = f.read()
    alignment_data = AlignIO.read(StringIO(alignment_data), "nexus")
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment_data)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    Phylo.draw(tree)
    # todo: different tree formats
    # todo: maybe different labels

def main():
    '''
    Main function.
    '''
    args = get_args()
    prot_seq_id = args['prot_seq_id']

    # Part 1
    get_protein_seq(prot_seq_id)

    # Part 2
    run_blast_search()

    # Part 3
    perform_msa()

    # Part 4
    build_phylotree()


if __name__ == "__main__":
    main()
