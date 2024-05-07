"""
Filename: build_phylotree.py

This script implements the phylogenetic pipeline. For a given protein sequence identifier
(UniProt ID), BLAST search, Multiple Sequence Alignment (MSA) and finally phylogenetic tree
construction are all performed consecutively.

Authors: (todo)
    - Mariana surname <email>
    - Jakub Duda <up202311235@edu.fc.up.pt>
    - Soul surname <email>

Citations:
    - Madeira F, Pearce M, Tivey ARN, et al. Search and sequence analysis tools services from
    EMBL-EBI in 2022. Nucleic Acids Research. 2022 Jul;50(W1):W276-W279. DOI: 10.1093/nar/gkac240.
    PMID: 35412617; PMCID: PMC9252731.

Usage:
    python build_phylotree.py <prot_seq_id>
"""

import re
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

ALN_DATA_FORMAT = "nexus"


def get_uniprot_id() -> str:
    '''
    Validate and get the protein sequence identifier from the command line arguments.

    Returns:
        str: Protein sequence identifier

    Raises:
        SystemExit: If the number of command line arguments is not equal to 2
            or if the protein sequence identifier is invalid
    '''
    if len(sys.argv) != 2:
        print("Usage: python build_phylotree.py <prot_seq_id>")
        sys.exit(1)
    uniprot_id = sys.argv[1]
    id_regex = r'^(?:[OPQ]\d[\dA-Z]{3}\d)|(?:[A-NR-Z]\d[A-Z][\dA-Z]{2}\d(?:[A-Z][A-Z\d]{2}\d)?)$'
    if not re.match(id_regex, uniprot_id):
        print("Invalid protein sequence identifier")
        sys.exit(1)
    return uniprot_id

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
        aln_data_format: str = ALN_DATA_FORMAT,
        filename: str = MSA_OUT_FN
    ) -> str:
    '''
    Perform Multiple Sequence Alignment (MSA) using the Clustal Omega algorithm.

    Args:
        sequences (list[SeqRecord]): List of protein sequences to be aligned
        aln_data_format (str): Format of the MSA results
        filename (str): Name of the file to save the MSA results to

    Returns:
        str: Result of the MSA

    Raises:
        ValueError: If the MSA fails to be performed or an error occurs when retrieving the results
    '''
    if not sequences:
        with open(UNIPROT_OUT_FN, encoding="utf-8") as f:
            sequences = [next(SeqIO.parse(f, "fasta"))]
        with open(BLAST_OUT_FN, encoding="utf-8") as f:
            sequences.extend(SeqIO.parse(f, "fasta"))
    url = 'https://www.ebi.ac.uk/Tools/services/rest/clustalo'
    form_data = {
        'email': 'up202311235@edu.fc.up.pt',
        'outfmt': aln_data_format,
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
    response = requests.get(f'{url}/result/{job_id}/aln-{aln_data_format}', timeout=100)
    if not response.ok:
        raise ValueError('Failed to get MSA results')
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(response.text)

def build_phylotree(alignment_data: Optional[str] = None, aln_data_format: str = ALN_DATA_FORMAT):
    '''
    Build a phylogenetic tree using the MSA results.

    Args:
        alignment_data (str): MSA results
        aln_data_format (str): Format of the MSA results
    '''
    if not alignment_data:
        with open(MSA_OUT_FN, encoding='utf-8') as f:
            alignment_data = f.read()
    alignment_data = AlignIO.read(StringIO(alignment_data), aln_data_format)
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
    prot_seq_id = get_uniprot_id()

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
