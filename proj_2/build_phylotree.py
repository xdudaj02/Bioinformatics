"""
Filename: build_phylotree.py

This script implements the phylogenetic pipeline. For a given protein sequence identifier
(UniProt ID), BLAST search, Multiple Sequence Alignment (MSA) and finally phylogenetic tree
construction are all performed consecutively.

Authors:
    - Mariana Lourenço <up201906985@edu.fc.up.pt>
    - Jakub Duda <up202311235@edu.fc.up.pt>
    - Soulaimane Salehddine <up202312271@edu.fc.up.pt>

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
import matplotlib.pyplot as plt
from Bio.Blast import NCBIWWW
from Bio.SeqIO import SeqRecord
from Bio import AlignIO, SeqIO, SearchIO, Phylo
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator, DistanceTreeConstructor,
    ParsimonyTreeConstructor, NNITreeSearcher, ParsimonyScorer
)


UNIPROT_OUT_FN = "sequence.fasta"
BLAST_OUT_FN = "sequences_to_analyse.fasta"
MSA_OUT_FN = "alignment.txt"

ALN_DATA_FORMAT = "nexus"
BLAST_LIMIT = 10 # should be over 10 to ensure that 10 distinct sequences are found


def get_uniprot_id() -> str:
    '''
    Get the protein sequence identifier from the command line arguments. Validate the identifier.

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


def get_protein_seq(prot_seq_id: str, filename: str = UNIPROT_OUT_FN) -> SeqRecord:
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


def _read_seqs_from_fasta(filename: str) -> list[SeqRecord]:
    '''
    Read protein sequences from a FASTA file.

    Args:
        filename (str): Name of the file containing the protein sequence

    Returns:
        list[SeqRecord]: List of protein sequences
    '''
    with open(filename, encoding="utf-8") as f:
        return list(SeqIO.parse(f, "fasta"))


def run_blast_search(
        sequence: Optional[SeqRecord] = None,
        filename: str = BLAST_OUT_FN
    ) -> list[SeqRecord]:
    '''
    Run a BLAST search for the protein sequence and get the 10 best distinct matches
    that are not human or synthetic constructs. Save the results to a file.

    Args:
        sequence (SeqRecord): Protein sequence
        filename (str): Name of the file to save the BLAST search results to

    Returns:
        list[SeqRecord]: List of protein sequences representing the BLAST search results

    Note:
        Some results may be labelled as "Homo sapiens". The reason for inclusion of these sequences
        is that they are found in multiple species and human is just one of them.
    '''
    if not sequence:
        sequence = _read_seqs_from_fasta(UNIPROT_OUT_FN)[0]

    entrez_query = 'NOT "Homo sapiens"[Organism] NOT "synthetic construct"[Organism]'
    with NCBIWWW.qblast("blastp", "nr", sequence.seq, entrez_query=entrez_query,
                        hitlist_size=BLAST_LIMIT) as result_handle:
        blast_results = SearchIO.parse(result_handle, "blast-xml")
        blast_records = []
        for query_result in blast_results:
            for hit in query_result.hits:
                for fragment in hit.fragments:
                    # todo: make sure to keep only one match for each organism
                    #           (need to query for more than 10)
                    blast_records.append(fragment.hit)
        SeqIO.write(blast_records, filename, "fasta")
    return blast_records


def perform_msa(
        sequences: Optional[list[SeqRecord]] = None,
        aln_data_format: str = ALN_DATA_FORMAT,
        filename: str = MSA_OUT_FN
    ) -> str:
    '''
    Perform Multiple Sequence Alignment (MSA) using the Clustal Omega algorithm and save the
    results to a file.

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
        sequences = _read_seqs_from_fasta(UNIPROT_OUT_FN) + _read_seqs_from_fasta(BLAST_OUT_FN)

    url = 'https://www.ebi.ac.uk/Tools/services/rest/clustalo'
    form_data = {
        'email': 'up202311235@edu.fc.up.pt',
        'outfmt': aln_data_format,
        'sequence': '\n'.join([seq.format('fasta') for seq in sequences]),
        'stype': 'protein',
    }
    response = requests.post(f'{url}/run', data=form_data, timeout=100)
    if not response.ok:
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


def _save_phylotree(tree: Phylo.BaseTree.Tree, out_filename: str, out_format: str):
    '''
    Save the phylogenetic tree to a file.

    Args:
        tree (Phylo.BaseTree.Tree): Phylogenetic tree
        out_format (str): Format of the phylogenetic tree
            Accepted values: 'png', 'pdf', 'svg', 'txt'

    Raises:
        ValueError: If the output format is invalid
    '''
    if out_format == 'txt':
        with open(f'{out_filename}.txt', 'w', encoding='utf-8') as f:
            Phylo.draw_ascii(tree, file=f)
    elif out_format in ['png', 'pdf', 'svg']:
        _, ax = plt.subplots(figsize=(12, 8))
        Phylo.draw(tree, do_show=False, axes=ax, branch_labels=lambda c: round(c.branch_length, 4))
        plt.savefig(f'{out_filename}.{out_format}')
    else:
        raise ValueError(f'Invalid output format: {out_format}')


def build_phylotree(
        alignment_data: Optional[str] = None,
        aln_data_format: str = ALN_DATA_FORMAT,
        constructor_method: str = 'upgma',
        out_filename: str = 'phylotree',
        out_format: str | list[str] = None
    ):
    '''
    Build a phylogenetic tree based on the MSA results using the specified constructor method.
    Plot and save the tree in the specified format(s).

    Args:
        alignment_data (str): MSA results
        aln_data_format (str): Format of the MSA results
        constructor_method (str): Method to construct the phylogenetic tree
            Accepted values: 'upgma', 'nj', 'parsimony', Default: 'upgma'
        out_format (str | list[str]): Format(s) of the phylogenetic tree
            Accepted values: 'png', 'pdf', 'svg', 'txt'. Default: all formats.
    '''
    if not out_format:
        out_format = ['png', 'pdf', 'svg', 'txt']
    if not alignment_data:
        with open(MSA_OUT_FN, encoding='utf-8') as f:
            alignment_data = f.read()
    alignment_data = AlignIO.read(StringIO(alignment_data), aln_data_format)

    match constructor_method:
        case 'upgma' | 'nj':
            calculator = DistanceCalculator("identity")
            dm = calculator.get_distance(alignment_data)
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm) if constructor_method == 'upgma' else constructor.nj(dm)
        case 'parsimony':
            constructor = ParsimonyTreeConstructor(NNITreeSearcher(ParsimonyScorer()))
            tree = constructor.build_tree(alignment_data)
        case _:
            raise ValueError(f'Invalid constructor method: {constructor_method}')

    if isinstance(out_format, str):
        _save_phylotree(tree, out_filename, out_format)
    elif isinstance(out_format, list):
        for fmt in out_format:
            _save_phylotree(tree, out_filename, fmt)


def main():
    '''
    Main function.
    '''
    # Step 1
    prot_seq_id = get_uniprot_id()
    get_protein_seq(prot_seq_id)

    # Step 2
    run_blast_search()

    # Step 3
    # perform_msa()

    # Step 4
    # build_phylotree()


if __name__ == "__main__":
    main()
