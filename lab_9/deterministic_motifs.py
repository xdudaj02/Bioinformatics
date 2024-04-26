'''Function to expand a motif with ambiguous characters into all possible motifs'''

from Bio.Data import IUPACData

def expand_motif(motif):
    '''Expand a motif with ambiguous characters into all possible motifs'''
    # Initialize list to store expanded motifs
    expanded = [""]

    # Iterate through each character in the motif
    for char in motif:
        # If the character is an IUPAC ambiguous character
        if char in IUPACData.ambiguous_dna_values:
            # Expand the current list of motifs with all possible nucleotides
            replacements = IUPACData.ambiguous_dna_values[char]
            expanded = [expanded + base for expanded in expanded for base in replacements]
        else:
            # If the character is not ambiguous, simply append it to each expanded motif
            expanded = [expanded + char for expanded in expanded]

    return expanded


def test_expand_motif():
    '''Test expanding motifs with ambiguous characters'''
    motif1 = "CCTKCCY"
    motif2 = "CCMCRCCC"

    # Expand motifs
    expanded_motif1 = expand_motif(motif1)
    expanded_motif2 = expand_motif(motif2)

    print(expanded_motif1)
    print(expanded_motif2)

    # Check that the motifs are expanded correctly
    assert expanded_motif1 == ['CCTGCCC', 'CCTGCCT', 'CCTTCCC', 'CCTTCCT']
    assert expanded_motif2 == ['CCACACCC', 'CCACGCCC', 'CCCCACCC', 'CCCCGCCC']


if __name__ == "__main__":
    test_expand_motif()
