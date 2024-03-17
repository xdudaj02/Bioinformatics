from sequence_alignments import *

SEQ_1 = 'TATTCG'
SEQ_2 = 'ATTTCC'
MATCH = 2
MISMATCH = -1
GAP = -3

def get_alignments(seq_1 : str = SEQ_1, seq_2 : str = SEQ_2, match : int = MATCH, mismatch : int = MISMATCH, gap : int = GAP):
    '''
    Given two sequences and the match, mismatch, and gap scores, print all the relevant information about the global and local alignments.

    Args:
    - seq_1 (str): the first sequence to be aligned
    - seq_2 (str): the second sequence to be aligned
    - match (int): the score for a match
    - mismatch (int): the score for a mismatch
    - gap (int): the score for a gap
    '''
    alphabet = set(seq_1 + seq_2)

    # Determine the score and traceback matrix. Determine the best score.
    sub_matrix = create_submat(match, mismatch, alphabet)
    global_score_matrix, global_traceback_matrix = needleman_Wunsch(seq_1, seq_2, sub_matrix, gap)
    global_best_score = global_score_matrix[-1][-1]
    local_score_matrix, local_traceback_matrix, local_best_score = smith_Waterman(seq_1, seq_2, sub_matrix, gap)

    print('Global score matrix:', global_score_matrix)
    print('Global traceback matrix:', global_traceback_matrix)
    print('Global alignment score:', global_best_score)
    print('Local score matrix:', local_score_matrix)
    print('Local traceback matrix:', local_traceback_matrix)
    print('Local alignment score:', local_best_score)
    
    # Retrieve the optimal sequence alignment.
    optimal_global_alignment = recover_align(global_traceback_matrix, seq_1, seq_2)
    optimal_local_alignment = recover_align_local(local_score_matrix, local_traceback_matrix, seq_1, seq_2)

    print('Optimal global alignment:', optimal_global_alignment)
    print('Optimal local alignment:', optimal_local_alignment)
    
    # Explain if there are multiple best alignments
    # TODO: ???
    # not sure what to do here as it cant really be done with the provided functions and the solution would probably be too complex to fit on one page


if __name__ == '__main__':
    get_alignments()
