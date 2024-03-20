# -*- coding: utf-8 -*-

## utility functions
def create_submat(match, mismatch, alphabet):
    '''Create a substitution matrix for the given alphabet, with the given match and mismatch scores. The matrix is a dictionary with keys of the form 'c1c2' and values the corresponding scores.'''
    sm = {}
    for c1 in alphabet:
        for c2 in alphabet:
            if (c1 == c2):
                sm[c1 + c2] = match
            else:
                sm[c1 + c2] = mismatch
    return sm

def score_pos(c1, c2, sm, g):
    '''Score of a position in the alignment. If c1 or c2 is a gap, return g. Otherwise, return the value in the substitution matrix.'''
    if c1 == "-" or c2 == "-":
        return g
    else:
        return sm[c1 + c2]


## global alignment functions
def needleman_wunsch(seq1, seq2, sm, g):
    '''Calculate the score and traceback matrices for the global alignment of seq1 and seq2.'''
    S = [[0]]
    T = [[0]]
    # initialize gaps in rows
    for j in range(1, len(seq2) + 1):
        S[0].append(g * j)
        T[0].append(3)  # horizontal move: 3
    # initialize gaps in cols
    for i in range(1, len(seq1) + 1):
        S.append([g * i])
        T.append([2])  # vertical move: 2
    # apply the recurrence to fill the matrices
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos(seq1[i], seq2[j], sm, g); # diagonal
            s2 = S[i][j+1] + g  # vertical
            s3 = S[i+1][j] + g # horizontal
            S[i+1].append(max(s1, s2, s3)) # na matrix score add max value
            T[i+1].append(max_3t(s1, s2, s3))
    return (S, T)

def max_3t(v1, v2, v3):
    '''Get direction of traceback. Return 1 if v1 is the maximum, 2 if v2 is the maximum, 3 if v3 is the maximum, 4 if v1 and v2 are equal maxima, 5 if v1 and v3 are equal maxima, 6 if v2 and v3 are equal maxima, 7 if all three are equal maxima.'''
    if v1 > v2 and v1 > v3:
        return 1
    elif v2 > v1 and v2 > v3:
        return 2
    elif v3 > v1 and v3 > v2:
        return 3
    elif v1 == v2 == v3:
        return 7
    elif v1 == v2:
        return 4
    elif v1 == v3:
        return 5
    elif v2 == v3:
        return 6

def recover_align(T, seq1, seq2, nseq1="", nseq2=""):
    '''Recover the optimal alignment from the traceback matrix T and the two sequences seq1 and seq2. The function returns a list of tuples, each containing a pair of aligned sequences.'''
    i = len(seq1)
    j = len(seq2)
    if i == 0 and j == 0: # end of alignment
        return [(nseq1, nseq2)]
    if T[i][j] == 1: # diagonal move
        return recover_align(T, seq1[0:i-1], seq2[0:j-1], seq1[i-1] + nseq1, seq2[j-1] + nseq2)
    elif T[i][j] == 2: # vertical move
        return recover_align(T, seq1[0:i-1], seq2, seq1[i-1] + nseq1, "-" + nseq2)
    elif T[i][j] == 3: # horizontal move
        return recover_align(T, seq1, seq2[0:j-1], "-" + nseq1, seq2[j-1] + nseq2)
    elif T[i][j] == 4: # diagonal or vertical move
        opt_1 = recover_align(T, seq1[0:i-1], seq2[0:j-1], seq1[i-1] + nseq1, seq2[j-1] + nseq2)
        opt_2 = recover_align(T, seq1[0:i-1], seq2, seq1[i-1] + nseq1, "-" + nseq2)
        return opt_1 + opt_2
    elif T[i][j] == 5: # diagonal or horizontal move
        opt_1 = recover_align(T, seq1[0:i-1], seq2[0:j-1], seq1[i-1] + nseq1, seq2[j-1] + nseq2)
        opt_2 = recover_align(T, seq1, seq2[0:j-1], "-" + nseq1, seq2[j-1] + nseq2)
        return opt_1 + opt_2
    elif T[i][j] == 6: # vertical or horizontal move
        opt_1 = recover_align(T, seq1[0:i-1], seq2, seq1[i-1] + nseq1, "-" + nseq2)
        opt_2 = recover_align(T, seq1, seq2[0:j-1], "-" + nseq1, seq2[j-1] + nseq2)
        return opt_1 + opt_2
    elif T[i][j] == 7: # diagonal or vertical or horizontal move
        opt_1 = recover_align(T, seq1[0:i-1], seq2[0:j-1], seq1[i-1] + nseq1, seq2[j-1] + nseq2)
        opt_2 = recover_align(T, seq1[0:i-1], seq2, seq1[i-1] + nseq1, "-" + nseq2)
        opt_3 = recover_align(T, seq1, seq2[0:j-1], "-" + nseq1, seq2[j-1] + nseq2)
        return opt_1 + opt_2 + opt_3
    
def get_global_alignment(seq1, seq2, sm, g):
    '''Get the global alignment of seq1 and seq2 using the Needleman-Wunsch algorithm. Return a list of tuples, each containing a pair of aligned sequences.'''
    _, T = needleman_wunsch(seq1, seq2, sm, g)
    return recover_align(T, seq1, seq2)

## local alignment functions
def smith_waterman(seq1, seq2, sm, g):
    '''Calculate the score and traceback matrices for the local alignment of seq1 and seq2. Return the score matrix, the traceback matrix, and the best score.'''
    S = [[0]]
    T = [[0]]
    maxscore = 0
    # first row filled with zero
    for j in range(1, len(seq2)+1):
        S[0].append(0)
        T[0].append(0)
    # first column filled with zero
    for i in range(1, len(seq1)+1):
        S.append([0])
        T.append([0])
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos(seq1[i], seq2[j], sm, g); 
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            b = max(s1, s2, s3)
            if b <= 0:
                S[i+1].append(0)
                T[i+1].append(0)
            else:
                S[i+1].append(b)
                T[i+1].append(max_3t(s1, s2, s3))
                if b > maxscore: 
                    maxscore = b
    return (S, T, maxscore)

def recover_align_local(S, T, seq1, seq2):
    '''Recover the optimal local alignment from the score and traceback matrices S and T, and the two sequences seq1 and seq2. The function returns a list of tuples, each containing a pair of aligned sequences.'''
    def recover_align_local_inner(S, T, seq1, seq2, nseq1, nseq2, i, j):
        if T[i][j] == 0: # end of alignment
            return [(nseq1, nseq2)]
        if T[i][j] == 1: # diagonal move
            return recover_align_local_inner(S, T, seq1, seq2, seq1[i-1] + nseq1, seq2[j-1] + nseq2, i-1, j-1)
        elif T[i][j] == 2: # vertical move
            return recover_align_local_inner(S, T, seq1, seq2, seq1[i-1] + nseq1, "-" + nseq2, i-1, j)
        elif T[i][j] == 3: # horizontal move
            return recover_align_local_inner(S, T, seq1, seq2, "-" + nseq1, seq2[j-1] + nseq2, i, j-1)
        elif T[i][j] == 4: # diagonal or vertical move
            opt_1 = recover_align_local_inner(S, T, seq1, seq2, seq1[i-1] + nseq1, seq2[j-1] + nseq2, i-1, j-1)
            opt_2 = recover_align_local_inner(S, T, seq1, seq2, seq1[i-1] + nseq1, "-" + nseq2, i-1, j)
            return opt_1 + opt_2
        elif T[i][j] == 5: # diagonal or horizontal move
            opt_1 = recover_align_local_inner(S, T, seq1, seq2, seq1[i-1] + nseq1, seq2[j-1] + nseq2, i-1, j-1)
            opt_2 = recover_align_local_inner(S, T, seq1, seq2, "-" + nseq1, seq2[j-1] + nseq2, i, j-1)
            return opt_1 + opt_2
        elif T[i][j] == 6: # vertical or horizontal move
            opt_1 = recover_align_local_inner(S, T, seq1, seq2, seq1[i-1] + nseq1, "-" + nseq2, i-1, j)
            opt_2 = recover_align_local_inner(S, T, seq1, seq2, "-" + nseq1, seq2[j-1] + nseq2, i, j-1)
            return opt_1 + opt_2
        elif T[i][j] == 7: # diagonal or vertical or horizontal move
            opt_1 = recover_align_local_inner(S, T, seq1, seq2, seq1[i-1] + nseq1, seq2[j-1] + nseq2, i-1, j-1)
            opt_2 = recover_align_local_inner(S, T, seq1, seq2, seq1[i-1] + nseq1, "-" + nseq2, i-1, j)
            opt_3 = recover_align_local_inner(S, T, seq1, seq2, "-" + nseq1, seq2[j-1] + nseq2, i, j-1)
            return opt_1 + opt_2 + opt_3
        
    res = []
    for (i, j) in max_mat(S):
        res += recover_align_local_inner(S, T, seq1, seq2, "", "", i, j)
    return res

def max_mat(mat):
    '''Find the maximum value in a matrix and return the coordinates of all the cells with that value.'''
    maxval = mat[0][0]
    max_coord = [(0, 0)]
    for i in range(0, len(mat)):
        for j in range(0, len(mat[i])):
            if mat[i][j] > maxval:
                maxval = mat[i][j]
                max_coord = [(i, j)]
            elif mat[i][j] == maxval:
                max_coord.append((i, j))
    return max_coord

def get_local_alignment(seq1, seq2, sm, g):
    '''Get the local alignment of seq1 and seq2 using the Smith-Waterman algorithm. Return a list of tuples, each containing a pair of aligned sequences.'''
    S, T, _ = smith_waterman(seq1, seq2, sm, g)
    return recover_align_local(S, T, seq1, seq2)
