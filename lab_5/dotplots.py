#############################
### DOTPLOT FUNCTIONS
#############################

# create matrix given dimensions: number of rows and columns filled with zeros
def create_mat(nrows, ncols):
    mat = []
    for i in range(nrows):
        mat.append([])
        for _ in range(ncols):
            mat[i].append(0)
    return mat

# basic dotplot algorithm: fills with ones coincident characters
def dotplot(seq1, seq2):
    ''' Create a matrix based on the the input sequences and fill the cells that correspond to a match
    '''
    mat = create_mat(len(seq1), len(seq2))
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]: mat[i][j] = 1
    return mat

# extended dotplot with window and stringency parameters
def extended_dotplot(seq1, seq2, window, stringency):
    mat = create_mat(len(seq1), len(seq2))
    start = int(window / 2)
    for i in range(start, len(seq1) - start):
        for j in range(start, len(seq2) - start):
            matches = 0
            l = j - start
            for k in range(i - start, i + start + 1):
                if seq1[k] == seq2[l]: matches += 1
                l += 1
                if matches >= stringency: mat[i][j] = 1
    return mat

# prints dotplot
def print_dotplot(mat, s1, s2):
    import sys
    sys.stdout.write("  " + " ".join(list(s2)) + "\n")
    for i in range(len(s1)):
        sys.stdout.write(s1[i])
        for j in range(len(s2)):
            sys.stdout.write(" *" if mat[i][j] == 1 else "  ")
        sys.stdout.write("\n")


def test_diagonal_length(mat, i, j):
    length = 0
    longest = 0
    while i < len(mat) and j < len(mat[0]):
        length += 1 if mat[i][j] == 1 else 0
        longest = length if length > longest else longest
        i += 1
        j += 1
    return length

def test():
    s1 = "CGATATAGATT"
    s2 = "TATATAGTAT"
    mat1 = dotplot(s1, s2)
    print_dotplot(mat1, s1, s2)
    print()
    print(test_diagonal_length(mat1, 2, 3 ))
    print()
    mat2 = extended_dotplot(s1, s2, 5, 4)
    print_dotplot(mat2, s1, s2)


if __name__ == "__main__":
    test()
