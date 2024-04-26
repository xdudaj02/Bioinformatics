from MyAlign import MyAlign
from MySeq import MySeq
from SubstMatrix import SubstMatrix

class PairwiseAlignment:
    def __init__(self, sm, g):
        self.g = g
        self.sm = sm
        self.s = None
        self.t = None
        self.seq1 = None
        self.seq2 = None

    def score_pos(self, c1, c2):
        if c1 == "-" and c2 == "-":
            return 0
        if c1 == "-" or c2 == "-":
            return self.g
        return self.sm[c1, c2]

    # def score_alin(self, alin):
    #     res = 0
    #     for i in range(len(alin)):
    #         res += self.score_pos(alin[0][i], alin[1][i])
    #     return res

    def needleman_Wunsch(self, seq1, seq2):
        if seq1.seq_type != seq2.seq_type:
            return None
        self.s = [[0]]
        self.t = [[0]]
        self.seq1 = seq1
        self.seq2 = seq2
        for j in range(1, len(seq2) + 1):
            self.s[0].append(self.g * j)
            self.t[0].append(3)
        for i in range(1, len(seq1) + 1):
            self.s.append([self.g * i])
            self.t.append([2])
        for i, _ in enumerate(seq1):
            for j, _ in enumerate(seq2):
                s1 = self.s[i][j] + self.score_pos(seq1[i], seq2[j])
                s2 = self.s[i][j + 1] + self.g
                s3 = self.s[i + 1][j] + self.g
                self.s[i + 1].append(max(s1, s2, s3))
                self.t[i + 1].append(max3t(s1, s2, s3))
        return self.s[len(seq1)][len(seq2)]

    def recover_align(self):
        res = ["", ""]
        i = len(self.seq1)
        j = len(self.seq2)
        while i > 0 or j > 0:
            if self.t[i][j]==1:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = self.seq2[j - 1] + res[1]
                i -= 1
                j -= 1
            elif self.t[i][j] == 3:
                res[0] = "-" + res[0]
                res[1] = self.seq2[j - 1] + res[1]
                j -= 1
            else:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = "-" + res[1]
                i -= 1
        return MyAlign(res, self.seq1.seq_type)

    def smith_Waterman(self, seq1, seq2):
        if seq1.seq_type != seq2.seq_type:
            return None
        self.s = [[0]]
        self.t = [[0]]
        self.seq1 = seq1
        self.seq2 = seq2
        maxscore = 0
        for j in range(1, len(seq2)+1):
            self.s[0].append(0)
            self.t[0].append(0)
        for i in range(1, len(seq1)+1):
            self.s.append([0])
            self.t.append([0])
        for i in enumerate(seq1):
            for j in enumerate(seq2):
                s1 = self.s[i][j] + self.score_pos(seq1[i], seq2[j])
                s2 = self.s[i][j + 1] + self.g
                s3 = self.s[i + 1][j] + self.g
                b = max(s1, s2, s3)
                if b <= 0:
                    self.s[i + 1].append(0)
                    self.t[i + 1].append(0)
                else:
                    self.s[i + 1].append(b)
                    self.t[i + 1].append(max3t(s1, s2, s3))
                    maxscore = b if b > maxscore else maxscore
        return maxscore

    def recover_align_local(self):
        res = ["", ""]
        maxscore = 0
        maxrow = 0
        maxcol = 0
        for i in range(1,len(self.s)):
            for j in range(1, len(self.s[i])):
                if self.s[i][j] > maxscore:
                    maxscore = self.s[i][j]
                    maxrow = i
                    maxcol = j
        i = maxrow
        j = maxcol
        while i > 0 or j > 0:
            if self.t[i][j] == 1:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = self.seq2[j - 1] + res[1]
                i -= 1
                j -= 1
            elif self.t[i][j] == 3:
                res[0] = "-" + res[0]
                res[1] = self.seq2[j - 1] + res[1]
                j -= 1
            elif self.t[i][j] == 2:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = "-" + res[1]
                i -= 1
            else:
                break
        return MyAlign(res, self.seq1.seq_type)


def max3t(v1, v2, v3):
    if v1 > v2:
        if v1 > v3:
            return 1
        return 3
    if v2 > v3:
        return 2
    return 3

def print_mat(mat):
    for i in mat:
        print(i)


#### TESTS #####

def test():
    seq1 = MySeq("ATGATATGATGATT")
    seq2 = MySeq("GATGAATAGATGTGT")
    sm = SubstMatrix()
    sm.create_submat(3, -1, "ACGT")
    alin = PairwiseAlignment(sm, -3)
    print(alin.smith_Waterman(seq1, seq2))
    print_mat(alin.s)
    print(alin.recover_align_local())

    print(alin.needleman_Wunsch(seq1,seq2))
    print_mat(alin.s)
    print(alin.recover_align())


if __name__ == "__main__":
    test()
