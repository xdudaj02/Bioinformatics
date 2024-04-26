from PairwiseAlignment import PairwiseAlignment
from MyAlign import MyAlign
from MySeq import MySeq
from SubstMatrix import SubstMatrix

class MultipleAlignment():
    def __init__(self, seqs, alignseq):
        self.seqs = seqs # list of MySeq objects
        self.alignpars = alignseq # PairwiseAlignment objects

    def num_seqs(self):
        return len(self.seqs)

    def add_seq_alignment(self, alignment, seq):
        res = []
        for i in range(len(alignment.listseqs) + 1):
            res.append("")
        cons = MySeq(alignment.consensus(), alignment.al_type)
        self.alignpars.needleman_Wunsch(cons, seq)
        align_2 = self.alignpars.recover_align()
        orig = 0
        for i in range(len(align_2)):
            if align_2[0, i]== '-':
                for k in range(len(alignment.listseqs)):
                    res[k] += "-"
            else:
                for k in range(len(alignment.listseqs)):
                    res[k] += alignment[k,orig]
                orig += 1
        res[len(alignment.listseqs)] = align_2.listseqs[1]
        return MyAlign(res, alignment.al_type)

    def align_consensus(self):
        self.alignpars.needleman_Wunsch(self.seqs[0], self.seqs[1])
        res = self.alignpars.recover_align()

        for i in range(2, len(self.seqs)):
            res = self.add_seq_alignment(res, self.seqs[i])
        return res

    def score_column(self, chars_col):
        score = 0
        for i, char1 in enumerate(chars_col):
            for char2 in chars_col[i + 1:]:
                part_score = self.alignpars.score_pos(char1, char2)
                score += part_score if part_score else 0
        return score

    def score_sp(self, alignment):
        score = 0
        for i in range(len(alignment)):
            score += self.score_column(alignment.column(i))
        return score


def print_matrix(mat):
    for x in mat:
        print(x)

def test_prot():
    s1 = MySeq("PHWAS","protein")
    s2 = MySeq("HWASW","protein")
    s3 = MySeq("HPHWA","protein")
    sm = SubstMatrix()
    sm.read_submat_file("data/blosum62.mat", "\t")
    aseq = PairwiseAlignment(sm, -8)
    ma = MultipleAlignment([s1, s2, s3], aseq)
    al = ma.align_consensus()
    print(al)
    print(ma.score_column(al.column(0)))
    print(ma.score_column(al.column(1)))
    print(ma.score_column(al.column(2)))
    print(ma.score_column(al.column(3)))
    print(ma.score_column(al.column(4)))
    print(ma.score_column(al.column(5)))
    print(ma.score_column(al.column(6)))
    print()
    print(ma.score_sp(al))

def test():
    s1 = MySeq("ATAGC")
    s2 = MySeq("AACC")
    s3 = MySeq("ATGAC")

    sm = SubstMatrix()
    sm.create_submat(1,-1,"ACGT")
    aseq = PairwiseAlignment(sm,-1)
    ma = MultipleAlignment([s1,s2,s3], aseq)
    al = ma.align_consensus()
    print(al)

def exercise1():
    s1 = MySeq("ACATATCAT")
    s2 = MySeq("AACAGATCT")
    s3 = MySeq("AGATATTAG")
    s4 = MySeq("GCATCGATT")

    sm = SubstMatrix()
    sm.create_submat(1,-1,"ACGT")
    aseq = PairwiseAlignment(sm,-1)
    ma = MultipleAlignment([s1, s2, s3, s4], aseq)
    al = ma.align_consensus()
    print(al)
    print(ma.score_sp(al))


if __name__ == "__main__":
    test_prot()
    print()
    test()
    print()
    exercise1()
