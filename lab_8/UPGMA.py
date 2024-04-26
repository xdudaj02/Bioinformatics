from NumMatrix import NumMatrix
from HierarchicalClustering import HierarchicalClustering
from MySeq import MySeq
from PairwiseAlignment import PairwiseAlignment
from SubstMatrix import SubstMatrix


class UPGMA:
    def __init__(self, seqs, alseq):
        self.seqs = seqs        # list of seqs of class MySeq
        self.alseq = alseq      # parameters for pairwise align.
        self.create_mat_dist()  # run this method on the construtor distance matrix

    def create_mat_dist(self):
        # create distance matrix with dim N x N sequences
        seq_count = len(self.seqs)
        self.matdist = NumMatrix(seq_count, seq_count)
        for i, s1 in enumerate(self.seqs):
            for j, s2 in enumerate(self.seqs):
                # align the sequences
                self.alseq.needleman_Wunsch(s1, s2)
                # recover the alignment
                alin = self.alseq.recover_align()
                ncd = 0
                # fill the matrix
                # distance defined as: number of different symbols in the alignment
                for k in range(len(alin)):
                    col = alin.column(k)
                    if col[0] != col[1]:
                        ncd += 1
                # set distance value in the matrix;
                # invoke here the method to update matdist ncd
                self.matdist.set_value(i, j, ncd)

    def run(self):
        # create an object of the class HierarchicalClustering
        ch = HierarchicalClustering(self.matdist)
        # execute the clustering algorithm
        t = ch.execute_clustering()
        return t

def test():
    seq1 = MySeq("ATAGC")
    seq2 = MySeq("ATGAC")
    seq3 = MySeq("AACG")
    seq4 = MySeq("AATCG")
    sm = SubstMatrix()
    sm.create_submat(1, -1, "ACGT")
    alseq = PairwiseAlignment(sm, -2)
    up  = UPGMA([seq1, seq2, seq3, seq4], alseq)
    arv = up.run()
    arv.print_tree()

if __name__ == '__main__':
    test()
