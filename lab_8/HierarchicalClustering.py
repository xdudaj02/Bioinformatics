from BinaryTree import BinaryTree
from NumMatrix import NumMatrix

class HierarchicalClustering:
    def __init__(self, matdists):
        self.matdists = matdists # sequence distance matrix

    def execute_clustering(self):
        '''
        returns a tree based on distance matrix
        '''
        # initialize the trees
        trees = []
        # traverse all seqs and create a tree for each sequence
        for i in range(self.matdists.num_rows()):
            # create a tree for each leaf
            BinaryTree(i)
            # add to list of trees
            trees.append(BinaryTree(i))
        # make a copy of the distance matrix to change it
        table_dist = self.matdists.copy()
        # iterations
        for k in range(self.matdists.num_rows(), 1, -1):
            # indices in the matrix for the minimum distance
            i, j = table_dist.min_dist_indexes()
            # create a new tree joining the clusters
            # this will be internal node; height will half of distance in the distance matrix
            # set left tree; set right tree
            n = BinaryTree(-1, table_dist.get_value(i, j) / 2, trees[i], trees[j])
            if k > 2:
                # remove trees being joined from the list
                ti = trees.pop(i)
                tj = trees.pop(j)
                si = len(ti.get_cluster()) #|si|
                sj = len(tj.get_cluster()) #|sj|
                dists = []
                # calculate the distance for the new cluster
                for x in range(table_dist.num_rows()):
                    if x not in [i, j]:
                        # use the weighted average to calculate the distances between the clusters
                        d = (si * table_dist.get_value(i, x) +
                             sj * table_dist.get_value(j, x)) / (si + sj)
                        dists.append(d)
                # update the distance matrix:
                # remove col corresponding to i and j
                table_dist.remove_col(i)
                table_dist.remove_col(j)
                # remove row corresponding to i and j
                table_dist.remove_row(i)
                table_dist.remove_row(j)
                # add row with new distances: dists
                table_dist.add_row(dists)
                # add col with zero distances: of len (|dists| + 1)
                table_dist.add_col([0] * (len(dists) + 1))
                trees.append(n)
            else:
                return n


def test():
    m = NumMatrix(5,5)
    m.set_value(0, 1, 2)
    m.set_value(0, 2, 5)
    m.set_value(0, 3, 7)
    m.set_value(0, 4, 9)
    m.set_value(1, 2, 4)
    m.set_value(1, 3, 6)
    m.set_value(1, 4, 7)
    m.set_value(2, 3, 4)
    m.set_value(2, 4, 6)
    m.set_value(3, 4, 3)
    hc = HierarchicalClustering(m)
    arv = hc.execute_clustering()
    arv.print_tree()

if __name__ == '__main__':
    test()
