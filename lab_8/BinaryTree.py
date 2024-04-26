class BinaryTree:

    def __init__(self, val, dist = 0, left = None, right = None):
        self.value = val
        self.distance = dist
        self.left = left
        self.right = right

    def get_cluster(self):
        # returns a list of all leaves for that tree
        if self.value >= 0:
            return [self.value]
        return self.left.get_cluster() + self.right.get_cluster()

    def print_tree(self):
        self.print_tree_rec(0, "Root")

    def print_tree_rec(self, level, side):
        tabs = "\t" * level
        if self.value >= 0:
            print(tabs, side, " - value:", self.value)
        else:
            print(tabs, side, "- Dist.: ", self.distance)
            if self.left:
                self.left.print_tree_rec(level + 1, "Left")
            if self.right:
                self.right.print_tree_rec(level + 1, "Right")

    def size(self):
        '''
        size of the tree: returns two values
        - number of internal nodes of the tree
        - number of leaves
        '''
        numleaves = 0
        numnodes = 0
        if self.value >= 0:
            numleaves = 1
        else:
            resl = self.left.size() if self.left else (0, 0)
            resr = self.right.size() if self.right else (0, 0)
            numnodes += (resl[0] + resr[0] + 1)
            numleaves += (resl[1] + resr[1])
        return numnodes, numleaves

    def exists_leaf(self, leafnum):
        '''
        returns true or false if leafnum appears in the leaves of the tree
        '''
        return leafnum in self.get_cluster()

    def common_ancestor(self, leaf1, leaf2):
        ''' Return simplest tree that contains leaf1, leaf2
        '''
        if not self.exists_leaf(leaf1) or not self.exists_leaf(leaf2):
            return None
        if self.left and self.left.exists_leaf(leaf1) and self.left.exists_leaf(leaf2):
            return self.left.common_ancestor(leaf1, leaf2)
        if self.right and self.right.exists_leaf(leaf1) and self.right.exists_leaf(leaf2):
            return self.right.common_ancestor(leaf1, leaf2)
        return self

    def distance_leaves(self, leafnum1, leafnum2):
        ''' distance between leafnum1 and leafnum2 using the common ancestor function.
        d(leaf1, leaf2) = 2 * height(common_ancest(leaf1, leaf2))
        '''
        common_ancestor = self.common_ancestor(leafnum1, leafnum2)
        if not common_ancestor:
            return None
        dist_1 = 0
        current = common_ancestor
        while current.value != leafnum1:
            dist_1 += 1
            is_in_left =  current.left and current.left.exists_leaf(leafnum1)
            current = current.left if is_in_left else current.right
        dist_2 = 0
        current = common_ancestor
        while current.value != leafnum2:
            dist_2 += 1
            is_in_left = current.left and current.left.exists_leaf(leafnum2)
            current = current.left if is_in_left else current.right
        return dist_1 + dist_2


def test():
    # leaf
    a = BinaryTree(1)
    b = BinaryTree(2)
    c = BinaryTree(3)
    d = BinaryTree(4)
    # internal nodes
    e = BinaryTree(-1, 2.0, b, c)
    f = BinaryTree(-1, 1.5, d, a)
    g = BinaryTree(-1, 4.5, e, f)
    g.print_tree()
    print(g.get_cluster())
    # testing exercise 3
    print(g.size())
    print(g.exists_leaf(1))
    print(g.exists_leaf(5))
    g.common_ancestor(1, 4).print_tree()
    print(g.distance_leaves(1, 4))
    print(g.distance_leaves(1, 2))

def test_aula():
    # leaves
    a = BinaryTree(2)
    b = BinaryTree(4)
    c = BinaryTree(3)
    d = BinaryTree(1)
    # internal nodes
    e = BinaryTree(-1, 2.0, a, b)
    f = BinaryTree(-1, 3.0, e, c)
    g = BinaryTree(-1, 4.0, f, d)
    g.print_tree()

if __name__ == '__main__':
    test()
