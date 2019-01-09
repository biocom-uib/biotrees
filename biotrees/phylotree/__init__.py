"""
For us, a Phylogenetic Tree (a `PhyloTree` instance) is a special case of Shape in which names of leaves can be
distinguished.
"""

from biotrees.shape import Shape, is_binary, count_leaves, get_depth


class PhyloTree(Shape):
    """
     Thus, the class `PhyloTree` must extend the class `Shape`, overriding those methods in `Shape` that do not take
     into account the distinguishability of leaves.
    """

    def __init__(self, leaf, children=None):
        """
        Create a new `PhyloTree` object.
        The boolean is_leaf is True if the object is a leaf; it is False otherwise.
        :param leaf: if is_leaf, its name; otherwise, `None`.
        :param children: `list` instance.
        :return: `PhyloTree` instance.
        """
        super(PhyloTree, self).__init__(children)
        self.leaf = leaf
        assert not self.is_leaf() or (leaf is not None)

    def clone(self):
        """
        Returns new `PhyloTree` instance which is exactly the same as self.
        :return: `PhyloTree` instance.
        """
        if self.is_leaf():
            return PhyloTree(self.leaf)
        else:
            return PhyloTree(None, [ch.clone() for ch in self.children])

    def compare_with_shape_lex(self, t2):
        l1, l2 = self.is_leaf(), t2.is_leaf()
        if l1 and l2:
            if self.leaf < t2.leaf:
                return 0, -1
            elif self.leaf == t2.leaf:
                return 0, 0
            else:
                return 0, 1
        elif l1:
            return -1, 0
        elif l2:
            return 1, 0

        cs = len(self.children) - len(t2.children)
        if cs != 0:
            return cs, 0
        cl_fst = 0

        for ch1, ch2 in zip(self.children, t2.children):
            cs, cl = ch1.compare_with_shape_lex(ch2)

            if cs != 0:
                return cs, 0
            if cl_fst == 0:
                cl_fst = cl

        return cs, cl_fst

    def compare(self, t2): # cambiar comentarios
        """
        Compare self with another `Shape` object. We use lexicographical order in order to compare two `Shape` instances.
        Leaves in this case are indistinguishable. It returns anint c, which is 0 if self and T2 are equal, < 0 if
        self < T2, and > 0 if self > T2.
        :param t2: the `Shape` object against which we compare self.
        :return: `int` instance.
        """

        c = self.compare_with_shape_lex(t2)
        # assert c[0] == self.shape().compare(t2.shape()), str(self) + '    ' + str(t2) + '     ' + str(c)

        if c < (0,0):
            return -1
        elif c > (0,0):
            return 1
        else:
            return 0

    def __str__(self):
        from biotrees.phylotree.newick import to_newick
        return to_newick(self)

    def shape(self):
        """
        Returns the `Shape` associated to self. Namely, it "forgets" the names of the leaves.
        :return: `Shape` instance.
        """
        if self.is_leaf():
            return Shape(None)
        else:
            return Shape([x.shape() for x in self.children])


def leaves(t):
    """
    Yields the leaves of t.
    :return: `PhyloTree` instance.
    """
    if t.is_leaf():
        yield t
    else:
        for x in t.children:
            for l in leaves(x):
                yield l


def get_leaves(t):
    """
    Returns a list with the leaves that appear in t, sorted in lexicographical order.
    Repetitions may arise if the user enters trees which are not phylogenetic.
    :return: `list` instance.
    """
    return sorted(leaves(t))


def get_leaves_names(t):
    """
    Returns a list with the leaves' names that appear in t, sorted in lexicographical order.
    Repetitions may arise if the user enters trees which are not phylogenetic.
    :return: `list` instance.
    """
    return [l.leaf for l in get_leaves(t)]


def get_leaves_names_set(t):
    """
    Returns a list with the leaves' names that appear in t, sorted in lexicographical order.
    Repetitions shall not be included.
    :return: `list` instance.
    """
    return sorted(set(get_leaves_names(t)))


def is_phylo(t):
    """
    Returns True if t is phylogenetic (namely, if it has no repeated leaves). Returns False otherwise.
    :return: `bool` instance.
    """
    l = get_leaves(t)
    for i in range(1, len(l)):
        if l[i] == l[i-1]:
            return False
    return True


def shape_to_phylotree(shape):
    """
    Returns a `PhyloTree` instance with the same shape as the one given as a parameter.
    :param shape: `Shape` instance.
    :return: `PhyloTree` instance.
    """
    i = [-1]

    def recurse(sh):
        if sh.is_leaf():
            i[0] += 1
            return PhyloTree(str(i[0]))
        else:
            return PhyloTree(None, [recurse(child) for child in sh.children])

    return recurse(shape)


def phylotree_to_shape(phylo):
    """
    Returns a `Shape` instance with the same shape as the `PhyloTree` given as a parameter.
    :param phylo: `PhyloTree` instance.
    :return: `Shape` instance.
    """
    return phylo.shape()
