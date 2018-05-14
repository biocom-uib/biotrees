"""
For us, a Phylogenetic Tree (a `PhyloTree` instance) is a special case of Shape in which names of leaves can be
distinguished.
We import modules Shape.py and _newick.py; the latter will be used for reading Newick code in string format and turning
it into trees; it can be found in https://github.com/glottobank/python-newick.
"""

from biotrees.shape import Shape


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
        assert not self.is_leaf or (leaf is not None)

    def clone(self):
        """
        Returns new `PhyloTree` instance which is exactly the same as self.
        :return: `PhyloTree` instance.
        """
        if self.is_leaf:
            return PhyloTree(self.leaf)
        else:
            return PhyloTree(None, [ch.clone() for ch in self.children])

    def compare(self, t2): # cambiar comentarios
        """
        Compare self with another `Shape` object. We use lexicographical order in order to compare two `Shape` instances.
        Leaves in this case are indistinguishable. It returns anint c, which is 0 if self and T2 are equal, < 0 if
        self < T2, and > 0 if self > T2.
        :param t2: the `Shape` object against which we compare self.
        :return: `int` instance.
        """
        if self.is_leaf != t2.is_leaf:
            return (not self.is_leaf) - (not t2.is_leaf)

        if self.is_leaf:
            if self.leaf < t2.leaf:
                return -1
            elif self.leaf > t2.leaf:
                return 1
            else:
                return 0

        if len(self.children) != len(t2.children):
            return len(self.children) - len(t2.children)

        for t1, t2 in zip(self.children, t2.children):
            c = t1.compare(t2)
            if c != 0:
                return c
        return 0

    def shape(self):
        """
        Returns the `Shape` associated to self. Namely, it "forgets" the names of the leaves.
        :return: `Shape` instance.
        """
        if self.is_leaf:
            return Shape(None)
        else:
            return Shape([x.shape() for x in self.children])

    def leaves(self):
        """
        Yields the (names of the) leaves of self.
        :return: `PhyloTree` instance.
        """
        if self.is_leaf:
            yield self
        else:
            for x in self.children:
                for l in x.leaves():
                    yield l

    def get_leaves(self):
        """
        Returns a list with the leaves that appear in self, sorted in lexicographical order.
        Repetitions may arise if the user enters trees which are not phylogenetic.
        :return: `list` instance.
        """
        return sorted(list(self.leaves()))

    def get_leaves_names(self):
        """
        Returns a list with the leaves' names that appear in self, sorted in lexicographical order.
        Repetitions may arise if the user enters trees which are not phylogenetic.
        :return: `list` instance.
        """
        return [l.leaf for l in self.get_leaves()]

    def leaves_names_set(self):
        """
        Returns a list with the leaves' names that appear in self, sorted in lexicographical order.
        Repetitions shall not be included.
        :return: `list` instance.
        """
        return sorted(list(set(self.get_leaves_names())))

    def is_phylo(self):
        """
        Returns True if self is phylogenetic (namely, if it has no repeated leaves). Returns False otherwise.
        :return: `bool` instance.
        """
        l = self.get_leaves()
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
        if sh.is_leaf:
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
