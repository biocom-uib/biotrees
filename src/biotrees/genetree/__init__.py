"""
For us, a Phylogenetic Tree (a `phylotree` instance) is a special case of shape in which labels of leaves can be distinguished.
We import modules shape.py and _newick.py; the latter will be used for reading Newick code in string format and turning it
into trees; it can be found in https://github.com/glottobank/python-newick.
"""

from biotrees.phylotree import PhyloTree, shape_to_phylotree


class GeneTree(PhyloTree):
    """
     Thus, the class `genetree` must extend the class `phylotree`, overriding those methods in `phylotree` that do not
     take into account the labelling morphism.
    """

    def __init__(self, leaf, label, children = None): #cambiar comentarios
        """
        Create a new `genetree` object.
        The boolean is_leaf is True if the object is a leaf; it is False otherwise.
        :param leaf: if is_leaf, its label; otherwise, `None`.
        :param children: if not is_leaf, the `phylotree` objects which are descendants of the considered object;
        otherwise, `None`.
        :return: `phylotree` instance.
        """
        super(GeneTree, self).__init__(leaf, children)
        self.label = label

    def __str__(self):
        from biotrees.genetree.newick import to_newick
        return to_newick(self)


    def compare(self, t2):
        """
        Compare self with another `shape` object. We use lexicographical order in order to compare two `shape` instances.
        Leaves in this case are indistinguishable. It returns anint c, which is 0 if self and T2 are equal, < 0 if self < T2,
        and > 0 if self > T2.
        :param t2: the `shape` object against which we compare self.
        :return: int instance.
        """
        if self.is_leaf() != t2.is_leaf():
            return (not self.is_leaf()) - (not t2.is_leaf())

        if self.is_leaf():
            if self.label < t2.label:
                return -1
            elif self.label > t2.label:
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


    def clone(self):
        """
        Returns new `phylotree` instance which is exactly the same as self.
        :return: `phylotree` instance.
        """
        if self.is_leaf:
            return GeneTree(self.leaf, self.label)
        else:
            return GeneTree(None, None, [ch.clone() for ch in self.children])


    def labels(self):
        """
        Yields the (labels of the) leaves of self.
        :return: `phylotree` instance.
        """
        if self.is_leaf():
            yield self.label
        else:
            for x in self.children:
                for l in x.labels():
                    yield l

    def label_set(self):
        return sorted(set(self.labels()))

    def count_labels(self):
        return len(set(self.labels()))

    def get_labels(self):
        return sorted(self.labels())

    def is_phylo(self):
        """
        Returns True if self is phylogenetic (namely, if it has no repeated leaves). Returns False otherwise.
        :return: bool instance.
        """
        L = self.get_labels()
        for i in range(1, len(L)):
            if L[i] == L[i-1]:
                return False
        return True


def phylotree_to_genetree(phylo):
    """
    Returns a `genetree` instance with the same labels as the leafs in the `phylotree` given as a parameter.
    :param phylo: `phylotree` instance.
    :return: `genetree` instance.
    """
    i = [-1]

    def recurse(shape):
        if shape.is_leaf:
            i[0] += 1
            return GeneTree(i[0], i[0])
        else:
            return GeneTree(None, None, [recurse(child) for child in shape.children])

    return recurse(phylo)


def genetree_to_phylo(genet):
    """
    Returns a `phylotree` instance with the same leaves as the leafs in the `genetree` given as a parameter. Namely,
    it "forgets" the labels of genet.
    :param genet: `genetree` instance.
    :return: `phylotree` instance.
    """
    if genet.is_leaf:
        return PhyloTree(genet.leaf)
    else:
        return PhyloTree(None, [genetree_to_phylo(child) for child in genet.children])


def shape_to_genetree(original_shape):
    phylo = shape_to_phylotree(original_shape)
    return phylotree_to_genetree(phylo)


def genetree_to_shape(genet):
    return genet.shape()


