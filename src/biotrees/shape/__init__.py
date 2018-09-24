from biotrees.util import iter_merge, skip_nth

"""
A `Shape` represents a topological tree. The data structure implemented here is of recursive type: a `Shape` can be
either a leaf or a list of `Shape` objects. Leaves are not distinguishable, but we know that they are leaves.
We choose a sorted shape to be the class representant of all shapes isomorphic to it.
"""

__all__ = ['Shape']


class Shape(object):
    LEAF = None # filled in after class def

    """
    A `Shape` instance is either a leaf or a list of `Shape` instances that hang from a root.
    """
    def __init__(self, children=None):
        """
        Create a new `Shape` object.
        The boolean is_leaf is True if the object is a leaf; it is False otherwise.
        :param children: `list` instance.
        :return: `Shape` instance.
        """
        assert children is None or len(children) > 0
        self.children = children

    def is_leaf(self):
        return self.children is None

    def clone(self):
        """
        Returns `Shape` instance which is exactly the same as self.
        :return: `Shape` instance.
        """
        if self.is_leaf():
            return Shape.LEAF
        else:
            return Shape([ch.clone() for ch in self.children])

    def _is_sorted(self):
        if self.is_leaf():
            return True

        children = self.children

        return all(ch._is_sorted() for ch in children) and \
               all(ch1 <= ch2 for ch1, ch2 in zip(children[:-1], children[1:]))

    def _sort(self):
        """
        Sorts self using lexicographical order.
        """
        if self.is_leaf():
            return

        for t in self.children:
            t._sort()

        self.children.sort()

    def compare(self, t2):
        """
        Compare self with another `Shape` object. We use lexicographical order in order to compare two `Shape` instances.
        Leaves in this case are indistinguishable. It returns anint c, which is 0 if self and T2 are equal, < 0 if
        self < T2, and > 0 if self > T2.
        :param t2: `Shape` instance.
        :return: `int` instance.
        """
        if self.is_leaf() and t2.is_leaf():
            return 0
        elif self.is_leaf():
            return -1
        elif t2.is_leaf():
            return 1
        else:
            c = len(self.children) - len(t2.children)
            if c != 0:
                return c

            for t1, t2 in zip(self.children, t2.children):
                c = t1.compare(t2)
                if c != 0:
                    return c

            return c

    def __lt__(self, t2):
        """
        Uses the comparing method above to decide if self is less than T2.
        :param t2: the `Shape` object against which we compare self.
        :return: `bool` instance.
        """
        return self.compare(t2) < 0

    def __le__(self, t2):
        """
        Uses the comparing method above to decide if self is less or equal than T2.
        :param t2: the `Shape` object against which we compare self.
        :return: `bool` instance.
        """
        return self.compare(t2) <= 0

    def __eq__(self, t2):
        """
        Uses the comparing method above to decide if self is equal to T2.
        :param t2: the `Shape` object against which we compare self.
        :return: `bool` instance.
        """
        return self.compare(t2) == 0

    def __ne__(self, t2):
        """
        Uses the comparing method above to decide if self is not equal to T2.
        :param t2: the `Shape` object against which we compare self.
        :return: `bool` instance.
        """
        return self.compare(t2) != 0

    def __ge__(self, t2):
        """
        Uses the comparing method above to decide if self is greater or equal than T2.
        :param t2: the `Shape` object against which we compare self.
        :return: `bool` instance.
        """
        return self.compare(t2) >= 0

    def __gt__(self, t2):
        """
        Uses the comparing method above to decide if self is greater than T2.
        :param t2: the `Shape` object against which we compare self.
        :return: `bool` instance.
        """
        return self.compare(t2) > 0

    def __str__(self):
        from .newick import to_newick
        return to_newick(self)

    def __repr__(self):
        return str(self)

    def shape(self):
        """
        Returns the `Shape` associated to self. Namely, it "forgets" the labels of the leafs.
        :return: `Shape` instance.
        """
        return self


Shape.LEAF = Shape()


def is_binary(t):
    """
    Returns True if t is a binary shape.
    :return: `bool` instance
    """
    return t.is_leaf() or \
        (len(t.children) == 2 and all(is_binary(ch) for t in t.children))


def count_leaves(t):
    """
    Returns the number of leaves in t.
    :return: `int` instance.
    """
    if t.is_leaf():
        return 1
    else:
        return sum(count_leaves(t) for t in t.children)


def get_depth(t):
    """
    Returns an integer representing the maximal depth of the shape, from the root to one
    of its furthest leaves.
    :return: `int` instance.
    """
    if t.is_leaf():
        return 0
    else:
        return max(get_depth(ch) for ch in t.children) + 1

def leaf_depths(t):
    """
    Returns a generator of integers representing the depth of each leaf in the tree
    :return: generator of integers
    """
    if t.is_leaf():
        yield 0
    else:
        for ch in t.children:
            for depth in leaf_depths(ch):
                yield depth+1

def get_leaf_depths(t):
    """
    Returns a list of integers representing the depth of each leaf in the tree
    :return: list of integers
    """
    return list(leaf_depths(t))

def count_nodes_by_depth(t):
    total_depth = get_depth(t)
    nodes_by_depth = [0]*(total_depth+1)

    def navigate(t2, d):
        if not t2.is_leaf():
            d1 = d+1
            nodes_by_depth[d1] += len(t2.children)

            for ch in t2.children:
                navigate(ch, d1)

    nodes_by_depth[0] += 1
    navigate(t, 0)
    return nodes_by_depth
