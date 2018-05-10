"""
A `Shape` represents a topological tree. The data structure implemented here is of recursive type: a `Shape` can be
either a leaf or a list of `Shape` objects. Leaves are not distinguishable, but we know that they are leaves.
We choose a sorted shape to be the class representant of all shapes isomorphic to it.
"""

__all__ = ['Shape', 'sorted_tree']


class Shape(object):
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
        self.is_leaf  = children is None
        self.children = children
        assert self.is_leaf or len(children) > 0

    def clone(self):
        """
        Returns `Shape` instance which is exactly the same as self.
        :return: `Shape` instance.
        """
        if self.is_leaf:
            return Shape()
        else:
            return Shape([ch.clone() for ch in self.children])

    def sort(self):
        """
        Sorts self using lexicographical order.
        """
        if self.is_leaf:
            return

        for t in self.children:
            t.sort()

        self.children.sort()

    def compare(self, t2):
        """
        Compare self with another `Shape` object. We use lexicographical order in order to compare two `Shape` instances. 
        Leaves in this case are indistinguishable. It returns anint c, which is 0 if self and T2 are equal, < 0 if
        self < T2, and > 0 if self > T2.
        :param t2: `Shape` instance.
        :return: `int` instance.
        """
        if self.is_leaf and t2.is_leaf:
            return 0
        elif self.is_leaf:
            return -1
        elif t2.is_leaf:
            return 1
        else:
            c = len(self.children) - len(t2.children)
            if c != 0:
                return c

            for i in range(0, len(self.children)):
                c = self.children[i].compare(t2.children[i])
                if c != 0:
                    return c
            return 0

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

    def iso(self, t2):
        """
        Since our `Shape` objects are sorted, to know whether two trees are isomorphic or not it suffices to know if
        they are equal or not.
        :param t2: the `Shape` object against which we compare self.
        :return: `bool` instance.
        """
        self.sort()
        t2.sort()
        return self == t2

    def is_symmetric(self):
        """
        Returns True if the root of self is a symmetric node, and False otherwise. If self is a leaf, it returns True:
        ex falso quodlibet.
        :return: `bool` instance.
        """
        if not self.is_leaf:
            return all(self.children[0].iso(x) for x in self.children)
        else:
            return True

    def is_binary(self):
        """
        Returns True if self is a binary shape.
        :return: `bool` instance
        """
        if self.is_leaf:
            return True
        else:
            return len(self.children) == 2 and all(t.is_binary() for t in self.children)

    def count_symmetries(self):
        """
        Returns the number of symmetric interior nodes in self.
        :return: `int` instance.
        """
        if self.is_leaf:
            return 0
        elif all(self.children[0].iso(x) for x in self.children):
            return 1 + sum(x.count_symmetries() for x in self.children)
        else:
            return sum(x.count_symmetries() for x in self.children)

    def count_leaves(self):
        """
        Returns the number of leaves in self.
        :return: `int` instance.
        """
        if self.is_leaf:
            return 1
        else:
            return sum(x.count_leaves() for x in self.children)

    def shape(self):
        """
        Returns the `Shape` associated to self. Namely, it "forgets" the labels of the leafs.
        :return: `Shape` instance.
        """
        return self

    def get_depth(self):
        """
        Returns an integer representing the maximal depth of the shape, from the root to one
        of its furthest leaves.
        :return: `int` instance.
        """
        if self.is_leaf:
            return 0
        else:
            return max(x.get_depth() for x in self.children) + 1


def delete_nodes_with_out_degree_one(t):
    """
    Deletes all nodes in the input `Shape` instance with only one child, for they are redundant
    :return: `Shape` instance
    """
    if t.is_leaf:
        return t
    else:
        if len(t.children) == 1:
            return t.children[0].delete_nodes_with_out_degree_one()
        else:
            return Shape([t.delete_nodes_with_out_degree_one() for t in t.children])


def sorted_by_shape(l):
    """
    Takes a list of lists or tuples whose first component is a Shape instance. It sorts the 
    list according only to its Shape-induced order in the first component of each element.
    :param: `list` instance.
    :return: `list` instance.
    """
    l1 = sorted([(l[i][0], i) for i in range(0, len(l))])
    return [l[tup[1]] for tup in l1]


def sorted_tree(t1):
    """
    Returns a sorted clone of t1.
    :param t1: `Shape` instance
    :return: `Shape` instance
    """
    t2 = t1.clone()
    t2.sort()
    return t2
