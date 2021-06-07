"""
This computes the extremal values of the Quadratic Colless index for any int n following the algorithm by Krzysztof
Bartoszek, Tomás M. Coronado, Arnau Mir and Francesc Rosselló in their paper "Squaring within the Colless index
yields a better balance index", as well as the value of the index in a given tree.
"""

from biotrees.shape.generator import comb, binary_max_balanced


def binary_qcolless_index(tree):
    """
    Returns the `int` value of the QColless index for a given `Shape` instance, tree.
    :param tree: `Shape` instance.
    :return: `int` instance.
    """
    def go(t):
        if t.is_leaf():
            return 0, 1

        left, right = t.children
        (cil, nl), (cir, nr) = go(left), go(right)
        return (nl - nr)**2 + cil + cir, nl + nr

    return go(tree)[0]


def min_qcolless(n):
    """
    Returns all the `Shape` instances that attain the minimum QColless index with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    return [binary_max_balanced(n)]


def max_qcolless(n):
    """
    Returns all the `Shape` instances that attain the maximum QColless index with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    return [comb(n)]
