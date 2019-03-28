from biotrees.util import binom2

from biotrees.shape.generator import comb, binary_max_balanced


def cophenetic_index(tree):
    def go(t):
        if t.is_leaf():
            return 0, 1

        cophs, kappas = zip(*map(go, t.children))
        kappa = sum(kappas)
        coph = binom2(kappa) + sum(cophs)
        return coph, kappa

    if tree.is_leaf():
        return 0
    else:
        return sum(go(ch)[0] for ch in tree.children)


def max_cophenetic(n):
    """
    Returns all the `Shape` instances that attain the maximum Cophenetic index with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    return [comb(n)]


def min_cophenetic(n):
    """
    Returns all the `Shape` instances that attain the minimum Cophenetic index with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    return [binary_max_balanced(n)]
