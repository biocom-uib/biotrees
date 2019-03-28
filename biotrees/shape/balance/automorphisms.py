from math import factorial

from biotrees.shape.iso import isomorphic
from biotrees.shape.generator import star, comb


def is_symmetric(t):
    """
    Returns True if the root of t is a symmetric node, and False otherwise. If t is a leaf, it returns True:
    ex falso quodlibet.
    :return: `bool` instance.
    """
    return t.is_leaf() or \
        all(isomorphic(t.children[0], ch) for ch in t.children[1:])


def count_symmetries(t):
    """
    Returns the number of symmetric interior nodes in t.
    :return: `int` instance.
    """
    if t.is_leaf():
        return 0
    else:
        return int(is_symmetric(t)) + sum(count_symmetries(ch) for ch in t.children)


def count_automorphisms(t):
    if t.is_leaf():
        return 1

    aut = 1
    cur_sym_class_rep = None
    cur_sym_class_aut = 1
    cur_sym_class_len = 1

    # def iter_len(it):
    #     return sum(1 for _ in it)

    # def node_aut(child, iso_class):
    #     iso_class_len = iter_len(iso_class)
    #     return count_automorphisms(child)**iso_class_len * factorial(iso_class_len)

    # return prod(starmap(node_aut, groupby(t.children)))

    for i in range(len(t.children)):
        ti = t.children[i]

        if cur_sym_class_rep is None or not isomorphic(cur_sym_class_rep, ti):
            aut *= cur_sym_class_aut**cur_sym_class_len * factorial(cur_sym_class_len)

            cur_sym_class_rep = ti
            cur_sym_class_aut = count_automorphisms(ti)
            cur_sym_class_len = 1
        else:
            cur_sym_class_len += 1

    aut *= cur_sym_class_aut**cur_sym_class_len * factorial(cur_sym_class_len)
    return aut


def min_automorphisms(n):
    """
    Returns all the `Shape` instances that attain the minimum number of automorphisms with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    return [comb(n)]


def max_automorphisms(n):
    """
    Returns all the `Shape` instances that attain the maximum number of automorphisms with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    return [star(n)]
