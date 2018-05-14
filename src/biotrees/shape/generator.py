"""
This file contains several functions that generate `Shape` instances.
"""

from itertools import groupby
from sympy import simplify

from biotrees.shape import Shape, sorted_tree, sorted_by_shape
from biotrees.shape.iso import iso

def add_leaf_to_edge(t):
    """
    Returns a `Shape` instance with a new root; both a new leaf and the input `Shape` pend from it.
    :param t: `Shape` instance.
    :return: `Shape` instance.
    """
    return Shape([Shape(), t])


def add_leaf_to_node(t):
    """
    Returns a `Shape` instance with the same root; we have added a pending leaf to it.
    :param t: `Shape` instance.
    :return: `Shape` instance.
    """
    if t.is_leaf:
        return add_leaf_to_edge(t)
    else:
        ch = t.children[:]
        ch[:0] = [Shape()]
        return Shape(ch)


def all_binary_trees_from_t(t):
    """
    Returns a list with all the binary `Shape` instances that can be obtained from the input tree.
    :param t: `Shape` instance.
    :return: `list` instance.
    """
    if t.is_leaf:
        return [add_leaf_to_edge(t)]
    else:
        lst = []
        for i in range(len(t.children)):
            for ti in all_binary_trees_from_t(t.children[i]):
                ts = t.children[:]
                ts[i] = ti
                lst.append(sorted_tree(Shape(ts)))

        lst.append(add_leaf_to_edge(t))

        return collapse_list(lst)


def all_trees_from_t(t):
    """
    Returns a list with all the `Shape` instances that can be obtained from the input tree.
    :param t: `Shape` instance.
    :return: `list` instance.
    """
    if t.is_leaf:
        return [add_leaf_to_edge(t)]
    else:
        lst = []
        for i in range(len(t.children)):
            for ti in all_trees_from_t(t.children[i]):
                ts = t.children[:]
                ts[i] = ti
                lst.append(sorted_tree(Shape(ts)))

        lst.append(add_leaf_to_edge(t))
        lst.append(add_leaf_to_node(t))

        return collapse_list(lst)


def all_binary_trees_with_n_leaves(n):
    """
    Returns a list with all the binary `Shape` instances that have n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    if n <= 0:
        return []
    elif n == 1:
        return [Shape()]
    elif n == 2:
        return [Shape([Shape(), Shape()])]
    else:
        ts = []

        for t in all_binary_trees_with_n_leaves(n-1):
            ts = ts + all_binary_trees_from_t(t)

        return collapse_list(ts)


def all_trees_with_n_leaves(n):
    """
    Returns a list with all the `Shape` instances that have n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    if n <= 0:
        return []
    elif n == 1:
        return [Shape()]
    elif n == 2:
        return [cherry()]
    else:
        ts = []

        for t in all_trees_with_n_leaves(n-1):
            ts = ts + all_trees_from_t(t)

        return collapse_list(ts)


def cherry():
    return Shape([Shape(), Shape()])


def collapse_tree_prob_list(tps, boolfunc):
    """
    Takes a list of tuples trees and probabilities and sums the probabilities of all equal trees. Then it returns a list
    in which each tree appears only once.
    :param tps: `list` instance.
    :param boolfunc: `function` instance.
    :return: `list` instance.
    """
    tps = sorted_by_shape(tps)
    i = 0
    while i < len(tps):
        ps = [tps[i][1]]
        while i + 1 < len(tps) and boolfunc(tps[i][0], tps[i + 1][0]):
            ps.append(tps[i + 1][1])
            tps.pop(i + 1)
        prob = lambda *args, ps=ps: simplify(sum(p(*args) for p in ps))
        tps[i] = (tps[i][0], prob)
        i += 1
    return tps


def filter_by_shape(lst, boolfunc):
    lst = sorted_by_shape(lst)
    i = 0
    while i < len(lst):
        while i + 1 < len(lst) and boolfunc(lst[i][0], lst[i + 1][0]):
            lst.pop(i + 1)

        i += 1
    return lst

def collapse_list(lst):
    return [k for k, _ in groupby(lst)]
