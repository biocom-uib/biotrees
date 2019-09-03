"""
This file contains several functions that generate `Shape` instances.
"""
from functools import lru_cache

from biotrees.shape import Shape
from biotrees.util import and_then, iter_merge, skip_nth, unique


def add_leaf_to_edge(t):
    """
    Returns a `Shape` instance with a new root; both a new leaf and the input `Shape` pend from it.
    :param t: `Shape` instance.
    :return: `Shape` instance.
    """
    return Shape([Shape.LEAF, t])


def add_leaf_to_node(t):
    """
    Returns a `Shape` instance with the same root; we have added a pending leaf to it.
    :param t: `Shape` instance.
    :return: `Shape` instance.
    """
    if t.is_leaf():
        return add_leaf_to_edge(t)
    else:
        return Shape([Shape.LEAF] + t.children)


def iter_insert_tree(ts, t):
    yield from iter_merge_forests_sorted(ts, [t])


def iter_merge_forests_sorted(ts1, ts2):
    yield from iter_merge(ts1, ts2)


def iter_replace_tree_at(ts, i, repl):
    yield from iter_insert_tree(skip_nth(ts, i), repl)


def delete_nodes_with_out_degree_one(t):
    """
    Deletes all nodes in the input `Shape` instance with only one child, for they are redundant
    :return: `Shape` instance
    """
    if t.is_leaf():
        return t
    else:
        if len(t.children) == 1:
            return delete_nodes_with_out_degree_one(t.children[0])
        else:
            return Shape([delete_nodes_with_out_degree_one(ch) for ch in t.children])


@and_then(unique)
def all_binary_trees_from_t(t):
    """
    Returns a list with all the binary `Shape` instances that can be obtained from the input tree.
    :param t: `Shape` instance.
    :return: `list` instance.
    """
    if not t.is_leaf():
        for i in range(len(t.children)):
            for ti in all_binary_trees_from_t(t.children[i]):
                ts2 = list(iter_replace_tree_at(t.children, i, ti))
                yield Shape(ts2)

    yield add_leaf_to_edge(t)


@and_then(unique)
def all_trees_from_t(t):
    """
    Returns a list with all the `Shape` instances that can be obtained from the input tree.
    :param t: `Shape` instance.
    :return: `list` instance.
    """
    if not t.is_leaf():
        for i in range(len(t.children)):
            for ti in all_trees_from_t(t.children[i]):
                ts2 = list(iter_replace_tree_at(t.children, i, ti))
                yield Shape(ts2)

        yield add_leaf_to_node(t)

    yield add_leaf_to_edge(t)


@lru_cache(maxsize=None)
@and_then(unique)
def all_binary_trees_with_n_leaves(n):
    """
    Returns a list with all the binary `Shape` instances that have n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    if n <= 0:
        pass
    elif n == 1:
        yield Shape.LEAF
    elif n == 2:
        yield Shape([Shape.LEAF, Shape.LEAF])
    else:
        for t in all_binary_trees_with_n_leaves(n-1):
            yield from all_binary_trees_from_t(t)


@lru_cache(maxsize=None)
@and_then(unique)
def all_trees_with_n_leaves(n):
    """
    Returns a list with all the `Shape` instances that have n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    if n <= 0:
        pass
    elif n == 1:
        yield Shape.LEAF
    elif n == 2:
        yield Shape.CHERRY
    else:
        for t in all_trees_with_n_leaves(n-1):
            yield from all_trees_from_t(t)


def star(n):
    """
    Returns a star `Shape` with n leaves.
    :param n: `int` instance.
    :return: `Shape` instance.
    """
    return Shape([Shape.LEAF for _ in range(n)])


@lru_cache(maxsize=None)
def binary_max_balanced(n):
    """
    Returns a binary maximum balanced `Shape` with n leaves.
    :param n: `int` instance.
    :return: `Shape` instance.
    """
    if n == 1:
        return Shape.LEAF
    elif n == 2:
        return Shape.CHERRY
    else:
        s = n % 2
        return Shape([binary_max_balanced((n - s) // 2), binary_max_balanced((n - s) // 2 + s)])


@lru_cache(maxsize=None)
def comb(n):
    """
    Returns a comb `Shape` with n leaves.
    :param n: `int` instance.
    :return: `Shape` instance.
    """
    if n == 1:
        return Shape.LEAF
    elif n == 2:
        return Shape.CHERRY
    else:
        return Shape([Shape.LEAF, comb(n-1)])
