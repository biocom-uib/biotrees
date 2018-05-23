"""
This file contains several functions that generate `PhyloTree` instances.
"""

from itertools import permutations

from biotrees.util import and_then, iter_merge, unique
from biotrees.phylotree import PhyloTree


def add_leaf_to_edge(t, leaf_id):
    """
    Returns a `PhyloTree` instance with a new root; both a new leaf and the input `PhyloTree` pend from it.
    :param t, leaf_id: `PhyloTree` instance, unicode instance.
    :return: `PhyloTree` instance.
    """
    new_leaf = PhyloTree(leaf_id)

    if t < new_leaf:
        return PhyloTree(None, [t, PhyloTree(leaf_id)])
    else:
        return PhyloTree(None, [PhyloTree(leaf_id), t])

def add_leaf_to_node(t, leaf_id):
    """
    Returns a `PhyloTree` instance with the same root; we have added a pending leaf to it.
    :param t, leaf_id: `PhyloTree` instance, unicode instance.
    :return: `PhyloTree` instance.
    """
    if t.is_leaf():
        return add_leaf_to_edge(t, leaf_id)
    else:
        return PhyloTree(None, list(iter_insert_tree(t.children, PhyloTree(leaf_id))))

def iter_insert_tree(ts, t):
    yield from iter_merge_forests_sorted(ts, [t])

def iter_merge_forests_sorted(ts, ts2):
    yield from iter_merge(ts, ts2)

def delete_nodes_with_out_degree_one(t):
    """
    Deletes all nodes in the input `PhyloTree` instance with only one leave, for they are redundant
    :return: `PhyloTree` instance
    """
    if t.is_leaf():
        return t
    else:
        if len(t.children) == 1:
            return delete_nodes_with_out_degree_one(t.children[0])
        else:
            return PhyloTree(None, [delete_nodes_with_out_degree_one(ch) for ch in t.children])

def duplicate_leaf(t, l, newl):
    """
    Adds a new leaf as the sibling of a given leaf.
    :param t: `PhyloTree` instance.
    :param l: `PhyloTree` instance.
    :param newl: `PhyloTree` instance.
    :return: `PhyloTree` instance.
    """
    l_leaf = l.leaf

    def recurse(t, l, newl):
        if t.is_leaf():
            return add_leaf_to_edge(t, newl) if t.leaf == l_leaf else t
        else:
            return PhyloTree(None, sorted([recurse(ch, l, newl) for ch in t.children]))

    return recurse(t, l, newl)


def relabel(t, rlbl):
    """
    Renames the leaves in the input tree with a given function.
    :param t: `PhyloTree` instance.
    :param rlbl: `list` instance.
    :return: `PhyloTree` instance.
    """
    if t.is_leaf():
        l2 = rlbl.get(t.leaf, None)
        if l2 is not None:
            return PhyloTree(l2)
        return t
    else:
        return PhyloTree(None, sorted([relabel(ch, rlbl) for ch in t.children]))


def relabellings(t):
    """
    Returns a list with all possible permutations of the names of the labels of t.
    Repetitions shall arise due to symmetries in the tree.
    This method has factorial complexity.
    :param t: `PhyloTree` instance.
    :return: `list` instance.
    """
    names = list(get_leaves_names_set(t))

    for perm in permutations(names):
        yield relabel(t, dict(zip(names, perm)))
