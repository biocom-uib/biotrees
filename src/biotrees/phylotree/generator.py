"""
This file contains several functions that generate `PhyloTree` instances.
"""

from phylotree import PhyloTree
from shape import sorted_tree
from combinatorics import permutations


def add_leaf_to_edge(t, leaf_id):
    """
    Returns a `PhyloTree` instance with a new root; both a new leaf and the input `PhyloTree` pend from it.
    :param t, leaf_id: `PhyloTree` instance, unicode instance.
    :return: `PhyloTree` instance.
    """
    return PhyloTree(None, [PhyloTree(leaf_id), t])


def add_leaf_to_node(t, leaf_id):
    """
    Returns a `PhyloTree` instance with the same root; we have added a pending leaf to it.
    :param t, leaf_id: `PhyloTree` instance, unicode instance.
    :return: `PhyloTree` instance.
    """
    if t.is_leaf:
        return add_leaf_to_edge(t, leaf_id)
    else:
        ch = t.children[:]
        ch[:0] = [PhyloTree(leaf_id)]
        return PhyloTree(None, ch)


def duplicate_leaf(t, l, newl):
    """
    Adds a new leaf as the sibling of a given leaf.
    :param t: `PhyloTree` instance.
    :param l: `PhyloTree` instance.
    :param newl: `PhyloTree` instance.
    :return: `PhyloTree` instance.
    """
    n = t.count_leaves()

    def recurse(t, l, newl, n):
        if t.leaf == l.leaf:
            return add_leaf_to_edge(t, newl)
        elif t.is_leaf and t.leaf != l.leaf:
            return t
        else:
            genet = PhyloTree(None, [recurse(ch, l, newl, n) for ch in t.children])
            return sorted_tree(genet)

    return recurse(t, l, newl, n)


def relabel(t, rlbl):
    """
    Renames the leaves in the input tree with a given function.
    :param t: `PhyloTree` instance.
    :param rlbl: `list` instance.
    :return: `PhyloTree` instance.
    """
    if t.is_leaf:
        for lbl1, lbl2 in rlbl:
            if t.leaf == lbl1:
                return PhyloTree(lbl2)
        return PhyloTree(t.leaf)
    else:
        return PhyloTree(None, [relabel(ch, rlbl) for ch in t.children])


def relabellings(t):
    """
    Returns a list with all possible permutations of the names of the labels of t.
    Repetitions shall arise due to symmetries in the tree.
    This method has factorial complexity.
    :param t: `PhyloTree` instance.
    :return: `list` instance.
    """
    perms = permutations(t.leaves_names_set())
    return [relabel(t, perm) for perm in perms]



def delete_nodes_with_out_degree_one(t):
    """
    Deletes all nodes in the input `PhyloTree` instance with only one leave, for they are redundant
    :return: `PhyloTree` instance
    """
    if t.is_leaf:
        return t
    else:
        if len(t.children) == 1:
            return delete_nodes_with_out_degree_one(t.children[0])
        else:
            return PhyloTree(None, [delete_nodes_with_out_degree_one(t) for t in t.children])