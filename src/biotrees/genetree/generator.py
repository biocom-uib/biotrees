"""
This file contains several functions that generate `genetree` instances.
"""
from itertools import combinations

from biotrees.genetree import GeneTree
#from biotrees.genetree.iso import equal

from biotrees.phylotree import count_leaves
from biotrees.phylotree.generator import iter_insert_tree, iter_merge_forests_sorted

from sympy import simplify
from biotrees.combinatorics import finite_bijections
from biotrees.util import and_then, unique_unsortable


def add_leaf_to_edge(t, leaf_id, lbl):
    """
    Returns a `genetree` instance with a new root; both a new leaf and the input `genetree` pend from it.
    :param t, leaf_id, lbl: `genetree` instance, unicode instance, unicode instance.
    :return: `genetree` instance.
    """
    new_leaf = GeneTree(leaf_id, lbl)

    if t.is_leaf and t.leaf < new_leaf.leaf:    # el orden en GeneTrees no acaba de ir muy allá
        return GeneTree(None, None, [t, new_leaf])
    else:
        return GeneTree(None, None, [new_leaf, t])


def add_leaf_to_node(t, leaf_id, lbl):
    """
    Returns a `genetree` instance with the same root; we have added a pending leaf to it.
    :param t, leaf_id, lbl: `genetree` instance, unicode instance, unicode instance.
    :return: `genetree` instance.
    """
    if t.is_leaf:
        return add_leaf_to_edge(t, leaf_id, lbl)
    else:
        return GeneTree(None, None, list(iter_insert_tree(t.children, GeneTree(leaf_id, lbl))))


def duplicate_leaf(t, l, newl):

    n = count_leaves(t)

    def recurse(t, l, newl, n):
        if t.is_leaf():
            return add_leaf_to_edge(t, str(n+1), newl) if t.leaf == l.leaf else t
        else:
            genet = GeneTree(None, None, sorted([recurse(ch, l, newl, n) for ch in t.children]))
            return genet

    return recurse(t, l, newl, n)


"""
def collapse_tree_prob_list(tps, boolfunc):
    i = 0
    while i < len(tps):
        j = i+1
        ps = [tps[i][1]]
        while j < len(tps):
            if boolfunc(tps[i][0], tps[j][0]):
                ps.append(tps[j][1])
                tps.pop(j)
            else:
                j += 1
        prob = lambda *args, ps=ps: simplify(sum(p(*args) for p in ps))
        tps[i] = (tps[i][0], prob)
        i += 1
    return tps
"""


def cherry(lf1, lf2, lbl1, lbl2):
    return GeneTree(None, None, [GeneTree(lf1, lbl1), GeneTree(lf2, lbl2)])


def relabel(t, rlbl):
    if t.is_leaf():
        lbl2 = rlbl.get(t.label, None)
        if lbl2 is not None:
            return GeneTree(t.leaf, lbl2)
        return GeneTree(t.leaf, t.label)
    else:
        return GeneTree(None, None, sorted([relabel(ch, rlbl) for ch in t.children]))

@and_then(unique_unsortable)
def relabellings(t, s):
    lbls = t.label_set()
    if len(s) < len(lbls):
        return []
    else:
        relbls = []
        for subset in combinations(s, len(lbls)):
            relbls = relbls + list(finite_bijections(lbls, subset))
        #return unique_unsortable([relabel(t, relbl) for relbl in relbls])
        return [relabel(t, relbl) for relbl in relbls]


def delete_nodes_with_out_degree_one(t):
    if t.is_leaf:
        return t
    else:
        if len(t.children) == 1:
            return t.children[0].delete_nodes_with_out_degree_one()
        else:
            return GeneTree(None, None, [delete_nodes_with_out_degree_one(t) for t in t.children])


