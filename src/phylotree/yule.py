from __future__ import division
from sympy import factorial
from shape import sorted_tree
from shape.generator import collapse_tree_prob_list
from phylotree.generator import duplicate_leaf, relabellings
from phylotree import PhyloTree


def yule_from_t(t, prob):

    lvs = t.get_leaves()

    n = len(lvs)

    tps = [(duplicate_leaf(t, l, str(n+1)), lambda: prob()/n) for l in lvs]

    return collapse_tree_prob_list(tps)


def pseudo_yule(n):
    if n < 1:
        return
    elif n == 1:
        return [(PhyloTree("1", None), lambda: 1)]
    elif n == 2:
        return [(PhyloTree(None, [PhyloTree("1", None), PhyloTree("2", None)]), lambda: 1)]
    else:
        tps = [(sorted_tree(t2), lambda prob2=prob2: prob2())
               for t1, prob1 in pseudo_yule(n-1)
               for t2, prob2 in yule_from_t(t1, prob1)]
        return collapse_tree_prob_list(tps)


def yule(n):
    tps = pseudo_yule(n)
    fact = factorial(n)
    return collapse_tree_prob_list([(sorted_tree(t), lambda prob=prob, t=t: prob() * 1 / fact)
                                    for t1, prob in tps
                                    for t in relabellings(t1)])



