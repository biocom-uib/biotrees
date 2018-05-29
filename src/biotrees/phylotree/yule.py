from functools import lru_cache
from math import factorial
from itertools import groupby


from biotrees.phylotree import PhyloTree, get_leaves
from biotrees.phylotree.generator import duplicate_leaf, relabellings
from biotrees.util import parametric_total_probabilities, and_then, groupby_sorted


@and_then(parametric_total_probabilities(grouper=groupby_sorted))
def yule_from_t(t, prob):
    lvs = get_leaves(t)
    n = len(lvs)

    for l in lvs:
        yield duplicate_leaf(t, l, str(n+1)), lambda: prob() / n


@lru_cache(maxsize=-1)
def pseudo_yule(n):
    if n < 1:
        pass
    elif n == 1:
        yield PhyloTree("1", None), lambda: 1
    elif n == 2:
        yield PhyloTree(None, [PhyloTree("1"), PhyloTree("2")]), lambda: 1
    else:
        for t1, prob1 in pseudo_yule(n-1):
            for t2, prob2 in yule_from_t(t1, prob1):
                yield t2, lambda prob2=prob2: prob2()


@lru_cache(maxsize=-1)
@and_then(parametric_total_probabilities(grouper=groupby_sorted))
def yule(n):
    n_factorial = factorial(n)

    for t1, prob in pseudo_yule(n):
        for t in relabellings(t1):
            yield t, lambda prob=prob, t=t: prob() * 1 / n_factorial
