from functools import lru_cache
from math import factorial
import random

from biotrees.phylotree import PhyloTree, count_leaves, get_leaves
from biotrees.phylotree.generator import duplicate_leaf, relabellings
from biotrees.util import parametric_total_probabilities, and_then


def sim_yule_from_t(t):
    lvs = get_leaves(t)
    n = len(lvs)

    return duplicate_leaf(t, random.choice(lvs), str(n+1))

def sim_yule(n):
    t = PhyloTree("1")

    for i in range(n-1):
        t = sim_yule_from_t(t)

    return t


@and_then(parametric_total_probabilities)
def yule_from_t(t, prob):
    lvs = get_leaves(t)
    n = len(lvs)

    for l in lvs:
        yield duplicate_leaf(t, l, str(n+1)), lambda *p: prob(*p)/n


@lru_cache(maxsize=-1)
def pseudo_yule(n):
    if n < 1:
        pass
    elif n == 1:
        yield PhyloTree("1"), lambda *p: 1
    elif n == 2:
        yield PhyloTree(None, [PhyloTree("1"), PhyloTree("2")]), lambda *p: 1
    else:
        for t1, prob1 in pseudo_yule(n-1):
            for t2, prob2 in yule_from_t(t1, prob1):
                yield t2, lambda *p, prob2=prob2: prob2(*p)


@lru_cache(maxsize=-1)
@and_then(parametric_total_probabilities)
def yule(n):
    n_factorial = factorial(n)

    for t1, prob in pseudo_yule(n):
        for t in relabellings(t1):
            yield t, lambda *p, prob=prob: prob(*p)/n_factorial

