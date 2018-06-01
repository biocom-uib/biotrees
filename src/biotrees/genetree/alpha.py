from sympy import var

from sys import setrecursionlimit

from biotrees.phylotree import get_leaves
from biotrees.genetree import GeneTree
from biotrees.genetree.generator import add_leaf_to_edge, duplicate_leaf, relabellings
from biotrees.genetree.newick import to_newick
from biotrees.genetree.iso import quotient
from biotrees.util import and_then, parametric_total_probabilities, groupby_sorted, iter_replace_at, skip_nth, iter_merge


setrecursionlimit(2000)


@and_then(parametric_total_probabilities(grouper=quotient))
def alpha_from_t(N, t, prob):
    n = len(get_leaves(t))

    def recurse(N, t, prob, n):
        if t.is_leaf():

            for i in range(N):
                yield add_leaf_to_edge(t, str(n + 1), str(i + 1)), lambda a: prob(a) * (1 - a) / (n - a) / N

        else:

            for i in range(N):
                yield add_leaf_to_edge(t, str(n + 1), str(i + 1)), lambda a: prob(a) * a / (n - a) / N

            i = 0

            for ch in t.children:
                children = list(skip_nth(t.children, i))
                chps = recurse(N, ch, prob, n)
                yield from [(GeneTree(None, None, list(iter_merge(children, [ch]))), p) for ch, p in chps]

                i += 1

    yield from recurse(N, t, prob, n)


@and_then(parametric_total_probabilities(grouper=quotient))
def pseudo_alpha(n, N):
        if n < 1:
            pass
        elif n == 1:
            for i in range(N):
                yield GeneTree("1", str(i+1)), lambda a: 1/N
        else:
            for t1, prob1 in pseudo_alpha(n-1, N):
                for t2, prob2 in alpha_from_t(N, t1, prob1):
                    yield t2, lambda a, prob=prob2: prob(a)


def alpha(n, N):
    return pseudo_alpha(n, N)





















def shape_invariant_list(tps):
    print(len(tps))
    i = 0
    while i < len(tps):
        print(i)
        j = i+1
        while j < len(tps):
            if iso(tps[i][0], tps[j][0]):
                p = var('p')
                if tps[i][1](p) != tps[j][1](p):
                    print(to_newick(tps[i][0]))
                    print(to_newick(tps[j][0]))
                    print(tps[i][1](p))
                    print(tps[j][1](p))
                    return False
                else:
                    j += 1
            else:
                j += 1
        i += 1

    return True


def shape_invariant(n, func):
    return shape_invariant_list(yule(n, func))



def sampling_consistent(n, marker):
    # la func debería ser el yule_from_t
    if marker == 1:
        func = pseudo_yule1
    elif marker == 2:
        func = pseudo_yule2
    else:
        return
    tpsn  = yule(n, func)
    #print(len(tpsn))

    tpsn1 = yule(n-1, func)

    #print(len(tpsn1))
    p = var('p')

    for t1, prob1 in tpsn1:
        if marker == 1:
            ts = [t for t, prob in yule_from_t1(t1, prob1)]
        elif marker == 2:
            ts = [t for t, prob in yule_from_t2(t1, prob1)]
        else:
            return
        print(len(ts))

        totprob = lambda p: 0

        for t, prob in tpsn:
            if t in ts:
                totprob = lambda p: totprob(p) + prob(p)

        if prob1(p) != totprob(p):
            return False
    return True


def debug(tps, n):
    lbls = [str(i+1) for i in range(n)]

    tps = [(t, lambda p, prob=prob2, t=t: prob2(p) * 1 / len(relabellings(t, lbls)))
           for t1, prob2 in tps
           for t in relabellings(t1, lbls)]

    ts = [tp[0] for tp in tps]
    reps = [(t, ts.count(t)) for t in ts]

    for t1, count1 in reps:
        for t2, count2 in reps:
            if iso(t1, t2) and count1 != count2:
                print(to_newick(t1))
                print(count1)
                print(to_newick(t2))
                print(count2)
                return True
    return False
# relabellings of ((1,2), (3,4)) have errors

"""
    lbls = [str(i+1) for i in range(n)]
    for t, prob in lst:
        relbls = [(ti, lambda p, prob=prob, t=t: prob(p) * 1 / len(relabellings(t, lbls)))
                          for ti in relabellings(t, lbls)]
        # check number of repetitions of each element
        reps = []
        for i in range(len(relbls)):
            t = relbls[i][0]
            if reps:
                rep = False
                for i in range(len(reps)):
                    numreps = reps[i][1]
                    print(reps[i][1])
                    if t == reps[i][0]:
                        numreps += 1
                        rep = True
                if rep:
                    reps[i] = ()
            else:
                reps.append((t, 1))
"""