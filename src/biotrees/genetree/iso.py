from itertools import groupby, chain

from biotrees.combinatorics import finite_bijections
from biotrees.genetree.generator import relabel
from biotrees.genetree import genetree_to_shape
from biotrees.genetree.newick import to_newick

import biotrees.shape.iso as shape


def equal(t1, t2):
    t1._sort()
    t2._sort()
    return t1 == t2


def isomorphic(t1, t2):  # horrible algorithm, O(n!)... but to be used only in small cases
    return get_iso_relabeling(t1, t2) is not None


def get_iso_relabeling(t1, t2):
    # horrible algorithm, O(n!)... but to be used only in small cases
    sh1 = t1.shape()
    sh2 = t2.shape()

    lbls1 = t1.label_set()
    lbls2 = t2.label_set()

    if shape.isomorphic(sh1, sh2) and len(lbls1) == len(lbls2):
        for bij in finite_bijections(lbls1, lbls2): # aqui esta el paso diabólico
            if equal(t2, relabel(t1, bij)):
                return bij

    return None


def quotient(ts, k):

    tps = [(t, p) for t, p in ts]

    ts = [list(lst) for _, lst in groupby(sorted(tps, key=lambda x: genetree_to_shape(k(x))),
                                    key=lambda x: genetree_to_shape(k(x)))]

    return chain.from_iterable(_quotient(lst, k) for lst in ts if lst is not None)
                                        #aquí era


def _quotient(lst, k):  # here we know that all trees have the same shape

    for i in range(len(lst)): # es una criba de Eratostenes
        if lst[i] is not None:
            rep = lst[i]

            klass = [rep]

            for j in range(i+1, len(lst)):  # o debo guardar la permutación?
                if lst[j] is not None:
                    if isomorphic(k(rep), k(lst[j])):
                        klass.append(lst[j])
                        lst[j] = None

            yield k(rep), klass
