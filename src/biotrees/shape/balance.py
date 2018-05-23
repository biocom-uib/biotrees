from functools import lru_cache
from math import factorial

from biotrees.shape import Shape
from biotrees.shape.iso import isomorphic

from biotrees.util import binom2



def is_symmetric(t):
    """
    Returns True if the root of t is a symmetric node, and False otherwise. If t is a leaf, it returns True:
    ex falso quodlibet.
    :return: `bool` instance.
    """
    return t.is_leaf() or \
        all(isomorphic(t.children[0], ch) for ch in t.children[1:])

def count_symmetries(t):
    """
    Returns the number of symmetric interior nodes in t.
    :return: `int` instance.
    """
    if t.is_leaf():
        return 0
    else:
        return int(is_symmetric(t)) + sum(count_symmetries(ch) for ch in t.children)

def count_automorphisms(t):
    if t.is_leaf():
        return 1

    aut = 1
    cur_sym_class_rep = None
    cur_sym_class_aut = 1
    cur_sym_class_len = 1

    # def iter_len(it):
    #     return sum(1 for _ in it)

    # def node_aut(child, iso_class):
    #     iso_class_len = iter_len(iso_class)
    #     return count_automorphisms(child)**iso_class_len * factorial(iso_class_len)

    # return prod(starmap(node_aut, groupby(t.children)))

    for i in range(len(t.children)):
        ti = t.children[i]

        if cur_sym_class_rep is None or not isomorphic(cur_sym_class_rep, ti):
            aut *= cur_sym_class_aut**cur_sym_class_len * factorial(cur_sym_class_len)

            cur_sym_class_rep = ti
            cur_sym_class_aut = count_automorphisms(ti)
            cur_sym_class_len = 1
        else:
            cur_sym_class_len += 1

    aut *= cur_sym_class_aut**cur_sym_class_len * factorial(cur_sym_class_len)
    return aut


def sackin_index(tree):
    def go(t):
        if t.is_leaf():
            return 0, 1

        sackins, kappas = zip(*map(go, t.children))
        node_kappa = sum(kappas)
        node_sackin = sum(sackins) + node_kappa
        return node_sackin, node_kappa

    return go(tree)[0]


def binary_colless_index(tree):
    def go(t):
        if t.is_leaf():
            return 0, 1

        left, right = t.children
        (cil, nl), (cir, nr) = go(left), go(right)
        return abs(nl - nr) + cil + cir, nl + nr

    return go(tree)[0]


def cophenetic_index(tree):
    def go(t):
        if t.is_leaf():
            return 0, 1

        cophs, kappas = zip(*map(go, t.children))
        kappa = sum(kappas)
        coph = binom2(kappa) + sum(cophs)
        return coph, kappa

    if tree.is_leaf():
        return 0
    else:
        return sum(go(ch)[0] for ch in tree.children)


def binary_quartet_index(tree):
    def go(t):
        if t.is_leaf():
            return 0, 1

        ts = t.children
        quartets, kappas = zip(*map(go, ts))
        kappa = sum(kappas)

        if kappa < 4:
            return 0, kappa

        s0 = sum(quartets)

        s3 = sum(binom2(kappas[i1]) * binom2(kappas[i2])
                    for i1 in range(len(ts))
                    for i2 in range(i1+1, len(ts)))

        return s0+s3, kappa

    return go(tree)[0]


@lru_cache(maxsize=1)
def get_quartets():
    q0 = Shape([Shape.LEAF, Shape([Shape.LEAF, Shape([Shape.LEAF, Shape.LEAF])])])
    q1 = Shape([Shape.LEAF, Shape.LEAF, Shape([Shape.LEAF, Shape.LEAF])])
    q2 = Shape([Shape.LEAF, Shape([Shape.LEAF, Shape.LEAF, Shape.LEAF])])
    q3 = Shape([Shape([Shape.LEAF, Shape.LEAF]), Shape([Shape.LEAF, Shape.LEAF])])
    q4 = Shape([Shape.LEAF, Shape.LEAF, Shape.LEAF, Shape.LEAF])
    return [q0, q1, q2, q3, q4]


def quartet_index(tree, vs = range(5)):
    Q = get_quartets()
    t0 = Shape([Shape.LEAF, Shape.LEAF, Shape.LEAF])

    def go(t):
        if t.is_leaf():
            return 0, 0, 1

        ts = t.children
        quartets, triples, kappas = zip(*map(go, ts))
        kappa = sum(kappas)

        if kappa < 3:
            triple = 0
        elif isomorphic(t, t0):
            triple = 1
        else:
            t_s0 = sum(triples)

            t_s1 = sum(kappas[i1] * kappas[i2] * kappas[i3]
                        for i1 in range(len(ts))
                        for i2 in range(i1+1, len(ts))
                        for i3 in range(i2+1, len(ts)))

            triple = t_s0 + t_s1

        if kappa < 4:
            return 0, triple, kappa

        for q, v in zip(Q, vs):
            if isomorphic(t, q):
                return v, triple, kappa

            s0 = sum(quartets)

            s1 = sum(binom2(kappas[i1]) * kappas[i2] * kappas[i3] +
                     binom2(kappas[i2]) * kappas[i1] * kappas[i3] +
                     binom2(kappas[i3]) * kappas[i1] * kappas[i2]
                        for i1 in range(len(ts))
                        for i2 in range(i1+1, len(ts))
                        for i3 in range(i2+1, len(ts)))

            s2 = sum(kappas[i1] * triples[i2] + kappas[i2] * triples[i1]
                        for i1 in range(len(ts))
                        for i2 in range(i1+1, len(ts)))

            s3 = sum(binom2(kappas[i1]) * binom2(kappas[i2])
                        for i1 in range(len(ts))
                        for i2 in range(i1+1, len(ts)))

            s4 = sum(kappas[i1] * kappas[i2] * kappas[i3] * kappas[i4]
                        for i1 in range(len(ts))
                        for i2 in range(i1+1, len(ts))
                        for i3 in range(i2+1, len(ts))
                        for i4 in range(i3+1, len(ts)))

            return s0 + vs[1]*s1 + vs[2]*s2 + vs[3]*s3 + vs[4]*s4, triple, kappa

    return go(tree)[0]
