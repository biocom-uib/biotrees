from functools import lru_cache

from biotrees.shape import Shape
from biotrees.shape.iso import isomorphic
from biotrees.shape.generator import star, comb, binary_max_balanced

from biotrees.util import binom2


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


def min_quartet(n):
    """
    Returns all the `Shape` instances that attain the minimum Quartet index with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    return [comb(n)]


def max_quartet(n):
    """
    Returns all the `Shape` instances that attain the maximum Quartet index with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    return [star(n)]


def max_binary_quartet(n):
    """
    Returns all the binary `Shape` instances that attain the minimum binary Quartet index with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    return [binary_max_balanced(n)]
