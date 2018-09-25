from itertools import combinations

from biotrees.shape.newick import from_newick
from biotrees.phylotree import shape_to_phylotree
from biotrees.phylotree.subtree import subtree

q0 = from_newick("(*,(*,(*,*)));")
q1 = from_newick("(*,*,(*,*));")
q2 = from_newick("(*,(*,*,*));")
q3 = from_newick("((*,*),(*,*));")
q4 = from_newick("(*,*,*,*);")
Q = [q0, q1, q2, q3, q4]


def binom2(n):
    return n*(n-1)/2


def qi(t, vs = list(range(len(Q)))):
    if t.count_leaves() < 4:
        return 0

    for q, v in zip(Q, vs):
        if t.iso(q):
            return v

    ts = t.children
    s0 = sum(qi(ch, vs) for ch in ts)
    s1 = sum(binom2(ts[i1].count_leaves()) * ts[i2].count_leaves() * ts[i3].count_leaves()
             + binom2(ts[i2].count_leaves()) * ts[i1].count_leaves() * ts[i3].count_leaves()
             + binom2(ts[i3].count_leaves()) * ts[i1].count_leaves() * ts[i2].count_leaves()
             for i1 in range(len(ts))
             for i2 in range(i1 + 1, len(ts))
             for i3 in range(i2 + 1, len(ts)))
    s2 = sum(ts[i1].count_leaves() * triples(ts[i2]) + ts[i2].count_leaves() * triples(ts[i1])
                      for i1 in range(len(ts))
                      for i2 in range(i1 + 1, len(ts)))
    s3 = sum(binom2(ts[i1].count_leaves()) * binom2(ts[i2].count_leaves())
             for i1 in range(len(ts))
             for i2 in range(i1 + 1, len(ts)))
    s4 = sum(ts[i1].count_leaves() * ts[i2].count_leaves() * ts[i3].count_leaves() * ts[i4].count_leaves()
                      for i1 in range(len(ts))
                      for i2 in range(i1 + 1, len(ts))
                      for i3 in range(i2 + 1, len(ts))
                      for i4 in range(i3 + 1, len(ts)))

    return s0 + qi(q1, vs) * s1 + qi(q2, vs) * s2 + qi(q3, vs) * s3 + qi(q4, vs) * s4


def qib2(t):
    if t.count_leaves() < 4:
        return 0
    else:
        ts = t.children
        return (binom2(ts[0].count_leaves()) * binom2(ts[1].count_leaves())
                + qib2(ts[0]) + qib2(ts[1]))

def qib(t):
    # for binary eyes only
    if t.count_leaves() < 4:
        return 0
    else:
        ts = t.children
        s0 = sum(qib(ch) for ch in ts)
        s3 = sum(binom2(ts[i1].count_leaves()) * binom2(ts[i2].count_leaves())
                 for i1 in range(len(ts))
                 for i2 in range(i1 + 1, len(ts)))

        return s0 + s3


def triples(t):
    if t.count_leaves() < 3:
        return 0
    elif t.iso(from_newick("(*,*,*);")):
        return 1
    else:
        ts = t.children
        s0 = sum(triples(ch) for ch in ts)
        s1 = sum(ts[i1].count_leaves() * ts[i2].count_leaves() * ts[i3].count_leaves()
                 for i1 in range(len(ts))
                 for i2 in range(i1 + 1, len(ts))
                 for i3 in range(i2 + 1, len(ts)))
        return s0 + s1


def qi_def(t):
    qs = combinations(range(t.count_leaves()), 4)
    phylot = shape_to_phylotree(t)
    return sum(qi(subtree(phylot, q)) for q in qs)
