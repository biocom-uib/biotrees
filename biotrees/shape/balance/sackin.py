"""
This computes the extremal values of the Sackin index for any int n following the algorithm by Mareike Fischer in their
paper "Extremal values of the Sackin balance index for rooted binary trees", as well as the value of the index in a
given tree.
"""

from biotrees.util import unique

from biotrees.shape import Shape, count_leaves
from biotrees.shape.generator import comb, binary_max_balanced

from biotrees.phylotree import PhyloTree, shape_to_phylotree, phylotree_to_shape


def log2(n):
    k = int.bit_length(n) - 1

    if n == 2**k:
        return k
    else:
        return k+1


def sackin_index(tree):
    def go(t):
        if t.is_leaf():
            return 0, 1

        sackins, kappas = zip(*map(go, t.children))
        node_kappa = sum(kappas)
        node_sackin = sum(sackins) + node_kappa
        return node_sackin, node_kappa

    return go(tree)[0]


def max_sackin(n):
    """
    Returns all the `Shape` instances that attain the maximum Sackin index with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    return [comb(n)]


def min_sackin(n):
    """
    Returns all the `Shape` instances that attain the minimum Sackin index with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """
    k = log2(n)
    nk = 2**k

    if n == nk:
        return [binary_max_balanced(nk)]
    else:
        ts = []

        for t in min_sackin(n + 1):
            ts = ts + cherry_picking(t, k)

        return unique(sorted(ts))


def cherry_picking(tree, k):
    if tree.is_leaf():
        return [tree]
    elif tree == Shape.CHERRY:
        return [Shape.LEAF]
    else:

        def go(t, i, d):    # for binary eyes only

            if d != k and not (t.is_leaf() or t.is_cherry()):
                chs = [go(ch, i, d + 1) for ch in t.children]
                if any(ch is None for ch in chs):
                    return None
                else:
                    return PhyloTree(None, sorted(chs))

            elif d != k and (t == PhyloTree(None, [PhyloTree(i), PhyloTree(i + 1)])
                             or t == PhyloTree(None, [PhyloTree(i - 1), PhyloTree(i)])):
                return None

            elif d == k and (t == PhyloTree(None, [PhyloTree(i), PhyloTree(i + 1)])
                             or t == PhyloTree(None, [PhyloTree(i - 1), PhyloTree(i)])):
                return PhyloTree(i)

            elif t == PhyloTree(i):
                return None

            else:
                return t

        n = count_leaves(tree)
        phyl = shape_to_phylotree(tree, gen=int)

        ts = []

        for i in range(n):
            t = go(phyl, i, 1)

            if t:
                ts.append(phylotree_to_shape(t))

        return unique(sorted(ts))
