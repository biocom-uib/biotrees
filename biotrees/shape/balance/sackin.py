from biotrees.util import unique

from biotrees.shape import Shape, count_leaves
from biotrees.shape.generator import comb, binary_max_balanced

from biotrees.phylotree import PhyloTree, shape_to_phylotree, phylotree_to_shape
from biotrees.phylotree.subtree import subtree


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


def minsackin(n):   # horrible e hiper costoso, esta Fischer...
    k = log2(n)
    nk = 2**k

    if n == nk:
        return [binary_max_balanced(nk)]
    else:
        ts = []

        for t in minsackin(n+1):
            ts = ts + cherrypicking(t, k)

        return ts


def cherrypicking(tree, k):     # y si encuentra la hoja i pero no es una cherry a max prof?
                                # devolver también un booleano?
                                # trabajar con maybe
    if tree.is_leaf():
        return [tree]
    elif tree == Shape.CHERRY:
        return [Shape.LEAF]
    else:

        def go(t, i, d):

            if d != k and not (t.is_leaf() or t.is_cherry()):
                chs = [go(ch, i, d + 1) for ch in t.children]
                if any(ch is None for ch in chs):
                    return None
                else:
                    return PhyloTree(None, chs)

            elif d != k and (t == PhyloTree(i)
                             or t == PhyloTree(None, [PhyloTree(i), PhyloTree(i + 1)])
                             or t == PhyloTree(None, [PhyloTree(i - 1), PhyloTree(i)])):
                return None

            elif d == k and (t == PhyloTree(None, [PhyloTree(i), PhyloTree(i + 1)])
                             or t == PhyloTree(None, [PhyloTree(i - 1), PhyloTree(i)])):
                return PhyloTree(i)

            else:
                return t

        n = count_leaves(tree)
        phyl = shape_to_phylotree(tree, gen=int)

        ts = []

        for i in range(n-1):
            t = go(phyl, i, 1)

            if t: ts.append(t)

        return unique(sorted(ts))


def go(t, i, d, k):        # I think we are assuming the entry is a binary tree
    print(i)

    if d != k and not (t.is_leaf() or t.is_cherry()):
        chs = [go(ch, i, d+1, k) for ch in t.children]
        if any(ch is None for ch in chs):
            return None
        else:
            return PhyloTree(None, chs)

    elif d != k and (t == PhyloTree(i)
            or t == PhyloTree(None, [PhyloTree(str(i)), PhyloTree(str(i+1))])
            or t == PhyloTree(None, [PhyloTree(str(i-1)), PhyloTree(str(i))])):
        print("ahoy")
        return None

    elif d == k and (t == PhyloTree(None, [PhyloTree(str(i)), PhyloTree(str(i+1))])
            or t == PhyloTree(None, [PhyloTree(str(i-1)), PhyloTree(str(i))])):
        return PhyloTree(i)

    else:
        print("nay")
        return t