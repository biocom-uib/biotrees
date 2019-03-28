
from biotrees.shape.generator import comb


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


