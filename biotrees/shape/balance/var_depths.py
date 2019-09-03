from statistics import variance
from math import log, floor

from biotrees.combinatorics import subsets_with_k_elements_that_contain_subset_s
from biotrees.shape import Shape, get_leaf_depths, get_depth
from biotrees.shape.generator import binary_max_balanced


def var_depths(t):
    return variance(get_leaf_depths(t))


def min_var_depths(n):
    m = int(floor(log(n, 2)))
    twotod = 2 ** (m + 1)
    k = n - twotod / 2
    if k == 0:
        return 0, binary_max_balanced(twotod / 2)
    else:
        var, ls = min_var_depths_vector(n)
        t = binary_max_balanced(twotod)
        t = _prune_lis(t, ls)
        t = _prune_m1(n, t, ls)

        return var, t


def min_var_depths_vector(n):
    m = int(floor(log(n, 2)))
    twotom = 2**m
    k = n - twotom

    if k == 0:
        return 0, []
    else:
        var = 2*k*(twotom - k) / n**2
        ls = []

        for lis in _liss(n):
            if sum(2**li - 1 for li in lis) <= twotom - k and max(lis) < m:
                var2 = (twotom - k - sum(2**li - li**2 - 1 for li in lis)) / n         \
                        - (twotom - k - sum(2**li - li - 1 for li in lis)**2) / n**2

                if var2 < var:
                    var = var2
                    ls = lis

            elif sum(2**(li-1) - 1 for li in lis) > twotom - k:
                var2 = (3 * twotom - k - sum(2**li - li**2 - 1 for li in lis)) / n        \
                        - (3 * twotom - k - sum(2**li - li - 1 for li in lis)**2) / n**2
                if var2 < var:
                    var = var2
                    ls = lis

        return var, ls


def _liss(n):     # meterlo en utils en combinatoria
    m = int(floor(log(n, 2)))
    liss = []
    for j in range(1, m+1):
        for lis in subsets_with_k_elements_that_contain_subset_s(range(5, m+1), j, []):
            liss.append(lis)
    return liss


def _prune_m1(n, t, lis):
    m = int(floor(log(n, 2)))
    m1 = 2**(m+1) - sum(2**li - 1 for li in lis) - n
    t2 = t.clone()
    t2.children[1] = binary_max_balanced(2**m - m1)

    return t2


def _prune_lis(t, lis):

    def go(t, dis):
        if dis:
            if dis[0] != 0:
                return Shape([go(t.children[0], [di - 1 for di in dis]), t.children[1]])
            else:
                return Shape([Shape.LEAF, go(t.children[1], [di - 1 for di in dis[1:]])])

        else:
            return t

    d = get_depth(t)
    dis = tuple(reversed(tuple(d - li - 1  for li in lis)))
    return go(t, dis)

