from functools import lru_cache

from biotrees.util import unique
from biotrees.shape import Shape


def binary(n):
    """
    Returns the binary representation of an `int` n in `String`.
    :param n: `int` instance.
    :return: `String` instance.
    """
    return "{0:b}".format(n)


def maxpower2(n):
    """
    Returns the maximum power of 2 that divides an `int` n.
    :param n: `int` instance.
    :return: `int` instance.
    """
    return last1(binary(n))


def last1(bn):
    """
    Returns the maximum power of 2 that divides an `int` n, given in binary `String` format.
    :param n: `String` instance.
    :return: `int` instance.
    """
    return len(bn) - 1 - max(filter(lambda i : bn[i] == '1', range(0, len(bn))))


@lru_cache(maxsize=None)
def mincolless(n):
    """
    Returns all the `Shape` instances that attain the minimum Colless index with `int` n leaves.
    :param n: `int` instance.
    :return: `list` instance.
    """

    if n == 0:
        return []
    elif n == 1:
        return [Shape.LEAF]
    elif n == 2:
        return [Shape.CHERRY]
    else:
        tss = []
        for n1, n2 in mincolless_root(n):
            ts = unique([Shape(sorted([t1, t2])) for t1 in mincolless(n1)
                         for t2 in mincolless(n2)])
            tss = tss + ts
        return tss


@lru_cache(maxsize=None)
def mincolless_root(n):
    """
    Returns all the possible distributions of the number of children in the root so that
    a tree with `int` n leaves can attain the minimum Colless index.
    :param n: `int` instance.
    :return: `list` instance.
    """

    ns = []

    if n % 2 == 0:
        ns.append((n // 2, n // 2))

    k = maxpower2(n)

    oddn = n // 2 ** k

    bn = binary(oddn)

    ls = list(map(lambda x : len(bn) - 1 - x, filter(lambda i : bn[i] == '1', range(0, len(bn)))))

    if oddn == 3:                   # this case must be treated apart
        ns.append((2**k, 2**(k+1)))

    elif len(ls) > 1:
        for i in range(1, len(ls)-1):
            if ls[i+1] < ls[i] - 1 or i == len(ls) - 2:
                l = ls[i]
                t = sum(2**ls[j] for j in range(i+1, len(ls)-1)) // 2
                p = (oddn - 2 ** l - 2 * t - 1) // 2 ** (l + 1)

                ns.append((2**k*(2**l*p + 2*t + 1), 2**(k + l)*(p + 1)))

        for i in range(1, len(ls)):
            if ls[i-1] > ls[i] + 1:
                l = ls[i] + 1
                p = (oddn - sum(2 ** ls[j] for j in range(len(ls) - 1, i, -1))) // 2 ** (l + 1)
                t = (2 ** (l+1) * p + 2 ** l - 1 - oddn) // 2

                ns.append((2**(k + l)*p, 2**k*(2**l*(p + 1) - 2*t - 1)))

    return ns
