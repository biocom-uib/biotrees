from collections import Hashable
from itertools import groupby
from functools import reduce
import operator


def binom2(m):
    return m * (m-1) // 2 if m >= 2 else 0

def skip_nth(iterable, n):
    for i, x in enumerate(iterable):
        if i != n:
            yield x

def last(iterable):
    return list(iterable)[-1]

def unique(lst, sort=False):
    """
    Checks if the first element in lst is hashable. If it is, lst is converted
    to a set and then back to a list. Otherwise it is first sorted and then
    linearly traversed skipping repeated items.
    """
    lst = list(lst)
    if lst and isinstance(lst[0], Hashable):
        return sorted(set(lst)) if sort else list(set(lst))
    else:
        return [k for k, _ in groupby(sorted(lst))]

def unique_unsortable(lst):
    ret = []
    for x in list(lst):
        if x not in ret:
            ret.append(x)
    return ret

def iter_merge(xs, ys):
    xs = iter(xs)
    ys = iter(ys)

    cur1 = next(xs, None)
    cur2 = next(ys, None)

    while cur1 is not None and cur2 is not None:
        if cur1 <= cur2:
            yield cur1
            cur1 = next(xs, None)
        else:
            yield cur2
            cur2 = next(ys, None)

    while cur1 is not None:
        yield cur1
        cur1 = next(xs, None)

    while cur2 is not None:
        yield cur2
        cur2 = next(ys, None)

def lifted_sum(fs):
    return lambda *args, **kwargs: sum(f(*args, **kwargs) for f in fs)

def lifted_prod(fs, *args, **kwargs):
    return reduce(operator.mul, (f(*args, **kwargs) for f in fs))

def parametric_total_probabilities(xps):
    """
    Takes a list of tuples (x, p) and sums the probabilities (p) of all equal x. Then it returns a list
    in which x appears only once.
    :param xps: `list` instance.
    :return: `list` instance.
    """
    xps = sorted(xps, key = lambda tp: tp[0])
    groups = groupby(xps, key = lambda tp: tp[0])

    return [(t, lifted_sum([p for _t, p in group]))
            for t, group in groups]

def and_then(then):
    def wrap(f):
        def wrapped(*args, **kwargs):
            return then(f(*args, **kwargs))
        return wrapped

    return wrap
