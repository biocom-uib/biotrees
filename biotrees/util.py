from itertools import groupby


def binom2(m):
    return m * (m-1) // 2 if m >= 2 else 0

def skip_nth(iterable, n):
    for i, x in enumerate(iterable):
        if i != n:
            yield x

def unique(lst):
    return [k for k, _ in groupby(sorted(lst))]

def unique_unsortable(lst):
    #lst = sorted_by_shape(lst)
    lst = lst[:]
    i = 0
    while i < len(lst):
        j = i + 1
        while j < len(lst):
            if lst[i] == lst[j]:
                lst.pop(j)
            else:
                j += 1

        i += 1
    return lst


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
    return lambda *args: sum(f(*args) for f in fs)


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
