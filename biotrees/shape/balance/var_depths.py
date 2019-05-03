from statistics import variance
from math import log, floor

from biotrees.shape import get_leaf_depths


def var_depths(t):
    return variance(get_leaf_depths(t))


def min_var_depths(n):
    k = int(floor(log(n, 2)))
    r = n - 2**k

    if r == 0:
        d = k
        m = [0 for _ in range(d)]
        m[0] = n
        return 0, m
    else:
        m1 = min_var_delta(n, k+1)
        m2 = min_var_delta(n, k+2)

        if m1[0] == m2[0]:
            return m1, m2
        elif m1[0] > m2[0]:
            return m2
        else:
            return m1


def min_var_delta(n, d):    # suponemos que d > 3

    def recurse(k, j):

        k1 = k[:]
        k2 = k[:]

        k1[j] = 2 * k1[j + 1] - 1
        k2[j] = 2 * k2[j + 1]

        if j > 0:
            ks = recurse(k1, j - 1)
            ks2 = recurse(k2, j - 1)

            for kii in ks2:
                ks.append(kii)

            return ks

        elif j == 0:
            return [k1, k2]

    k = [0 for _ in range(d-2)]
    k.append(1)

    vms2 = [(nvar(solve(n, ki)), solve(n, ki)) for ki in recurse(k, d-3)
            if solve(n, ki)[0] > 0 and solve(n, ki)[1] >= 0]

    return min(vms2)


def solve(n, k):

    d = len(k)

    m = [2*k[j] - k[j-1] for j in range(1, d)]
    ms = sum(mi for mi in m)
    m1 = 4*k[0] - n + ms  # revisar, como todo lo demas
    m0 = n - m1 - ms
    return [m0, m1] + m


def mean(m):
    n = sum(mi for mi in m)
    return sum(i*m[i] for i in range(len(m))) / n


def nvar(m):
    mm = mean(m)
    return sum(m[i]*(i - mm)**2 for i in range(len(m)))


def min_component(m, j):
    min = m[0]
    for i in range(1,len(m)):
        if m[i][j] < min[j]:
            min = m[i]
    return min


def feasiblity_proof(m): # las primeras entradas son las profundidades más altas
    ns = [0 for _ in range(len(m)-1)]
    ns.append(2)
    for i in reversed(range(0, len(m)-1)):
        ns[i] = 2*(ns[i+1] - m[i+1])
    return m[0] == ns[0] and all(mi >= 0 for mi in m)


def to_file(vars):
    f = open("vars", "w")

    for v in vars:
        f.write(str(v) + "\n")
    f.close()


def test(n):
    vars = [min_var_depths(i) for i in range(4, n + 1)]
    to_file(vars)


def arnau(n1, n2):
    for i in range(n1, n2+1):
        print(min_var_depths(i))


def max_bal_nvar(n):
    m = int(floor(log(n, 2)))
    twotom = 2 ** m
    k = n - twotom

    return 2*k*(twotom - k) / n


def cesc(n):
    m = int(floor(log(n, 2)))
    twotom = 2 ** m
    k = n - twotom

    ts = []

    for l in range(5, m+1):
        twotol = 2**l
        t = [0 for _ in range(m+1)]
        t[0] = twotol + 2 * (k-1)
        t[1] = twotom - twotol - k + 1
        t[l] = 1

        if nvar(t) < 2*k*(twotom - k) / n and feasiblity_proof(t):
            ts.append(l)

    return ts


def to_files(k1, k2):
    for k in range(k1, k2+1):
        f = open("vars " + str(k) + "-" + str(k+1), "w")
        for n in range(2**k, 2**(k+1)):
            f.write(str(min_var_depths(n)) + "\n")
        f.close()


def to_files2(k1, k2):
    for k in range(k1, k2+1):
        f = open("mins " + str(k) + "-" + str(k+1), "w")
        for n in range(2**k, 2**(k+1)):
            f.write(str(n) + " --> " + str(cesc(n)) + "\n")
        f.close()


def ratio_files(k):
    f = open("mins " + str(k) + "-" + str(k+1), "r")

    c = 0
    for l in f.readlines():
        if l[-3] == "[" and l[-2] == "]":
            c += 1

    f.close()

    return c / 2**k

