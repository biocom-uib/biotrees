from sympy import simplify
from functools import lru_cache

from biotrees.util import lifted_prod


def liu_polynomial(t):
    if t.is_leaf():
        return lambda y, x : x
    else:
        return lambda y, x : y + simplify(lifted_prod([liu_polynomial(ch) for ch in t.children], y, x))


def liu_polynomial_comb(n):
    return lambda y, x : x**n + sum(x**i * y for i in range(n-1))


@lru_cache(maxsize=None)
def liu_polynomial_max_balanced(n):
    if n <= 0:
        return 0
    elif n == 1:
        return lambda y, x : x
    elif n == 2:
        return lambda y, x : x**2 + y
    else:
        s = n % 2
        return lambda y, x : y + simplify(liu_polynomial_max_balanced((n-s)/2)(y, x)
                                          * liu_polynomial_max_balanced((n-s)/2 + s)(y, x))


def number_monomials(p):
    return p(1, 1)


def balance(t):
    return liu_polynomial(t)(1, 1)
