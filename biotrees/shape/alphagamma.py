from functools import lru_cache
from sympy import simplify

from biotrees.shape import Shape, count_leaves
from biotrees.shape.generator import add_leaf_to_edge, add_leaf_to_node, iter_replace_tree_at
from biotrees.util import and_then, parametric_total_probabilities


import sys
sys.setrecursionlimit(2000)


@and_then(parametric_total_probabilities)
def alphagamma_from_t(t, prob):
    """
    Returns a list of tuples containing each shape obtained from `Shape` sh (assuming sh has probability prob) with
    their associate probability under the Alpha-Gamma model.
    :param sh: `Shape` instance.
    :param prob: `function` instance.
    :return: `list` instance.
    """
    n = count_leaves(t)
    if t.is_leaf():
        return [(add_leaf_to_edge(t), lambda a, c: simplify(prob(a, c) * (1 - a) / (n - a)))]
    else:
        for i in range(len(t.children)):
            k = count_leaves(t.children[i])
            if n > 1:
                p_times_m = lambda a, c, k=k: simplify(prob(a, c) * (k - a) / (n - a))
            else:
                p_times_m = prob
            for ti, q in alphagamma_from_t(t.children[i], p_times_m):
                ts2 = list(iter_replace_tree_at(t.children, i, ti))
                yield (Shape(ts2), q)

        yield (add_leaf_to_edge(t), lambda a, c: simplify(prob(a, c) * c / (n - a)))

        k = len(t.children)
        if k > 1:
            prob2 = lambda a, c, k=k: simplify(prob(a, c) * ((k - 1) * a - c) / (n - a))
        else:
            prob2 = prob

        yield (add_leaf_to_node(t), prob2)


@lru_cache(maxsize=None)
@and_then(parametric_total_probabilities)
def alphagamma(n):
    """
    Returns a list of tuples containing all the shapes of n leaves that can be obtained under the Alpha-Gamma model and
    their corresponding probability.
    :param n: `int` instance.
    :return: `list` instance.
    """
    if n <= 0:
        pass
    elif n == 1:
        yield (Shape.LEAF, lambda a, c: 1)
    elif n == 2:
        yield (Shape([Shape.LEAF, Shape.LEAF]), lambda a, c: 1)
    else:
        for t, prob in alphagamma(n-1):
            yield from alphagamma_from_t(t, prob)
