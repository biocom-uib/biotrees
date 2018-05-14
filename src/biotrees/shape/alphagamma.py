from sys import setrecursionlimit
setrecursionlimit(2000)

from sympy import simplify

from biotrees.shape import Shape, sorted_tree
from biotrees.shape.generator import add_leaf_to_edge, add_leaf_to_node, collapse_tree_prob_list
from biotrees.shape.iso import iso



def alphagamma_from_t(t, prob):
    """
    Returns a list of tuples containing each shape obtained from `Shape` sh (assuming sh has probability prob) with
    their associate probability under the Alpha-Gamma model.
    :param sh: `Shape` instance.
    :param prob: `function` instance.
    :return: `list` instance.
    """
    n = t.count_leaves()
    if t.is_leaf:
        return [(add_leaf_to_edge(t), lambda a, c: simplify(prob(a, c) * (1 - a) / (n - a)))]
    else:
        tps = []
        for i in range(len(t.children)):
            k = t.children[i].count_leaves()
            if n > 1:
                p_times_m = lambda a, c, k = k: simplify(prob(a, c) * (k - a) / (n - a))
            else:
                p_times_m = prob
            for ti, q in alphagamma_from_t(t.children[i], p_times_m):
                ts = t.children[:]
                ts[i] = ti
                tps.append((sorted_tree(Shape(ts)), q))

        tps.append((sorted_tree(add_leaf_to_edge(t)), lambda a, c: simplify(prob(a, c) * c / (n - a))))
        k = len(t.children)
        if k > 1:
            prob2 = lambda a, c, k = k: simplify(prob(a, c) * ((k - 1) * a - c) / (n - a))
        else:
            prob2 = prob

        tps.append((sorted_tree(add_leaf_to_node(t)), prob2))

        return collapse_tree_prob_list(tps, iso)


def alphagamma(n):
    """
    Returns a list of tuples containing all the shapes of n leaves that can be obtained under the Alpha-Gamma model and
    their corresponding probability.
    :param n: `int` instance.
    :return: `list` instance.
    """
    if n <= 0:
        return []
    elif n == 1:
        return [(Shape(), lambda a, c: 1)]
    elif n == 2:
        return [(Shape([Shape(), Shape()]), lambda a, c: 1)]
    else:
        tps = []
        prev = alphagamma(n - 1)

        for t, prob in prev:
            tps.extend(alphagamma_from_t(t, prob))

        return collapse_tree_prob_list(tps, iso)
