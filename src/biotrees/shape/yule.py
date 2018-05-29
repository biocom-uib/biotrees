
from biotrees.phylotree import shape_to_phylotree
from biotrees.util import and_then, parametric_total_probabilities, groupby_sorted

import biotrees.phylotree.yule as phylo_yule


@and_then(parametric_total_probabilities(grouper=groupby_sorted))
def yule_from_t(sh, prob):
    """
    Returns a list of tuples containing each shape obtained from `Shape` sh (assuming sh has probability prob) with
    their associate probability under the Yule model.
    :param sh: `Shape` instance.
    :param prob: `function` instance.
    :return: `list` instance.
    """
    t = shape_to_phylotree(sh)
    for ti, p in phylo_yule.yule_from_t(t, prob):
        yield ti.shape(), p


@and_then(parametric_total_probabilities(grouper=groupby_sorted))
def yule(n):
    """
    Returns a list of tuples containing all the shapes of n leaves that can be obtained under the Yule model and their
    corresponding probability.
    :param n: `int` instance.
    :return: `list` instance.
    """
    for ti, p in phylo_yule.pseudo_yule(n):
        yield ti.shape(), p
