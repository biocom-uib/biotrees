from functools import lru_cache

from biotrees.phylotree import shape_to_phylotree, phylotree_to_shape
import biotrees.phylotree.yule as phylo_yule
from biotrees.util import and_then, parametric_total_probabilities



def sim_yule_from_t(t):
    return phylotree_to_shape(phylo_yule.sim_yule_from_t(shape_to_phylotree(t)))


def sim_yule(n):
    return phylotree_to_shape(phylo_yule.sim_yule(n))


@and_then(parametric_total_probabilities)
def yule_from_t(sh, prob):
    """
    Returns a list of tuples containing each shape obtained from `Shape` sh (assuming sh has probability prob) with
    their associate probability under the Yule model.
    :param sh: `Shape` instance.
    :param prob: `function` instance.
    :return: `list` instance.
    """
    for t, p in phylo_yule.yule_from_t(shape_to_phylotree(sh), prob):
        yield phylotree_to_shape(t), p


@lru_cache(maxsize=-1)
@and_then(parametric_total_probabilities)
def yule(n):
    """
    Returns a list of tuples containing all the shapes of n leaves that can be obtained under the Yule model and their
    corresponding probability.
    :param n: `int` instance.
    :return: `list` instance.
    """
    for t, p in phylo_yule.yule(n):
        yield phylotree_to_shape(t), p
