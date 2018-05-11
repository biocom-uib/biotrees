from phylotree import shape_to_phylotree
import phylotree.yule


def yule_from_t(sh, prob):
    """
    Returns a list of tuples containing each shape obtained from `Shape` sh (assuming sh has probability prob) with
    their associate probability under the Yule model.
    :param sh: `Shape` instance.
    :param prob: `function` instance.
    :return: `list` instance.
    """
    t = shape_to_phylotree(sh)
    tps = [(ti.shape(), p) for ti, p in phylotree.yule.yule_from_t(t, prob)]
    return phylotree.yule.collapse_tree_prob_list(tps)


def yule(n):
    """
    Returns a list of tuples containing all the shapes of n leaves that can be obtained under the Yule model and their
    corresponding probability.
    :param n: `int` instance.
    :return: `list` instance.
    """
    tps = [(ti.shape(), p) for ti, p in phylotree.yule.pseudo_yule(n)]
    return phylotree.yule.collapse_tree_prob_list(tps)
