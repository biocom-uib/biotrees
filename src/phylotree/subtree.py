from phylotree import PhyloTree, delete_nodes_with_out_degree_one

from shape.newick import from_newick
import combinatorics


def subtree(t, lvs):
    """
    Returns the subtree of t with the given leaves.
    :param t: `PhyloTree` instance.
    :param lvs: `list` instance.
    :return: `PhyloTree` instance.
    """
    children = []
    if t:
        if t.is_leaf and t.leaf in lvs:
            return t
        else:
            for ch in t.children:
                if ch.is_leaf and ch.leaf in lvs:
                    children.append(ch)
                elif not ch.is_leaf:
                    subt = subtree(ch, lvs)
                    if subt:
                        children.append(subt)
            if children:
                return delete_nodes_with_out_degree_one(PhyloTree(None, children))


def all_pairs_of_subtrees_of_m_leaves_that_share_k_leaves(t, m, k):
    """
    Returns a list of tuples with all pairs of subtrees of m leaves that share k leaves.
    :param t: `PhyloTree` instance.
    :param m: `int` instance.
    :param k: `int` instance.
    :return: `list` instance.
    """
    n = t.count_leaves()
    pairs = []
    if k == 0:
        ps = combinatorics.pairs_of_disjoint_subsets_with_k_elements(t.label_set(), m)
        return [(subtree(t, s1), subtree(t, s2)) for s1, s2 in ps]
    else:
        for s in combinatorics.subsets_with_k_elements(range(n), k):
            ps = combinatorics.pairs_of_subsets_with_k_elements_that_share_exactly_subset_s(range(n), m, s)
            forest = [(subtree(t, s1), subtree(t, s2)) for s1, s2 in ps]

            pairs = pairs + forest
        return pairs


def all_pairs_of_subtrees_of_m_leaves_that_share_k_leaves_with_shapes(t, m, k, s1, s2):
    """
    Returns a list of tuples with all pairs of subtrees of m leaves with shape given shapes that share k leaves.
    :param t: `PhyloTree` instance.
    :param m: `int` instance.
    :param k: `int` instance.
    :return: `list` instance.
    """
    p  = all_pairs_of_subtrees_of_m_leaves_that_share_k_leaves(t, m, k)
    P = []
    for pair in p:
        if (pair[0].shape().iso(s1) and pair[1].shape().iso(s2)) or (pair[0].shape().iso(s2) and pair[1].shape().iso(s1)):
            P.append(pair)
    return P


def all_pairs_of_subtrees_of_m_leaves_that_share_k_leaves_with_both_shapes_q3(t, m, k):
    """
    Returns a list of tuples with all pairs of subtrees of m leaves with shape Q3 that share k leaves.
    :param t: `PhyloTree` instance.
    :param m: `int` instance.
    :param k: `int` instance.
    :return: `list` instance.
    """
    q3 = from_newick("((*,*),(*,*));")
    return all_pairs_of_subtrees_of_m_leaves_that_share_k_leaves_with_shapes(t, m, k, q3, q3)
