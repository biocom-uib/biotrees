import numpy as np

from random import randint

from biotrees.shape import Shape, count_leaves, count_nodes_by_depth, get_depth
import biotrees.phylotree.taxo.generator as phylo_taxo_gen


def chain(height):
    if height == 0:
        return Shape.LEAF
    else:
        return Shape([chain(height-1)])


def random_adding_chains(nleaves, height, depth_dist=None):
    return phylo_taxo_gen.random_adding_chains(nleaves, height, depth_dist).shape()


def estimate_depth_dist(taxo):
    nodes_by_depth = count_nodes_by_depth(taxo)

    depth = len(nodes_by_depth) - 1
    n = nodes_by_depth[-1]

    return [(nodes_by_depth[i+1] - nodes_by_depth[i]) / (n-1) for i in range(0, depth)]


def tree_median_with(taxo, measure, nsample):
    nleaves = count_leaves(taxo)
    height = get_depth(taxo)

    depth_dist_probs = estimate_depth_dist(taxo)
    depth_dist = lambda: np.random.choice(np.arange(len(depth_dist_probs)), p=depth_dist_probs)

    samples = np.fromfunction(
            lambda i: measure(random_adding_chains(nleaves, height, depth_dist)),
            (nsample,))

    return np.median(samples)
