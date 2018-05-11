import shape.generator
import shape.newick
import phylotree

import random

cherry = shape.newick.from_newick("(*, *);")


def duplicate_leaf(t, l, newl):
    # notice that here t is a phylogenetic tree; elsewhere is just a shape
    if t.is_leaf and t.leaf == l:
        return phylotree.PhyloTree(None, [t, phylotree.PhyloTree(newl, None)])
    elif t.is_leaf and t.leaf != l:
        return t
    else:
        phylot = phylotree.PhyloTree(None, [duplicate_leaf(ch, l, newl) for ch in t.children])
        return phylot


def yule_from_tree(t):
    n = t.count_leaves()
    if n < 1:
        return shape.newick.from_newick("*;")
    elif n == 1:
        return shape.newick.from_newick("(*, *);")
    else:
        phylot = phylotree.shape_to_phylotree(t)    # I give my leaves names in order to be able
                                                    # to identify them; phylot is a phylogenetic
                                                    # tree
        l = random.randrange(n)
        result = duplicate_leaf(phylot, l, n).shape()
        result.sort()
        return result


def yule(n):
    if n < 1:
        return
    elif n == 1:
        return shape.newick.from_newick("*;")
    elif n == 2:
        return cherry
    elif n > 2:
        # print n
        return yule_from_tree(yule(n - 1))


def trees_from_yule(ntrees, nlabels):
    return [yule(nlabels) for _ in range(ntrees)]