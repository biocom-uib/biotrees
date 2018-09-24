from random import randint

from biotrees.phylotree import PhyloTree
from biotrees.shape import count_nodes_by_depth


def chain(x, height):
    if height == 0:
        return PhyloTree(x)
    else:
        return PhyloTree(None, [chain(x, height-1)])


def random_adding_chains(nleaves, height, depth_dist=None):
    if depth_dist is None:
        depth_dist = lambda: randint(0, height-1)

    def add_chain_at(x, t, depth, cur_level):
        assert not t.is_leaf()

        if depth == 0:
            # safe, still no sharing at this point
            t.children.append(chain(x, height-cur_level-1))
        else:
            child_index = randint(0, len(t.children)-1)
            add_chain_at(x, t.children[child_index], depth-1, cur_level+1)

    t = chain("0", height)

    for i in range(1, nleaves):
        add_chain_at(str(i), t, depth_dist(), 0)

    t._sort()
    return t

