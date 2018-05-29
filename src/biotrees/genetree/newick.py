from biotrees.genetree import GeneTree
import biotrees._newick as _newick


def to_newick(genet):
    return genetree_to_newick_node(genet).newick


def from_newick(nwk):
    """
    Create a `phylotree` object from a Newick code entered as a string.
    :param nwk: a string representing a Newick code.
    :return: `phylotree` instance.
    """
    return newick_node_to_genetree(_newick.loads(nwk)[0])


def from_newick_list(nwk):
    """
    Create a list of `phylotree` objects from a list of Newick codes entered as a string.
    :param nwk: a string representing a list of Newick codes.
    :return: [`phylotree`] instance.
    """
    return[newick_node_to_genetree(n) for n in _newick.loads(nwk)]


def newick_node_to_genetree(node):
    """
    #Create a `genetree` object from a `Node` object.
    #:param N: a Node.
    #:return: `genetree` instance.
    """
    i = [-1]
    def recurse(node):
        if not bool(node.descendants):
            i[0] += 1
            return GeneTree(str(i[0]+1), node.name, None)
        else:
            return GeneTree(None, None, sorted([recurse(ch) for ch in node.descendants]))

    return recurse(node)


def genetree_to_newick_node(genet):
    #problemas
    if genet.is_leaf():
        return _newick.Node.create(str(genet.label))
    return _newick.Node.create(descendants=[genetree_to_newick_node(child) for child in genet.children])


def trees_from_file(fname, encoding='utf8', strip_comments=False, **kw):
    """
    Load a list of trees from a Newick formatted file.
    :param fname: file path.
    :param strip_comments: Flag signaling whether to strip comments enclosed in square \
    brackets.
    :param kw: Keyword arguments are passed through to `Node.read`.
    :return: [`genetree`] instance.
    """
    l = _newick.read(fname, encoding, strip_comments, **kw)
    return [newick_node_to_genetree(n) for n in l]

#NEWICK ESTÁ ROTO