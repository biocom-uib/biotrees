from biotrees.phylotree import PhyloTree
import newick as _newick


def to_newick(phylo):
    """
    Returns a string representing the simplified Newick code of the input.
    :param: `PhyloTree` instance.
    :return: `str` instance.
    """
    return phylo_to_newick_node(phylo).newick


def from_newick(nwk):
    """
    Create a `PhyloTree` object from a Newick code entered as a string.
    :param nwk: a string representing a Newick code.
    :return: `PhyloTree` instance.
    """
    return newick_node_to_phylo(_newick.loads(nwk)[0])


def from_newick_list(nwk):
    """
    Create a list of `PhyloTree` objects from a list of Newick codes entered as a string.
    :param nwk: a string representing a list of Newick codes.
    :return: `list` instance.
    """
    return [newick_node_to_phylo(n) for n in _newick.loads(nwk)]


def newick_node_to_phylo(node):
    """
    Create a `PhyloTree` object from a `Node` object.
    :param node: `Node` instance.
    :return: `PhyloTree` instance.
    """
    if not bool(node.descendants):
        return PhyloTree(node.name, None)
    return PhyloTree(None, sorted([newick_node_to_phylo(x) for x in node.descendants]))


def phylo_to_newick_node(phylo):
    """
    Creates a `Node` instance with same shape as the shape of the input.
    :param phylo: `PhyloTree` instance.
    :return: `Node` instance.
    """
    if phylo.is_leaf():
        return _newick.Node.create(str(phylo.leaf))
    return _newick.Node.create(descendants=[phylo_to_newick_node(child) for child in phylo.children])


def trees_from_file(fname, encoding='utf8', strip_comments=False, **kw):
    """
    Load a list of trees from a Newick formatted file.
    :param fname: file path.
    :param strip_comments: Flag signaling whether to strip comments enclosed in square \
    brackets.
    :param kw: Keyword arguments are passed through to `Node.read`.
    :return: `list` instance.
    """
    l = _newick.read(fname, encoding, strip_comments, **kw)
    return [newick_node_to_phylo(n) for n in l]
