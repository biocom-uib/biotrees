"""
In order to read Newick codes we import the newick module from https://github.com/glottobank/python-newick.
"""

from shape import Shape

import _newick


def to_newick(shape):
    """
    Returns a string representing the simplified Newick code of the input.
    :param: `Shape` instance.
    :return: `str` instance.
    """
    return shape_to_newick_node(shape).newick


def from_newick(nwk):
    """
    Create a `Shape` object from a Newick code entered as a string.
    :param nwk: a string representing a Newick code.
    :return: `Shape` instance.
    """
    return newick_node_to_shape(_newick.loads(nwk)[0])


def from_newick_list(nwk):
    """
    Create a list of `Shape` objects from a list of Newick codes entered as a string.
    :param nwk: a string representing a list of Newick codes.
    :return: `list` instance.
    """
    return [newick_node_to_shape(n) for n in _newick.loads(nwk)]


def newick_node_to_shape(node):
    """
    Create a `Shape` object from a `Node` object.
    :param node: a `Node`.
    :return: `Shape` instance.
    """
    if not bool(node.descendants):
        return Shape(None)
    return Shape(sorted([newick_node_to_shape(x) for x in node.descendants]))


def shape_to_newick_node(shape):
    """
    Creates a `Node` instance with same shape as the input one.
    :param shape: `Shape` instance.
    :return: `Node` instance.
    """
    if shape.is_leaf:
        return _newick.Node.create("*")
    return _newick.Node.create(descendants = [shape_to_newick_node(child) for child in shape.children])


def shapes_from_file(fname, encoding='utf8', strip_comments=False, **kw):
    """
    Load a list of shapes from a Newick formatted file.
    :param fname: file path.
    :param strip_comments: Flag signaling whether to strip comments enclosed in square \
    brackets.
    :param kw: Keyword arguments are passed through to `Node.read`.
    :return: `list` instance.
    """
    l = _newick.read(fname, encoding, strip_comments, **kw)
    return [newick_node_to_shape(n) for n in l]