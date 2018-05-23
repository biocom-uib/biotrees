def equal(t1, t2):
    """
    Since our `PhyloTree` objects are sorted, to know whether two trees are equal or not it suffices to know if
    they are equal or not.
    :param t1, t2: the `PhyloTree` objects we want to compare.
    :return: `bool` instance.
    """
    assert t1._is_sorted()
    assert t2._is_sorted()
    return t1 == t2


def isomorphic(t1, t2):
    """
    Since our `PhyloTree` objects are sorted, to know whether two trees are isomorphic or not it suffices to know if
    they are equal or not.
    :param t1, t2: the `PhyloTree` objects we want to compare.
    :return: `bool` instance.
    """
    return equal(t1, t2)
