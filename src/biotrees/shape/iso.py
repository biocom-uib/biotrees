def equal(t1, t2):
    """
    Since our `Shape` objects are sorted, to know whether two trees are equal or not it suffices to know if
    they are equal or not.
    :param t1, t2: the `Shape` objects we want to compare.
    :return: `bool` instance.
    """
    t1._sort()  # esto a la larga tendrá que desaparecer
    t2._sort()
    assert t1._is_sorted()
    assert t2._is_sorted()
    return t1 == t2


def isomorphic(t1, t2):
    # arreglar
    """
    Since our `Shape` objects are sorted, to know whether two trees are isomorphic or not it suffices to know if
    they are equal or not.
    :param t1, t2: the `Shape` objects we want to compare.
    :return: `bool` instance.
    """
    return equal(t1, t2)
