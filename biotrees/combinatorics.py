from itertools import permutations, combinations


def subset(l, L):
    """
    Takes two lists and returns True if the first one is contained in the second one. If the lists could be sorted,
    it would be more efficient.
    :param l: `list` instance.
    :param L: `list` instance.
    :return: `bool` instance.
    """
    return all(x in L for x in l)


def subsets_with_k_elements_that_contain_subset_s(S, k, s):
    """
    Returns a list with all the sublists of S of length k that contain a given sublist s.
    :param S: `list` instance.
    :param k: `int` instance.
    :param s: `list` instance.
    :return: `list` instance.
    """
    s = tuple(s)

    if k < len(s):
        pass
    elif not subset(s, S):
        pass
    elif k == len(s):
        yield s
    else:
        for x in combinations(set(S) - set(s), k - len(s)):
            yield s + x


def pairs_of_disjoint_subsets_with_k_elements(S, k):  # los elementos son ordenables (por ejemplo, los arboles)
    """
    Returns a list of pairs of sublists with length k that share no elements. This assumes the elements to be
    ordinal.
    :param S: `list` instance.
    :param k: `int` instance.
    :return: `list` instance.
    """
    S = set(S)

    if len(S) >= 2*k:
        for s1 in combinations(S, k):
            for s2 in combinations(S - set(s1), k):
                if s1 <= s2:
                    yield s1, s2


def pairs_of_subsets_with_k_elements_that_share_exactly_subset_s(S, k, s):
    """
    Returns a list of pairs of sublists with length k that share a given subset s. This assumes the elements to be
    ordinal.
    :param S: `list` instance.
    :param k: `int` instance.
    :param s: `list` instance.
    :return: `list` instance.
    """
    s = tuple(s)
    if k >= len(s):
        for s1, s2 in pairs_of_disjoint_subsets_with_k_elements(set(S) - set(s), k - len(s)):
            yield s1+s, s2+s


def finite_bijections(A, B):
    """
    Gives a list with all finite bijections between lists A and B. A bijection is expressed as a list of pairs, whose
    first component is in A and its second component is the image of the first through the bijection, in B.
    :param A: `list` instance.
    :param B: `list` instance.
    :return: `list` instance.
    """
    A = tuple(A)
    B = tuple(B)
    if len(A) != len(B):
        pass
    else:
        for perm in permutations(B):
            yield {A[i]: b for i, b in enumerate(perm)}
