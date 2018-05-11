def subsets_with_k_elements(S, k):
    """
    Gives a list with all the sublists of S with k elements.
    :param S: `list` instance.
    :param k: `int` instance.
    :return: `list` instance.
    """
    if k == 0:
        return
    elif k == 1:
        for x in S:
            yield [x]
    else:
        for i in range(0, len(S)):
            s = S[i + 1:]
            for x in subsets_with_k_elements(s, k - 1):
                yield [S[i]] + x


def set_minus(L,l):
    """
    Takes two lists and substracts all the elements of the second list to the first one.
    :param L: `list` instance.
    :param l: `list` instance.
    :return: `list` instance.
    """
    return [x for x in L if x not in l]


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
    if k < len(s):
        return
    elif not subset(s, S):
        return
    elif k == len(s):
        yield s
    else:
        for x in subsets_with_k_elements(set_minus(S, s), k - len(s)):
            yield s + x


def pairs_of_disjoint_subsets_with_k_elements(S, k):  # los elementos son ordenables (por ejemplo, los arboles)
    """
    Returns a list of pairs of sublists with length k that share no elements. This assumes the elements to be
    ordinal.
    :param S: `list` instance.
    :param k: `int` instance.
    :return: `list` instance.
    """
    if len(S) < 2*k:
        return []
    else:
        pairs = []
        for s1 in subsets_with_k_elements(S, k):
            Sminuss1 = set_minus(S, s1)
            for s2 in subsets_with_k_elements(Sminuss1, k):
                if s1[0] < s2[0]:
                    pairs.append((s1, s2))
                elif s2[0] < s1[0]:
                    pairs.append((s2, s1))
        pairs.sort()
        return pairs


def pairs_of_subsets_with_k_elements_that_share_exactly_subset_s(S, k, s):
    """
    Returns a list of pairs of sublists with length k that share a given subset s. This assumes the elements to be
    ordinal.
    :param S: `list` instance.
    :param k: `int` instance.
    :param s: `list` instance.
    :return: `list` instance.
    """
    pairs = pairs_of_disjoint_subsets_with_k_elements(set_minus(S, s), k - len(s))
    if pairs:
        return [(s1 + s, s2 + s) for s1, s2 in pairs]
    else:
        return [(s, s)]


def finite_bijections(A, B):
    """
    Gives a list with all finite bijections between lists A and B. A bijection is expressed as a list of pairs, whose
    first component is in A and its second component is the image of the first through the bijection, in B.
    :param A: `list` instance.
    :param B: `list` instance.
    :return: `list` instance.
    """
    if len(A) != len(B):
        return []
    elif not A:
        return []
    elif len(A) == 1:
        return [[(A[0], B[0])]]
    else:
        return [[(A[0], b)] + bij for b in B
                for bij in finite_bijections(A[1:], set_minus(B, b))]


def permutations(S):
    """
    Gives a list with all bijections from S to itself.
    :param S: `list` instance.
    :return: `list` instance.
    """
    return finite_bijections(S, S)
