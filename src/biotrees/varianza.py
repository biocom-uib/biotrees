from scipy.special import binom

from shape.alphagamma import alphagamma
from phylotree import shape_to_phylotree
from phylotree.subtree import all_pairs_of_subtrees_of_m_leaves_that_share_k_leaves_with_shapes
from phylotree.newick import from_newick

q0 = from_newick("(0,(1,(2,3)));").shape()
q1 = from_newick("(0,1,(2,3));").shape()
q2 = from_newick("(0,(1,2,3));").shape()
q3 = from_newick("((0,1),(2,3));").shape()
q4 = from_newick("(0,1,2,3);").shape()
Q = [q0, q1, q2, q3, q4]

qaux = alphagamma(4)
q = [qaux[0], qaux[3], qaux[1], qaux[2], qaux[4]]


def theta(t, i, j):
    t = shape_to_phylotree(t)
    k = 8 - t.count_leaves()
    return len(all_pairs_of_subtrees_of_m_leaves_that_share_k_leaves_with_shapes(t, 4, k, Q[i], Q[j]))


def first_sum(n, a, c):
    bin = binom(n, 4)
    return bin * sum(i**2 * q[i][1](a, c) for i in range(len(q)))


def second_sum(n, a, c):
    bin = binom(n, 4)
    return bin**2 * sum(i * q[i][1](a, c) for i in range(len(q))) ** 2


def third_sum(n, a, c):
    return sum(i * j * binom(n, k) * theta(t, i, j) * p(a, c)
               for i in range(1, 5)
               for j in range(1, 5)
               for k in range(5, 9)
               for t, p in alphagamma(k))

def Var(n, a, c):
    # if thetas is None:
    #    thetas = compute_thetas(n)
    return first_sum(n, a, c) - second_sum(n, a, c) + third_sum(n, a, c)
# def Var(n, thetas = None)



def leading_coefficient(a,c,q):
    def first_sum(a, c, q):
        return sum(q[i] * q[j] * theta(t, i, j) * p(a, c)
                   for i in range(1, 5)
                   for j in range(1, 5)
                   for t, p in alphagamma(8))
    def second_sum(a, c, q):
        return (sum(q[i] * p(a, c)
                                      for i in range(1, 5)
                                      for t, p in alphagamma(4)))
    return (first_sum(a,c,q), second_sum(a,c,q))
    #return 1/factorial(8) * first_sum(a,c,q) - 1/factorial(4))**2 * second_sum(a,c,q)**2

