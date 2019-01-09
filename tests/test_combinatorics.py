import unittest

from math import factorial
from itertools import combinations
import biotrees.combinatorics as combinatorics

SETS = [set(), set(range(10)), {'a','b','c','d','e'}]
SUBSETS = [[set(s) for k in range(len(S)+1) for s in combinations(S, k)] for S in SETS]

LISTS = [list(S) for S in SETS]
SUBLISTS = [[list(s) for s in S] for S in SUBSETS]

TUPLES = [tuple(S) for S in SETS]
SUBTUPLES = [[tuple(s) for s in S] for S in SUBSETS]

SAMPLES = SETS + LISTS + TUPLES
SUBSAMPLES = SUBSETS + SUBLISTS + SUBTUPLES


class TestCombinatorics(unittest.TestCase):

    def test_subset(self):
        for i in range(len(SAMPLES)):
            S = SAMPLES[i]
            for s in SUBSAMPLES[i]:
                self.assertTrue(
                    combinatorics.subset(s, S))
                self.assertEqual(
                    combinatorics.subset(S, s), s == S)

        self.assertFalse(
            combinatorics.subset([1], {}))


    def test_subsets_with_k_elements_that_contain_subset_s(self):
        for i in range(len(SAMPLES)):
            S = SAMPLES[i]
            for k in range(len(S)+1):
                for s in SUBSAMPLES[i]:
                    r = list(combinatorics.subsets_with_k_elements_that_contain_subset_s(S, k, s))

                    if k < len(s):
                        self.assertEqual(r, [])

                    elif k == len(s):
                        self.assertEqual(r, [tuple(s)])

                    else:
                        self.assertEqual(
                            len(set(r)),
                            factorial(len(S)-len(s)) / factorial(len(S)-k) / factorial(k - len(s)))

                    self.assertTrue(
                        all(combinatorics.subset(s, sr) for sr in r))

                    self.assertTrue(
                        all(len(sr) == k for sr in r))


    def test_pairs_of_disjoint_subsets_with_k_elements(self):
        for i in range(len(SAMPLES)):
            S = SAMPLES[i]

            for k in range(len(S)+1):
                r = list(combinatorics.pairs_of_disjoint_subsets_with_k_elements(S, k))

                if k == 0:
                    self.assertEqual(r, [((), ())])

                elif len(S) < 2*k:
                    self.assertEqual(r, [])

                elif len(S) == 2*k:
                    self.assertEqual(
                        len(set(r)),
                        factorial(len(S)) / factorial(k) / factorial(len(S)-k)
                        *
                        factorial(len(S)-k) / factorial(k) / factorial(len(S)-2*k)
                        / 2)

                self.assertTrue(
                    all(combinatorics.subset(sr, S) for pair in r for sr in pair))

                self.assertTrue(
                    all(len(sr) == k for pair in r for sr in pair))


    def test_pairs_of_subsets_with_k_elements_that_share_exactly_subset_s(self):
        for i in range(len(SAMPLES)):
            S = SAMPLES[i]

            for k in range(len(S)+1):
                for s in SUBSAMPLES[i]:
                    r = list(combinatorics.pairs_of_subsets_with_k_elements_that_share_exactly_subset_s(S, k, s))

                    self.assertTrue(
                        all(set(sr1).intersection(set(sr2)) == set(s) for sr1, sr2 in r))

                    self.assertTrue(
                        all(len(sr) == k for pair in r for sr in pair))


    def test_finite_bijections(self):
        for A in SAMPLES:
            if len(A) > 8:
                break

            for B in SAMPLES:
                r = list(combinatorics.finite_bijections(A, B))

                for bij in r:
                    self.assertEqual(set(bij.keys())  , set(A))
                    self.assertEqual(set(bij.values()), set(B))

                if len(A) == len(B):
                    self.assertEqual(len(r), factorial(len(A)))
                else:
                    self.assertEqual(r, [])

        self.assertEqual(
            list(combinatorics.finite_bijections({1, 2}, {'a', 'b'})),
            [ {1: 'a', 2: 'b'},
              {1: 'b', 2: 'a'} ])
