import unittest

from math import factorial
import biotrees.util as util


class TestUtil(unittest.TestCase):

    def test_binom2(self):
        self.assertEqual(
            util.binom2(5),
            factorial(5) // (factorial(5-2) * factorial(2)))

        self.assertEqual(
            util.binom2(10),
            factorial(10) // (factorial(10-2) * factorial(2)))

        self.assertTrue(
            isinstance(util.binom2(3), int))

        self.assertTrue(
            isinstance(util.binom2(0), int))

        self.assertEqual(
            util.binom2(0),
            util.binom2(1),
            0)


    def test_skip_nth(self):
        self.assertEqual(
            list(util.skip_nth(range(1,10), 3)),
            [1,2,3,5,6,7,8,9])

        self.assertEqual(
            list(util.skip_nth(range(1,10), 20)),
            list(range(1,10)))

        self.assertEqual(
            list(util.skip_nth([], 0)),
            [])


    def test_unique(self):
        self.assertEqual(
            sorted(util.unique([1,2,2,4,2,2,3,4,3,2,1,2,4,1,3,2,1,4,1,2,4,1,2,1,1,1])),
            [1,2,3,4])

        self.assertEqual(
            util.unique([1,2,2,4,2,2,3,4,3,2,1,2,4,1,3,2,1,4,1,2,4,1,2,1,1,1], sort=True),
            [1,2,3,4])

        self.assertEqual(
            sorted(util.unique([[], [1,2,3], [], [1,2], [1,2,3]])),
            sorted([[], [1,2], [1,2,3]]))

        self.assertEqual(
            util.unique([]),
            [])

        self.assertEqual(
            util.unique([], sort=True),
            [])


    def test_unique_unsortable(self):
        self.assertEqual(
            util.unique_unsortable([1,2,2,4,2,2,3,4,3,2,1,2,4,1,3,2,1,4,1,2,4,1,2,1,1,1]),
            [1,2,4,3])

        self.assertEqual(
            util.unique_unsortable([[], {}, set(), {}, [], set(), [1], {'a': {}}]),
            [[], {}, set(), [1], {'a': {}}])


    def test_iter_merge(self):
        self.assertEqual(
            list(util.iter_merge(
                [1,3,5,7,9],
                [0,2,4,6,8])),
            list(range(10)))

        self.assertEqual(
            list(util.iter_merge(range(10), {})),
            list(range(10)))

        self.assertEqual(
            list(util.iter_merge([], range(10))),
            list(range(10)))

        self.assertEqual(
            list(util.iter_merge([1,4,6], [1,5,3,6])),
            [1,1,4,5,3,6,6])

    def test_lifted_sum(self):
        def f(*args, **kwargs):
            return len(args)

        def g(*args, **kwargs):
            return len(kwargs)

        test_args = (1,2,3,4)
        test_kwargs = {'d': 5}

        self.assertEqual(
            util.lifted_sum((f, g))(*test_args, **test_kwargs),
            len(test_args) + len(test_kwargs))

        self.assertEqual(
            util.lifted_sum((f, g))(),
            0)

        self.assertEqual(
            util.lifted_sum([lambda x: x]*2)(5),
            5*2)


    def test_parametric_total_probabilities(self):
        def eval_ps(tps, *args, **kwargs):
            return [(t, p(*args, **kwargs)) for t, p in tps]

        self.assertEqual(
            eval_ps(
                util.parametric_total_probabilities(
                    [('a', lambda: 1), ('b', lambda: 1), ('a', lambda: 2)])),
            [('a', 3), ('b', 1)])

        self.assertEqual(
            eval_ps(
                util.parametric_total_probabilities(
                    [('a', lambda: 1), ('b', lambda: 1), ('a', lambda: 2)])),
            [('a', 3), ('b', 1)])

        p=1
        r = util.parametric_total_probabilities(
                [('a', lambda a, b, p=p: a+b+p), ('b', lambda a, b, p=p: a*b*p), ('a', lambda a, b: a-b)])
        p=2

        self.assertEqual(
            eval_ps(r, 2, b=3),
            [('a', 2+3+1 + 2-3), ('b', 2*3*1)])


    def test_and_then(self):
        def times2(x):
            return x*2

        @util.and_then(times2)
        @util.and_then(lambda x: x+1)
        def f(x):
            return x

        self.assertEqual(
            f(1),
            (1+1)*2)

        self.assertEqual(
            util.and_then(times2)(f)(2),
            (2+1)*2*2)
