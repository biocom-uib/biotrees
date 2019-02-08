import unittest

from biotrees.shape import Shape
from biotrees.phylotree import PhyloTree

L1 = PhyloTree('1')
L2 = PhyloTree('2')
L3 = PhyloTree('3')

C12 = PhyloTree(None, [L1,L2])
C23 = PhyloTree(None, [L2,L3])
C13 = PhyloTree(None, [L1,L3])

C21 = PhyloTree(None, [L2,L1])
C32 = PhyloTree(None, [L3,L2])
C31 = PhyloTree(None, [L3,L1])


class TestPhyloTree(unittest.TestCase):

    def test_init(self):
        self.assertTrue(L1.is_leaf())
        self.assertEqual(L1.leaf, '1')

        self.assertFalse(C12.is_leaf())
        self.assertEqual(len(C12.children), 2)

    def test_clone(self):
        self.assertEqual(L1, L1.clone())

        L1.clone().leaf = '2'
        self.assertEqual(L1.leaf, '1')

        self.assertEqual(C12, C12.clone())

        C12.clone().children.append(L3)
        self.assertEqual(len(C12.children), 2)

    def test_compare_with_shape_lex(self):
        self.assertEqual(
            L1.compare_with_shape_lex(L1),
            (0, 0))

        self.assertEqual(
            L1.compare_with_shape_lex(L2),
            (0, -1))

        self.assertEqual(
            L2.compare_with_shape_lex(L1),
            (0, 1))

        t = PhyloTree(None, [L1, C23])

        self.assertEqual(
            C12.compare_with_shape_lex(t),
            L2.compare_with_shape_lex(C23))

        self.assertEqual(
            L2.compare_with_shape_lex(C23),
            (-1, 0))

        t123 = PhyloTree(None, [L1, L2, L3])

        self.assertEqual(
            t123.compare_with_shape_lex(C12),
            (1, 0))

    def test_shape(self):
        self.assertEqual(
            L1.shape(),
            L2.shape())

        self.assertEqual(
            C12.shape(),
            C31.shape())

        self.assertEqual(
            L1.shape(),
            Shape.LEAF)

        self.assertEqual(
            C23.shape(),
            Shape.CHERRY)

        self.assertEqual(
            PhyloTree(None, [L1, C31]).shape(),
            Shape([Shape.LEAF, Shape.CHERRY]))

    def test_compare_shape(self):
        self.assertEqual(
            L1.compare_with_shape_lex(L1)[0],
            L1.shape().compare(L1.shape()))

        t = PhyloTree(None, [L1, C23])

        self.assertEqual(
            C12.compare_with_shape_lex(t)[0],
            C12.shape().compare(t.shape()))

        self.assertEqual(
            L2.compare_with_shape_lex(C23)[0],
            L2.shape().compare(C23.shape()))

        t123 = PhyloTree(None, [L1, L2, L3])

        self.assertEqual(
            t123.compare_with_shape_lex(C12)[0],
            t123.shape().compare(C12.shape()))
