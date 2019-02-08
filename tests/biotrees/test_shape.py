import unittest

from biotrees.shape import Shape, LEAF, CHERRY, rooted_deg, unrooted_deg


class TestShape(unittest.TestCase):

    def test_init(self):
        self.assertTrue(LEAF.is_leaf())
        self.assertEqual(LEAF.leaf, '1')

        self.assertFalse(CHERRY.is_leaf())
        self.assertEqual(len(CHERRY.children), 2)

    def test_clone(self):
        self.assertEqual(LEAF, LEAF.clone())

        LEAF.clone().leaf = '2'
        self.assertEqual(LEAF.leaf, '1')

        self.assertEqual(CHERRY, CHERRY.clone())

        CHERRY.clone().children.append(LEAF)
        self.assertEqual(len(CHERRY.children), 2)

    def test_shape(self):
        self.assertEqual(
            LEAF,
            Shape())

        self.assertEqual(
            CHERRY,
            Shape([LEAF, LEAF]))

        self.assertEqual(
            LEAF.shape(),
            Shape.LEAF)

        self.assertEqual(
            CHERRY.shape(),
            Shape.CHERRY)

        self.assertEqual(
            Shape([LEAF, CHERRY]),
            Shape([Shape.LEAF, Shape.CHERRY]))

    def test_deg(self):
        self.assertTrue(LEAF)
        self.assertEqual(rooted_deg(self), 0)
        self.assertEqual(unrooted_deg(self), 1)
        self.assertEqual(unrooted_deg(self, True), 0)

        self.assertTrue(CHERRY)
        self.assertEqual(rooted_deg(self), 2)
        self.assertEqual(unrooted_deg(self), 3)
        self.assertEqual(unrooted_deg(self, True), 3)

