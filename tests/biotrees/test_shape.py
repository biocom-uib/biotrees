import unittest

from biotrees.shape import Shape, rooted_deg, unrooted_deg


class TestShape(unittest.TestCase):

    def test_init(self):
        self.assertTrue(Shape.LEAF.is_leaf())

        self.assertFalse(Shape.CHERRY.is_leaf())
        self.assertEqual(len(Shape.CHERRY.children), 2)

    def test_clone(self):
        self.assertEqual(Shape.LEAF, Shape.LEAF.clone())
        self.assertEqual(Shape.CHERRY, Shape.CHERRY.clone())

        Shape.CHERRY.clone().children[0].children = [1,2,3]
        self.assertTrue(Shape.CHERRY.children[0].is_leaf())

        Shape.CHERRY.clone().children.append(Shape.LEAF)
        self.assertEqual(len(Shape.CHERRY.children), 2)

    def test_shape(self):
        self.assertEqual(
            Shape.LEAF,
            Shape())

        self.assertEqual(
            Shape.CHERRY,
            Shape([Shape.LEAF, Shape.LEAF]))

        self.assertEqual(
            Shape.LEAF.shape(),
            Shape.LEAF)

        self.assertEqual(
            Shape.CHERRY.shape(),
            Shape.CHERRY)

        self.assertEqual(
            Shape([Shape.LEAF, Shape.CHERRY]),
            Shape([Shape.LEAF, Shape.CHERRY]))

    def test_deg(self):
        self.assertEqual(rooted_deg(Shape.LEAF), 0)
        self.assertEqual(rooted_deg(Shape.CHERRY), 2)

        self.assertEqual(unrooted_deg(Shape.LEAF), 1)
        self.assertEqual(unrooted_deg(Shape.CHERRY), 3)

        t = Shape([Shape.LEAF, Shape.CHERRY])

        self.assertEqual(rooted_deg(t), 2)
        self.assertEqual(unrooted_deg(t), 3)
        self.assertEqual(unrooted_deg(t, root=True), 2)

