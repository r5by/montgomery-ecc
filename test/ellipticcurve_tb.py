import unittest
from sage.all import GF, EllipticCurve as SageEllipticCurve

from ellipticcurve import EllipticCurve, EllipticCurveElt  # Adjust the import path as needed


class TestEllipticCurve(unittest.TestCase):
    def setUp(self):
        # Your elliptic curve
        self.ec = EllipticCurve(29, [4, 20])
        self.P = EllipticCurveElt(self.ec, [2, 6])
        self.Q = EllipticCurveElt(self.ec, [3, 28])

        # SageMath's elliptic curve for reference
        self.sage_ec = SageEllipticCurve(GF(29), [4, 20])
        self.sage_P = self.sage_ec(2, 6)
        self.sage_Q = self.sage_ec(3, 28)

    def test_point_addition(self):
        result = self.P + self.Q
        expected = self.sage_P + self.sage_Q
        self.assertEqual(result, EllipticCurveElt(self.ec, [expected[0], expected[1]]))

    def test_point_doubling(self):
        result = self.P + self.P
        expected = self.sage_P + self.sage_P
        self.assertEqual(result, EllipticCurveElt(self.ec, [expected[0], expected[1]]))

    def test_scalar_multiplication(self):
        result = 7 * self.P
        expected = 7 * self.sage_P
        self.assertEqual(result, EllipticCurveElt(self.ec, [expected[0], expected[1]]))

    def test_infinity_behavior(self):
        result = self.P - self.P
        expected = self.sage_P - self.sage_P
        self.assertEqual(result, EllipticCurveElt(self.ec, ("Infinity", "Infinity")))
        self.assertEqual(expected, self.sage_ec(0))  # Sage uses (0, 0) or specific representation for infinity

    def test_negation(self):
        result = -self.P
        expected = -self.sage_P
        self.assertEqual(result, EllipticCurveElt(self.ec, [expected[0], expected[1]]))


if __name__ == "__main__":
    unittest.main()
