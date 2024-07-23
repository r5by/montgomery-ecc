import os
import sys
import json
import unittest
import random
from mont.montgomery import Montgomery
from mont.utils import profiler
from mecc.short_weierstrass_ecc import ShortWeierstrassCurve
from sage_helper import sage_generate_points_on_curve, add_lut
from modp import GFmodp
from mecc.coordinate import PointAtInfinity, JacobianCoord
from nist_curves import F_256, SAGE_F256, p256, SAGE_p256, G256, SAGE_G256, G256_affine

import hypothesis.strategies as st
from hypothesis import given, assume, settings, example

USE_MONT = False  # montgomery domain is very slow, use it in caution
SLOW_SETTINGS = {}
if "--fast" in sys.argv:  # pragma: no cover
    SLOW_SETTINGS["max_examples"] = 2
else:
    SLOW_SETTINGS["max_examples"] = 10


class TestShortWeierstrassCurve(unittest.TestCase):

    def setUp(self):
        self.p = 29  # F = Z/pZ
        self.coeffs = [2, 20]  # coeffs [a, b] for the short weierstrass curve

        self.choices = 2  # if use mont, let us limit the test cases

        self.points_file = f"curve_points_{self.p}_{self.coeffs}.json"

        if USE_MONT:
            self.domain = Montgomery.factory(mod=self.p, mul_opt='real1').build()
            self.enable_profiler = True
        else:
            self.domain = GFmodp(self.p)

        self.curve = ShortWeierstrassCurve(self.domain, self.coeffs)
        self.affine_points, self.proj_points, self.double_points, self.addition_table = self.load_or_generate_points()
        self.cardinality = len(self.affine_points)

        if USE_MONT:
            picks = random.sample(range(0, self.cardinality), self.choices)

            self.affine_points = [self.affine_points[i] for i in picks]
            self.proj_points = [self.proj_points[i] for i in picks]
            self.double_points = [self.double_points[i] for i in picks]

        self.total = len(self.affine_points)
        print('done')

    def load_or_generate_points(self):
        if os.path.exists(self.points_file):
            with open(self.points_file, 'r') as f:
                affine_points, projective_points, double_points, addition_table = json.load(f)
                # Replace 'INF' with InfinitePoint() in affine_points
                affine_points = [PointAtInfinity() if pt == 'INF' else tuple(pt) for pt in affine_points]
                double_points = [PointAtInfinity() if pt == 'INF' else tuple(pt) for pt in double_points]
                addition_table = [PointAtInfinity() if pt == 'INF' else tuple(pt) for pt in addition_table]
        else:
            affine_points, projective_points, double_points, addition_table = sage_generate_points_on_curve(
                self.coeffs, self.p)
            with open(self.points_file, 'w') as f:
                json.dump([affine_points, projective_points, double_points, addition_table], f)
        return affine_points, projective_points, double_points, addition_table

    def test_init(self):

        # Test initialization with individual coefficients
        ec = ShortWeierstrassCurve(self.domain, self.coeffs[0], self.coeffs[1])
        self.assertEqual(ec.a, self.curve.a)
        self.assertEqual(ec.b, self.curve.b)
        self.assertEqual(ec._a.value, self.curve._a.value)
        self.assertEqual(ec._b.value, self.curve._b.value)

        # Test initialization with keyword arguments
        ec = ShortWeierstrassCurve(self.domain, a=self.coeffs[0], b=self.coeffs[1])
        self.assertEqual(ec.a, self.curve.a)
        self.assertEqual(ec.b, self.curve.b)
        self.assertEqual(ec._a.value, self.curve._a.value)
        self.assertEqual(ec._b.value, self.curve._b.value)

    def test_discriminant_zero(self):
        # A curve with discriminant zero, should raise ValueError
        with self.assertRaises(ValueError):
            ShortWeierstrassCurve(self.domain, validate=True, a=0, b=0)

    def test_type_error_on_invalid_domain(self):
        # Test that a TypeError is raised when the domain is neither int nor Montgomery
        with self.assertRaises(TypeError):
            ShortWeierstrassCurve("invalid", [2, 20])

    def test_value_error_on_incomplete_coeffs(self):
        # Test that a ValueError is raised when coefficients are incomplete
        with self.assertRaises(ValueError):
            ShortWeierstrassCurve(29, [2])  # Only one coefficient

        with self.assertRaises(ValueError):
            ShortWeierstrassCurve(29, a=2)  # Missing 'b'

    def test_affine_point_on_curve(self):
        for point in self.affine_points:
            self.assertTrue(self.curve.is_point_on_curve(point), f'point fail for {point}')

    def test_jacobian_point_on_curve(self):
        for point in self.affine_points:
            point_jaco = JacobianCoord.from_affine(point, self.domain)
            self.assertTrue(self.curve.is_point_on_curve(point_jaco))

    def test_point_off_curve(self):
        off_point = [1, 3]
        self.assertFalse(self.curve.is_point_on_curve(off_point))

        off_point_jaco = JacobianCoord.from_affine(off_point, self.domain)
        self.assertFalse(self.curve.is_point_on_curve(off_point_jaco))

    @profiler(1, enabled=USE_MONT)
    def test_affine_point_double(self):
        cnt, total = 0, self.total
        for i in range(self.total):
            p = self.affine_points[i]
            exp = self.double_points[i]

            act = self.curve.double_point_affine(p)
            self.assertEqual(exp, act)

            # print(f'done: {cnt + 1}/{total}%')
            cnt += 1

    @profiler(1, enabled=USE_MONT)
    def test_affine_point_add_normal(self):
        # g = add_lut(self.addition_table, 3, 7)

        cnt, total = 0, self.total - 1
        for i in range(self.total - 1):
            exp = add_lut(self.addition_table, i, i + 1)

            P, Q = self.affine_points[i], self.affine_points[i + 1]
            act = self.curve.add_points_affine(P, Q)
            self.assertEqual(exp, act)

            # print(f'done: {cnt + 1}/{total}%')
            cnt += 1

    def test_nist_single_double(self):
        act = p256.double_point_affine(G256_affine)
        exp = 2 * SAGE_G256
        exp = (int(exp[0]), int(exp[1]))

        self.assertEqual(act, exp)

    def test_nist_double_zeros(self):

        # test zero
        zero_point = JacobianCoord.from_int_coord(0, 0, 1, F_256)
        zero_point_sage = SAGE_p256.point([0, 0, 1], check=False)

        act = p256.double_point(zero_point).get_integer_coords()
        exp = 2 * zero_point_sage
        exp = (int(exp[0]), int(exp[1]), int(exp[2]))

        self.assertTrue(exp == (0, 1, 0))  # Inf for projective
        self.assertTrue(act == (1, 1, 0))  # Inf for jacobian

        # test zero-equivalent
        zero_point_ = JacobianCoord.from_int_coord(0, F_256.modulus, 1, F_256)
        zero_point_sage_ = SAGE_p256.point([0, F_256.modulus, 1], check=False)
        act = p256.double_point(zero_point_).get_integer_coords()
        exp = 2 * zero_point_sage_
        exp = (int(exp[0]), int(exp[1]), int(exp[2]))

        self.assertTrue(exp == (0, 1, 0))  # Inf for projective
        self.assertTrue(act == (1, 1, 0))  # Inf for jacobian

    @settings(**SLOW_SETTINGS)
    @given(
        st.integers(  # min_value is set above 0, simply get rid of k=0 and kP->Inf,
            # where for sage(projective): Inf = [0:1:0] !=
            # jacobian Inf = [1:1:0]
            min_value=1, max_value=int(F_256.modulus - 1)
        )
    )
    def test_nist_salar_mul_unfixed(self, k):
        exp = k * SAGE_G256
        exp = (int(exp[0]), int(exp[1]), int(exp[2]))

        act = p256.k_point(k, G256)
        act.to_affine()
        act = act.get_integer_coords()
        # print(f'success: k={k}')

        self.assertEqual(exp, act, f'Scalar multiplication on Curve: {p256} fails at k={k}')

    @settings(**SLOW_SETTINGS)
    @given(
        st.integers(
            min_value=1, max_value=int(F_256.modulus - 1)
        )
    )
    def test_nist_salar_mul_fixed_win(self, k):

        exp = k * SAGE_G256
        exp = (int(exp[0]), int(exp[1]), int(exp[2]))

        act1 = p256.k_point_fixed_win_naf(k, w=3, P=G256)  # 1) NAF
        act1.to_affine()
        act1 = act1.get_integer_coords()

        act2 = p256.k_point_fixed_win(k, w=5, P=G256)  # 2) windowing
        act2.to_affine()
        act2 = act2.get_integer_coords()

        # print(f'{i}-th test success: k={k}')
        self.assertEqual(exp, act1, f'Scalar multiplication on Curve: {p256} fails at k={k}')
        self.assertEqual(exp, act2, f'Scalar multiplication on Curve: {p256} fails at k={k}')

    @settings(**SLOW_SETTINGS)
    @given(
        st.integers(
            min_value=1, max_value=int(F_256.modulus - 1)
        )
    )
    def test_nist_salar_mul_fixed_comb(self, k):
        exp = k * SAGE_G256
        exp = (int(exp[0]), int(exp[1]), int(exp[2]))

        act = p256.k_point_fixed(k, w=3, P=G256)
        act.to_affine()
        act = act.get_integer_coords()
        self.assertEqual(exp, act, f'Scalar multiplication on Curve: {p256} fails at k={k}')

        # test maximum d value:
        max_t = F_256.modulus.bit_length()
        act1 = p256.k_point_fixed(k, w=4, P=G256, k_max_bits=max_t)
        act1.to_affine()
        act1 = act1.get_integer_coords()
        self.assertEqual(exp, act1, f'Scalar multiplication on Curve: {p256} fails at k={k}')


if __name__ == "__main__":
    unittest.main()
