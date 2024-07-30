import os
import sys
import json
import unittest
from mont.montgomery import Montgomery
from mont.utils import profiler
from mecc.twisted_edwards_ecc import TwistedEdwardsCurve
from sage_helper import sage_generate_points_on_curve, add_lut
from modp import GFmodp
from nist_curves import (F_448, edwards448, G448_affine, G448, SAGE_F448, SAGE_edwards448, F_25519, edwards25519,
                         G25519_affine, G25519, SAGE_F25519, SAGE_edwards25519)

import hypothesis.strategies as st
from hypothesis import given, assume, settings, example

USE_MONT = False  # montgomery domain is very slow, use it in caution
SLOW_SETTINGS = {}
if "--fast" in sys.argv:  # pragma: no cover
    SLOW_SETTINGS["max_examples"] = 2
else:
    SLOW_SETTINGS["max_examples"] = 10


class TestTwistedEdwardsCurve(unittest.TestCase):
    def setUp(self):
        self.p = 29  # F = Z/pZ

        # choose a from the following squars in the field, choose b not in the set
        # F = GF(29)
        # sqrts = [i for i in range(29) if F(i) ** 14 == 1]

        self.coeffs = [13, 21]  # coeffs [a, d] for the twisted edwards curve
        self.domain = GFmodp(self.p)

        self.points_file = f"adwards_points_{self.p}_{self.coeffs}.json"

        if USE_MONT:
            self.domain = Montgomery.factory(mod=self.p, mul_opt='real0').build()
            self.enable_profiler = True
        else:
            self.domain = GFmodp(self.p)

        self.curve = TwistedEdwardsCurve(self.domain, self.coeffs)
        self.affine_points, self.proj_points, self.double_points, self.addition_table = self.load_or_generate_points()
        self.cardinality = len(self.affine_points)

        self.total = len(self.affine_points)

    def load_or_generate_points(self):
        if os.path.exists(self.points_file):
            with open(self.points_file, 'r') as f:
                affine_points, projective_points, double_points, addition_table = json.load(f)
        else:
            affine_points, projective_points, double_points, addition_table = sage_generate_points_on_curve(
                self.coeffs, self.p, type='edwards')
            with open(self.points_file, 'w') as f:
                json.dump([affine_points, projective_points, double_points, addition_table], f)
        return affine_points, projective_points, double_points, addition_table

    def test_init_exceptions(self):
        domain = GFmodp(29)

        with self.assertRaises(ValueError) as context:
            coeff = [0, 0]
            TwistedEdwardsCurve(domain, coeff, validate=True)

        self.assertEqual("The coefficients 'a' and 'd' must be non-zero.", str(context.exception))

        with self.assertRaises(ValueError) as context:
            coeff = [13, 13]
            TwistedEdwardsCurve(domain, coeff, validate=True)

        self.assertEqual("The coefficients 'a' and 'd' must not be equal.", str(context.exception))

        with self.assertRaises(ValueError) as context:
            coeff = [12, 15]
            TwistedEdwardsCurve(domain, coeff, validate=True)

        self.assertEqual("Coefficient 'a' must be a square in the field.", str(context.exception))

        with self.assertRaises(ValueError) as context:
            coeff = [13, 20]
            TwistedEdwardsCurve(domain, coeff, validate=True)

        self.assertEqual("Coefficient 'd' must not be a square in the field.", str(context.exception))

    def test_affine_point_on_curve(self):
        for point in self.affine_points:
            self.assertTrue(self.curve.is_point_on_curve(point), f'point fail for {point}')

    def test_point_off_curve(self):
        off_point = [1, 3]
        self.assertFalse(self.curve.is_point_on_curve(off_point))

    @profiler(10, enabled=not USE_MONT)
    def test_affine_point_double(self):
        cnt, total = 0, self.total
        for i in range(self.total):
            p = self.affine_points[i]
            exp = tuple(self.double_points[i])

            act = self.curve.double_point_affine(p)
            self.assertEqual(exp, act)

            cnt += 1

    @profiler(1, enabled=not USE_MONT)
    def test_affine_point_add_normal(self):

        cnt, total = 0, self.total - 1
        for i in range(self.total - 1):
            exp = tuple(add_lut(self.addition_table, i, i + 1))

            P, Q = self.affine_points[i], self.affine_points[i + 1]

            # print(f'Adding points P={P} with Q={Q}')
            act = self.curve.add_points_affine(P, Q)
            self.assertEqual(exp, act)

            # print(f'i={i}-th test is done: {cnt + 1}/{total}%')
            cnt += 1

    def test_nist_single_double(self):

        # Double base point of goldilocks
        act = edwards448.double_point_affine(G448_affine)
        exp = SAGE_edwards448.double_point(G448_affine)
        exp = (int(exp[0]), int(exp[1]))

        self.assertEqual(act, exp)

        # Double base point of edwards25519
        act = edwards25519.double_point_affine(G25519_affine)
        exp = SAGE_edwards25519.double_point(G25519_affine)
        exp = (int(exp[0]), int(exp[1]))

        self.assertEqual(act, exp)

    @settings(**SLOW_SETTINGS)
    @given(
        st.integers(  # min_value is set above 0, simply get rid of k=0 and kP->Inf,
            # where for sage(projective): Inf = [0:1:0] !=
            # jacobian Inf = [1:1:0]
            min_value=1, max_value=int(F_448.modulus - 1)
        )
    )
    def test_nist_salar_mul_unfixed_edwards448(self, k):
        exp = SAGE_edwards448.k_points(k, G448_affine)

        act = edwards448.k_point(k, G448)
        act.to_affine()
        act = act.get_integer_coords()[:2]
        # print(f'success: k={k}')

        self.assertEqual(exp, act, f'Scalar multiplication on Curve: {edwards448} fails at k={k}')

    @settings(**SLOW_SETTINGS)
    @given(
        st.integers(  # min_value is set above 0, simply get rid of k=0 and kP->Inf,
            # where for sage(projective): Inf = [0:1:0] !=
            # jacobian Inf = [1:1:0]
            min_value=1, max_value=int(F_25519.modulus - 1)
        )
    )
    def test_nist_salar_mul_unfixed_edwards25519(self, k):
        exp = SAGE_edwards25519.k_points(k, G25519_affine)

        act = edwards25519.k_point(k, G25519)
        act.to_affine()
        act = act.get_integer_coords()[:2]
        # print(f'success: k={k}')

        self.assertEqual(exp, act, f'Scalar multiplication on Curve: {edwards25519} fails at k={k}')

    @settings(**SLOW_SETTINGS)
    @given(
        st.integers(
            min_value=1, max_value=int(F_448.modulus - 1)
        )
    )
    def test_nist_salar_mul_fixed_comb_448(self, k):
        exp = SAGE_edwards448.k_points(k, G448_affine)

        act = edwards448.k_point_fixed(k, w=3, P=G448)
        act.to_affine()
        act = act.get_integer_coords()[:2]
        self.assertEqual(exp, act, f'Scalar multiplication on Curve: {edwards448} fails at k={k}')

        # test maximum d value:
        max_t = F_448.modulus.bit_length()
        act1 = edwards448.k_point_fixed(k, w=4, P=G448, k_max_bits=max_t)
        act1.to_affine()
        act1 = act1.get_integer_coords()[:2]
        self.assertEqual(exp, act1, f'Scalar multiplication on Curve: {edwards448} fails at k={k}')

    @settings(**SLOW_SETTINGS)
    @given(
        st.integers(
            min_value=1, max_value=int(F_25519.modulus - 1)
        )
    )
    def test_nist_salar_mul_fixed_comb_25519(self, k):
        exp = SAGE_edwards25519.k_points(k, G25519_affine)

        act = edwards25519.k_point_fixed(k, w=3, P=G25519)
        act.to_affine()
        act = act.get_integer_coords()[:2]
        self.assertEqual(exp, act, f'Scalar multiplication on Curve: {edwards25519} fails at k={k}')

        # test maximum d value:
        max_t = F_25519.modulus.bit_length()
        act1 = edwards25519.k_point_fixed(k, w=4, P=G25519, k_max_bits=max_t)
        act1.to_affine()
        act1 = act1.get_integer_coords()[:2]
        self.assertEqual(exp, act1, f'Scalar multiplication on Curve: {edwards25519} fails at k={k}')
