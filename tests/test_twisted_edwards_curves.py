import os
import sys
import json
import unittest
import random
from sage.all import GF
from mont.montgomery import Montgomery
from mont.utils import profiler
from mecc.twisted_edwards_ecc import TwistedEdwards
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

        self.curve = TwistedEdwards(self.domain, self.coeffs)
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
            TwistedEdwards(domain, coeff, validate=True)

        self.assertEqual("The coefficients 'a' and 'd' must be non-zero.", str(context.exception))

        with self.assertRaises(ValueError) as context:
            coeff = [13, 13]
            TwistedEdwards(domain, coeff, validate=True)

        self.assertEqual("The coefficients 'a' and 'd' must not be equal.", str(context.exception))

        with self.assertRaises(ValueError) as context:
            coeff = [12, 15]
            TwistedEdwards(domain, coeff, validate=True)

        self.assertEqual("Coefficient 'a' must be a square in the field.", str(context.exception))

        with self.assertRaises(ValueError) as context:
            coeff = [13, 20]
            TwistedEdwards(domain, coeff, validate=True)

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
