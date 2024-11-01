import os
import sys
import json
import unittest
import random
from mecc.coordinate import PointAtInfinity, JacobianCoord
from nist_curves import F_256, SAGE_F256, p256, SAGE_p256, G256, SAGE_G256, G256_affine
from sm2_curve import  SAGE_SM2_G, SM2_G, SM2_F, SM2_c, MOD

import hypothesis.strategies as st
from hypothesis import given, assume, settings, example

# USE_MONT = True  # montgomery domain is very slow, use it in caution
SLOW_SETTINGS = {}
if "--fast" in sys.argv:  # pragma: no cover
    SLOW_SETTINGS["max_examples"] = 2
else:
    SLOW_SETTINGS["max_examples"] = 10


class TestNISTCurves(unittest.TestCase):

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
    def test_nist_salar_mul_fixed_comb(self):

        num_examples = SLOW_SETTINGS.get("max_examples", 10)
        ks = [random.getrandbits(random.randint(1, F_256.modulus.bit_length() - 1)) for _ in range(num_examples)]
        ks = [k if k < F_256.modulus else k % F_256.modulus for k in ks]  # Ensure k is within the modulus range

        for k in ks:
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
            print(f'Test case on k={k} succeed!')

    @settings(**SLOW_SETTINGS)
    def test_SM2_salar_mul_fixed_comb(self):

        ##region test single k
        # k = 486874
        # exp = k * SAGE_SM2_G
        # exp =  (int(exp[0]), int(exp[1]), int(exp[2]))
        #
        # act = SM2_c.k_point_fixed(k, w=4, P=SM2_G)
        # act.to_affine()
        # act = act.get_integer_coords()
        # self.assertEqual(exp, act)
        ##endregion

        num_examples = SLOW_SETTINGS.get("max_examples", 10)
        ks = [random.getrandbits(random.randint(1, MOD.bit_length() - 1)) for _ in range(num_examples)]
        ks = [k if k < MOD else k % MOD for k in ks]  # Ensure k is within the modulus range
        with open("../examples/sm2_golden.data", "a") as file:  # Open the file in append mode
            for k in ks:
                exp = k * SAGE_SM2_G
                exp = (int(exp[0]), int(exp[1]), int(exp[2]))

                act = SM2_c.k_point_fixed_256(k, w=4, P=SM2_G)
                act.to_affine()
                act = act.get_integer_coords()

                if exp[2] == 0 and act[2] == 0:
                    continue

                self.assertEqual(exp, act, f'Scalar multiplication on Curve: {SM2_c} fails at k={k}')

                print(f'Test case on k={k} succeed!')
                # Write the values of k, int(exp[0]), int(exp[1]) into the file
                # file.write(f"{k}, {exp[0]}, {exp[1]}\n")

if __name__ == "__main__":
    unittest.main()
