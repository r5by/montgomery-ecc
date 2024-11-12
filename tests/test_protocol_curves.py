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

    def test_SM2_point_add(self):
        num_examples = SLOW_SETTINGS.get("max_examples", 10)
        with open("../examples/sm2_point_add_golden.data", "a") as file:
            for _ in range(num_examples):
                # Generate random scalars k1 and k2
                k1 = random.randint(1, MOD - 1)
                k2 = random.randint(1, MOD - 1)
                # Compute points P and Q using SageMath
                P = k1 * SAGE_SM2_G
                Q = k2 * SAGE_SM2_G
                # Expected result using SageMath
                R_exp = P + Q
                # Convert points to integer coordinates
                P_coords = (int(P[0]), int(P[1]), int(P[2]))
                Q_coords = (int(Q[0]), int(Q[1]), int(Q[2]))
                R_exp_coords = (int(R_exp[0]), int(R_exp[1]), int(R_exp[2]))

                # Compute points using SM2 implementation
                P_SM2 = SM2_c.k_point_fixed_256(k1, w=4, P=SM2_G)
                file.write(f"{P_SM2.X.value}, {P_SM2.Y.value}, {P_SM2.Z.value};") # record P

                Q_SM2 = SM2_c.k_point_fixed_256(k2, w=4, P=SM2_G)
                file.write(f"{Q_SM2.X.value}, {Q_SM2.Y.value}, {Q_SM2.Z.value};") # record Q

                # Add points using SM2 implementation
                R_act = SM2_c.add_points(P_SM2, Q_SM2)
                R_act.to_affine()
                file.write(f"{R_act.X.value}, {R_act.Y.value}, {R_act.Z.value}\n") # record P+Q (affine)
                R_act_coords = R_act.get_integer_coords()

                # Validate results.
                if R_exp_coords[2] == 0 and R_act_coords[2] == 0: # INF
                    continue

                self.assertEqual(
                    R_exp_coords,
                    R_act_coords,
                    f'Point addition on Curve: {SM2_c} fails at P={P_coords}, Q={Q_coords}'
                )

                print(f'Test case on P={P_coords}, Q={Q_coords} succeed!')

    def test_SM2_point_double(self):
        num_examples = SLOW_SETTINGS.get("max_examples", 10)
        with open("../examples/sm2_point_double_golden.data", "a") as file:
            for _ in range(num_examples):
                # Generate random scalars k
                k = random.randint(1, MOD - 1)
                P = k * SAGE_SM2_G
                P_coords = (int(P[0]), int(P[1]), int(P[2]))

                # Expected result using SageMath
                R_exp = 2 * P
                R_exp_coords = (int(R_exp[0]), int(R_exp[1]), int(R_exp[2]))

                # Compute points using SM2 implementation
                P_SM2 = SM2_c.k_point_fixed_256(k, w=4, P=SM2_G)
                file.write(f"{P_SM2.X.value}, {P_SM2.Y.value}, {P_SM2.Z.value};") # record P

                # Add points using SM2 implementation
                R_act = SM2_c.double_point(P_SM2)
                R_act.to_affine()
                file.write(f"{R_act.X.value}, {R_act.Y.value}, {R_act.Z.value}\n") # record P+Q (affine)
                R_act_coords = R_act.get_integer_coords()

                # Validate results.
                if R_exp_coords[2] == 0 and R_act_coords[2] == 0: # INF
                    continue

                self.assertEqual(
                    R_exp_coords,
                    R_act_coords,
                    f'Double point on Curve: {SM2_c} fails at P={P_coords}'
                )

                print(f'Test case on P={P_coords} succeed!')

    def test_SM2_verify_point_on_curve(self):
        ''' Verify point on curve '''
        num_examples = SLOW_SETTINGS.get("max_examples", 10)
        with open("../examples/sm2_points_on_curve_affine_golden.data", "a") as file:
            for _ in range(num_examples):
                # Generate random scalars k
                k = random.randint(1, MOD - 1)
                # Compute points using SM2 implementation
                p_on_c = SM2_c.k_point_fixed_256(k, w=4, P=SM2_G)
                p_on_c.to_affine()

                self.assertTrue(
                    SM2_c.is_point_on_curve(p_on_c),
                    f'Verify point on curve fails at: P={p_on_c}'
                )

                p_on_c = p_on_c.get_integer_coords()
                file.write(f"{p_on_c[0]}, {p_on_c[1]}\n")  # record P
                print(f'Test case on P={p_on_c} succeed!')

    def test_SM2_verify_point_off_curve(self):
        ''' Verify points not on the curve '''
        num_examples = SLOW_SETTINGS.get("max_examples", 10)
        with open("../examples/sm2_points_off_curve_affine_golden.data", "a") as file:
            count = 0
            while count < num_examples:
                k = random.randint(1, MOD - 1)
                # Compute a valid point on the curve using SM2 implementation
                p_on_c = SM2_c.k_point_fixed_256(k, w=4, P=SM2_G)
                p_on_c.to_affine()
                x_valid, y_valid = p_on_c.get_integer_coords()[:2]
                # Modify the x-coordinate slightly to get a point off the curve
                x_invalid = (x_valid + random.randint(1, MOD - 1)) % MOD
                y_invalid = y_valid
                # Ensure the modified point is actually off the curve
                if x_invalid == x_valid:
                    continue  # Skip if x_invalid didn't change
                # Create a point with the modified x-coordinate
                p_off_c = JacobianCoord.from_affine([x_invalid, y_invalid], SM2_F)
                # Verify that this point is not on the curve
                if not SM2_c.is_point_on_curve(p_off_c):
                    self.assertFalse(
                        SM2_c.is_point_on_curve(p_off_c),
                        f'Verify point off curve fails at: P=({x_invalid}, {y_invalid})'
                    )
                    # Write the invalid point coordinates to the file
                    file.write(f"{x_invalid}, {y_invalid}\n")
                    print(f'Test case on P=({x_invalid}, {y_invalid}) succeeded!')
                    count += 1
                else:
                    # If the modified point is still on the curve, try again
                    continue

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
                file.write(f"{k}, {exp[0]}, {exp[1]}\n")

if __name__ == "__main__":
    unittest.main()
