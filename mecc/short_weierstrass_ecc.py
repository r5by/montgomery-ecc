from abc import ABC

from mecc.ecc import EllipticCurve
from typing import Union, List, Optional
from mecc.coordinate import JacobianCoord, PointAtInfinity
from mont.typing import GFType, GFElementType
from mont.common import ith_word
from mecc.utils import naf_prodinger, int_by_slider, msb, turn_off_msb
import math

from mecc.utils import fe_double, fe_const_scala


class ShortWeierstrassCurve(EllipticCurve, ABC):
    """
        A short weierstrass elliptic curve has the following form:
            y^2 = x^3 + ax + b  (affine)
            Y^2 = X^3 + aXZ^4 + bZ^6 (Jacobian)

            (not used)
            Y^2Z = X^3 + aXZ^2 + bZ^3 (projective)
        where a, b are coefficients in the domain
        References:
        [1] Guide to elliptic curve cryptography
        [2] https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-2007-bl
        [3] SP 800-186
    """

    def __init__(self, domain: GFType, *coeffs: Union[List[int], int, None], validate=False, **kwargs):
        # Determine if coefficients are provided directly, in a list, or as keyword arguments
        if len(coeffs) == 1 and isinstance(coeffs[0], list):  # case 1, like ec = EllipticCurve(29, [2, 20])
            a, b = coeffs[0]
        elif len(coeffs) == 2 and all(isinstance(x, int) for x in coeffs):  # case 2, like ec = EllipticCurve(29, 2, 20)
            a, b = coeffs
        elif 'a' in kwargs and 'b' in kwargs:  # case 3 like ec = EllipticCurve(29, a=2, b=20)
            a = kwargs['a']
            b = kwargs['b']
        else:
            raise ValueError(
                "Coefficients must be provided either as list [a, b], direct integers a, b, or as separate keyword "
                "arguments a, b")

        # domain associated with this elliptic curve, can be Montgomery domain or any others
        # validation? leaves it to the duck itself...
        self.domain = domain

        # Coefficients in Integer domain
        self.a = a
        self.b = b
        # Repr. of a, b in the domain
        self._a = self.domain(a)
        self._b = self.domain(b)

        if validate:
            self._validate_curve()

    def _validate_curve(self):
        """Private method to calculate and check the discriminant to ensure the curve is valid."""
        # Discriminant of the short Weierstrass curve: Δ = -16(4a³ + 27b²)
        discriminant = self.domain(4) * self._a ** 3 + self.domain(27) * self._b ** 2
        if discriminant == 0:
            raise ValueError("Invalid curve parameters; discriminant is zero.")

    def is_point_on_curve(self, point: Union[List[int], JacobianCoord, PointAtInfinity]):
        """
            Verify if given point is on curve, a "point" format is like
            1) [int, int]:  affine coord.
            2) [int, int, int]: jacobian coord.
            3) JacobianCoord: already in jacobian
        """
        if isinstance(point, JacobianCoord):
            if point.is_identity_point():
                return True

            if point.domain is None:
                point.enter_domain(self.domain)

            return self.verify_point_on_curve(point)

        if isinstance(point, PointAtInfinity):
            return True

        if len(point) == 2:
            point = JacobianCoord.from_affine(point, self.domain)
        elif len(point) == 3:
            X, Y, Z = point
            point = JacobianCoord.from_int_coord(X, Y, Z, self.domain)
        else:
            raise ValueError(f'Unrecognized coordinates: {point}')

        return self.verify_point_on_curve(point)

    def verify_point_on_curve(self, point: JacobianCoord):
        """
            Verify a point on the short weierstrass curve using Jacobian coordinates.
            *convert to mont domain first
        """
        X, Y, Z = point.X, point.Y, point.Z

        # Using the Jacobian coordinates curve equation: Y^2 = X^3 + aXZ^4 + bZ^6
        # Y2, X3 = Y ** 2, X ** 3
        Y2 = Y * Y
        X2 = X * X
        X3 = X2 * X

        Z2 = Z * Z
        Z4 = Z2 * Z2
        Z6 = Z4 * Z2

        AXZ4 = self._a * X * Z4
        BZ6 = self._b * Z6
        RHS = X3 + AXZ4 + BZ6

        return Y2 == RHS

    def add_points_affine(self, p1: List[int], p2: List[int]) -> Union[List[int], PointAtInfinity]:
        '''
            Add two points in affine coordinate, and return the result point in affine coordinate
            * We skip the `normal` way of calculating point addition on affine plan, directly apply jacobian addition
        '''
        P1 = JacobianCoord.from_affine(p1, domain=self.domain)
        P2 = JacobianCoord.from_affine(p2, domain=self.domain)
        R = self.add_points(P1, P2)
        return R.get_affine_coords()

    def add_points(self, p1: JacobianCoord, p2: JacobianCoord) -> JacobianCoord:
        '''
            Add p1 and p2 on the curve in Jacobian coord.
            ref: https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-mmadd-2007-bl
        '''

        # !! Important, or else this impl fails !!
        if p1.is_identity_point():
            return p2
        if p2.is_identity_point():
            return p1

        if p1 == p2:  # !! Important, or else this impl fails !!
            return self.double_point(p1)

        X1, Y1, Z1 = p1.X, p1.Y, p1.Z
        X2, Y2, Z2 = p2.X, p2.Y, p2.Z

        # Optimizations
        if Z1 == Z2:
            if Z1 == 1:
                return self._add_with_z_eq_1(X1, Y1, X2, Y2)

            return self._add_with_z_eq(X1, Y1, X2, Y2, Z1)

        if Z2 == 1:
            return self._add_with_z2_1(X1, Y1, Z1, X2, Y2)

        return self._add_with_z_ne(X1, Y1, Z1, X2, Y2, Z2)

    def _add_with_z_eq_1(self, X1, Y1, X2, Y2):
        """add points when both Z1 and Z2 equal 1"""
        # after:
        # http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-mmadd-2007-bl
        H = X2 - X1
        HH = H * H
        I = fe_const_scala(4, HH)
        J = H * I
        r = fe_double(Y2 - Y1)
        # ~~ H==r==0 indicates P==Q here ~~
        V = X1 * I

        rr = r * r
        X3 = rr - J - fe_double(V)
        Y3 = r * (V - X3) - fe_double(Y1) * J
        Z3 = fe_double(H)
        return JacobianCoord.copy(X3, Y3, Z3, self.domain)

    def _add_with_z_eq(self, X1, Y1, X2, Y2, Z1):
        """add points when Z1 == Z2"""
        # after:
        # http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-zadd-2007-m
        # A = (X2 - X1) ** 2
        _X21 = X2 - X1
        A = _X21 * _X21

        B = X1 * A
        C = X2 * A

        # D = (Y2 - Y1) ** 2
        _Y21 = Y2 - Y1
        D = _Y21 * _Y21

        # ~~ A == D == 0 indicates point double here ~~
        X3 = D - B - C
        Y3 = (Y2 - Y1) * (B - X3) - Y1 * (C - B)
        Z3 = Z1 * (X2 - X1)
        return JacobianCoord.copy(X3, Y3, Z3, self.domain)

    def _add_with_z2_1(self, X1, Y1, Z1, X2, Y2):
        """add points when Z2 == 1"""
        # after:
        # http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-madd-2007-bl
        Z1Z1 = Z1 * Z1
        U2 = X2 * Z1Z1
        S2 = Y2 * Z1 * Z1Z1
        H = U2 - X1
        HH = H * H
        I = fe_const_scala(4, HH)
        J = H * I
        r = fe_double(S2 - Y1)
        # ~~ r == H == 0 indicates point double here ~~
        V = X1 * I

        X3 = r * r - J - fe_double(V)
        Y3 = r * (V - X3) - fe_double(Y1) * J

        Z1_H = Z1 + H
        _Z1_H2 = Z1_H * Z1_H
        # Z3 = (Z1 + H) ** 2 - Z1Z1 - HH
        Z3 = _Z1_H2 - Z1Z1 - HH
        return JacobianCoord.copy(X3, Y3, Z3, self.domain)

    def _add_with_z_ne(self, X1, Y1, Z1, X2, Y2, Z2):
        #region: version(1) after: http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-2007-bl
        # # Z1Z1 = Z1 ** 2  # 1S
        # Z1Z1 = Z1 * Z1
        #
        # # Z2Z2 = Z2 ** 2  # 2S
        # Z2Z2 = Z2 * Z2  # 2S
        #
        # U1 = X1 * Z2Z2  # 1M
        # U2 = X2 * Z1Z1  # 2M
        # S1 = Y1 * Z2 * Z2Z2  # 4M
        # S2 = Y2 * Z1 * Z1Z1  # 6M
        # H = U2 - U1
        #
        # # I = fe_double(H) ** 2  # 1*2, 3S
        # _2H = fe_double(H)
        # I = _2H * _2H  # 1*2, 3S
        #
        # J = H * I  # 7M
        # r = fe_double(S2 - S1)  # 2*2
        #
        # # IMPORTANT NOTE:
        # # This add_point formula can't be use for adding two same point (doubling),
        # # as H will be Zero that cause incorrect result!!
        # # ~~ H == r == 0 indicates a point double ~~
        # V = U1 * I  # 8M
        #
        # rr = r * r
        # X3 = rr - J - fe_double(V)  # 4S, 3*2
        # Y3 = r * (V - X3) - fe_double(S1 * J)  # 10M, 4*2
        #
        # Z1_Z2 = Z1 + Z2
        # _Z1_Z2_2 = Z1_Z2 * Z1_Z2
        # # Z3 = ((Z1 + Z2) ** 2 - Z1Z1 - Z2Z2) * H  # 11M, 5S, 4*2
        # Z3 = (_Z1_Z2_2 - Z1Z1 - Z2Z2) * H  # 11M, 5S, 4*2
        #endregion

        #region: version(2) after: http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-1998-cmo-2
        Z1Z1 = Z1 * Z1
        Z2Z2 = Z2 * Z2
        U1 = X1 * Z2Z2
        U2 = X2 * Z1Z1
        S1 = Y1 * Z2 * Z2Z2
        S2 = Y2 * Z1 * Z1Z1

        H = U2 - U1
        HH = H * H
        HHH = H * HH
        r = S2 - S1
        V = U1 * HH

        rr = r * r
        _2V = fe_double(V)

        X3 = rr - HHH - _2V
        Y3 = r * (V - X3) - S1 * HHH
        Z3 = Z1 * Z2 * H
        #endregion

        return JacobianCoord.copy(X3, Y3, Z3, self.domain)

    def double_point_affine(self, p: Union[List[int], PointAtInfinity]) -> Union[List[int], PointAtInfinity]:
        if isinstance(p, PointAtInfinity):
            return p

        P = JacobianCoord.from_affine(p, domain=self.domain)
        return self.double_point(P).get_affine_coords()

    def double_point(self, p: JacobianCoord) -> JacobianCoord:
        '''
            Double p on the curve in Jacobian coord.
        '''
        X, Y, Z = p.X, p.Y, p.Z

        if Y == 0:
            return JacobianCoord.get_identity_point(self.domain)

        if Z == 1:
            return self._double_with_z_1(X, Y)

        return self._double(X, Y, Z)

    def _double(self, X1: GFElementType, Y1: GFElementType, Z1: GFElementType) -> JacobianCoord:
        #region version(2) http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-1998-cmo-2
        XX = X1 * X1
        YY = Y1 * Y1
        ZZ = Z1 * Z1

        t0 =  X1 * YY
        S = fe_const_scala(4, t0)

        t1 = ZZ * ZZ
        t2 = self._a * t1
        t3 = fe_const_scala(3, XX)
        M = t3 + t2

        t4 = M * M
        t5 = fe_double(S)
        T = t4 - t5

        X3 = T
        t6 = S - T
        t7 = YY * YY
        t8 = fe_const_scala(8, t7)
        t9 = M * t6
        Y3 = t9 - t8
        t10 = Y1 * Z1
        Z3 = fe_double(t10)

        #endregion

        #region: version(1) http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-2007-bl
        # XX = X1 * X1
        # YY = Y1 * Y1
        # ZZ = Z1 * Z1
        #
        # S = fe_const_scala(4, X1 * YY)
        #
        # ZZZZ = ZZ * ZZ
        # M = fe_const_scala(3, XX) + self._a * ZZZZ
        #
        # MM = M * M
        # T = MM - fe_double(S)
        #
        # X3 = T
        # YYYY = YY * YY
        # Y3 = M * (S - T) - fe_const_scala(8, YYYY)
        # Z3 = fe_double(Y1 * Z1)
        #endregion

        return JacobianCoord.copy(X3, Y3, Z3, self.domain)


    def _double_with_z_1(self, X1: GFElementType, Y1: GFElementType) -> JacobianCoord:
        """Add a point to itself with z == 1."""
        # after:
        # http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-mdbl-2007-bl
        XX = X1 * X1
        YY = Y1 * Y1
        YYYY = YY * YY
        S = fe_double((X1 + YY) ** 2 - XX - YYYY)
        M = fe_const_scala(3, XX) + self._a
        T = M * M - fe_double(S)

        # X3 = T
        Y3 = M * (S - T) - fe_const_scala(8, YYYY)
        Z3 = fe_double(Y1)
        return JacobianCoord.copy(T, Y3, Z3, self.domain)

    def k_point_affine(self, k, p):
        pass  # todo>

    def k_point_fixed_win(self, k: int, w: int, P: JacobianCoord) -> JacobianCoord:
        ''' Fixed-base windowing method for fixed-point multiplication
            ref: [1] Algorithm 3.41
        '''
        if k == 0:
            return JacobianCoord.get_identity_point(self.domain)

        # Prepare to perform the multiplication
        A = B = JacobianCoord.get_identity_point(self.domain)

        t = k.bit_length()
        d = math.ceil(t / w)
        precomputed = self._precompute_win(w, d, P)

        # [_.to_affine() for _ in precomputed]  # debug usage

        for j in range((1 << w) - 1, 0, -1):

            for i in range(d):
                Ki = ith_word(k, i, w)

                if Ki != j:
                    continue

                B = self.add_points(B, precomputed[i])

            A = self.add_points(A, B)

        #region naf
        # kp, kn = naf_prodinger(k)
        #
        # for i in range(d - 1, -1, -1):
        #     Q = self.double_point(Q)
        #     # ki = int_by_slider(k, d, i)
        #     kpi = int_by_slider(kp, d, i)
        #     kni = int_by_slider(kn, d, i)
        #     ki = kpi - kni
        #
        #     if ki > 0:
        #         Q = self.add_points(Q, precomputed[ki])
        #     else:
        #         Q = self.add_points(Q, -precomputed[-ki])
        #endregion

        return A

    def k_point_fixed_win_naf(self, k: int, w: int, P: JacobianCoord) -> JacobianCoord:
        ''' Fixed-base NAF windowing method for fixed-point multiplication
            ref: [1] Algorithm 3.42
        '''
        if k == 0:
            return JacobianCoord.get_identity_point(self.domain)

        A = B = JacobianCoord.get_identity_point(self.domain)

        kp, kn = naf_prodinger(k)
        t = max(kp.bit_length(), kn.bit_length())
        d = math.ceil(t / w)
        precomputed = self._precompute_win(w, d, P)

        # [_.to_affine() for _ in precomputed]  # debug usage

        for j in range((1 << w) - 1, 0, -1):

            for i in range(d):
                Kpi, Kni = ith_word(kp, i, w), ith_word(kn, i, w)

                if Kpi == j:
                    B = self.add_points(B, precomputed[i])

                if Kni == j:
                    B = self.add_points(B, -precomputed[i])

            A = self.add_points(A, B)

        return A
