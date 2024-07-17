from abc import ABC

from mecc.interface import EllipticCurve
from typing import Union, List
from mecc.coordinate import JacobianCoord, PointAtInfinity
from mont.typing import GFType


const_scala = lambda c, x: sum([x for _ in range(c)])  # c*x for c as a constant
double = lambda x: x + x


class ShortWeierstrassCurve(EllipticCurve, ABC):
    """
        A short weierstrass elliptic curve has the following form:
            y^2 = x^3 + ax + b  (affine)
            Y^2 = X^3 + aXZ^4 + bZ^6 (Jacobian)

            (not used)
            Y^2Z = X^3 + aXZ^2 + bZ^3 (projective)
        where a, b are coefficients in the domain
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
            if point.is_point_at_infinity():
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
        Y2, X3 = Y ** 2, X ** 3
        Z2 = Z ** 2
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
        return self.add_points(P1, P2).to_affine()

    def add_points(self, p1: JacobianCoord, p2: JacobianCoord) -> JacobianCoord:
        '''
            Add p1 and p2 on the curve in Jacobian coord.
            ref: https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-mmadd-2007-bl
        '''

        if p1.is_point_at_infinity():
            return p2
        if p2.is_point_at_infinity():
            return p1

        X1, Y1, Z1 = p1.X, p1.Y, p1.Z
        X2, Y2, Z2 = p2.X, p2.Y, p2.Z

        # TODO> Optimizations
        # if Z1 == Z2:
        #     if Z1 == 1:
        #         return self._add_points_Z1_eq_Z2_eq_1(X1, Y1, X2, Y2)
        #
        #     return self._add_points_Z1_eq_Z2(X1, Y1, X2, Y2, Z1)
        #
        # if Z2 == 1:
        #     return self._add_points_Z2_eq_1(X1, Y1, Z1, X2, Y2)

        return self._add_points(X1, Y1, Z1, X2, Y2, Z2)

    def _add_points_Z1_eq_Z2_eq_1(self, X1, Y1, X2, Y2):
        pass

    def _add_points_Z1_eq_Z2(self, X1, Y1, X2, Y2, Z1):
        pass

    def _add_points_Z2_eq_1(self, X1, Y1, Z1, X2, Y2):
        pass

    def _add_points(self, X1, Y1, Z1, X2, Y2, Z2):

        Z1Z1 = Z1 ** 2  # 1S
        Z2Z2 = Z2 ** 2  # 2S
        U1 = X1 * Z2Z2  # 1M
        U2 = X2 * Z1Z1  # 2M
        S1 = Y1 * Z2 * Z2Z2  # 4M
        S2 = Y2 * Z1 * Z1Z1  # 6M
        H = U2 - U1
        I = double(H) ** 2  # 1*2, 3S
        J = H * I  # 7M
        r = double(S2 - S1)  # 2*2
        V = U1 * I  # 8M

        X3 = r ** 2 - J - double(V)  # 4S, 3*2
        Y3 = r * (V - X3) - double(S1 * J)  # 10M, 4*2
        Z3 = ((Z1 + Z2) ** 2 - Z1Z1 - Z2Z2) * H  # 11M, 5S, 4*2

        return JacobianCoord.from_domain_coord(X3, Y3, Z3)

    def double_point_affine(self, p: Union[List[int], PointAtInfinity]) -> Union[List[int], PointAtInfinity]:
        if isinstance(p, PointAtInfinity):
            return p

        P = JacobianCoord.from_affine(p, domain=self.domain)
        return self.double_point(P).to_affine()

    def double_point(self, p: JacobianCoord) -> JacobianCoord:
        '''
            Double p on the curve in Jacobian coord.
        '''
        X, Y, Z = p.X, p.Y, p.Z

        if Y == 0:
            return JacobianCoord.point_at_infinity(self.domain)

        # S = self.domain(4) * X * Y ** 2
        S = const_scala(4, X * Y ** 2)
        # M = self.domain(3) * X ** 2 + self._a * Z ** 4
        M = const_scala(3, X ** 2) + self._a * Z ** 4

        # X_ = M ** 2 - self.domain(2) * S
        X_ = M ** 2 - double(S)
        # Y_ = M * (S - X_) - self.domain(8) * Y ** 4
        Y_ = M * (S - X_) - const_scala(8, Y ** 4)
        # Z_ = self.domain(2) * Y * Z
        Z_ = double(Y * Z)

        return JacobianCoord.from_domain_coord(X_, Y_, Z_)

    def k_point_affine(self, k, p):
        pass

    # todo> NAF for unfixed, COMB for fixed (P is known and pre-calculated)
    def k_point(self, k, p):
        # Naive implementation of scalar multiplication
        result = (None, None)
        addend = p
        while k:
            if k & 1:
                if result == (None, None):
                    result = addend
                else:
                    result = self.add_points(result, addend)
            addend = self.double_point(addend)
            k >>= 1
        return result
