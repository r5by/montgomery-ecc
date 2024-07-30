from abc import ABC

from mecc.ecc import EllipticCurve
from typing import Union, List
from mecc.coordinate import ExtendedCoord, PointAtInfinity
from mont.typing import GFType, GFElementType

const_scala = lambda c, x: sum([x for _ in range(c)])  # c*x for c as a constant
double = lambda x: x + x


class TwistedEdwardsCurve(EllipticCurve, ABC):
    """
        A twisted edwards curve has the following form:
            ax^2 + y^2 = 1 + dx^2y^2  (affine)

            by translates the coordinates into projective space by x=X/Z, y=Y/Z we have:
            (then multiply by Z^2)
            aX^2 + Y^2 = Z^2 + dX^2Y^2/Z^2  <== use this one

            by taking the extended coordinate T=XY/Z we have:
            aX^2 + Y^2 = Z^2 + dT^2  <== use this one

            Y^2 = X^3 + aXZ^4 + bZ^6 (Jacobian) <= not used!
        where a, d are coefficients in the domain

        An Edwards curve has the definition:
            x^2 + y^2 = c^2*(1 + dx^2y^2)
        where c and d are coefficients in the domain(if c==0, reduce the case back to twisted edwards)

        SP 800-186 specified:

            # Twisted Edward curve (a != 1)
            Edwards25519:
                p = 2^255 - 19,         a = -1,     d = -121665/121666


            # Edwards curves (a == 1):
            Edwards448:
                p = 2^448 - 2^224 - 1,  a = 1,      d = -39081

            E448:
                p = 2^448 - 2^224 - 1,  a = 1,      d = 39082/39081

        NOTE!!
            The addition law (affine) defines the group on the twisted edwards curves, that has
            the following properties (plugin these values to the formula shall check!):

            1) The identity point is (0, 1) which has order of 1  (always on the curve)
            2) The point (0, -1) is always on the curve (which has order of 2)
            3) the points (x, y) and (-x, y) are inverse in this group

            if a == 1, then this is an Edwards curve, additional to |(0,1)|=1 and |(0, -1)|=2
            it has |(+-1, 0)| = 4 that are on curve. Note for twisted edwards curve, it doesn't necessarily have
            an element of order 4, but its order shall be divided by 4 as well!

        References:
            [1] NIST SP 800-186
            [2] <<Twisted Edwards Curves Revisted 2008>>
    """

    def __init__(self, domain: GFType, *coeffs: Union[List[int], int, None], validate=False, **kwargs):
        # Determine if coefficients are provided directly, in a list, or as keyword arguments
        if len(coeffs) == 1 and isinstance(coeffs[0], list):  # case 1, like ec = TwistedEdwards(29, [2, 20])
            a, d = coeffs[0]
        elif len(coeffs) == 2 and all(
                isinstance(x, int) for x in coeffs):  # case 2, like ec = TwistedEdwards(29, 2, 20)
            a, d = coeffs
        elif 'a' in kwargs and 'd' in kwargs:  # case 3 like ec = TwistedEdwards(29, a=2, b=20)
            a = kwargs['a']
            d = kwargs['d']
        else:
            raise ValueError(
                "Coefficients must be provided either as list [a, d], direct integers a, d, or as separate keyword "
                "arguments a, d")

        self.domain = domain

        # Coefficients in Integer domain
        self.a = a
        self.d = d
        # Repr. of a, d in the domain
        self._a = self.domain(a)
        self._d = self.domain(d)

        if validate:
            self._validate_curve()

    def _validate_curve(self):
        """ Validate edwards's coefficients, following [1] 2.1.3:
          1) a, d != 0
          2) a != d
          3) a is a square in GF(p)
          4) b is not a square in GF(p)
         """

        def _is_square(x: GFElementType):
            # A number x is a square in GF(p) if it's a quadratic residue modulo p.
            # This can be checked by raising it to the power of (p-1)/2
            power = self.domain.order - 1 >> 1
            return x ** power == 1

        if self._a == 0 or self._d == 0:
            raise ValueError("The coefficients 'a' and 'd' must be non-zero.")
        if self._a == self._d:
            raise ValueError("The coefficients 'a' and 'd' must not be equal.")
        if not _is_square(self._a):
            raise ValueError("Coefficient 'a' must be a square in the field.")
        if _is_square(self._d):
            raise ValueError("Coefficient 'd' must not be a square in the field.")

    def is_point_on_curve(self, point: Union[List[int], ExtendedCoord]):
        """
            Verify if given point is on curve, a "point" format is like
            1) [int, int]:  affine coord.
            2) [int, int, int, int]: extended coord.
            3) ExtendedCoord: already in extended coord.
        """
        if isinstance(point, ExtendedCoord):
            if point.is_identity_point():
                return True

            if point.domain is None:
                point.enter_domain(self.domain)

            return self.verify_point_on_curve(point)

        if len(point) == 2:
            point = ExtendedCoord.from_affine(point, self.domain)
        elif len(point) == 3:
            X, Y, Z = point
            point = ExtendedCoord.from_int_coord(X, Y, z=Z, domain=self.domain)
        else:
            raise ValueError(f'Unrecognized coordinates: {point}')

        return self.verify_point_on_curve(point)

    def verify_point_on_curve(self, point: ExtendedCoord) -> bool:
        """
            Verify a point on the twisted edwards curve using extended coordinates.
            *convert to mont domain first
        """
        X, Y, Z, T = point.X, point.Y, point.Z, point.T

        # Using the extended coordinates curve equation aX^2 + Y^2 = Z^2 + dT^2
        X2 = X * X
        Y2 = Y * Y
        Z2 = Z * Z
        T2 = T * T

        LHS = self._a * X2 + Y2
        RHS = Z2 + self._d * T2
        return LHS == RHS

    def _add_with_z_eq_1(self, X1, Y1, T1, X2, Y2, T2):
        """
            <deprecated> add points when both Z1 and Z2 equal 1
            NOTE: WILL fail for (x,y)+(-x,-y)
        """

        #region Againd, dedicated addition fail for (x, y) + (-x, -y)
        A = X1 * X2
        B = Y1 * Y2
        C = T2
        D = T1

        E = D + C
        F = (X1 - Y1) * (X2 + Y2) + B - A
        G = B + self._a * A
        H = D - C

        X3 = E * F
        Y3 = G * H
        T3 = E * H
        Z3 = F * G

        return ExtendedCoord(X3, Y3, T3, Z3, domain=self.domain)

    def _add_with_z2_1(self, X1, Y1, Z1, T1, X2, Y2, T2):
        """Add two points if the second point is on affine plane"""

        #region Dedicated failed for (x,y)+(-x,-y)
        # A = X1 * X2
        # B = Y1 * Y2
        # C = Z1 * T2
        # D = T1
        #
        # E = D + C
        # F = (X1 - Y1) * (X2 + Y2) + B - A
        # G = B + self._a * A
        # H = D - C
        #
        # X3 = E * F
        # Y3 = G * H
        # T3 = E * H
        # Z3 = F * G
        #endregion

        A = X1 * X2
        B = Y1 * Y2
        C = T1 * self._d * T2
        D = Z1

        E = (X1 + Y1) * (X2 + Y2) - A - B
        F = D - C
        G = D + C
        H = B - self._a * A

        X3 = E * F
        Y3 = G * H
        T3 = E * H
        Z3 = F * G

        return ExtendedCoord(X3, Y3, T3, Z3, domain=self.domain)

    def _add_with_z_ne(self, X1, Y1, Z1, T1, X2, Y2, Z2, T2):
        """ Add points on twisted edwards curve (general case)"""
        #region the core idea:
        # Unified addition formula:
        # X3 = (X1 * Y2 + Y1 * X2) * (Z1 * Z2 - self._d * T1 * T2)
        # Y3 = (Y1 * Y2 - self._a * X1 * X2) * (Z1 * Z2 + self._d * T1 * T2)
        # T3 = (Y1 * Y2 - self._a * X1 * X2) * (X1 * Y2 + Y1 * X2)
        # Z3 = (Z1 * Z2 - self._d * T1 * T2) * (Z1 * Z2 + self._d * T1 * T2)

        # Dedicated addition formula: (NOTE, this fails for (x,y)+(-x, -y), but why?
        # Calculate the components for the result point
        # A = (X1 * Y2 - Y1 * X2)
        # B = (T1 * Z2 + Z1 * T2)
        # C = (Y1 * Y2 + self._a * X1 * X2)
        # D = (T1 * Z2 - Z1 * T2)
        #
        # X3 = A * B
        # Y3 = C * D
        # T3 = B * D
        # Z3 = C * A
        #endregion

        #region failed dedicated additioner
        # A = X1 * X2
        # B = Y1 * Y2
        # C = Z1 * T2
        # D = T1 * Z2
        #
        # E = D + C
        # F = (X1 - Y1) * (X2 + Y2) + B - A
        # G = B + self._a * A
        # H = D - C
        #
        # X3 = E * F
        # Y3 = G * H
        # T3 = E * H
        # Z3 = F * G
        #endregion

        # Calculations
        A = X1 * X2
        B = Y1 * Y2
        C = self._d * T1 * T2
        D = Z1 * Z2
        E = (X1 + Y1) * (X2 + Y2) - A - B

        F = D - C
        G = D + C
        H = B - self._a * A

        X3 = E * F
        Y3 = G * H
        T3 = E * H
        Z3 = F * G

        return ExtendedCoord(X3, Y3, T3, Z3, domain=self.domain)

    def add_points(self, p1, p2):

        if p1.is_identity_point():
            return p2
        if p2.is_identity_point():
            return p1

        X1, Y1, Z1, T1 = p1.X, p1.Y, p1.Z, p1.T
        X2, Y2, Z2, T2 = p2.X, p2.Y, p2.Z, p2.T


        # Optimizations: for one of Z is 1, use 8M optimized unified addition version
        if Z2 == 1:
            # if Z1 == Z2:  # 7M version that fails... doesn't worth the efforts
            #     return self._add_with_z_eq_1(X1, Y1, T1,
            #                                  X2, Y2, T2)

            return self._add_with_z2_1(X1, Y1, Z1, T1,
                                       X2, Y2, T2)

        if Z1 == 1:
            return self._add_with_z2_1(X2, Y2, Z2, T2,
                                       X1, Y1, T1)

        # 9M+2D unified addition formula, handles all cases
        return self._add_with_z_ne(X1, Y1, Z1, T1,
                                   X2, Y2, Z2, T2)

    def add_points_affine(self, p1: List[int], p2: List[int]) -> tuple[int, int]:
        '''
            Add two points in affine coordinate, and return the result point in affine coordinate
            * We skip the `normal` way of calculating point addition on affine plan, directly apply jacobian addition
        '''
        P1 = ExtendedCoord.from_affine(p1, domain=self.domain)
        P2 = ExtendedCoord.from_affine(p2, domain=self.domain)
        R = self.add_points(P1, P2)
        return R.get_affine_coords()

    def double_point_affine(self, p: List[int]) -> tuple[int, int]:

        P = ExtendedCoord.from_affine(p, domain=self.domain)
        return self.double_point(P).get_affine_coords()

    def _double_with_z_1(self, X1: GFElementType, Y1: GFElementType) -> ExtendedCoord:
        """Add a point to itself with z == 1."""
        A = X1 ** 2
        B = Y1 ** 2
        D = self._a * A
        E = (X1 + Y1) ** 2 - A - B

        G = D + B
        H = D - B

        X3 = E * (G - self.domain(2))
        Y3 = G * H
        T3 = E * H
        Z3 = G ** 2 - self.domain(2) * G
        return ExtendedCoord(X3, Y3, T3, Z3, self.domain)

    def double_point(self, p: ExtendedCoord) -> ExtendedCoord:
        '''
            Double p on the curve in extended coord.
            ref: [2] 3.3
        '''
        X1, Y1, Z1, T1 = p.X, p.Y, p.Z, p.T

        if Z1 == self.domain(1):
            return self._double_with_z_1(X1, Y1)

        # Perform the doubling using explicit formulas
        A = X1 ** 2
        B = Y1 ** 2
        C = self.domain(2) * Z1 ** 2
        D = self._a * A
        E = (X1 + Y1) ** 2 - A - B

        G = D + B
        F = G - C
        H = D - B

        X3 = E * F
        Y3 = G * H
        T3 = E * H
        Z3 = F * G

        return ExtendedCoord(X3, Y3, T3, Z3, self.domain)
