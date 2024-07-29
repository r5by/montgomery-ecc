from abc import ABC

from mecc.interface import EllipticCurve
from typing import Union, List
from mecc.coordinate import ExtendedCoord, PointAtInfinity
from mont.typing import GFType, GFElementType

const_scala = lambda c, x: sum([x for _ in range(c)])  # c*x for c as a constant
double = lambda x: x + x


class TwistedEdwards(EllipticCurve, ABC):
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

        References:
            [1] NIST SP 800-186
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

    def is_point_on_curve(self, point: Union[List[int], ExtendedCoord, PointAtInfinity]):
        """
            Verify if given point is on curve, a "point" format is like
            1) [int, int]:  affine coord.
            2) [int, int, int, int]: extended coord.
            3) ExtendedCoord: already in extended coord.
        """
        if isinstance(point, ExtendedCoord):
            if point.is_point_at_infinity():
                return True

            if point.domain is None:
                point.enter_domain(self.domain)

            return self.verify_point_on_curve(point)

        if isinstance(point, PointAtInfinity):
            return True

        if len(point) == 2:
            point = ExtendedCoord.from_affine(point, self.domain)
        elif len(point) == 3:
            X, Y, Z = point
            point = ExtendedCoord.from_int_coord(X, Y, Z, self.domain)
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

    def add_points(self, p1, p2):
        """Add two points on the curve."""
        pass

    def double_point(self, p):
        """Double a point on the curve."""
        pass

    def k_point(self, k, p):
        """Perform scalar multiplication of a point by k."""
        pass
