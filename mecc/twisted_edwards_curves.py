from abc import ABC

from mecc.interface import EllipticCurve
from typing import Union, List
from mecc.coordinate import JacobianCoord, PointAtInfinity
from mont.typing import GFType


const_scala = lambda c, x: sum([x for _ in range(c)])  # c*x for c as a constant
double = lambda x: x + x


class TwistedEdwards(EllipticCurve, ABC):
    """
        A twisted edwards curve has the following form:
            ax^2 + y^2 = 1 + dx^2y^2  (affine)
            Y^2 = X^3 + aXZ^4 + bZ^6 (Jacobian)
        where a, d are coefficients in the domain

        SP 800-186 specified:
            Edwards25519:
                p = 2^255 - 19, a = -1, d = -121665/121666

            Edwards448:
                p = 2^448 - 2^224 - 1, a = 1, d = -39081
    """