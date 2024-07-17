import re
from mecc.short_weierstrass_ecc import ShortWeierstrassCurve
from modp import GFmodp
from mecc.coordinate import JacobianCoord
from sage.all import GF, EllipticCurve


def remove_whitespace(text):
    """Removes all whitespace from passed in string"""
    return re.sub(r"\s+", "", text, flags=re.UNICODE)


# NIST Curve P-256:
_p = int(
    remove_whitespace(
        """
    1157920892103562487626974469494075735300861434152903141955
    33631308867097853951"""
    )
)
_r = int(
    remove_whitespace(
        """
    115792089210356248762697446949407573529996955224135760342
    422259061068512044369"""
    )
)
# s = 0xc49d360886e704936a6678e1139d26b7819f7e90L
# c = 0x7efba1662985be9403cb055c75d4f7e0ce8d84a9c5114abcaf3177680104fa0dL
_b = int(
    remove_whitespace(
        """
    5AC635D8 AA3A93E7 B3EBBD55 769886BC 651D06B0 CC53B0F6
    3BCE3C3E 27D2604B"""
    ),
    16,
)
_Gx = int(
    remove_whitespace(
        """
    6B17D1F2 E12C4247 F8BCE6E5 63A440F2 77037D81 2DEB33A0
    F4A13945 D898C296"""
    ),
    16,
)
_Gy = int(
    remove_whitespace(
        """
    4FE342E2 FE1A7F9B 8EE7EB4A 7C0F9E16 2BCE3357 6B315ECE
    CBB64068 37BF51F5"""
    ),
    16,
)


SAGE_F256 = GF(_p)
SAGE_p256 = EllipticCurve(SAGE_F256, [-3, _b])
SAGE_G256 = SAGE_p256(_Gx, _Gy)

F_256 = GFmodp(_p)
p256 = ShortWeierstrassCurve(domain=F_256, a=-3, b=_b)
G256_affine = [_Gx, _Gy]
G256 = JacobianCoord.from_affine(G256_affine, F_256)

