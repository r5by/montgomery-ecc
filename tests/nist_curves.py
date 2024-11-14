import re
from mecc.short_weierstrass_ecc import ShortWeierstrassCurve
from mecc.twisted_edwards_ecc import TwistedEdwardsCurve
from modp import GFmodp
from mecc.coordinate import JacobianCoord, ExtendedCoord
from sage.all import GF, EllipticCurve
from tests.sage_helper import CustomSageTwistedEdwardsCurve

# todo> remove
# _p = 29
# a, b = 2, 20
# x, y = 18, 28
# k = 129
# SAGE_F = GF(_p)
# SAGE_curve = EllipticCurve(SAGE_F, [a, b])
# SAGE_g = SAGE_curve(x, y)
# exp = k * SAGE_g
# exp = (int(exp[0]), int(exp[1]), int(exp[2]))
#
# F_256 = GFmodp(_p)
# p256 = ShortWeierstrassCurve(domain=F_256, a=a, b=b)
# G256_affine = [x, y]
# G256 = JacobianCoord.from_affine(G256_affine, F_256)
#
# # act = p256.k_point_fixed(k, w=4, P=G256, k_max_bits=64)
# act = p256.k_point_fixed(k, w=4, P=G256)
# act.to_affine()
# act = act.get_integer_coords()
#
# assert exp == act, f'Scalar multiplication on Curve: {p256} fails at k={k}'

# curve parameters are taken from:
# https://github.com/ijrsvt/python-ecdsa/blob/master/src/ecdsa/ecdsa.py
def remove_whitespace(text):
    """Removes all whitespace from passed in string"""
    return re.sub(r"\s+", "", text, flags=re.UNICODE)

#region NIST Curve P-521:
_p = int(
    "686479766013060971498190079908139321726943530014330540939"
    "446345918554318339765605212255964066145455497729631139148"
    "0858037121987999716643812574028291115057151"
)
_r = int(
    "686479766013060971498190079908139321726943530014330540939"
    "446345918554318339765539424505774633321719753296399637136"
    "3321113864768612440380340372808892707005449"
)
# s = 0xd09e8800291cb85396cc6717393284aaa0da64baL
# c = int(remove_whitespace(
#    """
#         0b4 8bfa5f42 0a349495 39d2bdfc 264eeeeb 077688e4
#    4fbf0ad8 f6d0edb3 7bd6b533 28100051 8e19f1b9 ffbe0fe9
#    ed8a3c22 00b8f875 e523868c 70c1e5bf 55bad637"""
# ), 16)
_b = int(
    remove_whitespace(
        """
         051 953EB961 8E1C9A1F 929A21A0 B68540EE A2DA725B
    99B315F3 B8B48991 8EF109E1 56193951 EC7E937B 1652C0BD
    3BB1BF07 3573DF88 3D2C34F1 EF451FD4 6B503F00"""
    ),
    16,
)
_Gx = int(
    remove_whitespace(
        """
          C6 858E06B7 0404E9CD 9E3ECB66 2395B442 9C648139
    053FB521 F828AF60 6B4D3DBA A14B5E77 EFE75928 FE1DC127
    A2FFA8DE 3348B3C1 856A429B F97E7E31 C2E5BD66"""
    ),
    16,
)
_Gy = int(
    remove_whitespace(
        """
         118 39296A78 9A3BC004 5C8A5FB4 2C7D1BD9 98F54449
    579B4468 17AFBD17 273E662C 97EE7299 5EF42640 C550B901
    3FAD0761 353C7086 A272C240 88BE9476 9FD16650"""
    ),
    16,
)

SAGE_F521 = GF(_p)
SAGE_p521 = EllipticCurve(SAGE_F521, [-3, _b])
SAGE_G521 = SAGE_p521(_Gx, _Gy)

from mont.montgomery import Montgomery
USE_MONT = True # toggle this to enable/disable montgomery

M8 = Montgomery.factory(mod=_p, mul_opt='real8').build(m=64, w=256)
M8 = M8.config(mul_opt='real0')
F_521 = M8 if USE_MONT else GFmodp(_p)
F_521.modulus = _p

p521 = ShortWeierstrassCurve(domain=F_521, a=-3, b=_b)
G521_affine = [_Gx, _Gy]
G521 = JacobianCoord.from_affine(G521_affine, F_521)
#endregion

#region NIST Curve P-256:
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
#endregion

#region NIST 3.2.3.2 Edwards448 (a.k.a the "Goldilocks")
# edwards448, defined in RFC7748
_p = 2**448 - 2**224 - 1
_a = 1
_d = -39081 % _p
# _h = 4

_Gx = int(
    remove_whitespace(
        "224580040295924300187604334099896036246789641632564134246125461"
        "686950415467406032909029192869357953282578032075146446173674602635"
        "247710"
    )
)
_Gy = int(
    remove_whitespace(
        "298819210078481492676017930443930673437544040154080242095928241"
        "372331506189835876003536878655418784733982303233503462500531545062"
        "832660"
    )
)
# _r = 2**446 - 0x8335DC163BB124B65129C96FDE933D8D723A70AADC873D6D54A7BB0D

F_448 = GFmodp(_p)
edwards448 = TwistedEdwardsCurve(domain=F_448, a=_a, d=_d)
G448_affine = [_Gx, _Gy]
G448 = ExtendedCoord.from_affine(G448_affine, F_448)

SAGE_F448 = GF(_p)
SAGE_edwards448 = CustomSageTwistedEdwardsCurve(SAGE_F448, [_a, _d])
#endregion

#region NIST 3.2.3.1 Edwards25519
# edwards25519, defined in RFC7748
_p = 2**255 - 19
_a = -1
_d = int(
    remove_whitespace(
        "370957059346694393431380835087545651895421138798432190163887855330"
        "85940283555"
    )
)
# _h = 8

_Gx = int(
    remove_whitespace(
        "151122213495354007725011514095885315114540126930418572060461132"
        "83949847762202"
    )
)
_Gy = int(
    remove_whitespace(
        "463168356949264781694283940034751631413079938662562256157830336"
        "03165251855960"
    )
)
# _r = 2**252 + 0x14DEF9DEA2F79CD65812631A5CF5D3ED

F_25519 = GFmodp(_p)
edwards25519 = TwistedEdwardsCurve(domain=F_25519, a=_a, d=_d)
G25519_affine = [_Gx, _Gy]
G25519 = ExtendedCoord.from_affine(G25519_affine, F_25519)

SAGE_F25519 = GF(_p)
SAGE_edwards25519 = CustomSageTwistedEdwardsCurve(SAGE_F25519, [_a, _d])

#endregion

