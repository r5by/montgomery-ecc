import re
from mecc.short_weierstrass_ecc import ShortWeierstrassCurve
from mecc.twisted_edwards_ecc import TwistedEdwardsCurve
from modp import GFmodp
from mecc.coordinate import JacobianCoord, ExtendedCoord
from sage.all import GF, EllipticCurve
from tests.sage_helper import CustomSageTwistedEdwardsCurve


def remove_whitespace(text):
    """Removes all whitespace from passed in string"""
    return re.sub(r"\s+", "", text, flags=re.UNICODE)


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

