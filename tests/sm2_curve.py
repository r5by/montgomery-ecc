from mecc.short_weierstrass_ecc import ShortWeierstrassCurve
from mecc.coordinate import JacobianCoord
from sage.all import GF, EllipticCurve
from mont.montgomery import Montgomery
from tests.modp import GFmodp


USE_MONT = True
# USE_MONT = False

#region SM2
_p = 0xfffffffeffffffffffffffffffffffffffffffff00000000ffffffffffffffff
_a = 0xfffffffeffffffffffffffffffffffffffffffff00000000fffffffffffffffc
_b = 0x28e9fa9e9d9f5e344d5a9e4bcf6509a7f39789f515ab8f92ddbcbd414d940e93
_Gx = 0x32c4ae2c1f1981195f9904466a39c9948fe30bbff2660be1715a4589334c74c7
_Gy = 0xbc3736a2f4f6779c59bdcee36b692153d0a9877cc62a474002df32e52139f0a0

SAGE_SM2_F = GF(_p)
SAGE_SM2_c = EllipticCurve(SAGE_SM2_F, [_a, _b])
SAGE_SM2_G = SAGE_SM2_c(_Gx, _Gy)

M8 = Montgomery.factory(mod=_p, mul_opt='real8').build(m=64, w=256)
M8 = M8.config(mul_opt='real0')
SM2_F = M8 if USE_MONT else GFmodp(_p)

SM2_c = ShortWeierstrassCurve(domain=SM2_F, a=_a, b=_b)
SM2_G_affine = [_Gx, _Gy]
SM2_G = JacobianCoord.from_affine(SM2_G_affine, SM2_F)
MOD = _p
#endregion