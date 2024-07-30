from abc import ABC, abstractmethod
from mecc.coordinate import ProjectiveCoord
from mecc.utils import naf_prodinger, int_by_slider, msb, turn_off_msb
from typing import Optional, List
import math


class EllipticCurve(ABC):
    """
        An interface for Elliptic Curves, including:
        * ShortWeierstrassCurve
        * MontgomeryCurve
        * TwistedEdwardCurve
    """

    @abstractmethod
    def is_point_on_curve(self, point):
        """Verify if a point is on the curve."""
        pass

    @abstractmethod
    def add_points(self, p1, p2):
        """Add two points on the curve."""
        pass

    @abstractmethod
    def double_point(self, p):
        """Double a point on the curve."""
        pass

    def k_point(self, k: int, P: ProjectiveCoord) -> ProjectiveCoord:
        ''' ECSM (elliptic curve scalar multiplication) of k*P for P is unfixed point
                use NAF(prodinger) to minimize the point add/sub operations
         '''

        if k == 0 or P.is_identity_point():
            return type(P).get_identity_point(self.domain)

        if k == 1:
            return P

        np, nm = naf_prodinger(k)

        # Project point p onto affine to safe the loop cost on addition and subtraction
        P.to_affine()
        Q = type(P).get_identity_point(self.domain)

        # Determine maximum bit length of np or nm to determine loop range
        max_bit_length = max(np.bit_length(), nm.bit_length())

        for i in range(max_bit_length - 1, -1, -1):
            Q = self.double_point(Q)
            if (np >> i) & 1:
                Q = self.add_points(P, Q)
            if (nm >> i) & 1:
                Q = self.add_points(-P, Q)

        return Q

    def _precompute_win(self, w: int, d: int, P: ProjectiveCoord) -> List[ProjectiveCoord]:
        '''
            Precompute 2^{wi}*P for i in [0 .. d-1] ({d} in total LUT entries)
                * is THE precomputation for windowing method
                * is Part of the precomputation for Comb method
        '''
        res = [P]  # length=d LUT

        for _ in range(1, d):
            nex_point = self.k_point(1 << w, res[-1])
            res.append(nex_point)

        return res

    def _precompute_comb(self, w: int, d: int, P: ProjectiveCoord) -> List[ProjectiveCoord]:
        ''' Precompute all combinations of w-bits repr of k with radix 2^d multiply the point P '''
        # NOTE: This precomputation depends on input (t) to generate (d) from (w), even though w is known!
        INF = type(P).get_identity_point(self.domain)
        res = [INF for _ in range(1 << w)]  # len(Lut) = 2^w (can optimized the first two elements [INF, P] out)

        # step 1) calculate power_dp[i] saves 2^(id) * P for i in [0 .. w-1]
        power_dp = self._precompute_win(d, w, P)

        # [_.to_affine() for _ in power_dp]  # debug usage

        # step 2) populate the precomputations suing DP
        for i in range(1, 1 << w):
            j, k = msb(i), turn_off_msb(i)
            tmp = self.add_points(res[k], power_dp[j])
            tmp.to_affine()
            res[i] = tmp

        return res

    def k_point_fixed(self, k: int, w: int, P: ProjectiveCoord, k_max_bits: Optional[int] = None) -> ProjectiveCoord:
        ''' Default fixed-point scalar multiplication using the Comb method
            ref [1] Algorithm 3.44
        '''
        # print(f'Starting scalar mul for k(={k})*P(={P})')
        if k == 0:
            return type(P).get_identity_point(self.domain)

        # Prepare to perform the multiplication
        Q = type(P).get_identity_point(self.domain)

        # simulate the maximum LUT (precomputed&load per curve)
        t = k_max_bits if k_max_bits else k.bit_length()
        d = math.ceil(t / w)
        precomputed = self._precompute_comb(w, d, P)  # NOTE: should load this LUT from HW impl. in the real practice

        for i in range(d - 1, -1, -1):
            Q = self.double_point(Q)
            ki = int_by_slider(k, d, i)
            Q = self.add_points(Q, precomputed[ki])

        return Q
