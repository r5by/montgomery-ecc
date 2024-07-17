from mont.montgomery import Montgomery
from typing import Union, List, Optional


class EllipticCurve:
    def __init__(self, domain: Union[int, Montgomery], *coeffs: Union[List[int], int, None], **kwargs):
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

        # Montgomery domain associated with this elliptic curve
        if isinstance(domain, Montgomery):
            self.mont = domain
        elif isinstance(domain, int):
            # default Montgomery, use a software-friendly one
            self.mont = Montgomery.factory(mod=domain, mul_opt='real1').build()
        else:
            raise TypeError("Montgomery domain or its moduli must be specified")

        # Coefficients in Integer domain
        self.a = a
        self.b = b
        # Repr. of a, b in underlying Montgomery domain
        self._a = self.mont(a)
        self._b = self.mont(b)

    @staticmethod
    def factory(domain: Union[int, Montgomery], a: int, b: int) -> 'EllipticCurve':
        # todo> use this for coordinate change?
        if isinstance(domain, int):
            # Create Montgomery domain first
            mont_domain = Montgomery.factory(domain)
            return EllipticCurve(mont_domain, a=a, b=b)
        elif isinstance(domain, Montgomery):
            return EllipticCurve(domain, a=a, b=b)
        else:
            raise TypeError("Domain must be either an integer or a Montgomery object")

    def __str__(self):
        return f"Elliptic Curve over Montgomery domain: {self.mont} with a={self.a}, b={self.b}"


# if __name__ == '__main__':
#     p = 61
#     # M = Montgomery.factory(mod=p, mul_opt='real3').build(w=8)
#
#     M = Montgomery.factory(mod=p).build(m=2, w=5)
#     x = 46
#
#     _x = M(x)
#
#     print('done')
