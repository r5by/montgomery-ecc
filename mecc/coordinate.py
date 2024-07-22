from mont.typing import GFElementType, GFType
from typing import Union, List, Optional, Tuple


class PointAtInfinity:
    '''Singleton repr. of the infinite point'''
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(PointAtInfinity, cls).__new__(cls)
        return cls._instance

    def __repr__(self):
        return "InfinitePoint"


class JacobianCoord:
    def __init__(self,
                 x: Union[int, GFElementType],
                 y: Union[int, GFElementType],
                 z: Union[int, GFElementType] = 1,
                 domain: Optional[GFType] = None):

        self.domain = domain
        # Verify that all instances are of the same type
        self.verify_instances(x, y, z)

        # (X, Y, Z) keeps the input (x, y) in the Z/pZ domain by default
        self.X = x
        self.Y = y
        self.Z = z

    def __neg__(self):
        """Return the additive inverse of the point."""
        if self.is_point_at_infinity():
            return self

        neg_y = self.domain(-self.Y) if isinstance(self.Y, int) else -self.Y
        return JacobianCoord(self.X, neg_y, self.Z, self.domain)

    def get_integer_coords(self) -> Tuple[int, int, int]:
        return int(self.X), int(self.Y), int(self.Z)

    def __repr__(self):
        return f"JacobianCoord(X={self.X}, Y={self.Y}, Z={self.Z})"

    def verify_instances(self, x, y, z):
        """
        Ensure that x, y, and z are of the same type.
        This is important to maintain consistency in the type of the coordinates.
        """
        # Get the types of x, y, z
        types = {type(x), type(y), type(z)}
        # Check if all types are the same or not
        if len(types) != 1:
            raise TypeError("All input coordinates must be of the same type.")
        # Further, check if they align with the expected types based on the domain
        if self.domain is None and not all(isinstance(v, GFElementType) for v in [x, y, z]):
            raise TypeError("All input coordinates must be instances of GFElementType when no domain is specified.")
        if all(isinstance(v, int) for v in [x, y, z]) and self.domain is None:
            raise TypeError("Arithmetic domain must be specified if all input coordinates are integers.")

    def enter_domain(self, domain: GFType):
        self.domain = domain

        self.X = domain(self.X)
        self.Y = domain(self.Y)
        self.Z = domain(self.Z)

    def is_point_at_infinity(self):
        return self.Z == 0

    def is_affine(self):
        return self.Z == 1

    @classmethod
    def point_at_infinity(cls, domain: GFType):
        # Note (1: 1: 0) for jacobian and (0: 1: 0) for projective, note the difference.
        return cls.from_int_coord(1, 1, 0, domain)

    @classmethod
    def from_affine(cls, point: Union[List[int], PointAtInfinity], domain: GFType) -> 'JacobianCoord':
        '''
            Transfer affine coordinate [x, y] into jacobian [X:Y:1] (using domain arithmetics if domain is specified)
            usage:
                $ INF = InfinitePoint()
                $ mont = Montgomery.factory(mod=29).build()
                $ p1 = JacobianCoord.from_affine([1, 2], mont)
                # p2 = JacobianCoord.from_affine(INF)
        '''
        if isinstance(point, PointAtInfinity):
            return cls.point_at_infinity(domain)

        if len(point) != 2:
            raise ValueError("Affine coordinates must be a list of two elements [x, y].")
        x, y = point
        return cls.from_int_coord(x, y, 1, domain)  # Z=1 for affine coordinates

    @classmethod
    def copy(cls, x: GFElementType, y: GFElementType, z: GFElementType, domain: GFType) -> 'JacobianCoord':
        return JacobianCoord(x, y, z, domain)

    @classmethod
    def from_int_coord(cls, x: int, y: int, z: int, domain: GFType) -> 'JacobianCoord':
        p = cls(x, y, z, domain)
        p.enter_domain(domain)
        return p

    def get_affine_coords(self) -> Union[Tuple[int, int], PointAtInfinity]:
        if self.is_point_at_infinity():
            return PointAtInfinity()

        if not self.is_affine():
            # todo> apply xgcd?
            self._to_affine()

        return int(self.X), int(self.Y)

    def to_affine(self):

        if self.is_point_at_infinity() or self.is_affine():
            return

        self._to_affine()

    def _to_affine(self):
        ''' Internal: project this point onto affine plane '''
        # Cost: 1I(inverse) + 1S(square) + 3M(multiplication)
        X, Y, Z = self.X, self.Y, self.Z
        _Z = 1 / Z
        _Z2 = _Z * _Z
        _Z3 = _Z * _Z2

        x = X * _Z2
        y = Y * _Z3

        self.X = x
        self.Y = y
        self.Z = self.domain(1)
