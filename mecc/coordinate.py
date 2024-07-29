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


class ProjectiveCoord:
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


class ExtendedCoord(ProjectiveCoord):

    def __init__(self,
                 x: Union[int, GFElementType],
                 y: Union[int, GFElementType],
                 z: Union[int, GFElementType] = 1,
                 t: Union[int, GFElementType] = None,
                 domain: Optional[GFType] = None):
        super().__init__(x, y, z, domain)
        self.T = t

    def __neg__(self):
        """Return the additive inverse of the point.
        for twisted edward curve, the negation of point (X,Y,Z,T)
        is conventionally chosen to be (-X,Y,Z,-T)
        """
        if self.is_point_at_infinity():
            return self

        neg_x = self.domain(-self.X) if isinstance(self.X, int) else -self.X
        neg_t = self.domain(-self.T) if isinstance(self.T, int) else -self.T
        return ExtendedCoord(neg_x, self.Y, self.Z, neg_t, self.domain)

    def __eq__(self, other):

        if not isinstance(other, ExtendedCoord):
            return NotImplemented

        if self.domain != other.domain:
            raise ValueError("Cannot compare points from different domains")

        if self.is_point_at_infinity() and other.is_point_at_infinity():
            return True
        if self.is_point_at_infinity() or other.is_point_at_infinity():
            return False

        # Use cross multiplication to compare without needing division
        # (X1/Z1 = X2/Z2) --> X1*Z2 = X2*Z1
        # (Y1/Z1 = Y2/Z2) --> Y1*Z2 = Y2*Z1
        norm_x1 = self.X * other.Z
        norm_x2 = other.X * self.Z
        norm_y1 = self.Y * other.Z
        norm_y2 = other.Y * self.Z

        return norm_x1 == norm_x2 and norm_y1 == norm_y2

    def get_integer_coords(self) -> Tuple[int, int, int, int]:
        return int(self.X), int(self.Y), int(self.Z), int(self.T)

    def __repr__(self):
        return f"ExtendedCoord(X={self.X}, Y={self.Y}, Z={self.Z}, T={self.T})"

    def enter_domain(self, domain: GFType):
        self.domain = domain

        self.X = domain(self.X)
        self.Y = domain(self.Y)
        self.Z = domain(self.Z)
        self.T = domain(self.T)

    def is_point_at_infinity(self):
        """
        Check if the point is the point at infinity in extended coordinates.
        The point at infinity is typically [0 : 1 : 0 : 0], but we only check
        X, Z, and T because Y becomes irrelevant when Z is zero.
        """
        return self.X == 0 and self.Z == 0 and self.T == 0

    def is_affine(self):
        """
        Check if the point is in affine coordinates.
        A point is in affine coordinates if Z == 1.
        """
        # Check if the Z-coordinate is 1, which indicates affine coordinates
        return self.Z == 1

    @classmethod
    def from_int_coord(cls, x: int, y: int, z: int = 1, t: Optional[int] = None, *, domain: GFType) -> \
            'ExtendedCoord':
        p = cls(x, y, z, t, domain) if t else cls(x, y, z, 0, domain)
        p.enter_domain(domain)

        # calc. T
        T = p.X * p.Y / p.Z
        if t and p.T != T:
            raise ValueError(f'Given T-coordinate:{t} does not match with given projective coordinate: (x={x},y={y},'
                             f'z={z})')

        p.T = T
        return p

    @classmethod
    def point_at_infinity(cls, domain: GFType):
        # (0: 1: 0: 0) for extended vs. (0: 1: 0) for projective, note the difference.
        return cls.from_int_coord(0, 1, 0, 0, domain=domain)

    @classmethod
    def from_affine(cls, point: Union[List[int], PointAtInfinity], domain: GFType) -> 'ExtendedCoord':
        '''
            Transfer affine coordinate [x, y] into ext [X:Y:1:T=x*y] (using domain arithmetics if domain is specified)
            usage:
                $ INF = InfinitePoint()
                $ mont = Montgomery.factory(mod=29).build()
                $ p1 = ExtendedCoord.from_affine([1, 2], mont)
                # p2 = ExtendedCoord.from_affine(INF)
        '''
        if isinstance(point, PointAtInfinity):
            return cls.point_at_infinity(domain)

        if len(point) != 2:
            raise ValueError("Affine coordinates must be a list of two elements [x, y].")
        x, y = point
        return cls.from_int_coord(x, y, 1, domain=domain)  # Z=1 for affine coordinates

    @classmethod
    def copy(cls, x: GFElementType, y: GFElementType, z: GFElementType, t: GFElementType, domain: GFType) -> 'ExtendedCoord':
        return ExtendedCoord(x, y, z, t, domain)  # a deep copy

    def get_affine_coords(self) -> Union[Tuple[int, int], PointAtInfinity]:
        if self.is_point_at_infinity():
            return PointAtInfinity()

        if not self.is_affine():
            self._to_affine()  # todo>

        return int(self.X), int(self.Y)

    def to_affine(self):

        if self.is_point_at_infinity() or self.is_affine():
            return

        self._to_affine()

    def _to_affine(self):
        ''' Internal: project this point onto affine plane '''
        # todo> what if Z == 0??
        X, Y, Z = self.X, self.Y, self.Z
        _Z = 1 / Z

        x = X * _Z
        y = Y * _Z

        self.X = x
        self.Y = y
        self.Z = self.domain(1)

        self.T = x * y


class JacobianCoord(ProjectiveCoord):

    def __neg__(self):
        """Return the additive inverse of the point."""
        if self.is_point_at_infinity():
            return self

        neg_y = self.domain(-self.Y) if isinstance(self.Y, int) else -self.Y
        return JacobianCoord(self.X, neg_y, self.Z, self.domain)

    def __eq__(self, other):
        if not isinstance(other, JacobianCoord):
            return NotImplemented

        if self.domain != other.domain:
            raise ValueError("Cannot compare points from different domains")

        if self.is_point_at_infinity() and other.is_point_at_infinity():
            return True
        if self.is_point_at_infinity() or other.is_point_at_infinity():
            return False

        # Cross-multiplication to avoid division
        # X1/Z1^2 == X2/Z2^2  --->  X1*Z2^2 == X2*Z1^2
        left_x = self.X * other.Z ** 2
        right_x = other.X * self.Z ** 2
        # Y1/Z1^3 == Y2/Z2^3  --->  Y1*Z2^3 == Y2*Z1^3
        left_y = self.Y * other.Z ** 3
        right_y = other.Y * self.Z ** 3

        return left_x == right_x and left_y == right_y

    def get_integer_coords(self) -> Tuple[int, int, int]:
        return int(self.X), int(self.Y), int(self.Z)

    def __repr__(self):
        return f"JacobianCoord(X={self.X}, Y={self.Y}, Z={self.Z})"

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
