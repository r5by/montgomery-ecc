from mont.typing import GFElementType, GFType
from typing import Union, List, Optional, Tuple
from abc import ABC, abstractmethod


class PointAtInfinity:
    '''Singleton repr. of the infinite point'''
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(PointAtInfinity, cls).__new__(cls)
        return cls._instance

    def __repr__(self):
        return "InfinitePoint"


class ProjectiveCoord(ABC):
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

    def is_affine(self):
        """Check if the point is in affine form."""
        return self.Z == 1

    def to_affine(self):

        if self.is_affine():
            return

        self._to_affine()

    @abstractmethod
    def __neg__(self):
        pass

    @abstractmethod
    def get_affine_coords(self):
        """Get the affine coordinates of the point in raw (integer) type"""
        pass

    @abstractmethod
    def is_identity_point(self):
        """Method to determine if the point is the identity element of the curve group."""
        pass

    @abstractmethod
    def _to_affine(self):
        """Convert point to affine coordinates."""
        pass

    @classmethod
    @abstractmethod
    def get_identity_point(cls, domain: GFType):
        pass


class ExtendedCoord(ProjectiveCoord):
    '''Refer to <<Twisted Edwards Curves Revisited>> Section 3.

        [X : Y: T: Z] as extended of projective coord [X: Y: Z] for T=XY/Z
    '''

    def __init__(self,
                 x: Union[int, GFElementType],
                 y: Union[int, GFElementType],
                 t: Union[int, GFElementType],
                 z: Union[int, GFElementType] = 1,
                 domain: Optional[GFType] = None):
        super().__init__(x, y, z, domain)
        self.T = t

    def __neg__(self):
        """Return the additive inverse of the point.
        for twisted edward curve, the negation of point (X,Y,Z,T)
        is conventionally chosen to be (-X,Y,Z,-T)
        """
        if self.is_identity_point():
            return self

        neg_x = self.domain(-self.X) if isinstance(self.X, int) else -self.X
        neg_t = self.domain(-self.T) if isinstance(self.T, int) else -self.T
        return ExtendedCoord(neg_x, self.Y, neg_t, self.Z, self.domain)

    def __eq__(self, other):

        if not isinstance(other, ExtendedCoord):
            return NotImplemented

        if self.domain != other.domain:
            raise ValueError("Cannot compare points from different domains")

        if self.is_identity_point() and other.is_identity_point():
            return True
        if self.is_identity_point() or other.is_identity_point():
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

    def is_identity_point(self):
        """
        Check if the point is the identity point in extended coordinates.
        The point at infinity is typically [0 : 1 : 0: 1] (extended of (0:1:1))
        """
        zero = self.domain(0)
        # Check the conditions for the identity point
        return self.X == zero and self.Y == self.Z and self.T == zero

    @classmethod
    def from_int_coord(cls, x: int, y: int, *, t: Optional[int] = None, z: int = 1, domain: GFType) -> \
            'ExtendedCoord':
        p = cls(x, y, t, z, domain) if t else cls(x, y, 0, z, domain)
        p.enter_domain(domain)

        # calc. T
        T = p.X * p.Y / p.Z
        if t and p.T != T:
            raise ValueError(f'Given T-coordinate:{t} does not match with given projective coordinate: (x={x},y={y},'
                             f'z={z})')

        p.T = T
        return p

    @classmethod
    def get_identity_point(cls, domain: GFType):
        # (0: 1: 1: 0) for extended coordinates of identity point (0,1) on affine plane of twisted edwards curves
        return cls.from_int_coord(0, 1, t=0, z=1, domain=domain)

    @classmethod
    def from_affine(cls, point: List[int], domain: GFType) -> 'ExtendedCoord':
        '''
            Transfer affine coordinate [x, y] into ext [X:Y: T=x*y: Z=1] (using domain arithmetics if domain is
            specified)
            usage:
                $ mont = Montgomery.factory(mod=29).build()
                $ p1 = ExtendedCoord.from_affine([1, 2], mont)
        '''
        if len(point) != 2:
            raise ValueError("Affine coordinates must be a list of two elements [x, y].")
        x, y = point
        return cls.from_int_coord(x=x, y=y, z=1, domain=domain)  # Z=1 for affine coordinates

    def get_affine_coords(self) -> Tuple[int, int]:

        if not self.is_affine():
            self._to_affine()

        return int(self.X), int(self.Y)

    def to_affine(self):

        if self.is_affine():
            return

        self._to_affine()

    def _to_affine(self):
        ''' Internal: project this point onto affine plane '''
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
        if self.is_identity_point():
            return self

        neg_y = self.domain(-self.Y) if isinstance(self.Y, int) else -self.Y
        return JacobianCoord(self.X, neg_y, self.Z, self.domain)

    def __eq__(self, other):
        if not isinstance(other, JacobianCoord):
            return NotImplemented

        if self.domain != other.domain:
            raise ValueError("Cannot compare points from different domains")

        if self.is_identity_point() and other.is_identity_point():
            return True
        if self.is_identity_point() or other.is_identity_point():
            return False

        # Cross-multiplication to avoid division
        # X1/Z1^2 == X2/Z2^2  --->  X1*Z2^2 == X2*Z1^2
        Z2 = other.Z
        Z2Z2 = Z2 * Z2
        # left_x = self.X * other.Z ** 2
        left_x = self.X * Z2Z2

        Z1 = self.Z
        Z1Z1 = Z1 * Z1
        # right_x = other.X * self.Z ** 2
        right_x = other.X * Z1Z1

        # Y1/Z1^3 == Y2/Z2^3  --->  Y1*Z2^3 == Y2*Z1^3
        Z2Z2Z2 = Z2Z2 * Z2
        # left_y = self.Y * other.Z ** 3
        left_y = self.Y * Z2Z2Z2

        Z1Z1Z1 = Z1Z1 * Z1
        # right_y = other.Y * self.Z ** 3
        right_y = other.Y * Z1Z1Z1

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

    def is_identity_point(self):
        return self.Z == 0

    def get_affine_coords(self):
        if self.is_identity_point():
            return PointAtInfinity()

        if not self.is_affine():
            self._to_affine()

        return self.X, self.Y

    @classmethod
    def get_identity_point(cls, domain: GFType):
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
            return cls.get_identity_point(domain)

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

    def _to_affine(self):
        ''' Internal: project this point onto affine plane '''
        if self.is_identity_point():
            # clear x, and y if z == 0 (identity)
            self.X = self.domain(1)
            self.Y = self.domain(1)
            return

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
