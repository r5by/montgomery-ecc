from abc import ABC, abstractmethod


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

    @abstractmethod
    def k_point(self, k, p):
        """Perform scalar multiplication of a point by k."""
        pass
