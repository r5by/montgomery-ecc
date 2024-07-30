from sage.all import GF, polygens, Curve, EllipticCurve as SageEllipticCurve
from typing import List


class CustomSageTwistedEdwardsCurve:
    # NOTE: Twisted Edwards Curves do NOT need a point at infinity!!
    def __init__(self, field, coeffs):
        self.field = field
        self.a = field(coeffs[0])
        self.d = field(coeffs[1])
        # Define the curve equation
        x, y = polygens(field, 'x,y')
        self.curve_eq = self.a * x ** 2 + y ** 2 - 1 - self.d * x ** 2 * y ** 2
        self.x = x
        self.y = y

    def is_on_curve(self, x_val, y_val):
        # Evaluate the curve's equation at (x, y)
        return self.curve_eq(x=x_val, y=y_val) == 0

    def points(self):
        # List all points in projective coord. on the curve (on affine plane)
        res = [(x_val, y_val, self.field(1)) for x_val in self.field for y_val in self.field
               if self.is_on_curve(x_val, y_val)]
        return res

    def double(self, x1, y1):
        if x1 == 0:   # short circuits for (0,1) and (0, -1)
            return 0, 1

        # Apply the doubling formula for twisted Edwards curves
        a, d = self.a, self.d
        x3 = (x1 * y1 + y1 * x1) / (self.field(1) + d * x1 * x1 * y1 * y1)
        y3 = (y1 * y1 - a * x1 * x1) / (self.field(1) - d * x1 * x1 * y1 * y1)
        return int(x3), int(y3)

    def add(self, x1, y1, x2, y2):
        if x1 == -x2 and y1 == y2:  # short circuits of inverse points
            return 0, 1

        a, d = self.a, self.d
        x3 = (x1 * y2 + y1 * x2) / (self.field(1) + d * x1 * x2 * y1 * y2)
        y3 = (y1 * y2 - a * x1 * x2) / (self.field(1) - d * x1 * x2 * y1 * y2)
        return int(x3), int(y3)


def sage_generate_points_on_curve(curve_params: List[int], field_size: int, type: str):
    """
    Generate a list of affine coordinates on an elliptic curve over a finite field.

    Args:
    curve_params (list or tuple): Coefficients of the curve (a, b) for y^2 = x^3 + ax + b.
    field_size (int): The size of the field (prime number).
    num_points (int): The number of points to generate

    Returns:
    list: A list of tuples representing the affine coordinates of points on the curve.
    """
    if field_size > (1 << 16):
        raise ValueError(f'Filed size is too large for this little toy test-tool to handle...')

    F = GF(field_size)
    if type == 'weierstrass':
        # Create the elliptic curve over the field
        curve = SageEllipticCurve(F, curve_params)
    elif type == 'edwards':
        curve = CustomSageTwistedEdwardsCurve(F, curve_params)
    else:
        raise ValueError(f'Undefined elliptic curve type: {type}')

    # Generate points (sage returns projective points by default)
    points = curve.points()
    cardinality = len(points)  # cardinality of the group over the ecc including point at infinity
    projective_points = [[int(point[0]), int(point[1]), int(point[2])] for point in points]

    if type == 'weierstrass':
        # map: proj -> affine
        proj2affine = {}

        # Filter out the point at infinity and convert to affine coordinates
        affine_points = []
        for p in points:

            if p[2] == 0:
                proj2affine[p] = 'INF'
            else:
                # if len(affine_points) >= num_points:
                #     break

                proj2affine[p] = proj_to_affine(p)

            affine_points.append(proj2affine[p])

        # 2) generate double points and addition table for testing point double
        double_points = [proj2affine[2 * p] for p in points]
        addition_table = generate_flat_addition_table_weierstrass(points, proj2affine)
    elif type == 'edwards':
        affine_points = [(int(point[0]), int(point[1])) for point in points]
        double_points = [curve.double(x[0], x[1]) for x in points]
        addition_table = generate_flat_addition_table_edwards(points, curve)

    return affine_points, projective_points, double_points, addition_table


def proj_to_affine(p):
    x = int(p[0] / p[2])
    y = int(p[1] / p[2])
    return (x, y)


def generate_flat_addition_table_weierstrass(points, proj2affine):
    size = len(points)
    flat_table = []

    for i in range(size):
        for j in range(i + 1):  # Only compute up to the diagonal
            p = proj2affine[points[i] + points[j]]
            flat_table.append(p)

    return flat_table


def generate_flat_addition_table_edwards(points, curve: CustomSageTwistedEdwardsCurve):
    size = len(points)
    flat_table = []

    for i in range(size):
        for j in range(i + 1):  # Only compute up to the diagonal
            pX, pY = points[i], points[j]
            p = curve.add(pX[0], pX[1], pY[0], pY[1])
            flat_table.append(p)

    return flat_table


def add_lut(flat_table, i, j):
    # Calculate the index in the flat table
    if j > i:
        i, j = j, i
    index = (i * (i + 1)) // 2 + j
    return flat_table[index]
