from sage.all import GF, polygens, Curve, EllipticCurve as SageEllipticCurve
from typing import List


class CustomSageTwistedEdwardsCurve:
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
        res.append((0, 1, 0))  # don't forget the point at infinity
        return res

    def double_affine(self, x1, y1):
        # Apply the doubling formula for twisted Edwards curves
        a, d = self.a, self.d
        x3 = (x1 * y1 + y1 * x1) / (1 + d * x1 * x1 * y1 * y1)
        y3 = (y1 * y1 - a * x1 * x1) / (1 - d * x1 * x1 * y1 * y1)
        return (int(x3), int(y3))

    def add_affine(self, x1, y1, x2, y2):
        # Apply the addition formula for twisted Edwards curves
        a, d = self.a, self.d
        x3 = (x1 * y2 + y1 * x2) / (1 + d * x1 * x2 * y1 * y2)
        y3 = (y1 * y2 - a * x1 * x2) / (1 - d * x1 * x2 * y1 * y2)
        return (int(x3), int(y3))

    def __neg__(self, x1, y1):
        # Apply the negation formula for twisted Edwards curves
        return (-x1, y1)


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

    # 2) generate double points and addition table for testing point doubling
    if type == 'weierstrass':
        double_points = [proj2affine[2 * p] for p in points]
        addition_table = generate_flat_addition_table_weierstrass(points, proj2affine)
    elif type == 'edwards':
        double_points = [(curve.double_affine(x[0], x[1]) if x != 'INF' else 'INF') for x in affine_points]
        addition_table = generate_flat_addition_table_edwards(affine_points, curve)

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


def generate_flat_addition_table_edwards(affine_points, curve: CustomSageTwistedEdwardsCurve):
    size = len(affine_points)
    flat_table = []

    for i in range(size):
        for j in range(i + 1):  # Only compute up to the diagonal
            pX, pY = affine_points[i], affine_points[j]

            if pX == 'INF' and pY == 'INF':
                p = 'INF'

            if pX == 'INF':
                p = pY
            elif pY == 'INF':
                p = pX
            else:
                p = curve.add_affine(pX[0], pX[1], pY[0], pY[1])

            flat_table.append(p)

    return flat_table


def add_lut(flat_table, i, j):
    # Calculate the index in the flat table
    if j > i:
        i, j = j, i
    index = (i * (i + 1)) // 2 + j
    return flat_table[index]
