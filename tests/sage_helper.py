from sage.all import GF, EllipticCurve as SageEllipticCurve


def sage_generate_points_on_curve(curve_params, field_size):
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

    # Create the elliptic curve over the field
    F = GF(field_size)
    curve = SageEllipticCurve(F, curve_params)

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

    # 2) generate double points for testing point doubling
    double_points = [proj2affine[2 * p] for p in points]

    # 3) generate addition table for testing point addition
    addition_table = generate_flat_addition_table(points, proj2affine)

    return affine_points, projective_points, double_points, addition_table


def proj_to_affine(p):
    x = int(p[0] / p[2])
    y = int(p[1] / p[2])
    return (x, y)


def generate_flat_addition_table(points, proj2affine):
    size = len(points)
    flat_table = []

    for i in range(size):
        for j in range(i + 1):  # Only compute up to the diagonal
            p = proj2affine[points[i] + points[j]]
            flat_table.append(p)

    return flat_table


def add_lut(flat_table, i, j):
    # Calculate the index in the flat table
    if j > i:
        i, j = j, i
    index = (i * (i + 1)) // 2 + j
    return flat_table[index]

