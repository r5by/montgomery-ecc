from dataclasses import dataclass

# Curve parameters for Elliptic Curve defined by y^2 = x^3 + 2*x + 20 over Finite Field of size 29
a = 2
b = 20
p = 29

@dataclass
class ProjectiveCoord:
    x: int
    y: int
    z: int

# Precomputed Look-Up Table (LUT)
precomputed = [
    ProjectiveCoord(1 % p, 1 % p, 0 % p),
    ProjectiveCoord(3 % p, 13 % p, 1 % p),
    ProjectiveCoord(3 % p, 13 % p, 1 % p),
    ProjectiveCoord(23 % p, 16 % p, 1 % p),
    ProjectiveCoord(3 % p, 13 % p, 1 % p),
    ProjectiveCoord(23 % p, 16 % p, 1 % p),
    ProjectiveCoord(23 % p, 16 % p, 1 % p),
    ProjectiveCoord(16 % p, 1 % p, 1 % p),
    ProjectiveCoord(3 % p, 13 % p, 1 % p),
    ProjectiveCoord(23 % p, 16 % p, 1 % p),
    ProjectiveCoord(23 % p, 16 % p, 1 % p),
    ProjectiveCoord(16 % p, 1 % p, 1 % p),
    ProjectiveCoord(23 % p, 16 % p, 1 % p),
    ProjectiveCoord(16 % p, 1 % p, 1 % p),
    ProjectiveCoord(16 % p, 1 % p, 1 % p),
    ProjectiveCoord(16 % p, 28 % p, 1 % p)
]

def mod_p(x):
    return (x % p + p) % p  # Handle negative values as well

def mod_p_inv(x):
    return pow(x, p - 2, p)

def double_point(p):
    X = p.x
    Y = p.y
    Z = p.z

    S = mod_p(4 * X * mod_p(Y * Y))
    Z2 = mod_p(Z * Z)
    Z4 = mod_p(Z2 * Z2)
    M = mod_p(3 * mod_p(X * X) + a * Z4)

    X_ = mod_p(mod_p(M * M) - mod_p(S + S))
    Y2 = mod_p(Y * Y)
    Y4 = mod_p(Y2 * Y2)
    Y_ = mod_p(mod_p(M * (S - X_)) - mod_p(8 * Y4))
    Z_ = mod_p(2 * Y * Z)

    return ProjectiveCoord(X_, Y_, Z_)

def add_point(p1, p2):
    X1 = p1.x
    Y1 = p1.y
    Z1 = p1.z
    X2 = p2.x
    Y2 = p2.y
    Z2 = p2.z

    if Z1 == 0:
        return p2

    if Z2 == 0:
        return p1

    if X1 == X2 and Y1 == Y2 and Z1 == Z2:
        return double_point(p1)

    Z1Z1 = mod_p(Z1 * Z1)
    Z2Z2 = mod_p(Z2 * Z2)
    U1 = mod_p(X1 * Z2Z2)
    U2 = mod_p(X2 * Z1Z1)
    S1 = mod_p(Y1 * Z2 * Z2Z2)
    S2 = mod_p(Y2 * Z1 * Z1Z1)

    H = mod_p(U2 - U1)
    _2H = mod_p(H + H)
    I = mod_p(_2H * _2H)

    J = mod_p(H * I)
    r = mod_p(2 * (S2 - S1))

    V = mod_p(U1 * I)
    _2V = mod_p(V + V)

    _r2 = mod_p(r * r)
    X3 = mod_p(_r2 - J - _2V)
    Y3 = mod_p(r * (V - X3) - 2 * mod_p(S1 * J))
    Z3 = mod_p((mod_p((Z1 + Z2) * (Z1 + Z2)) - Z1Z1 - Z2Z2) * H)

    return ProjectiveCoord(X3, Y3, Z3)

def int_by_slider(n, d, i):
    index = 0
    r = 0

    n >>= i  # Right shift n by i positions

    while n > 0:
        b = n & 1  # Get the least significant bit
        r |= (b << index)  # Set the corresponding bit in r

        n >>= d  # Right shift n by d positions
        index += 1  # Move to the next bit position in r

    return r

def test_slider(n, d, i, expected):
    result = int_by_slider(n, d, i)
    if result == expected:
        print(f"Test passed for n={n}, d={d}, i={i}: result={result}")
    else:
        print(f"Test failed for n={n}, d={d}, i={i}: expected={expected}, got={result}")

def k_point_fixed(k, px, py, pz):
    q = ProjectiveCoord(1, 1, 0)
    d = 3

    for i in range(d - 1, -1, -1):
        q = double_point(q)
        ki = int_by_slider(k, d, i)
        pre = precomputed[ki]
        q = add_point(q, pre)

    return q

def to_affine(X, Y, Z):
    if Z == 0:
        return ProjectiveCoord(1, 1, 0)

    _Z = mod_p_inv(Z)
    _Z2 = mod_p(_Z * _Z)
    _Z3 = mod_p(_Z * _Z2)

    x = mod_p(X * _Z2)
    y = mod_p(Y * _Z3)

    return ProjectiveCoord(x, y, 1)

def k_point_fixed_inlined(k, px, py, pz):
    q = ProjectiveCoord(1, 1, 0)
    d = 3  # load from HW

    for i in range(d - 1, -1, -1):
        # print(f'loop i={i}\n')
        # Inlined double_point function
        X = q.x
        Y = q.y
        Z = q.z

        _Y2 = mod_p(Y * Y)
        S = mod_p(4 * X * _Y2)

        _X2 = mod_p(X * X)  # X square
        _Z2 = mod_p(Z * Z)  # Z square
        _Z4 = mod_p(_Z2 * _Z2)  # Z**4
        M = mod_p(3 * _X2 + a * _Z4)

        _M2 = mod_p(M * M)  # M square
        X_ = mod_p(_M2 - mod_p(S + S))

        _Y4 = mod_p(_Y2 * _Y2)
        Y_ = mod_p(mod_p(M * (S - X_)) - mod_p(8 * _Y4))
        Z_ = mod_p(2 * Y * Z)

        q.x = X_
        q.y = Y_
        q.z = Z_

        # Inlined int_by_slider function
        n = k
        n >>= i
        index = 0
        r = 0
        while n > 0:
            b = n & 1  # Get the least significant bit
            r |= (b << index)  # Set the corresponding bit in r
            n >>= d  # Right shift n by d positions
            index += 1  # Move to the next bit position in r
        ki = r

        # Use precomputed points based on the value of ki
        pre = precomputed[ki]  # load from HW

        # Inlined add_point function
        X1 = q.x
        Y1 = q.y
        Z1 = q.z
        X2 = pre.x
        Y2 = pre.y
        Z2 = pre.z

        if Z1 == 0:
            q = ProjectiveCoord(X2, Y2, Z2)
            continue

        if Z2 == 0:
            continue

        if X1 == X2 and Y1 == Y2 and Z1 == Z2:
            # Inlined double_point for q = double_point(q)

            _Y2 = mod_p(Y * Y)
            S = mod_p(4 * X * _Y2)

            _X2 = mod_p(X * X)  # X square
            _Z2 = mod_p(Z * Z)  # Z square
            _Z4 = mod_p(_Z2 * _Z2)  # Z**4
            M = mod_p(3 * _X2 + a * _Z4)

            _M2 = mod_p(M * M)  # M square
            X_ = mod_p(_M2 - mod_p(S + S))

            _Y4 = mod_p(_Y2 * _Y2)
            Y_ = mod_p(mod_p(M * (S - X_)) - mod_p(8 * _Y4))
            Z_ = mod_p(2 * Y * Z)

            q.x = X_
            q.y = Y_
            q.z = Z_
            continue

        Z1Z1 = mod_p(Z1 * Z1)  # Z1^2
        Z2Z2 = mod_p(Z2 * Z2)  # Z2^2
        U1 = mod_p(X1 * Z2Z2)
        U2 = mod_p(X2 * Z1Z1)
        S1 = mod_p(Y1 * Z2 * Z2Z2)
        S2 = mod_p(Y2 * Z1 * Z1Z1)

        H = mod_p(U2 - U1)
        _2H = mod_p(H + H)
        I = mod_p(_2H * _2H)

        J = mod_p(H * I)
        r = mod_p(2 * (S2 - S1))

        V = mod_p(U1 * I)
        _2V = mod_p(V + V)

        _r2 = mod_p(r * r)
        X3 = mod_p(_r2 - J - _2V)
        Y3 = mod_p(r * (V - X3) - 2 * mod_p(S1 * J))
        Z3 = mod_p((mod_p((Z1 + Z2) * (Z1 + Z2)) - Z1Z1 - Z2Z2) * H)

        q.x = X3
        q.y = Y3
        q.z = Z3

    return q

def main():
    # Test k point
    k = 169
    Q = ProjectiveCoord(3, 13, 1)
    kQ = k_point_fixed(k, Q.x, Q.y, Q.z)
    kQ2 = k_point_fixed_inlined(k, Q.x, Q.y, Q.z)
    kQ = to_affine(kQ.x, kQ.y, kQ.z)
    kQ2 = to_affine(kQ2.x, kQ2.y, kQ2.z)
    print(f"k_point of Q is: X={kQ.x}, Y={kQ.y}, Z={kQ.z}")
    print(f"k_point of Q (inlined) is: X={kQ2.x}, Y={kQ2.y}, Z={kQ2.z}")

if __name__ == "__main__":
    main()
