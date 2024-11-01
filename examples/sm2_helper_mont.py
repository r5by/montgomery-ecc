from dataclasses import dataclass
from typing import List
from mecc.coordinate import JacobianCoord
from tests.sm2_curve import SM2_c, M8
from mecc.utils import fe_double, fe_const_scala

@dataclass
class SimpleProjCoord:
    x: int
    y: int
    z: int


# Curve parameters for SM2 ecc (in Mont domain)
p = 0xfffffffeffffffffffffffffffffffffffffffff00000000ffffffffffffffff

R = 2135987035920910082395021706169552114602704522356652769947041607822219725780640550022962086936576
_R = p - 188719626594729257751182451177044373699592459328094631476166851035125 # inverse of R
ONE = 26959946667150639796128516724350533591840829255256855500763837759489 # identity in mont domain

def _enter_domain(x: int) -> int:
    return x * R % p

def _exit_domain(x: int) -> int:
    return x * _R % p

# Montgomery multiplication
def MM(x, y):
    return (x * y * _R) % p

def MA(x, y):
    return (x + y) % p

def MS(x, y):
    return (x - y + p) % p

def MI(x: int):
    # inverse of x (in mont domain), satisfying: MM(inv_x, x) == ONE
    inv_x = M8._mont_inv(x)
    return inv_x


def point_exit_domain(point: SimpleProjCoord) -> SimpleProjCoord:
    return SimpleProjCoord(_exit_domain(point.x), _exit_domain(point.y), _exit_domain(point.z))

def to_affine(X: int, Y: int, Z: int) -> SimpleProjCoord:
    if Z == 0:
        return SimpleProjCoord(ONE, ONE, 0)

    _Z = MI(Z)
    _Z2 = MM(_Z, _Z)
    _Z3 = MM(_Z, _Z2)

    x = MM(X, _Z2)
    y = MM(Y, _Z3)

    return SimpleProjCoord(x, y, ONE)

# Curve parameters for SM2 ecc (over F_p)
_a = 0xfffffffeffffffffffffffffffffffffffffffff00000000fffffffffffffffc
_b = 0x28e9fa9e9d9f5e344d5a9e4bcf6509a7f39789f515ab8f92ddbcbd414d940e93
_Gx = 0x32c4ae2c1f1981195f9904466a39c9948fe30bbff2660be1715a4589334c74c7
_Gy = 0xbc3736a2f4f6779c59bdcee36b692153d0a9877cc62a474002df32e52139f0a0

# Load directly to save computation time
# a = _enter_domain(_a)
# b = _enter_domain(_b)
# Gx = _enter_domain(_Gx)
# Gy = _enter_domain(_Gy)
a = 115792089129476408754968425825635342593198753216401703688650627430998171713532
b = 22011861304025731791618878281070852302073234238990620437919231873049343083601
Gx = 33233751390280102669221384524433571339618901719018408594327814856354129491301
Gy = 51990514252569618031863371609527714278242372626802480674841694412900887993234

Gz = ONE
G = SimpleProjCoord(Gx, Gy, Gz) # base point (jacobian) over Montgomery domain

# Precomputed Look-Up Table (LUT) over montgomery domain
precomputed = [(26959946667150639796128516724350533591840829255256855500763837759489, 26959946667150639796128516724350533591840829255256855500763837759489), (33233751390280102669221384524433571339618901719018408594327814856354129491301, 51990514252569618031863371609527714278242372626802480674841694412900887993234), (36127976750173100903146405188765072942664794279968034595565245596706538631479, 16705954087254396076902395459485677519669634752816383431790635355072791026668), (113849852132070319862032403644217365649283713692366459803897028398633871935818, 89541874020243318879585961177230406817953682588886276904881336595194473260414), (61716983751213724254805231848712638613965673217358400838514465327279684820119, 6976980861615193897487328512143457069275727942954854828645272006983248852622), (112656745305372105887137422766706510415903008726284720275150745351373671734926, 97397822655754523875087786712097437270809706217298162047369834675827575477254), (69569175545051660252052483323074883473785361374982311437798602966529654934138, 37669341935196491112778075595021623915178172941339129494820344352480352174985), (105304122577785964684104609456907027886845313661457016607430982438923559448517, 100832269873530012789801144852368693594382870666351890108035981916345370455317), (43579333093403194840423304599181589915045473934913743744935574705007367471365, 105871265818870584607418339055614735641180378121436897527384582410143585539447), (16378933207214915342621136989715264587533298127026884545183881295898725879926, 75799883401937505152350606106936566506370616865149805780583874638405944138559), (59257525484990649564764034143230279823075657188593484357850979426291742764311, 103987251943312548389025119439365031889188770070954701873801830792439115359454), (93480391416923884533870500617211914431932426500261380249698110374940135160170, 28571273271495383148246712397686239275986914604035625988540019722175696247614), (46990075202872426116033196926967141657653317873018185663999878737306891471954, 2822525920645221648884689121985478750186665977012106059359828591379813835665), (28376707764488840387732819579971342710367919085647365419404143060101522844204, 6068052079241673593671652754709808091722623384824593298232364682212355594009), (87776473165800642513455083540391424318653662862792057785068327010467133241700, 40227679123314942194593183984066892088198639125266306073226301844358546348296), (58558914835608549663193756114592611134411871656620681949698435001208694136132, 74591717052139306186986153847722950594127039222898775495353115384919275233138)]


def _scale(c, x):
    r = fe_const_scala(c, x)
    return r % p

def _double(x):
    return MA(x, x)

# inline-free double point
def _double_point(X, Y, Z):
    YY = MM(Y, Y)
    XYY = MM(X, YY)
    S = _scale(4, XYY)

    Z2 = MM(Z, Z)
    Z4 = MM(Z2, Z2)
    XX = MM(X, X)
    _3XX = _scale(3, XX)
    _aZ4 = MM(a, Z4)
    M = MA(_3XX, _aZ4)

    _M2 = MM(M, M)
    _2S = MA(S, S)
    X_ = MS(_M2, _2S)

    Y2 = MM(Y, Y)
    Y4 = MM(Y2, Y2)
    _SX_ = MS(S, X_)
    _MSX = MM(M, _SX_)
    _8Y4 = _scale(8, Y4)
    Y_ = MS(_MSX, _8Y4)

    YZ = MM(Y, Z)
    Z_ = _double(YZ)
    return X_, Y_, Z_


def double_point(p):
    X_, Y_, Z_ = _double_point(p.x, p.y, p.z)
    return SimpleProjCoord(X_, Y_, Z_)

# def test_double_point():
#     G_j = JacobianCoord.from_affine([_Gx, _Gy], M8)
#     t_2G = SM2_c.double_point(G_j)
#     t_2G.to_affine()
#     exp = t_2G.get_integer_coords()
#
#     act = double_point(G)
#     act = to_affine(act.x, act.y, act.z)
#     act1 = point_exit_domain(act)
#
#     assert exp == act1

# inline-free point add
def _add_point(X1, Y1, Z1, X2, Y2, Z2):

    if Z1 == 0:
        return X2, Y2, Z2

    if Z2 == 0:
        return X1, Y1, Z1

    if X1 == X2 and Y1 == Y2 and Z1 == Z2:
        return _double_point(X1, Y1, Z1)

    Z1Z1 = MM(Z1, Z1)
    Z2Z2 = MM(Z2, Z2)
    U1 = MM(X1, Z2Z2)
    U2 = MM(X2, Z1Z1)

    Z2Z2Z2 = MM(Z2, Z2Z2)
    S1 = MM(Y1, Z2Z2Z2)

    Z1Z1Z1 = MM(Z1, Z1Z1)
    S2 = MM(Y2, Z1Z1Z1)

    H = MS(U2, U1)
    _2H = MA(H, H)
    I = MM(_2H, _2H)

    J = MM(H, I)
    _S2_S1 = MS(S2, S1)
    r = _double(_S2_S1)

    V = MM(U1, I)
    _2V = _double(V)
    rr = MM(r, r)
    _rr_j = MS(rr, J)
    X3 = MS(_rr_j, _2V)

    S1J = MM(S1, J)
    _2S1J = _double(S1J)
    _V_X3 = MS(V, X3)
    rV_X3 = MM(r, _V_X3)
    Y3 = MS(rV_X3, _2S1J)

    Z1_Z2 = MA(Z1, Z2)
    _Z1_Z1_2 = MM(Z1_Z2, Z1_Z2)
    t1 = MS(_Z1_Z1_2, Z1Z1)
    t2 = MS(t1, Z2Z2)
    Z3 = MM(t2, H)
    return X3, Y3, Z3

def add_point(p1, p2):

    X3, Y3, Z3 = _add_point(p1.x, p1.y, p1.z,
                            p2.x, p2.y, p2.z)

    return SimpleProjCoord(X3, Y3, Z3)

# def test_add_point():
#     G_j = JacobianCoord.from_affine([_Gx, _Gy], M8)
#     t_2G = SM2_c.double_point(G_j)
#     t_2G.to_affine()
#     exp = t_2G.get_integer_coords()
#
#     act = add_point(G, G)
#     act = to_affine(act.x, act.y, act.z)
#     act1 = point_exit_domain(act)
#
#     assert exp == act1


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

def k_point_fixed(k, px, py, pz):
    q = SimpleProjCoord(ONE, ONE, 0)
    d = 64

    # for i in range(d - 1, -1, -1):
    for i in range(d):
        q = double_point(q)
        ki = int_by_slider(k, d, d - i - 1)

        if ki == 0:
            continue

        pre_x, pre_y = precomputed[ki]
        pre = SimpleProjCoord(pre_x, pre_y, ONE) # note, only when ki > 0, Z <-- 1

        q = add_point(q, pre)

    return q

# def k_point_fixed_inlined(k, px, py, pz):
#     q = SimpleProjCoord(1, 1, 0)
#     d = 64  # load from HW
#
#     for i in range(d - 1, -1, -1):
#         # print(f'loop i={i}\n')
#         # Inlined double_point function
#         X = q.x
#         Y = q.y
#         Z = q.z
#
#         _Y2 = MM(Y * Y)
#         S = MM(4 * X * _Y2)
#
#         _X2 = MM(X * X)  # X square
#         _Z2 = MM(Z * Z)  # Z square
#         _Z4 = MM(_Z2 * _Z2)  # Z**4
#         M = MM(3 * _X2 + a * _Z4)
#
#         _M2 = MM(M * M)  # M square
#         X_ = MM(_M2 - MM(S + S))
#
#         _Y4 = MM(_Y2 * _Y2)
#         Y_ = MM(MM(M * (S - X_)) - MM(8 * _Y4))
#         Z_ = MM(2 * Y * Z)
#
#         q.x = X_
#         q.y = Y_
#         q.z = Z_
#
#         # Inlined int_by_slider function
#         n = k >> i
#         index = 0
#         ki = 0
#         tmp_n = n
#         while tmp_n > 0:
#             b = tmp_n & 1  # Get the least significant bit
#             ki |= (b << index)  # Set the corresponding bit in r
#             tmp_n >>= d  # Right shift n by d positions
#             index += 1  # Move to the next bit position in r
#
#         if ki == 0:
#             continue
#
#         # Use precomputed points based on the value of ki
#         new_x, new_y = precomputed[ki]  # load from HW
#         pre = SimpleProjCoord(new_x, new_y, 1)
#
#         # Inlined add_point function
#         X1 = q.x
#         Y1 = q.y
#         Z1 = q.z
#         X2 = pre.x
#         Y2 = pre.y
#         Z2 = pre.z
#
#         if Z1 == 0:
#             q = SimpleProjCoord(X2, Y2, Z2)
#             continue
#
#         if Z2 == 0:
#             continue
#
#         if X1 == X2 and Y1 == Y2 and Z1 == Z2:
#             # Inlined double_point for q = double_point(q)
#
#             _Y2 = MM(Y * Y)
#             S = MM(4 * X * _Y2)
#
#             _X2 = MM(X * X)  # X square
#             _Z2 = MM(Z * Z)  # Z square
#             _Z4 = MM(_Z2 * _Z2)  # Z**4
#             M = MM(3 * _X2 + a * _Z4)
#
#             _M2 = MM(M * M)  # M square
#             X_ = MM(_M2 - MM(S + S))
#
#             _Y4 = MM(_Y2 * _Y2)
#             Y_ = MM(MM(M * (S - X_)) - MM(8 * _Y4))
#             Z_ = MM(2 * Y * Z)
#
#             q.x = X_
#             q.y = Y_
#             q.z = Z_
#             continue
#
#         Z1Z1 = MM(Z1 * Z1)  # Z1^2
#         Z2Z2 = MM(Z2 * Z2)  # Z2^2
#         U1 = MM(X1 * Z2Z2)
#         U2 = MM(X2 * Z1Z1)
#         S1 = MM(Y1 * Z2 * Z2Z2)
#         S2 = MM(Y2 * Z1 * Z1Z1)
#
#         H = MM(U2 - U1)
#         _2H = MM(H + H)
#         I = MM(_2H * _2H)
#
#         J = MM(H * I)
#         r = MM(2 * (S2 - S1))
#
#         V = MM(U1 * I)
#         _2V = MM(V + V)
#
#         _r2 = MM(r * r)
#         X3 = MM(_r2 - J - _2V)
#         Y3 = MM(r * (V - X3) - 2 * MM(S1 * J))
#         Z3 = MM((MM((Z1 + Z2) * (Z1 + Z2)) - Z1Z1 - Z2Z2) * H)
#
#         q.x = X3
#         q.y = Y3
#         q.z = Z3
#
#     return q

def naf_prodinger(x):
    # https://en.wikipedia.org/wiki/Non-adjacent_form
    xh = x >> 1
    x3 = x + xh
    c = xh ^ x3
    np = x3 & c
    nm = xh & c
    return np, nm

def k_point(k: int, P: SimpleProjCoord) -> SimpleProjCoord:
    ''' ECSM (elliptic curve scalar multiplication) of k*P for P is unfixed point
            use NAF(prodinger) to minimize the point add/sub operations
     '''

    if k == 0 or (P.x == ONE and P.y == ONE and P.z == 0):
        return SimpleProjCoord(ONE, ONE, 0)

    if k == 1:
        return P

    np, nm = naf_prodinger(k)

    # Project point p onto affine to safe the loop cost on addition and subtraction
    P = to_affine(P.x, P.y, P.z)
    Q = SimpleProjCoord(ONE, ONE, 0)

    # Determine maximum bit length of np or nm to determine loop range
    max_bit_length = max(np.bit_length(), nm.bit_length())

    for i in range(max_bit_length - 1, -1, -1):

        Q = double_point(Q)

        # 1) original approach
        # if (np >> i) & 1:
        #     Q = add_point(P, Q)
        #
        # if (nm >> i) & 1:
        #     neg_P = ProjectiveCoord(P.x, p - P.y, P.z)
        #     Q = add_point(neg_P, Q)

        # 2) optimized for assembly, not f_np && f_nm is impossible for a valid NAF
        f_np = (np >> i) & 1
        f_nm = (nm >> i) & 1

        py = P.y
        if f_nm:
            py = p - P.y
        # elif f_np:
        #     py = P.y

        f_np_or_nm = f_nm or f_np
        if f_np_or_nm:
            Q.x, Q.y, Q.z = _add_point(P.x, py, P.z, Q.x, Q.y, Q.z)

    return Q

# def k_point_inlined(k: int, P: SimpleProjCoord) -> SimpleProjCoord:
#
#     if k == 0 or (P.x == 1 and P.y == 1 and P.z == 0):
#         return SimpleProjCoord(1, 1, 0)
#
#     if k == 1:
#         return P
#
#     # inline: np, nm = naf_prodinger(k)
#     xh = k >> 1
#     x3 = k + xh
#     c = xh ^ x3
#     np = x3 & c
#     nm = xh & c
#
#     # Project point p onto affine to safe the loop cost on addition and subtraction
#     # P = to_affine(P.x, P.y, P.z)
#     ## Inline replace P ==> PX, PY, PZ:
#     PX = P.x
#     PY = P.y
#     PZ = P.z
#
#     if PZ == 0:
#         PX = 1
#         PY = 1
#     else:
#         _Z = MI(PZ)
#         _Z2 = MM(_Z * _Z)
#         _Z3 = MM(_Z * _Z2)
#
#         PX = MM(PX * _Z2)
#         PY = MM(PY * _Z3)
#
#     # Q = ProjectiveCoord(1, 1, 0)
#     # Inline replace ==>
#     QX = 1
#     QY = 1
#     QZ = 0
#
#     # Determine maximum bit length of np or nm to determine loop range
#     # max_bit_length = max(np.bit_length(), nm.bit_length())
#     # Inline replace (for SM2) ==>
#     max_bit_length = 256
#
#     for i in range(max_bit_length - 1, -1, -1):
#         # Q = double_point(Q)
#         # inline replace ==>
#         # QX, QY, QZ = _double_point(QX, QY, QZ)
#         # inline replace ==>
#         S = MM(4 * QX * MM(QY * QY))
#         Z2 = MM(QZ * QZ)
#         Z4 = MM(Z2 * Z2)
#         M = MM(3 * MM(QX * QX) + a * Z4)
#
#         QX = MM(MM(M * M) - MM(S + S))
#         Y2 = MM(QY * QY)
#         QZ = MM(2 * QY * QZ)
#         Y4 = MM(Y2 * Y2)
#         QY = MM(MM(M * (S - QX)) - MM(8 * Y4))
#
#         f_np = (np >> i) & 1
#         f_nm = (nm >> i) & 1
#
#         py = PY
#         if f_nm:
#             py = p - PY
#
#         f_np_or_nm = f_nm or f_np
#         if f_np_or_nm:
#             # Q.x, Q.y, Q.z = _add_point(P.x, py, P.z, Q.x, Q.y, Q.z)
#             # inline replace ==>
#             X1 = PX
#             Y1 = py
#             Z1 = PZ
#             X2 = QX
#             Y2 = QY
#             Z2 = QZ
#
#             if Z1 == 0:
#                 QX = X2
#                 QY = Y2
#                 QZ = Z2
#
#             elif Z2 == 0:
#                 QX = X1
#                 QY = Y1
#                 QZ = Z1
#
#             elif X1 == X2 and Y1 == Y2 and Z1 == Z2:
#                 QX, QY, QZ = _double_point(X1, Y1, Z1)
#
#             else:
#                 Z1Z1 = MM(Z1 * Z1)
#                 Z2Z2 = MM(Z2 * Z2)
#                 U1 = MM(X1 * Z2Z2)
#                 U2 = MM(X2 * Z1Z1)
#                 S1 = MM(Y1 * Z2 * Z2Z2)
#                 S2 = MM(Y2 * Z1 * Z1Z1)
#
#                 H = MM(U2 - U1)
#                 _2H = MM(H + H)
#                 I = MM(_2H * _2H)
#
#                 J = MM(H * I)
#                 r = MM(2 * (S2 - S1))
#
#                 V = MM(U1 * I)
#                 _2V = MM(V + V)
#                 _r2 = MM(r * r)
#
#                 X3 = MM(_r2 - J - _2V)
#                 Y3 = MM(r * (V - X3) - 2 * MM(S1 * J))
#                 Z3 = MM((MM((Z1 + Z2) * (Z1 + Z2)) - Z1Z1 - Z2Z2) * H)
#
#                 QX = X3
#                 QY = Y3
#                 QZ = Z3
#
#     return SimpleProjCoord(QX, QY, QZ)


def FPM_test():
    golden = "sm2_golden.data"
    with open(golden, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            tokens = line.split(',')
            if len(tokens) != 3:
                print(f"Warning: Expected 3 tokens, got {len(tokens)} in line: {line}")
                continue
            k_str, var_x_str, var_y_str = tokens
            k = int(k_str.strip())
            exp_x = int(var_x_str.strip())
            exp_y = int(var_y_str.strip())

            # Test k point
            kG = k_point_fixed(k, G.x, G.y, G.z)
            kG = to_affine(kG.x, kG.y, kG.z)
            kG = point_exit_domain(kG)
            # print(f"k_point of G is: X={kG.x}, Y={kG.y}, Z={kG.z}")
            assert exp_x == kG.x and exp_y == kG.y


            # kG2 = k_point_fixed_inlined(k, G.x, G.y, G.z)
            # kG2 = to_affine(kG2.x, kG2.y, kG2.z)
            # # print(f"k_point of G (inlined) is: X={kG2.x}, Y={kG2.y}, Z={kG2.z}")
            # assert exp_x == kG2.x and exp_y == kG2.y

            print(f'Test case on FPM k={k} succeeded!\n')

def PM_test():
    golden = "sm2_golden.data"
    with open(golden, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            tokens = line.split(',')
            if len(tokens) != 3:
                print(f"Warning: Expected 3 tokens, got {len(tokens)} in line: {line}")
                continue
            k_str, var_x_str, var_y_str = tokens
            k = int(k_str.strip())
            exp_x = int(var_x_str.strip())
            exp_y = int(var_y_str.strip())

            # Test k point
            kG = k_point(k, G)
            kG = to_affine(kG.x, kG.y, kG.z)
            kG = point_exit_domain(kG)
            # print(f"k_point of G is: X={kG.x}, Y={kG.y}, Z={kG.z}")
            assert exp_x == kG.x and exp_y == kG.y

            # kG2 = k_point_inlined(k, G)
            # kG2 = to_affine(kG2.x, kG2.y, kG2.z)
            # # print(f"k_point of G (inlined) is: X={kG2.x}, Y={kG2.y}, Z={kG2.z}")
            # assert exp_x == kG2.x and exp_y == kG2.y

            print(f'Test case on PM k={k} succeeded!\n')

if __name__ == "__main__":
    # 1) fixed k point
    FPM_test()

    # 2) general k point
    PM_test()
