from dataclasses import dataclass

@dataclass
class ProjectiveCoord:
    x: int
    y: int
    z: int

# Curve parameters for SM2 ecc
a = 0xfffffffeffffffffffffffffffffffffffffffff00000000fffffffffffffffc
b = 0x28e9fa9e9d9f5e344d5a9e4bcf6509a7f39789f515ab8f92ddbcbd414d940e93
p = 0xfffffffeffffffffffffffffffffffffffffffff00000000ffffffffffffffff
Gx = 0x32c4ae2c1f1981195f9904466a39c9948fe30bbff2660be1715a4589334c74c7
Gy = 0xbc3736a2f4f6779c59bdcee36b692153d0a9877cc62a474002df32e52139f0a0
G = ProjectiveCoord(Gx, Gy, 1)


# Precomputed Look-Up Table (LUT)
precomputed = [(1, 1),
               (22963146547237050559479531362550074578802567295341616970375194840604139615431, 85132369209828568825618990617112496413088388631904505083283536607588877201568),
               (67705117572652837114471349525499634045614794997337890197046851983230464312599, 105231369456097397661386898492232980821817614166177179504943077059142844337611),
               (59085058512652262766744473634923081986036205385486894268141155278033732483669, 109091573048984472073089738700806795212297391679525613544153843136436839439305),
               (82580483505579888048248891498041876314786199024995590273835293449129297593069, 73029124966741147116568609485132862990712545849247788159814994149391419696555),
               (86848041686756723114584433308299608942636941729527524331542225209405118928267, 54921675316658680818388274484536285985620228771232775428701918326977374721529),
               (53209117963129755799540629530870949355596187090964143013579869346335355383829, 45310597926631255378856870098738145979682050257252233157641454365038168970041),
               (98567493925580762823069255045984118054797882661099345921668929848869459512328, 110336070580595011796555767631269242586466863542034033617751385756889574188365),
               (54842370261918212924648761054027838747829661419062947173569851166254391481878, 101426701570716438575530100482748050191536102803600044520927466495781527040827),
               (12269061587400982443539020516506469528566294741841090734900028672205580665746, 16529702793016095588534532837075850571549425155399206504289122598772051188685),
               (64434573173236480969138577404833282673523627011536882749755581097703792374498, 38739535235004306990435750688164009769761761056360970669716744546640331203293),
               (34527691700877729035117129267209078061476481175212747565893731231603176390121, 4057743830030335104756212222018194631112125931021813403712247080186416358025),
               (113019814001216926601167366326200148303289071208068658154223316262167530671662, 65771162203637495562352007069135952168431207244049062569251477870437195461279),
               (98806970682279766809788238647846303670311504174553780897122205688916923498009, 80258702523420980907564925663081481092346280892548128526246151985709225657322),
               (75015959162397789281735403142210955999780934098796070639937316793192444562854, 102622446084784894172630081175370402215434615166783993228004681223075583073513),
               (98290966315732394635804879472297552049863842953830735255423745352015426579195, 86045662416818120197860942815200917482442774209509364262885407182405740581028)]

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
    d = 64

    for i in range(d - 1, -1, -1):
        q = double_point(q)
        ki = int_by_slider(k, d, i)

        if ki == 0:
            continue

        pre_x, pre_y = precomputed[ki]
        pre = ProjectiveCoord(pre_x, pre_y, 1) # note, only when ki > 0, Z <-- 1

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
    d = 64  # load from HW

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

        if ki == 0:
            continue

        # Use precomputed points based on the value of ki
        new_x, new_y = precomputed[ki]  # load from HW
        pre = ProjectiveCoord(new_x, new_y, 1)

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
    golden = "sm2_fpm_golden.data"
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
            # print(f"k_point of G is: X={kG.x}, Y={kG.y}, Z={kG.z}")
            assert exp_x == kG.x and exp_y == kG.y


            kG2 = k_point_fixed_inlined(k, G.x, G.y, G.z)
            kG2 = to_affine(kG2.x, kG2.y, kG2.z)
            # print(f"k_point of G (inlined) is: X={kG2.x}, Y={kG2.y}, Z={kG2.z}")
            assert exp_x == kG2.x and exp_y == kG2.y

            print(f'Test case on k={k} succeeded!\n')

if __name__ == "__main__":
    main()
