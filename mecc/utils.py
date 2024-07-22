
msb = lambda n: n.bit_length() - 1
turn_off_msb = lambda n: n & ~(1 << (n.bit_length() - 1))


def _naf(mult):
    """Standard way of calculating non-adjacent form of number."""
    ret = []
    while mult:
        if mult % 2:
            nd = mult % 4
            if nd >= 2:
                nd -= 4
            ret.append(nd)
            mult -= nd
        else:
            ret.append(0)
        mult //= 2
    return ret


def naf(k):
    ''' Improved: Non-adjacent form of integer k'''
    naf = []
    while k > 0:
        if k & 1:  # Check if the least significant bit (LSB) is 1
            if k & 2:  # Check if the second LSB is also 1
                naf.append(-1)
                k += 1  # This effectively skips the next bit
            else:
                naf.append(1)
        else:
            naf.append(0)
        k >>= 1  # Right shift k by 1 bit (equivalent to k // 2)
    return naf


def naf_prodinger(x):
    # https://en.wikipedia.org/wiki/Non-adjacent_form
    xh = x >> 1
    x3 = x + xh
    c = xh ^ x3
    np = x3 & c
    nm = xh & c
    return np, nm


def wnaf(k):
    """Convert an integer to its Non-adjacent Form (NAF)."""
    naf_coeffs = []
    while k != 0:
        if k & 1:  # k is odd
            # Compute the remainder of k mod 4
            mod4 = k % 4
            if mod4 == 3:
                naf_coeffs.append(-1)
                k = (k + 1) // 2  # Increase k to make it even
            else:
                naf_coeffs.append(1)
                k = (k - 1) // 2  # Decrease k to make it even
        else:
            naf_coeffs.append(0)
            k //= 2  # Divide k by 2
    return naf_coeffs


# k = 1457
# naf_k = wnaf(k)
# kp, kn = naf_prodinger(1457)
# print(naf_k, kp, kn)


def int_by_slider(n: int, d: int, i: int) -> int:
    ''' From the position {i}, slide the given integer {n} by {d}-steps and take out the bit
        to form a new integer
        e.g.
            >>> n = 1457
            >>> x = int_by_slider(n, d=4, i=0)  # x=1
            >>> y = int_by_slider(n, d=4, i=1)  # x=14
    '''
    index = 0
    n >>= i
    r = 0
    while n > 0:
        b = n & 1
        r |= (b << index)

        n >>= d
        index += 1
    return r


# t = 27
# naf_t_correct = _naf(t)
# naf_t = naf(t)
#
# np_t, nm_t = naf_prodinger(t)
# assert naf_t == naf_t_correct
# print(naf_t)

# quick test...
# for i in range(111, 1000000):
#     t = i
#     naf_t_correct = _naf(t)
#     naf_t = naf(t)
#     assert naf_t == naf_t_correct
#
#     np_t, nm_t = naf_prodinger(t)
#     assert t == np_t - nm_t

