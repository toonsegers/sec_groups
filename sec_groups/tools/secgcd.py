from fractions import Fraction
from functools import reduce
from itertools import chain
import logging
from math import gcd, copysign, floor, log, log2
from operator import add, mul
import random
import sys
import time

from mpyc.runtime import mpc
from mpyc.sectypes import SecureInteger
from mpyc.finfields import GF
from mpyc.mpctools import reduce as mpyc_reduce


logger_sd = logging.getLogger("secure_division")
logger_sd.setLevel(logging.INFO)

logger_sm = logging.getLogger("secure_montgomery_exponentiation")
logger_sm.setLevel(logging.INFO)

logger_sx = logging.getLogger("secure_xgcd")
logger_sx.setLevel(logging.INFO)

logger_sb = logging.getLogger("secure_binary_xgcd")
logger_sb.setLevel(logging.DEBUG)


sign = lambda x: (1, -1)[x < 0] 


def extended_euclid_xgcd(a, b):
    """
    Returns d, u, v = xgcd(a,b)
    Where d = ua + vb = gcd(a, b)
    """
    s = 0
    old_s = 1
    t = 1
    old_t = 0
    r = b
    old_r = a

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    d, u, v = old_r, old_s, old_t
    return d, u, v


def binary_xgcd(x, y):
    """Binary extended GCD algorithm.

    Handbook of applied cryptography: [http://cacr.uwaterloo.ca/hac/about/chap14.pdf]
    by A. Menezes, P. van Oorschot, and S. Vanstone, CRC Press, 1996.
    Following Section 14.61 Algorithm Binary extended gcd algorithm.

    The numberof bits needed to represent either u or v decreases by (at least) 1,
    after at most two iterations of steps 4–7; thus,the algorithm takes at most
    2 (floor(lg(x))+floor(lg(y))+2) such iterations.
    """
    g = 1
    while (x | y) & 1 == 0:
        g = 2 * g
        x >>= 1
        y >>= 1
    u = x
    v = y

    A, B = 1, 0
    C, D = 0, 1

    while u != 0:
        while (u & 1) == 0:
            u >>= 1
            if (A | B) & 1 == 0:
                A >>= 1
                B >>= 1
            else:
                A = (A + y) >> 1
                B = (B - x) >> 1

        while (v & 1) == 0:
            v >>= 1
            if (C | D) & 1 == 0:
                C >>= 1
                D >>= 1
            else:
                C = (C + y) >> 1
                D = (D - x) >> 1
        if u >= v:
            u = u - v
            A = A - C
            B = B - D
        else:
            v = v - u
            C = C - A
            D = D - B

    a = C
    b = D
    return g * v, a, b


def shortgcd2(f, g):
    """Original shortgcd2 algorithm from [BY19].

    Invariant: f should be odd.
    See paper: https://eprint.iacr.org/2019/266
    """
    delta = 1
    assert f & 1
    m = 4 + 3 * max(f.bit_length(), g.bit_length())
    for _ in range(m):
        if delta > 0 and g & 1:
            delta, f, g = -delta, g, -f
        delta, g = 1 + delta, (g + (g & 1) * f) // 2
    return abs(f)


@mpc.coroutine
async def sec_shortgcd2(f, g):
    """Secure version of shortgcd2 from [BY19].

    Invariant: f should be odd at start.
        f remains odd throughout. Therefore, we only calculate lsb_g
        once per loop.
    See paper: https://eprint.iacr.org/2019/266
    """
    assert isinstance(f, SecureInteger)
    assert isinstance(g, SecureInteger)
    secint = type(g)
    await mpc.returnType(secint)
    delta = secint(1)
    assert await mpc.output(mpc.lsb(f))
    m = 4 + 3 * max(f.bit_length, g.bit_length)

    for i in range(m):
        lsb_g = mpc.lsb(g)
        # Optimized version of: check = (delta > 0) * lsb_g
        check = (1 - mpc.sgn(delta - 1, LT=True, l=(i + 1).bit_length() + 1)) * lsb_g
        delta, f, g = mpc.if_else(check, [-delta, g, -f], [delta, f, g])
        # truediv is less expensive and allowed here
        delta, g = 1 + delta, (g + lsb_g * f) / 2

    return mpc.abs(f)


@mpc.coroutine
async def to_bits_approx(a, l=None):
    """Secure extraction of l (or all) least significant bits of a,
    correct up to and including the least significant 1 (if any).
    """
    stype = type(a)  # secint
    if l is None:
        l = stype.bit_length
    await mpc.returnType(stype, l)
    field = stype.field

    r_bits = await mpc.random_bits(field, l)
    r_modl = 0
    for r_i in reversed(r_bits):
        r_modl <<= 1
        r_modl += r_i.value
    k = mpc.options.sec_param
    r_divl = mpc._random(field, 1<<(stype.bit_length + k - l)).value
    a = await mpc.gather(a)
    c = await mpc.output(a + ((1<<stype.bit_length) + (r_divl << l) + r_modl))
    c = c.value % (1<<l)
    return [1-r if (c >> i)&1 else r for i, r in enumerate(r_bits)]


def secure_even_gcd(a, b):
    x = to_bits_approx(a)  # low to high bits, correct up to and including the first 1 (if any)
    y = to_bits_approx(b)  # low to high bits
    z = mpc.vector_sub(mpc.vector_add(x, y), mpc.schur_prod(x, y))
    return secure_norm(z, msb=False, power=True)[0]  # TODO: clean up secure_norm()


@mpc.coroutine
async def secure_gcd(a, b):
    """Secure safegcd, generalized to both odd and even inputs.

    safegcd2 as defined in BY19 requires odd a. This protocol extends
    the functionality to even a 'by first finding the number of shared
    powers of 2 in a and b and then reducing to the odd case.' ([BY19])
    """
    secint = type(a)
    await mpc.returnType(secint)

    pow_of_2 = secure_even_gcd(a, b)
    a = a / pow_of_2
    b = b / pow_of_2

    # If a is even, swap a and b
    a, b = mpc.if_else(mpc.lsb(a), [a, b], [b, a])

    gcd_of_remainder = sec_shortgcd2(a, b)
    result = pow_of_2 * gcd_of_remainder

    return result


def truncate(f, t):
    """Truncate input f to t coefficients.

    Source: Bernstein & Yang 2019, 'Fast constant-time gcd computation and modular inversion'
    Link to paper: https://eprint.iacr.org/2019/266
    Link to original SAGE code: https://gcd.cr.yp.to/safegcd/11.sage
    """
    if t == 0:
        return 0
    twot = 1 << (t - 1)
    return ((f + twot) & (2 * twot - 1)) - twot


def divsteps2(n, t, delta, f, g):
    """Simple fast constant-time “division steps”.

    'When division steps are iterated a constant number of times, they reveal
    the gcd of the inputs, the modular reciprocal when the inputs are coprime, etc.'

    Source: Bernstein & Yang 2019, 'Fast constant-time gcd computation and modular inversion'
    Link to paper: https://eprint.iacr.org/2019/266
    Link to original SAGE code: https://gcd.cr.yp.to/safegcd/11.sage
    """
    assert t >= n and n >= 0
    f, g = truncate(f, t), truncate(g, t)
    u, v, q, r = Fraction(1, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)

    while n > 0:
        f = truncate(f, t)

        if delta > 0 and g & 1:
            delta, f, g, u, v, q, r = -delta, g, -f, q, r, -u, -v
        g0 = g & 1
        delta, g, q, r = (
            1 + delta,
            (g + g0 * f) // 2,
            (q + g0 * u) / 2,
            (r + g0 * v) / 2,
        )
        n, t = n - 1, t - 1

        g = truncate(int(g), t)
        # print("f=", f)
        # print("g=", g)

    return delta, f, g, [[u, v], [q, r]]


def jumpdivsteps2(n, t, delta, f, g):
    assert t >= n and n >= 0
    if n <= 1:
        return divsteps2(n, t, delta, f, g)
    j = n // 2
    delta, f1, g1, P1 = jumpdivsteps2(j, j, delta, f, g)

    f, g = P1[0][0] * f + P1[0][1] * g, P1[1][0] * f + P1[1][1] * g
    f, g = truncate(int(f), t - j), truncate(int(g), t - j)
    delta, f2, g2, P2 = jumpdivsteps2(n - j, n - j, delta, f, g)

    f, g = P2[0][0] * f + P2[0][1] * g, P2[1][0] * f + P2[1][1] * g
    f, g = truncate(int(f), t - n + 1), truncate(int(g), t - n)

    P3 = [
        [sum(a * b for a, b in zip(P2_row, P1_col)) for P1_col in zip(*P1)]
        for P2_row in P2
    ]

    return delta, f, g, P3


def iterations(d):
    """Bound on number of iterations for divstep.

    See Theorem 11.2 in BY19.
    Source: Bernstein & Yang 2019, 'Fast constant-time gcd computation and modular inversion'
    Link to paper: https://eprint.iacr.org/2019/266
    Link to original SAGE code: https://gcd.cr.yp.to/safegcd/11.sage
    """
    return (49*d + 80) // 17 if d < 46 else (49*d + 57) // 17


def gcd2(f, g, jumpdivsteps=False):
    """Compute gcd in constant time for odd f.

    Source: Bernstein & Yang 2019, 'Fast constant-time gcd computation and modular inversion'
    Link to paper: https://eprint.iacr.org/2019/266
    Link to original SAGE code: https://gcd.cr.yp.to/safegcd/11.sage
    """
    assert f & 1
    d = max(f.bit_length(), g.bit_length())
    m = iterations(d)

    if jumpdivsteps:
        delta, fm, gm, P = jumpdivsteps2(m, m + d, 1, f, g)
    else:
        delta, fm, gm, P = divsteps2(m, m + d, 1, f, g)

    assert gm == 0
    return abs(fm)


def recip2(f, g, jumpdivsteps=False):
    """Computes reciprocal of g modulo f.

    Requires f, g coprime; does not notice if not the case (see p. 26 of BY)
    """
    assert f & 1
    d = max(f.bit_length(), g.bit_length())
    m = iterations(d)
    precomp = (
        ((f + 1) // 2) ** (m - 1)
    ) % f  # // is allowed because f is odd (invariant)

    if jumpdivsteps:
        delta, fm, gm, P = jumpdivsteps2(m, m + 1, 1, f, g)
    else:
        delta, fm, gm, P = divsteps2(m, m + 1, 1, f, g)

    assert gm == 0
    p = P[0][1] * (2 ** (m - 1))
    V = int(copysign(1, fm)) * int(p)
    result = V * precomp
    return result % f


def extended_safegcd(f, g):
    v = gcd2(f, g)
    f_prime = f // v
    g_prime = g // v
    r = recip2(f // v, g // v)  # recip2 computes reciprocal of g mod f.
    v2 = gcd2(f_prime, g_prime)
    assert v2 == 1
    b = r % (f // v)
    a = (1 - (b * g_prime)) // (f_prime)
    assert a * f + b * g == v
    return v, a, b


def general_extended_safegcd(f, g):
    """Expand extended_safegcd for even f.
    """
    # Divide by common powers of 2.
    pow_of_2 = 1
    while (f | g) & 1 == 0:
        pow_of_2 = 2 * pow_of_2
        f >>= 1
        g >>= 1

    # If f is even, swap f and g
    swapped = 0
    if f & 1 == 0:
        swapped = 1
        f, g = g, f

    # Run extended gcd for odd f
    v, a, b = extended_safegcd(f, g)

    # Swap f and g back if needed
    if swapped:
        a, b = b, a
        f, g = g, f

    return v * pow_of_2, a, b


def divsteps22(a, b, l=None):
    """New divsteps without 2-adic computation."""
    delta, f, v, g, r = 1, a, 0, b, 1

    for i in range(iterations(l)):
        delta_gt0 = delta > 0
        g_0 = g%2
        if delta_gt0 * g_0:
            delta, f, v, g, r = -delta, g, r, -f, -v
        if g_0:
            g, r = g + f, r + v
        if r % 2:
            r = r - a
        delta, g, r = delta+1, g//2, r//2

    s = sign(f)
    f, v = s*f, s*v
    u = (f - v * b) // a
    return f, u, v


def even_gcd(x, y):
    g = 1
    while (x | y) & 1 == 0:
        g = 2 * g
        x >>= 1
        y >>= 1
    return g


def safe_xgcd(a, b, l=None):
    """Secure extended GCD based on BY.

    Bit length l for a and b, if known.
    Generalized to both odd and even inputs.
    """
    pow_of_2 = even_gcd(a, b)
    a, b = a//pow_of_2, b//pow_of_2 
    swap = 1 - a%2  # If a is even, swap a and b
    if swap:
        a, b = b, a
    
    if not l:
        l = int(log2(max(abs(a), abs(b))))+1

    g, u, v= divsteps22(a, b, l)
    if swap:
        u, v = v, u
    g = f * pow_of_2 
 
    return g, u, v



def secure_divsteps2(n, f, g):
    secint = type(g)

    delta = secint(1)
    u, v, q, r = secint(1), secint(0), secint(0), secint(1)

    for i in range(n):
        lsb_g = mpc.lsb(g)
        delta_gt0_g_odd = (
            1 - mpc.sgn(delta - 1, LT=True, l=(i + 1).bit_length() + 1)
        ) * lsb_g
        delta, f, g, u, v, q, r = mpc.if_else(
            delta_gt0_g_odd, [-delta, g, -f, q, r, -u, -v], [delta, f, g, u, v, q, r]
        )
        delta, g, q, r = (
            1 + delta,
            (g + lsb_g * f) / 2,
            (q + lsb_g * u) / 2,
            (r + lsb_g * v) / 2,
        )

    return f, g, [[u, v], [q, r]]


# def secure_jumpdivsteps2(n, t, delta, f, g):
    # assert t >= n and n >= 0
    # if n <= 1:
        # return secure_divsteps2(n, t, delta, f, g)
    # j = n // 2
    # delta, f1, g1, P1 = secure_jumpdivsteps2(j, j, delta, f, g)

    # f, g = P1[0][0] * f + P1[0][1] * g, P1[1][0] * f + P1[1][1] * g
    # # f, g = truncate(int(f), t-j), truncate(int(g), t-j)
    # delta, f2, g2, P2 = secure_jumpdivsteps2(n - j, n - j, delta, f, g)

    # f, g = P2[0][0] * f + P2[0][1] * g, P2[1][0] * f + P2[1][1] * g
    # # f, g = truncate(int(f), t-n+1), truncate(int(g), t-n)

    # P3 = [
        # [sum(a * b for a, b in zip(P2_row, P1_col)) for P1_col in zip(*P1)]
        # for P2_row in P2
    # ]

    # return delta, f, g, P3


async def secure_recip2(f, g, bit_length=None):
    """Computes reciprocal of g modulo f for secret f and g.

    Requires f, g coprime; does not notice if not the case (see p. 26 of BY)
    """
    # assert f & 1
    secint = type(g)

    if not bit_length:
        d = max(f.bit_length, g.bit_length)
    else:
        d = bit_length
    m = iterations(d)

    precomp = (1 - hensel_inverse_mod2k(f, m-1) * f) / 2**(m-1)
    fm, gm, P = secure_divsteps2(m, f, g)
    # assert gm == 0
    V = fm * P[0][1] * (2 ** (m - 1))  # NB: fm in {1, -1}
    _, recip = secure_division(V * precomp, f)  # reduction mod f;

    return recip


def secure_gcd2(f, g, bit_length=None):
    """Invariant: f is odd.

    d, m are public values.
    """
    secint = type(g)

    if not bit_length:
        d = max(f.bit_length, g.bit_length)
    else:
        d = bit_length
    m = iterations(d)

    fm, gm, P = secure_divsteps2(m, f, g)

    return mpc.abs(fm)


def secure_divsteps22(a, b, l=None):
    """New divsteps without 2-adic computation."""
    secint = type(a)
    l = l or secint.bit_length
    delta, f, v, g, r = secint(1), a, secint(0), b, secint(1)
    for i in range(iterations(l)):
        delta_gt0 = 1 - mpc.sgn(delta-1, LT=True, l=(i+1).bit_length()+1)
        g_0 = g%2
        delta, f, v, g, r = mpc.if_else(delta_gt0 * g_0,
                                       [-delta, g, r, -f, -v],
                                       [delta, f, v, g, r])
        g, r = mpc.if_else(g_0, [g + f, r + v], [g, r])
        r = mpc.if_else(r%2, r - a, r)
        delta, g, r = delta+1, g/2, r/2
    s = 1 - 2*(f < 0)  # sign of f
    f, v = mpc.scalar_mul(s, [f, v])
    u = (f - v * b) / a
    return f, u, v


def secure_xgcd(a, b, l=None):
    """Secure extended GCD based on BY.

    Bit length l for a and b, if known.
    Generalized to both odd and even inputs.
    safegcd2 as defined in BY19 requires odd a. This protocol extends
    the functionality to even a 'by first finding the common factors
    of 2 in a and b and then reducing to the odd case.' ([BY19])
    """
    # Divide by common powers of 2.
    pow_of_2 = secure_even_gcd(a, b)
    a, b = mpc.scalar_mul(1/pow_of_2, [a, b])
    swap = 1 - a%2  # If a is even, swap a and b
    a, b = mpc.if_else(swap, [b, a], [a, b])
    f, u, v = secure_divsteps22(a, b, l)
    u, v = mpc.if_else(swap, [v, u], [u, v])
    return f * pow_of_2, u, v


def flip(x):
    return 1 - x


def secure_norm(x, msb=True, power=True):
    """Determine most or least significant bit.

    Logarithmic round complexity and low computational overhead.

    Args:
        x (list): List of secure integers in {0, 1} (i.e. bits).
        msb (boolean): Flag to set most or least significant bit mode.
        power (boolean): Flag to return 2^pos or pos,
        where pos is the position of msb/lsb.
    """
    if msb:
        x = list(reversed(x))
    if power:
        nz, i = mpc.find(x, 1, e=None, cs_f=lambda b, i: (b+1) << i)
    else:
        nz, i = mpc.find(x, 1, e=None, cs_f=lambda b, i: b + i)
    if msb:
        if power:
            i = (1 << len(x)) / i  # proper division
        else:
            i = len(x) - i
    return i, 1 - nz


@mpc.coroutine
async def secure_binary_xgcd(x, y, progress_bar=False, allow_negative=True, bit_length = None):
    """Secure binary extended GCD algorithm.
    """
    secint = type(x)
    await mpc.returnType(secint, 3)

    # Securely calculate bit_length if not provided.
    if not bit_length:
        x_bl = secure_norm(mpc.to_bits(x), msb=True, power=False)[0]
        y_bl = secure_norm(mpc.to_bits(y), msb=True, power=False)[0]
        x_bl = int(await mpc.output(x_bl)) + 1
        y_bl = int(await mpc.output(y_bl)) + 1
        print("x_bl=", x_bl)
        print("y_bl=", y_bl)
        bit_length = max(x_bl, y_bl)
        logger_sx.debug(
            f"In secure_binary_xgcd, bit_length securely calculated (at performance penalty) = {bit_length}"
        )

    l = bit_length
    """The numberof bits needed to represent either u or v decreases by (at least) 1,
    after at most two iterations of steps 4–7; thus,the algorithm takes at most
    2 (lg(x)+lg(y)+2) such iterations. (See Handbook, p. 610)
    """
    # m = 2 * (x.bit_length + y.bit_length)
    m = 2*2*bit_length


    if progress_bar:
        toolbar_width = m
        sys.stdout.write("Constructing circuit: [%s]" % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write(
            "\b" * (toolbar_width + 1)
        )  # return to start of line, after '['

    if allow_negative:
        x_nonneg = x >= 0
        y_nonneg = y >= 0
        x = mpc.if_else(x_nonneg, x, -x)
        y = mpc.if_else(y_nonneg, y, -y)

    g = secure_even_gcd(x, y)
    x = x / g
    y = y / g
    u = x
    v = y

    A, B = secint(1), secint(0)
    C, D = secint(0), secint(1)

    # Page 610: "algorithm takes at most 2 (floor(lg(x))+floor(lg(y))+2) such iterations"
    for i in range(m):

        not_done = flip(mpc.is_zero(u))
        for i in range(l):
            not_done *= flip(mpc.lsb(u))  # TODO: avoid separate calls
            AB_both_even = flip(mpc.lsb(A)) * flip(mpc.lsb(B))
            u, A, B = mpc.if_else(
                not_done,
                [u / 2]
                + mpc.if_else(AB_both_even, [A / 2, B / 2], [(A + y) / 2, (B - x) / 2]),
                [u, A, B],
            )

        # TODO: avoid unnecessary execution of both loops by obliously swapping inputs
        not_done = secint(1)
        for i in range(l):
            not_done *= flip(mpc.lsb(v))
            CD_both_even = flip(mpc.lsb(C)) * flip(mpc.lsb(D))
            v, C, D = mpc.if_else(
                not_done,
                [v / 2]
                + mpc.if_else(CD_both_even, [C / 2, D / 2], [(C + y) / 2, (D - x) / 2]),
                [v, C, D],
            )

        # TODO: consider using bit representation and sgn?
        u, A, B, v, C, D = mpc.if_else(
            u >= v, [u - v, A - C, B - D, v, C, D], [u, A, B, v - u, C - A, D - B]
        )

        if progress_bar:
            sys.stdout.write("-")
            sys.stdout.flush()

    a = C
    b = D

    if progress_bar:
        sys.stdout.write("]\n")  # this ends the progress bar

    if allow_negative:
        a = mpc.if_else(x_nonneg, a, -a)
        b = mpc.if_else(y_nonneg, b, -b)

    return g * v, a, b


def pow_list(a, x, n):
    """Return [a,ax, ax^2, ..., ax^(n-1)].

    Runs in O(log n) rounds using minimal number of n-1 secure multiplications.
    NB: equivalent to list(mpyc.mpctools.accumulate([x] * (n-1), f=operator.mul, iv=a)),
    which also runs in O(log n) rounds but using O(n log n) secure multiplications.
    """
    if n == 1:
        powers = [a]
    elif n == 2:
        powers = [a, a * x]
    else:
        even_powers = pow_list(a, x * x, (n + 1) // 2)
        if n % 2:
            d = even_powers.pop()
        odd_powers = mpc.scalar_mul(x, even_powers)
        powers = [t for _ in zip(even_powers, odd_powers) for t in _]
        if n % 2:
            powers.append(d)
    return powers


@mpc.coroutine
async def prefix_mul(x):
    """WORK IN PROGRESS: Prefix multiplication for vector x of sectypes.

    Follows Protocol 4.17 from De Hoogh's PhD thesis,
    "Design of large scale applications of secure multiparty
    computation", TU Eindhoven, 2012.
    See: https://pure.tue.nl/ws/files/3430368/735328.pdf
    """

    sectype = type(x[0])
    k = len(x)
    await mpc.returnType(sectype, k)
    u = [None] * k
    r = [None] * k
    s = [None] * k
    # TODO: parallellize # generate all random at once and check if the vector
    for i in range(k):
        u[i] = 0
        while u[i] == 0:
            r[i] = mpc._random(sectype)
            s[i] = mpc._random(sectype)

            u[i] = await mpc.output(
                r[i] * s[i]
            )  # TODO: use Schur product to parallellize; avoid Schur; no resharing
    v = [None] * k
    for i in range(1, k):
        v[i] = r[i] * s[i - 1]  # TODO: use Schur product to parallellize
    w = [None] * k
    w[0] = r[0]
    for i in range(1, k):
        w[i] = v[i] * (
            1 / u[i - 1]
        )  # TODO: save 1/u for reuse below.  # TODO: use Schur product to parallellize
    m = [None] * k
    # TODO: parallellize
    for i in range(k):
        m[i] = await mpc.output(x[i] * w[i])  # TODO: use Schur product to parallellize
    y = [None] * k
    y[0] = x[0]
    for i in range(1, k):
        y[i] = s[i] * (1 / u[i]) * reduce(mul, m[0 : i + 1])  # note: typo in DH12.

    return y


@mpc.coroutine
async def secure_poly_A(omega, x):
    sectype = type(x)
    await mpc.returnType(sectype)
    p = prefix_mul([x] * omega)
    return sum(p) + 1


@mpc.coroutine
async def prefix_or(b):
    """Prefix or for vector b of secure field/int elements in {0,1}.

    Reversed to align it with definition of prefix-OR of DNT12.

    Follows Protocol 4.18 from De Hoogh's PhD thesis,
    "Design of large scale applications of secure multiparty
    computation", TU Eindhoven, 2012.
    See: https://pure.tue.nl/ws/files/3430368/735328.pdf
    """
    sectype = type(b[0])
    k = len(b)
    await mpc.returnType(sectype, k)

    b = b[::-1]
    z = prefix_mul([b_i + 1 for b_i in b])
    x = [None] * k
    x[0] = b[0]  # note: typo in DH12.
    for i in range(1, k):
        x[i] = 1 - mpc.lsb(z[i])
    x = x[::-1]
    return x


@mpc.coroutine
async def secure_2_pow_bit_length(a):
    sectype = type(a)
    await mpc.returnType(sectype)
    a_bits = mpc.to_bits(a)
    y = prefix_or(a_bits)
    return 1 + sum([y_i * 2 ** i for i, y_i in enumerate(y)])


@mpc.coroutine
async def secure_bit_length(a):
    sectype = type(a)
    await mpc.returnType(sectype)
    a_bits = mpc.to_bits(a)
    y = prefix_or(a_bits)
    return sum(y)


# @mpc.coroutine
# async def secure_division(n, d, l=None, l_d=None):
    # """Divides secint n by secint d.

    # Based on Dahl, Ning, Toft 2012: "On Secure Two-Party Integer Division"
    # See latest version: https://eprint.iacr.org/2012/164/20151016:230655

    # Args:
        # n (mpc.SecInt): numerator of the secure division.
        # d (mpc.SecInt): denominator of the secure division.
        # l (int): max bitlength of n and d. (optional, recommended for security)
        # l_d (int): bitlength of denominator. (optional, speed up, less secure)

    # Use sub-protocols with logarithmic round complexity and low computational
    # overhead (instead of constant round complexity, as in paper).

    # """
    # sectype = type(n)
    # original_sectype = sectype
    # await mpc.returnType(sectype, 2)

    # # If l is not provided, compute bitlength securely based on input n.
    # if not l:
        # # TODO: re-use bits between max and to_bits
        # logger_sd.debug(
            # "In secure_division(): Parameter l (max bit_length) not provided, calulcating with secure_norm()."
        # )
        # secure_max_n_d = mpc.max(n, d)
        # l = secure_norm(mpc.to_bits(secure_max_n_d), msb=True, power=False)[0]
        # l = int(await mpc.output(l)) + 1
        # logger_sd.debug(
            # f"In secure_division(): Parameter l (max bit_length) calculated: l = {l}"
        # )

    # # omega >= l_n - l_d sufficient to ensure approximation error of 2^k / d below 2^(k-l_d-ω).
    # omega = l

    # # Convert to larger sectype, required to fit terms of Taylor series and avoid wrap-around.
    # k = l ** 2 + l
    # safe_bit_length = (
        # l ** 2 + 2 * l
    # )  # Corresponds to l**2 + l in paper, but intermediate step n * a_tilde requires l**2 + 2l.
    # conversion_required = False
    # original_sectype_bit_length = original_sectype.bit_length
    # if original_sectype_bit_length < safe_bit_length:
        # logger_sd.debug(
            # f"In secure_division(): Sectype bit_length {original_sectype_bit_length} insufficient (l={l}), convert to sectype of bit_length {safe_bit_length} (excl. headroom). Required for Taylor series."
        # )
        # conversion_required = True
        # sectype = mpc.SecInt(
            # safe_bit_length
        # )  # kappa_s = 30 is implicit (via -K flag in MPYC)
        # d_original = d
        # n_original = n
        # n = mpc.convert(n, sectype)
        # d = mpc.convert(d, sectype)

    # # Calculate 2^l_d and its reciprocal.
    # if isinstance(l_d, int):
        # # If l_d is public, compute reciprocal directly.
        # two_to_ld = 2 ** l_d
        # two_to_ld_inv = type(n).field(two_to_ld).reciprocal()
    # else:
        # # If l_d is secret, follow paper. (log complexity)
        # two_to_ld = secure_norm(mpc.to_bits(d), msb=True, power=True)[0]
        # two_to_ld_inv = mpc.reciprocal(two_to_ld)

    # # Calculate a_tilde, the estimate of floor((2^k)/d). (log complexity)
    # p = (two_to_ld - d) * two_to_ld_inv
    # a_tilde_list = pow_list((2 ** k) * two_to_ld_inv, p, omega + 1)
    # a_tilde = mpyc_reduce(add, a_tilde_list, 0)

    # # Calculate q_tilde, the estimate of floor(n/d)
    # q_hat = n * a_tilde
    # q_tilde = mpc.trunc(q_hat, k)  # TODO: mpc.trunc is probabilistic. Consider performing a full bit-decomposition trunc.

    # # Convert back to original (smaller) sectype.
    # if conversion_required:
        # q_tilde = mpc.convert(q_tilde, original_sectype)
        # d = d_original
        # n = n_original

    # # Calculate remainder and correct q_tilde for approximation error.
    # r = n - d * q_tilde
    # epsilon_plus = mpc.ge(r + d, 2 * d)
    # epsilon_minus = mpc.ge(d - 1, r + d)
    # q = q_tilde + epsilon_plus - epsilon_minus

    # # Re-calculate remainder with correct q. (Figure 1 in DNT12 skips this step.)
    # r = n - d * q

    # return q, r


@mpc.coroutine
async def secure_division(a, b, l=None, l_d=None):
    """Integer division divmod(a, b) via NR. Ignores parameters l, l_d."""
    secint = type(a)
    await mpc.returnType(secint, 2)
#    assert await mpc.output(b>0)
    secfxp = mpc.SecFxp(2*secint.bit_length+2)
    a1, b1 = mpc.convert([a, b], secfxp)
    q = a1 / b1
    q = mpc.convert(q, secint)
    r = a - b * q
    q, r = mpc.if_else(r < 0, [q - 1, r + b], [q, r])  # correction using one <
#    q, r = mpc.if_else(r >= b, [q + 1, r - b], [q, r])  # correction using one <
#    assert await mpc.output(a == b * q + r), await mpc.output([q, r, a, b])
#    assert await mpc.output(0 <= r), await mpc.output([q, r, a, b])
#    assert await mpc.output(r < b), await mpc.output([q, r, a, b])
    return q, r


def hensel_inverse_mod2k(a, k):
    """Apply Hensel lifting to calculate modular inverse of a mod 2^k, k>=1.

    Solves f(x) = ax - 1 = 0 modulo 2^k iteratively, using
    a_k+1 = a_k - f(a_k) mod 2^k+1.
    """
    # TODO: fix doc string / function name to cover use of Newton-Raphson
    # Newton-Raphson with quadratic convergence:
    y = type(a)(1)
    # n = 2
    # for i in range(k.bit_length()-1):
        # # NR iteration y = (y * (2 - a*y)) % 2**(2**(i+1)) to double number of bits of y
        # y += (((y * (1 - a*y)) / n) % n) * n
        # n **= 2
    # if  n < 2**k:
        # y += (((y * (1 - a*y)) / n) % (2**k // n)) * n
    # return y

    # use Hensel lifting with linear convergence
    for i in range(1, k):
        # compute y one bit at a time, using one secure multiplication and one secure lsb
        y += (((a * y - 1)/(1<<i)) % 2) * (1<<i)  # NB: a*y - 1 equal to 0 or 2**i
    return y

    # b = a + 1
    # y = 1
    # for i in range(2, k + 1):
        # y = (b * y - 1) % 2 ** i  # NB: y + f(y) = y + ay - 1 = by - 1
    # return y


def mont(a, b, r, n, n_prime):
    """Montgomery multiplication.

    Compute Montgomery reduction of the product of two integers.
    Taken from: https://www.microsoft.com/en-us/research/wp-content/uploads/1996/01/j37acmon.pdf

    Args:
        a, b (int): Integers to multiply and apply Montgomery reduction to.
        r (int): Auxiliary modulus.
        n (int): original modulus.
        n_prime (int): -1 * modular inverse of n mod r.
    """
    t = a * b
    t_np_mod_r = (t * n_prime) % r
    u = (t + (t_np_mod_r) * n) // r
    if u >= n:
        return u - n
    else:
        return u


def montgomery_exponentiation(a, e, r, n, n_prime):
    """Modular exponentiation, a^e mod n, using Montgomery exponentiation.

    Algorithm follows Algorithm 14.94 from Handbook of Applied Cryptography, p. 620.
    """
    a_tilde = mont(a, r ** 2 % n, r, n, n_prime)

    result = r % n

    if e == 0:
        return result

    for i in range(e.bit_length() - 1, -1, -1):
        result = mont(result, result, r, n, n_prime)
        if (e >> i) & 1:
            result = mont(result, a_tilde, r, n, n_prime)

    result = mont(result, 1, r, n, n_prime)
    return result


async def secure_mont(a, b, r, n, n_prime):
    """Montgomery multiplication.

    See mont() docstring for details.
    """
    assert type(a).bit_length >= 3 * (r.bit_length())  # Ensure a*b*n_prime fits
    t = a * b
    t_np_mod_r = t*n_prime % r
    u = (t + t_np_mod_r * n) / r
    return mpc.if_else(u >= n, u - n, u)


async def secure_montgomery_exponentiation(a, e, r, n, n_prime):
    """Modular exponentiation, a^e mod n, using Montgomery exponentiation.

    Algorithm follows Algorithm 14.94 from Handbook of Applied Cryptography, p. 620.
    """
    secint = type(a)
    _, r2modn = secure_division(secint(r ** 2), n)
    _, result = secure_division(secint(r), n)

    # Convert to sectype with larger bit-length if needed.
    conversion_required = False
    safe_bit_length = 3 * (r.bit_length())  # for a*b*n_prime in secure_mont()
    if secint.bit_length < safe_bit_length:
        logger_sd.debug(
            f"In secure_montgomery_exponentiation(), original sectype bit_length {secint.bit_length} insufficient, convert to sectype of bit_length {safe_bit_length} (excl. headroom)."
        )
        conversion_required = True
        original_secint = type(a)
        secint = mpc.SecInt(safe_bit_length)
        r2modn = mpc.convert(r2modn, secint)
        result = mpc.convert(result, secint)
        a = mpc.convert(a, secint)
        n = mpc.convert(n, secint)
        n_prime = mpc.convert(n_prime, secint)

    a_tilde = await secure_mont(a, r2modn, r, n, n_prime)

    if await mpc.output(a == (n+1)/2):
        """The following shows that we can remove the computation for r2modn entirely.
        This saves the probably most expensive secure integet division to compute r**2 / n.
        Only the division to compute r / n to get result remains.

        This can be done when a = 1/2 mod n, which is the case of interest when we compute
        the secure reciprocal modulo n, see secure_recip2().

        We compute a_tilde = (a * r) mod n = (r / 2) mod n from result using only one lsb (%2).
        """
        new_a_tilde = (result + (result%2) * n)/2
        assert await mpc.output(a_tilde == new_a_tilde)  # confirms new approach

    if e == 0:
        return result

    for i in range(e.bit_length() - 1, -1, -1):
        result = await secure_mont(result, result, r, n, n_prime)
        if (e >> i) & 1:
            result = await secure_mont(result, a_tilde, r, n, n_prime)

    result = await secure_mont(result, 1, r, n, n_prime)

    # Convert back to original (smaller) sectype.
    if conversion_required:
        result = mpc.convert(result, original_secint)

    return result


async def secure_modular_exp(a, e, modulus, bit_length):
    """Wraps secure Montgomery exponentiation: a^e mod modulus.

    Calculated r and n_prime before calling Montgomery exponentiation.
    """
    r = 2 ** bit_length
    n_prime = -hensel_inverse_mod2k(modulus, bit_length)
    return await secure_montgomery_exponentiation(a, e, r, modulus, n_prime)


def secgcdKnuthAlgB(u, v):
    m = type(u).bit_length
    g = secure_even_gcd(u, v)
    u /= g
    v /= g
    t = mpc.if_else(mpc.lsb(u), -v, u)
    for _ in range(m):
        x = mpc.to_bits(t)
        s = x[-1]  # sign of t
        gg, nz = secure_norm(x, msb=False, power=True)
        t /= gg
        t, u, v = mpc.if_else(
            s,
            [t + u, u, -t],  # t < 0
            mpc.if_else(nz, [t - v, t, v], [t, u, v]),  # t > 0
        )  # t = 0
    return g * v


# TODO's
# 2: Remove unnecessary "async" and "await" statements (secure_recip2, _divsteps2, _gcd2, etc.)
