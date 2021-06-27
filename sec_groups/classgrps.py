"""Form class groups.

#TODO: add references to literature.
"""

import logging
from math import floor, sqrt, log, log2

from sec_groups.tools.secgcd import (
    extended_euclid_xgcd,
    secure_gcd,
    secure_xgcd,
    secure_binary_xgcd,
    secure_division,
)
from sec_groups.tools.bitlen import bit_length_integrated
from sec_groups.tools.repeat import secure_pow
from mpyc.runtime import mpc
import mpyc.gmpy as gmpy2
from sec_groups.tools.find_primes import find_primes_for_schnorr, _find_ike_prime

logger_cg = logging.getLogger("classgroups")
logger_cg.setLevel(logging.INFO)


def xgcd_(a, b):
    """Wraps extended euclid from secgcd module."""
    return extended_euclid_xgcd(a, b)


def discriminant(f):
    a, b, c = f[0], f[1], f[2]
    return b ** 2 - 4 * a * c


def lincong(a, b, m):
    """Solve ax = b mod m

    return mu, nu such that x = mu + nu n for all n in Z.
    Based on Lipa Long, "Binary Quadratic Forms", 2019.
    See: https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf
    """
    g, d, e = xgcd_(a, m)
    logger_cg.debug(f"In lincong, done xgcd: {g}, {d}, {e} = xgcd({a}, {m})")
    q, r = divmod(b, g)
    logger_cg.debug(f"In lincong, done {q}, {r} = division({b}, {g}).")

    # L19 Thm. 7.1: Congruence has a solution iff gcd(a,m) | b.
    if r != 0:    
        raise ValueError("The linear congruence has no solution")
    else:
        mu = (q * d) % m
        logger_cg.debug(f"In lincong, done _, {mu} = division({q}*{d}, {m}).")
        nu = m // g
        return mu, nu


def secure_lincong(a, b, m):
    """Solve ax = b mod m

    return mu, nu such that x = mu + nu n for all n in Z.
    """
    g, d, e = secure_xgcd(a, m)
    logger_cg.debug(f"In lincong, done secure_xgcd().")
    # q = floor(b/g)
    # q = b // g
    # r = b % g
    #    q, r = secure_division(b, g)
    q = b / g
    r = 0
    logger_cg.debug(f"In lincong, done secure_division(b, g).")

    if isinstance(r, int) and r != 0:
        raise ValueError("The congruence has no solution")
    else:
        # mu = (q * d) % m
        _, mu = secure_division(q * d, m)
        logger_cg.debug(f"In lincong, done secure_division(q*d, m).")
        # nu = m // g
        # nu, _ = secure_division(m, g)
        nu = m / g
        return mu, nu


def check_well_formed(f):
    a, b, c = f[0], f[1], f[2]
    disc = b ** 2 - 4 * a * c
    if a > 0 and disc < 0:
        pass
    else:
        raise ValueError(
            f"Form ({a}, {b}, {c}) does not have a > 0 and discriminant < 0: a={a}, disc={disc} "
        )


def check_reduced(f):
    a, b, c = f[0], f[1], f[2]

    if -a < b and b <= a:  # check normalized
        pass
    else:
        return False
    if a <= c:
        pass
    else:
        return False
    if a == c:
        if b >= 0:
            pass
        else:
            return False

    return True


def normalize(f):
    a, b, c = f[0], f[1], f[2]
    group = type(f)
    check_well_formed(f)
    r = (a - b) // (2 * a)
    eta = (a, b + 2 * r * a, a * r ** 2 + b * r + c)
    return group(eta)


def reduce_form(f):
    group = type(f)
    check_well_formed(f)
    f = normalize(f)
    while not check_reduced(f):
        a, b, c = f[0], f[1], f[2]
        s = (c + b) // (2 * c)
        f = group((c, -1 * b + 2 * s * c, c * s ** 2 - b * s + a))
    return f


@mpc.coroutine
async def secure_binary_reduce(f, size_b = None, leak_size_b = True):
    """Binary reduction algorithm by Agarwal and Frandsen.

    Based on Algorithm 3 from AF06: 'A New GCD Algorithm for Quadratic Number
    Rings with Unique Factorization' by Agarwal and Frandsen, 2006 (Aarhus)
    https://users-cs.au.dk/gudmund/Documents/38870030.pdf

    Requires:
        f is positive definite (iff discriminant < 0 and a > 0).

    NB: Option to open (leak) size(b) is default; to reduce number of 
    iterations of main loop. Alternative is to pass size_b bound.
    """

    def size(a):
        # Requires non-negative values
        return bit_length_integrated(mpc, a)

    def right_action_S_on_f(f):
        return [f[2], -f[1], f[0]]

    def right_action_Tm_on_f(m, f):
        fa, fb, fc = f[0], f[1], f[2]
        return [fa, fb + 2 * m * fa, (m ** 2) * fa + m * fb + fc]

    sec_grp = type(f)
    await mpc.returnType(sec_grp)

    secint = sec_grp.sectype_of_value
    a, b, c = f[0], f[1], f[2]

    if size_b:
        n = size_b
    elif not size_b and leak_size_b:
        n = await mpc.output(size(b))  # TODO: find good bound for for-loop
    else:
        raise NotImplementedError

    for i in range(n):
        sgn_b = 1 - 2 * mpc.sgn(
            b, l=n + 3, LT=True
        )  # TODO: check l; if n + 0, sgn_b produces inccorect values <-1
        abs_b_gt_abs_2a = sgn_b * b > 2 * a
        abs_a_gt_abs_c = a > c  # a always postive, because f positive definite
        ab_gt_0 = (sgn_b * sgn_b + sgn_b) // 2  # a always postive, because f positive definite

        size_abs_b = size(sgn_b * b)
        size_a = size(a)
        # TODO: find bound for (bit-length of) j.
        j = size_abs_b - size_a - 1
        # take |j| to avoid negative secint exponents. 2**j is used when |B|>2|A| and original j is positive
        sgn_j = 1 - 2 * mpc.sgn(j, l=n, LT=True)
        abs_j = sgn_j * j
        abs_j_bits = mpc.to_bits(abs_j, n)
        m = secure_pow(2, abs_j_bits, secint)

        m = mpc.if_else(ab_gt_0, -m, m)

        a, b, c = mpc.if_else(
            abs_b_gt_abs_2a,
            right_action_Tm_on_f(m, (a, b, c)),
            mpc.if_else(abs_a_gt_abs_c, right_action_S_on_f((a, b, c)), [a, b, c]),
        )

        print(f"Secure binary reduction: {round(100*i/n)}%", end="\r")

    assert f.group.discriminant < 0
    m = mpc.if_else(b > 0, secint(-1), secint(1))
    abs_b_gt_a = mpc.abs(b) > a

    a, b, c = mpc.if_else(abs_b_gt_a, right_action_Tm_on_f(m, (a, b, c)), [a, b, c])
    a_gt_c = a > c
    a, b, c = mpc.if_else(
        abs_b_gt_a * a_gt_c, right_action_S_on_f((a, b, c)), [a, b, c]
    )

    a, b, c = mpc.if_else((b < 0) * (a == c), right_action_S_on_f((a, b, c)), [a, b, c])
    a, b, c = mpc.if_else(
        (b < 0) * (a == -b), right_action_Tm_on_f(1, (a, b, c)), [a, b, c]
    )

    return sec_grp((a, b, c))


def parteucl(a, b, L):
    """Extended partial Euclides following Cohen Section 5.4.
    """
    # Step 1 Initialize
    v = 0
    d = a
    v2 = 1
    v3 = b
    z = 0

    while abs(v3) > L:
        # Step 3 Euclidean step
        q, t3 = d//v3, d%v3
        t2 = v - q*v2
        v = v2
        d = v3
        v2 = t2
        v3 = t3
        z = z+1

    # Step 2 Finished?
    if z % 2:
        v2 = -v2
        v3 = -v3

    return d, v, v2, v3, z


def nudupl(f):
    """Square(f) following Cohen, Alg. 5.4.8.
    """
    L = int(((abs(f.discriminant))/4)**(1/4))
    a, b, c = f[0], f[1], f[2]

    # Step 1 Euclidean step
    d1, u, v = extended_euclid_xgcd(b, a)
    A = a//d1
    B = b//d1
    C = (-c*u) % A
    C1 = A-C
    if C1 < C:
        C = -C1

    # Step 2 Partial reduction
    d, v, v2, v3, z = parteucl(A, C, L)

    # Step 3 Special case
    if z==0:
        g = (B*v3+c)//d
        a2 = d**2
        c2 = v3**2
        b2 = b + (d+v3)**2 - a2 - c2
        c2 = c2 + g*d1
    else:
        # Step 4 Final computations
        e = (c*v + B*d)//A
        g = (e*v2 - B)//v
        b2 = e*v2 + v*g
        if d1>1:
            b2 = d1*b2
            v = d1*v
            v2 = d1*v2
        a2 = d**2
        c2 = v3**2
        b2 = b2 + (d+v3)**2 - a2 - c2
        a2 = a2 + e*v
        c2 = c2 + g*v2

    f2 = type(f)((a2, b2, c2))
    return f2


def square(f):
    """Square form"""
    group = type(f)
    a, b, c = f[0], f[1], f[2]
    mu, _ = lincong(b, c, a)
    A = a ** 2
    B = b - 2 * a * mu
    C = mu ** 2 - (b * mu - c) // a
    return group((A, B, C))


def secure_square(f):
    sectype = type(f)
    """Square form"""
    a, b, c = f[0], f[1], f[2]
    mu, _ = secure_lincong(b, c, a)
    A = a ** 2
    B = b - 2 * a * mu
    # C = mu ** 2 - (b * mu - c) // a
    C = mu ** 2 - (b * mu - c) / a
    return sectype((A, B, C))


def repeat_square(f, n):
    new_f = f
    for i in range(n):
        new_f = reduce(square(new_f))
    return new_f


def nucomp(phi1, phi2):
    """Nucomp algorithm for composition of binary quadratic forms.

    Per Jacobson and Van der Poorten 'Computational Aspects of NUCOMP', 2002
    See: https://link.springer.com/chapter/10.1007%2F3-540-45455-1_10

    Alternatively see Cohen  'A course in computational number theory'
    All divisions are exact (see Cohen, p. 244)
    """
    delta = phi1.discriminant
    # JV02 uses the following L, which is different from Cohen's L
    L = int(abs(delta) ** (1 / 4))
    # L = int(((abs(delta))/4)**(1/4))  # L used in Cohen's book

    # Step 1
    if phi1[2] < phi2[2]:
        phi1, phi2 = phi2, phi1

    u1, v1, w1 = phi1[0], phi1[1], phi1[2]
    u2, v2, w2 = phi2[0], phi2[1], phi2[2]

    s = (v1 + v2) // 2
    m = v2 - s

    # Step 2
    F, b, c = extended_euclid_xgcd(u2, u1)
    if s % F == 0:  # F | s
        G = F
        A_x = G
        B_x = m * b
        B_y = u1 // G
        C_y = u2 // G
        D_y = s // G
        # go to Step 5
    else:  # F does not divide s
        # Step 3
        G, x, y = extended_euclid_xgcd(F, s)
        H = F // G
        B_y = u1 // G
        C_y = u2 // G
        D_y = s // G

        # Step 4
        l = y * (b * (w1 % H) + c * (w2 % H)) % H
        B_x = b * (m // H) + l * (B_y // H)

    # Step 5
    b_x = B_x % B_y
    b_y = B_y

    # Step 5a
    x, y, z = 1, 0, 0
    while abs(b_y) > L and b_x != 0:
        # Step 5c
        q, t = divmod(b_y, b_x)
        b_y = b_x
        b_x = t
        t = y - q * x
        y = x
        x = t
        z += 1

    # Step 5b if not abs(b_y) > L and b_x != 0
    if z % 2 == 1:
        b_y = -b_y
        y = -y
    a_x = G * x
    a_y = G * y

    # Step 6
    if z == 0:
        Q1 = C_y * b_x
        c_x = (Q1 - m) // B_y
        d_x = (b_x * D_y - w2) // B_y
        u3 = b_y * C_y
        w3 = b_x * c_x - G * d_x 
        v3 = v2 - 2 * Q1
    else:
        # Step 7
        c_x = (C_y * b_x - m * x) // B_y
        Q1 = b_y * c_x
        Q2 = Q1 + m
        d_x = (D_y * b_x - w2 * x) // B_y
        Q3 = y * d_x
        Q4 = Q3 + D_y
        d_y = Q4 // x
        if b_x != 0:
            c_y = Q2 // b_x
        else:
            c_y = (c_x * d_y - w1) // d_x
        u3 = b_y * c_y - a_y * d_y
        w3 = b_x * c_x - a_x * d_x
        v3 = G * (Q3 + Q4) - Q1 - Q2

    return type(phi1)((u3, v3, w3))


def compose(f1, f2):
    """Composition of binary quadratic forms.

    Based on Lipa Long, "Binary Quadratic Forms", 2019.
    See: https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf
    """

    group = type(f1)
    a1, b1, c1 = f1[0], f1[1], f1[2]
    a2, b2, c2 = f2[0], f2[1], f2[2]
    # step 1
    g = (b1 + b2) // 2
    h = -(b1 - b2) // 2
    w, _, _ = xgcd_(a1, a2)
    w, _, _ = xgcd_(w, g)
    logger_cg.debug("Done with 2 gcds in compose.")
    # step 2
    j = w
    s = a1 // w
    t = a2 // w
    u = g // w
    # step 3
    logger_cg.debug("Start lincong 1.")
    mu, nu = lincong(t * u, h * u + s * c1, s * t)
    # step 4
    logger_cg.debug("Start lincong 2.")
    lmb, _ = lincong(t * nu, h - t * mu, s)
    # step 5
    k = mu + nu * lmb
    l = (k * t - h) // s
    m = (t * u * k - h * u - c1 * s) // (s * t)
    # step 6
    A = s * t
    B = j * u - (k * t + l * s)
    C = k * l - j * m
    # step 7
    return group((A, B, C))


def secure_compose(f1, f2):
    sectype = type(f1)
    a1, b1, c1 = f1[0], f1[1], f1[2]
    a2, b2, c2 = f2[0], f2[1], f2[2]
    # step 1
    g = (b1 + b2) / 2
    h = -(b1 - b2) / 2
    w = secure_gcd(a1, a2)
    w = secure_gcd(w, g)
    logger_cg.debug("Done with 2 gcds in compose.")
    # step 2
    j = w
    s = a1 / w
    t = a2 / w
    u = g / w
    # step 3
    logger_cg.debug("Start secure_lincong 1.")
    mu, nu = secure_lincong(t * u, h * u + s * c1, s * t)
    # step 4
    logger_cg.debug("Start secure_lincong 2.")
    lmb, _ = secure_lincong(t * nu, h - t * mu, s)
    # step 5
    k = mu + nu * lmb
    l = (k * t - h) / s
    m = (t * u * k - h * u - c1 * s) / (s * t)
    # step 6
    A = s * t
    B = j * u - (k * t + l * s)
    C = k * l - j * m
    # step 7
    return sectype((A, B, C))


def shanks_compose(f1, f2):
    """Composition of positive definite forms.

    Originally by Shanks 'Class Number, a theory of factorization, 
    and genera', Proc. Symp. in Pure Maths, 1969.
    Taken from Coh93, Algorithm 5.4.7.
    """
    # Step 1
    if f1[0] > f2[0]:
        f1, f2 = f2, f1
    a1, b1, c1 = f1[0], f1[1], f1[2]
    a2, b2, c2 = f2[0], f2[1], f2[2]
    s = (b1 + b2)//2
    n = b2 - s

    # Step 2: First Euclidean step
    if a2 % a1 == 0:
        y1 = 0
        d = a1
    else:
        d, u, v = extended_euclid_xgcd(a2, a1)
        y1 = u

    # Step 3: Second Euclidean step
    if s % d == 0:
        y2 = -1
        x2 = 0
        d1 = d 
    else:
        d1, u, v = extended_euclid_xgcd(s, d)
        x2 = u 
        y2 = -v

    # Step 4: Compose
    v1 = a1//d1
    v2 = a2//d1
    r = (y1*y2*n - x2*c2) % v1
    b3 = b2 + 2*v2*r 
    a3 = v1*v2 
    c3 = (c2*d1 + r*(b2+v2*r))//v1 
    return type(f1)((a3, b3, c3))


def secure_shanks_compose(f1, f2):
    """Secure protocol for composition of positive definite forms.

    Originally by Shanks 'Class Number, a theory of factorization, 
    and genera', Proc. Symp. in Pure Maths, 1969.
    Taken from Coh93, Algorithm 5.4.7.
    """
    # Step 1
    a1, b1, c1, a2, b2, c2 = mpc.if_else(
        f1[0]>f2[0], 
        [f2[0], f2[1], f2[2], f1[0], f1[1], f1[2]], 
        [f1[0], f1[1], f1[2], f2[0], f2[1], f2[2]]
        )

    s = (b1 + b2)/2
    n = b2 - s

    # Step 2: First Euclidean step
    # Skip case distinction (if a1 | a2 as in Coh93) in the oblivious case.
    # Case distinction is pure for performance reasons.
    d, u, v = secure_xgcd(a2, a1)
    y1 = u

    # Step 3: Second Euclidean step
    # Skip case distinction (if d | s as in Coh93) in the oblivious case.
    d1, u, v = secure_xgcd(s, d)
    x2 = u 
    y2 = -v

    # Step 4: Compose
    v1 = a1/d1
    v2 = a2/d1
    # r = (y1*y2*n - x2*c2) % v1
    r_prep = (y1*y2*n - x2*c2)
    _, r = secure_division(r_prep, v1)
    b3 = b2 + 2*v2*r 
    a3 = v1*v2 
    c3 = (c2*d1 + r*(b2+v2*r))/v1 
    return type(f1)((a3, b3, c3))


def _compose_with_self(f, n):
    def composed(arg):
        for _ in range(n):
            arg = f(arg)
        return arg

    return composed


def number2ideal(n, D):
    def check_a(a, D):
        try:
            # print(legendre(D,a) ==1)
            # print(a % 8 in {3, 5, 7})
            return not (gmpy2.legendre(D, a) == 1 and a % 8 in {3, 5, 7})
        except:
            return True

    # print(f"{n=}")
    a = max(2, n - 1)
    # print(f"a at beginning {a=}")
    # while legendre(D, a) != 1 and a % 8 in {3, 5, 7}:
    while check_a(a, D):
        a = int(gmpy2.next_prime(a))
        # print(f"next prime {a=}")
    if a % 4 == 3:
        # b = D**int(((a+1)/4)) % a
        b = D ** ((a + 1) // 4) % a
    else:
        # if D**int(((a-1)/4)) % a == 1:
        #     b = D**int(((a+3)/8)) % a
        # else:
        #     b = 2*D*(4*D)**int(((a-5)/8)) % a
        if D ** ((a - 1) // 4) % a == 1:
            b = D ** ((a + 3) // 8) % a
        else:
            b = 2 * D * (4 * D) ** ((a - 5) // 8) % a

    if D % 2 != b:
        b = a - b
    return (a, b), n - a


def ideal2form(ideal, D):
    a = ideal[0]
    b = ideal[1]
    # f = (a, b, int((b**2 - D)/(4*a)))
    f = (a, b, (b ** 2 - D) // (4 * a))
    return reduce(f)


def number2form(n, D):
    ideal, distance = number2ideal(n, D)
    f = ideal2form(ideal, D)
    return f, distance


def ideal2number(ideal, distance):
    a, b = ideal
    return a + distance


def form2ideal(f):
    a, b = f[0], f[1]
    return (abs(a), b)


def forminverse(f):
    group = type(f)
    return group((f[0], -f[1], f[2]))


def idealinverse(ideal):
    return (ideal[0], -ideal[1])


def form2number(f, distance):
    ideal = form2ideal(f)
    return ideal2number(ideal, distance)


def principal_form(disc):
    """Construct principal form for given discriminant.

    Follows Def. 5.4 from `Binary quadratic forms` by Lipa Long, 2019:
    https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf
    """
    assert disc % 4 == 0 or disc % 4 == 1
    k = disc % 2
    f = (1, k, (k ** 2 - disc) // 4)
    return f


def create_generator_of_subgroup(discriminant):
    """'Generator' as per Chia VDF competition.

    See: https://www.chia.net/2018/11/07/chia-vdf-competition-guide.en.html
    Note: This element generates a cyclic subgroup for given discriminant,
    not per se the entire group.
    """
    if (1 - discriminant) % 8 == 0:
        g = (2, 1, (1 - discriminant) // 8)
    else:
        g = None
    return g


def find_fundamental_discriminant(bit_length):
    """Find fundamental discriminant with additional property: 8 | 1 - discriminant.

    Delta = 1 mod 4 and Delta is square-free, or,
    Delta = 0 mod 4, Delta/4 = 2 or 3 mod 4 and Delta/4 is square-free.
    Fundamental discriminants are those values which are discriminants of quadratic fields.
    Requirement 8 | (1 - discriminant) necessary for create_generator_subgroup method.
    """
    p = gmpy2.next_prime(1 << bit_length - 1)
    while (-p % 4) != 1 or (1 - -p) % 8 != 0:
        p = gmpy2.next_prime(p)
    return int(-p)


def prime_form(a, grp):
    """For prime a, take b square root of discriminant mod 4a.

    Algorithm 3.3 from Buchmann, Vollmer 'Binary Quadratic Forms' 2007.
    """
    assert gmpy2.is_prime(a)
    disc = grp.discriminant

    # Take b square root of disc mod 4a, a prime.
    if a % 4 == 3:
        b = disc ** ((a + 1) // 4) % a
    else:
        if disc ** ((a - 1) // 4) % a == 1:
            b = disc ** ((a + 3) // 8) % a
        else:
            b = 2 * disc * (4 * disc) ** ((a - 5) // 8) % a
    if disc % 2 != b:
        b = a - b

    return grp((a, b, (b ** 2 - disc) // (4 * a)))


def _kronecker(a, b):
    """Implements Legendre symbol and (m/p) for p = 2

    See BV07 Definition 3.4.3.
    """
    if b==2:
        if a % 2 == 0:
            return 0
        else:
            return (-1)**((a**2-1)//8)
    else:
        return gmpy2.legendre(a, b)


def generating_system(grp):
    """Based on Algorithm 9.1 from BV07.

    Time: O(|group.discriminant|^(1/2+o(1))), see Section 9.6.
    BV07: "Based on an idea of H.W. Lenstra. It is the fastest 
    known deterministic class number algorithm."

    Tested based on examples from:
        * BV07, Section 9.6.3
        * https://math.stackexchange.com/questions/2618232/class-group-of-a-field (requires group.nucomp=True)
    """
    disc = grp.discriminant

    def c2(disc):
        # See BV07 Prop. 9.5.1
        return int(sqrt(abs(disc) / 3))

    def c3(disc):
        # See BV07 Prop. 9.5.3
        return 6*int(abs(log2(abs(disc)))**2)

    def primes(start=2, stop=c2(disc)):
        n = start
        while gmpy2.next_prime(n) <= stop:
            n = gmpy2.next_prime(n)
            yield int(n)

    def update_skip_set(skip_set):
        # See BV07 Section 9.6.2, p. 199
        # Remove p if there exists...
        for p in p_set:
            # ... a prime form f = (a,b,c) with a == p
            skip_set.union(set([a[0] for a in prime_forms if a[0] == p]))
            # ... reduced form f = (a, b, p) with b <= 0
            skip_set.union(set([a[2] for a in grp_set if a[1] <= 0]))
            # ... reduced form f = (a, b, c) with p=a-b+c en 0<=2a-b<=p
            for a in grp_set:
                if a[0] - a[1] + a[2] == p:
                    skip_set.add(p)
                if 0 <= 2 * a[0] - a[1] and 2 * a[0] - a[1] <= p:
                    skip_set.add(p)
        return skip_set

    p_set = []
    
    # c = c2(disc)  # Does not work for tiny discriminant = -23
    c = c3(disc)
    # See BV07 Section 9.6.1 and eq (8.16) for definition of P set.
    p_set = [p for p in primes(1, c) if _kronecker(disc, p) != -1]

    skip_set = set()
    grp_set = set()
    gen_set = set()
    prime_forms = set()
    rotor = ["-", "\\", "|", "/", "-", "\\", "|", "/"]

    for p in p_set:
        if not p in skip_set:
            f = prime_form(p, grp)
            prime_forms.add(f.value)
            # Add form values to sets, because form instances are not hashable.
            if f.value not in grp_set:
                gen_set.add(f.value)
            fnew = f
            e = 1
            while not fnew.value in grp_set:
                grp_set.add(fnew.value)
                e += 1
                fnew = reduce_form(f ^ e)
                print(f"Calculating generating set for {grp}: {rotor[(e//1000) % 8]}", end="\r")
            # Update set P; P = p_set \ skip_set in this implementation.
            skip_set = update_skip_set(skip_set)
    return list(map(grp, grp_set)), list(map(grp, gen_set))
