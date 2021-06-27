import decimal
from math import log, floor
from mpyc import gmpy as gmpy2


def find_safe_primes(modulus_bitlength=256):
    """Find primes q = 2p + 1, where q has given bit length. """
    p, q = find_primes_for_schnorr(l=modulus_bitlength, r=2, q=5)
    return p, q


def find_primes_for_schnorr(l=256, r=2, q=5):
    """Find smallest prime q of bit length at least l, >q,
    with p=q r + 1, with both p, q prime.
    """
    # TODO: make sure p is a Blum prime (p mod 4 = 3); e.g., needed for QR groups
    if l <= 256:
        assert r % 2 == 0, "Parameter r should be even to find Schnorr primes."
        if l == 2:
            p = 3
        else:
            q = 2 ** l
            q = int(gmpy2.next_prime(q))
            p = (q - 1) // r
            while not (gmpy2.is_prime(p) and gmpy2.is_prime(q)):
                q = int(gmpy2.next_prime(q))
                p = (q - 1) // r
    elif r == 2:
        ike_b_k_options = {
            768: 149686,
            1024: 129093,
            2048: 124476,
            3072: 1690314,
            4096: 240904,
            6144: 929484,
            8192: 4743158,
        }

        if l in ike_b_k_options:
            b, k = l, ike_b_k_options[l]
        else:
            raise NotImplementedError(
                f"Specific prime with bitlength {l} not supported at the moment."
                f"Try l=2048 or 4096, for example."
            )
        q = _find_ike_prime(b, k)
        p = (q - 1) // 2
    else:
        raise NotImplementedError(
            f"Generating primes of bitlength l = {l} and r = {r} (as in q = p r + 1) not supported."
        )
    assert gmpy2.is_prime(p)
    assert gmpy2.is_prime(q)
    return p, q


def _make_pi_generator(max_i=1000):
    # From here: https://rosettacode.org/wiki/Pi#Python
    q, r, t, k, n, l = 1, 0, 1, 1, 3, 3
    i = 0
    while i < max_i:
        if 4 * q + r - t < n * t:
            i += 1
            yield n
            nr = 10 * (r - n * t)
            n = ((10 * (3 * q + r)) // t) - 10 * n
            q *= 10
            r = nr
        else:
            nr = (2 * q + r) * l
            nn = (q * (7 * k) + 2 + (r * l)) // (t * l)
            q *= k
            t *= l
            l += 2
            k += 1
            n = nn
            r = nr


def _make_pi(l=1024):
    # TODO: recipe voor pi in https://docs.python.org/3/library/decimal.html
    n = int(floor(log(2 ** l, 10)))
    context = decimal.Context(prec=n)
    decimal.setcontext(context)
    pi_array = []
    for i in _make_pi_generator(n):
        pi_array.append(str(i))
    pi_array = pi_array[:1] + ["."] + pi_array[1:]
    pi_string = "".join(pi_array)
    pi_big_decimal = decimal.Decimal(pi_string)
    return pi_big_decimal


def _find_ike_prime(b, k, fixedbits=64):
    # From here: https://kivinen.iki.fi/primes/
    # See also: https://link.springer.com/chapter/10.1007/978-3-540-24582-7_17
    pi = _make_pi(b)
    epi = floor(pi * 2 ** (b - 2 * fixedbits - 2)) + k
    prime_candidate = 2 ** (b) - 2 ** (b - fixedbits) - 1 + epi * 2 ** fixedbits
    return prime_candidate
