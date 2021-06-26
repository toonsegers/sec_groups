""" Helper functions to encode int or field element to group element.

Decoders return int or MPYC securet type (subclass of SecureObject)
"""

import logging

from mpyc.runtime import mpc
from mpyc.gmpy import legendre, next_prime, is_prime
from mpyc.sectypes import SecureObject

from sec_groups.secgroups import SecureQRElement, SecureEllipticCurveElement
from sec_groups.fingroups import (
    QuadraticResidueElement,
    EllipticCurveElement,
    FormClassGroupElement,
    affine_tuple_to_coord,
)
from sec_groups.classgrps import _kronecker, prime_form

logger_enc = logging.getLogger("Encoding")
logger_enc.setLevel(logging.DEBUG)


# WIP; consider jumping through range using hash(i)
def _next_candidate(candidate, group, a_shifted, k):
    if issubclass(group, QuadraticResidueElement):
        new_candidate = candidate + 1
        assert new_candidate < a_shifted + 2**k
        return new_candidate
    elif issubclass(group, EllipticCurveElement):
        new_candidate = candidate + 1
        assert new_candidate < a_shifted + 2**k
        return new_candidate
    elif issubclass(group, FormClassGroupElement):
        new_candidate = int(next_prime(candidate))
        assert new_candidate < a_shifted + 2**k
        return new_candidate
        # return int(next_prime(candidate))
    else:
        raise NotImplementedError


def _edwards_ysquared(x, group):
    fld = group.field
    return (group.curve.c * group.curve.c - group.curve.a * x * x) / (
        fld(1) - group.curve.c * group.curve.c * group.curve.d * x * x
    )

def _weierstrass_ysquared(x, group):
    fld = group.field 
    return x**3 + group.curve.a * x + group.curve.b


def _check(candidate, group):
    if issubclass(group, QuadraticResidueElement):
        return legendre(candidate, group._field.modulus)
    elif issubclass(group, EllipticCurveElement):
        if group.curve.name in ["ED25519", "ED448"]:  # TODO: check if Edwards by creating subclass for Edwards/Weierstrass/etc.
            ysquared = _edwards_ysquared(candidate, group)
            return ysquared.is_sqr()
        else:
            ysquared = _weierstrass_ysquared(candidate, group)
            return ysquared.is_sqr()
    elif issubclass(group, FormClassGroupElement):
        # See Buchmann, Vollmer 2007 "Binary Quadratic Forms", Proposition 3.4.5
        # Prime form (a, b, c) where a prime exists when following conditions are met.
        return is_prime(candidate) and _kronecker(group.discriminant, candidate) != -1
    else:
        raise NotImplementedError


def _map_to_group(candidate, group):
    if issubclass(group, QuadraticResidueElement):
        return group(candidate)
    elif issubclass(group, EllipticCurveElement):
        if group.curve.name in ["ED25519", "ED448"]:
            ysquared = _edwards_ysquared(candidate, group)
            y_coord = ysquared.sqrt()
            return affine_tuple_to_coord(group, (candidate, y_coord))            
        else:
            ysquared = _weierstrass_ysquared(candidate, group)
            logger_enc.debug("Calculating square root... (may take long for large groups)")
            y_coord = ysquared.sqrt()
            logger_enc.debug("Calculating square root... done.")
            return affine_tuple_to_coord(group, (candidate, y_coord))            
    elif issubclass(group, FormClassGroupElement):
        return prime_form(candidate, group)
    else:
        raise NotImplementedError


def encode_v3(a, group, k=8, encode_zero=False):  # TODO: add flag return_zero_encoding = False
    """k security parameter. (bitlength of range)
    """
    logger_enc.debug(f"Enter encode with k={k}, encode_zero={encode_zero}")
    if group.order is not None:
        assert a in range(group.order // 2**k), "Input out of range."
    else:
        # Assume group order is large enough.
        pass
    success = False
    a_shifted = a << k
    candidate = a_shifted

    while candidate < a_shifted + 2**k:
        # search
        candidate = _next_candidate(candidate, group, a_shifted, k)
        logger_enc.debug(f"Inside encode while loop; candidate={candidate}")

        # check candidate
        success = _check(candidate, group)
        if encode_zero:
            success &= _check(candidate & (2**k-1), group)

        # map to group
        if success:
            a_enc = _map_to_group(candidate, group)
            if encode_zero:
                z_enc = _map_to_group(candidate & (2**k-1), group)
            break
    if not success:
        raise ValueError("No encoding found. Consider increasing bitshift size.")
    
    if encode_zero:
        logger_enc.debug(f"Encoding returns zero encoding as well.")
        return a_enc, z_enc
    else:
        return a_enc


def decode_v3(a_enc, k=8):
    # If a_enc contains two elements, assume encoding of zero is passed along.
    encode_zero = False
    if isinstance(a_enc, list) and len(a_enc)==2:
        a_enc, z_enc = a_enc[0], a_enc[1]
        logger_enc.debug(f"z_enc={z_enc}.")
        encode_zero = True

    if isinstance(a_enc.value, (tuple, list)):
        a = a_enc.value[0]
        if encode_zero:
            z = z_enc.value[0]
    else:
        a = a_enc.value
        if encode_zero:
            z = z_enc.value

    if isinstance(a, SecureObject):
        if encode_zero:
            logger_enc.debug(f"Decoding secure object. Use encoding of zero for efficient decoding.")
            a = a - z
            a = a / (2 ** k)
        else:
            logger_enc.debug(f"Decoding secure object. Trunc last k={k} bits, using mpc.to_bits().")
            a = mpc.to_bits(a)[k:]
            a = mpc.from_bits(a)
    else:
        a = int(a) >> k
    return a


def encode_to_schnorr(m, schnorr_grp, bitshift=8):
    assert m in range(
        schnorr_grp.order // (2 ** bitshift)
    ), "Message out of encoding range."
    for i in range(1, 2 ** bitshift):
        m_enc, i = _encode_to_schnorr_return_i(
            m, schnorr_grp, i_start=i, bitshift=bitshift
        )
        z_enc, j = _encode_to_schnorr_return_i(
            0, schnorr_grp, i_start=i, bitshift=bitshift
        )
        if i == j:
            return m_enc, z_enc
    raise ValueError("Encoding parameter (bitshift={bitshift}) to small for successful encoding.")


def _encode_to_schnorr_return_i(m, schnorr_grp, i_start=1, bitshift=8):
    assert m in range(
        schnorr_grp.order // (2 ** bitshift)
    ), "Message out of encoding range."
    fld = schnorr_grp._field
    in_group = False
    m_shifted = m << bitshift
    for i in range(i_start, 2 ** bitshift):
        m_encode = m_shifted + fld(i)
        in_group = m_encode ** schnorr_grp.order == 1
        if in_group:
            break
    if not in_group:
        raise ValueError("No encoding found. Consider increasing bitshift size.")
    return schnorr_grp(m_encode), i


def decode_from_plain_schnorr(m_enc, bitshift=8):
    m = int(m_enc) >> bitshift
    return m


async def decode_from_schnorr(m_enc, z_enc=None, bitshift=8):
    if isinstance(m_enc, tuple):
        m_enc, z_enc = m_enc

    if isinstance(m_enc, SecureQRElement):
        assert (
            z_enc is not None
        ), f"No zero encoding passed to decoder. Required for decoding in MPC."
        m_type = type(m_enc)
        m_val = m_enc.value - z_enc.value
        m = m_type(m_val)
        m = m * (m_enc.group._field(2 ** bitshift)).reciprocal()
        m = m.value  # Convert to SecFld (instead of Secure QR element)
    elif isinstance(m_enc, (QuadraticResidueElement, int)):
        m = decode_from_plain_schnorr(m_enc, bitshift=8)
    else:
        raise NotImplementedError
    return m


def encode_on_curve(m, curve, bitshift=8):
    assert m in range(curve.order // (2 ** bitshift)), "Message out of encoding range."
    for i in range(1, 2 ** bitshift):
        M, i = _encode_on_curve_return_i(m, curve, i_start=i, bitshift=bitshift)
        Z, j = _encode_on_curve_return_i(0, curve, i_start=i, bitshift=bitshift)
        if i == j:
            return M, Z
        else:
            i = max(i, j)


def _encode_on_curve_return_i(m, curve, i_start=1, bitshift=8):
    """ Encode message m on curve (case: K = GF(q), q = p^n with n=1)

    From [Koblitz1987]: "Suppose that our plaintexts are integers m,
    0 < m < p/1000 - 1. We try appending three digits to m until we
    obtain an x, 1000m < x < 1000(m + 1) < p such that f(x) is a square
    in GF(p)."
    We look for x such that f(x,y) is a point on the Edwards curve.
    """
    def _edwards_ysquared(x, group):
        return (group.curve.c * group.curve.c - group.curve.a * x * x) / (
            fld(1) - group.curve.c * group.curve.c * group.curve.d * x * x
        )

    fld = curve.field

    if isinstance(m, int):
        m = fld(m)
    m_extra = m * (2 ** bitshift)
    for i in range(
        i_start, 2 ** bitshift
    ):  # skip 0 to avoid trivial encoding of message=0
        m_extra_i = m_extra + fld(i)
        ysquared_m = _edwards_ysquared(m_extra_i, curve)
        if legendre(int(ysquared_m), fld.modulus) == 1:
            y_coord = ysquared_m.sqrt()
            # encoded_pt = EdwardsGroupElement(curve, m_extra_i, y_coord)
            encoded_pt = affine_tuple_to_coord(curve, (m_extra_i, y_coord))
            break
    if not encoded_pt:
        raise ValueError(f"No encoding for message m: {m} on curve found")
    return encoded_pt, i


def decode_from_curve(M_enc, Z_enc=None, bitshift=8):
    # if isinstance(M_enc, SecureObject) or isinstance(M_enc, SecureEdwardsGroup):
    if isinstance(M_enc, SecureEllipticCurveElement):
        assert (
            Z_enc is not None
        ), f"No zero encoding passed to decoder. Required for decoding in MPC."
        M_x = M_enc.x - Z_enc.x
        m = M_x / M_enc.curve.field(2 ** bitshift)
    # elif isinstance(M_enc, Point):
    elif isinstance(M_enc, EllipticCurveElement):
        M_x = M_enc.x
        m = int(M_x) // 2 ** bitshift
    else:
        raise NotImplementedError
    return m


# TODOs
# 1: encode_to_schnorr / encode_to_curve: make uniform API
