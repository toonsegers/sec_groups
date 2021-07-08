import asyncio

from mpyc.runtime import mpc
from mpyc.finfields import FiniteFieldElement, GF
from mpyc.runtime import mpc
from mpyc.sectypes import SecureObject, SecureFiniteField, SecureInteger
from mpyc.thresha import _recombination_vector
import mpyc.mpctools as mpctools
import sec_groups.fingroups as fg
from operator import mul


def secure_mont_ladder(a, n, bitlen_n, sectype):
    """Compute a^[n]->[a^n] and return sectype type.

    Requires:
        n is a (secure) prime field element or (secure) int.
        n > 0
        a**x not allowed to overflow, i.e., > sectype.field.modulus

    Implements Montgomery's ladder https://en.wikipedia.org/wiki/Exponentiation_by_squaring
    """
    x1 = sectype(1)
    x2 = sectype(a)

    n_bits = mpc.to_bits(n, bitlen_n)
    l = len(n_bits)

    for i in range(l - 2, -1, -1):
        # print("in secure_mont_ladder, i=", i)
        x1, x2 = mpc.if_else(
            mpc.sgn(n_bits[i], l=1, EQ=True), [x1 * x1, x1 * x2], [x1 * x2, x2 * x2]
        )
    return x1


def exp_by_squaring(x, n):
    """Assumes n>=0

    See: https://en.wikipedia.org/wiki/Exponentiation_by_squaring
    """
    if n == 0:
        return 1
    if n % 2 == 0:
        return exp_by_squaring(x * x, n // 2)
    else:
        return x * exp_by_squaring(x * x, (n - 1) / 2)


def exp_by_squaring_bits(x, n):
    """Assumes n>=0

    Requires:
        n given in bits

    See: https://cp-algorithms.com/algebra/binary-exp.html
    """
    # if all(v == 0 for v in n):
    #     return 1

    res = 1
    for i in reversed(n):
        if i % 2 == 1:
            res = res * x
        x = x * x
    return res


def secure_pow(a, n, sectype):
    """Compute a^[n]->[a^n] and return sectype type.

    Requires:
        n   exponent >0, given as list of secure bits
        a**x not allowed to overflow, i.e., > sectype.field.modulus
        a   base is abelian.
    """
    if not isinstance(a, SecureObject):
        a = sectype(a)

    a_powlist = [a]
    for i, _ in enumerate(n[:-1]):
        a_powlist += [a_powlist[-1] ** 2]

    b = sectype(1)
    for i, _ in enumerate(n):
        b = mpc.if_else(n[i], b * a_powlist[i], b)
    return b


def secgrp_if_else(c, x, y):
    """Secure selection based on condition c between x and y.

    Requires: Value of c < order of new sectype for convert to work.

    This is the case when modulus of exponent field corresponds to
    the order of the group. (Then, the order of the group's base field
    is greater than the order of the exponent field.)
    """
    cc = mpc.convert(c, x.sectype_of_value)

    assert type(cc) == x.sectype_of_value
    sec_grp = type(x)
    if isinstance(x.value, (tuple, list)):
        z = [cc * (x.value[i] - y.value[i]) + y.value[i] for i, _ in enumerate(x.value)]
        z = tuple(z)
    else:
        z = cc * (x.value - y.value) + y.value
    return sec_grp(z)


def secure_repeat_secret_base_secret_output(a, n, sectype):
    """Compute a^[n]->[a^n] and return sectype type.

    Requires:
        n   exponent >0, given as list of secure bits
        a**x not allowed to overflow, i.e., > sectype.field.modulus
        a   base is abelian.
    """
    if hasattr(sectype, "group"):
        assert sectype.group.is_abelian
        op = sectype.operation
        unit = sectype(sectype.group.identity)
        if_else = secgrp_if_else
    else:
        unit = sectype(1)
        op = mul
        if_else = mpc.if_else

    if not isinstance(a, SecureObject):
        a = sectype(a)

    a_powlist = [a]
    for i, _ in enumerate(n[:-1]):
        a_powlist += [op(a_powlist[-1], a_powlist[-1])]

    # b = sectype(1)
    b = unit
    for i, _ in enumerate(n):
        b = if_else(n[i], op(b, a_powlist[i]), b)
    return b


@mpc.coroutine
async def repeat_public_base_secret_output(a, x, sec_grp):
    """Compute a^[x]->[a^x] and return sec_grp type.

    Requires:
        a in prime order group, x in Z.
        a is not a secure element.
        x is a (secure) prime field element or (secure) int.
    """
    assert isinstance(a, fg.FiniteGroupElement) and not isinstance(a, SecureObject)
    assert isinstance(x, (SecureFiniteField, SecureInteger))
    assert sec_grp.group.order == type(a).order

    await mpc.returnType(sec_grp)
    m = len(mpc.parties)
    lambda_i = _recombination_vector(x.field, range(1, m + 1), 0)[mpc.pid]
    x_i = await mpc.gather(x)
    c_i = sec_grp.group.repeat(a, int(lambda_i * x_i))
    c = mpc.input(sec_grp(c_i))
    return mpctools.reduce(sec_grp.operation, c)


@mpc.coroutine
async def secure_repeat_public_base_public_output(a, x) -> asyncio.Future:
    """Multi-exponentiation for given base(s) a and exponent(s) x.

    Argument a is a group element, or a list of group elements.
    Argument x is secure number (prime field element, or integer), or a list of secure numbers.
    The dimensions of a and x should match.

    Moreover, the order of the group should be prime and should match the order the field over
    which the secure numbers are shared. This ensures that "Shamir in the exponent" can be used.
    """
    if not isinstance(a, list):
        a, x = [a], [x]
    group = type(a[0])
    field = GF(modulus = group.order)

    m = len(mpc.parties)
    lambda_i = _recombination_vector(field, range(1, m + 1), 0)[mpc.pid]

    x_i = await mpc.gather(x)
    e_i = [int(lambda_i * s_i) for s_i in x_i]
    c_i = mpctools.reduce(group.operation, map(group.repeat, a, e_i))
    c = await mpc.transfer(c_i)
    return mpctools.reduce(group.operation, c)
