"""Create secure (secret-shared) types for finite groups.

This module provides methods and classes to construct secure groups
from finite groups.
Operators are inherited from non-secure finite groups and then
overloaded via secure group constructors.
SecureGroup() is the main entrypoint and exposes an optional sec_arithm
argument for additional secure arithmetic methods (or overloading).

"""

import asyncio
import logging

from mpyc.finfields import FiniteFieldElement
from mpyc.runtime import mpc
from mpyc.sectypes import SecureObject, SecureFiniteField, SecureInteger
from mpyc.thresha import _recombination_vector
from mpyc.seclists import seclist
import mpyc.mpctools as mpctools
import sec_groups.fingroups as fg
import sec_groups.classgrps as classgrps
from sec_groups.tools.repeat import (
    repeat_public_base_secret_output,
    secure_repeat_public_base_public_output,
    secure_repeat_secret_base_secret_output,
    secgrp_if_else,
)


logger_sg = logging.getLogger("SecureGroup")
logger_sg.setLevel(logging.INFO)


class SecureGroupElement(SecureObject):
    """Abstract base class for secure groups. """

    #    __slots__ = ()  # NB: to avoid "multiple bases have instance lay-out conflict" cf. SecureObject

    group = None
    sectype_of_value = None
    identity = None
    sec_arithm = None

    @property
    def share(self):
        return self.value

    @share.setter
    def share(self, value):
        self.value = value

    @classmethod
    def _input(cls, x, senders):
        """Input a list x of secure group elements. See mpc.input()."""
        is_tuple = isinstance(x[0].value, tuple)
        if is_tuple:
            r = len(x[0].value)
            x = [_ for a in x for _ in a.value]
        else:
            x = [a.value for a in x]
        shares = mpc.input(x, senders)
        if is_tuple:
            shares = [[s[i : i + r] for i in range(0, len(s), r)] for s in shares]
        return [[cls(a) for a in s] for s in shares]

    @classmethod
    async def _output(cls, x, receivers, threshold):
        """Output a list x of secure group elements.  See mpc.output()."""
        is_tuple = isinstance(x[0].value, tuple)
        if is_tuple:
            r = len(x[0].value)
            x = [_ for a in x for _ in a.value]
        else:
            x = [a.value for a in x]
        y = await mpc.output(x)
        if is_tuple:
            y = [y[i : i + r] for i in range(0, len(y), r)]
        return list(map(cls.group, y))


class SecureSymmetricGroupElement(SecureGroupElement, fg.SymmetricGroupElement):
    """Common base class for secure Symmetric Group elements."""

    def __init__(self, value=None):
        """Ensure all coefficients of value are of secure field type."""
        n = self.group.degree
        if value is None:
            value = [None] * n
        elif isinstance(value, self.group):
            value = value.value
        else:
            if not (isinstance(value, (tuple, list)) and len(value) == n):
                raise ValueError(f"tuple/list of length {n} required")

        secfld = self.sectype_of_value
        self.value = tuple(secfld(a) if not isinstance(a, secfld) else a for a in value)

    def set_share(self, value):
        for a, b in zip(self.value, value):
            a.set_share(b.share)

    def operation(p, q):
        """First p then q."""
        group = type(p)
        q = seclist(q.value)
        return group(tuple(q[j] for j in p.value))

    def inverse(p):
        group = type(p)
        n = len(p.value)
        q = seclist(p.value)  # use p.value as dummy of the right type
        for i in range(n):
            q[p.value[i]] = i
        return group(tuple(q))

    def equality(p, q):  # return type is self.sectype_of_value
        return seclist(p.value) == seclist(q.value)


class SecureQRElement(SecureGroupElement, fg.QuadraticResidueElement):
    """Common base class for secure Quadratic Residue elements."""

    def __init__(self, value=None):
        """Ensure value is of secure field type."""
        if isinstance(value, self.group):
            value = value.value
        secfld = self.sectype_of_value
        if not isinstance(value, secfld):
            value = secfld(value)
        self.value = value

    def set_share(self, v):
        self.share.set_share(v.share)


class SecureEllipticCurveElement(SecureGroupElement, fg.EllipticCurveElement):
    """Common base class for secure elliptic curve elements."""

    # __slots__ = "value"
    #    __slots__ = ()

    def __init__(self, value=None):
        """Ensure all coefficients are of secure field type.

        Enforce value is a tuple.
        """
        n = len(self.group.identity.value)
        if value is None:
            value = [None] * n
        elif isinstance(value, self.group):
            value = value.value
        else:
            if not (isinstance(value, (tuple, list)) and len(value) == n):
                raise ValueError(f"tuple/list of length {n} required")

        secfld = self.sectype_of_value
        self.value = tuple(secfld(a) if not isinstance(a, secfld) else a for a in value)

    def set_share(self, value):
        for a, b in zip(self.value, value):
            a.set_share(b.share)

    def operation(a, b):
        """Suggest to use secure arithmetic mixin at construction."""
        raise TypeError(
            "`operation` needs to be explicitly defined at construction. Pass elliptic curve arithmetic to constructor."
        )

    @property
    def curve(self):
        return self.group.curve

    @property
    def coord(self):
        return self.group.coord

    @property
    def arithm(self):
        raise TypeError(
            "Arithmetic of non-sectype group accessible via `.group.arithm` attribute. Secure arithmetic via `.sec_arithm`."
        )


class SecureFormClassGroupElement(SecureGroupElement, fg.FormClassGroupElement):
    """Common base class for secure elliptic curve elements."""

    #    __slots__ = ()

    def __init__(self, value=None):
        """Ensure all coefficients are of secure type.

        Enforce value is a tuple.
        """
        n = 3
        if value is None:
            value = [None] * n
        elif isinstance(value, self.group):
            value = value.value
        else:
            if not (isinstance(value, (tuple, list)) and len(value) == n):
                raise ValueError(f"tuple/list of length {n} required")

        sectype = self.sectype_of_value
        self.value = tuple(
            sectype(a) if not isinstance(a, sectype) else a for a in value
        )

    def set_share(self, value):
        for a, b in zip(self.value, value):
            a.set_share(b.share)

    def operation(a, b):
        c = classgrps.secure_shanks_compose(a, b)
        if a.group.auto_reduce:
            c = classgrps.secure_binary_reduce(c)
        return c

    def operation2(a):
        c = classgrps.secure_square(a)
        if a.group.auto_reduce:
            c = classgrps.secure_binary_reduce(c)
        return c

    def inverse(a):
        c = classgrps.forminverse(a)
        if a.group.auto_reduce:
            c = classgrps.secure_binary_reduce(c)
        return c

    def reduce_form(a):
        return classgrps.secure_binary_reduce(a)

    def equality(a, b):
        a_red = a.reduce_form()
        b_red = b.reduce_form()
        # Only test a, b coeffs, because c is redundant.
        v0 = a_red.value[0] == b_red.value[0]
        v1 = a_red.value[1] == b_red.value[1]
        return v0 * v1


def _parent_and_sectype(group):
    """Look up secure parent corresponding to group.

    Also return sectype used to represent coefficients.
    """
    if issubclass(group, fg.QuadraticResidueElement):
        parent = SecureQRElement
        field = group._field
        sectype_of_value = mpc.SecFld(field.order)
    elif issubclass(group, fg.EllipticCurveElement):
        parent = SecureEllipticCurveElement
        field = group.field
        sectype_of_value = mpc.SecFld(field.order)
    elif issubclass(group, fg.SymmetricGroupElement):
        parent = SecureSymmetricGroupElement
        field = group._field
        sectype_of_value = mpc.SecFld(field.order)
    elif issubclass(group, fg.FormClassGroupElement):
        parent = SecureFormClassGroupElement
        logger_sg.debug(
            "WARNING: crude bound for secint.bit_length, multiple of bit-length(discriminant)."
        )
        """TODO: check required secint length (see TODO list at bottom). 
        Bound depends on room required for secure gcd algorithm and
        on size of non-reduced forms. 
        """
        SAFETY = 4
        sectype_of_value = mpc.SecInt(
            SAFETY * group.bit_length
        )  # 4x set to enable secure_binary_reduce()
    else:
        raise NotImplementedError
    return parent, sectype_of_value


def SecureGroup(group, sec_arithm=None):
    """Secure group constructor.

    Args:
        group (FiniteGroupElement): Group to create secure group from.
        sec_arithm (class): Arithmetic for optional operator overloading.

    Returns:
        class: secure equivalent of group
    """
    parent, sectype_of_value = _parent_and_sectype(group)

    name = f"SecGrp({group.__name__})"
    #    body = {'__slots__': ()}
    body = {}
    # Mix in arithmetic if sec_arithm given, or group.arithm is oblivious.
    if sec_arithm:
        parents = (sec_arithm, parent)   
    elif hasattr(group, 'arithm') and group.arithm.oblivious:
        parents = (group.arithm, parent)
    else:
        parents = (parent,)
    sectype = type(name, parents, body)
    sectype.__doc__ = "Class of secret-shared finite group elements."
    sectype.group = group
    sectype.sectype_of_value = sectype_of_value
    sectype.sec_arithm = sec_arithm
    sectype.identity = sectype(group.identity)
    globals()[
        name
    ] = sectype  # NB: exploit (almost) unique name dynamic SecureGroup type
    return sectype


def secure_repeat(a, x, sec_grp=None, public_out=False, exp_as_bits=False):
    """Wraps several secure repeat methods.

    Methods maintained in module tools/repeat.

    Args
        a           base element
        x           exponent element
        sec_grp     defines output type

    Flags
        public_out  set if outputs is public
        exp_as_bits set if exponent is given as list of bits
    """
    assert isinstance(a, fg.FiniteGroupElement) or isinstance(
        a[0], fg.FiniteGroupElement
    )
    if not isinstance(a, list):
        group = type(a)
    else:
        group = type(a[0])

    # Case: Public exponent, secret output.
    if isinstance(x, (int, FiniteFieldElement)):
        return group.repeat(a, x)

    # Case: Public base, secret exponent, public output
    elif isinstance(x, (SecureFiniteField, SecureInteger)) and public_out == True:
        return secure_repeat_public_base_public_output(a, x)

    # Case: Public base, secret output.
    elif (
        not isinstance(a, SecureObject)
        and isinstance(x, (SecureFiniteField, SecureInteger))
        and not public_out
    ):

        if sec_grp is not None:
            assert issubclass(sec_grp, SecureGroupElement)
            return repeat_public_base_secret_output(a, x, sec_grp)
        else:
            logger_sg.debug(
                "WARNING/DEBUG: No secure group provided to secure_repeat. \
                Assuming secure arithmetic is well-defined by parent class."
            )
            if hasattr(group, "arithm"):
                sec_arithm = group.arithm
                sec_grp = SecureGroup(group, sec_arithm)
            else:
                sec_grp = SecureGroup(group)
            return repeat_public_base_secret_output(a, x, sec_grp)

    # Case: Other; Assume secret base, secret exponent, secret output
    else:
        if not exp_as_bits:
            x = mpc.to_bits(x)
        return secure_repeat_secret_base_secret_output(a, x, sec_grp)


# TODOs
# 1: Expand secure_repeat for non-prime order.
# 1: Define tighter bound for class group bitlength of secint (especially if automatic reduction is applied after composition); Consider converting sectypes temporarily during xgcd protocol.; Check Buchmann, Vollmer: 'Binary Quadratic Forms', Lemma 5.6.1 bound 2max(a,c)/sqrt(delta)
# 1: Use parent in the SecureGroup constructor for type checking? (if curve is an Edwards curve?)
# 2: Enforce that is_additive, is_multiplicative are checked at sectype.group.attribute level (not instance.attribute level)
# 1: Consider making .x, .y, .z attribute access part of SecureGroup constructor, instead of via (mandatory) passing of (secure) arithmetic.

