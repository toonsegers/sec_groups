"""This module implements finite groups in Python.

A finite group is a set together with a binary operation
optionally written additively or multiplicatively.
"""
# TODO: reconsider use of __slots__ (currently EC types still get a dict)
# TODO: decide whether operation(2), inverse, equality, repeat should be @staticmethod

import functools
import math
from typing import NamedTuple, Any

from mpyc import gmpy as gmpy2
import sec_groups.classgrps as classgrps
import sec_groups.rubiks as rubiks
from mpyc.finfields import GF, PrimeFieldElement, find_prime_root


class FiniteGroupElement:
    """Abstract base class for finite groups.

    Default: @, ~, ^
    Additive: +, -, *
    Multiplicative: *, 1/ (= **-1), **
    """

#    __slots__ = ()  # NB: to avoid "multiple bases have instance lay-out conflict" cf. SecureObject

    order = None
    is_additive = False
    is_multiplicative = False
    identity = None
    is_abelian = False
    is_cyclic = False
    generator = None
    generating_set = None

    def __matmul__(self, other):  # overload @
        group = type(self)
        #        return group(group.operation(self, other))
        return group.operation(self, other)

    def __invert__(self):  # overload ~
        group = type(self)
        #        return group(group.inverse(self))
        return group.inverse(self)

    def __xor__(self, other):  # overload ^
        group = type(self)
        return group.repeat(self, other)

    def __add__(self, other):
        group = type(self)
        if not group.is_additive:
            raise TypeError("group not additive")

        return group.__matmul__(self, other)

    def __neg__(self):
        group = type(self)
        if not group.is_additive:
            raise TypeError("group not additive")

        return group.__invert__(self)

    def __sub__(self, other):
        group = type(self)
        if not group.is_additive:
            raise TypeError("group not additive")

        return group.__matmul__(self, group.__invert__(other))

    def __mul__(self, other):
        group = type(self)
        if group.is_multiplicative:
            return group.__matmul__(self, other)

        if group.is_additive:
            return group.__xor__(self, other)

        raise TypeError("* not defined for group")

    def __rmul__(self, other):
        group = type(self)
        if group.is_multiplicative:
            if group.is_abelian:
                return group.__matmul__(self, other)

            return group.__matmul__(group(other), self)

        if group.is_additive:
            return group.__xor__(self, other)

        raise TypeError("* not defined for group")

    def __truediv__(self, other):
        group = type(self)
        if not group.is_multiplicative:
            raise TypeError("group not multiplicative")

        return group.__matmul__(self, group.__invert__(other))

    def __rtruediv__(self, other):
        group = type(self)
        if not group.is_multiplicative:
            raise TypeError("group not multiplicative")

        if other != 1:
            raise TypeError("only 1/ supported")

        return group.__invert__(self)

    def __pow__(self, other):
        group = type(self)
        if not group.is_multiplicative:
            raise TypeError("group not multiplicative")

        return group.__xor__(self, other)

    def __eq__(self, other):
        group = type(self)
        return group.equality(self, other)

    # instance method, maybe change to static or class method
    def operation(a, b):
        """Return a @ b."""
        raise NotImplementedError

    # instance method, maybe change to static or class method
    def operation2(a):
        """Return a @ a."""
        group = type(a)
        return group.operation(a, a)

    # instance method, maybe change to static or class method
    def inverse(a):
        """Return @-inverse of a (written ~a)."""
        raise NotImplementedError

    # instance method, maybe change to static or class method
    def equality(a, b):
        """Return a == b."""
        raise NotImplementedError

    # instance method, maybe change to static or class method
    def repeat(a, n):
        """Return nth @-power of a (written a^n), for any integer n."""
        if isinstance(n, PrimeFieldElement):
            n = int(n)

        group = type(a)
        if n == 0:
            return group.identity

        if n < 0:
            a = group.inverse(a)
            n = -n
        b = a
        for i in range(n.bit_length() - 2, -1, -1):
            b = group.operation2(b)
            if (n >> i) & 1:
                b = group.operation(b, a)
        return b


class SymmetricGroupElement(FiniteGroupElement):
    """Symmetric group of all permutations of a fixed length.

    Permutations represented as length-n tuples with entries in {0,...,n-1},
    following image map representation.
    """

#    __slots__ = "value"
    degree = None

    def __init__(self, value):
        if len(value) < self.degree:
            # Ensure value is an image map of full length.
            value = list(value) + list(range(len(value), self.degree))
        if isinstance(value, list):
            value = tuple(value)
        self.value = value

    def operation(p, q):
        """First p then q."""
        group = type(p)
        return group(tuple(q.value[int(j)] for j in p.value))

    def inverse(p):
        group = type(p)
        n = len(p.value)
        q = [None] * n
        for i in range(n):
            q[p.value[i]] = i
        return group(tuple(q))

    def equality(p, q):
        return p.value == q.value

    def __repr__(self):
        return str(self.value)


def S(n):
    """Create type for symmetric group Sym_n.

    Representation based on image maps, not cycle notation.
    """
    name = f'Sym({n})'
    Sym = type(name, (SymmetricGroupElement,), {})
#    Sym = type(name, (SymmetricGroupElement,), {'__slots__': ()})
    Sym.order = math.factorial(n)
    Sym.degree = n
    Sym._field = GF(int(gmpy2.next_prime(n*2)))  # TODO: Understand why factor *2 is needed (to patch issue after update of mpyc.to_bits())
    Sym.is_abelian = n <= 2
    Sym.is_cyclic = n <= 2
    Sym.identity = Sym(tuple(range(n)))
    # Generating set: {(0, 1), (0, .., n-1)} in cycle notation. We use image map notation.
    Sym.generating_set = [Sym([1, 0]), Sym(list(range(1, n))+[0])]
    globals()[name] = Sym  # NB: exploit (almost?) unique name dynamic Sym type
    return Sym


class QuadraticResidueElement(FiniteGroupElement):
    """Common base class for quadrtic residue group Z/pZ*, for prime modulus p.

    Quadratic residues represented by GF(p) elements.
    """

#    __slots__ = "value"

    is_multiplicative = True
    is_abelian = True
    is_cyclic = True

    def __init__(self, value):
        if isinstance(value, int):  # TODO: check if we want to catch more cases.
        # if not isinstance(value, self._field):
            value = self._field(value)
        self.value = value

    def operation(a, b):
        group = type(a)
        return group(a.value * b.value)

    def inverse(a):
        group = type(a)
        return group(1 / a.value)

    def equality(a, b):
        return a.value == b.value

    def repeat(a, n):  # override to get more speed
        group = type(a)
        return group(a.value ** int(n))

    def __repr__(self):
        return str(self.value)

    def __int__(self):
        assert isinstance(self.value, PrimeFieldElement)
        return int(self.value)


def QuadraticResidue(modulus):
    """Create a QR type for given modulus."""
    # assert gmpy2.is_prime(modulus >> 1), "Choose QR subgroup Z/pZ* of prime order."  # Assume modulus>>1 is prime to reduce overhead.
    name = f'QR({modulus})'
    QR = type(name, (QuadraticResidueElement,), {})
#    QR = type(name, (QuadraticResidueElement,), {'__slots__': ()})
    QR._field = GF(modulus)
    QR.order = modulus >> 1

    # Add identity and generator to class. Ensure rest of type gets defined above this line.
    QR.identity = QR(1)
    # Euler’s Criterion: a is a quadratic residue if and only if a^((p-1)/2)==1.
    for i in range(2, modulus):
        if gmpy2.legendre(i, modulus) == 1:
            generator = QR(i)
            break
    QR.generator = generator
    globals()[name] = QR  # NB: exploit (almost?) unique name dynamic QR type
    return QR


class EllCoordSys(NamedTuple):
    """Define coordinate system by identity and inverse/negative. """

    identity: tuple
    negative: tuple
    name: str


WEI_AFF = EllCoordSys((), (1, -1), "Weierstrass Affine")
WEI_HOM_PROJ = EllCoordSys((0, 1, 0), (1, -1, 1), "Weierstrass Homogeneous Projective")
WEI_JAC = EllCoordSys((1, 1, 0), (1, -1, 1), "Weierstrass Jacobian")
ED_AFF = EllCoordSys((0, 1), (-1, 1), "Edwards Affine")
ED_HOM_PROJ = EllCoordSys((0, 1, 1), (-1, 1, 1), "Edwards Homogeneous Projective")
ED_EXT_HOM_PROJ = EllCoordSys(
    (0, 1, 0, 1), (-1, 1, -1, 1), "Edwards Extended Homogeneous Projective"
)
TW_ED_EXT_HOM_PROJ = EllCoordSys((0, 1, 0, 1), (-1, 1, -1, 1), "Tw. Edwards Hom. Proj.")


def affine_tuple_to_coord(CurveElt_subtype, pt_tuple):
    """Convert tuple in affine notation to curve element.

    Invariant: pt_tuple should be defined as (x, y) tuple of field elements
    in affine notation.
    """
    assert len(pt_tuple) == 2
    gf = CurveElt_subtype.field
    target_coord = CurveElt_subtype.coord

    if target_coord == WEI_AFF or target_coord == ED_AFF:
        return CurveElt_subtype(pt_tuple)

    if target_coord in [WEI_HOM_PROJ, ED_HOM_PROJ, WEI_JAC, WEI_HOM_PROJ]:
        x, y = pt_tuple
        new_tuple = (x, y, gf(1))
        return CurveElt_subtype(new_tuple)

    if target_coord == ED_EXT_HOM_PROJ or target_coord == TW_ED_EXT_HOM_PROJ:
        x, y = pt_tuple
        new_tuple = (x, y, x * y, gf(1))
        return CurveElt_subtype(
            new_tuple
        )  # TODO: this said pt_tuple but that seems incorrect (corrected)

    raise NotImplementedError


class EllipticCurveElement(FiniteGroupElement):
    """Common base class for elliptic curve group elements.

    Note: Attribute access of x, y, z (and t) coordinates is defined in
    arithmetic.CurveArithmetic class.
    """

#    __slots__ = "value"

    is_additive = True
    is_multiplicative = False
    is_abelian = True

    def __init__(self, value):
        if isinstance(value, list):
            value = tuple(value)
        self.value = value

    def __repr__(self):
        return f"{self.value}"


@functools.lru_cache(maxsize=None)
def EllipticCurve(params, coord, arithm):
    """Create elliptic curve type for given curve parameters."""

    name = f'Curve({params.name})_{arithm.__name__}'
    EC = type(name, (arithm, EllipticCurveElement), {})
#    EC = type(name, (arithm, EllipticCurveElement), {'__slots__': ()}
    EC.curve = params
    EC.order = params.order
    EC.field = params.field
    EC.is_cyclic = params.is_cyclic
    EC.field.is_signed = False  # for consistency between sectypes and regular types
    EC.coord = coord
    EC.arithm = arithm

    # Add identity and generator to class. Ensure rest of type gets defined above this line.
    EC.identity = EC(coord.identity)
    EC.base_pt = affine_tuple_to_coord(EC, params.base_pt_tuple)
    EC.generator = EC.base_pt  # alias
    globals()[name] = EC  # NB: exploit (almost?) unique name dynamic EC type
    return EC


class FormClassGroupElement(FiniteGroupElement):
    """Class groups represented by binary quadratic form tuple (a, b, c).

    Requirements:
        Form is integral, i.e. a, b, c are integers
        Form is primitive, i.e. gcd(a, b, c) == 1
        Form is positive definite, i.e. Discriminant(f) < 0 and a > 0.

    """

#    __slots__ = "value"

    is_multiplicative = True
    is_abelian = True

    discriminant = None
    bit_length = None
    auto_reduce = True
    order = None 
    is_cyclic = None
    nucomp = True

    def __init__(self, value):
        if isinstance(value[0], PrimeFieldElement):
            value = [int(i) for i in value]
        if isinstance(value, list):
            value = tuple(value)
        # If c coefficient not given, calculate it.
        if len(value) == 2:
            value = list(value)
            assert (value[1]**2 - self.discriminant) % (4*value[0]) == 0
            value += [(value[1]**2 - self.discriminant)//(4*value[0])]
            value = tuple(value)
        self.value = value

    def __getitem__(self, key):
        return self.value[key]

    def __setitem__(self, key, value):
        self.value[key] = value

    def operation(a, b):
        if a.nucomp:
            c = classgrps.nucomp(a, b)
        else:
            # c = classgrps.compose(a, b)
            c = classgrps.shanks_compose(a, b)
        if a.auto_reduce:
            c = classgrps.reduce_form(c)
        return c

    def operation2(a):
        if a.nucomp:
            c = classgrps.nudupl(a)
        else:
            c = classgrps.square(a)
        if a.auto_reduce:
            c = classgrps.reduce_form(c)
        return c

    def inverse(a):
        c = classgrps.forminverse(a)
        if a.auto_reduce:
            c = classgrps.reduce_form(c)
        return c

    def reduce_form(a):
        return classgrps.reduce_form(a)

    def equality(a, b):
        a_red = a.reduce_form()
        b_red = b.reduce_form()
        return a_red.value == b_red.value

    def __repr__(self):
        return str(self.value)


def FormClassGroup(disc=None, bit_length=None, gen_guess=False, gen_set=False, auto_reduce=True, nucomp=True):
    """Create type for form class group for given discriminant (disc) or bit_length.

    Parameters
        disc        initialize group with given discriminant.
        bit_length  initialize group by generating discriminant of given bit_length.
        gen_guess   guess that group is cyclic (high probability), initialize attribute .generator.
        gen_set     create generating_set using generating_system, time: O(|group.discriminant|^(1/2+o(1)))
        auto_reduce enforce reduction after composition operation.
        nucomp      use NUCOMP and NUDUPL instead of compose 

    Requirements
        disc < 0 (we focus on imaginary quadratic fields) and
            disc = 1 mod 4 and Delta is square-free, or,
            disc = 0 mod 4, Delta/4 = 2 or 3 mod 4 and Delta/4 is square-free.

    """
    if disc:
        bit_length = disc.bit_length() + 1  # +1 To include sign bit.
    elif bit_length:
        disc = classgrps.find_fundamental_discriminant(bit_length)
    else:
        raise NotImplementedError
    name = f'FormClassGroup({bit_length})'
    FCG = type(name, (FormClassGroupElement,), {})
#    FCG = type(name, (FormClassGroupElement,), {'__slots__': ()}
    FCG.discriminant = disc
    FCG.bit_length = bit_length
    FCG.auto_reduce = auto_reduce
    FCG.nucomp = nucomp

    # See Cohen-Lenstra conjecture on probability that group is cyclic. 
    # See Cohen "A course in computational number theory", 
    # Conjecture 5.10.1: For any odd prime p, .. probability that Cl_o(D) 
    # (subgroup of elements of odd order) is cyclic .. approx. 0.977575
    # Buell 'Binary Quadratic Forms', Chapt. 8 p. 144: 
    # "class groups tend to be cyclic about 96% of the time [..]"
    if gen_guess:
        FCG.generator = FCG(classgrps.create_generator_of_subgroup(disc))
    if gen_set:
        group_elts, gen_set = classgrps.generating_system(FCG)
        FCG.order = len(group_elts)
        FCG.generating_set = gen_set
        FCG.is_cyclic = len(gen_set) == 1
    FCG.identity = FCG(classgrps.principal_form(disc))
    globals()[name] = FCG  # NB: exploit (almost?) unique name dynamic FCG type
    return FCG


class RubiksCubeGroupElement(SymmetricGroupElement):
    """Class for Rubik's Cube element.
    """
    is_abelian = False
    is_cyclic = False
    indices = rubiks.START_CUBE
    # Order from https://www.gap-system.org/Doc/Examples/rubik.html
    order = 43252003274489856000

    def indices_after_operation(self):
        return [self.value[c] for c in self.indices]

    def __repr__(self):
        new_indices = self.indices_after_operation()
        return rubiks.TEMPLATE.format(*new_indices)


def RubiksCube():
    """Create RubiksCube group element.

    Subgroup of S(n)
    Representation based on image maps, not cycle notation.
    See: https://en.wikipedia.org/wiki/Rubik%27s_Cube_group
    """
    name = f'RubiksCube'
    Rubiks = type(name, (RubiksCubeGroupElement,), {})
    # Super group is S(48) because we have 48 cube facets.
    n = 48
    Rubiks.degree = n
    Rubiks._field = GF(int(gmpy2.next_prime(n*2)))
    Rubiks.identity = Rubiks(tuple(range(n)))
    Rubiks.generating_set = rubiks.clockwise_rotations(Rubiks)

    # Alias to clockwise rotations of faces of the cube.
    Rubiks.top = Rubiks.generating_set[0]
    Rubiks.left = Rubiks.generating_set[1]
    Rubiks.front = Rubiks.generating_set[2]
    Rubiks.right = Rubiks.generating_set[3]
    Rubiks.rear = Rubiks.generating_set[4]
    Rubiks.bottom = Rubiks.generating_set[5]

    globals()[name] = Rubiks  # NB: exploit (almost?) unique name dynamic Rubiks type
    return Rubiks


