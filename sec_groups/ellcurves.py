from enum import Enum
import functools
from operator import add, mul
from typing import NamedTuple, Any

from mpyc.finfields import GF, PrimeFieldElement
from mpyc.gfpx import GFpX
from sec_groups.fingroups import (
    EllipticCurve,
    WEI_AFF,
    ED_AFF,
    ED_HOM_PROJ,
    ED_EXT_HOM_PROJ,
    WEI_HOM_PROJ,
    WEI_JAC,
)
from sec_groups.secgroups import SecureGroupElement


class CurveParams:
    """Contains curve parameters.

    Invariants:
        self.equation assumes affine coordinates
        self.base_pt_tuple assumes affine coordinates
    """

    def __init__(self, *, name, order, gf, is_cyclic):
        self.name = name
        self.order = order
        self.field = gf
        self.is_cyclic = is_cyclic

    def set_constants(self, *, a=1, b=1, c=1, d=1):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def set_equation(self, eq):
        self.equation = eq

    def set_base_pt(self, value):
        assert isinstance(value, tuple)
        assert len(value) == 2, "Base point should be (x, y) tuple, in affine notation."
        if isinstance(value[0], int):
            base_pt = map(self.field, value)
        elif isinstance(value[0], self.field):
            base_pt = value
        self.base_pt_tuple = base_pt


def set_edwards_eq(*, a=1, c=1, d=1):
    """Return equation that defines Edwards curve.

    To check if point is on curve.
    Edwards curve requires c, d params: x**2 + y**2 == c**2 * (1 + d * x**2 * y**2)
    Twisted curve requires a, d params: a * x**2 + y**2 == 1 + d * x**2 * y**2
    """

    def edwards_eq(x, y):
        return a * x ** 2 + y ** 2 == c ** 2 * (1 + d * x ** 2 * y ** 2)

    return edwards_eq


def set_weierstrass_eq(*, a=1, b=1):
    """Return equation that defines Weierstrass curve. """

    def wei_eq(x, y):
        return y ** 2 == x ** 3 + a * x + b

    return wei_eq


def _ed25519():
    """Edwards curve Ed25519.

    Link: https://en.wikipedia.org/wiki/EdDSA#Ed25519
    """
    q = 2 ** 255 - 19
    order = 2 ** 252 + 27742317777372353535851937790883648493
    gf = GF(q)
    ed = CurveParams(name="ED25519", order=order, gf=gf, is_cyclic = True)
    ed.set_constants(a=gf(-1), d=gf(-121665) / gf(121666))
    ed.set_equation(set_edwards_eq(a=ed.a, c=ed.c, d=ed.d))
    ed.set_base_pt(
        (
            gf(
                15112221349535400772501151409588531511454012693041857206046113283949847762202
            ),
            gf(4) / gf(5),
        )
    )
    return ed


def _ed448():
    """Edwards curve Ed448 'Goldilocks'.

    Link: https://en.wikipedia.org/wiki/Curve448
    """
    q = 2 ** 448 - 2 ** 224 - 1
    order = 2 ** 446 - int(
        "0x8335dc163bb124b65129c96fde933d8d723a70aadc873d6d54a7bb0d", 16
    )  # from: https://eprint.iacr.org/2015/625.pdf
    gf = GF(q)
    ed = CurveParams(name="ED448", order=order, gf=gf, is_cyclic = True)
    ed.set_constants(a=gf(1), c=gf(1), d=gf(-39081))
    ed.set_equation(set_edwards_eq(a=ed.a, c=ed.c, d=ed.d))
    ed.set_base_pt(
        (
            gf(
                117812161263436946737282484343310064665180535357016373416879082147939404277809514858788439644911793978499419995990477371552926308078495
            ),
            gf(19),
        )
    )
    return ed


def _bn256():
    """Define Barreto-Naehrig 256 curve.

    Curve equation: y^2 = x^3 + 3 over F_p
    """
    u = 1868033 ** 3
    p = (((u + 1) * 6 * u + 4) * u + 1) * 6 * u + 1
    order = p - 6 * u ** 2
    gf = GF(p)
    bn = CurveParams(name="BN256", order=order, gf=gf, is_cyclic = True)
    bn.set_constants(a=gf(0), b=gf(3))
    bn.set_equation(set_weierstrass_eq(a=bn.a, b=bn.b))
    bn.set_base_pt((gf(1), gf(-2)))
    return bn


def _bn256_twist():
    """Define sextic twist Barreto-Naehrig 256 curve.

    Curve equation: y^2 = x^3 + b/xi, with xi in F_p^2 following [NNS10],
    Section 3.1. For [NNS10], see: https://eprint.iacr.org/2010/186
    """
    u = 1868033 ** 3
    p = 36 * u ** 4 + 36 * u ** 3 + 24 * u ** 2 + 6 * u + 1
    assert p % 4 == 3
    order = p - 6 * u ** 2
    irred_poly = GFpX(p)([1, 0, 1])  # i^2 == -1 from Section 3.1
    gf = GF(irred_poly)
    b = gf([3, 0, 0])
    i = gf([0, 1, 0])
    xi = i + gf([3, 0, 0])  # xi = i + 3

    bn = CurveParams(name="BN256_twist", order=order, gf=gf, is_cyclic = True)
    bn.set_constants(a=gf([0, 0, 0]), b=b / xi)
    bn.set_equation(set_weierstrass_eq(a=bn.a, b=bn.b))
    base_x = gf(
        [
            64746500191241794695844075326670126197795977525365406531717464316923369116492,
            21167961636542580255011770066570541300993051739349375019639421053990175267184,
            0,
        ]
    )
    base_y = gf(
        [
            17778617556404439934652658462602675281523610326338642107814333856843981424549,
            20666913350058776956210519119118544732556678129809273996262322366050359951122,
            0,
        ]
    )
    bn.set_base_pt((base_x, base_y))

    return bn


ED25519 = _ed25519()
ED448 = _ed448()
BN256 = _bn256()
BN256_TWIST = _bn256_twist()


def negative(pt):
    neg = tuple([z1 * z2 for z1, z2 in zip(pt.coord.negative, pt.value)])
    return type(pt)(neg)


class _Point_xyz(NamedTuple):
    """Lightweight helper container to represent point.

    Note: This is not a curve element with all attributes.
    Attributes:
        x, y, z: Access coordinate of point (x, y, z)
    """

    x: Any
    y: Any
    z: Any


def _add_hwcd08_twisted_edw_extended_coord(pt1, pt2):
    """Addition formula for Twisted Edwards points in extended projective coordinates.

    Requires: Twisted Edwards curve (curve constant a = -1) and extended coordinates.
    See Hisil et al. (2008) [HWCD08], Section 4.2 (p. 11) for 4 processors
    Link: https://eprint.iacr.org/2008/522
    """

    x1, y1, t1, z1 = pt1.x, pt1.y, pt1.t, pt1.z
    x2, y2, t2, z2 = pt2.x, pt2.y, pt2.t, pt2.z

    d = pt1.curve.d

    r1 = y1 - x1
    r2 = y2 - x2
    r3 = y1 + x1
    r4 = y2 + x2

    r5 = r1 * r2
    r6 = r3 * r4
    r7 = t1 * t2
    r8 = z1 * z2

    r7 = (
        2 * d * r7
    )  # k = 2d' or k=2d in ed25519, as per hisil et al., p. 6; see also: http://hyperelliptic.org/efd/g1p/auto-twisted-extended-1.html
    r8 = 2 * r8

    r1 = r6 - r5
    r2 = r8 - r7
    r3 = r8 + r7
    r4 = r6 + r5

    x3 = r1 * r2
    y3 = r3 * r4
    t3 = r1 * r4
    z3 = r2 * r3

    return (x3, y3, t3, z3)


def add_edwards_extended(pt1, pt2):
    pt3 = _add_hwcd08_twisted_edw_extended_coord(pt1, pt2)
    return type(pt1)(pt3)


def add_rcb16_hom_proj(pt1, pt2):
    """Implementation of complete addition formula from [RCB16].

    Formula is complete on every prime order short Weierstrass 
    curve defined over a field k with char(k) != 2, 3.
    Algorithm 7.
    Link: https://www.microsoft.com/en-us/research/wp-content/uploads/2016/06/complete-2.pdf
    """
    b = pt1.curve.b

    b3 = b * 3

    x1 = pt1.x
    y1 = pt1.y
    z1 = pt1.z
    x2 = pt2.x
    y2 = pt2.y
    z2 = pt2.z

    t0 = x1 * x2
    t1 = y1 * y2
    t2 = z1 * z2

    t3 = x1 + y1
    t4 = x2 + y2
    t3 = t3 * t4

    t4 = t0 + t1
    t3 = t3 - t4
    t4 = y1 + z1

    x3 = y2 + z2
    t4 = t4 * x3
    x3 = t1 + t2

    t4 = t4 - x3
    x3 = x1 + z1
    y3 = x2 + z2

    x3 = x3 * y3
    y3 = t0 + t2
    y3 = x3 - y3

    x3 = t0 + t0
    t0 = x3 + t0
    t2 = b3 * t2

    z3 = t1 + t2
    t1 = t1 - t2
    y3 = b3 * y3

    x3 = t4 * y3
    t2 = t3 * t1
    x3 = t2 - x3

    y3 = y3 * t0
    t1 = t1 * z3
    y3 = t1 + y3

    t0 = t0 * t3
    z3 = z3 * t4
    z3 = z3 + t0

    return type(pt1)((x3, y3, z3))


def double_wei_hom_proj(pt1):
    # TODO: Add optimized formula for doubling.
    return add_rcb16_hom_proj(pt1, pt1)


def add_2007_bl(pt1, pt2):
    """Add points, formula add_2007_bl by Bernstein and Lange (2007). 

    Link: http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
    """
    x1, y1, z1 = pt1.x, pt1.y, pt1.z
    x2, y2, z2 = pt2.x, pt2.y, pt2.z

    z1z1 = z1 ** 2
    z2z2 = z2 ** 2
    u1 = x1 * z2z2
    u2 = x2 * z1z1
    s1 = y1 * z2 * z2z2
    s2 = y2 * z1 * z1z1
    h = u2 - u1
    i = (2 * h) ** 2
    j = h * i
    r = 2 * (s2 - s1)
    v = u1 * i
    x3 = r ** 2 - j - 2 * v
    y3 = r * (v - x3) - 2 * s1 * j
    z3 = ((z1 + z2) ** 2 - z1z1 - z2z2) * h

    return _Point_xyz(x3, y3, z3)


def dbl_2009_l(pt1):
    """Double point, formula dbl-2009-l by Lange (2009). 

    Link: http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
    """
    x1, y1, z1 = pt1.x, pt1.y, pt1.z

    a = x1 ** 2
    b = y1 ** 2
    c = b ** 2
    d = 2 * ((x1 + b) ** 2 - a - c)
    e = 3 * a
    f = e ** 2
    x3 = f - 2 * d
    y3 = e * (d - x3) - 8 * c
    z3 = 2 * y1 * z1

    return _Point_xyz(x3, y3, z3)


def add_weierstrass_jacobian(pt1, pt2):
    """Add Weierstrass points with Jacobian coordinates.

    Requires short Weierstrass form and Jacobian coordinates.    
    Does not assume z1 = 1 or z2 = 1.

    Args:
        pt1, pt2 (CurveElement): Points to apply group law to.

    Returns:
        type(pt1)
    """
    assert pt1.coord == WEI_JAC
    assert pt2.coord == WEI_JAC

    if pt1 == type(pt1).identity:
        return pt2
    elif pt2 == type(pt1).identity:
        return pt1
    elif pt1 == pt2:
        pt3 = dbl_2009_l(pt1)
        return type(pt1)((pt3.x, pt3.y, pt3.z))
    else:
        pt3 = add_2007_bl(pt1, pt2)
        return type(pt1)((pt3.x, pt3.y, pt3.z))


def double_weierstrass_jacobian(pt1, pt2):
    """Double Weierstrass point with Jacobian coordinates.

    Requires short Weierstrass form and Jacobian coordinates.    
    Does not assume z1 = 1.

    Args:
        pt1 (CurveElement): Point to double.

    Returns:
        type(pt1)
    """
    pt3 = dbl_2009_l(pt1)
    return type(pt1)((pt3.x, pt3.y, pt3.z))


def add_weierstrass_affine(pt1, pt2):
    """Add Weierstrass points with affine coordinates.

    Requires short Weierstrass form and affine coordinates.    

    Args:
        pt1, pt2 (CurveElement): Points to apply group law to.

    Returns:
        type(pt1)

    Algorithm documented in Silverman, The arithmetic of elliptic curves(1994) 
    and thesis Hisil (2010), for example.
    Link: https://core.ac.uk/download/pdf/10898289.pdf (Algorithm 4.1.1)
    """
    assert pt1.coord == WEI_AFF
    assert pt2.coord == WEI_AFF

    if pt1 == type(pt1).identity:
        return pt2
    elif pt2 == type(pt1).identity:
        return pt1
    elif pt1.x == pt2.x:
        if pt1.y != pt2.y:
            return type(pt1).identity
        elif pt1.y == 0:
            return type(pt1).identity
        else:
            a = pt1.curve.a
            x3 = ((3 * pt1.x ** 2 + a) / (2 * pt1.y)) ** 2 - 2 * pt1.x
            y3 = ((3 * pt1.x ** 2 + a) / (2 * pt1.y)) * (pt1.x - x3) - pt1.y
            return type(pt1)((x3, y3))
    else:
        x3 = ((pt1.y - pt2.y) / (pt1.x - pt2.x)) ** 2 - pt1.x - pt2.x
        y3 = ((pt1.y - pt2.y) / (pt1.x - pt2.x)) * (pt1.x - x3) - pt1.y
        return type(pt1)((x3, y3))


def add_mmadd_2007_bl(pt1, pt2):
    """Apply addition law 'mmadd-2007-bl'.

    Args:
        pt1 (CurveElement): First point.
        pt2 (CurveElement): Second point.

    Returns:
        tuple: Coordinates of third point.

    See: http://hyperelliptic.org/EFD/g1p/data/edwards/projective/addition/mmadd-2007-bl
    """
    x1, y1 = pt1.x, pt1.y
    x2, y2 = pt2.x, pt2.y

    a = pt1.curve.a
    c = pt1.curve.c
    d = pt1.curve.d

    C = x1 * x2
    D = y1 * y2
    E = d * C * D
    x3 = (1 - E) * ((x1 + y1) * (x2 + y2) - C - D)
    y3 = (1 + E) * (D - a * C)
    z3 = c * (1 - E * E)
    return _Point_xyz(x3, y3, z3)


def add_2008_bbjlp(pt1, pt2):
    """Apply addition law "add_2008_bbjlp".

    Requires twisted Edwards curve in hom. proj. coordinates.
    Does not assume z1 = 1 and z2 = 1.
    See: https://www.hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#addition-add-2008-bbjlp
    """
    x1, y1, z1 = pt1.x, pt1.y, pt1.z
    x2, y2, z2 = pt2.x, pt2.y, pt2.z

    a = pt1.curve.a
    d = pt1.curve.d

    A = z1 * z2
    B = A ** 2
    C = x1 * x2
    D = y1 * y2
    E = d * C * D
    F = B - E
    G = B + E
    new_x3 = A * F * ((x1 + y1) * (x2 + y2) - C - D)
    new_y3 = A * G * (D - a * C)
    new_z3 = F * G
    return _Point_xyz(new_x3, new_y3, new_z3)


def set_returntype_for_operation(pt1, pt2):
    # If one of the inputs is a SecureGroupElement, then ensure return type is a sectype..
    if isinstance(pt1, SecureGroupElement):
        newtype = type(pt1)
    else:
        newtype = type(pt2)
    return newtype


def add_edwards_hom_proj(pt1, pt2):
    """Add Edwards points with homogeneous projective coordinates.

    Requires twisted Edwards curve in hom. proj. coordinates.
    Does not assume z1 = 1 and z2 = 1.

    Args:
        pt1, pt2 (CurveElement): Points to apply group law to.

    Returns:
        type(pt1)

    Addition law 'add_2008_bbjlp'.
    See: https://www.hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#addition-add-2008-bbjlp
    """
    assert pt1.coord == ED_HOM_PROJ
    assert pt2.coord == ED_HOM_PROJ
    pt3 = add_2008_bbjlp(pt1, pt2)
    newtype = set_returntype_for_operation(pt1, pt2)
    return newtype((pt3.x, pt3.y, pt3.z))


def add_edwards_hom_proj_z1(pt1, pt2):
    """Add Edwards points with homogeneous projective coordinates.

    Requires z1 == 1 and z2 == 1.

    Args:
        pt1, pt2 (CurveElement): Points to apply group law to.

    Returns:
        type(pt1)

    Addition law 'mmadd-2007-bl'.
    See: http://hyperelliptic.org/EFD/g1p/data/edwards/projective/addition/mmadd-2007-bl
    """
    assert pt1.coord == ED_HOM_PROJ
    assert pt2.coord == ED_HOM_PROJ
    assert pt1.z == 1
    assert pt2.z == 1
    pt3 = add_mmadd_2007_bl(pt1, pt2)
    newtype = set_returntype_for_operation(pt1, pt2)
    return newtype((pt3.x, pt3.y, pt3.z))


def add_edwards_affine(pt1, pt2):
    """ Add Edwards points using 'mmadd-2007-bl', convert to affine. """
    assert pt1.coord == ED_AFF
    assert pt2.coord == ED_AFF
    pt3 = add_mmadd_2007_bl(pt1, pt2)
    x3 = pt3.x / pt3.z
    y3 = pt3.y / pt3.z
    newtype = set_returntype_for_operation(pt1, pt2)
    return newtype((x3, y3))


def ed_affine_to_extended(pt):
    """Map (x, y) to (x : y : x*y :  1)."""
    new_curve = EllipticCurve(pt.curve, ED_EXT_HOM_PROJ, Edwards_ExtProj_Arithm)
    return new_curve((pt.x, pt.y, pt.x * pt.y, new_curve.field(1)))


def extended_to_affine(pt):
    """Convert extended Edwards point to affine point."""

    x, y, t, z, = pt.x, pt.y, pt.t, pt.z
    z_inv = z.reciprocal()
    x = x * z_inv
    y = y * z_inv
    new_curve = EllipticCurve(pt.curve, ED_AFF, Edwards_Affine_Arithm)
    return new_curve((x, y))


def ed_hom_proj_to_extended(pt):
    """Convert hom. proj. coordinates to extended hom. proj. coordinates."""
    a_pt = ed_hom_proj_to_affine(pt)
    e_pt = ed_affine_to_extended(a_pt)
    return e_pt


def ed_affine_to_hom_proj(pt):
    """Convert affine coordinates to extended hom. proj. coordinates."""
    new_curve = EllipticCurve(pt.curve, ED_HOM_PROJ, Edwards_HomProj_Arithm)
    return new_curve((pt.x, pt.y, pt.field(1)))


def ed_hom_proj_to_affine(pt):
    """Map (X:Y:Z) to (X/Z, Y/Z, 1), assumes hom. projective coordinates."""
    ptz = int(pt.z)
    new_curve = EllipticCurve(pt.curve, ED_AFF, Edwards_Affine_Arithm)
    if ptz == 1:
        # pt is affine
        return new_curve((pt.x, pt.y))

    if ptz == 0:
        return new_curve.identity
    else:
        z_inv = pt.z.reciprocal()
        new_x = pt.x * z_inv
        new_y = pt.y * z_inv
        return new_curve((new_x, new_y))


def wei_hom_proj_to_affine(pt):
    """Map (X:Y:Z) to (X/Z, Y/Z, 1), assumes hom. projective coordinates."""
    ptz = int(pt.z)
    new_curve = EllipticCurve(pt.curve, WEI_AFF, Weierstr_Affine_Arithm)
    if ptz == 1:
        # pt is affine
        return new_curve((pt.x, pt.y))

    if ptz == 0:
        # TODO: check: Does this covers all cases for identity element?
        return new_curve.identity
    else:
        z_inv = pt.z.reciprocal()
        new_x = pt.x * z_inv
        new_y = pt.y * z_inv
        return new_curve((new_x, new_y))


def wei_jacobian_to_affine(pt):
    """Map (X:Y:Z) to (X/Z^2, Y/Z^3, 1), assumes Jacobian coordinates."""
    ptz = int(pt.z)
    new_curve = EllipticCurve(pt.curve, WEI_AFF, Weierstr_Affine_Arithm)
    if ptz == 1:
        # pt is affine
        return new_curve((pt.x, pt.y))

    if ptz == 0:
        # TODO: check: Does this covers all cases for identity element?
        return new_curve.identity
    else:
        z_inv = pt.z.reciprocal()
        new_x = pt.x * (z_inv ** 2)
        new_y = pt.y * (z_inv ** 3)
        return new_curve((new_x, new_y))


def wei_affine_to_jacobian(pt):
    """Map (x, y) to (x : y :  1)."""
    new_curve = EllipticCurve(pt.curve, WEI_JAC, Weierstr_Jacobian_Arithm)
    return new_curve((pt.x, pt.y, new_curve.field(1)))


def wei_hom_proj_to_jacobian(pt):
    """Convert hom. proj. coordinates to Jacobian coordinates."""
    a_pt = wei_hom_proj_to_affine(pt)
    j_pt = wei_affine_to_jacobian(a_pt)
    return j_pt


def check_equal_in_affine(pt1, pt2):
    pt1_aff = pt1.to_affine()
    pt2_aff = pt2.to_affine()
    return pt1_aff.value == pt2_aff.value


def on_affine_curve(pt):
    affine_pt = pt.to_affine()
    return pt.curve.equation(affine_pt.x, affine_pt.y)


def on_curve_affine(pt):
    return pt.curve.equation(pt.x, pt.y)


class CurveArithmetic:
    """Abstract base class for curve arithmetic.

    Defaults are defined via mixins. Subtype factory EllipticCurve()
    consumes mixin to add default operators to curve() instance.
    """

    __slots__ = ()

    @property
    def x(self):
        return self.value[0]

    @property
    def y(self):
        return self.value[1]

    @property
    def z(self):
        return self.value[2]

    @x.setter
    def x(self, value):
        self.value[0] = value

    @y.setter
    def y(self, value):
        self.value[1] = value

    @z.setter
    def z(self, value):
        self.value[2] = value

    inverse = negative
    # Set oblivious = True if arithmetic works in MPC setting.
    oblivious = False


class Edwards_Affine_Arithm(CurveArithmetic):
    """Implement Edwards curve arithmetic for affine coordinates."""

    operation = add_edwards_affine
    on_curve = on_curve_affine
    to_affine = _ = lambda x: x  # identity function
    to_extended = ed_affine_to_extended
    to_homproj = ed_affine_to_hom_proj
    equality = check_equal_in_affine
    on_curve = on_affine_curve
    # TODO: verify
    oblivious = True


class Edwards_HomProj_Arithm(CurveArithmetic):
    """Implement Edwards curve arithmetic for affine coordinates."""

    operation = add_edwards_hom_proj
    to_affine = ed_hom_proj_to_affine
    to_extended = ed_hom_proj_to_extended
    equality = check_equal_in_affine
    on_curve = on_affine_curve
    oblivious = True


class Edwards_ExtProj_Arithm(CurveArithmetic):
    """Implement Edwards curve arithmetic for extended coordinates."""

    operation = add_edwards_extended
    to_affine = extended_to_affine
    equality = check_equal_in_affine
    on_curve = on_affine_curve
    oblivious = True

    @property
    def t(self):
        return self.value[2]

    @t.setter
    def t(self, value):
        self.value[2] = value

    @property
    def z(self):
        return self.value[3]

    @z.setter
    def z(self, value):
        self.value[3] = value


class Weierstr_Affine_Arithm(CurveArithmetic):
    """Implement Weierstrass curve arithmetic for affine coordinates."""

    operation = add_weierstrass_affine
    to_affine = _ = lambda x: x  # identity function
    to_jacobian = wei_affine_to_jacobian
    equality = check_equal_in_affine
    on_curve = on_affine_curve
    oblivious = False


class Weierstr_HomProj_Arithm(CurveArithmetic):
    """Implement Weierstrass curve arithmetic for hom. proj. coordinates."""

    operation = add_rcb16_hom_proj
    to_affine = wei_hom_proj_to_affine
    to_jacobian = wei_hom_proj_to_jacobian
    equality = check_equal_in_affine
    on_curve = on_affine_curve
    oblivious = True


class Weierstr_Jacobian_Arithm(CurveArithmetic):
    """Implement Weierstrass curve arithmetic for Jacobian coordinates."""

    operation = add_weierstrass_jacobian
    to_affine = wei_jacobian_to_affine
    equality = check_equal_in_affine
    on_curve = on_affine_curve
    oblivious = False
