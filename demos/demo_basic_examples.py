from math import log
import random
import os
import sys

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from mpyc.runtime import mpc

from tools.find_primes import find_safe_primes
from sec_groups.fingroups import S, QuadraticResidue, EllipticCurve, FormClassGroup
import sec_groups.ellcurves as ell
import sec_groups.pairing as pairing
import sec_groups.secgroups as sfg
from sec_groups.classgrps import secure_binary_reduce


def generate_dixon_cube_w(generating_set, k):
    """Generate random group element with k iterations.

    See: https://people.math.carleton.ca/~jdixon/RandomGroupElts.pdf
    """
    from functools import reduce

    s = type(generating_set[0])
    d = len(generating_set)
    # dixon_w = copy.deepcopy(generating_set)
    dixon_w = generating_set.copy()
    for i in range(d + 1, k + 1):
        x_rhs = reduce(
            s.operation, [x for x in dixon_w if bool(random.getrandbits(1))], s.identity
        )
        x_lhs = reduce(
            s.operation,
            [x.inverse() for x in dixon_w if bool(random.getrandbits(1))][::-1],
            s.identity,
        )
        x_next = x_lhs @ x_rhs
        dixon_w.append(x_next)
    return dixon_w


def random_group_element(dixon_cube):
    dix = generate_dixon_cube_w(dixon_cube, len(dixon_cube) + 1)
    return dix[-1]


async def suite1():
    await mpc.start()

    print("Test very small class group (discriminant -23)")
    grp = FormClassGroup(disc=-23, gen_guess=True)
    g = grp.generator
    assert g == grp((2, 1, 3))
    assert g * g == grp((2, -1, 3))
    assert g @ g == grp((2, -1, 3))
    assert g @ g @ g == grp.identity
    grp.order = 3
    assert g ^ grp.order == grp.identity

    print("Test secure class group")
    sec_grp = sfg.SecureGroup(grp)
    secfld = sec_grp.sectype_of_value
    out = sfg.repeat_public_base_secret_output(grp.generator, secfld(3), sec_grp)
    out = await mpc.output(out)
    assert out == grp.identity
    g = sec_grp(grp.generator)
    out = await mpc.output(g * g)
    assert out == grp((2, -1, 3))

    print("Test (secure) class group with larger bit-length.")
    grp = FormClassGroup(bit_length=32, gen_guess=True)
    g = grp.generator
    assert g * g.inverse() == grp.identity
    assert (g ^ 100) * (g.inverse() ^ 100) == grp.identity
    g2 = g * g

    sec_grp = sfg.SecureGroup(grp)
    g = sec_grp(grp.generator)
    out = await mpc.output(g * g)
    assert g2 == out

    print("Test secure_binary_reduce")
    out = await mpc.output(secure_binary_reduce(g * g))
    assert g2 == out

    f = g.inverse()
    out = await mpc.output(f * g)
    assert out == grp.identity

    print("Test class group with very large discriminant (not in MPC).")
    grp = FormClassGroup(bit_length=256, gen_guess=True)  # works for 512 and 4096
    g = grp.generator
    assert g * g.inverse() == grp.identity
    assert (g ^ 100) * (g.inverse() ^ 100) == grp.identity

    print("Test finite_grps in non-MPC setting.")
    print("Test creation of Edwards base point with affine coordinates.")
    curve = EllipticCurve(ell.ED25519, ell.ED_AFF, ell.Edwards_Affine_Arithm)
    base_pt = curve.base_pt
    assert base_pt.on_curve()
    assert base_pt.to_homproj() == base_pt

    print("Test arithmetic for Edwards curve with homogeneous projective coordinates.")
    curve = EllipticCurve(ell.ED25519, ell.ED_HOM_PROJ, ell.Edwards_HomProj_Arithm)
    base_pt = curve.base_pt
    assert base_pt * -1 == -base_pt
    bp4 = base_pt + base_pt + base_pt + base_pt
    assert bp4 == (base_pt * 4)
    assert bp4.on_curve()

    print("Test conversion to extended for Edwards")
    e_bpt = base_pt.to_extended()
    e_bpt4 = e_bpt + e_bpt + e_bpt + e_bpt
    assert e_bpt4 == bp4
    assert e_bpt4.to_affine() == bp4.to_affine()

    print("Test multiplicative arithmetic for same Edwards curve.")
    curve.is_additive = False
    curve.is_multiplicative = True
    bp4 = base_pt * base_pt * base_pt * base_pt
    assert bp4 == (base_pt ** 4)
    # restore values two flags for later use
    curve.is_additive = True
    curve.is_multiplicative = False

    print("Test arithmetic for BN256 curve")
    curve = EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
    bp = curve.base_pt
    bp4 = bp + bp + bp + bp
    assert bp4 == bp * 4
    assert (bp * 4).on_curve()

    print("Test arithmetic for BN256_TWIST curve")
    curve = EllipticCurve(
        ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm
    )
    bp = curve.base_pt
    bp4 = bp + bp + bp + bp
    assert bp4 == bp * 4
    assert (bp * 4).on_curve()

    print("Test multiplicative arithmetic for BN256_TWIST curve")
    curve = EllipticCurve(
        ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm
    )
    curve.is_additive = False
    curve.is_multiplicative = True
    bp = curve.base_pt
    bp4 = bp * bp * bp * bp
    assert bp4 == bp ** 4
    assert (bp ** 4).on_curve()
    # restore values two flags for later use
    curve.is_additive = True
    curve.is_multiplicative = False

    print("Test optimal ate pairing for BN256 curve.")
    # import pairing

    bn_curve = EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
    bn_base_pt = bn_curve.base_pt
    bn_twist = EllipticCurve(
        ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm
    )
    bn_twist_base_pt = bn_twist.base_pt

    a = bn_twist_base_pt
    b = bn_base_pt
    out = pairing.optimal_ate(a, b)
    assert (
        str(out)
        == "((21196725492326442667627602847448291017588785091108454726553939930507871042194x+43377286344385619261366633786771510260308325523330699049673437042931316243080,6519643025219345838017222895845262421933106748812739147102715685220942156402x+49167428805016148020211692433963935169227827998815966523759923767649702416286,20859153989312102139134197789193625833054352751342763983124765227328906753159x+54487609867103086000679472636440811221594782878893833473117980998935745956453),(39332187218762343173097683922586548248512461497033773500717078587710862589062x+53094021597187285527531149248961924798783924165048746910730430368152798446315,30733062774817315099333283633560206070773769397463591953591634872711994123707x+13560812839206871407210482171294630929511117297628119163695762035076749250365,57080396271372976186541856910255122494067079575964633601781325258931774013251x+60034081832369659851990276215766505463071420460830584339216728009871686767851))"
    )

    print("Tests with MPyC")

    print("Test secure Edwards.")
    print("Test Edwards curve with hom. proj. coordinates (non-MPC)")
    ed_curve = EllipticCurve(ell.ED25519, ell.ED_HOM_PROJ, ell.Edwards_HomProj_Arithm)
    base_pt = ed_curve.base_pt
    assert base_pt * -1 == -base_pt
    bp4 = base_pt + base_pt + base_pt + base_pt
    assert bp4 == (base_pt * 4)
    assert bp4.on_curve()

    print("Test Edwards curve with hom. proj. coordinates (secure group)")
    sec_grp = sfg.SecureGroup(ed_curve, ell.Edwards_HomProj_Arithm)
    b = sec_grp(ed_curve.base_pt.value)
    c = sec_grp(ed_curve.base_pt)
    b, c = await mpc.output([b, c])
    assert b == c
    secfld = mpc.SecFld(modulus=sec_grp.group.order)
    out = sfg.repeat_public_base_secret_output(ed_curve.base_pt, secfld(2), sec_grp)
    out = await mpc.output(out)
    assert out == ed_curve.base_pt ^ 2
    out = await sfg.secure_repeat_public_base_public_output(ed_curve.base_pt, secfld(2))
    assert out == ed_curve.base_pt ^ 2
    print("Test Shamir in exponent with secret vector")
    a = [ed_curve.base_pt, ed_curve.base_pt]
    x = [secfld(2), secfld(2)]
    out = await sfg.secure_repeat_public_base_public_output(a, x)
    assert out == ed_curve.base_pt ^ 4

    print("Test secure quadratic residue group.")
    order, modulus = find_safe_primes(64)
    group = QuadraticResidue(modulus=modulus)
    sec_grp = sfg.SecureGroup(group)
    sec_g4 = sec_grp(2) * sec_grp(2)
    g4 = await mpc.output(sec_g4)
    assert int(group(2)) == int(
        group.identity * group(2)
    )  # test identity and __int__ method
    assert g4 == group(2) * group(2)
    secfld = mpc.SecFld(modulus=sec_grp.group.order)
    sec_g4 = sfg.secure_repeat(group(2), secfld(2), sec_grp)
    new_g4 = await mpc.output(sec_g4)
    assert new_g4 == g4

    print("Test secure Symmetric group")
    n = 5
    group = S(n)
    a = group([3, 4, 2, 1, 0])
    b = a @ a
    sec_group = sfg.SecureGroup(group)
    c = sec_group(a)
    d = c @ c
    d_out = await mpc.output(d)
    assert b == d_out
    e = ~c
    f = e @ d
    assert await mpc.output(f == c)
    assert a == await mpc.output(f)

    print("Test secure repeat with S(11) and an element of prime order 11.")
    n = 11
    group = S(n)
    sec_group = sfg.SecureGroup(group)
    a = group([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0])  # ord(a) = 11
    b = a ^ 7
    secfld = mpc.SecFld(modulus=11)
    out = sfg.repeat_public_base_secret_output(a, secfld(7), sec_group)
    out = await mpc.output(out)
    assert out == b

    # WIP
    print("Test generating set")
    # Sn is generated by cycles (0, 1) and (0, .., n-1)
    g1 = group([1, 0])  # cycle (0, 1) in impage map notation.
    g2 = group(list(range(1, n)) + [0])  # cycle (0, .., n-1) in image map notation.
    assert (g2 ^ 2) @ g1 @ (g2 ^ 3) @ g1 @ (g2 ^ 2) == group(
        (7, 8, 9, 10, 0, 1, 3, 2, 4, 6, 5)
    )
    k = int(log(group.order, 2)) ** 2
    dix = generate_dixon_cube_w([g1, g2], k)
    r = random_group_element(dix)
    print("WIP: random based on Dixon cube=", r)

    print("Test more of S(11).")
    group = S(
        n
    )  # TODO: check required size of base field (e.g. for seclists/to_bits to work)
    sec_group = sfg.SecureGroup(group)
    g = group([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0])  # ord(a) = 11
    g4 = g ^ 4
    sec_g = sec_group(g)
    sec_id = sec_group.identity
    id_out = await mpc.output(sec_id)
    sec_g4 = sec_g ^ 4
    sec_g4 = sec_g4 @ sec_id  # include identity in test
    new_g4 = await mpc.output(sec_g4)
    assert new_g4 == g4
    secfld = mpc.SecFld(
        modulus=11
    )  # group order is not prime, but element order is; hence, don't do: secfld = mpc.SecFld(modulus=sec_group.group.order)
    sec_g8 = sfg.repeat_public_base_secret_output(new_g4, secfld(2), sec_group)
    new_g8 = await mpc.output(sec_g8)
    assert new_g8 == g4 @ g4
    sec_g8 = sfg.secure_repeat(new_g4, secfld(2), sec_group)
    new_g8 = await mpc.output(sec_g8)
    assert new_g8 == g4 @ g4
    sec_g8 = sfg.secure_repeat(new_g4, secfld(2))
    new_g8 = await mpc.output(sec_g8)
    assert new_g8 == g4 @ g4

    print("Test secure Ed25519 and repeat of secure point with public exponent.")
    curve = EllipticCurve(ell.ED25519, ell.ED_HOM_PROJ, ell.Edwards_HomProj_Arithm)
    sec_grp = sfg.SecureGroup(
        curve, ell.Edwards_HomProj_Arithm
    )  # secure curve api requires explicit passing of (secure) arithmetic as argument
    bp4 = curve.base_pt * 4
    sec_bp = sec_grp(curve.base_pt)
    sec_id = sec_grp.identity
    id_out = await mpc.output(sec_id)
    sec_bp4 = sec_bp * 4
    sec_bp4 = sec_bp4 + sec_id  # include identity in test
    new_bp4 = await mpc.output(sec_bp4)
    assert new_bp4 == bp4
    secfld = mpc.SecFld(modulus=sec_grp.group.order)
    sec_bp8 = sfg.repeat_public_base_secret_output(new_bp4, secfld(2), sec_grp)
    new_bp8 = await mpc.output(sec_bp8)
    assert new_bp8 == bp4 + bp4
    # WIP
    print("Test secure_repeat(.. public_out=True)")
    new_bp8 = await sfg.secure_repeat(new_bp4, secfld(2), sec_grp, public_out=True)
    assert new_bp8 == bp4 + bp4
    # TODO: understand who secure_repeat needs to be awaited
    # sec_bp8 = sfg.secure_repeat(new_bp4, secfld(2), sec_grp, public_out = True)
    # new_bp8 = await mpc.output(sec_bp8)
    # assert new_bp8 == bp4 + bp4
    print(
        "Test secure_repeat_secret_base_secret_output (slow if exponent not short list of bits)"
    )
    sec_bp8 = sfg.secure_repeat(
        sec_bp4, mpc.to_bits(secfld(2), l=3), sec_grp, exp_as_bits=True
    )
    new_bp8 = await mpc.output(sec_bp8)
    assert new_bp8 == new_bp4 + new_bp4

    sec_bp8 = sfg.secure_repeat(new_bp4, secfld(2), sec_grp)
    new_bp8 = await mpc.output(sec_bp8)
    assert new_bp8 == bp4 + bp4
    sec_bp8 = sfg.secure_repeat(new_bp4, secfld(2))
    new_bp8 = await mpc.output(sec_bp8)
    assert new_bp8 == bp4 + bp4

    print("Test creation of secure group for BN256 curve.")
    curve = EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
    sec_grp = sfg.SecureGroup(curve, ell.Weierstr_HomProj_Arithm)
    bp4 = curve.base_pt * 4
    sec_bp = sec_grp(curve.base_pt)
    sec_id = sec_grp.identity
    id_out = await mpc.output(sec_id)
    sec_bp4 = sec_bp * 4
    sec_bp4 = sec_bp4 + sec_id  # include identity in test
    new_bp4 = await mpc.output(sec_bp4)
    assert new_bp4 == bp4
    secfld = mpc.SecFld(modulus=sec_grp.group.order)
    sec_bp8 = sfg.repeat_public_base_secret_output(new_bp4, secfld(2), sec_grp)
    new_bp8 = await mpc.output(sec_bp8)
    assert new_bp8 == bp4 + bp4
    sec_bp8 = sfg.secure_repeat(new_bp4, secfld(2), sec_grp)
    new_bp8 = await mpc.output(sec_bp8)
    assert new_bp8 == bp4 + bp4
    sec_bp8 = sfg.secure_repeat(new_bp4, secfld(2))
    new_bp8 = await mpc.output(sec_bp8)
    assert new_bp8 == bp4 + bp4

    print("Test creation of secure group for BN256_TWIST curve.")
    curve = EllipticCurve(
        ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm
    )
    # print(curve.curve.name == "BN256_twist")
    sec_grp = sfg.SecureGroup(curve, ell.Weierstr_HomProj_Arithm)
    bp4 = curve.base_pt * 4
    sec_bp = sec_grp(curve.base_pt)
    sec_id = sec_grp.identity
    id_out = await mpc.output(sec_id)
    sec_bp4 = sec_bp * 4
    sec_bp4 = sec_bp4 + sec_id  # include identity in test
    new_bp4 = await mpc.output(sec_bp4)
    assert new_bp4 == bp4
    secfld = mpc.SecFld(modulus=sec_grp.group.order)
    sec_bp8 = sfg.repeat_public_base_secret_output(new_bp4, secfld(2), sec_grp)
    new_bp8 = await mpc.output(sec_bp8)
    assert new_bp8 == bp4 + bp4
    sec_bp8 = sfg.secure_repeat(new_bp4, secfld(2), sec_grp)
    new_bp8 = await mpc.output(sec_bp8)
    assert new_bp8 == bp4 + bp4
    sec_bp8 = sfg.secure_repeat(new_bp4, secfld(2))
    new_bp8 = await mpc.output(sec_bp8)
    assert new_bp8 == bp4 + bp4

    print("Test affine arithmetic for short Weierstrass form (BN256)")
    curve = EllipticCurve(ell.BN256, ell.WEI_AFF, ell.Weierstr_Affine_Arithm)
    bp2 = curve.base_pt + curve.base_pt
    bp2 = bp2 + curve.identity
    bp3 = bp2 + curve.base_pt
    bp4 = curve.base_pt * 4
    assert bp4 == bp3 + curve.base_pt

    # TODO: Apply Jacobian arithmetic
    print("Test Jacobian arithmetic for short Weierstrass form (BN256)")
    curve = EllipticCurve(ell.BN256, ell.WEI_AFF, ell.Weierstr_Affine_Arithm)
    bp2 = curve.base_pt + curve.base_pt
    bp2 = bp2 + curve.identity
    bp3 = bp2 + curve.base_pt
    bp4 = curve.base_pt * 4
    assert bp4 == bp3 + curve.base_pt

    print("Test Jacobian arithmetic for short Weierstrass form (BN256_TWIST)")
    curve = EllipticCurve(ell.BN256_TWIST, ell.WEI_AFF, ell.Weierstr_Affine_Arithm)
    bp2 = curve.base_pt + curve.base_pt
    bp2 = bp2 + curve.identity
    bp3 = bp2 + curve.base_pt
    bp4 = curve.base_pt * 4
    assert bp4 == bp3 + curve.base_pt

    print("Test arithmetic conversion for short Weierstrass form (BN256_TWIST)")
    curve = EllipticCurve(
        ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm
    )
    bp = curve.base_pt
    bp = curve.base_pt.to_jacobian()
    bp2 = bp + bp
    bp2 = bp2 + type(bp).identity
    bp3 = bp2 + bp
    bp4 = bp * 4
    assert bp4 == bp3 + bp

    # # Work in progress: requires secure exponentiation with secret output
    # print("Test optimal ate pairing for secure BN256 curve elements, requires secure exponentiation with secret output.")
    # bn = EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
    # bn_twist = EllipticCurve(ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
    # sec_bn = sfg.SecureGroup(bn, ell.Weierstr_HomProj_Arithm)
    # sec_bn_twist = sfg.SecureGroup(bn_twist, ell.Weierstr_HomProj_Arithm)
    # sec_bn_base_pt = sec_bn(bn.base_pt)
    # sec_bn_twist_base_pt = sec_bn_twist(bn_twist.base_pt)
    # # secfld = mpc.SecFld(modulus=sec_bn.group.order)
    # # sec_bits = mpc.to_bits(secfld(5))
    # # bits = await mpc.output(sec_bits)
    # # print(bits)
    # a = sec_bn_twist_base_pt
    # b = sec_bn_base_pt
    # out = pairing.optimal_ate(a, b)

    await mpc.shutdown()
    return True


if __name__ == "__main__":
    mpc.run(suite1())


# TODO:
# 1: Reduce is_prime checks if redundant (for QR group)
# 2: Consider unifying group.field (curve group) and group._field interface (QR group)
# 1: Allow QR element @ secure QR element? define operation() in SecureQRElement?
# 1: Ensure a.to_affine() and a.to_jacobian() (etc) methods result in inplace replacements of the variable a
# 1: Prover side of AC20 to support Weierstrass homproj coord's in MPC. (Done? check)
# 1: For pairing check: Public AC20 proof (curve point) to be converted to jacobian (or affine?) as input for pairing?
# 3: Verify if affine coordinates are a requirement for this optimal_ate implementation (hyp: projective coords with z=1 work, check if z==1 holds in Miller loop)
# 3: Consider representing affine identity as (0,0) instead of ()
# 3: Als numpy arrays ondersteund worden met Share, dan sec_grps laten inheriten van Share?
# 3: Schnorr toevoegen; QR is een subclass van Schnorr; consider: encoding voor generieke Schnorr vinden en implementeren?
# 4: Support extension fields for quadratic residues (currently, the int-casting: outp = int(mpc.run(func(element))) in _mpc_func breaks compatibility)
