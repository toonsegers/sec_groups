import random
import sys, os

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import mpyc.gmpy as gmpy2
from mpyc.runtime import mpc
from sec_groups.fingroups import (
    S,
    QuadraticResidue,
    EllipticCurve,
    FormClassGroup,
    QuadraticResidueElement,
    EllipticCurveElement,
    FormClassGroupElement,
    RubiksCube,
)
import sec_groups.ellcurves as ell
from sec_groups.secgroups import SecureGroup, secure_repeat
from sec_groups.tools import sampling, group_encode
from sec_groups.tools.find_primes import find_safe_primes


async def suite1():
    await mpc.start()

    print("Setting up groups...")
    groups = {
        "QR64": QuadraticResidue(modulus=find_safe_primes(64)[1]),
        "S5": S(5),
        "S11": S(11),
        "ED25519_AFF": EllipticCurve(ell.ED25519, ell.ED_AFF, ell.Edwards_Affine_Arithm),
        "ED25519_PROJ": EllipticCurve(ell.ED25519, ell.ED_HOM_PROJ, ell.Edwards_HomProj_Arithm),
        "BN256_PROJ": EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm),
        "BN256_TWIST": EllipticCurve(ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm),
        "ED25519_PROJ": EllipticCurve(ell.ED25519, ell.ED_HOM_PROJ, ell.Edwards_HomProj_Arithm),
        "BN256_AFF": EllipticCurve(ell.BN256, ell.WEI_AFF, ell.Weierstr_Affine_Arithm),  # No secure arithmetic.
        "Rubiks": RubiksCube(),
        "FCG6": FormClassGroup(disc=-23, gen_guess=True, gen_set=True),
        "FCG32": FormClassGroup(bit_length=32, gen_guess=True, gen_set=True),
        "FCG64": FormClassGroup(bit_length=64, gen_guess=True),
    }  
    print("\nDone.")

    for key, group in groups.items():
        print(f"# Group = {group}")
        if group.generator:
            g = group.generator
        elif group.generating_set:
            g = group.generating_set[0]
        else:
            raise NotImplementedError
        print("Group abelian: ", group.is_abelian)
        print("Group cyclic: ", group.is_cyclic)

        print("Demonstrate inverse, group law, exponentiation.")
        g_inv = g.inverse()
        assert g @ g_inv == group.identity
        g4 = g @ g @ g @ g
        assert g4 == g ^ 4
        x = random.randint(1, 10)
        assert (g ^ x) @ (g.inverse() ^ x) == group.identity
        assert g ^ -1 == ~g

        print("Demonstrate secure arithmetic.")
        # Exclude BN256_AFF because it does not have secure arithmetic.
        if not key in ["BN256_AFF"]:
            sec_group = SecureGroup(group)
            sec_g = sec_group(g)
            sec_id = sec_group.identity
            sec_g4 = sec_id @ sec_g ^ 4
            assert g4 == await mpc.output(sec_g4)

            print("Demo encoding for QR, elliptic curve and form class groups.")
            if issubclass(group, (QuadraticResidueElement, EllipticCurveElement, FormClassGroupElement)):
                # Exclude BN256_TWIST because calculating root of extension field element is inefficient.
                # Exclude FCG6 because encoding codomain too small. 
                if not key in ["BN256_TWIST", "FCG6", "FCG32"]:
                    print("Encode x=", x)
                    encode_zero = True
                    k = 16  # encoding parameter 
                    e = group_encode.encode_v3(x, group, k=k, encode_zero=encode_zero)
                    print(f"Encoding({x=}): {e=}")
                    d = group_encode.decode_v3(e, k=k)
                    print(f"Decoding({e=}): {d=}")
                    assert d == x
                    if encode_zero:
                        d_sec = group_encode.decode_v3((sec_group(e[0]), sec_group(e[1])), k=k)
                    else:
                        d_sec = group_encode.decode_v3(sec_group(e), k=k)
                    d = await mpc.output(d_sec)
                    print("After secure decode: d=", d)
                    assert d == x

            #TODO: add exponentiation for groups of unknown order based on BIB

            print("Demonstrate secure exponentiation if group is of prime order.")
            # TODO: expand to non-prime order groups
            if group.order and gmpy2.is_prime(group.order):
                print("Group is of prime order =", group.order)
                secfld = mpc.SecFld(modulus=sec_group.group.order)
                sec_g8 = secure_repeat(g4, secfld(2), sec_group)
                new_g8 = await mpc.output(sec_g8)
                assert new_g8 == g4 @ g4
                new_g8 = await secure_repeat(g4, secfld(2), sec_group, public_out=True)
                assert new_g8 == g4 @ g4

                # TODO: fix BN256_TWIST; secfld does not work because base fld is an extension field?
                # TODO: FCG6 does not work; why? 
                if not key in ["BN256_TWIST", "FCG6"]:
                    sec_g8 = secure_repeat(sec_g4, mpc.to_bits(secfld(2), l=3), sec_group, exp_as_bits=True)
                    new_g8 = await mpc.output(sec_g8)
                    assert new_g8 == g4 @ g4

        print("Demonstrate random sampling in the group if generating set is available.")
        if group.generating_set != None:
            print(group.generating_set)
            k = 25
            dix = sampling.generate_dixon_cube_lemma13b(group.generating_set, k)
            r = sampling.random_group_element(dix)
            print("Random element: ", r)
        print("===")

    await mpc.shutdown()
    return True


if __name__ == "__main__":
    mpc.run(suite1())


# TODO:
# 1. Add random sampling of element; group.random()
# 1. Debug secure_repeat with bits
# 3. Random exponents in MPC
# # r1 = await mpc.transfer(random.randrange(0, 100), senders=0)
