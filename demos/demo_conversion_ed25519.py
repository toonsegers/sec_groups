"""Demonstrate protocol to distribute shares with one public key.
"""
import os
import sys

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from mpyc.runtime import mpc
from tools.elgamal import encrypt, threshold_decrypt
from tools.group_encode import encode_on_curve, decode_from_curve, encode_v3, decode_v3
from sec_groups.fingroups import EllipticCurve
from sec_groups.secgroups import SecureGroup, repeat_public_base_secret_output
import sec_groups.ellcurves as ell

async def main():
    # group = EllipticCurve(ell.ED25519, ell.ED_AFF, ell.Edwards_Affine_Arithm)
    group = EllipticCurve(ell.ED448, ell.ED_AFF, ell.Edwards_Affine_Arithm)
    sec_grp = SecureGroup(group, group.arithm)  # Passing of (secure) arithmetic is mandatory with current API (also required to receive .x, .y, .z getters/setters).
    g = group.generator

    await mpc.start()
    # Create ElGamal key pair
    sec_a = mpc._random(mpc.SecFld(modulus=group.order))
    # sec_aG = g ^ sec_a
    sec_aG = repeat_public_base_secret_output(g, sec_a, sec_grp)
    pubkey = await mpc.output(sec_aG)

    # Encode input m to curve point; Also encode 0 (to Z), required for decoding in MPC
    m = 100
    # M_encoded, Z_encoded = encode_on_curve(m, group)
    M_encoded = encode_v3(m, group)

    # Encrypt m using ElGamal under pubkey (k is dummy random)
    k = 3809836274351126414438016410872729431457283748709071912489406741072364988984
    elgamal_M = encrypt(M_encoded, k, g, pubkey)
    # elgamal_Z = encrypt(Z_encoded, k, g, pubkey)

    # Decrypt with MPC to shares of m (sec_m)
    sec_M = await threshold_decrypt(elgamal_M, sec_a)
    # sec_Z = await threshold_decrypt(elgamal_Z, sec_a)
    # sec_m = decode_from_curve(sec_M, sec_Z)
    sec_m = decode_v3(sec_M)

    # Show that decrypted ciphertext == plaintext (i.e. [m] recombines to m)
    plain_m = await mpc.output(sec_m)
    assert plain_m == m

    share_m = await mpc.gather(sec_m)
    print(f"MPC party {mpc.pid} holds share of m: {share_m}")

    await mpc.shutdown()


if __name__ == "__main__":
    mpc.run(main())


# TODOs
# 1: Test mpc.gather with sec_grps2
# 2: Make threshold_decrypt an mpc_coro