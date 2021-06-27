""" Demonstrate protocol to distribute secret shares with one public key

Work-in-progress file.
"""
import os
import sys

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from mpyc.runtime import mpc
from sec_groups.tools.group_encode import encode_to_schnorr, decode_from_schnorr, encode_v3, decode_v3
from sec_groups.tools.elgamal import encrypt, threshold_decrypt
from sec_groups.tools.find_primes import find_safe_primes
from sec_groups.fingroups import QuadraticResidue
from sec_groups.secgroups import SecureGroup, repeat_public_base_secret_output


async def main():
    order, modulus = find_safe_primes(1024)
    group = QuadraticResidue(modulus=modulus)
    sec_grp = SecureGroup(group)
    g = group.generator

    await mpc.start()

    sec_a = mpc._random(mpc.SecFld(modulus=group.order))
    sec_aG = repeat_public_base_secret_output(g, sec_a, sec_grp)
    pubkey = await mpc.output(sec_aG)

    # Encode input m to curve point; Also encode 0 (to Z), required for decoding in MPC
    m = 301
    print("Plaintext: ", m)
    print("Binary representation of secret input: ", bin(m))
    # m_encoded, z_encoded = encode_to_schnorr(m, group)
    m_encoded = encode_v3(m, group)

    # Encrypt m using ElGamal under pubkey (k is dummy random)
    k = 3809836274351126414438016410872729431457283748709071912489406741072364988984
    elgamal_m = encrypt(m_encoded, k, g, pubkey)
    # elgamal_z = encrypt(z_encoded, k, g, pubkey)

    # # Decrypt with MPC to shares of m (sec_m)
    sec_m = await threshold_decrypt(elgamal_m, sec_a)
    # sec_z = await threshold_decrypt(elgamal_z, sec_a)
    # sec_m = await decode_from_schnorr(sec_m, sec_z)
    sec_m = decode_v3(sec_m)

    # Show that decrypted ciphertext == plaintext (i.e. [m] recombines to m)
    plain_m = await mpc.output(sec_m)
    print("Plaintext after recombining shares from threshold decryption: ", plain_m)
    assert int(plain_m) == m
    share_m = await mpc.gather(sec_m)
    print(f"MPC party {mpc.pid} holds share of m: {share_m}")

    await mpc.shutdown()


if __name__ == "__main__":
    mpc.run(main())


# TODO:
# 1: encode_to_schnorr / encode_to_curve: make uniform API
# 1. Enable all regular MPC methods on secgrp types
