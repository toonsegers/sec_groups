from math import gcd
import random
import os
import sys
import time

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from mpyc.runtime import mpc
from sec_groups.tools.secgcd import (
    secure_binary_xgcd,
    secure_xgcd,
    secure_division,
    secgcdKnuthAlgB,
    secure_gcd,
    hensel_inverse_mod2k,
    secure_montgomery_exponentiation,
    montgomery_exponentiation,
)


async def demo_secure_xgcd(instances=5):
    verifications = {}

    await mpc.start()
    secint = mpc.SecInt()

    t0 = time.time()
    print("Demo secure_binary_xgcd().")
    for i in range(instances):
        # f = await mpc.transfer(random.randrange(-10000, 10000, 1), senders=0)
        # g = await mpc.transfer(random.randrange(-10000, 10000, 1), senders=0)
        f = await mpc.transfer(random.randrange(0, 1000), senders=0)
        g = await mpc.transfer(random.randrange(0, 1000), senders=0)
        print("f=", f)
        print("g=", g)
        gv, a, b = secure_binary_xgcd(secint(f), secint(g), progress_bar=True)
        print("Evaluating circuit, running MPC rounds... (no progress bar shown)")
        v, a, b = await mpc.output([gv, a, b])
        assert gcd(f, g) == v
        print("a*f + b*g=", a * f + b * g)
        assert a * f + b * g == v, (a, b, secint.field.modulus)
        verifications[f"secure_binary_xgcd(); instance:{i}; inputs:{f},{g}"] = (
            a * f + b * g == v
        )
        print("Secure extended GCD computed correctly!")

    t1 = time.time()
    print(f"Time to run {instances} instances of secure_binary_xgcd: {t1-t0}")

    print("Demo secure_xgcd().")
    for i in range(instances):
        f = await mpc.transfer(random.randrange(-2**30, 2**30), senders=0)
        g = await mpc.transfer(random.randrange(-2**30, 2**30), senders=0)        
        # f = await mpc.transfer(random.randrange(0, 10000, 1), senders=0)
        # g = await mpc.transfer(random.randrange(0, 10000, 1), senders=0)
        print("f=", f)
        print("g=", g)
        v, a, b = secure_xgcd(secint(f), secint(g))
        vab_out = await mpc.output([v, a, b])
        v, a, b = vab_out
        print("secure v, a, b=")
        print(v, a, b)
        assert gcd(f, g) == v
        print("a*f + b*g=", a * f + b * g)
        assert a * f + b * g == v
        assert int(a) * f + int(b) * g == int(v)
        print("Secure extended GCD computed correctly!")
        verifications[f"secure_xgcd(); instance:{i}; inputs: {f},{g}"] = int(
            a
        ) * f + int(b) * g == int(v)
    t2 = time.time()

    print("Total time secure_binary_xgcd:", t1 - t0)    
    print("Total time secure_xgcd       :", t2 - t1)

    await mpc.shutdown()
    return verifications


async def demo_secure_division(instances=2):
    verifications = {}
    await mpc.start()
    secint = mpc.SecInt()

    print("Demo Dahl2012 integer division approximation.")
    rounding_errors = 0
    for i in range(instances):
        d_max = 2 ** 14
        n_max = 2 ** 16

        d = await mpc.transfer(random.randrange(1, d_max, 1), senders=0)
        n = await mpc.transfer(random.randrange(1, n_max, 1), senders=0)
        print("Inputs (n = numerator, d = denominator)", n)
        print("n=", n)
        print("d=", d)

        q, r = secure_division(secint(n), secint(d))
        qo, ro = await mpc.output([q, r])
        print("q=", qo)
        print("r=", ro)
        print("n // d=", n // d)
        assert n // d == qo
        print("Secure division computed correctly!")
        verifications[f"secure_division(); instance:{i}; inputs: {n},{d}"] = (
            n // d == qo
        )

    await mpc.shutdown()
    return verifications


async def demo_gcd_other(instances=2):
    verifications = {}
    await mpc.start()
    secint = mpc.SecInt()
    print("Demo Knuth-based gcds:")
    for i in range(instances):
        f = await mpc.transfer(random.randrange(1, 10000), senders=0)
        g = await mpc.transfer(random.randrange(1, 10000), senders=0)
        v = await mpc.output(secgcdKnuthAlgB(secint(f), secint(g)))
        assert gcd(f, g) == v, (f, g)
        print("Computed correctly!")
        verifications[f"secgcdKnuthAlgB(); instance:{i}; inputs: {f},{g}"] = (
            gcd(f, g) == v
        )
    print()

    print("Demo secure_gcd(), i.e. not _extended_ gcd")
    for i in range(instances):
        f = await mpc.transfer(random.randrange(1, 10000), senders=0)
        g = await mpc.transfer(random.randrange(1, 10000), senders=0)
        v = secure_gcd(secint(f), secint(g))
        vo = await mpc.output(v)
        assert gcd(f, g) == vo
        print("Computed correctly!")
        verifications[f"secure_gcd(); instance:{i}; inputs: {f},{g}"] = gcd(f, g) == vo

    await mpc.shutdown()
    return verifications


async def demo_secure_montgomery_exponentiation(instances=2):
    verifications = {}
    await mpc.start()
    print("Test secure Montgomery exponentiation.")
    for i in range(instances):
        a = await mpc.transfer(random.randrange(1, 10000, 1), senders=0)
        n = await mpc.transfer(random.randrange(1, 10000, 1), senders=0)
        N = await mpc.transfer(random.randrange(1, 10000, 2), senders=0)

        l = max(a.bit_length(), N.bit_length())
        R = 2 ** l
        N_prime = hensel_inverse_mod2k(N, l) * -1
        print("a=", a)
        print("n=", n)
        print("N=", N)
        print("R=", R)
        assert R > N
        assert (0 <= a) and (a < N * R)

        print(f"N_prime (-({N}^-1) mod {R}) =", N_prime)
        print()

        print("Test non-secure montgomery_exponentiation:")
        r = montgomery_exponentiation(a, n, R, N, N_prime)
        print("Result of Montgomery exponentiation=", r)
        assert r % N == (a ** n) % N
        print("Direct calculation of a^n mod N=", a ** n % N)
        print("Correct.")
        print()

        print("Test secure_montgomery_exponentiation:")
        secint = mpc.SecInt()
        seca = secint(a)
        secN = secint(N)
        secN_prime = -hensel_inverse_mod2k(secN, l)
        secN_primeo = await mpc.output(secN_prime)
        print("secN_prime=", secN_primeo)
        assert N_prime == secN_primeo
        print("Secure computation of N_prime correct.")
        print(f"a = {a}, n = {n}, modulus N = {N}, aux modulus R = {R}")
        r = await secure_montgomery_exponentiation(seca, n, R, secN, secN_prime)
        out = await mpc.output(r)
        print("Result of secure Montgomery exponentiation=", out)
        assert out % N == (a ** n) % N
        print(f"Secure Montgomery exponentiation computed correctly.")
        verifications[
            f"secure_montgomery_exponentiation(); instance:{i}; inputs: {a},{n},{N}"
        ] = (out % N == (a ** n) % N)

    await mpc.shutdown()
    return verifications


if __name__ == "__main__":
    verifications = mpc.run(demo_secure_xgcd())
    print(verifications)
    verifications = mpc.run(demo_secure_division())
    print(verifications)
    verifications = mpc.run(demo_gcd_other())
    print(verifications)
    verifications = mpc.run(demo_secure_montgomery_exponentiation())
    print(verifications)
