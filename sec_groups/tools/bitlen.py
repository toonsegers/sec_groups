"""Secure norm protocols by Thijs Veugen, adapted by Berry Schoenmakers.

See: https://www.researchgate.net/profile/Thijs_Veugen
"""

import itertools
from mpyc.runtime import mpc


def norm(self, x):
    """Recursive norm (adapted from mpc._norm())."""

    def _norm(x):
        n = len(x)
        if n == 1:
            t = x[0]
            return 1 - t, t

        i0, nz0 = _norm(x[:n//2])  # low bits
        i1, nz1 = _norm(x[n//2:])  # high bits
        i0 += (((n+1)//2))
        return self.if_else(nz1, [i1, nz1], [i0, nz0])

    l = type(x[0]).bit_length
    i, _ = _norm(x)
    return l - i


def bit_length_mpc_norm(self, a):
    """Bit length currently in MPyC."""
    x = self.to_bits(a)  # low to high bits
    return norm(self, x)


def bit_length_mpc_find(self, a):
    """Bit length currently in MPyC."""
    l = type(a).bit_length
    x = self.to_bits(a)  # low to high bits
    x.reverse()
    nf = self.find(x, 1)
    return l - nf


@mpc.coroutine
async def bit_length_new(self, a):
    stype = type(a)
    await self.returnType(stype)

    Zp = stype.field
    l = stype.bit_length
    r_bits = await self.random_bits(Zp, l)
    r_modl = 0
    for r_i in reversed(r_bits):
        r_modl <<= 1
        r_modl += r_i.value
    r_divl = self._random(Zp, 1<<self.options.sec_param).value
    a = await self.gather(a)
    c = await self.output(a + ((1<<l) + (r_divl << l) + r_modl))
    c = c.value % (1<<l)

    c_bits = [(c >> i) & 1 for i in range(l)]
    r_bits = [stype(r.value) for r in r_bits]
    d_bits = [1-r_bits[i] if c_bits[i] else r_bits[i] for i in range(l)]
    h_bits = mpc.schur_prod([1-r_bits[i] for i in range(l-1) if not c_bits[i]],
                            [d_bits[i+1] for i in range(l-1) if not c_bits[i]])
    for i in range(l-2, -1, -1):
        if not c_bits[i]:
            d_bits[i+1] = h_bits.pop()

    k = norm(self, d_bits) - 1
    k_u = self.unit_vector(k, l)  # 0<=k<l assumed
    k = mpc.in_prod(k_u, list(map(stype, range(l+1))))
    k2 = mpc.in_prod(k_u, list(map(stype, list(map(lambda a: 2**a, range(l+1))))))

    return k - (a < k2) + 1  # position + 1 is bit length


@mpc.coroutine
async def bit_length_integrated(self, a):
    stype = type(a)
    await self.returnType(stype)

    Zp = stype.field
    l = stype.bit_length
    r_bits = await self.random_bits(Zp, l)
    r_modl = 0
    for r_i in reversed(r_bits):
        r_modl <<= 1
        r_modl += r_i.value
    r_divl = self._random(Zp, 1<<self.options.sec_param).value
    a = await self.gather(a)
    c = await self.output(a + ((1<<l) + (r_divl << l) + r_modl))
    c = c.value % (1<<l)

    c_bits = [(c >> i) & 1 for i in range(l)]
    d_bits = [stype((1 - r_bits[i] if c_bits[i] else r_bits[i]).value) for i in range(l)]
    h_bits = mpc.schur_prod([stype(1 - r_bits[i]) for i in range(l-1) if not c_bits[i]],
                            [d_bits[i+1] for i in range(l-1) if not c_bits[i]])
    for i in range(l-2, -1, -1):
        if not c_bits[i]:
            d_bits[i+1] = h_bits.pop()

    k = norm(self, d_bits) - 1
    k_u = self.unit_vector(k, l)  # 0<=k<l assumed
    k_u = await mpc.gather(k_u)
    psums = list(itertools.accumulate(k_u))
    pp = await mpc.schur_prod(psums, [c_bits[i] - r_bits[i] for i in range(l)])
    for i in range(l):
        r_bits[i] += pp[i]

    s_sign = (await self.random_bits(Zp, 1, signed=True))[0].value
    e = [None] * (l+1)
    sumXors = 0
    for i in range(l-1, -1, -1):
        c_i = c_bits[i]
        r_i = r_bits[i].value
        e[i] = Zp(s_sign + r_i - c_i + 3*sumXors)
        sumXors += 1 - r_i if c_i else r_i
    e[l] = Zp(s_sign - 1 + 3*sumXors)
    g = await self.is_zero_public(stype(self.prod(e)))
    z = Zp(1 - s_sign if g else 1 + s_sign)/2
    return k - z + 1  # position + 1 is bit length


async def test(text, bit_length):
    print(f'Secure bit length: using {text}.')
    async with mpc:
        for i in range(1, 128):  # TODO: case i=0, case i<0
            a = i
            n = await mpc.output(bit_length(mpc, secint(a)))
            print(f'{round(100*i/127)}%', end='\r')
            assert a.bit_length() == n, (a.bit_length(), n, i)


# secint = mpc.SecInt()
# print(f'Using secure {secint.bit_length}-bit integers: {secint.__name__}')

# mpc.run(test('MPyC bit decomposition (norm)', bit_length_mpc_norm))
# mpc.run(test('MPyC bit decomposition (find)', bit_length_mpc_find))
# mpc.run(test('new O(m) approach', bit_length_new))
# mpc.run(test('integrated O(m) approach', bit_length_integrated))
