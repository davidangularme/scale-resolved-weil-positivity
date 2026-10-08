"""Compute the first 500 zeta zeros (Arb) and three control point sets with the same smooth density."""
import numpy as np, mpmath as mp
from flint import acb, ctx
ctx.dps = 60
with open('zeros500.txt', 'w') as f:
    for z in acb.zeta_zeros(1, 500):
        f.write(z.imag.str(55, radius=False) + '\n')
rng = np.random.default_rng(20261007)
mp.mp.dps = 40
zeros = [mp.mpf(l) for l in open('zeros500.txt')]
Nsmooth = lambda T: T / (2 * mp.pi) * mp.log(T / (2 * mp.pi * mp.e)) + mp.mpf(7) / 8
Ninv = lambda x: mp.findroot(lambda T: Nsmooth(T) - x, 14 + x * 1.5)
def write(name, xs):
    with open(name, 'w') as f:
        for x in xs:
            f.write(mp.nstr(Ninv(mp.mpf(x)), 45) + '\n')
start, n = float(Nsmooth(zeros[0])), 500
write('pts_fence.txt', [start + j for j in range(n)])
sp = rng.exponential(1.0, n)
write('pts_poisson.txt', list(start + np.concatenate([[0], np.cumsum(sp[:-1])])))
M = 1400
A = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M)); H = (A + A.conj().T) / 2
x = np.linalg.eigvalsh(H); x = x / np.max(np.abs(x))
uf = M * (0.5 + (x * np.sqrt(1 - x ** 2) + np.arcsin(x)) / np.pi)
mid = uf[M // 2 - 250:M // 2 + 250]; mid = (mid - mid[0]) * (499 / (mid[-1] - mid[0]))
write('pts_gue.txt', list(start + mid))
