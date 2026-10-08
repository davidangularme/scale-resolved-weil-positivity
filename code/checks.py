"""Corrections 1-2: rigidity-matrix zero mode, transverse Hessian sign, generic Fisher pole identity,
and a check of the Weil explicit formula (time-domain archimedean term) on a Gaussian."""
import mpmath as mp
mp.mp.dps = 60
g = [mp.mpf(l) for l in open('zeros500.txt')][:45]; N = len(g)
R = mp.matrix(N, N)
for n in range(N):
    for m in range(N):
        if n != m: R[n, m] = -1 / (g[n] - g[m]) ** 2
    R[n, n] = -sum(R[n, m] for m in range(N) if m != n)
ev = sorted(mp.eigsy(R, eigvals_only=True))
print("rigidity matrix: mu0 =", mp.nstr(ev[0], 5), " mu1 =", mp.nstr(ev[1], 8))
def E(x, y):
    return -mp.fsum(mp.log(mp.sqrt((x[n]-x[m])**2 + (g[n]+y[n]-g[m]-y[m])**2)) for n in range(N) for m in range(n+1, N))
h = mp.mpf('1e-12'); z = [mp.mpf(0)] * N
def hess(d, i, j):
    def sh(a, b):
        v = [mp.mpf(0)] * N; v[i] += a; v[j] += b
        return E(v, z) if d == 'x' else E(z, v)
    return (sh(h, h) - sh(h, -h) - sh(-h, h) + sh(-h, -h)) / (4 * h * h)
print("H_transverse[0,0] =", mp.nstr(hess('x', 0, 0), 8), " H_longitudinal[0,0] =", mp.nstr(hess('y', 0, 0), 8), " R[0,0] =", mp.nstr(R[0, 0], 8))
a, b = mp.mpc(0.8, 20), mp.mpc(0.2, 20)
f = lambda s: (s - a) * (s - b) * mp.exp(s)
s = a + mp.mpf('1e-8') * mp.expjpi(0.3)
print("off-line toy zero: |I|*|s-rho|^2 =", mp.nstr(abs(mp.diff(lambda t: mp.log(f(t)), s, 2)) * abs(s - a) ** 2, 12))
mp.mp.dps = 40
sg = mp.mpf('0.1'); hf = lambda x: mp.exp(-x**2 / (2 * sg**2)); hh = lambda r: sg * mp.sqrt(2 * mp.pi) * mp.exp(-sg**2 * r**2 / 2)
zeros = [mp.mpf(l) for l in open('zeros500.txt')]
lhs = 2 * mp.fsum(hh(t) for t in zeros)
arch = -mp.euler * hf(0) + mp.quad(lambda t: (mp.exp(-t) * hf(0) - mp.exp(-t / 4) * hf(t / 2)) / (1 - mp.exp(-t)), [0, 0.5, 2, 10, mp.inf])
primes = [p for p in range(2, 200) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
ps = mp.fsum(mp.log(p) / mp.sqrt(p ** k) * hf(mp.log(p ** k)) for p in primes for k in range(1, 20) if p ** k < 10 ** 6)
rhs = mp.re(hh(mp.j / 2) + hh(-mp.j / 2)) - mp.log(mp.pi) * hf(0) + arch - 2 * ps
print("explicit formula: zero side", mp.nstr(lhs, 25), " prime side", mp.nstr(rhs, 25))
