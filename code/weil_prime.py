"""Prime-side (non-circular) Weil quadratic form in the archive's cosine basis.

Basis: g_j(t) = cos(j t) on [-L, L] (t = log u), j = 0..N, L = log(lambda).
For h = g_j * g_k (autocorrelation; g real even), Weil's explicit formula gives

  sum_gamma hhat(gamma) = hhat(i/2) + hhat(-i/2) - log(pi) h(0)
                          + (1/2pi) int hhat(r) Re psi(1/4 + i r/2) dr
                          - 2 sum_n Lambda(n) n^{-1/2} h(log n)

The right-hand side uses only primes and the archimedean place. No zeta zero is used.
W_jk := right-hand side.  Under RH, W = sum over ALL zeros of ghat_j ghat_k (PSD).
The archive's matrix M is the same sum truncated to the first 500 zeros.
"""
import sys, json, time
import mpmath as mp
from flint import arb, arb_mat, acb_mat, ctx


def von_mangoldt_upto(X):
    out = []
    for n in range(2, X + 1):
        m, p = n, None
        for q in range(2, n + 1):
            if m % q == 0:
                p = q
                break
        while m % p == 0:
            m //= p
        if m == 1:
            out.append((n, mp.log(p)))
    return out


def h_jk(j, k, x, L):
    """(g_j * g_k)(x) = int g_j(t) g_k(t - x) dt for 0 <= x <= 2L."""
    lo, hi = x - L, L
    def F(a, b):
        if a == 0:
            return (hi - lo) * mp.cos(b)
        return (mp.sin(a * hi + b) - mp.sin(a * lo + b)) / a
    return (F(j - k, k * x) + F(j + k, -k * x)) / 2


def build_W(N, lam, dps, gl_degree):
    mp.mp.dps = dps
    L = mp.log(lam)
    twoL = 2 * L
    # pole terms
    c = mp.mpf(1) / 2
    a = []
    for j in range(N + 1):
        a.append(2 * (c * mp.sinh(c * L) * mp.cos(j * L) + j * mp.cosh(c * L) * mp.sin(j * L)) / (c * c + j * j))
    # prime terms
    lam2 = int(mp.floor(mp.e ** twoL + mp.mpf('1e-30')))
    vm = [(n, Lam) for (n, Lam) in von_mangoldt_upto(lam2) if mp.log(n) < twoL]
    # Gauss-Legendre nodes on [0, 4L] for the archimedean integral
    gl = mp.calculus.quadrature.GaussLegendre(mp.mp)
    nodes = gl.calc_nodes(gl_degree, mp.mp.prec)
    T = 4 * L
    pts = [((x + 1) * T / 2, w * T / 2) for (x, w) in nodes]
    pref = [(mp.exp(-t) / (1 - mp.exp(-t)), mp.exp(-t / 4) / (1 - mp.exp(-t)), w, t) for (t, w) in pts]
    tail = -mp.log(1 - mp.exp(-T))
    W = mp.matrix(N + 1, N + 1)
    for j in range(N + 1):
        for k in range(j, N + 1):
            h0 = h_jk(j, k, mp.mpf(0), L)
            arch = -mp.euler * h0 + tail * h0
            s = mp.mpf(0)
            for (e1, e4, w, t) in pref:
                s += w * (e1 * h0 - e4 * h_jk(j, k, t / 2, L))
            arch += s
            primes = mp.fsum(Lam / mp.sqrt(n) * h_jk(j, k, mp.log(n), L) for (n, Lam) in vm)
            val = 2 * a[j] * a[k] - mp.log(mp.pi) * h0 + arch - 2 * primes
            W[j, k] = val
            W[k, j] = val
    return W, vm


def ghat(j, r, L):
    if j == 0:
        return 2 * mp.sin(r * L) / r
    return mp.sin((j - r) * L) / (j - r) + mp.sin((j + r) * L) / (j + r)


def build_M(N, lam, zeros):
    L = mp.log(lam)
    Phi = [[ghat(j, g, L) for g in zeros] for j in range(N + 1)]
    M = mp.matrix(N + 1, N + 1)
    for j in range(N + 1):
        for k in range(j, N + 1):
            v = 2 * mp.fsum(Phi[j][i] * Phi[k][i] for i in range(len(zeros)))
            M[j, k] = v
            M[k, j] = v
    return M


def eig_sym(A):
    n = A.rows
    B = acb_mat([[mp.nstr(A[i, j], mp.mp.dps) for j in range(n)] for i in range(n)])
    ev = B.eig(algorithm="approx")
    vals = sorted([mp.mpf(e.real.mid().str(mp.mp.dps, radius=False)) for e in ev])
    return vals
