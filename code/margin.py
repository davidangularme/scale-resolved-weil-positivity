"""Positivity margin m(lambda) = min_f W_lambda(f)/||f||^2 (Galerkin upper bound, L2-normalised),
and first-order rigidity of each prime weight: eps_n = m / |T_n(f0)| where f0 is the minimiser and
T_n(f0) the n-th prime term evaluated on f0 (so W(f0) changes sign when Lambda(n) -> (1+eps_n) Lambda(n)).
Basis: cos(kt) on [-L,L]; Gram G_jk = int cos(jt)cos(kt) dt.  Generalised eigenproblem W a = mu G a
solved via Cholesky of G.
"""
import sys, json, time
import mpmath as mp
from flint import ctx
import weil_prime as wp


def gram(N, L):
    G = mp.matrix(N + 1, N + 1)
    for j in range(N + 1):
        for k in range(j, N + 1):
            if j == k == 0:
                v = 2 * L
            elif j == k:
                v = L + mp.sin(2 * j * L) / (2 * j)
            else:
                v = mp.sin((j - k) * L) / (j - k) + mp.sin((j + k) * L) / (j + k)
            G[j, k] = v; G[k, j] = v
    return G


def gen_min(W, G):
    C = mp.cholesky(G)
    Ci = mp.inverse(C)
    A = Ci * W * Ci.T
    E, Q = mp.eigsy(A)
    i = min(range(len(E)), key=lambda t: E[t])
    y = Q[:, i]
    a = Ci.T * y          # eigenvector in coefficient space, with a^T G a = 1
    return E[i], a


if __name__ == "__main__":
    lam = float(sys.argv[1]); N = int(sys.argv[2]); dps = int(2.6 * N) + 80
    ctx.dps = dps; mp.mp.dps = dps
    t0 = time.time()
    W, vm = wp.build_W(N, lam, dps, 8)
    L = mp.log(lam)
    G = gram(N, L)
    m, a = gen_min(W, G)
    out = dict(lam=lam, N=N, L=float(L), margin=float(m), log10_margin=float(mp.log10(abs(m))), sign=('+' if m > 0 else '-'), secs=round(time.time() - t0, 1), primes=[])
    for n, Lam in vm:
        w = Lam / mp.sqrt(n)
        Tn = mp.fsum(a[j] * a[k] * (-2 * w * wp.h_jk(j, k, mp.log(n), L)) for j in range(N + 1) for k in range(N + 1))
        out['primes'].append(dict(n=n, T_f0=float(Tn), eps_star=float(m / abs(Tn)) if Tn != 0 else None, log10_eps=float(mp.log10(abs(m / Tn))) if Tn != 0 else None))
    print(json.dumps(out), flush=True)
