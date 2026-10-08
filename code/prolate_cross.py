"""Low/high coupling of the prime operator in the prolate basis at lambda = 2.

Prolate basis psi_n (bandwidth Omega, interval [-L, L], L = log 2), orthonormal on [-L, L].
Prime operator  (P f)(t) = -2 sum_n Lambda(n)/sqrt(n) * [f(t - u_n) + f(t + u_n)]/2 restricted to [-L, L], u_n = log n,
so that <f, P f> is the prime term of the Weil form.  We compute P_mn = <psi_m, P psi_n> and the spectral norm of the
low/high block P_LH for cuts K, together with the prolate-basis Galerkin lower bound mu_min(K) of the full form on
span(psi_0..psi_{K-1}) (even modes) obtained by re-expanding the certified cosine matrix is NOT needed here: we only
quantify the coupling, which decides whether a block + tail (Schur complement) argument can work at the margin m(2).
Double precision throughout (the coupling is O(1), no cancellation issue).
"""
import sys
import numpy as np
from numpy.polynomial import legendre as leg

L = np.log(2.0)


def prolate_basis(c, parity, M):
    ks = np.array([parity + 2 * i for i in range(M)])
    A = np.zeros((M, M))
    for i, k in enumerate(ks):
        A[i, i] = k * (k + 1) + c ** 2 * (2 * k * (k + 1) - 1) / ((2 * k + 3) * (2 * k - 1))
        if i + 1 < M:
            v = c ** 2 * (k + 2) * (k + 1) / ((2 * k + 3) * np.sqrt((2 * k + 1) * (2 * k + 5)))
            A[i, i + 1] = A[i + 1, i] = v
    chi, D = np.linalg.eigh(A)
    # coefficient vectors in normalised Legendre on [-1,1]: psi(x) = sum d_k sqrt(k+1/2) P_k(x)
    coefs = []
    for col in range(M):
        cvec = np.zeros(ks[-1] + 1)
        cvec[ks] = D[:, col] * np.sqrt(ks + 0.5)
        coefs.append(cvec)
    return chi, coefs


def evaluate(coefs, t):
    """psi_n(t) on [-L, L] orthonormal: psi_n(t) = (1/sqrt L) * series(t / L)."""
    x = t / L
    out = np.zeros((len(coefs), len(t)))
    inside = np.abs(x) <= 1
    for n, cv in enumerate(coefs):
        out[n, inside] = leg.legval(x[inside], cv) / np.sqrt(L)
    return out


def translation_matrix(coefs, u, nq=1200):
    """T[m, n] = int psi_m(t) psi_n(t - u) dt over the overlap [u - L, L] (u >= 0)."""
    a, b = u - L, L
    if b <= a:
        return np.zeros((len(coefs), len(coefs)))
    x, w = leg.leggauss(nq)
    t = (b - a) / 2 * x + (b + a) / 2
    w = (b - a) / 2 * w
    Pm = evaluate(coefs, t)
    Pn = evaluate(coefs, t - u)
    return (Pm * w) @ Pn.T


def von_mangoldt(nmax):
    out = []
    for n in range(2, nmax + 1):
        m, p = n, None
        for q in range(2, n + 1):
            if m % q == 0:
                p = q; break
        while m % p == 0:
            m //= p
        if m == 1:
            out.append((n, np.log(p)))
    return out


if __name__ == "__main__":
    Omega = float(sys.argv[1]) if len(sys.argv) > 1 else 80.0
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 90
    c = Omega * L
    chi, coefs = prolate_basis(c, 0, M)      # even modes only (the form acts on even functions)
    nmodes = len(coefs)
    P = np.zeros((nmodes, nmodes))
    for n, Lam in von_mangoldt(4):
        u = np.log(n)
        if u >= 2 * L:
            continue
        T = translation_matrix(coefs, u)
        P += -2 * Lam / np.sqrt(n) * (T + T.T) / 2
    # sanity: orthonormality
    x, w = leg.leggauss(1500); t = L * x; w = L * w
    Psi = evaluate(coefs, t)
    G = (Psi * w) @ Psi.T
    print(f"Omega = {Omega}, c = {c:.3f}, even modes M = {M}; max |G - I| = {np.abs(G - np.eye(nmodes)).max():.1e}")
    print(f"||P|| = {np.linalg.norm(P, 2):.4f}   (bound 2*sum Lambda(n)/sqrt(n) over log n < 2L = "
          f"{2 * sum(Lam / np.sqrt(n) for n, Lam in von_mangoldt(4) if np.log(n) < 2 * L):.4f})")
    print(" K (even-mode index, cut after mode K-2) : ||P_LH||_2   max|P_LH|")
    for K in [10, 20, 28, 36, 44, 56, 64, 72, 80]:
        kk = K // 2
        if kk >= nmodes:
            break
        B = P[:kk, kk:]
        print(f"   K = {K:3d} : {np.linalg.norm(B, 2):.4e}   {np.abs(B).max():.3e}")
    # decay of the coupling of the SMOOTH modes (n <= 2c/pi) to high modes: row norms
    k0 = int(2 * c / np.pi) // 2
    print(f" row norms of P[m, high] for high = modes with index >= {2 * (k0 + 10)} (well past the plunge):")
    for m in [0, 2, 4, 8, 12, 2 * k0 - 4, 2 * k0, 2 * k0 + 4, 2 * k0 + 8, 2 * k0 + 16]:
        mm = m // 2
        if mm < k0 + 10:
            print(f"   mode {m:3d}: {np.linalg.norm(P[mm, k0 + 10:]):.3e}")
