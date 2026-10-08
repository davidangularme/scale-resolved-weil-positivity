"""Prolate spheroidal wave function eigenvalues lambda_n(c) for the sinc kernel on [-1,1]
(Slepian): int_{-1}^{1} sin(c(x-t))/(pi(x-t)) psi_n(t) dt = lambda_n psi_n(x),  sum_n lambda_n = 2c/pi.

Method (Bouwkamp / Xiao-Rokhlin-Yarvin): expand psi_n in normalised Legendre polynomials; the
differential operator is a symmetric tridiagonal matrix on each parity class.  Then
mu_n psi_n(0) = int psi_n = sqrt(2) d_0                 (even n)
mu_n psi_n'(0) = i c int t psi_n = i c sqrt(2/3) d_1     (odd n)
and lambda_n = c |mu_n|^2 / (2 pi).

Use: prolate.py Omega [Omega ...]  with c = Omega * log 2 (the lambda = 2 window, L = log 2).
Reports K(1e-20), K(1e-30): the first index with lambda_K below the threshold (even modes only,
since the Weil form acts on even test functions), and a(Omega) = Re psi(1/4 + i Omega/2) - log pi.
"""
import sys
import mpmath as mp

mp.mp.dps = 50


def legendre0(kmax):
    """P_k(0) and P_k'(0) for k = 0..kmax (exact recurrences)."""
    P = [mp.mpf(1), mp.mpf(0)]
    dP = [mp.mpf(0), mp.mpf(1)]
    for k in range(1, kmax):
        P.append(-k * P[k - 1] / (k + 1))                 # (k+1)P_{k+1}(0) = -k P_{k-1}(0)
        dP.append(((2 * k + 1) * P[k] - k * dP[k - 1]) / (k + 1) + 0)  # placeholder, fixed below
    # P'_{k+1}(0) = (2k+1) P_k(0) + P'_{k-1}(0)
    dP = [mp.mpf(0), mp.mpf(1)]
    for k in range(1, kmax):
        dP.append((2 * k + 1) * P[k] + dP[k - 1])
    return P, dP


def prolate_eigs(c, parity, M):
    """Eigenvalues lambda_n (n of given parity, ascending) from an M x M tridiagonal block."""
    c = mp.mpf(c)
    ks = [parity + 2 * i for i in range(M)]
    A = mp.matrix(M, M)
    for i, k in enumerate(ks):
        A[i, i] = k * (k + 1) + c ** 2 * (2 * k * (k + 1) - 1) / ((2 * k + 3) * (2 * k - 1))
        if i + 1 < M:
            v = c ** 2 * (k + 2) * (k + 1) / ((2 * k + 3) * mp.sqrt((2 * k + 1) * (2 * k + 5)))
            A[i, i + 1] = v; A[i + 1, i] = v
    chi, D = mp.eigsy(A)
    order = sorted(range(M), key=lambda i: chi[i])
    P0, dP0 = legendre0(ks[-1] + 1)
    lams = []
    for col in order:
        d = [D[i, col] for i in range(M)]
        if parity == 0:
            psi0 = mp.fsum(d[i] * mp.sqrt(k + mp.mpf(1) / 2) * P0[k] for i, k in enumerate(ks))
            mu = mp.sqrt(2) * d[0] / psi0
        else:
            dpsi0 = mp.fsum(d[i] * mp.sqrt(k + mp.mpf(1) / 2) * dP0[k] for i, k in enumerate(ks))
            mu = c * mp.sqrt(mp.mpf(2) / 3) * d[0] / dpsi0
        lams.append(c * mu ** 2 / (2 * mp.pi))
    return lams


def run(Omega, L=mp.log(2), M=None):
    c = Omega * L
    M = M or int(2 * c) + 60
    ev = prolate_eigs(c, 0, M)
    od = prolate_eigs(c, 1, M)
    total = mp.fsum(ev) + mp.fsum(od)
    allv = sorted(ev + od, reverse=True)
    K20 = next(2 * i for i, v in enumerate(ev) if v < mp.mpf('1e-20'))
    K30 = next(2 * i for i, v in enumerate(ev) if v < mp.mpf('1e-30'))
    Kall20 = next(i for i, v in enumerate(allv) if v < mp.mpf('1e-20'))
    a = mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * mp.mpf(Omega) / 2)) - mp.log(mp.pi)
    print(f"Omega = {Omega}: c = Omega*L = {mp.nstr(c, 8)}, 2c/pi = {mp.nstr(2 * c / mp.pi, 8)}, "
          f"sum lambda_n = {mp.nstr(total, 8)} (M = {M})")
    print(f"  even modes: K(1e-20) = {K20}  (lambda_K = {mp.nstr(ev[K20 // 2], 4)}),  "
          f"K(1e-30) = {K30}  (lambda_K = {mp.nstr(ev[K30 // 2], 4)})")
    print(f"  all modes:  number with lambda_n >= 1e-20: {Kall20}")
    print(f"  a(Omega) = Re psi(1/4 + i Omega/2) - log pi = {mp.nstr(a, 8)}   (log(Omega/2) - log pi = {mp.nstr(mp.log(Omega / 2) - mp.log(mp.pi), 8)})")
    print("  lambda_n (even n) around the plunge:")
    for i in range(max(0, K20 // 2 - 12), min(len(ev), K30 // 2 + 3), 2):
        print(f"    n = {2 * i:3d}: {mp.nstr(ev[i], 6)}")
    return dict(Omega=Omega, c=float(c), K20=K20, K30=K30, a=float(a))


if __name__ == "__main__":
    # sanity: c = 1, lambda_0 = 0.57258... ; c = 10, 2c/pi = 6.366
    e1 = prolate_eigs(1, 0, 20)
    print("check c=1: lambda_0 =", mp.nstr(e1[0], 10), "(ref 0.5725804...)")
    for arg in sys.argv[1:]:
        run(int(arg))
