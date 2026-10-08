"""Modified-symbol lower bound for the Weil form at scale lambda (exploratory, mpmath).

Frequency form of the prime-side Weil form on even f supported in [-L, L]:
    W(f) = (1/2pi) int |f^(r)|^2 sigma(r) dr + 2 f^(i/2)^2,
    sigma(r) = Re psi(1/4 + i r/2) - log pi - 2 sum_{log n < 2L} Lambda(n)/sqrt(n) cos(r log n).
For |r| >= Omega:  sigma(r) >= m_out := Re psi(1/4 + i Omega/2) - log pi - 2 sum Lambda(n)/sqrt(n)  (> 0 for Omega large).
Replace sigma by sigma~ = sigma on |r| < Omega, m_out outside:  W >= B,
    B = Pi M_{sigma 1_in} Pi + m_out (I - G_in) + pole,   G_in = Pi M_{1_in} Pi  (Slepian operator).
In an orthonormal basis of L^2[-L, L] the in-band parts have entries bounded by (max|sigma| + m_out) sqrt(l_j l_k),
l_k = in-band energy of the k-th basis function, so B is "diagonal = m_out" on all basis functions whose in-band
energy is negligible, and B >= 0 reduces (Schur complement) to positivity of a finite block with margin.

Basis: orthonormal even Legendre polynomials p_k(t) = sqrt((2k+1)/(2L)) P_k(t/L), k = 0, 2, ..., 2K-2;
p_k^(r) = sqrt((2k+1) L / 2) * 2 i^k j_k(r L).
Usage: symbol_B.py lam Omega K [dps] [nodes]
"""
import sys, time
import mpmath as mp


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
            out.append((n, mp.log(p)))
    return out


def run(lam, Omega, K, dps=30, nodes=None):
    mp.mp.dps = dps
    L = mp.log(lam)
    vm = [(n, Lam) for n, Lam in von_mangoldt(int(lam * lam)) if mp.log(n) < 2 * L]
    S = mp.fsum(Lam / mp.sqrt(n) for n, Lam in vm)
    a_Om = mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * mp.mpf(Omega) / 2)) - mp.log(mp.pi)
    m_out = a_Om - 2 * S
    sigma = lambda r: mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * r / 2)) - mp.log(mp.pi) - 2 * mp.fsum(Lam / mp.sqrt(n) * mp.cos(r * mp.log(n)) for n, Lam in vm)
    # Gauss-Legendre on [0, Omega]
    nodes = nodes or int(6 * Omega) + 200
    xs, ws = mp.gauss_legendre(nodes) if hasattr(mp, 'gauss_legendre') else (None, None)
    if xs is None:
        from mpmath.calculus.quadrature import GaussLegendre
        gl = GaussLegendre(mp.mp)
        # degree m gives 3*2^(m-1) nodes
        m = 1
        while 3 * 2 ** (m - 1) < nodes:
            m += 1
        pts = gl.calc_nodes(m, mp.mp.prec)
        xs = [p[0] for p in pts]; ws = [p[1] for p in pts]
    rs = [Omega / 2 * (x + 1) for x in xs]; wr = [Omega / 2 * w for w in ws]
    ks = [2 * i for i in range(K)]
    t0 = time.time()
    # phat[k][q] = p_k^(r_q) (real, even): sqrt((2k+1) L/2) * 2 * (-1)^(k/2) * j_k(r L)
    phat = []
    for k in ks:
        pref = mp.sqrt((2 * k + 1) * L / 2) * 2 * (-1) ** (k // 2)
        row = []
        for r in rs:
            x = r * L
            jk = mp.sqrt(mp.pi / (2 * x)) * mp.besselj(k + mp.mpf(1) / 2, x)
            row.append(pref * jk)
        phat.append(row)
    sig = [sigma(r) for r in rs]
    T = mp.matrix(K, K); G = mp.matrix(K, K)
    for a in range(K):
        for b in range(a, K):
            t = mp.fsum(wr[q] * phat[a][q] * phat[b][q] * sig[q] for q in range(len(rs))) / mp.pi   # (1/2pi) * 2 * int_0^Omega
            g = mp.fsum(wr[q] * phat[a][q] * phat[b][q] for q in range(len(rs))) / mp.pi
            T[a, b] = T[b, a] = t; G[a, b] = G[b, a] = g
    # pole: 2 <p_j, cosh(t/2)> <p_k, cosh(t/2)>  (f^(i/2) = int f cosh(t/2))
    chi = [mp.quad(lambda t: mp.sqrt((2 * k + 1) / (2 * L)) * mp.legendre(k, t / L) * mp.cosh(t / 2), [-L, L]) for k in ks]
    B = mp.matrix(K, K)
    for a in range(K):
        for b in range(K):
            B[a, b] = T[a, b] + m_out * ((1 if a == b else 0) - G[a, b]) + 2 * chi[a] * chi[b]
    E = mp.eigsy(B, eigvals_only=True)
    Emin = min(E)
    # time-domain check of the full W on the same basis is not done here; report the in-band energies
    ell = [G[a, a] for a in range(K)]
    print(f"lambda = {lam}, Omega = {Omega}, K = {K} even Legendre modes (degrees 0..{ks[-1]}), dps = {dps}, nodes = {len(rs)}, {time.time()-t0:.0f}s")
    print(f"  primes: {[n for n, _ in vm]}, 2*sum Lambda/sqrt n = {mp.nstr(2*S, 6)}, a(Omega) = {mp.nstr(a_Om, 6)}, m_out = {mp.nstr(m_out, 6)}")
    print(f"  sigma(0) = {mp.nstr(sigma(mp.mpf(0)), 6)}, max|sigma| on [0,Omega] ~ {mp.nstr(max(abs(s) for s in sig), 5)}")
    print(f"  smallest eigenvalue of B_PP: {mp.nstr(Emin, 8)}  (log10 |.| = {mp.nstr(mp.log10(abs(Emin)), 5)}),  next: {mp.nstr(sorted(E)[1], 6)}, {mp.nstr(sorted(E)[2], 6)}")
    print(f"  sum of in-band energies (trace G_in on block) = {mp.nstr(mp.fsum(ell), 6)}  vs 2 Omega L / pi = {mp.nstr(2*Omega*L/mp.pi, 6)}")
    print("  in-band energy l_k by degree k:")
    for a in range(0, K, max(1, K // 12)):
        print(f"    k = {ks[a]:3d}: {mp.nstr(ell[a], 4)}")
    print(f"    k = {ks[-1]:3d}: {mp.nstr(ell[-1], 4)}")
    return Emin


if __name__ == "__main__":
    lam = int(sys.argv[1]); Omega = int(sys.argv[2]); K = int(sys.argv[3])
    dps = int(sys.argv[4]) if len(sys.argv) > 4 else 30
    nodes = int(sys.argv[5]) if len(sys.argv) > 5 else None
    run(lam, Omega, K, dps, nodes)
