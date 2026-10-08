"""Rigorous certificate of Weil positivity at scale lambda via the modified symbol (ball arithmetic, Arb).

Setting.  f even, supported in [-L, L], L = log lambda.  With sigma(r) = Re psi(1/4 + i r/2) - log pi
- 2 sum_{log n < 2L} Lambda(n) n^{-1/2} cos(r log n),
    W(f) = (1/2pi) int |f^(r)|^2 sigma(r) dr + 2 f^(i/2)^2.
Since Re psi(1/4 + i y) is increasing in |y|, sigma(r) >= m_out := Re psi(1/4 + i Omega/2) - log pi - 2 S for |r| >= Omega,
S = sum Lambda(n)/sqrt(n).  Hence W >= B with
    B(f) = (1/2pi) int_{|r|<Omega} |f^|^2 (sigma - m_out) dr + m_out ||f||^2 + 2 f^(i/2)^2.
Basis: orthonormal even Legendre polynomials p_k, k = 0, 2, ..., D, on [-L, L];  P = their span, Q = complement.
    p_k^(r) = sqrt((2k+1) L / 2) * 2 (-1)^{k/2} j_k(rL),   chi_k = <p_k, cosh(t/2)> = sqrt((2k+1) L/2) * 2 i_k(L/2).
B_PP computed rigorously (acb.integral).  Tail: the in-band energy of p_k is
    l_k = (2/pi) (2k+1) int_0^{Omega L} j_k(x)^2 dx <= (2/pi) (Omega L)^{2k+1} / ((2k+1)!!)^2    (|J_nu(x)| <= (x/2)^nu / Gamma(nu+1)),
so eps_Q := ||Q G_in Q|| <= sum_{k > D, even} l_k, and with sigma_max = max_{|r|<=Omega} |sigma|:
    ||B_PQ|| <= (sigma_max + m_out) sqrt(eps_Q) + 2 ||chi|| ||Q chi||,   B_QQ >= m_out - (sigma_max + m_out) eps_Q,
and (Schur)  B >= 0  <=  lambda_min(B_PP) >= ||B_PQ||^2 / lambda_min(B_QQ).
Usage: certify_symbol.py lam Omega D [dps] [tol_digits] [--pairs a:b]  (pairs subset for parallel runs, results cached in json)
"""
import sys, time, json, os
from flint import arb, acb, arb_mat, ctx
import mpmath as mp


def primes_below(lam, L):
    vm = []
    for n in range(2, int(lam * lam) + 1):
        m, p = n, None
        for q in range(2, n + 1):
            if m % q == 0:
                p = q; break
        while m % p == 0:
            m //= p
        if m == 1 and arb(n).log() < 2 * L:
            vm.append((n, arb(p).log()))
    return vm


def dfact_odd(n):   # n!! for odd n, exact integer
    r = 1
    for i in range(1, n + 1, 2):
        r *= i
    return r


def main():
    lam = int(sys.argv[1]); Omega = int(sys.argv[2]); D = int(sys.argv[3])
    dps = int(sys.argv[4]) if len(sys.argv) > 4 else 45
    told = int(sys.argv[5]) if len(sys.argv) > 5 else 32
    pairs_arg = None
    if '--pairs' in sys.argv:
        a, b = sys.argv[sys.argv.index('--pairs') + 1].split(':'); pairs_arg = (int(a), int(b))
    ctx.dps = dps
    L = arb(lam).log()
    ks = list(range(0, D + 1, 2)); K = len(ks)
    vm = primes_below(lam, L)
    S = sum((Lam / arb(n).sqrt() for n, Lam in vm), arb(0))
    quarter = arb(1) / 4
    a_Om = (acb(quarter, arb(Omega) / 2).digamma()).real - arb.pi().log()
    m_out = a_Om - 2 * S
    # sigma_max on [0, Omega]: |Re psi(1/4 + i r/2)| <= max(|psi(1/4)|, Re psi(1/4 + i Omega/2)) by monotonicity
    psi_q = abs(arb(quarter).digamma())
    sigma_max = max(psi_q, a_Om + arb.pi().log()) + arb.pi().log() + 2 * S
    # ---- chi_k = <p_k, cosh(t/2)> and ||chi||^2 = L + sinh L
    chi = []
    for k in ks:
        ik = (arb.pi() / (2 * (L / 2))).sqrt() * arb(L / 2).bessel_i(arb(k) + arb(1) / 2)   # modified spherical Bessel i_k(L/2)
        chi.append(((2 * k + 1) * L / 2).sqrt() * 2 * ik)
    chi_norm2 = L + L.sinh()
    Qchi2 = chi_norm2 - sum((c * c for c in chi), arb(0))
    # ---- tail bound eps_Q
    OL = Omega * L
    def ell_bound(k):
        return 2 / arb.pi() * OL ** (2 * k + 1) / arb(dfact_odd(2 * k + 1)) ** 2
    k1 = D + 2
    ratio = OL ** 2 / ((2 * k1 + 3) * (2 * k1 + 5))     # successive ratio bound (decreasing in k)
    eps_Q = ell_bound(k1) / (1 - ratio) if ratio < 1 else arb("1e300")   # D too small for the tail bound
    # ---- block B_PP
    cache = f"symbol_cache_l{lam}_O{Omega}_D{D}_d{dps}.json"
    done = json.load(open(cache)) if os.path.exists(cache) else {}
    tol = arb(10) ** (-told)
    eps0 = arb(10) ** (-(told + 10))
    half = arb(1) / 2
    lnn = [(Lam / arb(n).sqrt(), arb(n).log()) for n, Lam in vm]

    def sigma_c(z):   # analytic continuation of sigma(r) for complex r: Re psi -> (psi(1/4+iz/2)+psi(1/4-iz/2))/2
        iz = acb(0, 1) * z / 2
        ps = ((acb(quarter) + iz).digamma() + (acb(quarter) - iz).digamma()) / 2
        s = ps - arb.pi().log()
        for w, u in lnn:
            s -= 2 * w * (z * u).cos()
        return s

    def phat(k, z, sq):   # p_k^(z), sq = sqrt(pi / (2 z L)) precomputed
        pref = ((2 * k + 1) * L / 2).sqrt() * 2 * (1 if k % 4 == 0 else -1)
        return pref * sq * (z * L).bessel_j(acb(k) + half)

    pairs = [(a, b) for a in range(K) for b in range(a, K)]
    if pairs_arg:
        pairs = pairs[pairs_arg[0]:pairs_arg[1]]
    t0 = time.time(); n_new = 0
    for (a, b) in pairs:
        key = f"{a},{b}"
        if key in done:
            continue
        ja, jb = ks[a], ks[b]
        def integrand(z, analytic=False):
            if analytic and not (z.real > 0):
                return acb("nan")
            sq = (arb.pi() / (2 * z * L)).sqrt()
            return phat(ja, z, sq) * phat(jb, z, sq) * (sigma_c(z) - m_out)
        I = acb.integral(integrand, eps0, arb(Omega), rel_tol=tol, abs_tol=tol).real
        # [0, eps0] piece: |integrand| <= p_ja(0)^2-type bound: |p_k^(r)| <= sqrt((2k+1) L/2) * 2 * 1 (|j_k| <= 1), |sigma - m_out| <= sigma_max + m_out
        bnd = ((2 * ja + 1) * L / 2).sqrt() * ((2 * jb + 1) * L / 2).sqrt() * 4 * (sigma_max + m_out) * eps0
        I = I + arb(0).union(bnd).union(-bnd)
        val = I / arb.pi()
        done[key] = [val.mid().str(dps + 5, radius=False), val.rad().str(12, radius=False)]
        n_new += 1
        if n_new % 20 == 0:
            json.dump(done, open(cache, 'w'))
            print(f"  {len(done)}/{len(pairs) if not pairs_arg else K*(K+1)//2} entries, {time.time()-t0:.0f}s", flush=True)
    json.dump(done, open(cache, 'w'))
    if len(done) < K * (K + 1) // 2:
        print(f"partial: {len(done)} of {K*(K+1)//2} entries done ({time.time()-t0:.0f}s)"); return
    # assemble
    B = arb_mat(K, K)
    for a in range(K):
        for b in range(a, K):
            m, r = done[f"{a},{b}"]
            v = arb(m) + arb(0).union(2 * arb(r)).union(-2 * arb(r))
            if a == b:
                v = v + m_out
            v = v + 2 * chi[a] * chi[b]
            B[a, b] = v; B[b, a] = v
    # certify lambda_min(B_PP) by Gershgorin on Q^T B Q
    mp.mp.dps = dps
    Bm = mp.matrix([[mp.mpf(B[i, j].mid().str(dps, radius=False)) for j in range(K)] for i in range(K)])
    E, Qm = mp.eigsy(Bm)
    Qa = arb_mat([[arb(Qm[i, j].__str__()) for j in range(K)] for i in range(K)])
    Bt = Qa.transpose() * B * Qa
    gaps = []
    for i in range(K):
        off = arb(0)
        for j in range(K):
            if j != i:
                off += Bt[i, j].abs_upper()
        gaps.append(Bt[i, i].lower() - off.upper())
    lam_min_PP = min(gaps)        # lower bound for the smallest eigenvalue of Q^T B Q, hence of B_PP (Q invertible)
    # Q invertible and B_PP similar: use the Gershgorin bound divided by the Rayleigh factor: eig(B) >= lambda_min(Q^T B Q)/||Q||^2
    QtQ = Qa.transpose() * Qa
    normQ2 = max(sum(QtQ[i, j].abs_upper() for j in range(K)) for i in range(K))   # row-sum bound on ||Q^T Q||
    lam_min_B = lam_min_PP / normQ2
    # Schur budget
    cross = (sigma_max + m_out) * eps_Q.sqrt() + 2 * chi_norm2.sqrt() * Qchi2.abs_upper().sqrt()
    lam_min_QQ = m_out - (sigma_max + m_out) * eps_Q
    need = cross ** 2 / lam_min_QQ
    maxrad = max(B[i, j].rad() for i in range(K) for j in range(K))
    print(f"lambda = {lam}, Omega = {Omega}, D = {D} (K = {K} even Legendre modes), dps = {dps}, primes {[n for n,_ in vm]}")
    print(f"  m_out = {m_out.str(10)},  sigma_max = {sigma_max.str(8)},  a(Omega) = {a_Om.str(8)}")
    print(f"  eps_Q <= {eps_Q.str(5)},  ||Q chi||^2 = {Qchi2.str(5)},  max entry radius {maxrad.str(3)}")
    print(f"  approx eigenvalues of B_PP: {[mp.nstr(e, 5) for e in sorted(E)[:4]]}")
    print(f"  certified lambda_min(B_PP) >= {lam_min_B.str(8)}   (Gershgorin gap {lam_min_PP.str(6)}, ||Q||^2 <= {normQ2.str(6)})")
    print(f"  Schur requirement ||B_PQ||^2 / lambda_min(B_QQ) <= {need.str(5)}   (||B_PQ|| <= {cross.str(5)}, lambda_min(B_QQ) >= {lam_min_QQ.str(6)})")
    ok = (lam_min_B > need) and (lam_min_QQ > 0)
    print("CERTIFIED: W_lambda >= B >= 0 on all even f supported in [-L, L]" if ok else "NOT CERTIFIED")


if __name__ == "__main__":
    main()
