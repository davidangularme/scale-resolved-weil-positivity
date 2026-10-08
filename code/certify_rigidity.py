"""Certified rigidity of the prime data under Weil positivity at scale lambda (ball arithmetic).

For a perturbed form W' (weight Lambda(n) -> (1+eps) Lambda(n), or position log n -> log n + eta) we certify
INDEFINITENESS by exhibiting an explicit even trigonometric polynomial f0 (coefficients in the cosine basis) with
W'(f0) < 0, every quantity enclosed in an Arb ball.  Since the cosine polynomials are in L^2[-L, L], this proves that
the perturbed Weil form is indefinite on the whole space; the unperturbed form is positive there (Zhu 2026, Thm 1.2,
for L <= 0.8; independently certified for L = log 2 in certify_symbol.py).  Candidate boundaries come from a bisection
on the floating-point Galerkin matrix; the certified statement uses the boundary enlarged by 5 %.
Usage: certify_rigidity.py lam N [dps]
"""
import sys, time, json
from flint import arb, arb_mat, ctx
import mpmath as mp
import certify as C


def ball_H(N, L, x, weight):
    """prime-term matrix T = -2 weight [h_jk(x)] as balls"""
    T = arb_mat(N + 1, N + 1)
    for j in range(N + 1):
        for k in range(j, N + 1):
            v = -2 * weight * C.h_jk_arb(j, k, x, L)
            T[j, k] = v; T[k, j] = v
    return T


def mid(M, N, dps):
    return mp.matrix([[mp.mpf(M[i, j].mid().str(dps, radius=False)) for j in range(N + 1)] for i in range(N + 1)])


def min_vec(A):
    E, Q = mp.eigsy(A)
    i = min(range(len(E)), key=lambda t: E[t])
    return E[i], [Q[k, i] for k in range(len(E))]


def quad_ball(M, a, N):
    av = arb_mat(N + 1, 1)
    for k in range(N + 1):
        av[k, 0] = arb(a[k].__str__())
    return (av.transpose() * M * av)[0, 0]


def bisect(Wm, Tm, lo, hi, tol):
    """largest s in [lo, hi] with Wm + s Tm positive definite (lo feasible), bisection on mid matrices"""
    f = lambda s: min(mp.eigsy(Wm + s * Tm, eigvals_only=True)) > 0
    lo = mp.mpf('1e-40')
    assert f(lo)
    if f(hi):
        return hi, True
    while hi / lo > 1 + tol:          # geometric bisection: the boundary may be tiny
        m = mp.sqrt(lo * hi)
        if f(m): lo = m
        else: hi = m
    return lo, False


if __name__ == "__main__":
    lam = int(sys.argv[1]); N = int(sys.argv[2]); dps = int(sys.argv[3]) if len(sys.argv) > 3 else int(2.6 * N) + 80
    ctx.dps = dps; mp.mp.dps = dps
    t0 = time.time()
    W, vm = C.build_W_rigorous(N, lam, arb(10) ** (-(dps - 20)))
    L = arb(lam).log()
    print(f"lambda = {lam}, N = {N}, dps = {dps}, primes {[n for n,_ in vm]}, build {time.time()-t0:.0f}s, max radius {max(W[i,j].rad() for i in range(N+1) for j in range(N+1)).str(3)}", flush=True)
    Wm = mid(W, N, dps)
    mu0, a0 = min_vec(Wm)
    print(f"  Galerkin min eigenvalue (coefficient space) {mp.nstr(mu0, 6)}; W(f0) ball = {quad_ball(W, a0, N).str(6)}")
    out = {}
    for n, Lam in vm:
        w = Lam / arb(n).sqrt()
        T = ball_H(N, L, arb(n).log(), w)
        Tm = mid(T, N, dps)
        res = {}
        for sgn, name in [(1, 'eps_plus'), (-1, 'eps_minus')]:
            b, free = bisect(Wm, sgn * Tm, mp.mpf(0), mp.mpf(10), mp.mpf('1e-3'))
            if free:
                res[name] = None; continue
            s = b * mp.mpf('1.05')
            _, a = min_vec(Wm + sgn * s * Tm)
            val = quad_ball(W + arb((sgn * s).__str__()) * T, a, N)
            res[name] = dict(boundary=float(b), certified_eps=float(sgn * s), form_value=val.str(6), certified=bool(val.upper() < 0))
            print(f"  n = {n}: {name}: boundary {mp.nstr(b, 6)}, at eps = {mp.nstr(sgn*s, 6)} the form takes value {val.str(6)} on f0 -> {'INDEFINITE CERTIFIED' if val.upper() < 0 else 'not certified'}", flush=True)
        out[n] = res
    # position shift of each prime power: log n -> log n + eta
    for n, Lam in vm:
        w = Lam / arb(n).sqrt()
        T0 = ball_H(N, L, arb(n).log(), w); T0m = mid(T0, N, dps)
        def Weta_m(eta):
            return Wm - T0m + mid(ball_H(N, L, arb(n).log() + arb(eta.__str__()), w), N, dps)
        res = {}
        for sgn, name in [(1, 'eta_plus'), (-1, 'eta_minus')]:
            lo, hi = mp.mpf('1e-40'), mp.mpf('0.3')
            ok = lambda e: min(mp.eigsy(Weta_m(sgn * e), eigvals_only=True)) > 0
            if ok(hi):
                res[name] = None; continue
            while hi / lo > mp.mpf('1.001'):
                m = mp.sqrt(lo * hi)
                if ok(m): lo = m
                else: hi = m
            e = lo * mp.mpf('1.05')
            _, a = min_vec(Weta_m(sgn * e))
            Wb = W - T0 + ball_H(N, L, arb(n).log() + arb((sgn * e).__str__()), w)
            val = quad_ball(Wb, a, N)
            res[name] = dict(boundary=float(sgn * lo), certified_eta=float(sgn * e), form_value=val.str(6), certified=bool(val.upper() < 0))
            print(f"  n = {n}: {name}: boundary {mp.nstr(sgn*lo, 6)}, at eta = {mp.nstr(sgn*e, 6)} the form takes value {val.str(6)} on f0 -> {'INDEFINITE CERTIFIED' if val.upper() < 0 else 'not certified'}", flush=True)
        out[f"shift_{n}"] = res
    json.dump(out, open(f"rigidity_cert_l{lam}_n{N}.json", "w"), indent=1)
    print(f"total {time.time()-t0:.0f}s")
