"""Rigidity of the prime weights under Weil positivity at scale lambda.

W(eps) = W_lambda + eps * T_n,   T_n = -2 Lambda(n)/sqrt(n) * [h_jk(log n)]   (the n-th prime term of W)
i.e. the weight Lambda(n) is replaced by (1+eps) Lambda(n).  Galerkin indefiniteness of W(eps) is a NECESSARY
condition for true indefiniteness to fail... precisely: if the Galerkin matrix is indefinite then the true form
is indefinite.  So the Galerkin tolerance interval CONTAINS the true one: it is an upper bound on how far the
weight can be moved.  Also: shift the position log n -> log n + eta.
"""
import sys, json, time
import mpmath as mp
from flint import ctx
import weil_prime as wp


def min_eig(A):
    return min(mp.eigsy(A, eigvals_only=True))


def pert_matrix(N, L, u, weight):
    H = mp.matrix(N + 1, N + 1)
    for j in range(N + 1):
        for k in range(j, N + 1):
            v = -2 * weight * wp.h_jk(j, k, u, L)
            H[j, k] = v; H[k, j] = v
    return H


def bisect_boundary(W, T, lo, hi, tol):
    """largest eps in [lo,hi] (lo feasible) with W + eps*T positive definite; assumes monotone."""
    assert min_eig(W + lo * T) > 0
    if min_eig(W + hi * T) > 0:
        return hi, True
    while hi - lo > tol * max(1, abs(hi)):
        mid = (lo + hi) / 2
        if min_eig(W + mid * T) > 0:
            lo = mid
        else:
            hi = mid
    return lo, False


if __name__ == "__main__":
    lam = int(sys.argv[1]); N = int(sys.argv[2]); dps = int(2.6 * N) + 80
    ctx.dps = dps; mp.mp.dps = dps
    t0 = time.time()
    W, vm = wp.build_W(N, lam, dps, 8)
    L = mp.log(lam)
    print(json.dumps(dict(lam=lam, N=N, build_s=round(time.time() - t0, 1), mu0=float(min_eig(W)))), flush=True)
    for n, Lam in vm:
        T = pert_matrix(N, L, mp.log(n), Lam / mp.sqrt(n))
        up, upfree = bisect_boundary(W, T, mp.mpf(0), mp.mpf(10), mp.mpf('1e-12'))
        dn, dnfree = bisect_boundary(W, -T, mp.mpf(0), mp.mpf(10), mp.mpf('1e-12'))
        print(json.dumps(dict(n=n, weight=float(Lam / mp.sqrt(n)), eps_plus=float(up), eps_minus=float(-dn),
                              unbounded_plus=upfree, unbounded_minus=dnfree)), flush=True)
    # position shift of the prime 2 and of the largest prime power present
    for n in [vm[0][0], vm[-1][0]]:
        Lam = dict(vm)[n]
        T0 = pert_matrix(N, L, mp.log(n), Lam / mp.sqrt(n))
        def ok(eta):
            Tn = pert_matrix(N, L, mp.log(n) + eta, Lam / mp.sqrt(n))
            return min_eig(W - T0 + Tn) > 0
        res = {}
        for sgn, name in [(1, 'eta_plus'), (-1, 'eta_minus')]:
            lo, hi = mp.mpf(0), mp.mpf('0.5')
            if ok(sgn * hi):
                res[name] = float(sgn * hi); res[name + '_unbounded'] = True; continue
            while hi - lo > mp.mpf('1e-12'):
                mid = (lo + hi) / 2
                if ok(sgn * mid): lo = mid
                else: hi = mid
            res[name] = float(sgn * lo); res[name + '_unbounded'] = False
        print(json.dumps(dict(n=n, shift=res)), flush=True)
