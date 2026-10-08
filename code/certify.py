"""Rigorous (ball-arithmetic) certificate: the prime-side Weil form W_lambda is positive definite on the
span of cos(k t), k = 0..N, on [-L, L].  Every matrix entry is an Arb ball enclosing the true value
(closed forms for the pole and prime terms, rigorous integration for the archimedean term), and positivity
is certified by Gershgorin on Q^T W Q with Q an approximate eigenbasis.

This certifies the FINITE BLOCK only; the tail (orthogonal complement) is the remaining step of the theorem.
"""
import sys, time
from flint import arb, acb, arb_mat, ctx
import mpmath as mp


def h_jk_arb(j, k, x, L):
    """(g_j * g_k)(x) for 0 <= x <= 2L, g_j = cos(j t) on [-L, L], as a closed form (ball arithmetic)."""
    lo, hi = x - L, L
    def F(a, b):   # int_{lo}^{hi} cos(a t + b) dt
        if a == 0:
            return (hi - lo) * b.cos()
        return ((a * hi + b).sin() - (a * lo + b).sin()) / a
    return (F(arb(j - k), k * x) + F(arb(j + k), -k * x)) / 2


def build_W_rigorous(N, lam, tol):
    L = arb(lam).log()
    half = arb(1) / 2
    a = [(half.sinh() * (L / 2).cos() if False else None) for _ in range(N + 1)]
    for j in range(N + 1):
        a[j] = ((L / 2).sinh() * (j * L).cos() + 2 * j * (L / 2).cosh() * (j * L).sin()) / (arb(1) / 4 + j * j)
    # prime powers n <= lam^2 with log n < 2L
    vm = []
    for n in range(2, int(lam * lam) + 1):
        m, p = n, None
        for q in range(2, n + 1):
            if m % q == 0:
                p = q
                break
        while m % p == 0:
            m //= p
        if m == 1 and arb(n).log() < 2 * L:
            vm.append((n, arb(p).log()))
    eps = arb(10) ** (-120)
    T = 4 * L
    tail = -(1 - (-T).exp()).log()
    W = arb_mat(N + 1, N + 1)
    for j in range(N + 1):
        for k in range(j, N + 1):
            h0 = h_jk_arb(j, k, arb(0), L)
            def integrand(t, analytic=False):
                tt = t
                return ((-tt).exp() * h0 - (-tt / 4).exp() * h_jk_acb(j, k, tt / 2, L)) / (1 - (-tt).exp())
            I = acb.integral(integrand, eps, T, rel_tol=tol, abs_tol=tol).real
            # |integral over [0, eps]| <= 2 * M * eps where M bounds the derivative of the numerator on [0, eps]
            xb = arb(0).union(eps / 2)
            dh = dh_jk_arb(j, k, xb, L)
            M = abs(h0) * (1 + arb(1) / 4) + abs(dh) / 2
            I = I + arb(0).union(2 * M * eps).union(-2 * M * eps)
            arch = -arb.const_euler() * h0 + tail * h0 + I
            primes = arb(0)
            for n, Lam in vm:
                primes += Lam / arb(n).sqrt() * h_jk_arb(j, k, arb(n).log(), L)
            val = 2 * a[j] * a[k] - arb.pi().log() * h0 + arch - 2 * primes
            W[j, k] = val
            W[k, j] = val
    return W, vm


def h_jk_acb(j, k, x, L):
    lo, hi = x - L, L
    def F(a, b):
        if a == 0:
            return (hi - lo) * b.cos()
        return ((a * hi + b).sin() - (a * lo + b).sin()) / a
    return (F(acb(j - k), k * x) + F(acb(j + k), -k * x)) / 2


def dh_jk_arb(j, k, x, L):
    """d/dx (g_j * g_k)(x) enclosure: h'(x) = -g_j(x - L) g_k(L) ... use the closed form derivative."""
    lo, hi = x - L, L
    def dF(a, b):   # d/dx of int_{x-L}^{L} cos(a t + b(x)) dt with b = +-k x  (b' = +-k)
        # F = (sin(a hi + b) - sin(a lo + b)) / a ; dF/dx = (cos(a hi + b) b' - cos(a lo + b)(a + b')) / a
        return None
    # simpler: bound |h'| by the max of |d/dx| via direct formula h'(x) = -(cos(j(x-L)) cos(k L)) + int ... ; use a crude bound:
    # |h'(x)| <= 2L * (j + k) + 1  (derivative of a product of unit-amplitude cosines integrated over length <= 2L, plus boundary term)
    return arb(2 * L * (j + k) + 1)


def certify(W, N, dps):
    mp.mp.dps = dps
    Wm = mp.matrix([[mp.mpf(W[i, j].mid().str(dps, radius=False)) for j in range(N + 1)] for i in range(N + 1)])
    E, Q = mp.eigsy(Wm)
    Qa = arb_mat([[arb(Q[i, j].__str__()) for j in range(N + 1)] for i in range(N + 1)])
    B = Qa.transpose() * W * Qa
    worst = None
    for i in range(N + 1):
        lower = B[i, i].lower()
        off = arb(0)
        for j in range(N + 1):
            if j != i:
                off += B[i, j].abs_upper()
        gap = lower - off.upper()
        if worst is None or gap < worst[0]:
            worst = (gap, i)
    detQ = Qa.det()
    return worst, detQ, min(E)


if __name__ == "__main__":
    lam = int(sys.argv[1]); N = int(sys.argv[2]); dps = int(sys.argv[3]) if len(sys.argv) > 3 else int(2.6 * N) + 80
    ctx.dps = dps
    t0 = time.time()
    W, vm = build_W_rigorous(N, lam, arb(10) ** (-(dps - 20)))
    tb = time.time() - t0
    maxrad = max(W[i, j].rad() for i in range(N + 1) for j in range(N + 1))
    worst, detQ, mu0 = certify(W, N, dps)
    print("lambda =", lam, " N =", N, " dps =", dps, " primes:", [n for n, _ in vm], " build %.0fs" % tb)
    print("max entry radius:", maxrad.str(3))
    print("approx smallest eigenvalue (coefficient space):", mp.nstr(mu0, 6))
    print("det Q (ball):", detQ.str(5))
    print("Gershgorin worst gap (lower bound of diag minus off-diagonals):", worst[0].str(6), " at row", worst[1])
    print("CERTIFIED POSITIVE DEFINITE" if worst[0] > 0 and not detQ.contains(0) else "NOT CERTIFIED")
