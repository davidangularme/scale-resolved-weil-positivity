"""Detection experiment: does Weil positivity on [1/lambda, lambda] see one hypothetical off-line zero?

Zero-side quadratic form in the cosine basis g_k(t)=cos(kt), t in [-L,L], L=log(lambda):
  Q[a] = sum over zeros rho of F(rho) conj F(1-conj(rho)),   F = sum_k a_k ghat_k.
On-line pair (+-gamma):      2 * ghat_j(gamma) ghat_k(gamma)
Off-line quadruplet at height gamma, distance delta from the line:
      4 * Re[ ghat_j(r) ghat_k(r) ],  r = gamma - i*delta
      = 4 (Re ghat_j Re ghat_k - Im ghat_j Im ghat_k).
We take the first NZ true zeros, move the one nearest `target` off the line by delta,
and report the smallest eigenvalue of the (N+1)x(N+1) matrix.  Negative => positivity fails.
"""
import sys, json
from flint import arb, acb, arb_mat, acb_mat, ctx


def ghat(k, r, L):
    if k == 0:
        return 2 * (r * L).sin() / r
    a, b = k - r, k + r
    return (a * L).sin() / a + (b * L).sin() / b


def base_matrix(zeros, N, L):
    rows = [[ghat(k, g, L) for g in zeros] for k in range(N + 1)]
    Phi = arb_mat(rows)
    return 2 * Phi * Phi.transpose(), Phi


def min_eig(M):
    ev = acb_mat(M).eig(algorithm="approx")
    return min(e.real.mid() for e in ev)


def run(lam, N, NZ, target, deltas, dps):
    ctx.dps = dps
    zeros = [arb(l.strip()) for l in open('zeros3000.txt')][:NZ]
    L = arb(lam).log()
    M0, Phi = base_matrix(zeros, N, L)
    m = min(range(NZ), key=lambda i: abs(zeros[i] - target))
    gam = zeros[m]
    col = arb_mat([[Phi[k, m]] for k in range(N + 1)])
    Mbase = M0 - 2 * col * col.transpose()
    out = dict(lam=lam, L=float(L), N=N, NZ=NZ, gamma=float(gam), mu0_online=float(min_eig(M0)), rows=[])
    for d in deltas:
        r = acb(gam, -arb(d))
        psi = [ghat(k, r, L) for k in range(N + 1)]
        Re = arb_mat([[p.real] for p in psi]); Im = arb_mat([[p.imag] for p in psi])
        M = Mbase + 4 * (Re * Re.transpose() - Im * Im.transpose())
        mu = min_eig(M)
        out['rows'].append((d, float(mu)))
        print(json.dumps(dict(lam=lam, N=N, gamma=round(float(gam), 3), delta=d, min_eig=float(mu))), flush=True)
    return out


if __name__ == "__main__":
    lam = float(sys.argv[1]); N = int(sys.argv[2]); NZ = int(sys.argv[3]); target = float(sys.argv[4])
    deltas = [float(x) for x in sys.argv[5].split(',')]
    dps = int(sys.argv[6]) if len(sys.argv) > 6 else 120
    run(lam, N, NZ, target, deltas, dps)
