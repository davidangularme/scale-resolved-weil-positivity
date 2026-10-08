"""Reproduce and extend the 'Weil minimiser' spectral table of the CST/FDBC archive.

M = 2 * Phi Phi^T with
  Phi_0(g) = 2 sin(g L)/g
  Phi_k(g) = sin((k-g)L)/(k-g) + sin((k+g)L)/(k+g)   (k >= 1)
summed over a point set {g} (Riemann zero ordinates, or control point sets).
Prints log10 mu0, log10 mu1, log10(mu1/mu0), alpha(N)=log10(mu1/mu0)/N.
"""
import sys, time, json
import flint
from flint import arb, arb_mat, acb_mat, ctx


def load_points(path, n):
    with open(path) as f:
        return [l.strip() for l in f if l.strip()][:n]


def build_M(points, N, L):
    rows = []
    for k in range(N + 1):
        row = []
        for g in points:
            g = arb(g)
            if k == 0:
                row.append(2 * (g * L).sin() / g)
            else:
                a = (k - g)
                b = (k + g)
                row.append((a * L).sin() / a + (b * L).sin() / b)
        rows.append(row)
    Phi = arb_mat(rows)
    return 2 * Phi * Phi.transpose()


def two_smallest(M):
    # eigenvalues of a real symmetric matrix via acb_mat.eig (approximate, then sorted)
    A = acb_mat(M)
    ev = A.eig(algorithm="approx")
    vals = sorted([e.real.mid() for e in ev], key=lambda x: x)
    return vals[0], vals[1]


def main():
    path = sys.argv[1]
    lam = float(sys.argv[2])
    Ns = [int(x) for x in sys.argv[3].split(",")]
    npts = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    pts = load_points(path, npts)
    out = []
    for N in Ns:
        ctx.dps = int(2.6 * N) + 80
        L = arb(lam).log()
        t = time.time()
        M = build_M(pts, N, L)
        m0, m1 = two_smallest(M)
        l0 = float(abs(m0).log() / arb(10).log()) if m0 != 0 else float("-inf")
        l1 = float(abs(m1).log() / arb(10).log())
        r = l1 - l0
        rec = dict(N=N, log10_mu0=round(l0, 2), log10_mu1=round(l1, 2),
                   log10_ratio=round(r, 2), alpha=round(r / N, 4),
                   sign_mu0=("+" if m0 > 0 else "-"), secs=round(time.time() - t, 1))
        print(json.dumps(rec), flush=True)
        out.append(rec)


if __name__ == "__main__":
    main()
