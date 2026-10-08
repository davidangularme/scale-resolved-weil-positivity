"""Reconstruction test: how far above the resolution height H = 2 pi lambda^2 do the zeros of the
prime-side Weil minimiser's transform keep tracking the true zeta zeros?

The maximum-entropy (Krein/Burg) extension of the prime kernel has spectral density ~ 1/|Psi(r)|^2 with
Psi the transform of the minimiser, so its atoms are the real zeros of Psi.  We locate all real zeros of
Psi on (0, N] (N = basis frequency limit), match each to the nearest true zero, and report digits of agreement,
spurious zeros (no true zero within 0.5) and missed true zeros.
"""
import sys, json, time
import mpmath as mp
from flint import ctx
import weil_prime as wp

lam = int(sys.argv[1]); N = int(sys.argv[2]); dps = int(2.6 * N) + 80
ctx.dps = dps; mp.mp.dps = dps
t0 = time.time()
W, vm = wp.build_W(N, lam, dps, 8)
tb = time.time() - t0
E, Q = mp.eigsy(W)
i0 = min(range(N + 1), key=lambda i: E[i])
v = [Q[k, i0] for k in range(N + 1)]
L = mp.log(lam)
Psi = lambda r: mp.fsum(v[k] * wp.ghat(k, r, L) for k in range(N + 1))
zeros = [mp.mpf(l) for l in open('zeros500_120d.txt')]
Rmax = N
step = mp.mpf('0.05'); off = mp.mpf('0.0123456789')
grid = [off + step * i for i in range(1, int(Rmax / step))]
vals = [Psi(r) for r in grid]
found = []
for i in range(len(grid) - 1):
    if vals[i] * vals[i + 1] < 0:
        found.append(mp.findroot(Psi, (grid[i], grid[i + 1]), solver='anderson'))
true_below = [g for g in zeros if g < Rmax]
rows = []; used = set()
for z in found:
    j = min(range(len(zeros)), key=lambda t: abs(zeros[t] - z))
    d = abs(zeros[j] - z)
    rows.append(dict(r=float(z), nearest=float(zeros[j]), k=j + 1, diff=float(d), digits=float(-mp.log10(d)) if d > 0 else 99.0, spurious=bool(d > 0.5)))
    if d <= 0.5: used.add(j)
missed = [float(g) for t, g in enumerate(true_below) if t not in used]
H = float(2 * mp.pi * lam ** 2)
print(json.dumps(dict(lam=lam, N=N, H=H, build_s=round(tb, 1), n_true_below=len(true_below), n_found=len(found),
                      n_matched=len(used), n_spurious=sum(r['spurious'] for r in rows), missed=missed[:20])))
for r in rows:
    print(json.dumps(r))
