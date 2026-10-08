"""Landau-Widom check: number of prolate eigenvalues above eps vs 2c/pi + (1/pi^2) ln c ln((1-eps)/eps)."""
import sys, json
import mpmath as mp
import prolate
mp.mp.dps = 30
out = []
for c in [10, 30, 100, 300]:
    M = int(2 * c / mp.pi) // 2 + 40
    ev = prolate.prolate_eigs(c, 0, M); od = prolate.prolate_eigs(c, 1, M)
    allv = sorted(ev + od, reverse=True)
    row = dict(c=c, landau=float(2 * c / mp.pi), total=float(mp.fsum(allv)), counts={})
    for e in ['0.5', '1e-2', '1e-4', '1e-6', '1e-10']:
        eps = mp.mpf(e)
        n = sum(1 for v in allv if v > eps)
        pred = 2 * c / mp.pi + mp.log(c) * mp.log((1 - eps) / eps) / mp.pi ** 2
        row['counts'][e] = dict(count=n, lw=float(pred))
    out.append(row); print(json.dumps(row), flush=True)
json.dump(out, open('lw_check.json', 'w'), indent=1)
