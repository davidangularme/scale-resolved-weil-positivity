import sys, time, json, mpmath as mp
from flint import ctx
import weil_prime as wp
N=int(sys.argv[1]); lam=int(sys.argv[2]); dps=int(2.6*N)+80
ctx.dps=dps
t=time.time()
W,vm=wp.build_W(N,lam,dps,8)
tb=time.time()-t
E,Q=mp.eigsy(W)
order=sorted(range(N+1),key=lambda i:E[i])
mu0,mu1=E[order[0]],E[order[1]]
v=[Q[k,order[0]] for k in range(N+1)]
L=mp.log(lam)
G=lambda r: mp.fsum(v[k]*wp.ghat(k,r,L) for k in range(N+1))
zeros=[mp.mpf(l) for l in open('zeros500.txt')]
# locate real zeros of ghat_v on (0, 70]
grid=[mp.mpf(i)/10+mp.mpf("0.0123456789") for i in range(1,701)]
vals=[G(r) for r in grid]
found=[]
for i in range(len(grid)-1):
    if vals[i]==0 or vals[i]*vals[i+1]<0:
        found.append(mp.findroot(G,(grid[i],grid[i+1]),solver='anderson'))
print(json.dumps(dict(N=N,lam=lam,build_s=round(tb,1),log10_mu0=float(mp.log10(abs(mu0))),log10_mu1=float(mp.log10(abs(mu1))),
      log10_ratio=float(mp.log10(mu1/mu0)),sign_mu0='+' if mu0>0 else '-')))
for z in found:
    near=min(zeros,key=lambda g:abs(g-z))
    print(f"  zero of minimiser at r = {mp.nstr(z,14):>18}   nearest zeta zero = {mp.nstr(near,14):>18}   |diff| = {mp.nstr(abs(z-near),3)}")
