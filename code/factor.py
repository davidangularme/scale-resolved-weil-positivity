"""Test the factorisation Psi_lambda(r) = Xi(r) * G_lambda(r): is G smooth and zero-free below the horizon?
Xi(r) = xi(1/2 + i r), xi(s) = (1/2) s (s-1) pi^{-s/2} Gamma(s/2) zeta(s)   (real on the line)."""
import sys, json, mpmath as mp
from flint import ctx
import weil_prime as wp
lam=int(sys.argv[1]); N=int(sys.argv[2]); dps=int(2.6*N)+80
ctx.dps=dps; mp.mp.dps=dps
W,vm=wp.build_W(N,lam,dps,8)
E,Q=mp.eigsy(W); i0=min(range(N+1),key=lambda i:E[i]); v=[Q[k,i0] for k in range(N+1)]
L=mp.log(lam)
Psi=lambda r: mp.fsum(v[k]*wp.ghat(k,r,L) for k in range(N+1))
mp.mp.dps=60
def Xi(r):
    s=mp.mpc(0.5,r)
    return mp.re(mp.mpf(1)/2*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))
H=2*mp.pi*lam**2
out=[]
r=mp.mpf('0.3123')
while r<=1.6*H:
    x=Xi(r); p=Psi(r)
    if abs(x)>mp.mpf('1e-40'):
        out.append((float(r), float(mp.log10(abs(p))), float(mp.log10(abs(x))), float(mp.log10(abs(p/x))), int(mp.sign(p/x))))
    r+=mp.mpf('0.75')
print(json.dumps(dict(lam=lam,N=N,H=float(H),mu0=float(E[i0]))))
for o in out: print(json.dumps(o))
