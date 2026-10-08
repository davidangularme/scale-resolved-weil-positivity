"""Prime-to-zero Jacobian at scale lambda: S_jn = d gamma_j / d eps_n where Lambda(n) -> (1+eps) Lambda(n).
Finite difference with eps = 1e-20 at full precision.  Prediction to test: S_jn ∝ (Lambda(n)/sqrt n) cos(gamma_j log n) x envelope."""
import sys, json, mpmath as mp
from flint import ctx
import weil_prime as wp
from rigidity import pert_matrix
lam=int(sys.argv[1]); N=int(sys.argv[2]); dps=int(2.6*N)+80
ctx.dps=dps; mp.mp.dps=dps
W,vm=wp.build_W(N,lam,dps,8)
L=mp.log(lam)
zeros=[mp.mpf(l) for l in open('zeros500_120d.txt')]
def minimiser(M):
    E,Q=mp.eigsy(M); i0=min(range(N+1),key=lambda i:E[i]); return [Q[k,i0] for k in range(N+1)]
def zeros_of(v, targets):
    Psi=lambda r: mp.fsum(v[k]*wp.ghat(k,r,L) for k in range(N+1))
    out=[]
    for g in targets:
        try: out.append(mp.findroot(Psi,(g-mp.mpf('0.3'),g+mp.mpf('0.3')),solver='anderson'))
        except Exception: out.append(mp.findroot(Psi,g))
    return out
Hh=2*mp.pi*lam**2
targets=[g for g in zeros if g< min(N, Hh)]
v0=minimiser(W); z0=zeros_of(v0,targets)
eps=mp.mpf(sys.argv[3]) if len(sys.argv)>3 else mp.mpf('1e-45')
print(json.dumps(dict(lam=lam,N=N,n_zeros=len(targets),primes=[n for n,_ in vm])))
for n,Lam in vm:
    T=pert_matrix(N,L,mp.log(n),Lam/mp.sqrt(n))
    v1=minimiser(W+eps*T); z1=zeros_of(v1,targets)
    S=[(z1[j]-z0[j])/eps for j in range(len(targets))]
    print(json.dumps(dict(n=n, w=float(Lam/mp.sqrt(n)), S=[float(s) for s in S], cosg=[float(mp.cos(targets[j]*mp.log(n))) for j in range(len(targets))], gam=[float(g) for g in targets])), flush=True)
