"""Effective prime mass: top eigenvalue (L2-normalised) of the compressed prime shift operator on [-L, L]."""
import mpmath as mp, math
import weil_prime as wp
from margin import gram
def vm(N):
    out=[]
    for n in range(2,N+1):
        m=n;p=None
        for q in range(2,n+1):
            if m%q==0: p=q;break
        while m%p==0: m//=p
        if m==1: out.append((n,mp.log(p)))
    return out
for lam, N in [(2,30),(3,30),(4,40),(5,40),(10,50)]:
    mp.mp.dps=int(2.6*N)+80
    L=mp.log(lam); G=gram(N,L)
    P=mp.matrix(N+1,N+1); S=mp.mpf(0)
    for n,Lam in vm(int(lam**2+1e-9)):
        if mp.log(n) >= 2*L: continue
        w=Lam/mp.sqrt(n); S+=w
        for j in range(N+1):
            for k in range(j,N+1):
                v=w*wp.h_jk(j,k,mp.log(n),L)
                P[j,k]+=v
                if j!=k: P[k,j]+=v
    C=mp.cholesky(G); Ci=mp.inverse(C); A=Ci*P*Ci.T
    ev=mp.eigsy(A,eigvals_only=True); top=max(ev)
    print(f"lam={lam} N={N} S={float(S):.3f} S_eff={float(top):.3f} K_new~{4*float(L)*math.exp(2*float(top)):.3g}", flush=True)
