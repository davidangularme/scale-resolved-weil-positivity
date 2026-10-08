"""Window w_lambda(t) = f_lambda(t) / (c Phi(t)) on [0, L], normalised w(0)=1, sampled at s = t/L.
Then fit against (a) Gaussian exp(-b s^2), (b) bump exp(-b s^2/(1-s^2)), (c) prolate ground function psi_0(c, s)."""
import sys, json, mpmath as mp, numpy as np
from flint import ctx
import weil_prime as wp
from scipy.special import pro_ang1
from scipy.optimize import minimize_scalar
def Phi(u):
    return mp.fsum((2*mp.pi**2*n**4*mp.exp(9*u/2) - 3*mp.pi*n**2*mp.exp(5*u/2))*mp.exp(-mp.pi*n**2*mp.exp(2*u)) for n in range(1,40))
lam=float(sys.argv[1]); N=int(sys.argv[2]); dps=int(2.6*N)+80
ctx.dps=dps; mp.mp.dps=dps
W,vm=wp.build_W(N,lam,dps,8)
E,Q=mp.eigsy(W); i0=min(range(N+1),key=lambda i:E[i]); v=[Q[k,i0] for k in range(N+1)]
L=mp.log(lam)
f=lambda t: mp.fsum(v[k]*mp.cos(k*t) for k in range(N+1))
c0=f(0)/Phi(0)
S=np.linspace(0,0.98,50)
w=np.array([float(f(s*L)/(c0*Phi(s*L))) for s in S])
def fit(model):
    m=(S<=0.9)&(w>1e-12)
    r=minimize_scalar(lambda b: np.sum((np.log(np.abs(w[m]))-np.log(np.abs(model(b,S[m]))))**2), bounds=(0.05,60), method='bounded')
    return r.x, r.fun
gauss=lambda b,s: np.exp(-b*s**2)
bump=lambda b,s: np.exp(-b*s**2/(1-s**2))
prol=lambda c,s: pro_ang1(0,0,c,s)[0]/pro_ang1(0,0,c,0.0)[0]
res={}
for name,m in [('gauss',gauss),('bump',bump),('prolate',prol)]:
    p,err=fit(m); res[name]=dict(param=float(p), rms_log=float(np.sqrt(err/len(S))))
print(json.dumps(dict(lam=lam,N=N,L=float(L),fits=res)))
for s,ww in zip(S,w): print(json.dumps(dict(s=round(float(s),3), w=float(ww))))
