"""Compare the Weil minimiser f_lambda(t) with Riemann's Phi(t) on [-L, L]:
Xi(r) = xi(1/2+ir) = int_{-inf}^{inf} Phi(u) e^{iru} du,  Phi(u) = sum_n (2 pi^2 n^4 e^{9u/2} - 3 pi n^2 e^{5u/2}) e^{-pi n^2 e^{2u}}."""
import sys, json, mpmath as mp
from flint import ctx
import weil_prime as wp
mp.mp.dps=60
def Phi(u):
    return mp.fsum((2*mp.pi**2*n**4*mp.exp(9*u/2) - 3*mp.pi*n**2*mp.exp(5*u/2))*mp.exp(-mp.pi*n**2*mp.exp(2*u)) for n in range(1,40))
def xi(r):
    s=mp.mpc(0.5,r); return mp.re(mp.mpf(1)/2*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))
# sanity: Fourier transform of Phi vs xi
for r in [0, 5, 14.134725141734693]:
    I=2*mp.quad(lambda u: Phi(u)*mp.cos(r*u), [0, 1, 2, 4])
    print("r=%6.3f  int Phi cos = %s   xi(1/2+ir) = %s" % (r, mp.nstr(I,12), mp.nstr(xi(r),12)))
lam=int(sys.argv[1]); N=int(sys.argv[2]); dps=int(2.6*N)+80
ctx.dps=dps; mp.mp.dps=dps
W,vm=wp.build_W(N,lam,dps,8)
E,Q=mp.eigsy(W); i0=min(range(N+1),key=lambda i:E[i]); v=[Q[k,i0] for k in range(N+1)]
L=mp.log(lam)
mp.mp.dps=60
f=lambda t: mp.fsum(v[k]*mp.cos(k*t) for k in range(N+1))
# normalise so that f(0)/Phi(0) = 1
c=f(0)/Phi(0)
print("lambda=%d N=%d L=%.4f  normalisation c=%s" % (lam,N,float(L),mp.nstr(c,6)))
print("   t      f(t)/(c Phi(t))      Phi(t)        log10 Phi")
t=mp.mpf(0)
while t<=L+mp.mpf('1e-9'):
    ratio=f(t)/(c*Phi(t))
    print("%6.3f   %s   %s   %6.2f" % (float(t), mp.nstr(ratio,10), mp.nstr(Phi(t),6), float(mp.log10(abs(Phi(t))))))
    t+=L/12
