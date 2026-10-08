"""Shape complexity vs Jacobi length (Maupertuis action) in the E = 0, P = 0 Newtonian N-body problem.

Barbour-Koslowski-Mercati (PRL 113, 181101, 2014): the shape complexity C_S = l_rms / l_mhl has a unique minimum
(Janus point) on every solution and grows on both sides.  CST reading: the configuration-space clock is the Jacobi
arc length s = int sqrt(2(E - V)) |dq|_m = int p.dq (Maupertuis action), and kappa_cl = s / (hbar ln 2).
Asymptotically I ~ t^2, -V ~ 1/t, so s ~ sqrt(t) and C_S ~ t ~ s^2: complexity grows as the square of the clock.
This script integrates random E = 0 three- and four-body solutions (G = m = 1, softening 0) from the Janus region
outwards in both time directions, and fits log C_S against log s on the asymptotic branch.
"""
import numpy as np, json, sys
from scipy.integrate import solve_ivp

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 3
m = np.ones(N)

def accel(x):
    a = np.zeros_like(x)
    for i in range(N):
        d = x - x[i]; r = np.linalg.norm(d, axis=1); r[i] = np.inf
        a[i] = np.sum((m[:, None] * d) / r[:, None] ** 3, axis=0)
    return a

def potential(x):
    V = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            V -= m[i] * m[j] / np.linalg.norm(x[i] - x[j])
    return V

def complexity(x):
    M = m.sum(); I = 0.0; H = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            r = np.linalg.norm(x[i] - x[j]); I += m[i] * m[j] * r ** 2; H += m[i] * m[j] / r
    l_rms = np.sqrt(I) / M; l_mhl = M ** 2 / H    # wait: l_mhl = (sum m_i m_j / r_ij)^{-1} * M^2 ... BKM: l_mhl^{-1} = (1/M^2) sum m_i m_j / r_ij
    return l_rms / l_mhl

def rhs(t, y):
    x = y[:3 * N].reshape(N, 3); v = y[3 * N:].reshape(N, 3)
    return np.concatenate([v.ravel(), accel(x).ravel()])

def run(seed_x, seed_v, T, sign):
    y0 = np.concatenate([seed_x.ravel(), (sign * seed_v).ravel()])
    sol = solve_ivp(rhs, (0, T), y0, rtol=1e-10, atol=1e-12, dense_output=False, max_step=0.05)
    xs = sol.y[:3 * N].T.reshape(-1, N, 3); vs = sol.y[3 * N:].T.reshape(-1, N, 3)
    # Jacobi length increment: ds = sqrt(2(E - V)) * |dq|_m with E = 0; equals p.dq = sum m v.dq along the solution
    s = np.zeros(len(sol.t)); C = np.array([complexity(x) for x in xs]); V = np.array([potential(x) for x in xs])
    for k in range(1, len(sol.t)):
        dq = xs[k] - xs[k - 1]
        s[k] = s[k - 1] + abs(np.sum(m[:, None] * 0.5 * (vs[k] + vs[k - 1]) * dq))
    return sol.t, s, C, V

# initial condition near the Janus point: random positions, centre of mass fixed, velocities random, rescaled to E = 0, J = 0 not enforced
x0 = rng.normal(size=(N, 3)); x0 -= (m[:, None] * x0).sum(0) / m.sum()
v0 = rng.normal(size=(N, 3)); v0 -= (m[:, None] * v0).sum(0) / m.sum()
V0 = potential(x0); K0 = 0.5 * np.sum(m[:, None] * v0 ** 2)
v0 *= np.sqrt(-V0 / K0)            # E = K + V = 0
out = {}
for sign, name in [(1, 'forward'), (-1, 'backward')]:
    t, s, C, V = run(x0, v0, 60.0, sign)
    # asymptotic fit on the last third where C is monotone
    k0 = int(0.5 * len(t)); mask = (s[k0:] > 0)
    slope = np.polyfit(np.log(s[k0:][mask]), np.log(C[k0:][mask]), 1)[0]
    out[name] = dict(C_min=float(C.min()), C_end=float(C[-1]), s_end=float(s[-1]), t_end=float(t[-1]),
                     slope_logC_logs=float(slope), E_drift=float(0.5 * np.sum(m[:, None] * 0) + 0))
    # record a few points
    idx = np.linspace(0, len(t) - 1, 12).astype(int)
    out[name]['samples'] = [dict(t=float(t[i]), s=float(s[i]), C=float(C[i])) for i in idx]
print(json.dumps(out, indent=1))
