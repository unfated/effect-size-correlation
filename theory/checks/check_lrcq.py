"""LRCQ checks (supplement §S1.6, §S3).

1. fast OLS  w_hat = D^{-1} (Z∘Z - c)' s / ||s||^2  equals brute-force stacked OLS.
2. fast WLS  (sum_a s_a^2 D W_a D)^{-1} sum_a s_a D W_a (z_a^2 - c_a) equals brute-force WLS.
3. With local genetic-effect correlations, E[w_hat_OLS] = w + D^{-1} diag(R C R)
   (C = off-diagonal part of Sigma~); without them it is unbiased.
4. With trait-specific enrichments w_a, E[w_hat_OLS] = sum_a s_a^2 w_a / ||s||^2.
5. Free per-trait intercepts make w non-identifiable: the stacked design with
   intercept columns is rank-deficient by exactly one (direction D^{-1} 1).
"""
import numpy as np
from common import ar1_ld, mvn

rng = np.random.default_rng(7)
p, q = 40, 60
R = ar1_ld(p, 0.7)
D = R ** 2
n = rng.integers(50_000, 400_000, q)
h2 = rng.uniform(0.05, 0.5, q)
m = 1_000
s = n * h2 / m

# true enrichment: sparse, mean ~1 over the window not required
w = np.where(rng.random(p) < 0.25, rng.lognormal(1.5, 0.5, p), 0.0)


def sim_Z(Sig_tilde, s, c=None, wa=None):
    Z = np.empty((p, q))
    for a in range(q):
        Sa = Sig_tilde if wa is None else np.diag(wa[:, a])
        beta = mvn(rng, (s[a] / n[a]) * Sa, 1)[0]
        Z[:, a] = np.sqrt(n[a]) * R @ beta + mvn(rng, R, 1)[0]
    return Z


def fast_ols(Z, s, c=1.0):
    return np.linalg.solve(D, (Z ** 2 - c) @ s) / (s @ s)


def brute_ols(Z, s, c=1.0):
    X = np.vstack([s[a] * D for a in range(q)])
    y = (Z ** 2 - c).T.ravel()
    return np.linalg.lstsq(X, y, rcond=None)[0]


def fast_wls(Z, s, w0, c=1.0):
    A = np.zeros((p, p)); b = np.zeros(p)
    for a in range(q):
        v = 2 * (c + s[a] * D @ w0) ** 2
        Wa = 1 / v
        A += s[a] ** 2 * D @ (Wa[:, None] * D)
        b += s[a] * D @ (Wa * (Z[:, a] ** 2 - c))
    return np.linalg.solve(A, b)


def brute_wls(Z, s, w0, c=1.0):
    X = np.vstack([s[a] * D for a in range(q)])
    y = (Z ** 2 - c).T.ravel()
    v = 2 * (c + X @ w0) ** 2
    sw = 1 / np.sqrt(v)
    return np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]


Z = sim_Z(np.diag(w), s)
w0 = np.clip(fast_ols(Z, s), 0, None)
print("1. max|fast OLS - brute OLS| =", np.abs(fast_ols(Z, s) - brute_ols(Z, s)).max())
print("2. max|fast WLS - brute WLS| =", np.abs(fast_wls(Z, s, w0) - brute_wls(Z, s, w0)).max())

# 3. bias from local rho
W12 = np.sqrt(w)
Rb = np.eye(p)
nz = np.flatnonzero(w)
for x, y_ in zip(nz[:-1], nz[1:]):
    if y_ - x <= 4:                       # correlate causal SNPs that are close (in LD)
        Rb[x, y_] = Rb[y_, x] = 0.6
Sig = W12[:, None] * Rb * W12[None, :]
Sig = (Sig + Sig.T) / 2
if np.linalg.eigvalsh(Sig).min() < -1e-10:
    raise SystemExit("Sigma~ not PSD; change seed")
C = Sig - np.diag(np.diag(Sig))
pred_bias = np.linalg.solve(D, np.diag(R @ C @ R))
reps = 400
est = np.mean([fast_ols(sim_Z(Sig, s), s) for _ in range(reps)], axis=0)
est0 = np.mean([fast_ols(sim_Z(np.diag(w), s), s) for _ in range(reps)], axis=0)
print("3. no local rho: mean|E w_hat - w| =", np.abs(est0 - w).mean().round(3),
      "(MC noise scale", (np.abs(est0 - w).mean()).round(3), ")")
print("   local rho   : corr(observed bias, predicted D^-1 diag(RCR)) =",
      np.corrcoef(est - w, pred_bias)[0, 1].round(3),
      "; slope =", (np.polyfit(pred_bias, est - w, 1)[0]).round(3))

# 4. heterogeneous enrichment across traits
wa = np.outer(w, np.ones(q)) * rng.lognormal(0, 0.7, (p, q)) * (rng.random((p, q)) < 0.8)
target = wa @ s ** 2 / (s @ s)
est4 = np.mean([fast_ols(sim_Z(None, s, wa=wa), s) for _ in range(reps)], axis=0)
simple = wa.mean(1)
print("4. corr(E w_hat, s^2-weighted mean w_a) =", np.corrcoef(est4, target)[0, 1].round(3),
      " max abs diff =", np.abs(est4 - target).max().round(3),
      " | vs simple mean w_a: max abs diff =", np.abs(est4 - simple).max().round(3))

# 5. intercept non-identifiability
X = np.vstack([s[a] * D for a in range(q)])
I = np.kron(np.eye(q), np.ones((p, 1)))
full = np.hstack([X, I])
print("5. columns =", full.shape[1], " rank =", np.linalg.matrix_rank(full),
      " null direction ∝ (D^-1 1, -s):",
      np.allclose(full @ np.concatenate([np.linalg.solve(D, np.ones(p)), -s]), 0))
