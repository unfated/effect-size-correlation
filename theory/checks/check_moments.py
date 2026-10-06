"""Fourth moments and mediator models (supplement §S2.4, §S3.4, §S5).

1. Normal effects: Var(z_k^2) = 2 sigma_k^4, Cov(z_k^2, z_l^2) = 2 sigma_kl^2,
   Var(z_k z_l) = sigma_k^2 sigma_l^2 + sigma_kl^2.
2. Point-normal (spike-and-slab) independent effects with causal prob pi:
   Cov(z_k^2, z_l^2) = 2 sigma_kl^2 + n^2 sum_s r_ks^2 r_ls^2 kappa_s,
   kappa_s = 3 v_s^2 (1 - pi) / pi  (fourth cumulant of beta_s, v_s = Var beta_s).
3. One-mediator model B = p q' (p ~ N(0,U), q ~ N(0,V)) has Cov(vec B) = V ⊗ U
   (same second moments as MN(0,U,V)) but is not matrix normal: E B_ia^4 = 9 U_ii^2 V_aa^2.
"""
import numpy as np
from common import ar1_ld

rng = np.random.default_rng(5)
m, n = 5, 1.0
R = ar1_ld(m, 0.7); L = np.linalg.cholesky(R)
v = np.full(m, 2.0)                       # n * Var(beta_s) (n absorbed)
reps = 2_000_000
for name, pi in (("normal", 1.0), ("point-normal pi=0.1", 0.1)):
    g = rng.standard_normal((reps, m)) * np.sqrt(v / pi) * (rng.random((reps, m)) < pi)
    z = g @ R.T + rng.standard_normal((reps, m)) @ L.T
    S = R + R @ np.diag(v) @ R
    kappa = 3 * v ** 2 * (1 - pi) / pi
    pred = 2 * S ** 2 + (R ** 2) @ np.diag(kappa) @ (R ** 2).T
    emp = np.cov(z ** 2, rowvar=False)
    pred_x = S[0, 0] * S[1, 1] + S[0, 1] ** 2 + ((R[0] * R[1]) ** 2 * kappa).sum()
    print(f"{name:22s} Cov(z^2) max rel err = {np.abs(emp / pred - 1).max():.3f}; "
          f"Var(z0 z1) emp {np.var(z[:, 0] * z[:, 1]):.2f} pred {pred_x:.2f}; "
          f"slide formula 2*sigma^2 = {2 * S[0, 0]:.2f} vs emp Var(z0^2) {emp[0, 0]:.2f}")

U = np.array([[1, .5], [.5, 1]]); V = np.array([[2, -.6], [-.6, 1]])
p = rng.standard_normal((reps, 2)) @ np.linalg.cholesky(U).T
qv = rng.standard_normal((reps, 2)) @ np.linalg.cholesky(V).T
Bv = np.einsum('ri,ra->ria', p, qv).reshape(reps, 4)        # vec over (i,a) i-major
print("mediator: max|Cov(vec B) - U⊗V| =", np.abs(np.cov(Bv, rowvar=False) - np.kron(U, V)).max().round(3),
      "| E B11^4 =", np.mean(Bv[:, 0] ** 4).round(2), "(MN would give", 3 * (U[0, 0] * V[0, 0]) ** 2,
      "; one-mediator theory", 9 * (U[0, 0] * V[0, 0]) ** 2, ")")
