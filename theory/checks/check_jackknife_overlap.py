"""Why the delete-one-trait-cluster jackknife fails under sample overlap (supplement §S3.10).

One null SNP (no LD), q = 300 traits in 20 r_g clusters; full sample overlap with phenotypic
correlation within (r_w) and BETWEEN (r_b) clusters. w_hat = sum_a s_a (z_a^2 - 1)/||s||^2.
The jackknife treats clusters as independent replicates, so it estimates only the within-cluster
part of Var(w_hat); the between-cluster covariance 2 sum_{k<l} Cov(cluster k, cluster l) is missed.
"""
import numpy as np

rng = np.random.default_rng(71)
q, K = 300, 20
sizes = rng.multinomial(q - K, np.ones(K) / K) + 1
cl = np.repeat(np.arange(K), sizes)
s = rng.uniform(0.02, 0.06, q); ss = s @ s
for rw, rb in [(0.4, 0.0), (0.4, 0.1), (0.4, 0.2)]:
    Gam = np.where(cl[:, None] == cl[None, :], rw, rb); np.fill_diagonal(Gam, 1)
    L = np.linalg.cholesky(Gam)
    model = np.sqrt(2 * s @ (Gam ** 2) @ s) / ss
    within = np.sqrt(2 * s @ (Gam ** 2 * (cl[:, None] == cl[None, :])) @ s) / ss
    est, jk = [], []
    for _ in range(3000):
        z = L @ rng.standard_normal(q)
        y = s * (z ** 2 - 1)
        est.append(y.sum() / ss)
        loo = np.array([y[cl != k].sum() / (s[cl != k] @ s[cl != k]) for k in range(K)])
        jk.append(np.sqrt((K - 1) / K * np.sum((loo - loo.mean()) ** 2)))
    est = np.array(est)
    print(f"r_w={rw} r_b={rb}: SD emp {est.std():.2f}, model (3.5) {model:.2f}, "
          f"within-cluster only {within:.2f}, jackknife median {np.median(jk):.2f}")
