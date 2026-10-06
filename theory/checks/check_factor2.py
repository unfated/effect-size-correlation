"""Factor-of-2 question (supplement §S1.5).

Re-runs the 2022 cross-block toy simulation (N=10,000, h2=0.8, two blocks of
50 SNPs, q=200 GWAS, all SNPs causal, one non-zero genetic-effect correlation
rho between the top-LD-score SNP of each block) and regresses the cross-block
Z products on three candidate regressors:

  correct  : g = N h2/M * (r_ki r_lj + r_kj r_li)   (pair counted once)
  chapter2 : g = N h2/M * 2 r_ki r_lj               (pair counted twice)

Expected: slope(correct) ~ rho, slope(chapter2) ~ rho/2.
"""
import numpy as np
from common import ar1_ld, block_diag, mvn

rng = np.random.default_rng(2022)
N, h2, b, q = 10_000, 0.8, 50, 200
NREP = 400
M = 2 * b
R1, R2 = ar1_ld(b, 0.85), ar1_ld(b, 0.8)
R = block_diag(R1, R2)
ld = (R ** 2).sum(0)
i = int(np.argmax(ld[:b]))
j = b + int(np.argmax(ld[b:]))

res = []
for rho in (0.8, 0.5, 0.2, -0.5):
    slopes = {"correct": [], "chapter2": []}
    for rep in range(NREP):
        S = np.eye(M)
        S[i, j] = S[j, i] = rho
        B = mvn(rng, h2 / M * S, q).T          # M x q
        U = mvn(rng, R, q).T
        Z = np.sqrt(N) * R @ B + U
        y = (Z[:b] @ Z[b:].T / q).ravel()      # mean cross-block product, k in block1, l in block2
        base = N * h2 / M
        g_ok = base * (np.outer(R[:b, i], R[b:, j]) + np.outer(R[:b, j], R[b:, i])).ravel()
        g_ch = base * 2 * np.outer(R[:b, i], R[b:, j]).ravel()
        for key, g in (("correct", g_ok), ("chapter2", g_ch)):
            X = np.column_stack([np.ones_like(g), g])
            slopes[key].append(np.linalg.lstsq(X, y, rcond=None)[0][1])
    res.append((rho, np.mean(slopes["correct"]), np.std(slopes["correct"]) / np.sqrt(NREP),
                np.mean(slopes["chapter2"])))

print("rho_true  slope_correct (MC se)  slope_chapter_form")
for r in res:
    print(f"{r[0]:8.2f}  {r[1]:8.3f} ({r[2]:.3f})      {r[3]:8.3f}")
