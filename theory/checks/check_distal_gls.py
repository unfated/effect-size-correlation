"""Collapsed GLS for distal SNP sets (supplement §S4.6).

Ybar = sum_a s_a z_Aa z_Ba' / ||s||^2,  E Ybar = R_A Sigma~_AB R_B.
Under Cov(vec Ybar) ∝ R_B ⊗ R_A, GLS for Sigma~ on candidate sets P_A x P_B is
  C_hat = R_A[P,P]^{-1} Ybar[P_A,P_B] R_B[P,P]^{-1}.
Checks: (1) equals brute-force GLS; (2) E C_hat = R_A[P,P]^{-1} R_A[P,:] Sigma~_AB R_B[:,P] R_B[P,P]^{-1}
(projection bias when non-candidates carry signal); (3) data-driven prescreening on w_hat:
null calibration and bias under the alternative, with and without trait splitting.
"""
import numpy as np
from common import ar1_ld

rng = np.random.default_rng(23)
b, q = 20, 300
RA, RB = ar1_ld(b, 0.7), ar1_ld(b, 0.7)
s = rng.uniform(0.03, 0.12, q)

# (1) algebra
PA, PB = [3, 4, 10], [7, 15]
Y = rng.standard_normal((b, b))
X = np.kron(RB[:, PB], RA[:, PA])
Vi = np.linalg.inv(np.kron(RB, RA))
gls = np.linalg.solve(X.T @ Vi @ X, X.T @ Vi @ Y.T.ravel(order="C").reshape(-1)) if False else None
y = Y.ravel(order="F")
gls = np.linalg.solve(X.T @ Vi @ X, X.T @ Vi @ y).reshape(len(PA), len(PB), order="F")
coll = np.linalg.solve(RA[np.ix_(PA, PA)], Y[np.ix_(PA, PB)]) @ np.linalg.inv(RB[np.ix_(PB, PB)])
print("(1) max|collapsed - brute GLS| =", np.abs(gls - coll).max())


def simulate(wA, wB, CAB):
    Sig = np.block([[np.diag(wA), CAB], [CAB.T, np.diag(wB)]])
    L = np.linalg.cholesky(Sig + 1e-10 * np.eye(2 * b))
    beta = (L @ rng.standard_normal((2 * b, q))) * np.sqrt(s)
    ZA = RA @ beta[:b] + np.linalg.cholesky(RA) @ rng.standard_normal((b, q))
    ZB = RB @ beta[b:] + np.linalg.cholesky(RB) @ rng.standard_normal((b, q))
    return ZA, ZB


def collapsed(ZA, ZB, PA, PB, idx=None):
    idx = np.arange(q) if idx is None else idx
    Yb = (ZA[:, idx] * s[idx]) @ ZB[:, idx].T / (s[idx] @ s[idx])
    return np.linalg.solve(RA[np.ix_(PA, PA)], Yb[np.ix_(PA, PB)]) @ np.linalg.inv(RB[np.ix_(PB, PB)])


def what(Z, R, idx):
    D = R ** 2
    return np.linalg.solve(D, (Z[:, idx] ** 2 - 1) @ s[idx]) / (s[idx] @ s[idx])


# (2) off-candidate signal: true pair (5, 9) but candidates {6} x {9}
wA = np.ones(b); wB = np.ones(b); wA[5] = wB[9] = 20; wA[6] = 20
CAB = np.zeros((b, b)); CAB[5, 9] = 0.5 * 20
est = np.mean([collapsed(*simulate(wA, wB, CAB), [6], [9])[0, 0] for _ in range(1500)])
pred = (np.linalg.solve(RA[np.ix_([6], [6])], RA[[6], :]) @ CAB @ RB[:, [9]] @ np.linalg.inv(RB[np.ix_([9], [9])]))[0, 0]
print(f"(2) off-candidate signal: mean C_hat(6,9) = {est:.2f}, predicted projection {pred:.2f}, true C(6,9) = 0")

# (3) data-driven prescreen: top-3 w_hat per block
def run(rho, split, reps=800):
    out = []
    for _ in range(reps):
        CAB = np.zeros((b, b)); CAB[5, 9] = rho * 20
        ZA, ZB = simulate(wA, wB, CAB)
        perm = rng.permutation(q)
        i1, i2 = (perm[: q // 2], perm[q // 2:]) if split else (np.arange(q), np.arange(q))
        PA = list(np.argsort(-what(ZA, RA, i1))[:3]); PB = list(np.argsort(-what(ZB, RB, i1))[:3])
        if 5 not in PA or 9 not in PB:
            continue
        C = collapsed(ZA, ZB, PA, PB, i2)
        out.append(C[PA.index(5), PB.index(9)])
    return np.mean(out), np.std(out), len(out)
for rho in (0.0, 0.5):
    for split in (False, True):
        mu, sd, n = run(rho, split)
        print(f"(3) rho={rho} split={split}: mean C_hat = {mu:.2f} (true {rho*20:.1f}), sd {sd:.2f}, kept {n}")
