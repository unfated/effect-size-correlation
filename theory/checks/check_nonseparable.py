"""LRCP pair tests when locus effects are non-separable (supplement §S8.5).

Mediator model with FIXED mediator->trait profiles q_r (r = 1..L) and random locus->mediator effects p.
Genome-wide architecture: loci load on all L mediators -> trait covariance T = sum_r q_r q_r'.
"Lipid-class" loci load only on mediators 1-2 (a narrow, correlated trait domain): V_class = sum_{r<=2} q_r q_r'.
Single tag per locus (no LD); traits fully overlapping (Gamma), q = 100.
Statistic: C_hat = x'Sy/||s||^2 with x, y the two loci's Z profiles.
Null variances (conditional on y, §S8.2):
  separable   Sigma_x = Gamma + kappa_A * T_hat          (T_hat = genome-wide, trace-matched to locus A)
  class       Sigma_x = Gamma + kappa_A * V_hat_class    (V_hat_class from other class loci, held out)
  own-block   Sigma_x = x x' (locus covariance from its own realised profile): self-normalising, |z| ~ 1
H0: A and B are independent lipid-class loci (random-effect rho = 0, realised profiles nearly collinear).
H1: B's mediator loadings are proportional to A's (rho = 1 within the class space).
"""
import numpy as np

rng = np.random.default_rng(17)
q, L = 100, 25
cl = np.repeat(np.arange(10), 10)
Gam = np.where(cl[:, None] == cl[None, :], 0.4, 0.05); np.fill_diagonal(Gam, 1.0)
LGam = np.linalg.cholesky(Gam)
Q = rng.standard_normal((q, L)) * (rng.random((q, L)) < 0.15)          # sparse mediator->trait profiles
Q[:, 0] = 0; Q[cl == 0, 0] = rng.uniform(0.8, 1.2, 10)                 # mediator 1: "LDL" cluster
Q[:, 1] = 0; Q[cl == 0, 1] = rng.normal(0, .3, 10); Q[cl == 1, 1] = rng.uniform(.5, 1, 10)   # mediator 2: related
s = np.full(q, 1.0)
amp = 8.0                                                               # per-trait noncentrality scale of a class locus


def locus(cls, p=None):
    if p is None:
        p = rng.standard_normal(2 if cls else L)
    beta = Q[:, :2] @ p if cls else Q @ p / np.sqrt(L / 2)
    return amp * beta, p


T = Q @ Q.T / (L / 2)
V_true = Q[:, :2] @ Q[:, :2].T


def noise():
    return LGam @ rng.standard_normal(q)


def v_class_hat(k=40):
    """Average of (x x' - Gamma) over k held-out class loci, i.e. amp^2 V_class + noise."""
    acc = np.zeros((q, q))
    for _ in range(k):
        mu, _ = locus(True); x = mu + noise(); acc += np.outer(x, x) - Gam
    return acc / k


def ztest(x, y, V):
    """Conditional-on-y test; V is the locus trait covariance scaled to x's realised signal size."""
    kap = max((x @ x - np.trace(Gam)) / np.trace(V), 1e-6)
    Sx = Gam + kap * V
    c = x @ (s * y)
    return c / np.sqrt((s * y) @ Sx @ (s * y))


def qeff(V):
    return np.trace(V) ** 2 / np.trace(V @ V)


print(f"q_eff: genome-wide T {qeff(T):.1f}, lipid class {qeff(V_true):.1f}")
reps = 2000
out = {}
for hyp in ["H0", "H1 (rho=1)", "H1 (rho=0.7)"]:
    zs, zc, zo, rr = [], [], [], []
    Vc = v_class_hat()
    for _ in range(reps):
        muA, pA = locus(True)
        if hyp == "H0":
            muB, _ = locus(True)
        else:
            rho = 1.0 if "rho=1)" in hyp else 0.7
            pB = rho * pA * rng.choice([-1, 1]) + np.sqrt(1 - rho ** 2) * rng.standard_normal(2)
            muB, _ = locus(True, pB)
        x, y = muA + noise(), muB + noise()
        zs.append(ztest(x, y, T)); zc.append(ztest(x, y, Vc))
        zo.append((x @ (s * y)) / np.sqrt((s * y) @ (Gam + np.outer(x, x)) @ (s * y)))
        rr.append(np.corrcoef(muA, muB)[0, 1])
    zs, zc, zo = np.abs(zs), np.abs(zc), np.abs(zo)
    print(f"{hyp:13s} realised |profile corr| median {np.median(np.abs(rr)):.2f} | "
          f"reject 5%: separable {np.mean(zs > 1.96):.3f}, class null {np.mean(zc > 1.96):.3f} | "
          f"own-block {np.mean(zo > 1.96):.3f} (max |z| {zo.max():.2f}) | median |z| {np.median(zs):.1f} vs {np.median(zc):.2f}")
