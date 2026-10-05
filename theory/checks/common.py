"""Shared helpers for the theory-supplement verification scripts."""
import numpy as np


def ar1_ld(p, phi):
    """AR(1) LD correlation matrix, r_kl = phi^|k-l|."""
    idx = np.arange(p)
    return phi ** np.abs(idx[:, None] - idx[None, :])


def random_ld(p, rng, n_ref=None, phi=0.9):
    """A realistic-looking LD block: sample correlation of AR(1) latent haplotypes."""
    if n_ref is None:
        return ar1_ld(p, phi)
    L = np.linalg.cholesky(ar1_ld(p, phi))
    G = rng.standard_normal((n_ref, p)) @ L.T
    G = (G > 0).astype(float) + (rng.standard_normal((n_ref, p)) @ L.T > 0)
    G = (G - G.mean(0)) / G.std(0)
    return np.corrcoef(G, rowvar=False)


def block_diag(*blocks):
    p = sum(b.shape[0] for b in blocks)
    out = np.zeros((p, p))
    i = 0
    for b in blocks:
        k = b.shape[0]
        out[i:i + k, i:i + k] = b
        i += k
    return out


def mvn(rng, cov, size):
    """Draws from N(0, cov) via eigen-decomposition (tolerates PSD)."""
    vals, vecs = np.linalg.eigh(cov)
    vals = np.clip(vals, 0, None)
    A = vecs * np.sqrt(vals)
    return rng.standard_normal((size, cov.shape[0])) @ A.T


def summary_z(rng, R, beta, n, c=1.0):
    """Summary-level Z for one trait: z = sqrt(n) R beta + u, u ~ N(0, c R)."""
    return np.sqrt(n) * R @ beta + np.sqrt(c) * mvn(rng, R, 1)[0]
