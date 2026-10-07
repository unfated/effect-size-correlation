"""Null behaviour of naive phenome-wide correlations between two distal SNPs.

Two SNPs k and l in different LD blocks, genetic-effect correlation rho = 0.
Across q traits their Z vectors are independent Gaussians with covariance
M_k = n*tau_k*G + C and M_l = n*tau_l*G + C, where G is the (scaled) trait
genetic covariance, C the inter-GWAS null correlation (sample overlap x
phenotypic correlation, the LDSC cross-trait intercepts) and tau_k the
per-SNP expected marginal variance (LD score x enrichment x h2/m).

We compare three naive statistics across traits:
  * Pearson correlation of signed Z           (what LRCP targets)
  * Pearson correlation of chi-square / -log10 P (what Figure 6.1 used)
and show the effective number of traits q_eff = tr(M_k)tr(M_l)/tr(M_k M_l)
that sets the null spread of the signed correlation.
"""
import numpy as np, pandas as pd, sys
from scipy import stats

rng = np.random.default_rng(2026)

def trait_cov(q, n_clusters, within_rg, overlap_r):
    """Clustered genetic correlation G (unit diagonal) and intercept matrix C."""
    lab = np.arange(q) % n_clusters
    same = lab[:, None] == lab[None, :]
    G = np.where(same, within_rg, 0.0); np.fill_diagonal(G, 1.0)
    C = np.where(same, overlap_r, 0.0); np.fill_diagonal(C, 1.0)
    return G, C

def simulate(q=300, n_clusters=30, within_rg=0.6, overlap_r=0.3, power_sd=1.0,
             signal=(2.0, 2.0), reps=4000):
    """signal = expected mean chi-square excess of SNP k and l (n*tau*avg power)."""
    G, C = trait_cov(q, n_clusters, within_rg, overlap_r)
    # trait-level power heterogeneity: n_a*h2_a varies log-normally across traits
    pw = np.exp(rng.normal(0, power_sd, q)); pw /= pw.mean()
    S = np.sqrt(pw)[:, None] * G * np.sqrt(pw)[None, :]
    Mk, Ml = signal[0] * S + C, signal[1] * S + C
    Lk, Ll = np.linalg.cholesky(Mk), np.linalg.cholesky(Ml)
    zk = rng.standard_normal((reps, q)) @ Lk.T
    zl = rng.standard_normal((reps, q)) @ Ll.T
    def rowcor(a, b):
        a = a - a.mean(1, keepdims=True); b = b - b.mean(1, keepdims=True)
        return (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1))
    r_z = rowcor(zk, zl)
    r_chi = rowcor(zk ** 2, zl ** 2)
    lp_k = -stats.norm.logsf(np.abs(zk)) / np.log(10); lp_l = -stats.norm.logsf(np.abs(zl)) / np.log(10)
    r_lp = rowcor(lp_k, lp_l)
    q_eff = np.trace(Mk) * np.trace(Ml) / np.trace(Mk @ Ml)
    return dict(q=q, n_clusters=n_clusters, within_rg=within_rg, overlap_r=overlap_r, power_sd=power_sd,
                signal=signal[0], q_eff=q_eff,
                sd_rz=r_z.std(), sd_rz_theory=1 / np.sqrt(q_eff), sd_rz_naive=1 / np.sqrt(q - 3),
                typeI_naive=np.mean(np.abs(r_z) > 1.96 / np.sqrt(q - 3)),
                typeI_qeff=np.mean(np.abs(r_z) > 1.96 / np.sqrt(q_eff)),
                mean_rchi=r_chi.mean(), mean_rlogp=r_lp.mean(), mean_rz=r_z.mean())

if __name__ == "__main__":
    out = sys.argv[1]
    rows = []
    for q in (30, 100, 300):
        for nc, rg in ((q, 0.0), (max(q // 10, 3), 0.3), (max(q // 10, 3), 0.6), (max(q // 30, 1), 0.8)):
            for ov in (0.0, 0.3):
                for sig in (0.0, 0.5, 2.0, 10.0):
                    for psd in (0.0, 1.0):
                        rows.append(simulate(q, nc, rg, ov, psd, (sig, sig), reps=2000))
                        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False, float_format="%.4g")
