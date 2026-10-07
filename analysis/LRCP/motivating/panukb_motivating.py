"""Motivating analysis in Pan-UKB (EUR): phenome-wide Z correlations of the
TP53/MDM2 SNPs and lipid control SNPs, judged against an empirical null of
distal SNP pairs from random 20-kb windows, and against q_eff.

Inputs: per-trait region TSVs from fetch_panukb_regions.py.
Outputs (shared folder): tables + figures in results/motivating/.
"""
import os, sys, glob, itertools, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

IN, OUT = sys.argv[1], sys.argv[2]
traits = pd.read_csv(f"{IN}/traits.tsv", sep="\t", low_memory=False)
traits["key"] = traits.filename.str.replace(".tsv.bgz", "", regex=False)
regions = pd.read_csv(f"{IN}/regions.tsv", sep="\t")

# ---- build Z matrix (variant x trait) ----
Zs, AF = {}, None
for f in sorted(glob.glob(f"{IN}/*.tsv")):
    k = os.path.basename(f)[:-4]
    if k in ("traits", "regions"): continue
    d = pd.read_csv(f, sep="\t", dtype={"chr": str})
    d["id"] = d.region + ":" + d.chr + ":" + d.pos.astype(str) + ":" + d.ref + ":" + d.alt
    d = d.drop_duplicates("id").set_index("id")
    Zs[k] = d.beta / d.se
    if k == "continuous-50-both_sexes-irnt" or AF is None:
        AF = d.af if AF is None or k == "continuous-50-both_sexes-irnt" else AF
Z = pd.DataFrame(Zs)
keep_traits = [k for k in traits.key if k in Z.columns]
Z = Z[keep_traits]
ok = Z.notna().mean(1) > 0.98
af = pd.to_numeric(AF.reindex(Z.index), errors="coerce")
EXACT = ("TP53_rs78378222:17:7571752:", "MDM2_rs3730556:12:69216521:", "LDLR_rs6511720:19:11202306:",
         "PCSK9_rs11591147:1:55505647:", "HMGCR_rs12916:5:74656539:", "APOB_rs1367117:2:21263900:")
is_target = pd.Series(Z.index, index=Z.index).str.startswith(EXACT)  # named SNPs kept even if rare (TP53, PCSK9 ~1%)
Z = Z[ok & (af.between(0.05, 0.95) | (is_target & af.between(0.005, 0.995)))]
Z = Z.apply(lambda c: c.fillna(0.0))  # rare missing -> 0 (no information)
print("Z matrix", Z.shape)
reg_of = Z.index.str.split(":").str[0]
chrom = Z.index.str.split(":").str[1].astype(int); pos = Z.index.str.split(":").str[2].astype(int)
q = Z.shape[1]

# ---- trait-level null covariance C from random windows (mostly null SNPs) ----
rand = reg_of.str.startswith("rand")
chi_mean = (Z ** 2).mean(1)
null_snps = rand & (chi_mean < 1.5)
Zn = Z[null_snps].values
C_hat = np.corrcoef(Zn.T)
qeff_null = q ** 2 / np.sum(C_hat ** 2)
print("q =", q, " q_eff from null-SNP trait correlation:", round(qeff_null, 1))

# ---- one SNP per random window: random and lead (max mean chi2) ----
rng = np.random.default_rng(11)
pick_rand, pick_lead = [], []
for r, idx in pd.Series(np.arange(len(Z)), index=Z.index).groupby(reg_of.values):
    if not r.startswith("rand"): continue
    ids = idx.values
    pick_rand.append(rng.choice(ids)); pick_lead.append(ids[np.argmax(chi_mean.values[ids])])

def pair_stats(ix, label, Zm):
    rows = []
    X = Zm[ix]; lp = -stats.norm.logsf(np.abs(X)) / np.log(10) - np.log10(2)
    Xc = X - X.mean(1, keepdims=True); Xc /= np.sqrt((Xc ** 2).sum(1, keepdims=True))
    L = lp - lp.mean(1, keepdims=True); L /= np.sqrt((L ** 2).sum(1, keepdims=True))
    Rz, Rl = Xc @ Xc.T, L @ L.T
    for a, b in itertools.combinations(range(len(ix)), 2):
        i, j = ix[a], ix[b]
        if chrom[i] == chrom[j] and abs(pos[i] - pos[j]) < 5_000_000: continue
        rows.append((label, Z.index[i], Z.index[j], Rz[a, b], Rl[a, b], chi_mean.values[i], chi_mean.values[j]))
    return rows

def run(cols, tag):
    Zm = Z[cols].values; qq = len(cols)
    C = np.corrcoef(Zm[null_snps.values].T); qe = qq ** 2 / np.sum(C ** 2)
    rows = pair_stats(np.array(pick_rand), "random", Zm) + pair_stats(np.array(pick_lead), "lead", Zm)
    nul = pd.DataFrame(rows, columns=["pool", "snp1", "snp2", "r_z", "r_neglog10p", "chi1", "chi2"])
    nul["min_chi"] = nul[["chi1", "chi2"]].min(axis=1)
    # targets
    def lead_of(region):
        idx = np.where(reg_of == region)[0]
        return idx[np.argmax(chi_mean.values[idx])] if len(idx) else None
    named = {}
    for name, exact in [("TP53_rs78378222", "17:7571752"), ("MDM2_rs3730556", "12:69216521"),
                        ("LDLR_rs6511720", "19:11202306"), ("PCSK9_rs11591147", "1:55505647"),
                        ("HMGCR_rs12916", "5:74656539"), ("APOB_rs1367117", "2:21263900")]:
        hits = [i for i, s in enumerate(Z.index) if s.startswith(name + ":" + exact + ":")]
        named[name] = hits[0] if hits else lead_of(name)
    tgt = []
    for a, b in itertools.combinations(named, 2):
        i, j = named[a], named[b]
        if i is None or j is None: continue
        x, y = Zm[i], Zm[j]
        r = np.corrcoef(x, y)[0, 1]
        lx = -stats.norm.logsf(np.abs(x)); ly = -stats.norm.logsf(np.abs(y))
        mc = min(chi_mean.values[i], chi_mean.values[j])
        # empirical null: distal pairs whose weaker SNP has similar signal (within factor 2), pooled
        ref = nul[(nul.min_chi > mc / 2) & (nul.min_chi < mc * 2)]
        if len(ref) < 200: ref = nul.reindex((nul.min_chi - mc).abs().sort_values().index[:500])
        p_emp = (np.sum(np.abs(ref.r_z) >= abs(r)) + 1) / (len(ref) + 1)
        p_naive = 2 * stats.norm.sf(abs(np.arctanh(r)) * np.sqrt(qq - 3))
        p_qeff = 2 * stats.norm.sf(abs(np.arctanh(r)) * np.sqrt(max(qe - 3, 1)))
        tgt.append(dict(set=tag, snp1=a, snp2=b, id1=Z.index[i], id2=Z.index[j], q=qq, q_eff_null=qe,
                        r_z=r, r_neglog10p=np.corrcoef(lx, ly)[0, 1], mean_chi1=chi_mean.values[i], mean_chi2=chi_mean.values[j],
                        p_naive=p_naive, p_qeff=p_qeff, n_null_ref=len(ref), null_sd_ref=ref.r_z.std(), p_empirical=p_emp))
    nsum = nul.groupby("pool").agg(n=("r_z", "size"), mean_r_z=("r_z", "mean"), sd_r_z=("r_z", "std"),
                                   mean_r_neglog10p=("r_neglog10p", "mean"), median_min_chi=("min_chi", "median")).reset_index()
    nsum["set"] = tag; nsum["q"] = qq; nsum["sd_expected_naive"] = 1 / np.sqrt(qq - 3); nsum["sd_expected_qeff_null"] = 1 / np.sqrt(qe)
    nsum["typeI_naive"] = [np.mean(np.abs(np.arctanh(nul[nul.pool == p].r_z)) * np.sqrt(qq - 3) > 1.96) for p in nsum.pool]
    return pd.DataFrame(tgt), nsum, nul

all_cols = list(Z.columns)
indep = [k for k in all_cols if bool(traits.set_index("key").loc[k, "in_max_independent_set"])]
t1, s1, n1 = run(all_cols, "all_EUR_PASS")
t2, s2, n2 = run(indep, "max_independent_set")
T = pd.concat([t1, t2]); S = pd.concat([s1, s2])
T.to_csv(f"{OUT}/panukb_target_pairs.tsv", sep="\t", index=False, float_format="%.4g")
S.to_csv(f"{OUT}/panukb_null_summary.tsv", sep="\t", index=False, float_format="%.4g")
n1.to_csv(f"{OUT}/panukb_null_pairs_all.tsv.gz", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250)
print(S.to_string()); print(T[["set", "snp1", "snp2", "q", "q_eff_null", "r_z", "r_neglog10p", "mean_chi1", "mean_chi2", "p_naive", "p_qeff", "p_empirical", "null_sd_ref"]].to_string())

# ---- figure: null distributions ----
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
for pool, c in (("random", "#999999"), ("lead", "#C44E52")):
    d = n1[n1.pool == pool]
    ax[0].hist(d.r_z, bins=60, density=True, alpha=.6, color=c, label=f"{pool} SNPs (SD {d.r_z.std():.3f})")
    ax[1].hist(d.r_neglog10p, bins=60, density=True, alpha=.6, color=c, label=f"{pool} (mean {d.r_neglog10p.mean():.3f})")
xx = np.linspace(-.4, .4, 200)
ax[0].plot(xx, stats.norm.pdf(xx, 0, 1 / np.sqrt(len(all_cols) - 3)), "k--", lw=1, label=f"naive N(0,1/(q-3)), q={len(all_cols)}")
ax[0].set_title("a  signed Z: distal null pairs"); ax[0].set_xlabel("r across traits"); ax[0].legend(fontsize=7, frameon=False)
ax[1].axvline(0, c="k", lw=.6); ax[1].set_title("b  −log10 P: distal null pairs"); ax[1].set_xlabel("r across traits"); ax[1].legend(fontsize=7, frameon=False)
b = n1.copy(); b["bin"] = pd.qcut(np.log(b.min_chi), 8)
g = b.groupby("bin", observed=True).agg(x=("min_chi", "median"), sd=("r_z", "std"), ml=("r_neglog10p", "mean"))
ax[2].plot(g.x, g.sd, "o-", label="SD of signed-Z r"); ax[2].plot(g.x, g.ml, "s-", label="mean of −log10P r")
ax[2].axhline(1 / np.sqrt(len(all_cols) - 3), c="k", ls="--", lw=.8, label="naive SD")
ax[2].set_xscale("log"); ax[2].set_xlabel("weaker SNP's mean χ² across traits"); ax[2].set_title("c  null spread grows with SNP signal"); ax[2].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/fig_panukb_null.png", dpi=200)
