"""Fig 1a-b: TP53 rs78378222 vs MDM2 rs3730556 Z scores across the CMD Knowledge
Portal PheWAS traits (meta-analysed), overall and by trait group."""
import sys, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D = sys.argv[1]
a = pd.read_csv(f"{D}/phewas_SNP1_meta.tsv", sep="\t"); b = pd.read_csv(f"{D}/phewas_SNP2_meta.tsv", sep="\t")
for d in (a, b): d["z"] = d.beta / d.se
m = a.merge(b, on=["pheno", "group", "desc"], suffixes=("_tp53", "_mdm2"))
m.to_csv(f"{D}/tp53_mdm2_matched.tsv", sep="\t", index=False)

def summ(d, label):
    r, p = stats.pearsonr(d.z_tp53, d.z_mdm2); s, sp = stats.spearmanr(d.z_tp53, d.z_mdm2)
    # unsigned statistic, as in the MAGMA figure
    rl, _ = stats.pearsonr(-np.log10(d.p_tp53), -np.log10(d.p_mdm2))
    return dict(group=label, n_traits=len(d), pearson_z=r, p_naive=p, spearman_z=s, p_spearman=sp, pearson_neglog10p=rl)
rows = [summ(m, "ALL")] + [summ(d, g) for g, d in m.groupby("group") if len(d) >= 8]
# leave-top-trait-out sensitivity: drop the 5 traits with largest |z_tp53|
top = m.reindex(m.z_tp53.abs().sort_values(ascending=False).index).iloc[5:]
rows.append(summ(top, "ALL minus top-5 |z_TP53| traits"))
res = pd.DataFrame(rows); res.to_csv(f"{D}/tp53_mdm2_correlations.tsv", sep="\t", index=False, float_format="%.4g")
print(res.to_string())

groups = m.group.value_counts().index.tolist()
cmap = plt.get_cmap("tab10")
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw=dict(width_ratios=[1.2, 1]))
for i, g in enumerate(groups):
    d = m[m.group == g]; ax[0].scatter(d.z_tp53, d.z_mdm2, s=18, color=cmap(i % 10), label=g.title(), alpha=.85)
ax[0].axhline(0, c="grey", lw=.5); ax[0].axvline(0, c="grey", lw=.5)
ax[0].set_xlabel("Z, rs78378222 (TP53, chr17)"); ax[0].set_ylabel("Z, rs3730556 (MDM2, chr12)")
ax[0].set_title(f"a  {len(m)} traits: r = {res.pearson_z[0]:.2f}")
ax[0].legend(fontsize=6.5, frameon=False, loc="lower right")
s = res[(res.group != "ALL") & ~res.group.str.startswith("ALL")].sort_values("pearson_z")
ax[1].barh(s.group.str.title(), s.pearson_z, color="#4C72B0")
for y, (r, n) in enumerate(zip(s.pearson_z, s.n_traits)): ax[1].text(r + (0.02 if r >= 0 else -0.02), y, f"n={n}", va="center", ha="left" if r >= 0 else "right", fontsize=7)
ax[1].axvline(0, c="k", lw=.6); ax[1].set_xlabel("Pearson r of Z across traits"); ax[1].set_title("b  by trait group")
ax[1].set_xlim(-1, 1)
fig.tight_layout(); fig.savefig(f"{D}/fig1ab_tp53_mdm2.png", dpi=200)
