"""Fig. 3 draft: targeted LRCP for the 15 named-variant pairs."""
import sys, pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D = sys.argv[1]
s = pd.read_csv(f"{D}/stage4_pairs_z3.tsv", sep="\t"); s = s[s.target_pair].reset_index(drop=True)
b = pd.read_csv(f"{D}/blockpair_z3.tsv", sep="\t")
d = s.merge(b, on=["block1", "block2"])
lip = {"PCSK9", "APOB", "HMGCR", "LDLR"}
d["kind"] = np.where(d.block1.isin(lip) & d.block2.isin(lip), "lipid", np.where({*()} == set(), "", ""))
d.loc[(d.block1 + d.block2).isin(["MDM2TP53"]), "kind"] = "TP53-MDM2"
d.loc[d.kind == "", "kind"] = "cross (null)"
d["label"] = d.block1 + " × " + d.block2
d = d.sort_values(["kind", "rho"]).reset_index(drop=True)
col = {"lipid": "#1f6f8b", "TP53-MDM2": "#d1495b", "cross (null)": "#999"}
fig, ax = plt.subplots(1, 2, figsize=(11, 5.5), gridspec_kw={"width_ratios": [1.4, 1]})
y = np.arange(len(d))
lo = np.where(np.isfinite(d.se_plugin), d.rho - 1.96 * d.se_plugin, d.rho - 1.96 * d.se_null)
hi = np.where(np.isfinite(d.se_plugin), d.rho + 1.96 * d.se_plugin, d.rho + 1.96 * d.se_null)
for k, c in col.items():
    m = d.kind == k
    ax[0].errorbar(d.rho[m], y[m], xerr=[d.rho[m] - lo[m], hi[m] - d.rho[m]], fmt="o", color=c, label=k, capsize=2)
    ax[0].scatter(d.r_naive[m], y[m], marker="x", color="k", s=18)
ax[0].axvline(0, color="k", lw=.5); ax[0].set_yticks(y); ax[0].set_yticklabels(d.label, fontsize=8)
ax[0].set(xlabel="ρ̂ (LRCP-GLS, 95% CI); × = naive r", title="a  Named-variant pairs, 444 Pan-UKB traits", xlim=(-2.6, 2.6))
ax[0].legend(fontsize=7, loc="lower right")
ax[1].barh(y, -np.log10(d.p_Q), color=[col[k] for k in d.kind])
ax[1].axvline(-np.log10(0.05), color="k", lw=.5, ls="--"); ax[1].set_yticks(y); ax[1].set_yticklabels([])
ax[1].set(xlabel="−log10 P, block-pair Q test", title="b  Block-pair test")
plt.tight_layout(); plt.savefig(f"{D}/fig3_targeted.png", dpi=150)
