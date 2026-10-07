"""Figure 4: genome-wide LRCP. (a) QQ of the conditional screen over 30.2 M distal tag pairs;
(b) modules at screen P < 1e-5 with their lead trait axis.
usage: fig4.py <genomewide_results_dir> <out.png>"""
import sys, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D, OUT = sys.argv[1:3]
qq = pd.read_csv(f"{D}/screen_qq.tsv", sep="\t"); M = pd.read_csv(f"{D}/modules_p1e-05.tsv", sep="\t")
fig, ax = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1, 1.6]})
x = -np.log10(qq.exp); y = -np.log10(qq.p)
ax[0].plot(x, y, ".", ms=3, color="C0", label="LRCP conditional test")
ax[0].plot([0, x.max()], [0, x.max()], "k--", lw=0.8)
ax[0].axhline(-np.log10(0.05 / 30215516), color="grey", lw=0.6, ls=":")
ax[0].set_xlabel("expected -log10 P"); ax[0].set_ylabel("observed -log10 P")
ax[0].set_title("a  30.2 M distal pairs of 7,793 enriched tags", loc="left", fontsize=10)
ax[0].text(0.3, 12, "naive Z correlation:\n77% of pairs P < 0.05\n11.1 M Bonferroni", fontsize=8)
M = M.head(12)[::-1]
lab = [f"{r.module}. " + ",".join(dict.fromkeys(r.genes.split(","))) for r in M.itertuples()]
lab = [l if len(l) < 60 else l[:57] + "..." for l in lab]
ax[1].barh(range(len(M)), M.n_loci, color=["C3" if q < 5 else "C2" for q in M.min_q_eff_class.fillna(99)])
ax[1].set_yticks(range(len(M))); ax[1].set_yticklabels(lab, fontsize=7)
for i, r in enumerate(M.itertuples()):
    ax[1].text(r.n_loci + 0.2, i, r.top_traits.split(";")[0], va="center", fontsize=7)
ax[1].set_xlabel("loci in module"); ax[1].set_xlim(0, 24)
ax[1].set_title("b  modules (P < 1e-5 edges); red = narrow trait class (descriptive)", loc="left", fontsize=10)
plt.tight_layout(); plt.savefig(OUT, dpi=150)
