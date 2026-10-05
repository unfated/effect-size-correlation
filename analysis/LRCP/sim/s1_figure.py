"""Fig. 2 draft from S1 summary tables."""
import sys, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D = sys.argv[1]
S = pd.read_csv(f"{D}/s1_estimation.tsv", sep="\t"); T = pd.read_csv(f"{D}/s1_typeI.tsv", sep="\t")
S.columns = S.columns.str.strip(); T.columns = T.columns.str.strip()
for c in S.columns:
    if S[c].dtype == object: S[c] = S[c].str.strip()
for c in T.columns:
    if T[c].dtype == object: T[c] = T[c].str.strip()
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
d = S[(S.ld_r == 0.6)]
col = {"gls": "#1f6f8b", "ols": "#7fb3c8", "naive": "#d1495b", "oracle": "#888"}
for tr, ls in [("indep", "-"), ("clustered", "--")]:
    x = d[(d.traits == tr) & (d.rho == 0.5)]
    for m in ["gls", "ols", "naive"]:
        ax[0].plot(x.q, x[f"{m}_bias"], ls, marker="o", color=col[m], label=f"{m.upper() if m!='naive' else 'naive r'} ({tr})")
    ax[1].plot(x.q, x.gls_sd, ls, marker="o", color=col["gls"], label=f"GLS SD ({tr})")
    ax[1].plot(x.q, x.gls_pse, ls, marker="x", color="k", label=f"GLS mean SE ({tr})")
    ax[1].plot(x.q, x.oracle_sd, ls, marker="s", color=col["oracle"], label=f"oracle SD ({tr})")
    t = T[(T.traits == tr) & (T.ld_r == 0.6)]
    ax[2].plot(t.q, t.gls, ls, marker="o", color=col["gls"], label=f"LRCP-GLS ({tr})")
    ax[2].plot(t.q, t.naive, ls, marker="o", color=col["naive"], label=f"naive Fisher ({tr})")
ax[0].axhline(0, color="k", lw=.5); ax[0].set(xscale="log", xlabel="traits q", ylabel="bias of ρ̂ (true ρ = 0.5)", title="a  Bias")
ax[1].set(xscale="log", xlabel="traits q", ylabel="SD / SE", title="b  Precision and SE calibration (ρ = 0.5)")
ax[2].axhline(.05, color="k", lw=.5); ax[2].set(xscale="log", xlabel="traits q", ylabel="type-I error at 5%", title="c  Null distal pairs")
for a in ax:
    a.legend(fontsize=6.5); a.set_xticks([30, 100, 300]); a.set_xticklabels(["30", "100", "300"]); a.minorticks_off()
plt.tight_layout(); plt.savefig(f"{D}/fig2_s1.png", dpi=150)
