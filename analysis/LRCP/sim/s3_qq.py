"""S3: QQ plots of null distal pairs from the S1 grid (about 116k null pairs per cell)."""
import sys, subprocess, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
rds, outdir = sys.argv[1], sys.argv[2]
tsv = "/tmp/claude-0/s1_nulls.tsv"
subprocess.run(["Rscript", "-e", f'x<-readRDS("{rds}")$nulls; write.table(x[x$ld_r==0.6,],"{tsv}",sep="\\t",quote=F,row.names=F)'], check=True)
N = pd.read_csv(tsv, sep="\t")
fig, ax = plt.subplots(1, 2, figsize=(9, 4.3)); rows = []
for k, tr in enumerate(["indep", "clustered"]):
    for q, c in zip([30, 100, 300], ["#9ecae1", "#4292c6", "#08519c"]):
        d = N[(N.traits == tr) & (N.q == q)]
        for m, ls, cc in [("z_gls", "-", c), ("z_naive", ":", "#d1495b" if q == 300 else "#f4a3a8")]:
            p = np.sort(2 * stats.norm.sf(np.abs(d[m].values))); e = (np.arange(len(p)) + 0.5) / len(p)
            sel = np.unique(np.r_[np.arange(0, min(2000, len(p))), np.linspace(0, len(p) - 1, 3000).astype(int)])
            ax[k].plot(-np.log10(e[sel]), -np.log10(p[sel]), ls, color=cc, label=f"{'LRCP-GLS' if m=='z_gls' else 'naive'} q={q}")
            chi = d[m].values ** 2
            rows.append(dict(traits=tr, q=q, method=m, n=len(p), lambda_gc=np.median(chi) / stats.chi2.ppf(.5, 1),
                             p_lt_1e3=np.mean(p < 1e-3), p_lt_1e4=np.mean(p < 1e-4)))
    lim = 5.5; ax[k].plot([0, lim], [0, lim], "k", lw=.5); ax[k].set(xlim=(0, lim), ylim=(0, 12), title=f"{'abcd'[k]}  {tr} traits",
        xlabel="expected −log10 P", ylabel="observed −log10 P"); ax[k].legend(fontsize=6.5)
plt.tight_layout(); plt.savefig(f"{outdir}/figS_s3_qq.png", dpi=150)
pd.DataFrame(rows).to_csv(f"{outdir}/s3_null_qq_summary.tsv", sep="\t", index=False, float_format="%.4g")
print(pd.DataFrame(rows).to_string())
