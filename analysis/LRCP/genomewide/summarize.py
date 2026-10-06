# Numbers quoted in the genome-wide Results and Robustness sections, from one set of screen
# and confirmation outputs. usage: summarize.py <main_dir> <indep_dir> <q01_dir> <top30_dir> <cand_tags.tsv>
import sys, numpy as np, pandas as pd
D, DI, DQ, DT, CAND = sys.argv[1:6]
rd = lambda p: pd.read_csv(p, sep="\t")
key = lambda d: d.ID1 + "|" + d.ID2
sm = rd(f"{D}/screen_summary.tsv").iloc[0]
sc = rd(f"{D}/screen_pairs_p1e-3.tsv.gz"); sc["k"] = key(sc)
cf = rd(f"{D}/confirm_top.tsv"); cf["k"] = key(cf)
cu = rd(f"{D}/confirm_top_ukbld.tsv"); cu["k"] = key(cu)
hits = sc[sc.q_BH < 0.05]
print(f"screen: FDR5% {int(sm.n_q05)}, Bonferroni {int(sm.n_bonf)}, lambda {sm.lambda_median:.2f}, P<0.05 {sm.frac_p05:.4f}, P<1e-3 {sm.frac_p1e3:.2e}")
print(f"pairs P<1e-5: {(sc.p < 1e-5).sum()}, P<1e-4: {(sc.p < 1e-4).sum()}")
print(f"top100: bootstrap P max {cf.p_boot.max():.5f}; boot SD mean {np.nanmean(np.r_[cf.boot_sd_A, cf.boot_sd_B]):.3f}")
nar = lambda d: (d.q_eff_class1 < 5) | (d.q_eff_class2 < 5)
both = (cf.q_eff_class1 < 5) & (cf.q_eff_class2 < 5)
qmin = np.minimum(cf.q_eff_class1, cf.q_eff_class2)
t = cf.testable.astype(bool)
print(f"class null: >=1 narrow {nar(cf).sum()} (q_eff {qmin[nar(cf)].min():.2f}-{qmin[nar(cf)].max():.2f}); both narrow {both.sum()} (max |z_class| {cf.z_class[both].abs().max():.2f}); testable {t.sum()} (class P max {cf.p_class[t].max():.2e}, below 5e-8: {(cf.p_class[t] < 5e-8).sum()})")
for g in [("GSDMC", "CX3CR1"), ("CX3CR1", "GSDMC")]:
    r = cf[(cf.gene1 == g[0]) & (cf.gene2 == g[1])]
    if len(r): print(f"  {g}: z_class {r.z_class.iloc[0]:.2f}")
print(f"naive |r| of FDR hits: median {hits.r_naive.abs().median():.2f}")
m = cf.merge(cu, on="k", suffixes=("", "_u"))
print(f"UKB LD: n {len(m)}, z_gene r {np.corrcoef(m.z_gene, m.z_gene_u)[0,1]:.3f}, signs agree {(np.sign(m.z_gene) == np.sign(m.z_gene_u)).mean():.3f}, median |z| ratio {np.median(m.z_gene_u.abs() / m.z_gene.abs()):.3f}, min |z| {m.z_gene_u.abs().min():.2f}, boot P max {cu.p_boot.max():.5f}")
for lab, d2 in [("indep150", DI), ("q01", DQ), ("top30", DT)]:
    s2 = rd(f"{d2}/screen_summary.tsv").iloc[0]; p2 = rd(f"{d2}/screen_pairs_p1e-3.tsv.gz"); p2["k"] = key(p2)
    b = sc.merge(p2, on="k", suffixes=("", "_2"))
    hk = set(hits.k); h2 = p2[p2.k.isin(hk)]
    print(f"{lab}: FDR {int(s2.n_q05)}, Bonf {int(s2.n_bonf)}, lambda {s2.lambda_median:.2f}; pairs P<1e-3 in both {len(b)}, z r {np.corrcoef(b.z, b.z_2)[0,1]:.3f}, rho r {np.corrcoef(b.rho, b.rho_2)[0,1]:.3f}, signs agree {(np.sign(b.z) == np.sign(b.z_2)).mean():.3f}; main hits with P<1e-3 here {len(h2)}, at FDR5% here {(h2.q_BH < 0.05).sum()}")
    if lab == "q01":
        c = rd(CAND); ok = set(c.ID[c.q_gt0 < 0.01]); e = hits[hits.ID1.isin(ok) & hits.ID2.isin(ok)]
        print(f"  main hits with both tags q<0.01: {len(e)}, of which FDR5% in q01 screen: {p2[p2.k.isin(set(e.k))].q_BH.lt(0.05).sum()}")
# main hits under the top30 panel: |z| from the full-pair screen is not stored below P 1e-3, so read profiles
print(f"FDR hits: {len(hits)}")
