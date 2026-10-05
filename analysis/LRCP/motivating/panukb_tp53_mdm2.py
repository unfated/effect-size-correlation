"""TP53 rs78378222 x MDM2 rs3730556 across 452 Pan-UKB EUR traits: scatter,
per-category correlations, leave-one-category-out, and trait contributions."""
import sys, glob, os, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
IN, OUT = sys.argv[1], sys.argv[2]
tr = pd.read_csv(f"{IN}/traits.tsv", sep="\t", low_memory=False)
tr["key"] = tr.filename.str.replace(".tsv.bgz", "", regex=False)
rows = []
for k in tr.key:
    f = f"{IN}/{k}.tsv"
    if not os.path.exists(f): continue
    d = pd.read_csv(f, sep="\t", dtype={"chr": str})
    a = d[(d.chr == "17") & (d.pos == 7571752) & (d.alt == "G")]
    b = d[(d.chr == "12") & (d.pos == 69216521) & (d.alt == "G")]
    if len(a) and len(b):
        rows.append((k, a.beta.iloc[0] / a.se.iloc[0], b.beta.iloc[0] / b.se.iloc[0]))
m = pd.DataFrame(rows, columns=["key", "z_tp53", "z_mdm2"]).merge(
    tr[["key", "description", "trait_type", "category", "in_max_independent_set"]], on="key")
m = m.dropna(subset=["z_tp53", "z_mdm2"])
def dom(row):
    c = str(row.category); d = str(row.description).lower(); t = row.trait_type
    if t == "biomarkers" or "blood count" in c.lower() or "haematology" in c.lower(): return "blood/biomarker"
    if any(s in d for s in ("body mass", "weight", "fat", "height", "impedance", "waist", "hip", "sitting")): return "anthropometric"
    if any(s in d for s in ("blood pressure", "pulse", "heart", "cardi")): return "cardiovascular"
    if t in ("phecode", "icd10"): return "disease"
    if t == "prescriptions": return "medication"
    if "diet" in c.lower() or "food" in d or "intake" in d: return "diet"
    return "other"
m["domain"] = m.apply(dom, axis=1)
m["contrib"] = (m.z_tp53 - m.z_tp53.mean()) * (m.z_mdm2 - m.z_mdm2.mean())
m.to_csv(f"{OUT}/panukb_tp53_mdm2_by_trait.tsv", sep="\t", index=False, float_format="%.4g")
r_all = np.corrcoef(m.z_tp53, m.z_mdm2)[0, 1]
res = [dict(subset="all", n=len(m), r=r_all)]
for g, d in m.groupby("domain"):
    res.append(dict(subset=f"only {g}", n=len(d), r=np.corrcoef(d.z_tp53, d.z_mdm2)[0, 1] if len(d) > 4 else np.nan))
    o = m[m.domain != g]; res.append(dict(subset=f"without {g}", n=len(o), r=np.corrcoef(o.z_tp53, o.z_mdm2)[0, 1]))
mm = m.sort_values("contrib", ascending=False)
for k in (5, 10, 20):
    o = mm.iloc[k:]; res.append(dict(subset=f"without top-{k} contributing traits", n=len(o), r=np.corrcoef(o.z_tp53, o.z_mdm2)[0, 1]))
R = pd.DataFrame(res); R.to_csv(f"{OUT}/panukb_tp53_mdm2_subsets.tsv", sep="\t", index=False, float_format="%.4g")
print(R.to_string()); print(mm[["description", "domain", "z_tp53", "z_mdm2", "contrib"]].head(20).to_string())
print(m.domain.value_counts())
# figure
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw=dict(width_ratios=[1.2, 1]))
cmap = plt.get_cmap("tab10")
for i, (g, d) in enumerate(m.groupby("domain")):
    ax[0].scatter(d.z_tp53, d.z_mdm2, s=14, color=cmap(i), alpha=.8, label=f"{g} ({len(d)})")
ax[0].axhline(0, c="grey", lw=.5); ax[0].axvline(0, c="grey", lw=.5)
ax[0].set_xlabel("Z, rs78378222 (TP53)"); ax[0].set_ylabel("Z, rs3730556 (MDM2)")
ax[0].set_title(f"c  Pan-UKB EUR, {len(m)} traits: r = {r_all:.2f}"); ax[0].legend(fontsize=7, frameon=False)
s = R[R.subset.str.startswith("without")].sort_values("r")
ax[1].barh(s.subset, s.r, color="#55A868"); ax[1].axvline(r_all, c="k", ls="--", lw=.8)
ax[1].set_xlabel("r after removing traits"); ax[1].set_title("d  robustness"); ax[1].tick_params(axis="y", labelsize=7)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1cd_panukb_tp53_mdm2.png", dpi=200)
