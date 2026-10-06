"""Select Pan-UKB EUR traits for the LRCQ real-data application.

Inputs: Pan-UKB phenotype manifest and h2 manifest (bgzipped TSVs).
Output: trait-list.tsv with one row per selected trait.

Rules (primary set): EUR phenotype QC == PASS (Pan-UKB's own QC: GWAS run,
reasonable N, defined and in-bounds h2, h2 z > 4, normal lambda and ratio),
both sexes only. The 'indep' flag marks Pan-UKB's maximally independent set
(pairwise phenotypic correlation < 0.1), used as a sensitivity set.
"""
import sys
import pandas as pd

pheno_path, h2_path, out = sys.argv[1:4]
key = ["trait_type", "phenocode", "pheno_sex", "coding", "modifier"]
p = pd.read_csv(pheno_path, sep="\t", compression="gzip", low_memory=False)
h = pd.read_csv(h2_path, sep="\t", compression="gzip", low_memory=False)
h = h[h["pop"] == "EUR"].copy()
for k in key:
    p[k] = p[k].astype(str)
    h[k] = h[k].astype(str)
m = p.merge(h, on=key, how="left", suffixes=("", "_h2"))
s = m[(m.phenotype_qc_EUR == "PASS") & (m.pheno_sex == "both_sexes")].copy()

s["n_cases"] = s.n_cases_EUR
s["n_controls"] = s.n_controls_EUR
s["N"] = s.n_cases_EUR + s.n_controls_EUR.fillna(0)
s["binary"] = s.n_controls_EUR.notna()
# effective N for binary traits (4 / (1/Ncase + 1/Ncontrol)), used in sensitivity only
s["N_eff"] = s.N
b = s.binary
s.loc[b, "N_eff"] = 4 / (1 / s.loc[b, "n_cases"] + 1 / s.loc[b, "n_controls"])
s["prevalence"] = s.n_cases / s.N
s["category_top"] = s.category.astype(str).str.split(" > ").str[0]
s["trait_id"] = s[key].fillna("").astype(str).agg("-".join, axis=1)

cols = ["trait_id"] + key + ["description", "category", "category_top", "binary",
        "n_cases", "n_controls", "N", "N_eff", "prevalence",
        "estimates.final.h2_observed", "estimates.final.h2_observed_se",
        "estimates.final.h2_liability", "estimates.final.h2_z",
        "estimates.ldsc.h2_observed", "estimates.ldsc.intercept",
        "estimates.ldsc.intercept_se", "estimates.ldsc.ratio", "lambda_gc_EUR",
        "in_max_independent_set", "aws_path", "size_in_bytes"]
s = s[cols].rename(columns={"in_max_independent_set": "indep",
                            "estimates.final.h2_observed": "h2_obs",
                            "estimates.final.h2_observed_se": "h2_obs_se",
                            "estimates.final.h2_liability": "h2_liab",
                            "estimates.final.h2_z": "h2_z",
                            "estimates.ldsc.h2_observed": "h2_ldsc_obs",
                            "estimates.ldsc.intercept": "ldsc_intercept",
                            "estimates.ldsc.intercept_se": "ldsc_intercept_se",
                            "estimates.ldsc.ratio": "ldsc_ratio"})
s = s.sort_values(["trait_type", "phenocode"]).reset_index(drop=True)
s.to_csv(out, sep="\t", index=False)
print(len(s), "traits;", int(s.indep.sum()), "in max independent set")
print(s.trait_type.value_counts().to_string())
print(s.category_top.value_counts().to_string())
