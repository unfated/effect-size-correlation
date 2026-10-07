# LRCP paper analyses

Analysis scripts for the LRCP paper ("Identifying functionally related SNPs and genes from phenome-wide genetic-effect correlations"). Plans, notes, drafts and result tables are in the project's shared folder (`papers/LRCP/`). Estimators come from the `lrcpq` package in this repository.

| Folder | What |
|---|---|
| `motivating/` | TP53/MDM2 PheWAS parsing; Pan-UKB region fetcher (tabix linear index + HTTP range requests) for the motivating example and its empirical null |
| `naive_null/` | Simulation of the null behaviour of naive phenome-wide SNP-pair correlations (effective number of traits; bias of unsigned statistics) |
| `sim/` | LRCP simulations; `proto_distal.py` is a numpy prototype of the distal-pair GLS estimator used to cross-check the package |
