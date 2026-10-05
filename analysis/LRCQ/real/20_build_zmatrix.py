"""Assemble per-trait Z files into one SNP x trait float32 matrix.

Usage: 20_build_zmatrix.py <trait-list.tsv> <z_dir> <out_prefix> [--indep-only]
Writes <out_prefix>.Z.npy and <out_prefix>.Z.f32 (row-major, m x q, NaN = missing),
<out_prefix>.traits.tsv (trait rows in column order) and prints coverage.
Only traits whose Z file exists are included.
"""
import sys, gzip, os
import numpy as np
import pandas as pd

lst, zdir, out = sys.argv[1:4]
indep_only = "--indep-only" in sys.argv
t = pd.read_csv(lst, sep="\t")
if indep_only:
    t = t[t.indep]
t = t[[os.path.exists(os.path.join(zdir, x.replace("/", "_") + ".z.gz")) for x in t.trait_id]].reset_index(drop=True)
m = sum(1 for _ in open("/home/user/data/ref/snp_universe.txt"))
Z = np.lib.format.open_memmap(out + ".Z.npy", mode="w+", dtype=np.float32, shape=(m, len(t)))
for j, tid in enumerate(t.trait_id):
    v = pd.read_csv(os.path.join(zdir, tid.replace("/", "_") + ".z.gz"), header=None, na_values="NA",
                    dtype=np.float32).values[:, 0]
    assert len(v) == m, tid
    Z[:, j] = v
Z.flush()
# raw row-major float32 copy for readers without .npy support (R readBin + seek)
np.asarray(Z).tofile(out + ".Z.f32")
t.to_csv(out + ".traits.tsv", sep="\t", index=False)
miss = np.isnan(Z).mean(0)
print(f"{len(t)} traits x {m} SNPs; missing per trait: median {np.median(miss):.4f}, max {miss.max():.4f}")
