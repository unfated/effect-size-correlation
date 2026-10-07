"""Fetch Pan-UKB EUR Z scores for whole LD blocks (Berisa-Pickrell EUR) across
the 452 traits, keeping only HapMap3 SNPs in the Pan-UKB EUR LD-score file.
Reuses the range-request reader in ../motivating/fetch_panukb_regions.py.

usage: fetch_panukb_blocks.py <ldscore.gz> <out_dir> <name:chr:start:end> ...
Writes <out_dir>/<trait>.tsv with block, id (chr:pos:ref:alt), beta, se.
"""
import os, sys, gzip, zlib
from multiprocessing import Pool
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "motivating"))
import fetch_panukb_regions as F

ldscore, OUT = sys.argv[1], sys.argv[2]
BLOCKS = [(b.split(":")[0], int(b.split(":")[1]), int(b.split(":")[2]), int(b.split(":")[3])) for b in sys.argv[3:]]
ls = pd.read_csv(ldscore, sep="\t", compression="gzip", usecols=["CHR", "BP", "SNP"] if False else None)
ids = set()
for name, c, s, e in BLOCKS:
    x = ls[(ls.iloc[:, 0].astype(str) == str(c)) & (ls.iloc[:, 2] >= s) & (ls.iloc[:, 2] <= e)]
    ids |= set(x.iloc[:, 1].astype(str))   # column 2 = chr:pos:ref:alt id

def work(fname):
    out = os.path.join(OUT, fname.replace(".tsv.bgz", ".tsv"))
    if os.path.exists(out): return fname, "cached"
    url = F.BASE + "sumstats_flat_files/" + fname
    try:
        cols = zlib.decompressobj(31).decompress(F.get(url, (0, 65535))).split(b"\n")[0].decode().split("\t")
        ix = [cols.index(k) for k in ("chr", "pos", "ref", "alt", "beta_EUR", "se_EUR")]
        lin = F.read_tbi(F.get(F.BASE + "sumstats_flat_files_tabix/" + fname + ".tbi"))
        rows = []
        for name, c, s, e in BLOCKS:
            for line in F.fetch_region(url, lin, c, s, e):
                f = line.split("\t"); sid = ":".join(f[i] for i in ix[:4])
                if sid in ids: rows.append((name, sid, f[ix[4]], f[ix[5]]))
        pd.DataFrame(rows, columns=["block", "id", "beta", "se"]).to_csv(out + ".tmp", sep="\t", index=False)
        os.rename(out + ".tmp", out)
        return fname, len(rows)
    except Exception as ex:
        return fname, "FAILED " + repr(ex)

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    print("HM3 SNPs in blocks:", len(ids), flush=True)
    tr = pd.read_csv("/mnt/project-files/papers/LRCQ/results/real/trait-list.tsv", sep="\t", low_memory=False)
    files = [p.split("/")[-1] for p in tr.aws_path]
    with Pool(int(os.environ.get("NPROC", 12))) as pool:
        for i, (f, r) in enumerate(pool.imap_unordered(work, files)):
            print(i, f, r, flush=True)
