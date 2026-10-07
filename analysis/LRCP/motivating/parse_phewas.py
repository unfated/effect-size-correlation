"""Parse the Common Metabolic Diseases Knowledge Portal PheWAS exports in
`observed data hints/` into tidy tables (one row per SNP x phenotype, plus one
row per SNP x phenotype x dataset)."""
import json, re, sys, os
import pandas as pd

SRC = "/mnt/project-files/Methodology papers/observed data hints"
LINE = re.compile(r'^"(?P<var>[^"]+)","[^"]*",\d+,(?P<ph>\{.*?\}),(?P<p>[^,]+),(?P<b>[^,]+),"(?P<se>[^"]*)",(?P<n>[^,]+),(?P<ds>\[.*\])\s*$')

def parse(fn):
    meta, per = [], []
    with open(os.path.join(SRC, fn), encoding="utf-8-sig") as f:
        next(f)
        for line in f:
            m = LINE.match(line.strip())
            if not m:
                print("skip:", line[:80], file=sys.stderr); continue
            ph = json.loads(m["ph"])
            row = dict(var=m["var"], pheno=ph["name"], desc=ph["description"], group=ph["group"],
                       dichotomous=ph["dichotomous"], p=float(m["p"]), beta=float(m["b"]),
                       se=float(m["se"]) if m["se"] else float("nan"), n=float(m["n"]))
            meta.append(row)
            for d in json.loads(m["ds"]):
                per.append(dict(var=m["var"], pheno=ph["name"], group=ph["group"], dataset=d["dataset"],
                                p=d["pValue"], beta=d["beta"], n=d.get("n")))
    return pd.DataFrame(meta), pd.DataFrame(per)

if __name__ == "__main__":
    out = sys.argv[1]
    for fn, tag in [("PheWAS_associations SNP1.csv", "SNP1"), ("PheWAS_associations SNP2.csv", "SNP2")]:
        a, b = parse(fn)
        a.to_csv(f"{out}/phewas_{tag}_meta.tsv", sep="\t", index=False)
        b.to_csv(f"{out}/phewas_{tag}_datasets.tsv", sep="\t", index=False)
        print(tag, a.shape, b.shape, a["var"].unique(), a["p"].max())
