"""Pathway validation of genome-wide LRCP links (work plan step 13).
Tag pairs are mapped to nearest genes (LRCQ gene table); a pair is 'co-annotated' if both genes share
at least one gene set of size 10-500 in a collection (MSigDB v7.5.1 Reactome, KEGG, GO:BP; via
msigdbr's GitHub data). Observed co-annotation of the LRCP edges is compared with
 (a) random distal pairs of LRCQ candidate tags (same tag universe),
 (b) degree-preserving rewiring of the LRCP edges (same nodes, controls for annotation-rich genes),
 (c) the same number of top pairs by naive |Z correlation| (does LRCP add over naive similarity?).
usage: pathway_enrichment.py <genomewide_dir> <tagZ.tsv.gz> <genesets.tsv.gz> <p_edge> <out.tsv>"""
import sys, numpy as np, pandas as pd
D, ZF, GS, PE, OUT = sys.argv[1:6]; PE = float(PE); rng = np.random.default_rng(13)
cand = pd.read_csv(f"{D}/candidate_tags_q05.tsv", sep="\t")
genes = pd.read_csv("/mnt/project-files/papers/LRCQ/results/real/genes_all.genes.tsv.gz", sep="\t")
genes["chr"] = genes.chr.astype(str)
near = {}
for c, cc in cand.groupby("CHR"):
    g = genes[genes.chr == str(c)]; s = g.start.values; e = g.end.values
    for i, b in zip(cc.ID, cc.BP): near[i] = g.gene.values[np.argmin(np.maximum(np.maximum(s - b, b - e), 0))]
cand["gene"] = cand.ID.map(near)
gs = pd.read_csv(GS, sep="\t")
sz = gs.groupby("gs_id").symbol.transform("size"); gs = gs[(sz >= 10) & (sz <= 500)]
coll = {k: {g: set(v.gs_id) for g, v in d.groupby("symbol")} for k, d in gs.groupby("gs_subcat")}
def coann(a, b, k):
    A = coll[k].get(a); B = coll[k].get(b); return bool(A and B and (A & B))
sc = pd.read_csv(f"{D}/screen_pairs_p1e-3.tsv.gz", sep="\t")
e = sc[sc.p < PE].copy(); e["g1"] = e.ID1.map(near); e["g2"] = e.ID2.map(near); e = e[e.g1 != e.g2]
pairs = list(zip(e.g1, e.g2)); n = len(pairs)
chrom = dict(zip(cand.ID, cand.CHR)); pos = dict(zip(cand.ID, cand.BP))
def distal(a, b): return chrom[a] != chrom[b] or abs(pos[a] - pos[b]) >= 5e6
ids = cand.ID.values
def rand_pairs(k):
    out = []
    while len(out) < k:
        a, b = rng.choice(ids, 2, replace=False)
        if distal(a, b) and near[a] != near[b]: out.append((near[a], near[b]))
    return out
def rewire(pr, nsw=2000):
    pr = [list(p) for p in pr]
    for _ in range(nsw):
        i, j = rng.choice(len(pr), 2, replace=False)
        a, b = pr[i]; c, d = pr[j]
        if len({a, b, c, d}) == 4: pr[i], pr[j] = [a, d], [c, b]
    return [tuple(p) for p in pr]
# naive comparison: top-n distal pairs by |naive r| among the same candidate tags
Z = pd.read_csv(ZF, sep="\t", index_col=0).fillna(0.0).values
Zs = (Z - Z.mean(1, keepdims=True)) / Z.std(1, keepdims=True); Rn = Zs @ Zs.T / Z.shape[1]
iu = np.triu_indices(len(ids), 1); ch = cand.CHR.values; bp = cand.BP.values
dm = (ch[iu[0]] != ch[iu[1]]) | (np.abs(bp[iu[0]] - bp[iu[1]]) >= 5e6)
r = np.abs(Rn[iu]); r[~dm] = -1; top = np.argsort(-r)
naive = []
for t in top:
    a, b = ids[iu[0][t]], ids[iu[1][t]]
    if near[a] != near[b]: naive.append((near[a], near[b]))
    if len(naive) == n: break
rows = []
for k in coll:
    obs = np.mean([coann(a, b, k) for a, b in pairs])
    nr = [np.mean([coann(a, b, k) for a, b in rand_pairs(n)]) for _ in range(500)]
    nw = [np.mean([coann(a, b, k) for a, b in rewire(pairs)]) for _ in range(500)]
    nv = np.mean([coann(a, b, k) for a, b in naive])
    rows.append(dict(collection=k, n_pairs=n, obs=obs, random_mean=np.mean(nr), fold_random=obs / np.mean(nr),
                     p_random=(1 + np.sum(np.array(nr) >= obs)) / 501, rewired_mean=np.mean(nw),
                     p_rewired=(1 + np.sum(np.array(nw) >= obs)) / 501, naive_top=nv))
R = pd.DataFrame(rows); R.to_csv(OUT, sep="\t", index=False); print(R.round(4).to_string())
