"""Genome-wide LRCP modules (work plan step 13): connected components of the tag-pair
graph at screen P < P_EDGE, annotated with nearest genes, profile class and q_eff (from
confirm.R), and the traits that load most on each module's sign-aligned mean profile.
usage: modules.py <out_dir> <candidate_tags.tsv> <tagZ.tsv.gz> <traits.tsv> <p_edge>"""
import sys, numpy as np, pandas as pd
OUT, CAND, ZF, TR, PE = sys.argv[1:6]; PE = float(PE)
sc = pd.read_csv(f"{OUT}/screen_pairs_p1e-3.tsv.gz", sep="\t")
cf = pd.read_csv(f"{OUT}/confirm_top.tsv", sep="\t")
e = sc[sc.p < PE]
genes = pd.read_csv("/mnt/project-files/papers/LRCQ/results/real/genes_all.genes.tsv.gz", sep="\t")
def near(i):
    c, b = int(i.split(":")[0]), int(i.split(":")[1]); g = genes[genes.chr.astype(str) == str(c)]
    d = np.maximum(np.maximum(g.start - b, b - g.end), 0); return g.gene.values[np.argmin(d)]
# union-find
par = {}
def f(x):
    par.setdefault(x, x)
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for a, b in zip(e.ID1, e.ID2): par[f(a)] = f(b)
nodes = sorted(set(e.ID1) | set(e.ID2)); comp = {n: f(n) for n in nodes}
Z = pd.read_csv(ZF, sep="\t", index_col=0).fillna(0.0); tr = pd.read_csv(TR, sep="\t")
desc = dict(zip(tr.trait_id, tr.description))
rows, mods = [], []
for k, (root, mem) in enumerate(sorted(pd.Series(comp).groupby(lambda x: comp[x]).groups.items(), key=lambda t: -len(t[1]))):
    mem = list(mem); ee = e[e.ID1.isin(mem) & e.ID2.isin(mem)]
    # sign-align members to the first by the sign of rho along a spanning order
    sign = {mem[0]: 1.0}; changed = True
    while changed:
        changed = False
        for a, b, r in zip(ee.ID1, ee.ID2, ee.rho):
            if a in sign and b not in sign: sign[b] = sign[a] * np.sign(r); changed = True
            if b in sign and a not in sign: sign[a] = sign[b] * np.sign(r); changed = True
    X = np.vstack([Z.loc[m].values * sign[m] / np.linalg.norm(Z.loc[m].values) for m in mem]).mean(0)
    top = np.argsort(-np.abs(X))[:6]
    traits = "; ".join(f"{desc[Z.columns[t]]} ({'+' if X[t] > 0 else '-'})" for t in top)
    cls = cf.set_index("ID1").class1.to_dict(); cls.update(cf.set_index("ID2").class2.to_dict())
    qe = cf.set_index("ID1").q_eff_class1.to_dict(); qe.update(cf.set_index("ID2").q_eff_class2.to_dict())
    gs = [near(m) for m in mem]
    mods.append(dict(module=k + 1, n_loci=len(mem), n_edges=len(ee), genes=",".join(gs),
                     min_q_eff_class=np.nanmin([qe.get(m, np.nan) for m in mem]) if any(m in qe for m in mem) else np.nan,
                     top_traits=traits))
    for m, g in zip(mem, gs): rows.append(dict(module=k + 1, ID=m, gene=g, sign=sign.get(m, np.nan), cls=cls.get(m, np.nan)))
M = pd.DataFrame(mods); M.to_csv(f"{OUT}/modules_p{PE:g}.tsv", sep="\t", index=False)
pd.DataFrame(rows).to_csv(f"{OUT}/module_members_p{PE:g}.tsv", sep="\t", index=False)
print("edges", len(e), "nodes", len(nodes), "modules", len(M))
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print(M.head(15).to_string())
