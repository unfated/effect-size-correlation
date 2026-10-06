#!/usr/bin/env bash
# Genome-wide LRCQ on Pan-UKB: Z matrix -> matrix LDSC -> stage 3 (sharded) ->
# post-processing, annotation regression, motivation figure.
# Usage: run_real.sh [label] [method] [trait_ids_file|all]
# Env: D (data root, default /home/user/data), NSHARD (default 4), RES (shared results dir)
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
D=${D:-/home/user/data}
lab=${1:-all}; method=${2:-ols}; subset=${3:-all}
NSHARD=${NSHARD:-4}
RES=${RES:-/mnt/project-files/papers/LRCQ/results/real}
lst=/mnt/project-files/papers/LRCQ/results/real/trait-list.tsv
zp=$D/real/panukb; lp=$D/real/panukb_ldsc
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p $D/real $RES

if [ ! -f $zp.Z.f32 ]; then
  python3 $here/20_build_zmatrix.py $lst $D/z $zp
fi
if [ ! -f $lp.h2.tsv ]; then
  OPENBLAS_NUM_THREADS=4 python3 $here/21_ldsc_matrix.py $zp $D/ref/UKBB.EUR.l2.ldscore.gz $(cat $D/ref/UKBB.EUR.l2.M_5_50) $lp
fi

# stage 1-2 outputs shared with the LRCP thread
for f in h2 intercept intercept_se gcov gcov_se rg; do
  [ -f $RES/panukb_ldsc.$f.tsv ] || cp $lp.$f.tsv $RES/panukb_ldsc.$f.tsv
done
[ -f $RES/panukb.traits.tsv ] || cp $zp.traits.tsv $RES/panukb.traits.tsv

# stage 3, sharded by chromosome groups balanced by SNP count
out=$D/real/s3_${lab}_${method}
mkdir -p $out
for k in $(seq 1 $NSHARD); do mkdir -p $out/ld$k; done
python3 - "$D/ld/kg_hm3" "$out" "$NSHARD" <<'EOF'
import os, sys, glob, re
ld, out, ns = sys.argv[1], sys.argv[2], int(sys.argv[3])
files = sorted(glob.glob(os.path.join(ld, "*.snps.tsv")))
by = {}
for f in files:
    c = int(re.search(r"chr(\d+)_", os.path.basename(f)).group(1))
    by.setdefault(c, []).append(f)
size = {c: sum(sum(1 for _ in open(f)) for f in fs) for c, fs in by.items()}
load = [0] * ns; asg = {}
for c in sorted(size, key=lambda c: -size[c]):
    k = load.index(min(load)); asg[c] = k; load[k] += size[c]
for c, fs in by.items():
    for f in fs:
        for ext in (".snps.tsv", ".R.f32"):
            src = f[:-len(".snps.tsv")] + ext
            dst = os.path.join(out, f"ld{asg[c]+1}", os.path.basename(src))
            if not os.path.exists(dst): os.symlink(src, dst)
print("shard loads:", load)
EOF
pids=()
for k in $(seq 1 $NSHARD); do
  nice -n 5 Rscript $here/30_stage3_lrcq.R $zp $lp $out/ld$k $out/shard$k.tsv $method 503 $subset 20 > $out/shard$k.log 2>&1 &
  pids+=($!)
done
for p in "${pids[@]}"; do wait $p; done

# merge shards (tsv rows and jackknife rows stay aligned)
python3 - "$out" "$NSHARD" "$lab" "$method" <<'EOF'
import sys, os, pandas as pd, numpy as np
out, ns, lab, meth = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
t = [pd.read_csv(f"{out}/shard{k}.tsv", sep="\t") for k in range(1, ns + 1)]
m = [pd.read_csv(f"{out}/shard{k}.snp2tag.tsv", sep="\t") for k in range(1, ns + 1)]
pre = f"{out}/lrcq_{lab}_{meth}"
pd.concat(t).to_csv(pre + ".tsv", sep="\t", index=False)
pd.concat(m).to_csv(pre + ".snp2tag.tsv", sep="\t", index=False)
if all(os.path.exists(f"{out}/shard{k}.jk.f32") for k in range(1, ns + 1)):
    with open(pre + ".jk.f32", "wb") as fo:
        for k in range(1, ns + 1):
            fo.write(open(f"{out}/shard{k}.jk.f32", "rb").read())
os.replace(f"{out}/shard1.clusters.tsv", pre + ".clusters.tsv")
print("merged", sum(len(x) for x in t), "tags")
EOF

# annotation LD scores + pooled response for this trait set
an=$D/real/annot_${lab}
[ -f $an.npz ] || nice -n 5 python3 $here/41_annot_ld.py $D/ld/kg_hm3 $D/annot/baselineLF_v2.2.UKB $zp $lp 503 $an $([ "$subset" != all ] && echo $subset) > $an.log 2>&1
python3 $here/40_postprocess.py $out/lrcq_${lab}_${method} $RES $an.npz
python3 $here/42_annot_regression.py $an.npz $RES/annot_regression_${lab}.tsv none > /dev/null
python3 $here/42_annot_regression.py $an.npz $RES/annot_regression_${lab}_l2w.tsv l2 > /dev/null
cp $out/lrcq_${lab}_${method}.clusters.tsv $RES/
echo "done $lab $method"
