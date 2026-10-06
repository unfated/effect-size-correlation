#!/usr/bin/env bash
# Extract Z for every trait in trait-list.tsv in parallel (NPAR, default 3).
# Pan-UKB's maximally independent traits go first so a 150-trait analysis can
# start before the full set is done. Restartable: finished traits are skipped.
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
list=${1:-/mnt/project-files/papers/LRCQ/results/real/trait-list.tsv}
uni=${2:-/home/user/data/ref/snp_universe.txt}
out=${3:-/home/user/data/z}
mkdir -p "$out"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}{print ($c["indep"]=="True"?0:1)"\t"NR"\t"$c["trait_id"]"\t"$c["aws_path"]}' "$list" |
  sort -k1,1n -k2,2n | cut -f3,4 |
  xargs -P ${NPAR:-3} -L 1 bash -c '"$0" "$1" "$3" "$4" "$2"' "$here/02_extract_z.sh" "$uni" "$out" 2>&1
