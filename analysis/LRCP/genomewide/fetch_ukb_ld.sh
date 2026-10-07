#!/bin/bash
# Robustness (step 14): UKB in-sample LD (PolyFun release) among the candidate tags of
# the windows holding the top genome-wide pairs; one 1 GB window at a time, deleted after use.
# usage: fetch_ukb_ld.sh <windows.txt> <universe.tsv> <tmp_dir> <out_dir>
W=$1; U=$2; T=$3; O=$4; S=https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD
DIR=$(dirname "$0")
while read -r w; do
  [ -s "$O/w$w.R.f32" ] && continue
  for ext in gz npz; do
    for i in 1 2 3 4; do curl -sf -o "$T/$w.$ext" "$S/$w.$ext" && break; sleep $((2**i)); done
  done
  python3 "$DIR/ukb_ld_subset.py" "$T/$w" "$U" "$O/w$w" || echo "FAILED $w"
  rm -f "$T/$w.gz" "$T/$w.npz"
done < "$W"
echo DONE
