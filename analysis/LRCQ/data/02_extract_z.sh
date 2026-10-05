#!/usr/bin/env bash
# Stream one Pan-UKB flat sumstats file and keep EUR Z-scores for the SNP universe.
# Usage: 02_extract_z.sh <snp_universe.txt> <trait_id> <s3_path> <out_dir>
# snp_universe.txt: one "chr:pos:ref:alt" per line (Pan-UKB EUR LD-score SNPs, HapMap3).
# Output: <out_dir>/<trait_id>.z.gz, one line per universe SNP, "NA" when missing,
# low-confidence in EUR, or EUR AF (controls AF for binary traits) outside [0.01, 0.99].
# Parsing is done by src/extract_z.c (columns resolved from the header, so
# quantitative and binary layouts both work); it is built on first use.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
uni=$1; tid=$2; s3=$3; out=$4
bin=${EXTRACT_Z_BIN:-$here/src/extract_z}
[ -x "$bin" ] || gcc -O2 -o "$bin" "$here/src/extract_z.c" -lz
url="https://pan-ukb-us-east-1.s3.amazonaws.com/${s3#s3://pan-ukb-us-east-1/}"
dest="$out/$tid.z.gz"
[ -s "$dest" ] && exit 0
for attempt in 1 2 3 4; do
  if curl -sSf --retry 3 "$url" | "$bin" "$uni" 2>/dev/null | gzip > "$dest.tmp"; then
    mv "$dest.tmp" "$dest"; echo "$(date +%T) done $tid"; exit 0
  fi
  echo "$(date +%T) retry $attempt $tid" >&2; sleep $((2**attempt))
done
echo "FAILED $tid" >&2; exit 1
