#!/usr/bin/env bash
# Stream one Pan-UKB flat sumstats file and keep EUR Z-scores for the SNP universe.
# Usage: 02_extract_z.sh <snp_universe.txt> <trait_id> <s3_path> <out_dir>
# snp_universe.txt: one "chr:pos:ref:alt" per line (Pan-UKB EUR LD-score SNPs, HapMap3).
# Output: <out_dir>/<trait_id>.z.gz, one line per universe SNP, "NA" when missing,
# low-confidence in EUR, or EUR AF (controls AF for binary traits) outside [0.01, 0.99].
set -euo pipefail
uni=$1; tid=$2; s3=$3; out=$4
url="https://pan-ukb-us-east-1.s3.amazonaws.com/${s3#s3://pan-ukb-us-east-1/}"
dest="$out/$tid.z.gz"
[ -s "$dest" ] && exit 0
# Column positions differ between quantitative (af_EUR) and binary
# (af_cases_EUR / af_controls_EUR) traits, so resolve them from the header.
hdr=""
for attempt in 1 2 3 4; do
  hdr=$(curl -sSf -r 0-65535 "$url" | zcat 2>/dev/null | head -1 || true)
  [ -n "$hdr" ] && break; sleep $((2**attempt))
done
fields=$(echo "$hdr" | tr '\t' '\n' | awk '{c[$1]=NR} END{af=("af_EUR" in c)? c["af_EUR"] : c["af_controls_EUR"];
  if(!af || !c["beta_EUR"] || !c["se_EUR"] || !c["low_confidence_EUR"]) exit 1;
  print "1,2,3,4,"af","c["beta_EUR"]","c["se_EUR"]","c["low_confidence_EUR"]}') || { echo "BADHEADER $tid" >&2; exit 1; }
for attempt in 1 2 3 4; do
  if curl -sSf --retry 3 "$url" | zcat | cut -f"$fields" | ${AWK:-mawk} -v U="$uni" '
      BEGIN{FS=OFS="\t"; while((getline l < U)>0){n++; idx[l]=n}}
      NR==1{next}
      { k=$1":"$2":"$3":"$4; if(k in idx){ if($6!="NA" && $7!="NA" && $7>0 && $8=="false" && $5!="NA" && $5>=0.01 && $5<=0.99) z[idx[k]]=$6/$7 } }
      END{for(i=1;i<=n;i++) print ((i in z)? sprintf("%.5g", z[i]) : "NA")}' | gzip > "$dest.tmp"; then
    mv "$dest.tmp" "$dest"; echo "$(date +%T) done $tid"; exit 0
  fi
  echo "$(date +%T) retry $attempt $tid" >&2; sleep $((2**attempt))
done
echo "FAILED $tid" >&2; exit 1
