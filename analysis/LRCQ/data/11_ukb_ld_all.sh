#!/usr/bin/env bash
# Download UKB in-sample LD windows (3 Mb, one every 2 Mb: files starting at 0, 2, 4, ... Mb) and
# keep only the LRCQ SNP universe. Consecutive windows overlap by 1 Mb; LRCQ uses
# the middle 2 Mb of each window ([start+0.5, start+2.5) Mb) as its core, so cores
# tile the genome with >= 0.5 Mb of LD context on each side.
# Restartable; raw .npz files are deleted after subsetting.
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
uni=${1:-/home/user/data/ref/snp_universe.tsv}
out=${2:-/home/user/data/ld/hm3}
tmp=${3:-/home/user/data/ld/tmp}
B=https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD
mkdir -p "$out" "$tmp"
list="$out/../window_list.txt"
if [ ! -s "$list" ]; then
  tok=""; : > "$list.all"
  while :; do
    url="https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/?list-type=2&prefix=UKBB_LD/chr&max-keys=1000"
    [ -n "$tok" ] && url="$url&continuation-token=$(python3 -c "import urllib.parse,sys;print(urllib.parse.quote(sys.argv[1]))" "$tok")"
    page=$(curl -sS "$url")
    echo "$page" | grep -o '<Key>UKBB_LD/chr[^<]*\.npz</Key>' | sed 's/<Key>UKBB_LD\///; s/\.npz<\/Key>//' >> "$list.all"
    tok=$(echo "$page" | grep -o '<NextContinuationToken>[^<]*' | sed 's/<NextContinuationToken>//')
    [ -z "$tok" ] && break
  done
  # keep windows starting at an even Mb (file start = 1 + k*2e6)
  awk -F_ '{s=$2; if(((s-1)/1000000)%2==0) print}' "$list.all" | sort -V > "$list"
fi
echo "$(wc -l < "$list") windows"
while read -r w; do
  [ -s "$out/$w.R.f32" ] && continue
  for attempt in 1 2 3 4; do
    curl -sSf -o "$tmp/$w.gz" "$B/$w.gz" && curl -sSf -o "$tmp/$w.npz" "$B/$w.npz" && break
    sleep $((2**attempt))
  done
  python3 "$here/10_ukb_ld_subset.py" "$tmp/$w" "$uni" "$out/$w" || echo "FAILED $w"
  rm -f "$tmp/$w.gz" "$tmp/$w.npz"
done < "$list"
