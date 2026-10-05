"""Fetch Pan-UKB EUR summary statistics for a set of small genomic regions across
many phenotypes from the public S3 bucket.

Uses the tabix linear index to issue one small HTTP Range request per region
(htslib's remote reader pulled far more data than needed). Regions: the
TP53 / MDM2 motivating SNPs, a few positive-control loci, and random 20-kb
windows that give an empirical null for distal SNP pairs. Output (outside
the shared folder): one TSV per phenotype with region, chr, pos, ref, alt,
af, beta_EUR, se_EUR.
"""
import os, sys, io, gzip, struct, random, zlib, time, urllib.request
from multiprocessing import Pool
import pandas as pd

MANIFEST = "/mnt/project-files/Methodology papers/UKBB data/phenotype_manifest.tsv.bgz"
BASE = "https://pan-ukb-us-east-1.s3.amazonaws.com/"
OUT = sys.argv[1] if len(sys.argv) > 1 else "/home/user/data/panukb_regions"

CHRLEN = {1:249250621,2:243199373,3:198022430,4:191154276,5:180915260,6:171115067,7:159138663,
          8:146364022,9:141213431,10:135534747,11:135006516,12:133851895,13:115169878,14:107349540,
          15:102531392,16:90354753,17:81195210,18:78077248,19:59128983,20:63025520,21:48129895,22:51304566}

TARGETS = [  # name, chr, pos (GRCh37)
    ("TP53_rs78378222", 17, 7571752),
    ("MDM2_rs3730556", 12, 69216521),
    ("LDLR_rs6511720", 19, 11202306),
    ("PCSK9_rs11591147", 1, 55505647),
    ("HMGCR_rs12916", 5, 74656539),
    ("APOB_rs1367117", 2, 21263900),
]

def regions(n_random=150, half=10000, seed=1):
    rng = random.Random(seed)
    regs = [(name, c, p - half, p + half) for name, c, p in TARGETS]
    tot = sum(CHRLEN.values())
    while len(regs) < len(TARGETS) + n_random:
        x = rng.randrange(tot)
        for c, L in CHRLEN.items():
            if x < L: break
            x -= L
        if c == 6 and 25_000_000 < x < 34_000_000:  # skip MHC
            continue
        if x < 1_000_000 or x > CHRLEN[c] - 1_000_000:
            continue
        regs.append((f"rand{len(regs)-len(TARGETS):03d}", c, x - half, x + half))
    return regs

def get(url, rng=None, tries=5):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={rng[0]}-{rng[1]}"} if rng else {})
            return urllib.request.urlopen(req, timeout=120).read()
        except Exception:
            if i == tries - 1: raise
            time.sleep(2 ** i)

def read_tbi(raw):
    """Return {seqname: linear_index(list of virtual offsets)}."""
    b = gzip.decompress(raw)
    assert b[:4] == b"TBI\x01"
    n_ref, = struct.unpack_from("<i", b, 4)
    l_nm, = struct.unpack_from("<i", b, 32)
    names = b[36:36 + l_nm].split(b"\x00")[:n_ref]
    off = 36 + l_nm
    lin = {}
    for r in range(n_ref):
        n_bin, = struct.unpack_from("<i", b, off); off += 4
        for _ in range(n_bin):
            _, n_chunk = struct.unpack_from("<Ii", b, off); off += 8 + 16 * n_chunk
        n_intv, = struct.unpack_from("<i", b, off); off += 4
        lin[names[r].decode()] = struct.unpack_from(f"<{n_intv}Q", b, off); off += 8 * n_intv
    return lin

def bgzf_text(raw):
    """Decompress a run of whole BGZF blocks (drop a trailing partial block)."""
    out, pos = [], 0
    while pos + 18 <= len(raw):
        bsize = struct.unpack_from("<H", raw, pos + 16)[0] + 1
        if pos + bsize > len(raw): break
        out.append(zlib.decompress(raw[pos:pos + bsize], 31)); pos += bsize
    return b"".join(out)

def fetch_region(url, lin, c, s, e):
    li = lin.get(str(c))
    if not li: return []
    i0, i1 = s >> 14, (e >> 14) + 1
    if i0 >= len(li): return []
    v0 = li[i0]
    nxt = [v for v in li[i1:] if v > v0]
    c0 = v0 >> 16; u0 = v0 & 0xFFFF
    c1 = (nxt[0] >> 16) + 70000 if nxt else c0 + 2_000_000
    txt = bgzf_text(get(url, (c0, c1)))[u0:]
    rows = []
    for line in txt.split(b"\n"):
        f = line.split(b"\t")
        if len(f) < 3 or f[0] != str(c).encode(): continue
        p = int(f[1])
        if p > e: break
        if p >= s: rows.append(line.decode())
    return rows

def work(args):
    fname, regs = args
    out = os.path.join(OUT, fname.replace(".tsv.bgz", ".tsv"))
    if os.path.exists(out):
        return fname, "cached"
    url = BASE + "sumstats_flat_files/" + fname
    try:
        cols = zlib.decompressobj(31).decompress(get(url, (0, 65535))).split(b"\n")[0].decode().split("\t")
        af = "af_EUR" if "af_EUR" in cols else "af_controls_EUR"  # binary traits report case/control AF
        ix = [cols.index(k) for k in ("chr", "pos", "ref", "alt", af, "beta_EUR", "se_EUR")]
        lin = read_tbi(get(BASE + "sumstats_flat_files_tabix/" + fname + ".tbi"))
        rows = []
        for name, c, s, e in regs:
            for line in fetch_region(url, lin, c, max(1, s), e):
                f = line.split("\t")
                rows.append([name] + [f[i] for i in ix])
        tmp = out + ".tmp"
        pd.DataFrame(rows, columns=["region", "chr", "pos", "ref", "alt", "af", "beta", "se"]).to_csv(tmp, sep="\t", index=False)
        os.rename(tmp, out)
        return fname, len(rows)
    except Exception as ex:
        return fname, "FAILED " + repr(ex)

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    man = pd.read_csv(MANIFEST, sep="\t", compression="gzip", low_memory=False)
    man = man[man.phenotype_qc_EUR == "PASS"]
    man.to_csv(os.path.join(OUT, "traits.tsv"), sep="\t", index=False)
    regs = regions()
    pd.DataFrame(regs, columns=["region", "chr", "start", "end"]).to_csv(os.path.join(OUT, "regions.tsv"), sep="\t", index=False)
    files = list(man.filename) if len(sys.argv) < 3 else sys.argv[2:]
    with Pool(int(os.environ.get("NPROC", 16))) as pool:
        for i, (f, r) in enumerate(pool.imap_unordered(work, [(f, regs) for f in files])):
            print(i, f, r, flush=True)
