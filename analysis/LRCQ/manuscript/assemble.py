#!/usr/bin/env python3
"""Assemble papers/LRCQ/LRCQ-manuscript.md from the section drafts and add a
reference list built from notes/literature.md (entries cited in the text).

Usage: assemble.py [papers/LRCQ dir]
"""
import re
import sys
from pathlib import Path

root = Path(sys.argv[1] if len(sys.argv) > 1 else "/mnt/project-files/papers/LRCQ")
parts = ["01-abstract.md", "02-introduction.md", "03-results.md", "05-discussion.md", "04-methods.md", "06-figures.md"]
body = []
for p in parts:
    txt = (root / "draft" / p).read_text()
    txt = re.sub(r"^\*Draft[^\n]*\*\n\n?", "", txt, flags=re.M)   # drop per-section draft notes
    if p != "01-abstract.md":
        txt = re.sub(r"^(#+) ", r"#\1 ", txt, flags=re.M)   # demote every heading one level
    body.append(txt.strip())
text = "\n\n".join(body)

# reference list
lit = (root / "notes" / "literature.md").read_text().splitlines()
refs = {}
for line in lit:
    if "doi:" not in line or line.startswith("@"):
        continue
    entry = re.sub(r"^\s*-\s*(\*\*[^*]+\*\*\s*)?", "", line)
    entry = re.sub(r"\s*\*\*\[[^\]]*\]\*\*.*$", "", entry).strip()
    m = re.match(r"([^\s,;]+(?: [a-z]+)?)[^.]*\.\s.*?(\d{4})[;:]", entry)
    if not m:
        continue
    sur, yr = m.group(1).split()[0], m.group(2)
    if entry.startswith("1000 Genomes"):
        sur = "Genomes"
    refs.setdefault((sur, yr), [])
    if entry not in refs[(sur, yr)]:
        refs[(sur, yr)].append(entry)
cited = set()
for (sur, yr) in refs:
    if re.search(re.escape(sur) + r"[^()\[\]]{0,40}?" + yr, text):
        cited.add((sur, yr))
reflist = sorted(e for k in cited for e in refs[k])
# cited author-year strings with no entry
mentions = set(re.findall(r"([A-Z][A-Za-z'\-]+)(?: et al\.| & [A-Z][A-Za-z'\-]+)? (\d{4})[a-z]?", text))
missing = sorted(f"{a} {y}" for a, y in mentions if (a, y) not in refs and a != "Consortium" and 1990 < int(y) < 2030)
out = text + "\n\n## References\n\n" + "\n".join(f"{i}. {r}" for i, r in enumerate(reflist, 1)) + "\n"
hdr = ("<!-- Assembled by analysis/LRCQ/manuscript/assemble.py from draft/*.md; edit the section files, not this one. -->\n\n"
       "*Working draft. Numbers in [square brackets] are placeholders for analyses still running.*\n\n")
(root / "LRCQ-manuscript.md").write_text(hdr + out)
print(f"{len(reflist)} references; citations without an entry: {missing}")
