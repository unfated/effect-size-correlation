#!/usr/bin/env python3
"""Fig. 1a of the LRCP paper: how two unlinked SNPs acquire a phenome-wide
genetic-effect correlation through a shared mediator. Writes fig1a_mediator_schematic.png next to this script."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

INK, GREY = "#222222", "#6b6b6b"
SNP, GENE, MED, PHE = "#e07b39", "#3b6ea8", "#4f9a55", "#e9b949"

fig, ax = plt.subplots(figsize=(7.4, 4.3), dpi=200)
ax.set_xlim(0, 10); ax.set_ylim(0, 6.2); ax.axis("off")

def box(x, y, w, h, text, fc, tc="white", fs=8.5, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=fc, ec="none"))
    ax.text(x, y, text, ha="center", va="center", color=tc, fontsize=fs, weight="bold" if bold else "normal")

def arrow(x0, y0, x1, y1, label=None, lx=0, ly=0, color=INK, ls="-", rad=0.0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11, lw=1.2,
                                 color=color, ls=ls, connectionstyle=f"arc3,rad={rad}"))
    if label:
        ax.text((x0 + x1) / 2 + lx, (y0 + y1) / 2 + ly, label, ha="center", va="center", fontsize=7.5,
                color=GREY, style="italic")

# chromosomes and SNPs
for x0, name in [(0.3, "Chr 1"), (4.5, "Chr 2")]:
    ax.add_patch(FancyBboxPatch((x0, 5.55), 3.0, 0.22, boxstyle="round,pad=0,rounding_size=0.11", fc="#d9dee6", ec="none"))
    ax.text(x0 + 3.05, 5.66, name, va="center", fontsize=7.5, color=GREY)
for x, lab in [(1.8, "SNP 1"), (6.0, "SNP 2")]:
    ax.plot([x, x], [5.5, 5.82], color=SNP, lw=3)
    box(x, 5.05, 1.05, 0.42, lab, SNP, bold=True)
ax.annotate("", xy=(5.4, 5.05), xytext=(2.4, 5.05), arrowprops=dict(arrowstyle="-", ls=(0, (3, 3)), color=GREY))
ax.text(3.9, 4.88, "unlinked (r = 0)", ha="center", va="top", fontsize=7.5, color=GREY)

# genes
box(1.8, 3.65, 2.3, 0.62, "Gene A\n(transcription factor)", GENE)
box(6.0, 3.65, 2.3, 0.62, "Gene B", GENE)
arrow(1.8, 4.82, 1.8, 3.98, "cis-eQTL", lx=0.48)
arrow(6.0, 4.82, 6.0, 3.98, "cis-eQTL", lx=0.48)
arrow(2.97, 3.65, 4.83, 3.65, "trans-regulation", ly=0.17)

# mediator and phenotypes
box(6.2, 2.35, 2.2, 0.62, "Mediator m:\nexpression of gene B", MED)
arrow(6.0, 3.32, 6.0, 2.68)
ys = [5.05, 4.25, 3.45, 2.55, 1.6]
labs = ["Phenotype 1", "Phenotype 2", "Phenotype 3", "\u22ee", "Phenotype K"]
for y, lab in zip(ys, labs):
    if lab == "\u22ee":
        ax.text(9.05, y, lab, ha="center", va="center", fontsize=12, color=INK)
        continue
    box(9.05, y, 1.6, 0.48, lab, PHE, tc=INK)
    arrow(7.32, 2.35, 8.23, y, rad=0.0)
ax.text(7.55, 1.55, "loadings Q", fontsize=7.5, color=GREY, style="italic", ha="center")

# model inset
ax.add_patch(FancyBboxPatch((0.25, 0.15), 4.6, 2.55, boxstyle="round,pad=0.02,rounding_size=0.12",
                            fc="#f4f5f7", ec="#c9ced6", lw=0.8))
lines = [
    (r"$\bf{Mediator\ model}$", 8.5),
    (r"$\mathbf{m} = \mathbf{P}\mathbf{x} + \mathbf{e},\quad \mathbf{y} = \mathbf{Q}\mathbf{m} + \mathbf{u}$", 9),
    (r"$\Rightarrow\ \mathbf{B} = \mathbf{P}^{\top}\mathbf{Q}^{\top},\quad \Sigma_\beta \propto \mathbf{P}^{\top} E[\mathbf{Q}^{\top}\mathbf{Q}]\,\mathbf{P}$", 9),
    ("Both SNPs act on the K traits only through m,", 7.5),
    ("so their effect rows are proportional to Q:", 7.5),
    (r"one mediator gives $|\rho_{12}| = 1$, several $|\rho_{12}| \leq 1$.", 7.5),
    (r"Hence $\rho_{12} \neq 0$ without any LD between the SNPs.", 7.5),
]
y = 2.42
for txt, fs in lines:
    ax.text(0.45, y, txt, fontsize=fs, color=INK, va="center")
    y -= 0.34

out = Path(__file__).with_suffix(".png")
fig.savefig(out, bbox_inches="tight", facecolor="white")
print(out)
