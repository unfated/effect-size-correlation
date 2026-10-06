#!/usr/bin/env python3
"""Fig 6: fold enrichment of phenome-wide heritability in baseline-LF annotations
(joint annotation regression and model-free clump-level estimate).
Usage: 46_annot_figure.py <results_real_dir> <fig_dir> [label]"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

rd, fd = sys.argv[1:3]
lab = sys.argv[3] if len(sys.argv) > 3 else "all"
r = pd.read_csv(os.path.join(rd, f"annot_regression_{lab}.tsv"), sep="\t")
c = pd.read_csv(os.path.join(rd, f"lrcq_{lab}_ols.annot_clump.tsv"), sep="\t").set_index("annotation")
main = ["Coding_UCSC", "UTR_5_UCSC", "UTR_3_UCSC", "TSS_Hoffman", "Promoter_UCSC", "Conserved_LindbladToh",
        "Conserved_Mammal_phastCons46way", "Enhancer_Hoffman", "H3K9ac_Trynka", "H3K4me3_Trynka", "DGF_ENCODE",
        "TFBS_ENCODE", "FetalDHS_Trynka", "SuperEnhancer_Hnisz", "DHS_Trynka", "H3K4me1_Trynka", "H3K27ac_Hnisz",
        "Transcr_Hoffman", "Intron_UCSC", "Repressed_Hoffman"]
j = r[r["model"] == "joint"].set_index("annotation").reindex(main)
cc = c.reindex(main)
tab = pd.DataFrame(dict(prop=j["prop"], E_joint=j["E"], E_joint_se=j["E_se"], p_joint=j["p_E"],
                        E_clump=cc["E_clump"], E_clump_se=cc["E_clump_se"]))
tab.round(4).to_csv(os.path.join(rd, f"table_annotation_{lab}.tsv"), sep="\t")
fig, ax = plt.subplots(figsize=(7, 6))
y = range(len(main))
ax.errorbar(tab["E_joint"], [v + 0.15 for v in y], xerr=1.96 * tab["E_joint_se"], fmt="o", label="annotation regression (joint)")
ax.errorbar(tab["E_clump"], [v - 0.15 for v in y], xerr=1.96 * tab["E_clump_se"], fmt="s", ms=4, label="clump-level (model-free)")
ax.axvline(1, c="k", lw=0.7)
ax.set_yticks(list(y)); ax.set_yticklabels([m.replace("_", " ") for m in main], fontsize=8); ax.invert_yaxis()
ax.set_xlabel("fold enrichment of phenome-wide heritability"); ax.legend(fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(fd, "fig6_annotations.png"), dpi=150)
print(tab.round(2).to_string())
