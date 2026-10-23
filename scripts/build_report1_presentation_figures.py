#!/usr/bin/env python3
"""Presentation-only figures from Report 1's validated genus matrix and PCoA scores.

No ordination, inferential test, filtering or scientific result is recomputed.
Run from any directory with the alfalfa-snakemake Python interpreter.
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/figures/report1-presentation"
OUT.mkdir(parents=True, exist_ok=True)
COLORS = {"CK": "#0072B2", "LS": "#E69F00", "MS": "#009E73"}
INK = "#263A47"
MUTED = "#62727D"
LIGHT = "#E8EEF2"
SAMPLES = [f"{t}-{i}" for t in ("CK", "LS", "MS") for i in (1, 2, 3)]
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": INK,
    "svg.fonttype": "none", "figure.facecolor": "white", "axes.facecolor": "white",
})


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", pad_inches=.18)
    fig.savefig(OUT / f"{stem}.png", dpi=210, bbox_inches="tight", pad_inches=.18)
    plt.close(fig)


with (ROOT / "data/processed/genus_abundance.tsv").open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
assert len(rows) == 15805 and list(rows[0]) == ["genus", *SAMPLES]
totals = {sample: sum(float(r[sample]) for r in rows) for sample in SAMPLES}
ranked = sorted(rows, key=lambda r: -sum(float(r[s]) / totals[s] for s in SAMPLES) / len(SAMPLES))
selected = ranked[:6]
assert [r["genus"] for r in selected] == [
    "SCGC-AG-212-J23", "Gp6-AA40", "Sediminibacterium",
    "PSRF01", "12-FULL-67-14b", "Bradyrhizobium",
]


# Compact sample-level composition view, replacing two crowded older figures.
fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.7), sharex=True, sharey=True)
offsets = [-.13, 0, .13]
for ax, row in zip(axes.flat, selected):
    for pos, treatment in enumerate(("CK", "LS", "MS")):
        values = [100 * float(row[f"{treatment}-{i}"]) / totals[f"{treatment}-{i}"] for i in (1, 2, 3)]
        ax.scatter([pos + d for d in offsets], values, s=55, color=COLORS[treatment], alpha=.88,
                   edgecolor="white", linewidth=.6, zorder=3)
        ax.plot([pos-.19, pos+.19], [np.mean(values)]*2, color=COLORS[treatment], lw=2.1, zorder=4)
    ax.set_title(row["genus"], fontsize=11, fontweight="bold", loc="left", pad=10)
    ax.set_xticks([0, 1, 2], ["CK", "LS", "MS"])
    ax.set_xlim(-.40, 2.40)
    ax.set_ylim(0, 4.25)
    ax.grid(axis="y", color=LIGHT)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", length=0)
    ax.spines["left"].set_color(LIGHT)
    ax.spines["bottom"].set_color(LIGHT)
for ax in axes[:, 0]:
    ax.set_ylabel("Within-sample abundance (%)")
fig.suptitle("Leading genus labels, sample by sample", x=.075, y=.98, ha="left", fontsize=16, weight="bold")
fig.text(.075, .93, "Six highest mean proportions across nine samples · points are individual samples · bars mark treatment means",
         color=MUTED, fontsize=10, ha="left")
fig.subplots_adjust(left=.075, right=.985, top=.85, bottom=.08, wspace=.23, hspace=.45)
save(fig, "leading_genera_sample_patterns")


# Re-present the saved ordination coordinates; axis direction and values unchanged.
with (ROOT / "results/tables/pcoa_scores.tsv").open(newline="", encoding="utf-8") as handle:
    scores = list(csv.DictReader(handle, delimiter="\t"))
assert [r["sample"] for r in scores] == SAMPLES
fig, ax = plt.subplots(figsize=(9.6, 6.4))
ax.axhline(0, color=LIGHT, lw=1, zorder=0)
ax.axvline(0, color=LIGHT, lw=1, zorder=0)
offsets = {
    "CK-1": (-.006, .005), "CK-2": (-.006, -.007), "CK-3": (.004, .006),
    "LS-1": (.004, .004), "LS-2": (.004, .001), "LS-3": (-.005, -.004),
    "MS-1": (.004, -.005), "MS-2": (.004, .005), "MS-3": (-.013, -.006),
}
for treatment in ("CK", "LS", "MS"):
    block = [r for r in scores if r["treatment"] == treatment]
    xs = [float(r["PCoA1"]) for r in block]
    ys = [float(r["PCoA2"]) for r in block]
    ax.scatter(xs, ys, s=155, color=COLORS[treatment], edgecolor="white", linewidth=1.4,
               zorder=3)
    ax.scatter(np.mean(xs), np.mean(ys), s=125, color=COLORS[treatment], marker="D",
               edgecolor="white", linewidth=1, zorder=4)
    for r, x, y in zip(block, xs, ys):
        dx, dy = offsets[r["sample"]]
        ax.text(x+dx, y+dy, r["sample"], fontsize=10.5, color=COLORS[treatment], weight="bold")
ax.set_xlim(-.129, .078)
ax.set_ylim(-.069, .069)
ax.set_xlabel("PCoA1 (64.8% of positive-eigenvalue variation)")
ax.set_ylabel("PCoA2 (16.1%)")
ax.set_title("Genus-level community structure", loc="left", fontsize=16, weight="bold", pad=14)
ax.spines["left"].set_color(LIGHT)
ax.spines["bottom"].set_color(LIGHT)
fig.subplots_adjust(left=.12, right=.97, top=.91, bottom=.13)
save(fig, "genus_bray_curtis_pcoa_refined")

print("Created two presentation figure pairs in", OUT)
