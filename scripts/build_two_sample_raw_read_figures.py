#!/usr/bin/env python3
"""Make the three descriptive figures for the two-sample raw-read report.

Inputs are validated stage counts, Bracken genus tables, and the published
supplementary genus matrix. This script does not perform inferential analysis.
Run from the repository root with the alfalfa-snakemake Python interpreter.
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import braycurtis


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/figures/raw-read-report"
OUT.mkdir(parents=True, exist_ok=True)
with (ROOT / "data/metadata/sample_manifest.tsv").open(newline="", encoding="utf-8") as handle:
    MANIFEST = {r["sample_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
assert (MANIFEST["CK-1"]["run_accession"], MANIFEST["CK-1"]["treatment"]) == ("SRR36739085", "Control")
assert (MANIFEST["MS-1"]["run_accession"], MANIFEST["MS-1"]["treatment"]) == ("SRR36739079", "Moderate salinity")
CK = "#2673A8"
MS = "#C96A48"
INK = "#263A47"
MUTED = "#62727D"
LIGHT = "#E8EEF2"
PALE = "#B8C9D3"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": INK,
    "svg.fonttype": "none", "figure.facecolor": "white", "axes.facecolor": "white",
})


def metric_table(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return {r["metric"]: r["value"] for r in csv.DictReader(handle, delimiter="\t")}


def bracken(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_taxid = {r["taxonomy_id"]: (r["name"], int(r["new_est_reads"])) for r in rows}
    assert len(by_taxid) == len(rows)
    total = sum(value for _, value in by_taxid.values())
    by_name = {name: count / total for name, count in by_taxid.values()}
    assert len(by_name) == len(rows)
    return {taxid: count / total for taxid, (_, count) in by_taxid.items()}, by_name


def published(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 15805
    result = {}
    for sample in ("CK-1", "MS-1"):
        values = {r["genus"]: float(r[sample]) for r in rows}
        assert len(values) == len(rows)
        total = sum(values.values())
        result[sample] = {name: value / total for name, value in values.items()}
    return result


ck_tax, ck_name = bracken(ROOT / "analysis/raw_read_comparison/CK1_pluspf8_c0_genus.bracken")
ms_tax, ms_name = bracken(ROOT / "results/taxonomy/raw-read-workflow/MS-1/MS-1_pluspf8_c0_genus.bracken")
pub = published(ROOT / "data/processed/genus_abundance.tsv")

ck_raw = 34283431
ck_clean = 34271117
ck_nonhost = 33808912
ck_classified = 1833747
ms_raw = int(metric_table(ROOT / "results/validation/raw-read-workflow/MS-1/raw_read_validation.tsv")["r1_reads"])
ms_clean = int(metric_table(ROOT / "results/validation/raw-read-workflow/MS-1/clean_read_validation.tsv")["r1_reads"])
ms_nonhost = int(metric_table(ROOT / "results/validation/raw-read-workflow/MS-1/nonhost_validation.tsv")["r1_reads"])
ms_classified = int(metric_table(ROOT / "results/validation/raw-read-workflow/MS-1/kraken2_validation.tsv")["classified_pairs"])
assert (ms_raw, ms_clean, ms_nonhost, ms_classified) == (34492089, 34479851, 33983612, 1756206)


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", pad_inches=0.16)
    fig.savefig(OUT / f"{stem}.png", dpi=210, bbox_inches="tight", pad_inches=0.16)
    plt.close(fig)


# Figure 1: account for each original paired fragment; enlarge classified fraction.
fig, (ax, zoom) = plt.subplots(1, 2, figsize=(12.4, 4.25), gridspec_kw={"width_ratios": [1.25, 0.8]})
samples = ["CK-1", "MS-1"]
values = [(ck_raw, ck_clean, ck_nonhost, ck_classified), (ms_raw, ms_clean, ms_nonhost, ms_classified)]
for i, (raw, clean, nonhost, classified) in enumerate(values):
    parts = [(raw-clean, PALE), (clean-nonhost, "#D5A044"),
             (nonhost-classified, LIGHT), (classified, "#337E89")]
    left = 0
    for count, color in parts:
        ax.barh(i, count / raw * 100, left=left, color=color, height=.46, edgecolor="white", linewidth=.4)
        left += count / raw * 100
    zoom.barh(i, classified / nonhost * 100, color="#337E89", height=.46)
    zoom.text(classified / nonhost * 100 + .11, i, f"{classified/nonhost*100:.2f}%", va="center", fontsize=12, weight="bold")
    ax.text(0, i-.37, f"{raw/1e6:.2f} M raw pairs  →  {nonhost/1e6:.2f} M non-host pairs", color=MUTED, fontsize=10)
for axis in (ax, zoom):
    axis.set_yticks([0, 1], samples)
    axis.set_ylim(1.6, -.85)
    axis.spines["left"].set_visible(False)
    axis.spines["bottom"].set_color(LIGHT)
    axis.tick_params(axis="y", length=0)
    axis.grid(axis="x", color=LIGHT, linewidth=.8)
    axis.set_axisbelow(True)
ax.set_xlim(0, 100)
ax.set_xlabel("Share of original paired fragments (%)")
ax.set_title("Read fate from archive to taxonomy", loc="left", weight="bold", pad=14)
zoom.set_xlim(0, 7)
zoom.set_xlabel("Classified / non-host fragments (%)")
zoom.set_title("Classification, enlarged scale", loc="left", weight="bold", pad=14)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(facecolor=PALE, label="fastp removed"), Patch(facecolor="#D5A044", label="host-associated"),
                    Patch(facecolor=LIGHT, label="non-host, unclassified"), Patch(facecolor="#337E89", label="classified")],
           loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(.5, -.035), fontsize=10)
fig.subplots_adjust(bottom=.23, left=.10, right=.97, top=.83, wspace=.36)
save(fig, "read_fate_ck1_ms1")


# Figure 2: identical genus set and within-Bracken denominator in both samples.
names = sorted(set(ck_name) | set(ms_name), key=lambda n: -(ck_name.get(n, 0) + ms_name.get(n, 0)))[:13]
names = sorted(names, key=lambda n: max(ck_name.get(n, 0), ms_name.get(n, 0)))
fig, ax = plt.subplots(figsize=(10.8, 6.8))
y = np.arange(len(names))
x1 = np.array([100 * ck_name.get(n, 0) for n in names])
x2 = np.array([100 * ms_name.get(n, 0) for n in names])
for yi, a, b in zip(y, x1, x2):
    ax.plot([a, b], [yi, yi], color="#C5D1D7", linewidth=2, zorder=1)
ax.scatter(x1, y, s=75, color=CK, label="CK-1", zorder=3)
ax.scatter(x2, y, s=75, color=MS, label="MS-1", zorder=3)
ax.set_yticks(y, names)
ax.tick_params(axis="y", length=0, pad=10)
ax.set_xlim(0, 10.5)
ax.set_xlabel("Abundance within each exported Bracken genus profile (%)")
ax.grid(axis="x", color=LIGHT)
ax.set_axisbelow(True)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_color(LIGHT)
ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.10), frameon=False, ncol=2)
ax.set_title("A shared view of the reconstructed genera", loc="left", weight="bold", pad=18)
taxids = sorted(set(ck_tax) | set(ms_tax))
distance = braycurtis([ck_tax.get(t, 0) for t in taxids], [ms_tax.get(t, 0) for t in taxids])
assert abs(distance - 0.17388273087616) < 1e-9
fig.text(.97, .035, f"Bray–Curtis = {distance:.3f}   ·   1,839 shared TaxIDs",
         ha="right", va="bottom", color=MUTED, fontsize=10)
fig.subplots_adjust(left=.24, right=.97, top=.88, bottom=.14)
save(fig, "shared_genus_composition_ck1_ms1")


# Figure 3: exact-label overlap, without implying unmatched taxa are absent.
fig, (ax, note) = plt.subplots(1, 2, figsize=(12.2, 4.6), gridspec_kw={"width_ratios": [1.3, .9]})
shares = []
distances = []
for sample, reconstructed in [("CK-1", ck_name), ("MS-1", ms_name)]:
    original = pub[sample]
    shared = set(original) & set(reconstructed)
    shares.append((100 * sum(original[n] for n in shared), 100 * sum(reconstructed[n] for n in shared)))
    labels = sorted(set(original) | set(reconstructed))
    distances.append(braycurtis([original.get(n, 0) for n in labels], [reconstructed.get(n, 0) for n in labels]))
assert abs(distances[0] - 0.8489889179026905) < 1e-9
assert abs(distances[1] - 0.8556155574) < 1e-8
y = np.arange(2)
ax.barh(y-.16, [x[0] for x in shares], height=.28, color="#8997A1", label="Published profile")
ax.barh(y+.16, [x[1] for x in shares], height=.28, color="#337E89", label="Reconstructed profile")
for i, (p, b) in enumerate(shares):
    ax.text(p+1, i-.16, f"{p:.1f}%", va="center", fontsize=10)
    ax.text(b+1, i+.16, f"{b:.1f}%", va="center", fontsize=10)
ax.set_yticks(y, samples)
ax.set_ylim(1.55, -.65)
ax.set_xlim(0, 102)
ax.set_xlabel("Profile abundance carried by exact shared genus labels (%)")
ax.set_title("The same sample, different taxonomic vocabulary", loc="left", weight="bold", pad=15)
ax.legend(loc="upper right", frameon=False, fontsize=10, ncol=2)
ax.grid(axis="x", color=LIGHT)
ax.set_axisbelow(True)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
note.axis("off")
note.text(.04, .87, "Exact-label profile distance", weight="bold", fontsize=13, transform=note.transAxes)
for j, sample in enumerate(samples):
    ypos = .62 - j*.29
    note.text(.04, ypos, sample, weight="bold", fontsize=12, transform=note.transAxes)
    note.text(.35, ypos, f"Bray–Curtis  {distances[j]:.3f}", fontsize=12, transform=note.transAxes)
    note.text(.35, ypos-.10, "Top-10 overlap  2 / 10", color=MUTED, fontsize=10, transform=note.transAxes)
note.text(.04, .03, "Unmatched labels reflect reference and nomenclature differences;\nexact-name matching is not a synonym or TaxID crosswalk.",
          color=MUTED, fontsize=9.5, transform=note.transAxes, va="bottom")
fig.subplots_adjust(left=.10, right=.98, top=.83, bottom=.16, wspace=.22)
save(fig, "published_vs_reconstructed_ck1_ms1")

print("created three figure pairs (SVG and PNG) in", OUT)
print("Bray CK-1/MS-1", distance)
print("published/reconstructed shares", shares)
print("exact-label Bray", distances)
