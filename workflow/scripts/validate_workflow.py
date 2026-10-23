#!/usr/bin/env python3
"""Lightweight validations for the raw-read workflow."""

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


def write_metrics(path, metrics):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(metrics.items())


def read_metrics(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def read_seqkit(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 2 or "num_seqs" not in rows[0]:
        raise ValueError(f"Expected two SeqKit rows with num_seqs in {path}")
    counts = [int(float(row["num_seqs"])) for row in rows]
    if counts[0] != counts[1]:
        raise ValueError(f"Paired read counts differ in {path}: {counts}")
    return counts[0]


def command_version(command):
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    text = "\n".join(part for part in (completed.stdout, completed.stderr) if part).strip()
    return text.splitlines()[0] if text else f"exit status {completed.returncode}"


def versions(args):
    commands = {
        "prefetch": ["prefetch", "--version"],
        "fasterq-dump": ["fasterq-dump", "--version"],
        "FastQC": ["fastqc", "--version"],
        "fastp": ["fastp", "--version"],
        "MultiQC": ["multiqc", "--version"],
        "SeqKit": ["seqkit", "version"],
        "Bowtie2": ["bowtie2", "--version"],
        "Samtools": ["samtools", "--version"],
        "Kraken2": ["kraken2", "--version"],
        "Bracken": ["bracken", "-v"],
    }
    bio_prefix = Path(args.bio_prefix).resolve()
    rows = []
    for name, command in commands.items():
        executable = shutil.which(command[0])
        if executable is None:
            raise ValueError(f"Missing executable: {command[0]}")
        Path(executable).resolve().relative_to(bio_prefix)
        rows.append((name, executable, command_version(command)))
    rows.append(("workflow Python", sys.executable, sys.version.split()[0]))
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["tool", "executable", "version"])
        writer.writerows(rows)


def paired_stats(args):
    count = read_seqkit(args.stats)
    write_metrics(args.output, {
        "stage": args.stage,
        "r1_reads": count,
        "r2_reads": count,
        "paired_counts_equal": "TRUE",
    })


def host_counts(args):
    with open(args.counts, newline="", encoding="utf-8") as handle:
        counts = {row["category"]: int(row["pairs"]) for row in csv.DictReader(handle, delimiter="\t")}
    partition_names = ["both_unmapped", "r1_mapped_only", "r2_mapped_only", "both_mapped"]
    if "total" not in counts or any(name not in counts for name in partition_names):
        raise ValueError("Host pair-count table is incomplete")
    partition = sum(counts[name] for name in partition_names)
    if partition != counts["total"]:
        raise ValueError(f"Host categories sum to {partition}, expected {counts['total']}")
    clean = int(read_metrics(args.clean_validation)["r1_reads"])
    if counts["total"] != clean:
        raise ValueError(f"BAM pair total {counts['total']} differs from clean FASTQ count {clean}")
    one_mate = counts["r1_mapped_only"] + counts["r2_mapped_only"]
    strict_removed = counts["both_mapped"] + one_mate
    write_metrics(args.output, {
        "clean_pairs": clean,
        "categories_sum": partition,
        "both_unmapped_pairs": counts["both_unmapped"],
        "one_mate_only_pairs": one_mate,
        "both_mapped_pairs": counts["both_mapped"],
        "strict_removed_pairs": strict_removed,
        "strict_removed_percent": f"{100 * strict_removed / clean:.6f}",
        "status": "validated",
    })


def nonhost(args):
    count = read_seqkit(args.stats)
    expected = int(read_metrics(args.host_validation)["both_unmapped_pairs"])
    if count != expected:
        raise ValueError(f"Non-host FASTQ count {count} differs from expected {expected}")
    text = Path(args.extraction_log).read_text(encoding="utf-8", errors="replace")
    singleton = re.search(r"discarded\s+(\d+)\s+singletons", text)
    processed = re.search(r"processed\s+(\d+)\s+reads", text)
    if not singleton or int(singleton.group(1)) != 0:
        raise ValueError("Samtools extraction did not report zero discarded singletons")
    if not processed or int(processed.group(1)) != 2 * count:
        raise ValueError("Samtools processed-read total does not equal twice the paired count")
    write_metrics(args.output, {
        "r1_reads": count,
        "r2_reads": count,
        "paired_counts_equal": "TRUE",
        "discarded_singletons": 0,
        "status": "validated",
    })


def parse_kraken_log(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    matches = [
        re.search(r"(\d+) sequences .* processed", text),
        re.search(r"(\d+) sequences classified", text),
        re.search(r"(\d+) sequences unclassified", text),
    ]
    if not all(matches):
        raise ValueError(f"Could not parse Kraken2 totals from {path}")
    return tuple(int(match.group(1)) for match in matches)


def kraken(args):
    total, classified, unclassified = parse_kraken_log(args.log)
    expected = int(read_metrics(args.nonhost_validation)["r1_reads"])
    if total != expected or classified + unclassified != total:
        raise ValueError("Kraken2 totals are inconsistent with non-host input")
    with open(args.report, encoding="utf-8") as handle:
        first = handle.readline().rstrip("\n").split("\t")
    if len(first) < 3 or int(first[1]) != unclassified:
        raise ValueError("Kraken2 report unclassified count differs from execution log")
    write_metrics(args.output, {
        "input_pairs": total,
        "classified_pairs": classified,
        "unclassified_pairs": unclassified,
        "classified_percent": f"{100 * classified / total:.6f}",
        "status": "validated",
    })


def bracken(args):
    with open(args.abundance, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = {"name", "taxonomy_id", "taxonomy_lvl", "kraken_assigned_reads", "added_reads", "new_est_reads", "fraction_total_reads"}
    if not rows or set(rows[0]) != expected:
        raise ValueError("Unexpected or empty Bracken abundance table")
    taxids = [row["taxonomy_id"] for row in rows]
    values = [int(row["new_est_reads"]) for row in rows]
    if len(taxids) != len(set(taxids)):
        raise ValueError("Duplicate TaxIDs in Bracken abundance table")
    if any(value < 0 for value in values) or Path(args.report).stat().st_size == 0:
        raise ValueError("Invalid Bracken values or empty Bracken report")
    write_metrics(args.output, {
        "retained_genera": len(rows),
        "sum_new_est_reads": sum(values),
        "duplicate_taxids": 0,
        "negative_values": 0,
        "status": "validated",
    })


def ck1(args):
    paths = [args.provenance, args.r1, args.r2, args.kraken_report, args.bracken, args.metrics]
    for path in paths:
        if Path(path).stat().st_size == 0:
            raise ValueError(f"Missing or empty CK-1 checkpoint artifact: {path}")
    total, classified, unclassified = parse_kraken_log(args.kraken_log)
    if (total, classified, unclassified) != (33808912, 1833747, 31975165):
        raise ValueError("CK-1 Kraken2 checkpoint differs from the verified result")
    with open(args.bracken, newline="", encoding="utf-8") as handle:
        genera = list(csv.DictReader(handle, delimiter="\t"))
    if len(genera) != 1897:
        raise ValueError(f"Expected 1,897 CK-1 genera, observed {len(genera)}")
    metrics = json.loads(Path(args.metrics).read_text(encoding="utf-8"))
    if metrics["kraken2_total_pairs_from_session_log"] != total:
        raise ValueError("CK-1 comparison metrics disagree with Kraken2 log")
    provenance = Path(args.provenance).read_text(encoding="utf-8")
    phrases = [
        "extract strict paired non-host reads",
        "classify CK-1 non-host pairs with PlusPF-8",
        "estimate CK-1 genus abundances",
    ]
    if any(phrase not in provenance for phrase in phrases):
        raise ValueError("CK-1 executed-command provenance is incomplete")
    write_metrics(args.output, {
        "sample": "CK-1",
        "accession": "SRR36739085",
        "nonhost_pairs": total,
        "classified_pairs": classified,
        "classification_percent": f"{100 * classified / total:.6f}",
        "bracken_genera": len(genera),
        "comparison_metrics_consistent": "TRUE",
        "executed_command_provenance_present": "TRUE",
        "status": "validated",
    })


def build_parser():
    main = argparse.ArgumentParser()
    sub = main.add_subparsers(dest="command", required=True)
    definitions = {
        "versions": (["bio-prefix", "output"], versions),
        "paired-stats": (["stats", "stage", "output"], paired_stats),
        "host-counts": (["counts", "clean-validation", "output"], host_counts),
        "nonhost": (["stats", "host-validation", "extraction-log", "output"], nonhost),
        "kraken": (["report", "log", "nonhost-validation", "output"], kraken),
        "bracken": (["abundance", "report", "output"], bracken),
        "ck1": (["provenance", "r1", "r2", "kraken-log", "kraken-report", "bracken", "metrics", "output"], ck1),
    }
    for name, (options, function) in definitions.items():
        item = sub.add_parser(name)
        for option in options:
            item.add_argument(f"--{option}", required=True)
        item.set_defaults(function=function)
    return main


if __name__ == "__main__":
    arguments = build_parser().parse_args()
    arguments.function(arguments)
