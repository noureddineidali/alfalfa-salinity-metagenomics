#!/usr/bin/env python3
"""Spot-check the first N synchronized paired FASTQ records; do not edit reads."""

import argparse
import gzip
from contextlib import ExitStack


def open_fastq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def record(handle):
    return [handle.readline().rstrip("\n") for _ in range(4)]


def identifier(header):
    return header.split()[0].removesuffix("/1").removesuffix("/2")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("r1")
    parser.add_argument("r2")
    parser.add_argument("--pairs", type=int, default=10000)
    args = parser.parse_args()
    if args.pairs < 1:
        parser.error("--pairs must be positive")
    with ExitStack() as stack:
        r1 = stack.enter_context(open_fastq(args.r1))
        r2 = stack.enter_context(open_fastq(args.r2))
        for i in range(1, args.pairs + 1):
            a, b = record(r1), record(r2)
            if not all(a) or not all(b):
                raise SystemExit(f"Incomplete record at pair {i}")
            if not a[0].startswith("@") or not b[0].startswith("@"):
                raise SystemExit(f"Invalid read header at pair {i}")
            if identifier(a[0]) != identifier(b[0]):
                raise SystemExit(f"Mate identifier mismatch at pair {i}: {a[0]} != {b[0]}")
            if a[2] != "+" or b[2] != "+" or len(a[1]) != len(a[3]) or len(b[1]) != len(b[3]):
                raise SystemExit(f"Invalid FASTQ sequence/quality record at pair {i}")
    print(f"PASS: first {args.pairs:,} paired FASTQ records have matching IDs and valid record structure")


if __name__ == "__main__":
    main()
