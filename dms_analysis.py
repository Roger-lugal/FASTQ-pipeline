#!/usr/bin/env python3
"""
Deep Mutational Scanning (DMS) analysis from FASTQ libraries.

This script:
1. parses FASTQ files from an input library and a selected library,
2. extracts a mutant sequence per read,
3. counts variant abundance,
4. estimates enrichment as log2 fold-change,
5. exports a TSV summary for downstream analysis.

Usage:
    python dms_analysis.py --input input.fastq --selected selected.fastq --output dms_scores.tsv
    python dms_analysis.py --input input.fastq --selected selected.fastq --variant-start 10 --variant-end 200
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path


def open_fastq(path: str):
    """Open a FASTQ or FASTQ.gz file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_fastq(path: str):
    """Yield (header, sequence) pairs from a FASTQ file."""
    with open_fastq(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            handle.readline()  # +
            handle.readline()  # quality
            if sequence:
                yield header.strip(), sequence.upper()


def extract_variant(sequence: str, start: int = 0, end: int | None = None) -> str:
    """Extract the mutant region from a read sequence."""
    if end is None:
        end = len(sequence)
    return sequence[start:end]


def count_variants(fastq_path: str, variant_start: int = 0, variant_end: int | None = None, min_length: int = 1):
    """Count variants from a FASTQ file."""
    counts = defaultdict(int)

    for _, sequence in parse_fastq(fastq_path):
        variant = extract_variant(sequence, start=variant_start, end=variant_end)
        if len(variant) < min_length:
            continue
        if "N" in variant:
            continue
        counts[variant] += 1

    return counts


def compute_enrichment(input_counts, selected_counts, pseudocount: float = 0.5):
    """Compute log2 enrichment score for each variant."""
    all_variants = sorted(set(input_counts) | set(selected_counts))
    total_input = sum(input_counts.values())
    total_selected = sum(selected_counts.values())

    rows = []
    for variant in all_variants:
        input_count = input_counts.get(variant, 0)
        selected_count = selected_counts.get(variant, 0)

        input_fraction = (input_count + pseudocount) / (total_input + pseudocount * len(all_variants))
        selected_fraction = (selected_count + pseudocount) / (total_selected + pseudocount * len(all_variants))

        log2fc = math.log2(selected_fraction / input_fraction)
        activity = 1.0 / (1.0 + math.exp(-log2fc))

        rows.append(
            {
                "variant": variant,
                "input_count": input_count,
                "selected_count": selected_count,
                "log2fc": round(log2fc, 6),
                "activity": round(activity, 6),
            }
        )

    rows.sort(key=lambda x: x["log2fc"], reverse=True)
    return rows


def write_summary(rows, output_path: str):
    """Write DMS scores to a TSV file."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["variant", "input_count", "selected_count", "log2fc", "activity"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Summary saved to: {output}")


def main():
    parser = argparse.ArgumentParser(description="Deep Mutational Scanning enrichment analysis")
    parser.add_argument("--input", required=True, help="FASTQ file for the input library")
    parser.add_argument("--selected", required=True, help="FASTQ file for the selected library")
    parser.add_argument("--variant-start", type=int, default=0, help="Start index of variant region in each read")
    parser.add_argument("--variant-end", type=int, default=None, help="End index of variant region in each read")
    parser.add_argument("--min-length", type=int, default=1, help="Minimum sequence length accepted")
    parser.add_argument("--output", default="dms_scores.tsv", help="Output TSV file")
    args = parser.parse_args()

    input_counts = count_variants(
        args.input,
        variant_start=args.variant_start,
        variant_end=args.variant_end,
        min_length=args.min_length,
    )
    selected_counts = count_variants(
        args.selected,
        variant_start=args.variant_start,
        variant_end=args.variant_end,
        min_length=args.min_length,
    )

    rows = compute_enrichment(input_counts, selected_counts)
    write_summary(rows, args.output)

    print(f"Input variants: {len(input_counts)}")
    print(f"Selected variants: {len(selected_counts)}")
    print("Top variants by log2FC:")
    for row in rows[:10]:
        print(f"  {row['variant']}: log2FC={row['log2fc']}, activity={row['activity']}")


if __name__ == "__main__":
    main()
