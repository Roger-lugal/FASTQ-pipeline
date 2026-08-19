#!/usr/bin/env python3
"""Extract and translate the top variants from a DMS TSV result."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_dna(sequence: str, frame: int = 0, stop_at_stop: bool = False) -> str:
    """Translate a DNA sequence; ambiguous codons become X."""
    sequence = sequence.strip().upper().replace("U", "T")
    protein = []
    for index in range(frame, len(sequence) - 2, 3):
        amino_acid = CODON_TABLE.get(sequence[index:index + 3], "X")
        if stop_at_stop and amino_acid == "*":
            break
        protein.append(amino_acid)
    return "".join(protein)


def read_top_variants(input_path: str, top_n: int) -> list[dict[str, str]]:
    """Read and rank variants by descending log2 fold-change."""
    with open(input_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"variant", "log2fc"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required TSV columns: {', '.join(sorted(missing))}")
        rows = list(reader)

    rows.sort(key=lambda row: float(row["log2fc"]), reverse=True)
    return rows[:top_n]


def write_outputs(rows: list[dict[str, str]], output_tsv: str, output_fasta: str, frame: int, stop_at_stop: bool) -> None:
    """Write translated variants to TSV and protein FASTA files."""
    output_tsv_path = Path(output_tsv)
    output_fasta_path = Path(output_fasta)
    output_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    output_fasta_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(rows[0].keys()) + ["rank", "protein_sequence"] if rows else ["rank", "protein_sequence"]
    with output_tsv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            row = dict(row)
            row["rank"] = rank
            row["protein_sequence"] = translate_dna(row["variant"], frame, stop_at_stop)
            writer.writerow(row)

    with output_fasta_path.open("w", encoding="utf-8") as handle:
        for rank, row in enumerate(rows, start=1):
            protein = translate_dna(row["variant"], frame, stop_at_stop)
            score = row["log2fc"]
            handle.write(f">variant_{rank}|log2fc={score}\n{protein}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract and translate top DMS variants")
    parser.add_argument("input_tsv", help="DMS TSV produced by dms_analysis.py")
    parser.add_argument("--top", type=int, default=5, help="Number of top variants to extract (default: 5)")
    parser.add_argument("--frame", type=int, choices=(0, 1, 2), default=0, help="Reading frame (default: 0)")
    parser.add_argument("--stop-at-stop", action="store_true", help="Stop translation at the first stop codon")
    parser.add_argument("--output-tsv", default="top_variants_translated.tsv", help="Translated TSV output")
    parser.add_argument("--output-fasta", default="top_variants.fasta", help="Protein FASTA output")
    args = parser.parse_args()

    if args.top < 1:
        parser.error("--top must be at least 1")

    rows = read_top_variants(args.input_tsv, args.top)
    write_outputs(rows, args.output_tsv, args.output_fasta, args.frame, args.stop_at_stop)
    print(f"Extracted {len(rows)} variants")
    print(f"Translated TSV: {args.output_tsv}")
    print(f"Protein FASTA: {args.output_fasta}")


if __name__ == "__main__":
    main()
