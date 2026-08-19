#!/usr/bin/env python3
"""Predict a protein structure with ESMFold and compute SASA and contacts."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

def read_fasta(path: str) -> tuple[str, str]:
    """Read one protein sequence from a FASTA file."""
    name = Path(path).stem
    sequence_parts: list[str] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence_parts:
                    break
                name = line[1:].split()[0] or name
            else:
                sequence_parts.append(line)

    sequence = "".join(sequence_parts).upper().replace(" ", "")
    if not sequence:
        raise ValueError(f"No sequence found in FASTA file: {path}")
    invalid = sorted(set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"))
    if invalid:
        raise ValueError(f"Invalid amino-acid symbols: {', '.join(invalid)}")
    return name, sequence


def predict_pdb(sequence: str, output_pdb: str, device: str | None = None) -> None:
    """Run ESMFold and write the predicted PDB file."""
    try:
        import torch
    except ImportError as error:
        raise RuntimeError("ESMFold requires PyTorch. Install requirements-structure.txt first.") from error

    try:
        import esm
    except ImportError as error:
        raise RuntimeError(
            "ESMFold requires the 'fair-esm' package. Install requirements-structure.txt first."
        ) from error

    model = esm.pretrained.esmfold_v1()
    selected_device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    model = model.eval().to(selected_device)
    if selected_device == "cpu":
        model.set_chunk_size(64)

    with torch.inference_mode():
        pdb_text = model.infer_pdb(sequence)

    output = Path(output_pdb)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(pdb_text, encoding="utf-8")


def residue_label(residue) -> str:
    """Return a stable chain/residue identifier."""
    chain = residue.get_parent().id
    number = residue.id[1]
    insertion = residue.id[2].strip()
    return f"{chain}:{number}{insertion}"


def calculate_sasa(pdb_path: str, output_tsv: str, probe_radius: float) -> None:
    """Calculate total and per-residue SASA with Biopython."""
    try:
        from Bio.PDB import PDBParser, ShrakeRupley
        from Bio.PDB.Polypeptide import is_aa
    except ImportError as error:
        raise RuntimeError("SASA calculation requires Biopython. Install requirements-structure.txt first.") from error

    structure = PDBParser(QUIET=True).get_structure("prediction", pdb_path)
    ShrakeRupley(probe_radius=probe_radius, n_points=960).compute(structure, level="R")

    output = Path(output_tsv)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["residue", "chain", "residue_number", "resname", "sasa_angstrom2"])
        for model in structure:
            for chain in model:
                for residue in chain:
                    if not is_aa(residue, standard=False):
                        continue
                    writer.writerow(
                        [
                            residue_label(residue),
                            chain.id,
                            residue.id[1],
                            residue.resname,
                            f"{residue.sasa:.3f}",
                        ]
                    )


def calculate_contacts(pdb_path: str, output_tsv: str, cutoff: float, min_sequence_separation: int) -> None:
    """Find residue contacts using the minimum heavy-atom distance."""
    try:
        from Bio.PDB import NeighborSearch, PDBParser
        from Bio.PDB.Polypeptide import is_aa
    except ImportError as error:
        raise RuntimeError("Contact calculation requires Biopython. Install requirements-structure.txt first.") from error

    structure = PDBParser(QUIET=True).get_structure("prediction", pdb_path)
    model = next(structure.get_models())
    residues = [
        residue
        for chain in model
        for residue in chain
        if is_aa(residue, standard=False)
    ]
    atoms = [
        atom
        for residue in residues
        for atom in residue
        if atom.element != "H"
    ]
    neighbors = NeighborSearch(atoms)
    contacts: dict[tuple[str, str], float] = {}

    for atom in atoms:
        for partner in neighbors.search(atom.coord, cutoff, level="A"):
            residue_a = atom.get_parent()
            residue_b = partner.get_parent()
            if residue_a is residue_b:
                continue
            chain_a = residue_a.get_parent().id
            chain_b = residue_b.get_parent().id
            number_a = residue_a.id[1]
            number_b = residue_b.id[1]
            if chain_a == chain_b and abs(number_a - number_b) < min_sequence_separation:
                continue

            label_a = residue_label(residue_a)
            label_b = residue_label(residue_b)
            key = tuple(sorted((label_a, label_b)))
            distance = float(atom - partner)
            contacts[key] = min(distance, contacts.get(key, float("inf")))

    output = Path(output_tsv)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["residue_a", "residue_b", "minimum_distance_angstrom"])
        for (residue_a, residue_b), distance in sorted(contacts.items()):
            writer.writerow([residue_a, residue_b, f"{distance:.3f}"])


def main() -> None:
    parser = argparse.ArgumentParser(description="ESMFold structure prediction, SASA and contact analysis")
    parser.add_argument("fasta", help="FASTA file containing one protein sequence")
    parser.add_argument("--output-dir", default="structure_results", help="Output directory")
    parser.add_argument("--device", choices=("cpu", "cuda"), help="Inference device; default: CUDA when available")
    parser.add_argument("--sasa-probe", type=float, default=1.4, help="SASA probe radius in Angstrom (default: 1.4)")
    parser.add_argument("--contact-cutoff", type=float, default=4.5, help="Heavy-atom contact cutoff in Angstrom (default: 4.5)")
    parser.add_argument("--min-sequence-separation", type=int, default=3, help="Exclude same-chain contacts closer than this separation")
    args = parser.parse_args()

    if args.device == "cuda":
        try:
            import torch
        except ImportError:
            parser.error("--device cuda requires PyTorch. Install requirements-structure.txt first.")
        if not torch.cuda.is_available():
            parser.error("--device cuda was requested but CUDA is unavailable")
    if args.sasa_probe <= 0 or args.contact_cutoff <= 0:
        parser.error("--sasa-probe and --contact-cutoff must be positive")

    name, sequence = read_fasta(args.fasta)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = output_dir / f"{name}.pdb"
    predict_pdb(sequence, str(pdb_path), args.device)
    calculate_sasa(str(pdb_path), str(output_dir / f"{name}.sasa.tsv"), args.sasa_probe)
    calculate_contacts(
        str(pdb_path),
        str(output_dir / f"{name}.contacts.tsv"),
        args.contact_cutoff,
        args.min_sequence_separation,
    )
    print(f"PDB: {pdb_path}")
    print(f"SASA: {output_dir / f'{name}.sasa.tsv'}")
    print(f"Contacts: {output_dir / f'{name}.contacts.tsv'}")


if __name__ == "__main__":
    main()
