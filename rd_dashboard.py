#!/usr/bin/env python3
"""Interactive R&D dashboard for ESMFold structures, SASA and contacts."""

from __future__ import annotations

import io
import tempfile
from pathlib import Path

import pandas as pd
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

from predict_structure_esmfold import (
    calculate_contacts,
    calculate_sasa,
    predict_pdb,
    read_fasta,
)


st.set_page_config(
    page_title="DMS Structure Lab",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


@st.cache_data(show_spinner=False)
def run_structure_metrics(pdb_text: str, probe_radius: float, contact_cutoff: float, min_separation: int):
    """Run Biopython metrics once per PDB/settings combination."""
    with tempfile.TemporaryDirectory() as directory:
        directory_path = Path(directory)
        pdb_path = directory_path / "structure.pdb"
        sasa_path = directory_path / "sasa.tsv"
        contacts_path = directory_path / "contacts.tsv"
        pdb_path.write_text(pdb_text, encoding="utf-8")
        calculate_sasa(str(pdb_path), str(sasa_path), probe_radius)
        calculate_contacts(str(pdb_path), str(contacts_path), contact_cutoff, min_separation)
        return sasa_path.read_text(encoding="utf-8"), contacts_path.read_text(encoding="utf-8")


@st.cache_data(show_spinner=False)
def predict_structure(fasta_bytes: bytes, device: str | None) -> tuple[str, str]:
    """Predict a PDB from FASTA and return its name and content."""
    with tempfile.TemporaryDirectory() as directory:
        fasta_path = Path(directory) / "input.fasta"
        pdb_path = Path(directory) / "prediction.pdb"
        fasta_path.write_bytes(fasta_bytes)
        _, sequence = read_fasta(str(fasta_path))
        predict_pdb(sequence, str(pdb_path), device)
        return pdb_path.name, pdb_path.read_text(encoding="utf-8")


def render_structure(pdb_text: str, highlighted_residues: list[str], representation: str) -> None:
    """Render the PDB with py3Dmol inside Streamlit."""
    viewer = py3Dmol.view(width=900, height=620)
    viewer.addModel(pdb_text, "pdb")
    viewer.setBackgroundColor("#101820")
    viewer.setStyle({"cartoon": {"color": "spectrum"}})

    if representation == "Surface":
        viewer.addSurface(
            py3Dmol.VDW,
            {"opacity": 0.72, "color": "white"},
            {"hetflag": False},
        )
    elif representation == "B-factor":
        viewer.setStyle({"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 0, "max": 100}}})

    for residue in highlighted_residues:
        chain, number = residue.split(":", 1)
        viewer.setStyle(
            {"chain": chain, "resi": int(number.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))},
            {"stick": {"colorscheme": "yellowCarbon", "radius": 0.22}, "cartoon": {"color": "yellow"}},
        )
    viewer.zoomTo()
    components.html(viewer._make_html(), height=650, scrolling=False)


def load_pdb(uploaded_file, fasta_file, device: str | None, predict_requested: bool) -> tuple[str, str] | None:
    """Return a structure name and PDB text from the selected input."""
    if uploaded_file is not None:
        return Path(uploaded_file.name).stem, uploaded_file.getvalue().decode("utf-8")
    if fasta_file is not None:
        fasta_bytes = fasta_file.getvalue()
        fasta_key = hash(fasta_bytes)
        if st.session_state.get("fasta_key") != fasta_key:
            st.session_state["fasta_key"] = fasta_key
            st.session_state.pop("predicted_structure", None)
        if not predict_requested:
            return st.session_state.get("predicted_structure")
        with st.spinner("ESMFold calcule la structure 3D..."):
            try:
                structure = predict_structure(fasta_bytes, device)
                st.session_state["predicted_structure"] = structure
                return structure
            except Exception as error:
                st.error(f"Échec de la prédiction ESMFold : {error}")
    return None


def main() -> None:
    st.title("DMS Structure Lab")
    st.caption("Exploration interactive des structures, surfaces accessibles et contacts moléculaires")

    with st.sidebar:
        st.header("Entrée")
        uploaded_pdb = st.file_uploader("Structure PDB", type=("pdb", "ent"))
        uploaded_fasta = st.file_uploader("Ou séquence FASTA pour ESMFold", type=("fasta", "fa", "faa"))
        device = st.selectbox("Accélérateur ESMFold", ("Auto", "CPU", "CUDA"))
        device_value = {"Auto": None, "CPU": "cpu", "CUDA": "cuda"}[device]
        predict_requested = st.button("Prédire avec ESMFold", type="primary", disabled=uploaded_fasta is None)
        st.divider()
        st.header("Paramètres d'analyse")
        probe_radius = st.number_input("Rayon de sonde SASA (Å)", min_value=0.1, value=1.4, step=0.1)
        contact_cutoff = st.number_input("Seuil de contact (Å)", min_value=0.1, value=4.5, step=0.1)
        min_separation = st.number_input("Séparation minimale de séquence", min_value=0, value=3, step=1)
        representation = st.radio("Vue 3D", ("Cartoon", "Surface", "B-factor"))

    if uploaded_pdb is not None and uploaded_fasta is not None:
        st.warning("Le PDB est prioritaire. Retire-le pour lancer une prédiction FASTA avec ESMFold.")

    structure = load_pdb(uploaded_pdb, uploaded_fasta, device_value, predict_requested)
    if structure is None:
        st.info("Charge un fichier PDB ou un FASTA pour commencer l'analyse.")
        st.stop()

    structure_name, pdb_text = structure
    try:
        sasa_text, contacts_text = run_structure_metrics(
            pdb_text, probe_radius, contact_cutoff, min_separation
        )
    except Exception as error:
        st.error(f"Échec de l'analyse Biopython : {error}")
        st.stop()

    sasa = pd.read_csv(io.StringIO(sasa_text), sep="\t")
    contacts = pd.read_csv(io.StringIO(contacts_text), sep="\t")
    sasa = sasa.sort_values("sasa_angstrom2", ascending=False).reset_index(drop=True)

    st.subheader(structure_name)
    summary_1, summary_2, summary_3 = st.columns(3)
    summary_1.metric("Résidus analysés", f"{len(sasa):,}")
    summary_2.metric("SASA totale (Å²)", f"{sasa['sasa_angstrom2'].sum():,.1f}")
    summary_3.metric("Contacts détectés", f"{len(contacts):,}")

    viewer_col, data_col = st.columns([1.35, 1])
    with viewer_col:
        st.subheader("Structure 3D")
        top_contact_residues: list[str] = []
        if not contacts.empty:
            top_contact_residues = sorted(
                set(contacts.head(20)["residue_a"]).union(contacts.head(20)["residue_b"])
            )
        render_structure(pdb_text, top_contact_residues, representation)

    with data_col:
        st.subheader("SASA par résidu")
        chart_data = sasa.set_index("residue")["sasa_angstrom2"].head(20)
        st.bar_chart(chart_data, horizontal=True, height=420)
        st.download_button(
            "Télécharger SASA TSV",
            sasa.to_csv(sep="\t", index=False).encode("utf-8"),
            file_name=f"{structure_name}.sasa.tsv",
            mime="text/tab-separated-values",
        )

    st.subheader("Contacts inter-résidus")
    st.dataframe(contacts, use_container_width=True, hide_index=True)
    st.download_button(
        "Télécharger contacts TSV",
        contacts.to_csv(sep="\t", index=False).encode("utf-8"),
        file_name=f"{structure_name}.contacts.tsv",
        mime="text/tab-separated-values",
    )
    st.download_button(
        "Télécharger structure PDB",
        pdb_text.encode("utf-8"),
        file_name=f"{structure_name}.pdb",
        mime="chemical/x-pdb",
    )


if __name__ == "__main__":
    main()
