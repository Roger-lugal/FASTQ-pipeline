# FASTQ-pipeline

Pipeline reproductible de Deep Mutational Scanning (DMS), d'analyse structurale
et d'exploration interactive pour la recherche en biologie moléculaire.

## Fonctionnalités

- comptage de variants à partir de FASTQ/FASTQ.gz et calcul du `log2fc` ;
- exécution Nextflow en local, sur Slurm, AWS Batch ou Google Batch ;
- extraction et traduction des meilleurs variants ;
- prédiction de structure 3D avec ESMFold ;
- calcul de la SASA et des contacts avec Biopython ;
- dashboard Streamlit et py3Dmol pour l'exploration R&D.

## Installation

Python 3.10+ est recommandé. Pour la partie structurelle et le dashboard :

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-structure.txt
```

Pour une exécution GPU, installer une version de PyTorch compatible avec le
driver CUDA de la machine avant `fair-esm`.

## Démarrage rapide DMS

```bash
python3 dms_analysis.py \
	--input input.fastq.gz \
	--selected selected.fastq.gz \
	--output results/dms_scores.tsv
```

Extraire et traduire les cinq meilleurs variants :

```bash
python3 extract_top_variants.py results/dms_scores.tsv \
	--output-tsv results/top5.translated.tsv \
	--output-fasta results/top5.fasta
```

## Pipeline Nextflow

```bash
nextflow run main.nf \
	--samplesheet example_samplesheet.csv \
	--output results \
	-profile local
```

Profils disponibles : `local`, `slurm`, `aws` et `gcp`. Les chemins de stockage,
la région cloud, le projet GCP et la queue Slurm doivent être adaptés dans
`nextflow.config` avant une exécution distante. Le détail des commandes se
trouve dans [NEXTFLOW.md](NEXTFLOW.md).

## Structure 3D et dashboard

Prédire une structure et calculer ses métriques :

```bash
python3 predict_structure_esmfold.py protein.fasta \
	--output-dir structure_results
```

Lancer le dashboard :

```bash
streamlit run rd_dashboard.py
```

Le dashboard accepte un FASTA ou un PDB, affiche la structure 3D interactive,
la SASA par résidu et les contacts détectés, et permet de télécharger les
résultats PDB/TSV.

## Validation

```bash
python3 -m py_compile \
	dms_analysis.py extract_top_variants.py \
	predict_structure_esmfold.py rd_dashboard.py
```

## Licence

Ce projet est distribué sous licence MIT. Voir [LICENSE](LICENSE).