# Pipeline Nextflow pour Deep Mutational Scanning

## 1. Prérequis
- Nextflow installé
- Python 3 avec le script `dms_analysis.py`
- Deux FASTQ ou FASTQ.gz par échantillon : input et selected

## 2. Exécution locale
```bash
nextflow run main.nf \
  --samplesheet example_samplesheet.csv \
  -profile local
```

## 3. Exécution sur un cluster HPC (Slurm)
```bash
nextflow run main.nf \
  --samplesheet example_samplesheet.csv \
  -profile slurm
```

## 4. Exécution sur AWS Batch
```bash
nextflow run main.nf \
  --samplesheet example_samplesheet.csv \
  -profile aws \
  -work-dir s3://my-bucket/fastq-dms/work
```

## 5. Exécution sur GCP Batch
```bash
nextflow run main.nf \
  --samplesheet example_samplesheet.csv \
  -profile gcp \
  -work-dir gs://my-bucket/fastq-dms/work
```

## 6. Exécution simple sans samplesheet
```bash
nextflow run main.nf \
  --input input.fastq.gz \
  --selected selected.fastq.gz \
  --output results \
  -profile local
```

## 7. Paramètres utiles
- `--variant-start` : indice du début de la région mutante
- `--variant-end` : indice de fin de la région mutante
- `--min-length` : longueur minimale acceptée
- `--output` : dossier de sortie

## 8. Résultats
Le pipeline produit un fichier `.dms.tsv` par échantillon dans le dossier `results/`.

## 9. Extraire et traduire les 5 meilleurs variants
Après l'analyse DMS, exécuter :

```bash
python3 extract_top_variants.py results/sample_1.dms.tsv \
  --output-tsv results/sample_1.top5.translated.tsv \
  --output-fasta results/sample_1.top5.fasta
```

Le classement utilise `log2fc` décroissant. Le fichier TSV contient la séquence
ADN originale, son rang et la séquence protéique traduite. Le fichier FASTA est
directement utilisable comme entrée pour une étape de modélisation structurale.

Pour arrêter la traduction au premier codon stop :

```bash
python3 extract_top_variants.py results/sample_1.dms.tsv --stop-at-stop
```

Le cadre de lecture peut être changé avec `--frame 1` ou `--frame 2`.
