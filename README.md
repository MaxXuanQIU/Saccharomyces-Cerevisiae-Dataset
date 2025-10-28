# Saccharomyces-Cerevisiae-Dataset

A ready-to-run multimodal dataset and processing code for protein function prediction in budding yeast (Saccharomyces cerevisiae). It ties together:

- Protein–protein interaction (PPI) network (STRING)
- Gene Ontology (GO) annotations (SGD GAF)
- Protein sequences (UniProt FASTA)
- Optional learned sequence embeddings (ESM-2)

The processing scripts build a consistent subset of proteins with both network context and GO labels, then optionally extract ESM-2 embeddings and align all artifacts into a single pickle for downstream ML.


## Repository layout

```
Saccharomyces-Cerevisiae-Dataset/
├── data/
│   ├── original_data/
│   │   ├── idmapping_2025_10_27.fasta         # UniProt FASTA for S. cerevisiae
│   │   ├── sgd.gaf                            # GO annotations from SGD (GAF format)
│   │   ├── protein_network_data.txt           # STRING network edges
│   │   ├── unique_STRING_proteins.txt         # Unique STRING protein IDs
│   │   └── list_of_STRING_proteins.txt        # Supplemental list (optional)
│   └── processed_data/                        # Outputs are written here
├── src/
│   ├── common.py                              # Paths, helpers, ESM loader
│   ├── data_processor.py                      # Cleans/merges network + GO; builds labels
│   ├── esm_feature_extractor.py               # Extracts ESM-2 embeddings and re-aligns dataset
│   └── data_inspection.py                     # Prints shapes/stats of the processed pickle
├── scripts/
│   ├── run_pipeline.ps1                       # Run preprocess → ESM extract → inspection
│   ├── preprocess.ps1                         # Only preprocessing
│   ├── extract_esm.ps1                        # Only ESM feature extraction
│   └── inspect.ps1                            # Only dataset inspection
├── requirements.txt                           # Core deps (install PyTorch separately)
├── LICENSE                                    # Apache-2.0
└── README.md
```


## What the pipeline does

1) Parse STRING PPI network and unique protein list.

2) Parse SGD GAF and keep high-confidence experimental GO evidence codes:
	 {EXP, IDA, IPI, IMP, IGI, IEP}.

3) Map GO annotations to STRING protein IDs using yeast systematic names
	 (e.g., YAL001C) found in GAF synonyms and symbols; STRING IDs use the
	 4932.<SYSTEMATIC_NAME> convention.

4) Select the top 50 most frequent GO terms among the high-confidence set.

5) Build a dataset of proteins that both appear in the network and have at
	 least one of the selected GO terms. Produce:
	 - labels: a binary matrix of shape (N_proteins, N_go_terms)
	 - network_edges: network edges filtered to the selected proteins
	 - protein_features: placeholder vectors (replaced by ESM features later)

6) Optionally extract ESM-2 embeddings (layer 33 mean-pooled, 1280-d) from
	 sequences and replace placeholders. The extractor aligns proteins, labels,
	 and edges to only those proteins with available sequences/features.


## Setup

Requirements:

- Python 3.9–3.11
- pip
- NumPy, pandas
- PyTorch (CPU or CUDA, per your hardware)
- fair-esm (ESM-2 models)

Install basics (Windows PowerShell):

```powershell
python -m venv .venv
. .\.venv\Scripts\Activate
pip install -U pip
pip install numpy pandas fair-esm
# Install PyTorch following the selector at https://pytorch.org/get-started/locally/
# Example (CPU-only):
# pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
```

Or install from the provided requirements file (PyTorch still installed separately):

```powershell
python -m venv .venv
. .\.venv\Scripts\Activate
pip install -U pip
pip install -r requirements.txt
# Then install PyTorch per your platform:
# pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
```

Notes on ESM:

- This repo loads ESM-2 650M (esm2_t33_650M_UR50D) via fair-esm.
- If you previously installed a different `esm` package and import fails, try:

```powershell
pip uninstall -y esm
pip install -U fair-esm
```


## Running the pipeline

The repo already includes the required inputs under `data/original_data`. Outputs are written to `data/processed_data`.

### Quick start (PowerShell helper scripts)

From the repo root:

```powershell
# Full pipeline
./scripts/run_pipeline.ps1

# Or run steps individually
./scripts/preprocess.ps1
./scripts/extract_esm.ps1
./scripts/inspect.ps1
```

If you encounter an execution policy warning, you can temporarily allow local scripts for the current user:

```powershell
Set-ExecutionPolicy RemoteSigned -Scope CurrentUser
```

### Alternative: run modules directly

1) Preprocess network + GO to build labels and split subset:

```powershell
python -m src.data_processor
```

This writes:

- `data/processed_data/processed_protein_data.pkl`
- `data/processed_data/data_statistics.txt`

2) Extract ESM-2 features and re-align the dataset to proteins with sequences/features:

```powershell
python -m src.esm_feature_extractor
```

This writes:

- `data/processed_data/processed_protein_data_with_esm.pkl`

3) Inspect the processed pickle (prints keys, shapes, and a small preview):

```powershell
python -m src.data_inspection
```


## Output contents

Both pickles use simple Python/Numpy containers:

- `proteins`: list[str] — STRING IDs (e.g., `4932.YAL001C`).
- `go_terms`: list[str] — selected GO term IDs (top 50 by frequency in the filtered set).
- `labels`: np.ndarray[float32] with shape (N, K) — multi-hot GO labels in the same order as `proteins` and `go_terms`.
- `network_edges`: list[tuple[str, str]] — filtered PPI edges between STRING IDs; only edges where both endpoints are in `proteins` are kept.
- `protein_features`:
	- In the first pickle: placeholder vectors (float32, 1280-d) for development.
	- In the ESM pickle: mean-pooled ESM-2 layer-33 embeddings (float32, 1280-d) for proteins with available sequences.


## Data sources and attribution

- STRING: protein–protein association data and protein ID namespace.
	- Site: https://string-db.org/
	- In this repo: `data/original_data/protein_network_data.txt`, `unique_STRING_proteins.txt`

- SGD (Saccharomyces Genome Database): GO annotations for S. cerevisiae (GAF format) and gene symbols/synonyms.
	- Site: https://www.yeastgenome.org/
	- In this repo: `data/original_data/sgd.gaf`

- UniProt: FASTA sequences with GN fields used to recover gene symbols.
	- Site: https://www.uniprot.org/
	- In this repo: `data/original_data/idmapping_2025_10_27.fasta`

Please respect the licenses/terms of the original providers when reusing these data.


## How GO mapping works (brief)

- Yeast systematic names (regex: `^Y[A-P][LR][0-9]{3}[CW](?:-[A-Z])?$`) are extracted from GAF synonyms/symbols.
- Each systematic name is prefixed with taxon code `4932.` to match STRING IDs.
- GO evidence is filtered to experimental codes {EXP, IDA, IPI, IMP, IGI, IEP}.
- The 50 most frequent GO terms in the filtered set are kept as targets.


## Troubleshooting

- ESM import error: uninstall any conflicting `esm` and install `fair-esm` (see above).
- CUDA OOM: reduce ESM batch size; the extractor uses 8 (GPU) or 2 (CPU) by default.
- File not found: ensure files exist under `data/original_data` with the exact filenames expected in `src/common.py`.


## Citation

If you use this repo or derived artifacts, please consider citing the original resources:

- STRING: Szklarczyk et al. (Latest release); https://string-db.org/
- SGD: The Saccharomyces Genome Database; https://www.yeastgenome.org/
- ESM-2: Rives et al., Lin et al.; https://github.com/facebookresearch/esm

You may also cite this repository (Zenodo DOI if available).


## License

This repository is released under the Apache License 2.0 (see `LICENSE`).


## Contributing

Issues and pull requests are welcome. For significant changes, please open an issue first to discuss what you would like to change.

