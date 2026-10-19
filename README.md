# plate-truth

**Does having a cultured relative mean a microbe will actually grow?**
A cross-study test, against real plate data, of the assumption behind every estimate of the "uncultured majority."

## The question

Estimates of how much of the microbial world can be cultured often rest on a shortcut: if an environmental 16S sequence closely matches a strain that someone has grown in a lab, the organism is treated as culturable. That shortcut is disputed (Lloyd et al. 2018; Martiny 2019; Steen et al. 2019), but it has not been tested directly against what actually grows on plates, across many studies.

This project measures how well 16S identity to the nearest cultured relative predicts whether an environmental taxon grows on plates, and how that relationship changes with habitat, taxon, medium, temperature, oxygen and incubation time.

Throughout, **"predicted to grow"** (a model or database output) is kept strictly separate from **"detected" or "enriched on plates"** (an observation).

## Status

Early setup. The repository structure and environment are in place; no results yet. Progress will be reported here as phases complete.

## Setup

Requires [conda](https://docs.conda.io/) (Miniconda or Anaconda).

```bash
conda env create -f environment.yml
conda activate plate-truth
pip install -e . --no-deps
python scripts/00_get_vsearch.py     # Windows only; on Linux/macOS: conda install -c bioconda vsearch
pytest
```

Copy `.env.example` to `.env` and fill in only the keys you need. Never commit `.env`.

## Repository layout

```
plate-truth/
├── README.md
├── environment.yml
├── pyproject.toml
├── config/                 # thresholds, paths, study list (YAML)
├── data/
│   ├── raw/                # downloads (not in git)
│   ├── interim/            # intermediate files (not in git)
│   └── processed/          # small harmonized tables (in git)
├── src/plate_truth/        # importable package
├── scripts/                # numbered pipeline entry points, in run order
├── tests/                  # pytest
├── notebooks/              # exploration only, never imported
├── results/                # model outputs and tables
├── figures/                # figures and captions
├── app/                    # Streamlit tool
└── manuscript/             # preprint draft
```

`data/plate_truth.sqlite` (the harmonized database) is rebuilt by the pipeline and is not kept in git.

## Author

Owen Swanson ([@oswan-1225](https://github.com/oswan-1225)).

## Citing

Not yet available. A bioRxiv preprint and a Zenodo data release are planned.
