# plate-truth

**Does having a cultured relative mean a microbe will actually grow?**
A cross-study test, against real plate data, of the assumption behind every estimate of the "uncultured majority."

## The question

Estimates of how much of the microbial world can be cultured often rest on a shortcut: if an environmental 16S sequence closely matches a strain that someone has grown in a lab, the organism is treated as culturable. That shortcut is disputed (Lloyd et al. 2018; Martiny 2019; Steen et al. 2019), but it had not been tested directly against what actually grows on plates, across many studies.

This project measures how well 16S identity to the nearest cultured relative predicts whether an environmental taxon grows on plates, and how that relationship changes with habitat, taxon, medium, temperature, oxygen and incubation time.

Throughout, **"predicted to grow"** (a model or database output) is kept strictly separate from **"detected" or "enriched on plates"** (an observation).

## What was found

I reprocessed 11 public datasets from 10 independent studies that sequenced both an environmental sample and the communities grown from it on plates (soil, activated sludge, a marine sponge, oak, broiler carcasses, rumen): 665,023 combinations of sequence variant and plating condition. The hypotheses, models and criteria were registered in [ANALYSIS_PLAN.md](ANALYSIS_PLAN.md) before any model was fitted, and every performance number comes from a study held out of the fit.

![Growth on plates against identity to the nearest cultured relative, held out one study at a time](figures/fig2_calibration.png)

- **Identity ranks taxa, but does not set the level.** In a study the model had not seen, identity alone ranked which taxa grew (mean held-out AUC 0.801, 95% interval 0.715-0.870). But the overall rate of growth differed between studies more than tenfold (observed / expected 0.20 to 2.81) and could not be predicted from the plating conditions.
- **Most taxa with a close cultured relative did not grow.** For a new study, the predicted share enriched was 5.3% at 99% identity (95% interval 1.2-18.8%). No sequence without a cultured relative at 75% identity or above was enriched (0 of 3,059).
- **The shortcut overstates recovery in every study.** The share of reads with a cultured relative at 94.5% or above exceeded the share actually enriched on plates in all 10 studies, by a median factor of 5.91 (3.88-33.00).
- **Most families grow less often than relatedness predicts.** 48 of 80 evaluated families were systematically under-recovered relative to their identity to cultured strains, as were 8 of 14 phyla, including Acidobacteriota, Planctomycetota and Verrucomicrobiota; recovery concentrated in a few Pseudomonadota families.
- **Culture-collection media records are not shown to carry medium-specific information.** The registered test was positive, but the effect disappeared once each study had its own baseline.
- **In public surveys the shortcut exceeds what one plating could recover.** Across 28 MGnify survey studies from seven habitats, the shortcut counted 44% to 98% of reads as culturable. That is above what one plating is predicted to recover even in a high-recovery study, in every habitat and in both amplicon and shotgun data.

A cultured relative is useful for ranking which taxa are likely to grow, not as a measure of how much of a community can be grown.

The full draft is [manuscript/preprint.md](manuscript/preprint.md). Every number in it is filled in from `results/` by `scripts/47_manuscript.py`, and every reference is built from its Crossref record. Figures and their legends are in `figures/` and in the manuscript.

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

## Running the pipeline

The numbered scripts in `scripts/` run in order; each one's docstring gives its inputs, outputs and options. Downloads and sequence searches are cached, so every step can be resumed.

| Scripts | Step |
|---|---|
| 00-10 | Pilot on one soil dataset: reads to ASVs, cultured reference, identity, calibration |
| 11-14 | Literature and archive census of paired studies; sample sheet |
| 15-25 | Region references, reads to ASVs, identity, taxonomy, database, quality checks |
| 26-31 | Cross-study calibration model, held-out evaluation, sensitivity analyses |
| 32-38 | Media classes, culture-collection records, medium and taxon analyses |
| 39-45 | Public surveys: selection, processing, application of the model |
| 46-49 | Figures, manuscript, data release |

`scripts/21_verify_references.py` checks every study's DOI, first author, year and accession against Crossref and the paper's full text; it must pass before the manuscript is built.

## Repository layout

```
plate-truth/
├── README.md
├── ANALYSIS_PLAN.md        # pre-registered plan, with dated deviations
├── environment.yml
├── pyproject.toml
├── config/                 # thresholds, study lists, model settings (YAML)
├── data/
│   ├── raw/                # downloads (not in git)
│   ├── interim/            # intermediate files (not in git)
│   └── processed/          # census tables and small harmonized tables (in git)
├── src/plate_truth/        # importable package
├── scripts/                # numbered pipeline entry points, in run order
├── tests/                  # pytest
├── results/                # model outputs and tables
├── figures/                # figures
└── manuscript/             # preprint draft and its template
```

`data/plate_truth.sqlite` (the harmonized database) is rebuilt by the pipeline and is not kept in git. `scripts/49_data_release.py` exports it, with the census and the analysis table, as a documented data package for Zenodo.

## Author

Owen Swanson ([@oswan-1225](https://github.com/oswan-1225)).

## Citing

A bioRxiv preprint and a Zenodo data release are in preparation. Until then, please cite this repository.
