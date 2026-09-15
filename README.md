# Dietary readout from an MS/MS food reference library

Code and food-enriched MS/MS feature tables accompanying the manuscript describing **food composition scores**: a reference-based strategy that matches untargeted LC-MS/MS data from biological samples against MS/MS spectra that are enriched in specific foods, and summarises the matches as a per-food score.

- **Web app:** https://foodreadouts.gnps2.org/
- **App documentation / SOP:** https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/metaboapp_Food_Readout/
- **Supplementary data (Zenodo):** [10.5281/zenodo.18203282](https://doi.org/10.5281/zenodo.18203282)

## Repository contents

```
Code/
  01_food_enriched_features_500foods.Rmd      Food-enriched features from the 500-food reference set
  02_food_enriched_features_CDTREAT_foods.Rmd Food-enriched features from 92 locally sourced UK foods
  03_NIST_stool_food_readout.Rmd              Food readout of NIST vegetarian/omnivore stool samples
  04_ITIS_mouse_diet_study.R                  Food readout of ITIS vs Western diet-fed mice
Files/
  food_enriched_features_level3.csv           Food-enriched features, food ontology level 3
  food_enriched_features_level3_annotated.csv   ...subset with a GNPS library annotation
  food_enriched_features_level5.csv           Food-enriched features, food ontology level 5
  food_enriched_features_level5_annotated.csv   ...subset with a GNPS library annotation
```

See [Files/README.md](Files/README.md) for column descriptions.

### Scripts, figures and data sources

| Script | Manuscript | Inputs | Public data |
|---|---|---|---|
| `01_food_enriched_features_500foods.Rmd` | Figure 1 | MZmine feature table (`featuretable_reformated.csv`), GNPS2 annotations (`merged_results_with_gnps.tsv`), food ontology metadata (`Metadata_500food.csv`) | MassIVE MSV000096584 |
| `02_food_enriched_features_CDTREAT_foods.Rmd` | Figure 5 (reference library for the Crohn's urine cohort) | MZmine feature table (`UK_foods_iimn_gnps_quant.csv`), food classes (`Food_classes.csv`) | MassIVE MSV000097255 |
| `03_NIST_stool_food_readout.Rmd` | Figure 3 | Food-enriched features (`Biomarkers.csv`), GNPS2 library matches (`library_matches.tsv`), MZmine feature table (`NIST_Vegan_Ominvore_iimn_gnps_quant.csv`), sample metadata (`metadata.csv`) | MassIVE MSV000098937; [FBMN](https://gnps2.org/status?task=640fe3da659b44b6a7bac0bee1985cb6), [library matching](https://gnps2.org/status?task=2a20536c11874dcaa4e5753143e07140) |
| `04_ITIS_mouse_diet_study.R` | Figure 2b-c | Sample metadata (`metadata_RA_mice_2025.csv`), food readout output (`food_metadata.tsv`) from the web app | MassIVE MSV000097333; [FBMN](https://gnps2.org/status?task=057283987ac946b9b220e0094faef7fb), [library matching](https://gnps2.org/status?task=b00468d72ac642f8b1109aa627f1649c) |

Input files are not stored in this repository. Feature tables are produced by MZmine from the MassIVE datasets listed above (processing parameters are given in the manuscript Methods), and annotation/library-matching tables are downloadable from the linked GNPS2 tasks. Place the inputs next to the script before running it.

## Requirements

Tested on Windows 11 with R 4.4.1 (RStudio). Packages used:

| Package | Version tested |
|---|---|
| tidyverse | 2.0.0 |
| mixOmics (Bioconductor) | 6.28.0 |
| ggpubr | 0.6.0 |
| ggrepel | 0.9.6 |
| vegan | 2.7.2 |
| CoDaSeq | 0.99.7 |
| knitr | 1.50 |

```r
install.packages(c("tidyverse", "ggpubr", "ggrepel", "vegan", "knitr"))
BiocManager::install("mixOmics")
remotes::install_github("ggloor/CoDaSeq/CoDaSeq")
```

## Obtaining a food readout for your own data

1. Run the GNPS2 library search workflow on your `.mgf` file against the food MS/MS reference library.
2. Enter the library search task ID and your feature quantification table in the food readout app (https://foodreadouts.gnps2.org/) and download the food composition scores.
3. Analyse the scores in the app or in the statistical software of your choice (the scripts in `Code/` are worked examples).
