# Food-enriched MS/MS feature tables

Food-enriched MS/MS features identified from the 500-food reference set (MassIVE MSV000096584) with `Code/01_food_enriched_features_500foods.Rmd`. The same tables are deposited on Zenodo ([10.5281/zenodo.18203282](https://doi.org/10.5281/zenodo.18203282)).

| File | Features | Food categories |
|---|---|---|
| `food_enriched_features_level3.csv` | 7,282 | 11 (food ontology level 3, e.g. dairy, legume, fleshy fruit) |
| `food_enriched_features_level3_annotated.csv` | 503 | 11 |
| `food_enriched_features_level5.csv` | 6,127 | 120 (food ontology level 5, individual foods, e.g. cheddar cheese, salmon, walnut) |
| `food_enriched_features_level5_annotated.csv` | 235 | 115 |

## Columns

| Column | Description |
|---|---|
| `feature` | Feature ID (MZmine `row ID`, equal to the GNPS2 feature-based molecular networking scan number) |
| `category` | Comma-separated food categories in which the feature is enriched |
| `Compound_Name` | GNPS2 spectral library annotation (`*_annotated.csv` files only) |
