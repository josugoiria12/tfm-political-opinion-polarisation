# Measuring Political Opinion Polarisation from Ordinal Survey Data

This repository contains the unified reproduction code for a Master's thesis on measuring political opinion polarisation from ordinal survey data in Spain, the United Kingdom, and Germany.

The empirical workflow estimates latent political-opinion scores from the POL-AXES survey and applies Duclos-Esteban-Ray-style polarisation summaries to the resulting score distributions. The main reproducibility script is:

```text
scripts/reproduce_current_results.R
```

## Repository contents

- `scripts/reproduce_current_results.R`: unified R script for the main empirical reproduction workflow.

## Data availability

The raw POL-AXES data file and source PDFs are not included in this repository. The scripts expect the data file at:

```text
raw/32020284/POL-AXES data.csv
```

To reproduce the empirical outputs, place the POL-AXES CSV at that path and run the scripts from the project root.

## Software

The reproduction script is written in R and uses `lavaan` for ordinal latent-variable modelling.

## Citation

Goiria Kortajarena, J. (2026). *Measuring Political Opinion Polarisation from Ordinal Survey Data*.
