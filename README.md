# Measuring Political Opinion Polarisation from Ordinal Survey Data

This repository contains the unified reproduction code for a Master's thesis on measuring political opinion polarisation from ordinal survey data in Spain, the United Kingdom, and Germany.

The empirical workflow estimates latent political-opinion scores from the POL-AXES survey and applies Duclos-Esteban-Ray-style polarisation summaries to the resulting score distributions. Kernel density estimation uses the alpha-specific Duclos bandwidth rule, `h_alpha = 4.7 * n^(-1/5) * sigma * alpha^(1/5)`. The main reproducibility script is:

```text
scripts/reproduce_current_results.R
```

## Repository contents

- `scripts/reproduce_current_results.R`: unified R script for the main empirical reproduction workflow.

## Data availability

The raw POL-AXES data file and source PDFs are not included in this repository. The script expects the data file at:

```text
TFM/Scripts/POL-AXES data.csv
```

To reproduce the empirical outputs, place the POL-AXES CSV in the same folder as the reproduction script and run the script from the TFM project folder.

## Software

The reproduction script is written in R and uses `lavaan` for ordinal latent-variable modelling.

## Citation

Goiria Kortajarena, J. (2026). *Measuring Political Opinion Polarisation from Ordinal Survey Data*.
