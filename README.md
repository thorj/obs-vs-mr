# Title of the Paper (Short)

This repository contains the R code used to generate the figures, tables, and Bayesian model results for the paper:

> **Systematic comparison of observational and Mendelian Randomization estimates for cardiometabolic proteomic signatures**  
> Thorarinn Jonmundsson, Valur Emilsson, Elisabet A Frick, Heida Bjarnadottir, Eva Jacobsen, Thor Aspelund, Lenore J Launer, Joseph J Loureiro, Anthony P Orth, Nancy Finkel, Vilmundur Gudnason, Valborg Gudmundsdottir

## Overview

The analyses in this repository reproduce all main and supplementary results reported in the manuscript, including:

- Observational proteomic analyses
- Mendelian randomization results
- Hierarchical Bayesian models implemented in Stan
- Figures and tables appearing in the main text and supplement

The code is organized to allow exact reproduction of reported results from precomputed model outputs.

---

## Repository Structure

...


---

## Data Availability

Individual-level genetic and proteomic data used in this study cannot be shared publicly due to ethical and legal restrictions associated with the AGES–Reykjavik Study.

Specifically:
- Genotype data
- Individual-level proteomic measurements
- Linked phenotype data

### Provided Alternatives

To support reproducibility and transparency, this repository includes:

- Full Stan model source code
- Fitted model objects
- Model summaries used in the manuscript
- R scripts that recreate all figures

Researchers with approved access to the underlying data may adapt the provided scripts to rerun the analyses.

See `data/README.md` for details.

---

## Bayesian Models

Bayesian models were implemented in **Stan** and fitted using:

- R version: `x`
- Stan version: `x.y.z`
- Interface: `cmdstanr` / `rstan`

All reported results are based on the final fitted models stored in `results/`.

Model files:
- `models/model_name.stan` — description

Posterior samples were generated with fixed random seeds.

---

## Reproducing Figures and Tables

To reproduce all figures and tables from the manuscript