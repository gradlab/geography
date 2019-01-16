# Outpatient antibiotic stewardship interventions: geographic scale and associations between use and resistance

*author*: Scott Olesen <olesen@hsph.harvard.edu>

This repository contains code to reproduce the simulations and a subset of the empirical analyses in the manuscript.

## Simulations

The simulations are self-contained in the R script `simulations.R`, which
produces figures in the folder `fig/`. The script reproduces:

- Figures 1bc and 2 from this manuscript
- Supplemental Figures 1 and 2 from this manuscript
- Figure 2g with $k=1$ from Davies *et al*.
- Figure 3 from Lehtinen *et al*.

## Empirical analyses

The `data/` folder contains a readme and an R script `clean_data.R` that will
download and clean the CDC (IMS Quintiles and NHSN) antibiotic use and
resistance data.

The manuscript also includes an analysis of two other datasets. The
MarketScan/ResistanceOpen dataset is not publicly available and we do not have
license to post that data. The ECDC dataset is publicly available but we do not
have license to post the cleaned version used in this analysis. The scraping
and cleaning process is slow and would not make for a good demonstration of the
methods used for analysis.

The analyses of the CDC data are in `analyze_cdc.R`, which reproduces one panel
each of Supplemental Figures 3, 4, and 5 in the folder `fig/`. (Figures 3 and
4, and the other panels of Supplemental Figures 3, 4, and 5 are the same
analyses applied to the MarketScan/ResistanceOpen and ECDC datasets.) The
script also reproduces one of the lines in Supplemental Table 2.

The confidence intervals reported in the empirical results may differ slightly
from those shown in the manuscript because the bootstrapping process is
inherently random.
