# The potential for "spillover" in outpatient antibiotic stewardship interventions among US states

*author*: Scott Olesen <olesen@hsph.harvard.edu>

This repository contains code to reproduce the simulations and a subset of the
empirical analyses in the manuscript.

## Simulations

The simulations are self-contained in the R script `simulations.R`, which
produces figures in the folder `fig/`. The script reproduces:

- Figures 1b-e from this manuscript
- Supplemental Figure 2 from this manuscript
- Figure 2g with *k = 1* from [Davies *et al*](https://dx.doi.org/10.1038/s41559-018-0786-x).
- Figure 3 from [Lehtinen *et al*.](https://dx.doi.org/10.1073/pnas.1617849114)

## Empirical analyses

The `data/` folder contains readme's that explain the sources of two of the
datasets, Xponent/NHSN and ECDC, as well as the metadata for US states and
European countries.

The manuscript also includes an analysis of the MarketScan/ResistanceOpen
dataset, which is not publicly available. We do not have license to post that
data.

The R script `analyze_observational.R` reproduces (parts of) Figures 2-4 and
Supplemental Tables 2-5 in the folder `fig/`.
