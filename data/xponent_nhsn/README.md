# CDC (Xponent/NHSN) data

QuintilesIMS Xponent antibiotic use data were downloaded from the [CDC
website](https://gis.cdc.gov/grasp/PSA/indexAU.html).  NHSN antibiotic
resistance data were downloaded from the [CDC
website](https://gis.cdc.gov/grasp/PSA/indexAU.html).

The R script `raw/clean_data.R` cleans the data from the CDC website and
combines it into `use_resistance_data.tsv`.

## Data dictionary

- `unit`: Name of state
- `rx_person_year`: Prescriptions (treatments) per person per year
- `n_isolates`: Number of total *E. coli* isolates
- `n_resistant`: Number of quinolone-resistant isolates
