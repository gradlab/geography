# ECDC data

ESAC-Net antibiotic use data were scraped from the [CDC
website](https://ecdc.europa.eu/en/antimicrobial-consumption/database/distribution-by-antimicrobial-group)
into `esac.tsv`.

EARS-Net antibiotic resistance were downloaded from the [CDC
website](https://atlas.ecdc.europa.eu/public/index.aspx) as a single csv file.

The R script `raw/clean_data.R` cleans the downloaded data to create
`ecdc_data.tsv`.

## Attribution and disclaimer

Data from The European Surveillance System -- TESSy released by ECDC.

The views and opinions of the authors expressed herein do not necessarily state
or reflect those of ECDC. The accuracy of the authors’ statistical analysis and
the findings they report are not the responsibility of ECDC. ECDC is not
responsible for conclusions or opinions drawn from the data provided. ECDC is
not responsible for the correctness of the data and for data management, data
merging and data collation after provision of the data. ECDC shall not be held
liable for improper or incorrect use of the data.
