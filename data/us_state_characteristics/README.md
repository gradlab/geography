# US State characteristics

The file `state_adjacency.tsv` encodes which US states are adjacent. Two states
are considered adjacent if they have any adjacent counties. Adjacent counties
were determined using the US Census Bureau's County Adjacency File.

The file `state_characteristics.tsv` includes metadata about each state used in
the adjacency analysis:

- `unit`: US state
- `region`: US Census region
- `division`: US Census division
- `population`: 2010 Census population estimates
- `area`: total state area (square kilometers) from US Census Bureau's State Area Measurements
- `temperature`: average temperature (degrees Fahrenheit, averaged over all NOAA stations in the state over 1981-2010)
- `income`: median income (dollars) from 2011-2015 American Community Survey
