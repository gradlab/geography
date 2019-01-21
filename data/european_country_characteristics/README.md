# European country characteristics

The file `europe_adjacency.tsv` encodes which European countries are adjacent.
This data was drawn from the [Correlates of War
project](http://wwww.correlatesofwar.org/data-sets/direct-contiguity) using
their definition of adjacency.

The file `europe_characteristics.tsv` includes metadata about each country used
in the adjacency analysis:

- `unit`: Name of the country
- `eu_region`: UN geoscheme sub-region
- `area`: Land and inland water area, in square kilometers, from the [UN Statistics Division](http://unstats.un.org/unsd/environment/totalarea.htm)
- `population`: 2010 population estimates from [UN World Population Prospect](http://population.un.org/wpp)
- `temperature`: Mean annual temperature in degrees Celsius, averaged over 1961-1999, from [World Bank Climate Change Knowledge Portal](http://datacatalog.worldbank.org/dataset/climate-change-knowledge-portal-historical-data). Malta's value, which is not in the World Bank data, was set at 19.5 degrees Celsius.
- `income`: Per capita GDP at constant 2010 prices in US dollars, averaged over 2011-2015, from [UN Statistics Division](http://unstats.un.org/unsd/snaama/introduction.asp)
