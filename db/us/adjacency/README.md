# Adjacency

https://www.census.gov/geo/reference/county-adjacency.html

County Adjacency File
Overview:
The county adjacency file lists each county, or county equivalent, and which county, or counties, are neighboring.  The file includes all 50 states, the District of Columbia, Puerto Rico and the Island Areas (American Samoa, the Commonwealth of the Northern Mariana Islands, Guam, and the U.S. Virgin Islands).

Download:
Download the County Adjacency File [txt]

File Layout:
The county adjacency file is a national, tab delimited text file.  The list is sorted by County GEOID and then neighboring counties are sorted by Neighbor GEOID within each county.  The Island Areas appear at the end of the list.

The table below provides the layout of the county adjacency file:

Column	Column Name	Column Description
1
County Name	2010 State and County Name
2
County GEOID	2010 State and County FIPS Codes
3
Neighbor Name	2010 State and County name of each neighboring county or county equivalent
4
Neighbor GEOID	2010 State and County FIPS Codes of each neighboring county or county equivalent
Notes:
In some instances the boundary between two counties may fall within a body of water, so it seems as if the two counties do not physically touch.  These counties are included on the list as neighbors.
Every county has itself listed as a neighboring county.

# FIPS

https://www.census.gov/geo/reference/codes/cou.html

2010 FIPS Codes for Counties and County Equivalent Entities
National and state files containing FIPS codes for counties and county equivalent entities are available for download.

Download
Select a state, territory, or United States from this list to view the file:


File Format and Record Layout
These text files contain comma-delimited records for each county. The records are of the format:

Field Name	Field Description	Example
FIPS Class Codes
H1:  identifies an active county or statistically equivalent entity that does not qualify under subclass C7 or H6.
H4:  identifies a legally defined inactive or nonfunctioning county or statistically equivalent entity that does not qualify under subclass H6.
H5:  identifies census areas in Alaska, a statistical county equivalent entity.
H6:  identifies a county or statistically equivalent entity that is areally coextensive or governmentally consolidated with an incorporated place, part of an incorporated place, or a consolidated city.
C7:  identifies an incorporated place that is an independent city; that is, it also serves as a county equivalent because it is not part of any county, and a minor civil division (MCD) equivalent because it is not part of any MCD.
STATE	State Postal Code	FL
STATEFP	State FIPS Code	12
COUNTYFP	County FIPS Code	011
COUNTYNAME	County Name and Legal/Statistical Area Description	Broward County
CLASSFP	FIPS Class Code	H1


County Changes
The inventory, names, and codes for counties and equivalent areas change periodically. Please review Substantial Changes to Counties and County Equivalent Entities: 1970-Present.
