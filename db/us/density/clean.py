#!/usr/bin/env python3

from bs4 import BeautifulSoup
import csv

with open('2010 Census_ Population Density Data (Text Version).htm') as f:
    contents = f.read()

soup = BeautifulSoup(contents, 'html.parser')

table = soup.find(class_="datatable")

# the columns in the table are years
years = [x.text for x in table.contents[1].find_all('th')[1:]]

# extract rows from the table
body = table.contents[2]
rows = body.find_all("tr")

# parse the rows
data = []

for row in rows:
    # rows come in two flavors, depending on whether a th is present
    th = row.find("th")

    if th is not None:
        # this is a title row
        assert "class" in row.attrs
        assert row["class"][0] == "grouphead"
        state = row.text.strip()

        assert "id" in th.attrs
        group_id = th["id"]
    else:
        # this is a data row
        tds = row.find_all("td")
        assert "class" in tds[0].attrs
        assert tds[0]["class"][0] == "rowhead"

        row_title = tds[0].text
        row_id = tds[0]["id"]

        values = [x.text for x in tds[1:]]

        for year, value in zip(years, values):
            data.append({"state": state, "group_id": group_id, "row_title": row_title, "row_id": row_id, "year": year, "value": value})

# output only the 2010 densities
non_states = ["United States", "Puerto Rico", "District of Columbia"]
density2010 = [x for x in data if x["row_title"] == "People per sq. mile" and x["year"] == "2010" and not (x["state"] in non_states)]
assert len(density2010) == 50

density2010 = [{"state": x["state"], "density": x["value"]} for x in density2010]

with open("../density.tsv", "w") as f:
    fieldnames = ["state", "density"]
    writer = csv.DictWriter(f, fieldnames=fieldnames, dialect="excel-tab")
    writer.writeheader()

    for row in density2010:
        writer.writerow(row)
