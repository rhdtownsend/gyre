#!/usr/bin/env python

import sys
import csv
import ads
import pickle

if len(sys.argv) > 1:
    infile = sys.argv[1]
    outfile = sys.argv[2]
else:
    infile = 'ads_refs.tsv'
    outfile = 'ads_refs.dat'

# Read the infile

ads_data = {}

with open(infile, 'r') as f:

    reader = csv.reader(f, delimiter='\t')
    next(reader, None)

    for line in reader:
        ref = line[0]
        bibcode = line[1]
        ads_data[ref] = bibcode

# Replace the bibcodes with the ADS articles

for key, value in ads_data.items():
    print('Processing key:', key)
    ads_data[key] = list(ads.SearchQuery(bibcode=value))[0]

# Write the data

with open(outfile, 'wb') as f:

    pickle.dump(ads_data, f)
