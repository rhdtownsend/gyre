#!/usr/bin/env python3
#
# Simple script for stress-testing GYRE

import sys
import os
import argparse
import subprocess
import tempfile
import astropy.table as at
import astropy.io as ai
import numpy as np
import time

# Parse arguments

parser = argparse.ArgumentParser(description='GYRE stress test')

parser.add_argument('template_file', help='template GYRE inlist file')
parser.add_argument('output_file', help='output table file')
parser.add_argument('-t', '--threads', help='number of threads (can be a comma-separated list')
parser.add_argument('-c', '--count', type=int, help='count of runs to average over')
parser.add_argument('-n', '--name', help='name of template variable')
parser.add_argument('-v', '--values', help='values of template variable (comma-separated list)')

args = parser.parse_args()

# Process arguments

if args.threads is not None:
    threads = args.threads.split(',')
else:
    threads = (1,)

if args.count is not None:
    count = args.count
else:
    runs = 1

if args.name is not None:
    name = args.name
else:
    raise Exception('--name must be specified')

if args.values is not None:
    values = args.values.split(',')
else:
    raise Exception('--values must be specified')

# Read the template file

with open(args.template_file, 'r') as f:
    template = f.read()

# Timing routine

def measure_time(f):

    times = []

    for i in range(count):

        # Run GYRE, capturing output

        start_time = time.perf_counter()

        try:
            result = subprocess.run(('./gyre', f'{f.name}'), capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            print('gyre failed with error', e.stderr)
            
        end_time = time.perf_counter()

        times += [end_time - start_time]

    return np.mean(times)

# Set up the table

if len(threads) > 1:
    tbl = at.Table(names=(name, 'threads', 'exec. time (s)'),
                   dtype=(str, int, float))
else:
    tbl = at.Table(names=(name, 'exec. time (s)'),
                   dtype=(str, float))

# Perform the run
            
for v in values:

    with tempfile.NamedTemporaryFile(mode='w+t') as f:

        # Create the input file

        f.write(template.format_map({name: v}))
        f.flush()

        # Do the timing

        if len(threads) > 1:

            first = True

            for t in threads:

                os.environ['OMP_NUM_THREADS'] = f'{t}'

                if first:
                    tbl.add_row([v, t, measure_time(f)])
                    first = False
                else:
                    tbl.add_row(['', t, measure_time(f)])

        else:

            os.environ['OMP_NUM_THREADS'] = f'{threads[0]}'

            tbl.add_row([v, measure_time(f)])

# Print out results
        
for col in tbl.itercols():
    if col.info.dtype.kind == 'f':
        col.info.format = '.2f'

tbl.write(args.output_file, overwrite=True)
