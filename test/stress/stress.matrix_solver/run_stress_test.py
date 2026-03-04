#!/usr/bin/env python

import sys
import re
import subprocess
import tempfile
import astropy.table as at
import astropy.io as ai
import numpy as np
import time

# Parameters

tmpl_file = sys.argv[1]
n_reps = int(sys.argv[2])

# Read the template file

with open(tmpl_file, 'r') as f:
    template = f.read()

# Define run variables

MATRIX_SOLVERS = ['BANDED', 'CYCLIC', 'ROWPP']

# Loop through run variables, storing timing data in table

tbl = at.Table(names=('matrix_solver', 'runtime'),
               dtype=(str, float))

for matrix_solver in MATRIX_SOLVERS:

    with tempfile.NamedTemporaryFile(mode='w+t') as f:

        # Create the input file

        f.write(template.format(matrix_solver=matrix_solver))
        f.flush()

        # Loop over repititions

        times = []

        for i in range(n_reps):

            # Run GYRE, capturing output

            start_time = time.perf_counter()

            try:
                result = subprocess.run(('./gyre', f'{f.name}'), capture_output=True, text=True, check=True)
            except subprocess.CalledProcessError as e:
                print('gyre failed with error', e.stderr)

            end_time = time.perf_counter()

            times += [end_time - start_time]

        tbl.add_row([matrix_solver, np.mean(times)])

# Print out results
        
for col in tbl.itercols():
    if col.info.dtype.kind == 'f':
        col.info.format = '.2f'
        
ai.ascii.write(tbl, sys.stdout, format='rst')
