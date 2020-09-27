#!/usr/bin/env python
#
# Generate tar_fit files for the data/tar directory

import os

# Parameters

K_MIN = -5
K_MAX = 5

M_MIN = -5
M_MAX = 5

TOLER = 1E-10

# Loop over parameter combinations

for k in range(K_MIN, K_MAX+1):
    for m in range(M_MIN, M_MAX+1):

        if k < 0 and m == 0: continue

        print("Processing m={:d}, k={:d}".format(m, k))

        infix = 'm{:+d}.k{:+d}'.format(m, k)

        # Run build_tar_fit

        os.system('./build_tar_fit {:d} {:d} {:e} ../../data/tar/tar_fit.{:s}.h5'.format(m, k, TOLER, infix))
