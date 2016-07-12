#!/usr/bin/env python

import os

# Parameters

K_MIN = -5
K_MAX = 5

M_MIN = -5
M_MAX = 5

TOLER = 0.

# Loop over parameter combinations

for k in range(K_MIN, K_MAX+1):
    for m in range(M_MIN, M_MAX+1):

        print("Processing m={:d}, k={:d}".format(m, k))

        infix = 'm{:+d}.k{:+d}'.format(m, k)

        # Run build_trad_fit

        os.system('./build_trad_fit {:d} {:d} {:e} ../../data/trad/trad_fit.{:s}.h5'.format(m, k, TOLER, infix))
