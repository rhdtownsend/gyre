#!/usr/bin/env python
#
# Plot trad_fit data

import sys
import numpy as np
import matplotlib.pyplot as plt

import gyre_tar_fit
import gyre_cheb_fit

# Arguments

n = len(sys.argv)

if n < 6:
    raise Exception("Syntax: plot_tar_fit nu_min nu_max n_nu infile1 [infile2...] outfile")

nu_min = float(sys.argv[1])
nu_max = float(sys.argv[2])
n_nu = int(sys.argv[3])
infiles = sys.argv[4:n-1]
outfile = sys.argv[n-1]

# Set up the plot

fig, ax = plt.subplots()

ax.set_xlabel('nu')
ax.set_ylabel('lambda')

ax.set_yscale('log')

# Loop over files

for infile in infiles:

    # Load the file

    tf = gyre_tar_fit.TarFit.load(infile)

    # Calculate data

    nu = np.linspace(nu_min, nu_max, n_nu)
    lam = np.vectorize(tf.lam)

    # Plot the data

    ax.plot(nu, lam(nu), label='m={:d}, k={:d}'.format(tf.m, tf.k))

# Write out the figure

ax.legend(prop={'size':6})

fig.savefig(outfile)


