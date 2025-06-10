#!/usr/bin/env python3
#
# Plot tar_fit data

import sys
import numpy as np
import matplotlib.pyplot as plt

import tar_fit
import cheb_fit

# Arguments

n = len(sys.argv)

if n < 6:
    raise Exception("Syntax: plot_tar_fit q_min q_max n_q infile1 [infile2...] outfile")

q_min = float(sys.argv[1])
q_max = float(sys.argv[2])
n_q = int(sys.argv[3])
infiles = sys.argv[4:n-1]
outfile = sys.argv[n-1]

# Set up the plot

fig, ax = plt.subplots()

ax.set_xlabel('q')
ax.set_ylabel('lambda')

ax.set_yscale('log')

# Loop over files

for infile in infiles:

    # Load the file

    tf = tar_fit.TarFit.load(infile)

    # Calculate data

    q = np.linspace(q_min, q_max, n_q)
    lam = np.vectorize(tf.lam)

    # Plot the data

    ax.plot(q, lam(q), label='m={:d}, k={:d}'.format(tf.m, tf.k))

# Write out the figure

ax.legend(prop={'size':6})

fig.savefig(outfile)


