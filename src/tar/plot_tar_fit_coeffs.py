#!/usr/bin/env python
#
# Plot tar_fit coefficient data

import sys
import numpy as np
import matplotlib.pyplot as plt

import gyre_tar_fit
import gyre_cheb_fit

# Arguments

n = len(sys.argv)

if n != 3:
    raise Exception("Syntax: plot_tar_fit_coeffs infile outfile")

infile = sys.argv[1]
outfile = sys.argv[2]

# Set up the plot

fig, ax = plt.subplots()

ax.set_xlabel('k')
ax.set_ylabel('|c|')
ax.set_title(infile)

ax.set_yscale('log')

# Load the file

tf = gyre_tar_fit.TarFit.load(infile)

# Plot coefficients

if hasattr(tf, 'cf_neg'):
    ax.plot(np.abs(tf.cf_neg.coeffs())/np.max(np.abs(tf.cf_neg.coeffs())), label='cf_neg')
if hasattr(tf, 'cf_ctr'):
    ax.plot(np.abs(tf.cf_ctr.coeffs())/np.max(np.abs(tf.cf_ctr.coeffs())), label='cf_ctr')
if hasattr(tf, 'cf_pos'):
    ax.plot(np.abs(tf.cf_pos.coeffs())/np.max(np.abs(tf.cf_pos.coeffs())), label='cf_pos')

# Write out the figure

ax.legend(prop={'size':6})

fig.savefig(outfile)


