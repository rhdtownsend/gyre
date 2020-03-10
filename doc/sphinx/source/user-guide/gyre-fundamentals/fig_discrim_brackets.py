#!/usr/bin/env python
#
# Discriminant brackets plot

# Imports

import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import numpy as np
import numpy.linalg as la
import scipy.optimize as op

from colors import *

# Plot settings

plt.style.use('web.mplstyle')

# Calculation & plot parameters

N = 50

sigma_min = 0
sigma_max = 5.5
delta_sigma = 0.17

D_min = -0.3
D_max = 1.1

# System matrix

def M (sigma):

    tau = np.pi/(N-1)

    M = np.zeros([N, N])

    M[0,0] = 1.

    for i in range(1, N-1):
        M[i,i] = tau**2*sigma**2 - 2.
        M[i,i-1] = 1.
        M[i,i+1] = 1.

    M[N-1,N-1] = 1.

    return M

# Discriminant function

def discrim (sigma):

    return la.det(M(sigma))/((-1)**(N-2)*(N-1))

# Evaluate data

sigma = np.concatenate([np.arange(0, sigma_min, -delta_sigma)[1:], np.arange(0, sigma_max, delta_sigma)])
n_sigma = len(sigma)
D = np.empty(n_sigma)

for i in range(n_sigma):
    D[i] = discrim(sigma[i])

# Bracket roots

i_bracket = np.where(D[1:]*D[:-1] <= 0.)[0]

# Do the plot

fig, ax = plt.subplots()

ax.plot(sigma, D, 'o', color=SKY_BLUE, zorder=0)

ax.plot(sigma[i_bracket], D[i_bracket], 'o', mfc='None', mec=ORANGE, mew=2, ms=7, zorder=1)
ax.plot(sigma[i_bracket+1], D[i_bracket+1], 'o', mfc='None', mec=ORANGE, mew=2, ms=7, zorder=1)

ax.set_xlabel(r'$\sigma\ [\pi c/L]$')
ax.set_ylabel(r'$\mathcal{D}(\sigma)$')

ax.set_xlim(sigma_min, sigma_max)
ax.set_ylim(D_min, D_max)

ax.grid(True, which='both')

ax.xaxis.set_major_locator(tkr.MultipleLocator(1))
ax.xaxis.set_minor_locator(tkr.MultipleLocator(0.25))

ax.yaxis.set_major_locator(tkr.MultipleLocator(0.2))
ax.yaxis.set_minor_locator(tkr.MultipleLocator(0.1))

# Write out the figure

fig.tight_layout()
fig.savefig('fig_discrim_brackets.svg')
