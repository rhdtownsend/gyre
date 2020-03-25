#!/usr/bin/env python
#
# Discriminant function plot

# Imports

import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import numpy as np
import numpy.linalg as la

from colors import *

# Plot settings

plt.style.use('web.mplstyle')

# Calculation & plot parameters

N = 50
M = 1000

sigma_min = 0
sigma_max = 5.5
n_sigma = M

D_min = -0.3
D_max = 1.1

# System matrix

def S (sigma):

    tau = np.pi/(N-1)

    S = np.zeros([N, N])

    S[0,0] = 1.

    for i in range(1, N-1):
        S[i,i] = tau**2*sigma**2 - 2.
        S[i,i-1] = 1.
        S[i,i+1] = 1.

    S[N-1,N-1] = 1.

    return S

# Discriminant function

def discrim (sigma):

    return la.det(S(sigma))/((-1)**(N-2)*(N-1))

# Evaluate data

sigma = np.linspace(sigma_min, sigma_max, n_sigma)
D = np.empty(n_sigma)

for i in range(n_sigma):
    D[i] = discrim(sigma[i])

# Do the plot

fig, ax = plt.subplots()

ax.plot(sigma, D, color=SKY_BLUE, zorder=0)

ax.set_xlabel(r'$\sigma\ [\pi c/L]$')
ax.set_ylabel(r'$\mathcal{D}(\sigma)$')

ax.set_xlim(sigma_min, sigma_max)
ax.set_ylim(D_min, D_max)

ax.grid(True, which='both')

ax.xaxis.set_major_locator(tkr.MultipleLocator(1))
ax.xaxis.set_minor_locator(tkr.MultipleLocator(0.25))

ax.yaxis.set_major_locator(tkr.MultipleLocator(0.2))
ax.yaxis.set_minor_locator(tkr.MultipleLocator(0.1))

# Add analytic roots

tau = np.pi/(N-1)

for n in range(1, N):
    sigma_ana = np.sqrt(2-2*np.cos(n*tau))/tau
    ax.scatter([sigma_ana,-sigma_ana], [0.,0.], 25, color=ORANGE, zorder=1)

# Write out the figure

fig.tight_layout()
fig.savefig('fig_discrim_func.svg')
