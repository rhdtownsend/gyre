#!/usr/bin/env python
#
# Discriminant brackets plot (low-res grid)

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
M = 5
M_LINE = 1000

sigma_min = 0
sigma_max = 5.5
n_sigma = M
n_sigma_line = M_LINE

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

sigma_line = np.linspace(sigma_min, sigma_max, n_sigma_line)
D_line = np.empty(n_sigma_line)

for i in range(n_sigma_line):
    D_line[i] = discrim(sigma_line[i])

# Bracket roots

i_bracket = np.where(D[1:]*D[:-1] <= 0.)[0]

# Construct line segments for each interval

n_seg = n_sigma-1

sigma_seg = np.ndarray([n_seg], dtype=np.ndarray)
D_seg = np.ndarray([n_seg], dtype=np.ndarray)
color_seg = np.ndarray([n_seg], dtype=list)

for i in range(n_seg):
    j = np.where(np.logical_and(sigma_line >= sigma[i],sigma_line <= sigma[i+1]))[0]
    sigma_seg[i] = sigma_line[j]
    D_seg[i] = D_line[j]
    if D[i]*D[i+1] < 0.:
        color_seg[i] = ORANGE
    else:
        color_seg[i] = SKY_BLUE

# Do the plot

fig, ax = plt.subplots()

for i in range(n_seg):
    ax.plot(sigma_seg[i], D_seg[i], color=color_seg[i], alpha=0.4, zorder=0)

ax.plot(sigma, D, 'o', color=SKY_BLUE, zorder=0)

ax.plot(sigma[i_bracket], D[i_bracket], 'o', mfc='None', mec=ORANGE, mew=2, ms=7, zorder=1)
ax.plot(sigma[i_bracket+1], D[i_bracket+1], 'o', mfc='None', mec=ORANGE, mew=2, ms=7, zorder=1)

ax.set_xlabel(r'$\sigma\ [\pi c/L]$')
ax.set_ylabel(r'$\mathcal{D}(\sigma)$')

ax.set_xlim(0., sigma_max)
ax.set_ylim(D_min, D_max)

ax.grid(True, which='both', zorder=-1)

ax.xaxis.set_major_locator(tkr.MultipleLocator(1))
ax.xaxis.set_minor_locator(tkr.MultipleLocator(0.25))

ax.yaxis.set_major_locator(tkr.MultipleLocator(0.2))
ax.yaxis.set_minor_locator(tkr.MultipleLocator(0.1))

# Write out the figure

fig.tight_layout()
fig.savefig('fig_discrim_brackets_low.svg')
