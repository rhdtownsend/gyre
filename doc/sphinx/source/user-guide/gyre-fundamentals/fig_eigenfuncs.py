#!/usr/bin/env python
#
# Eigenfunction plot

# Imports

import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import scipy.optimize as op

from colors import *

# Plot settings

plt.style.use('web.mplstyle')

# Calculation & plot parameters

N = 50
M = 32

sigma_min = 0
sigma_max = 5.5
n_sigma = M

n_min = 1
n_max = 3

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

# Evaluate discriminant data

sigma = np.linspace(sigma_min, sigma_max, n_sigma)
D = np.empty(n_sigma)

for i in range(n_sigma):
    D[i] = discrim(sigma[i])

# Bracket and find roots

i_bracket = np.where(D[1:]*D[:-1] <= 0.)[0]
n_bracket = len(i_bracket)

sigma_root = np.empty(n_bracket)

for j in range(n_bracket):
    sigma_root[j] = op.brentq(discrim, sigma[i_bracket[j]], sigma[i_bracket[j]+1])

# Do the plot

fig, ax = plt.subplots()

x = np.linspace(0., 1., N)
x_ana = np.linspace(0., 1., 1000)

# Select the n'th root, and plot the eigenvector

colors = [SKY_BLUE, ORANGE, BLUE_GREEN, VERMILLION]

for n in range(n_min, n_max+1):

    y = sla.null_space(S(sigma_root[n-1]), rcond=1.)[:,-1]
    y *= np.sin(n*np.pi*x[1])/y[1]

    y_ana = np.sin(n*np.pi*x_ana)

    ax.plot(x_ana, y_ana, '-', color=colors[n-n_min], zorder=n)
    ax.plot(x, y, 'o', color=colors[n-n_min], zorder=n, ms=5, label=r'$n={:d}$'.format(n))

ax.set_xlabel(r'$x/L$')
ax.set_ylabel(r'$\tilde{y}$')

ax.legend()

ax.set_xlim(0., 1.)
ax.set_ylim(-1.1, 1.1)

ax.grid(True, which='both')

ax.xaxis.set_major_locator(tkr.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(tkr.MultipleLocator(0.1))

ax.yaxis.set_major_locator(tkr.MultipleLocator(0.5))
ax.yaxis.set_minor_locator(tkr.MultipleLocator(0.25))

# Write out the figure

fig.tight_layout()
fig.savefig('fig_eigenfuncs.svg')
