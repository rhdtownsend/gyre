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

omega_min = 0
omega_max = 5.5
delta_omega = 0.17

n_min = 1
n_max = 3

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

# Evaluate discriminant data

omega = np.concatenate([np.arange(0, omega_min, -delta_omega)[1:], np.arange(0, omega_max, delta_omega)])
n_omega = len(omega)
D = np.empty(n_omega)

for i in range(n_omega):
    D[i] = discrim(omega[i])

# Bracket and find roots

i_bracket = np.where(D[1:]*D[:-1] <= 0.)[0]
n_bracket = len(i_bracket)

omega_root = np.empty(n_bracket)

for j in range(n_bracket):
    omega_root[j] = op.brentq(discrim, omega[i_bracket[j]], omega[i_bracket[j]+1])

# Do the plot

fig, ax = plt.subplots()

x = np.linspace(0., 1., N)

# Select the n'th root, and plot the eigenvector

colors = [SKY_BLUE, ORANGE, BLUE_GREEN, VERMILLION]

for n in range(n_min, n_max+1):

    u = sla.null_space(M(omega_root[n-1]))[:,0]
    u /= np.max(np.abs(u))

    ax.plot(x, u, 'o', color=colors[n-n_min], zorder=n, label=r'$n={:d}$'.format(n))

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
