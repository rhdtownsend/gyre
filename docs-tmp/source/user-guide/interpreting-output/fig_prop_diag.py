#!/usr/bin/env python
#
# Displacement eigenfunctions

# Imports

import matplotlib.pyplot as plt
import pygyre as pg

from colors import *

# Plot settings

plt.style.use('web.mplstyle')

# Read data

d = pg.read_output('detail.l1.n-7.h5')

# Evaluate characteristic frequencies

l = d.meta['l']
omega = d.meta['omega'].real
n_pg = d.meta['n_pg']

x = d['x']
V = d['V_2']*x**2
As = d['As']
c_1 = d['c_1']
Gamma_1 = d['Gamma_1']

N2 = d['As']/d['c_1']
Sl2 = l*(l+1)*Gamma_1/(V*c_1)

# Do the plot

fig, ax = plt.subplots()

ax.plot(x, N2, color=SKY_BLUE, label=r'$N^{2}$', zorder=1)
ax.plot(x, Sl2, color=ORANGE, label=r'$S_{\ell}^{2}$', zorder=0)

ax.axhline(omega**2, color=BLUE, dashes=(4,2), label=r'$n_{{\rm pg}}={:d}$'.format(n_pg))

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$\omega^{2}$')

ax.set_xlim(0,1)

ax.set_ylim(5e-2, 5e2)
ax.set_yscale('log')

ax.legend(loc=1)

fig.tight_layout()
fig.savefig('fig_prop_diag.svg')
