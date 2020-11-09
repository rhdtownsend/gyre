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

# Do the plot

fig, ax = plt.subplots()

ax.plot(d['x'], d['xi_r'].real, color=SKY_BLUE, label=r'$\tilde{\xi}_{r}/R$', zorder=1)
ax.plot(d['x'], d['xi_h'].real, color=ORANGE, label=r'$\tilde{\xi}_{\rm h}/R$', zorder=0)

ax.set_xlabel(r'$x$')

ax.set_xlim(0, 1)

ax.legend(loc=4)

fig.tight_layout()
fig.savefig('fig_disp_eigfunc.svg')
