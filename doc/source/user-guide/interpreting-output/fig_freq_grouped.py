#!/usr/bin/env python
#
# Frequency vs radial order (grouped)

# Imports

import matplotlib.pyplot as plt
import pygyre as pg

from colors import *

# Plot settings

plt.style.use('web.mplstyle')

# Read data

s = pg.read_output('summary.h5')

sg = s.group_by('l')

# Do the plot

fig, ax = plt.subplots()

ax.plot(sg.groups[0]['n_pg'], sg.groups[0]['freq'].real, color=SKY_BLUE, label=r'$\ell=1$')
ax.plot(sg.groups[1]['n_pg'], sg.groups[1]['freq'].real, color=ORANGE, label=r'$\ell=2$')

ax.set_xlabel(r'$n_{\rm pg}$')
ax.set_ylabel(r'$\nu\,[{\rm cyc\,d^{-1}}]$')

ax.xaxis.get_major_locator().set_params(integer=True)

ax.legend()

fig.tight_layout()
fig.savefig('fig_freq_grouped.svg')
