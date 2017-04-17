#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# Read a discriminant map

def read_map (filename):
    
    import numpy as np
    import h5py as h5
    
    f = h5.File(filename, 'r')
    
    D_e = f['discrim_map_e'][...]
    D_f = f['discrim_map_f'][...]
    
    log_D = 0.5*np.log10(D_f['re']**2 + D_f['im']**2) + D_e*np.log10(2.)
    
    D_re = D_f['re']*2.**(D_e-np.min(D_e))
    D_im = D_f['im']*2.**(D_e-np.min(D_e))
    
    omega_re = f['omega_re'][...]
    omega_im = f['omega_im'][...]
    
    f.close()
    
    return {'omega_re': omega_re,
            'omega_im': omega_im,
            'D_re': D_re,
            'D_im': D_im,
            'log_D': log_D}

# Calculate figure dimensions

fig_width = 5
fig_height = 4

# Set up display environment

params = {'backend': 'pdf',
          'figure.figsize': [fig_width, fig_height],
          'font.family':'serif',
          'font.size':11,
          'font.serif': 'Times Roman',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'medium',
          'legend.fontsize': 11,
          'legend.frameon' : False,
          'text.usetex': True,
          'figure.dpi': 600,
          'lines.markersize': 4,
          'lines.linewidth': 1,
          'lines.antialiased': False,
          'path.simplify': False }

mpl.rcParams.update(params)

#

D = read_map('map.h5')

fig, ax = plt.subplots()

c_re = ax.contour(D['omega_re'], D['omega_im'], D['D_re'], colors='blue', levels=[0])
c_im = ax.contour(D['omega_re'], D['omega_im'], D['D_im'], colors='red', levels=[0])

c_re.collections[0].set_label(r'${\rm Re}(\mathcal{D})=0$')
c_im.collections[0].set_label(r'${\rm Im}(\mathcal{D})=0$')

#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.02))

ax.set_xlabel(r'${\rm Re}(\omega)$')
ax.set_ylabel(r'${\rm Im}(\omega)$')

ax.legend(loc=2, frameon=True)

fig.tight_layout(pad=0.2)
fig.savefig('map.pdf')
