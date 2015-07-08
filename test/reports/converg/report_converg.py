#!/usr/bin/env python
#
# Program : report_converg.py
# Purpose : Prepare a convergence report for GYRE

import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import re

import gyre

from colors import *

from matplotlib.backends.backend_pdf import PdfPages
 
# Run parameters

IVP_SOLVER = ['MAGNUS_GL2', 
              'MAGNUS_GL4',
              'MAGNUS_GL6',
              'COLLOC_GL2',
              'COLLOC_GL4',
              'COLLOC_GL6']

VARS_SET = ['DZIEM',
            'JCD',
            'MIX',
            'LAGP']

N_GRID_MIN = 10
N_GRID_MAX = 10000
N_N_GRID = 31

N_GRID = (10.**np.linspace(np.log10(N_GRID_MIN), np.log10(N_GRID_MAX), N_N_GRID)).astype(int)

L = 0
N_PG = 1

### PROCS ###

def run_gyre (n_grid, l, ivp_solver, vars_set, summ_file) :

    # Create the input file

    with open('gyre_ad.in', 'w') as fout :
        fout.write('''
&constants
/
&model
	model_type = 'HOM'
/
&mode
	l = {l:d}
/
&osc
        variables_set = '{vars_set:s}'
/
&num
	ivp_solver = '{ivp_solver:s}'
/
&scan
        grid_type = 'LINEAR'
        freq_min = 0.5
        freq_max = 20
        n_freq = 100
/
&shoot_grid
	op_type = 'CREATE_GEOM'
        n = {n_grid:d}
        s = 1000
/
&recon_grid
/
&output
	summary_file = '{summ_file:s}'
	summary_item_list = 'l,n_pg,n_p,n_g,omega'
/
'''.format(l=l, vars_set=vars_set, ivp_solver=ivp_solver, n_grid=n_grid, summ_file=summ_file))

    # Run gyre
         
    os.system('./gyre_ad gyre_ad.in > /dev/null')
#    os.remove('gyre_ad.in')

#

def read_summary (summ_file, n_pg) :

    # Read the summary file

    d, r = gyre.read_output(summ_file)

    omega_n = dict(zip(r['n_pg'], r['omega']))

    return omega_n[n_pg]

#

def eval_omega_ana (l, n_pg) :

    # Evaluate analytical frequencies

    n_ana = n_pg - 1

    Gamma_1 = 5./3.

    D_n = 0.5*(-4. + Gamma_1*(n_ana*(2*l + 2*n_ana + 5) + 2*l + 3))

    omega_p = np.sqrt(D_n + np.sqrt(D_n**2 + l*(l+1)))
    omega_f = np.sqrt(2.*l*(l-1.)/(2.*l+1.))

    omega_ana = np.where(n_ana >= 0, omega_p, omega_f)

    return omega_ana

#

def plot_defaults () :

    # Calculate figure dimensions

    fig_width = 7
    fig_height = 10

    # Set up display environment

    params = {'backend': 'pdf',
              'figure.figsize': [fig_width, fig_height],
              'font.family':'serif',
              'font.size':12,
              'font.serif': 'Times Roman',
              'legend.fontsize': 12,
              'legend.frameon' : False,
              'text.usetex': True,
              'figure.dpi': 600,
              'lines.markersize': 4,
              'lines.linewidth': 1,
              'lines.antialiased': False
              }

    mpl.rcParams.update(params)

#

def plot_data (pp, n_grid, domega, title) :

    # Plot the data

    fig, ax = plt.subplots()

    ax.plot(n_grid, np.abs(domega), color=ORANGE)
    ax.plot(n_grid, (N_GRID_MIN/n_grid.astype(float))**2, color=BLACK, ls='--', label=r'$\propto n^{-2}$')
    ax.plot(n_grid, (N_GRID_MIN/n_grid.astype(float))**4, color=BLACK, ls=':', label=r'$\propto n^{-4}$')
    ax.plot(n_grid, (N_GRID_MIN/n_grid.astype(float))**6, color=BLACK, ls='-.', label=r'$\propto n^{-6}$')

    ax.set_ylim(1E-16,10)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(r'$n$')
    ax.set_ylabel(r'$|\Delta \omega|$')

    ax.set_title(title)

    ax.legend(loc=1)

    pp.savefig(fig)

    plt.close(fig)
    
#                   

def tex_escape (text) :

    conv = {
        '&': r'\&',
        '%': r'\%',
        '$': r'\$',
        '#': r'\#',
        '_': r'\_',
        '{': r'\{',
        '}': r'\}',
        '~': r'\textasciitilde{}',
        '^': r'\^{}',
        '\\': r'\textbackslash{}',
        '<': r'\textless',
        '>': r'\textgreater',
    }

    pattern = re.compile('|'.join(re.escape(key) for key in conv.keys()))

    return pattern.sub(lambda x: conv[x.group()], text)

#

### MAIN ###

if __name__ == '__main__' :

    plot_defaults()

    pp = PdfPages('report_converg.pdf')

    # Loop through run parameters

    for ivp_solver in IVP_SOLVER :
        for vars_set in VARS_SET :

            # Loop over grid sizes

            n_grid = []
            domega = []

            for i in range(0, len(N_GRID)) :

                # Run gyre

                run_gyre(N_GRID[i], L, ivp_solver, vars_set, 'summary.h5')

                # Read the summary results

                if os.path.isfile('summary.h5') : 

                    omega = read_summary('summary.h5', N_PG).real
                    os.remove('summary.h5')

                    # Calculate the frequency error

                    n_grid.append(N_GRID[i])
                    domega.append(np.abs(omega - eval_omega_ana(L, N_PG)))

            # Plot the results
                    
            title = r'{{\tt ivp\_solver = {ivp_solver:s}}}, {{\tt variables\_set = {vars_set:s}}}'.format(ivp_solver=tex_escape(ivp_solver),
                                                                                                          vars_set=tex_escape(vars_set))

            plot_data(pp, np.array(n_grid), np.array(domega), title)

    # Finish

    pp.close()

#
