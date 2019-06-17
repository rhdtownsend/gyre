#!/usr/bin/env python

import sys
import h5py
import numpy as np
import scipy.special as ss
import matplotlib.pyplot as plt

#

def read_response (filename):

    # Read data from gyre_response

    f = h5py.File(filename, 'r')

    xi_r_re = f['xi_r']['re'][...]
    xi_r_im = f['xi_r']['im'][...]

    lag_L_re = f['lag_L']['re'][...]
    lag_L_im = f['lag_L']['im'][...]

    k_max = f.attrs['k_max']
    l_max = f.attrs['l_max']

    Omega_rot = f.attrs['Omega_rot']
    Omega_orb = f.attrs['Omega_orb']

    f.close()

    # Return it

    return {'xi_r': xi_r_re + 1j*xi_r_im,
            'lag_L': lag_L_re + 1j*lag_L_im,
            'k_max': k_max,
            'l_max': l_max,
            'Omega_rot': Omega_rot,
            'Omega_orb': Omega_orb}

#

def eval_fourier (data, theta_0, phi_0):

    # Initialize the fourier amplitudes array

    A = np.zeros(2*d['k_max']+1, dtype=np.complex)

    # Loop over l, m and k

    for l in range(2, d['l_max']+1):
        for m in range(-l, l+1):
            for k in range(-d['k_max'], d['k_max']+1):

                i_l = l
                i_m = m + d['l_max']
                i_k = k + d['k_max']

                # Townsend (2003), eqns. 17 & 18

                Del_R = np.sqrt(4.*np.pi) * d['xi_r'][i_k,i_m,i_l]
                Del_L = np.sqrt(4.*np.pi) * d['lag_L'][i_k,i_m,i_l]
                Del_T = 0.25*(Del_L - 2*Del_R)

                # Townsend (2003), eqns. (15) & (16) assuming I_x
                # = const. (no limb darkening)

                I_0 = 1/2

                if l == 0:
                    I_l = I_0
                elif l == 1:
                    I_l = 1/3
                elif l == 2:
                    I_l = 1/8
                elif l == 3:
                    I_l = 0
                elif l == 4:
                    I_l = -1/48
                elif l == 5:
                    I_l = 0
                elif l == 6:
                    I_l = 1/128
                else:
                    raise Exception('No I_l defined')

                # Townsend (2003), eqns. (12) & (13) assuming
                # dlnI_x/dlnTeff = 4 (black body)
                
                Yml = ss.sph_harm(m, l, phi_0, theta_0)

                Rml = (2 + l)*(1 - l)*I_l/I_0 * Yml
                Tml = 4*I_l/I_0 * Yml

                # Add the Fourier contribution

                A[i_k] += Del_R*Rml + Del_T*Tml

    # Return data

    return A

####

if __name__ == "__main__":

    # Read arguments

    response_file = sys.argv[1]
    inc = float(sys.argv[2])
    omega = float(sys.argv[3])
    light_fig = sys.argv[4]
    fourier_fig = sys.argv[5]

    # Set the sub-observer co-latitude and azimuth

    theta_0 = inc/180 * np.pi
    phi_0 = (90-omega)/180 * np.pi

    # Read tidal response data

    d = read_response(response_file)

    # Calculate fourier amplitudes (frequencies in units of Omega_orb)

    freq = np.linspace(-d['k_max'], d['k_max'], 2*d['k_max']+1)

    A = eval_fourier(d, theta_0, phi_0)

    # Use amplitudes to construct light curve

    n = 1001

    phase = np.linspace(-0.2, 1.2, n)

    dF_F = np.zeros(n)

    for k in range(-d['k_max'], d['k_max']+1):

        i_k = k + d['k_max']

        # Add the flux contribution from harmonic k

        dF_F += np.real(A[i_k] * np.exp(1j*k*2.*np.pi*phase))

    # Plot the lightcurve

    fig, ax = plt.subplots()

    ax.plot(phase, 1 + dF_F)

    ax.set_xlabel('Orbital Phase')
    ax.set_ylabel(r'Relative Flux')

    ax.set_ylim(0.9, 1.1)

    fig.tight_layout()
    fig.savefig(light_fig)

    # Plot the flux amplitude spectrum

    amp = np.abs(A)

    amp_min = 1E-5
    amp_max = 1E-2

    fig, ax = plt.subplots()

    for i, f in enumerate(freq):
        ax.plot([f,f], [amp_min,amp[i]], color='g')

    ax.scatter(freq, amp, color='k')

    ax.grid()

    ax.set_xlabel('Orbital Harmonic')
    ax.set_ylabel('Relative Flux Amplitude')

    ax.set_xlim(0, np.max(freq))

    ax.set_yscale('log')
    ax.set_ylim(amp_min, amp_max)

    fig.tight_layout()
    fig.savefig(fourier_fig)
