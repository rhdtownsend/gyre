import pygyre as pg
import numpy as np

# Read summary file from gyre_tides

s = pg.read_output('summary.h5')

# Extract the first set of responses

sg = s.group_by('id').groups[0]

Omega_orb = sg['Omega_orb']               
R_a = sg['R_a']
q = sg['q']

eps_T = R_a**3*q
    
l = sg['l']
m = sg['m']
k = sg['k']

cbar = sg['cbar']
        
Fbar = -0.5*sg['eul_Phi_ref']/(cbar*eps_T)

x = sg['x_ref']

Gbar_1 = sg['Gbar_1']
Gbar_2 = sg['Gbar_2']
Gbar_3 = sg['Gbar_3']
Gbar_4 = sg['Gbar_4']

# Evaluate secular rates-of-change

kap = np.empty(len(l))

for i in range(len(kap)):
    if k[i] == 0:
        if m[i] == 0:
            kap[i] = 0.5
        elif m[i] > 0 and m[i] <= l[i]:
            kap[i] = 1.
        else:
            kap[i] = 0.
    elif k[i] > 0:
        kap[i] = 1.
    else:
        kap[i] = 0.

# Semi-major axis (units of R per dynamical timescale)

a_dot = np.sum(4. * Omega_orb * q / R_a * 
    (R_a)**(l+3) * (x)**(l+1) * kap * Fbar.imag * Gbar_2)

# Eccentricity (units of per dynamical timescale)

e_dot = np.sum(4. * Omega_orb * q * 
    (R_a)**(l+3) * (x)**(l+1) * kap * Fbar.imag * Gbar_3)

# Argument of periastron (units of radians per dynamical timescale)

pom_dot = np.sum(4. * Omega_orb * q * 
    (R_a)**(l+3) * (x)**(l+1) * kap * Fbar.real * Gbar_1)

# Angular momentum (units of GM^2/R)

J_dot = np.sun(4. * q**2 * R_a *
    (R_a)**(l+3) * (x)**(l+1) * kap * Fbar.imag * Gbar_4)
