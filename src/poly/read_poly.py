# Program  : read_poly
# Purpose  : read composite polytropes in POLY format
#
# Copyright 2020 Rich Townsend & The GYRE Team
#
# This file is part of GYRE. GYRE is free software: you can
# redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation, version 3.
#
# GYRE is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import h5py
import numpy as np

def read_poly (file) :

    # Read data from the file

    file = h5py.File(file, 'r')

    z = file['z'][...]

    theta = file['theta'][...]
    dtheta = file['dtheta'][...]

    n_poly = file.attrs['n_poly']
    Gamma_1 = file.attrs['Gamma_1']

    n_r = file.attrs['n_r']

    if n_r > 1:
        Delta_b = file.attrs['Delta_b']
    else:
        Delta_b = np.empty([0])

    # Reconstruct other data

    t = np.empty(n_r)
    B = np.empty(n_r)
    z_b = np.empty(n_r-1)

    t[0] = 1.
    B[0] = 1.

    n_z = len(z)

    n_poly_z = np.empty([n_z])
    t_z = np.empty([n_z])
    B_z = np.empty([n_z])

    n_poly_z[0] = n_poly[0]
    t_z[0] = t[0]
    B_z[0] = B[0]

    i = 0

    for k in range(1, n_z):

        if z[k] == z[k-1]:
            z_b[i] = z[k]
            i += 1
            t[i] = t[i-1]*np.exp(n_poly[i-1]*np.log(theta[k-1]) - n_poly[i]*np.log(theta[k]) + Delta_b[i-1])
            B[i] = B[i-1]*(n_poly[i-1]+1)/(n_poly[i]+1)*theta[k]**(n_poly[i]+1)/theta[k-1]**(n_poly[i-1]+1)*t[i]**2/t[i-1]**2
 
        n_poly_z[k] = n_poly[i]
        t_z[k] = t[i]
        B_z[k] = B[i]

    # Return data

    return {'z': z,
            'z_b': z_b,
            'z_s': z[-1],
            'theta': theta,
            'dtheta': dtheta,
            'n_poly': n_poly,
            'n_poly_z': n_poly_z,
            't': t,
            't_z': t_z,
            'B': B,
            'B_z': B_z,
            'Delta_b': Delta_b,
            'Gamma_1': Gamma_1, 
            'n_z': len(z),
            'n_r': n_r}
