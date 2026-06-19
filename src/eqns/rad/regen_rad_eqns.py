#!/usr/bin/env python3
# Program : regen_ad_eqns.py
# Purpose : regenerate equations using sympy (radial adiabatic)
#
# Copyright 2013-2026 Rich Townsend & The GYRE Team
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

import sympy as sp
import sympy.printing.fortran as spf

from symbols import *

# Routines to generate inhomogeneous terms

def generate_tidal_r(A):

    r = sp.zeros(2, 1)

    return r


def generate_tidal_s(B):

    s = sp.zeros(1, 1)

    return s


def generate_tidal_t(C):

    t = sp.zeros(2, 1)

    return t


# Declare equation matrices

# Differential Jacobian matrix and inhomogeneous vector

A = sp.Matrix([
    [V_g - 1, -V_g],
    [c_1(x)*alpha_omg*omega_c**2 + U(x) - As, As - U(x) + 3]
])

r = generate_tidal_r(A)

# Inner boundary condition matrices and inhomogeneous vectors

IB_regular = sp.Matrix([
    [c_1(x)*alpha_omg*omega_c**2, 0]
])

Is_regular = generate_tidal_s(IB_regular)

IB_zero_r = sp.Matrix([
    [1, 0]
])

Is_zero_r = generate_tidal_s(IB_zero_r)

# Outer boundary condition matrices

OB_vacuum = sp.Matrix([
    [1, -1]
])

Os_vacuum = generate_tidal_s(OB_vacuum)

OB_zero_r = sp.Matrix([
    [1, 0]
])

Os_zero_r = generate_tidal_s(OB_zero_r)

OB_dziem = sp.Matrix([
    [1 - (4 + c_1(x)*alpha_omg*omega_c**2)/V, -1]
])

Os_dziem = generate_tidal_s(OB_dziem)

OB_decomp = sp.Matrix([
    [-(chi - a_11), a_12]
])

Os_decomp = generate_tidal_s(OB_decomp)

OB_jcd = sp.Matrix([
    [chi - b_11, -b_12]
])

Os_jcd = generate_tidal_s(OB_jcd)

# Match condition matrix and inhomogeneous vector

C = sp.Matrix([
    [1, 0],
    [-U(x), U(x)]
])

t = generate_tidal_t(C)

# Main program

if __name__ == '__main__':

    # Define variable transformation matrices

    T_gyre = sp.eye(2)

    T_dziem = sp.eye(2)

    T_jcd = sp.Matrix([
        [1, 0],
        [0, 1/(c_1(x)*alpha_omg*omega_c**2)]
    ])

    T_mix = sp.eye(2)

    T_lagp = sp.Matrix([
        [1, 0],
        [-V_2(x), V_2(x)]
    ])

    # Regenerate equation include files

    for vars, T in zip(('gyre', 'dziem', 'jcd', 'mix', 'lagp'), (T_gyre, T_dziem, T_jcd, T_mix, T_lagp)):

        with open(f'{vars}/A.inc', 'w') as file:
            file.write(generate_A(A, r, T)+'\n')

        with open(f'{vars}/IB_regular.inc', 'w') as file:
            file.write(generate_IB(IB_regular, Is_regular, T)+'\n')

        with open(f'{vars}/IB_zero_r.inc', 'w') as file:
            file.write(generate_IB(IB_zero_r, Is_zero_r, T)+'\n')

        with open(f'{vars}/OB_vacuum.inc', 'w') as file:
            file.write(generate_OB(OB_vacuum, Os_vacuum, T)+'\n')

        with open(f'{vars}/OB_zero_r.inc', 'w') as file:
            file.write(generate_OB(OB_zero_r, Os_zero_r, T)+'\n')

        with open(f'{vars}/OB_dziem.inc', 'w') as file:
            file.write(generate_OB(OB_dziem, Os_dziem, T)+'\n')

        with open(f'{vars}/OB_decomp.inc', 'w') as file:
            file.write(generate_OB(OB_decomp, Os_decomp, T)+'\n')

        with open(f'{vars}/OB_jcd.inc', 'w') as file:
            file.write(generate_OB(OB_jcd, Os_jcd, T)+'\n')

        with open(f'{vars}/C.inc', 'w') as file:
            file.write(generate_C(C, t, T)+'\n')

        with open(f'{vars}/R.inc', 'w') as file:
            file.write(generate_R(T)+'\n')
