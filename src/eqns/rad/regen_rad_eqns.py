#!/usr/bin/env python3
# Program : regen_ad_eqns.py
# Purpose : regenerate equations using sympy (radial adiabatic)
#
# Copyright 2013-2025 Rich Townsend & The GYRE Team
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

# Declare equation matrices

# Differential Jacobian matrix

A = sp.Matrix([
    [V_g - 1, -V_g],
    [c_1(x)*alpha_omg*omega_c**2 + U(x) - As, As - U(x) + 3]
])

# Match condition matrix

C = sp.Matrix([
    [1, 0],
    [-U(x), U(x)]
])

# Inner boundary condition matrices

B_i_regular = sp.Matrix([
    [c_1(x)*alpha_omg*omega_c**2, 0]
])

B_i_zero_r = sp.Matrix([
    [1, 0]
])

# Outer boundary condition matrices

B_o_vacuum = sp.Matrix([
    [1, -1]
])

B_o_zero_r = sp.Matrix([
    [1, 0]
])

B_o_dziem = sp.Matrix([
    [1 - (4 + c_1(x)*alpha_omg*omega_c**2)/V, -1]
])

B_o_decomp = sp.Matrix([
    [-(chi - a_11), a_12]
])

B_o_jcd = sp.Matrix([
    [chi - b_11, -b_12]
])

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

        with open(f'{vars}/A.inc', 'w') as f:
            f.write(generate_A(A, T)+'\n')

        with open(f'{vars}/B_i_regular.inc', 'w') as f:
            f.write(generate_B_i(B_i_regular, T)+'\n')

        with open(f'{vars}/B_i_zero_r.inc', 'w') as f:
            f.write(generate_B_i(B_i_zero_r, T)+'\n')

        with open(f'{vars}/B_o_vacuum.inc', 'w') as f:
            f.write(generate_B_o(B_o_vacuum, T)+'\n')

        with open(f'{vars}/B_o_zero_r.inc', 'w') as f:
            f.write(generate_B_o(B_o_zero_r, T)+'\n')

        with open(f'{vars}/B_o_dziem.inc', 'w') as f:
            f.write(generate_B_o(B_o_dziem, T)+'\n')

        with open(f'{vars}/B_o_decomp.inc', 'w') as f:
            f.write(generate_B_o(B_o_decomp, T)+'\n')

        with open(f'{vars}/B_o_jcd.inc', 'w') as f:
            f.write(generate_B_o(B_o_jcd, T)+'\n')

        with open(f'{vars}/C.inc', 'w') as f:
            f.write(generate_C(C, T)+'\n')

        with open(f'{vars}/R.inc', 'w') as f:
            f.write(generate_R(T)+'\n')
