#!/usr/bin/env python3
# Program : regen_ad_eqns.py
# Purpose : regenerate equations using sympy (static adiabatic)
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
    [3 - U(x) - l_e, 1],
    [lamda - As*U(x) - U(x)*V_g, -U(x) - l_e + 2]
])

# Match condition matrix

C = sp.Matrix([
    [1, 0],
    [-U(x), 1]
])

# Inner boundary condition matrices

IB_regular = sp.Matrix([
    [l_e, -1]
])

# Outer boundary condition matrices

OB_vacuum = sp.Matrix([
    [l_e + 1 - U(x), 1]
])

# Main program

if __name__ == '__main__':

    # Define variable transformation matrices

    T_gyre = sp.eye(2)

    # Regenerate equation include files

    for vars, T in zip(('gyre', ), (T_gyre, )):

        with open(f'{vars}/A.inc', 'w') as f:
            f.write(generate_A(A, T)+'\n')

        with open(f'{vars}/IB_regular.inc', 'w') as f:
            f.write(generate_IB(IB_regular, T)+'\n')

        with open(f'{vars}/OB_vacuum.inc', 'w') as f:
            f.write(generate_OB(OB_vacuum, T)+'\n')

        with open(f'{vars}/C.inc', 'w') as f:
            f.write(generate_C(C, T)+'\n')

        with open(f'{vars}/R.inc', 'w') as f:
            f.write(generate_R(T)+'\n')
