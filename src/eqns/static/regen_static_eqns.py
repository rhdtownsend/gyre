#!/usr/bin/env python3
# Program : regen_static_eqns.py
# Purpose : regenerate equations using sympy (static-tide)
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

# Declare equation matrices

# Differential Jacobian matrix and inhomogeneous vector

A = sp.Matrix([
    [3 - U(x) - l_e, 1],
    [lamda - V_g*U(x) - As*U(x), -U(x) - l_e + 2]
])

f = sp.Matrix([
    0,
    (-V_g*U(x)-As*U(x))*y_T_1*c_1(x)
])

# Match condition matrix

C = sp.Matrix([
    [1, 0],
    [-U(x), 1]
])

# Inner boundary condition matrices and inhomogeneous vectors

IB_regular = sp.Matrix([
    [l_e, -1]
])

IB_regular_g = sp.Matrix([
    0
])

# Outer boundary condition matrices and inhomogeneous vectors

OB_vacuum = sp.Matrix([
    [l_e + 1 - U(x), 1]
])

OB_vacuum_g = sp.Matrix([
    -U(x)*y_T_1
])

# Main program

if __name__ == '__main__':

    # Define the variable transformation matrix

    T = sp.eye(2)

    # Regenerate equation include files

    with open(f'A_t.inc', 'w') as file:
        file.write(generate_E(A, f, T, transpose=True)+'\n')

    with open(f'IB_regular.inc', 'w') as file:
        file.write(generate_IB(IB_regular, IB_regular_g, T)+'\n')

    with open(f'OB_vacuum.inc', 'w') as file:
        file.write(generate_OB(OB_vacuum, OB_vacuum_g, T)+'\n')

    with open(f'C_t.inc', 'w') as file:
        file.write(generate_C(C, T, transpose=True)+'\n')
