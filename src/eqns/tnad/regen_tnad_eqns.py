#!/usr/bin/env python3
# Program : regen_nad_eqns.fypp
# Purpose : regenerate equations using sympy (non-adiabatic)
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
from regen_nad_eqns import A as A_nad

# Declare equation matrices

# Differential Jacobian matrix. This is formed from the non-adiabatic
# matrix with a turbulent correction to the radial momentum equation:
#
# x dy/dx = A_nad y - G_trb y_trb - x d/dx(y_trb)
#
# with
#
# y_trb = F_trb [x dy_1/dx + (l-1) y_1] e_2
#
# where e_2 is the basis vector associated with the radial momentum
# equation, and F_trb = i omega nu_trb / (g r)

e_1 = sp.eye(6)[:,0]
e_2 = sp.eye(6)[:,1]

P = (-A_nad + (As + V_g - U(x) + 1 - l_i)*sp.eye(6)) @ e_2
Q = -F_trb/(F_trb*A_nad[0,1] + 1) * e_1.T @ (A_nad + (l_i-1)*sp.eye(6))

A = A_nad - P @ Q

# Match condition matrix

C = sp.Matrix([
    [1, 0, 0, 0, 0, 0],
    [-U(x), U(x), 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [U(x), 0, 0, 1, 0, 0],
    [-V*nabla_ad, V*nabla_ad, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1]
])

# Main program

if __name__ == '__main__':

    # Define variable transformation matrices

    T_gyre = sp.eye(6)

    # Regenerate equation include files

    for vars, T in zip(('gyre', ), (T_gyre, )):

        with open(f'{vars}/A.inc', 'w') as f:
            f.write(generate_A(A, T)+'\n')

        with open(f'{vars}/C.inc', 'w') as f:
            f.write(generate_C(C, T)+'\n')

        with open(f'{vars}/Q.inc', 'w') as f:
            f.write(generate(Q.row(0), 'Q')+'\n')
