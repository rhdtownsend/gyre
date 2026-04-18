#!/usr/bin/env python3
# Program : regen_nad_eqns.fypp
# Purpose : regenerate equations using sympy (non-adiabatic)
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
from regen_nad_eqns import A as A_nad

# Declare equation matrices

# Differential Jacobian matrix. This is formed from the non-adiabatic
# matrix with a turbulent correction to the radial momentum equation:
#
#  x dy/dx = A_nad y + G y_trb - x d/dx(y_trb)
#
# with
#
#  G = As + V_g - U + 1 - l
#
# and
#
#  y_trb = F_trb [x dy_1/dx + (l-1) y_1] e_2
#
# where e_2 is the basis vector associated with the radial momentum
# equation, and F_trb = i omega nu_trb / (g r).

# H is an intermediate matrix defined so that y_trb = H y

H = F_trb * sp.Matrix([
    [0, 0, 0, 0, 0, 0],
    [A_nad[0,0] + l_i - 1, A_nad[0,1], A_nad[0,2], A_nad[0,3], A_nad[0,4], A_nad[0,5]],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
    ])

# Use A_nad, G and H to rewrite the equations as
#
#  x dz/dx = A z
#
# where z = y + y_trb. The matrix A is given by
#
#  A_nad + (G - A) H Q
#
# where Q = (I + H)^-1. Note that y can be reconstructed
# from z via y = Q z

G = As + V_g - U(x) + 1 - l_i

Q = sp.Inverse(sp.eye(6) + H)

A = A_nad + (G*sp.eye(6) - A_nad) @ H @ Q

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

        with open(f'{vars}/A_t.inc', 'w') as f:
            f.write(generate_A(A, T, transpose=True)+'\n')

        with open(f'{vars}/C_t.inc', 'w') as f:
            f.write(generate_C(C, T, transpose=True)+'\n')

        with open(f'{vars}/Q.inc', 'w') as f:
            f.write(generate(Q, 'Q')+'\n')
