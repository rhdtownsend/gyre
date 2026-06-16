#!/usr/bin/env python3
# Program : regen_ad_eqns.py
# Purpose : regenerate equations using sympy (adiabatic)
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

def generate_tidal_A_F(A):

    A_F = (A[:,2] + l_e*A[:,3])*y_T_1*c_1(x)
    A_F[2] = 0
    A_F[3] = 0

    A_F = A_F.subs(alpha_grv, 1)

    return A_F


def generate_tidal_B_F(B):

    B_F = -(B[:,2] + l_e*B[:,3])*y_T_1*c_1(x)
    B_F[1] = 0

    B_F = B_F.subs(alpha_grv, 1)

    return B_F


def generate_tidal_C_F(C):

    C_F = sp.zeros(4, 1)

    return C_F


# Declare equation matrices

# Differential Jacobian matrix and inhomogeneous vector

A = sp.Matrix([
    [V_g - 1 - l_i,
     lamda/(c_1(x)*alpha_omg*omega_c**2) - V_g*alpha_gam,
     alpha_grv*lamda/(c_1(x)*alpha_omg*omega_c**2),
     0],
    [c_1(x)*alpha_omg*omega_c**2 - As_iso,
     As - U(x) + 3 - l_i,
     0,
     -alpha_grv],
    [0,
     0,
     alpha_grv*(3 - U(x) - l_i),
     alpha_grv],
    [alpha_grv*U(x)*As,
     alpha_grv*U(x)*V_g,
     alpha_grv*lamda,
     alpha_grv*(-U(x) - l_i + 2)]
])

A_F = generate_tidal_A_F(A)

# Match condition matrix and inhomogeneous vector

C = sp.Matrix([
    [1, 0, 0, 0],
    [-U(x), U(x), 0, 0],
    [0, 0, 1, 0],
    [U(x), 0, 0, 1]
])

C_F = generate_tidal_C_F(C)

# Inner boundary condition matrices and inhomogeneous vectors

IB_regular = sp.Matrix([
    [c_1(x)*alpha_omg*omega_c**2, -l_i, -alpha_grv*l_i, 0],
    [0, 0, alpha_grv*l_i, 1-2*alpha_grv]
])

IB_F_regular = generate_tidal_B_F(IB_regular)

IB_zero_r = sp.Matrix([
    [1, 0, 0, 0],
    [0, 0, 0, 1]
])

IB_F_zero_r = generate_tidal_B_F(IB_zero_r)

IB_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0],
    [0, 0, 0, 1]
])

IB_F_zero_h = generate_tidal_B_F(IB_zero_h)

# Outer boundary condition matrices

OB_vacuum = sp.Matrix([
    [1, -1, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_F_vacuum = generate_tidal_B_F(OB_vacuum)

OB_zero_r = sp.Matrix([
    [1, 0, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_F_zero_r = generate_tidal_B_F(OB_zero_r)

OB_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_F_zero_h = generate_tidal_B_F(OB_zero_h)

OB_dziem = sp.Matrix([
    [1 + (lamda/(c_1(x)*alpha_omg*omega_c**2) - 4 - c_1(x)*alpha_omg*omega_c**2)/V,
     -1,
     alpha_grv*(lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)/V,
     0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_F_dziem = generate_tidal_B_F(OB_dziem)

OB_decomp = sp.Matrix([
    [-(chi-a_11), a_12, -alpha_grv*G_1, alpha_grv*G_2],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_F_decomp = generate_tidal_B_F(OB_decomp)

OB_jcd = sp.Matrix([
    [chi-b_11,
     -b_12,
     alpha_grv*((lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)*b_12/(V_g + As)),
     0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_F_jcd = generate_tidal_B_F(OB_jcd)

OB_gamma1 = sp.Matrix([
    [1, 0, 0, 0],
    [0, 1, alpha_grv, 0]
])

OB_F_gamma1 = generate_tidal_B_F(OB_gamma1)

OB_gamma2 = sp.Matrix([
    [V_g+l_e, lamda/(c_1(x)*alpha_omg*omega_c**2), lamda/(c_1(x)*alpha_omg*omega_c**2), 0],
    [0, 0, alpha_grv*l_e+1, 1]
])

OB_F_gamma2 = generate_tidal_B_F(OB_gamma2)

# Main program

if __name__ == '__main__':

    # Define variable transformation matrices

    T_gyre = sp.eye(4)

    T_dziem = sp.Matrix([
        [1, 0, 0, 0],
        [0, 1, 1, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])

    T_jcd = sp.Matrix([
        [1, 0, 0, 0],
        [0, lamda/(c_1(x)*alpha_omg*omega_c**2), lamda/(c_1(x)*alpha_omg*omega_c**2), 0],
        [0, 0, -1, 0],
        [0, 0, U(x)-1, -1]
    ])

    T_rjcd = sp.Matrix([
        [1, 0, 0, 0],
        [0, 1/(c_1(x)*alpha_omg*omega_c**2), 1/(c_1(x)*alpha_omg*omega_c**2), 0],
        [0, 0, -1, 0],
        [0, 0, U(x)-1, -1]
    ])

    T_mix = sp.Matrix([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, -1, 0],
        [0, 0, U(x)-1, -1]
    ])

    T_lagp = sp.Matrix([
        [1, 0, 0, 0],
        [-V_2(x), V_2(x), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])

    # Regenerate equation include files

    for vars, T in zip(('gyre', 'dziem', 'jcd', 'rjcd', 'mix', 'lagp'), (T_gyre, T_dziem, T_jcd, T_rjcd, T_mix, T_lagp)):

        with open(f'{vars}/A.inc', 'w') as file:
            file.write(generate_A(A, A_F, T)+'\n')

        with open(f'{vars}/C.inc', 'w') as file:
            file.write(generate_C(C, C_F, T)+'\n')

        with open(f'{vars}/IB_regular.inc', 'w') as file:
            file.write(generate_IB(IB_regular, IB_F_regular, T)+'\n')

        with open(f'{vars}/IB_zero_r.inc', 'w') as file:
            file.write(generate_IB(IB_zero_r, IB_F_zero_r, T)+'\n')

        with open(f'{vars}/IB_zero_h.inc', 'w') as file:
            file.write(generate_IB(IB_zero_h, IB_F_zero_h, T)+'\n')

        with open(f'{vars}/OB_vacuum.inc', 'w') as file:
            file.write(generate_OB(OB_vacuum, OB_F_vacuum, T)+'\n')

        with open(f'{vars}/OB_zero_r.inc', 'w') as file:
            file.write(generate_OB(OB_zero_r, OB_F_zero_r, T)+'\n')

        with open(f'{vars}/OB_zero_h.inc', 'w') as file:
            file.write(generate_OB(OB_zero_h, OB_F_zero_h, T)+'\n')

        with open(f'{vars}/OB_dziem.inc', 'w') as file:
            file.write(generate_OB(OB_dziem, OB_F_dziem, T)+'\n')

        with open(f'{vars}/OB_decomp.inc', 'w') as file:
            file.write(generate_OB(OB_decomp, OB_F_decomp, T)+'\n')

        with open(f'{vars}/OB_jcd.inc', 'w') as file:
            file.write(generate_OB(OB_jcd, OB_F_jcd, T)+'\n')

        with open(f'{vars}/OB_gamma1.inc', 'w') as file:
            file.write(generate_OB(OB_gamma1, OB_F_gamma1, T)+'\n')

        with open(f'{vars}/OB_gamma2.inc', 'w') as file:
            file.write(generate_OB(OB_gamma2, OB_F_gamma2, T)+'\n')

        with open(f'{vars}/R.inc', 'w') as file:
            file.write(generate_R(T)+'\n')
