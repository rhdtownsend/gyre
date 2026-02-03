#!/usr/bin/env python3
# Program : regen_ad_eqns.py
# Purpose : regenerate equations using sympy (adiabatic)
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

# A_subs = []

# for sym in (V_g, As, sp.Symbol('U'), sp.Symbol('c_1')):
#     A_subs += [(sym, fcode_symbol(sym.name, f'this%{sym.name}'))]

# for sym in (alpha_grv, alpha_omg):
#     A_subs += [(sym, fcode_symbol(sym.name, f'this%cx%{sym.name}'))]

# Match condition matrix

C = sp.Matrix([
    [1, 0, 0, 0],
    [-U(x), U(x), 0, 0],
    [0, 0, 1, 0],
    [U(x), 0, 0, 1]
])

# Inner boundary condition matrices

IB_regular = sp.Matrix([
    [c_1(x)*alpha_omg*omega_c**2, -l_i, -alpha_grv*l_i, 0],
    [0, 0, alpha_grv*l_i, 1-2*alpha_grv]
])

IB_zero_r = sp.Matrix([
    [1, 0, 0, 0],
    [0, 0, 0, 1]
])

IB_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0],
    [0, 0, 0, 1]
])

# Outer boundary condition matrices

OB_vacuum = sp.Matrix([
    [1, -1, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_zero_r = sp.Matrix([
    [1, 0, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_dziem = sp.Matrix([
    [1 + (lamda/(c_1(x)*alpha_omg*omega_c**2) - 4 - c_1(x)*alpha_omg*omega_c**2)/V,
     -1,
     alpha_grv*(lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)/V,
     0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_decomp = sp.Matrix([
    [-(chi-a_11), a_12, -alpha_grv*G_1, alpha_grv*G_2],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_jcd = sp.Matrix([
    [chi-b_11,
     -b_12,
     alpha_grv*((lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)*b_12/(V_g + As)),
     0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_gamma1 = sp.Matrix([
    [1, 0, 0, 0],
    [0, 1, alpha_grv, 0]
])

OB_gamma2 = sp.Matrix([
    [V_g+l_e, lamda/(c_1(x)*alpha_omg*omega_c**2), lamda/(c_1(x)*alpha_omg*omega_c**2), 0],
    [0, 0, alpha_grv*l_e+1, 1]
])

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

        with open(f'{vars}/A.inc', 'w') as f:
            f.write(generate_A(A, T)+'\n')

        with open(f'{vars}/IB_regular.inc', 'w') as f:
            f.write(generate_IB(IB_regular, T)+'\n')

        with open(f'{vars}/IB_zero_r.inc', 'w') as f:
            f.write(generate_IB(IB_zero_r, T)+'\n')

        with open(f'{vars}/IB_zero_h.inc', 'w') as f:
            f.write(generate_IB(IB_zero_h, T)+'\n')

        with open(f'{vars}/OB_vacuum.inc', 'w') as f:
            f.write(generate_OB(OB_vacuum, T)+'\n')

        with open(f'{vars}/OB_zero_r.inc', 'w') as f:
            f.write(generate_OB(OB_zero_r, T)+'\n')

        with open(f'{vars}/OB_zero_h.inc', 'w') as f:
            f.write(generate_OB(OB_zero_h, T)+'\n')

        with open(f'{vars}/OB_dziem.inc', 'w') as f:
            f.write(generate_OB(OB_dziem, T)+'\n')

        with open(f'{vars}/OB_decomp.inc', 'w') as f:
            f.write(generate_OB(OB_decomp, T)+'\n')

        with open(f'{vars}/OB_jcd.inc', 'w') as f:
            f.write(generate_OB(OB_jcd, T)+'\n')

        with open(f'{vars}/OB_gamma1.inc', 'w') as f:
            f.write(generate_OB(OB_gamma1, T)+'\n')

        with open(f'{vars}/OB_gamma2.inc', 'w') as f:
            f.write(generate_OB(OB_gamma2, T)+'\n')

        with open(f'{vars}/C.inc', 'w') as f:
            f.write(generate_C(C, T)+'\n')

        with open(f'{vars}/R.inc', 'w') as f:
            f.write(generate_R(T)+'\n')
