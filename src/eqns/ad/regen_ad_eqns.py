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

def generate_tidal_f(A):

    f = (A[:,2] + l_e*A[:,3])*y_T_1*c_1(x)
    f[2] = 0
    f[3] = 0

    f = f.subs(alpha_grv, 1)

    return f

def generate_tidal_g(B):

    g = -(B[:,2] + l_e*B[:,3])*y_T_1*c_1(x)
    g[1] = 0

    g = g.subs(alpha_grv, 1)

    return g

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

f = generate_tidal_f(A)

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

# Inner boundary condition matrices and inhomogeneous vectors

IB_regular = sp.Matrix([
    [c_1(x)*alpha_omg*omega_c**2, -l_i, -alpha_grv*l_i, 0],
    [0, 0, alpha_grv*l_i, 1-2*alpha_grv]
])

IB_regular_g = generate_tidal_g(IB_regular)

IB_zero_r = sp.Matrix([
    [1, 0, 0, 0],
    [0, 0, 0, 1]
])

IB_zero_r_g = generate_tidal_g(IB_zero_r)

IB_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0],
    [0, 0, 0, 1]
])

IB_zero_h_g = generate_tidal_g(IB_zero_h)

# Outer boundary condition matrices

OB_vacuum = sp.Matrix([
    [1, -1, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_vacuum_g = generate_tidal_g(OB_vacuum)

OB_zero_r = sp.Matrix([
    [1, 0, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_zero_r_g = generate_tidal_g(OB_zero_r)

OB_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_zero_h_g = generate_tidal_g(OB_zero_h)

OB_dziem = sp.Matrix([
    [1 + (lamda/(c_1(x)*alpha_omg*omega_c**2) - 4 - c_1(x)*alpha_omg*omega_c**2)/V,
     -1,
     alpha_grv*(lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)/V,
     0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_dziem_g = generate_tidal_g(OB_dziem)

OB_decomp = sp.Matrix([
    [-(chi-a_11), a_12, -alpha_grv*G_1, alpha_grv*G_2],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_decomp_g = generate_tidal_g(OB_decomp)

OB_jcd = sp.Matrix([
    [chi-b_11,
     -b_12,
     alpha_grv*((lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)*b_12/(V_g + As)),
     0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e+1, alpha_grv]
])

OB_jcd_g = generate_tidal_g(OB_jcd)

OB_gamma1 = sp.Matrix([
    [1, 0, 0, 0],
    [0, 1, alpha_grv, 0]
])

OB_gamma1_g = generate_tidal_g(OB_gamma1)

OB_gamma2 = sp.Matrix([
    [V_g+l_e, lamda/(c_1(x)*alpha_omg*omega_c**2), lamda/(c_1(x)*alpha_omg*omega_c**2), 0],
    [0, 0, alpha_grv*l_e+1, 1]
])

OB_gamma2_g = generate_tidal_g(OB_gamma2)

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

        with open(f'{vars}/A_t.inc', 'w') as file:
            file.write(generate_E(A, f, T, transpose=True)+'\n')

        with open(f'{vars}/IB_regular.inc', 'w') as file:
            file.write(generate_IB(IB_regular, IB_regular_g, T)+'\n')

        with open(f'{vars}/IB_zero_r.inc', 'w') as file:
            file.write(generate_IB(IB_zero_r, IB_zero_r_g, T)+'\n')

        with open(f'{vars}/IB_zero_h.inc', 'w') as file:
            file.write(generate_IB(IB_zero_h, IB_zero_h_g, T)+'\n')

        with open(f'{vars}/OB_vacuum.inc', 'w') as file:
            file.write(generate_OB(OB_vacuum, OB_vacuum_g, T)+'\n')

        with open(f'{vars}/OB_zero_r.inc', 'w') as file:
            file.write(generate_OB(OB_zero_r, OB_zero_r_g, T)+'\n')

        with open(f'{vars}/OB_zero_h.inc', 'w') as file:
            file.write(generate_OB(OB_zero_h, OB_zero_h_g, T)+'\n')

        with open(f'{vars}/OB_dziem.inc', 'w') as file:
            file.write(generate_OB(OB_dziem, OB_dziem_g, T)+'\n')

        with open(f'{vars}/OB_decomp.inc', 'w') as file:
            file.write(generate_OB(OB_decomp, OB_decomp_g, T)+'\n')

        with open(f'{vars}/OB_jcd.inc', 'w') as file:
            file.write(generate_OB(OB_jcd, OB_jcd_g, T)+'\n')

        with open(f'{vars}/OB_gamma1.inc', 'w') as file:
            file.write(generate_OB(OB_gamma1, OB_gamma1_g, T)+'\n')

        with open(f'{vars}/OB_gamma2.inc', 'w') as file:
            file.write(generate_OB(OB_gamma2, OB_gamma2_g, T)+'\n')

        with open(f'{vars}/C_t.inc', 'w') as file:
            file.write(generate_C(C, T, transpose=True)+'\n')

        with open(f'{vars}/R.inc', 'w') as file:
            file.write(generate_R(T)+'\n')
