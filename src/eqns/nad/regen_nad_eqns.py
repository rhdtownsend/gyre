#!/usr/bin/env python3
# Program : regen_nad_eqns.py
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

# Declare equation matrices

# Differential Jacobian matrix

A = sp.Matrix([
    [V_g - 1 - l_i,
     lamda/(c_1(x)*alpha_omg*omega_c**2) - V_g*alpha_gam,
     alpha_grv*lamda/(c_1(x)*alpha_omg*omega_c**2),
     0,
     ups_T,
     0],
    [c_1(x)*alpha_omg*omega_c**2 - As_iso,
     As - U(x) + 3 - l_i,
     0,
     -alpha_grv,
     ups_T,
     0],
    [0,
     0,
     alpha_grv*(3 - U(x) - l_i),
     alpha_grv,
     0,
     0],
    [alpha_grv*U(x)*As,
     alpha_grv*U(x)*V_g,
     alpha_grv*lamda,
     alpha_grv*(-U(x) - l_i + 2),
     -alpha_grv*U(x)*ups_T,
     0],
    [V*(nabla_ad*(U(x) - c_1(x)*alpha_omg*omega_c**2) - 4*(nabla_ad - nabla) + c_kap_ad*V*nabla + c_dif)/f_rht,
     V*(lamda/(c_1(x)*alpha_omg*omega_c**2)*(nabla_ad - nabla) - c_kap_ad*V*nabla - c_dif)/f_rht,
     alpha_grv*(V*lamda/(c_1(x)*alpha_omg*omega_c**2)*(nabla_ad - nabla))/f_rht,
     alpha_grv*(V*nabla_ad)/f_rht,
     V*nabla*(4*f_rht - c_kap_S)/f_rht - df_rht - (l_i - 2),
     -V*nabla/(c_rad*f_rht)],
    [alpha_hfl*lamda*(nabla_ad/nabla - 1)*c_rad - V*c_eps_ad - alpha_egv*c_egv*nabla_ad*V,
     V*c_eps_ad - lamda*c_rad*alpha_hfl*nabla_ad/nabla + lamda*f_conv/(c_1(x)*alpha_omg*omega_c**2) + alpha_egv*c_egv*nabla_ad*V,
     alpha_grv*lamda*f_conv/(c_1(x)*alpha_omg*omega_c**2),
     0,
     c_eps_S - alpha_hfl*lamda*c_rad/(nabla*V) + alpha_thm*i_omega_c*c_thk + alpha_egv*c_egv,
     -1 - l_i]
])

# Match condition matrix

C = sp.Matrix([
    [1, 0, 0, 0, 0, 0],
    [-U(x), U(x), 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [U(x), 0, 0, 1, 0, 0],
    [-V*nabla_ad, V*nabla_ad, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1]
])

# Inner boundary condition matrices

IB_regular = sp.Matrix([
    [c_1(x)*alpha_omg*omega_c**2, -l_i, -alpha_grv*l_i, 0, 0, 0],
    [0, 0, alpha_grv*l_i, 1-2*alpha_grv, 0, 0],
    [0, 0, 0, 0, 1, 0]
])

IB_zero_r = sp.Matrix([
    [1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0]
])

IB_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0]
])

# Outer boundary condition matrices

OB_vacuum = sp.Matrix([
    [1, -1, 0, 0, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

OB_zero_r = sp.Matrix([
    [1, 0, 0, 0, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

OB_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

OB_dziem = sp.Matrix([
    [1 + (lamda/(c_1(x)*alpha_omg*omega_c**2) - 4 - c_1(x)*alpha_omg*omega_c**2)/V,
     -1,
     alpha_grv*(lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)/V,
     0,
     0,
     0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

OB_decomp = sp.Matrix([
    [-(chi-a_11), a_12, -alpha_grv*G_1, alpha_grv*G_2, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

OB_jcd = sp.Matrix([
    [chi-b_11,
     -b_12,
     alpha_grv*((lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)*b_12/(V_g + As)),
     0,
     0,
     0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

# Main program

if __name__ == '__main__':

    # Define variable transformation matrices

    T_gyre = sp.eye(6)

    T_dziem = sp.Matrix([
        [1, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1]
    ])

    T_jcd = sp.Matrix([
        [1, 0, 0, 0, 0, 0],
        [0, lamda/(c_1(x)*alpha_omg*omega_c**2), lamda/(c_1(x)*alpha_omg*omega_c**2), 0, 0, 0],
        [0, 0, -1, 0, 0, 0],
        [0, 0, U(x)-1, -1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1]
    ])

    T_rjcd = sp.Matrix([
        [1, 0, 0, 0, 0, 0],
        [0, 1/(c_1(x)*alpha_omg*omega_c**2), 1/(c_1(x)*alpha_omg*omega_c**2), 0, 0, 0],
        [0, 0, -1, 0, 0, 0],
        [0, 0, U(x)-1, -1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1]
    ])

    T_mix = sp.Matrix([
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [0, 0, -1, 0, 0, 0],
        [0, 0, U(x)-1, -1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1]
    ])

    T_lagp = sp.Matrix([
        [1, 0, 0, 0, 0, 0],
        [-V_2(x), V_2(x), 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1]
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

        with open(f'{vars}/C.inc', 'w') as f:
            f.write(generate_C(C, T)+'\n')

        with open(f'{vars}/R.inc', 'w') as f:
            f.write(generate_R(T)+'\n')
