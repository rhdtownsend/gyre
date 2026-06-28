#!/usr/bin/env python3
# Program : regen_nad_eqns.py
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

# Routines to generate inhomogeneous terms

def generate_tidal_r(A):

    r = (A[:,2] + l_e*A[:,3])*y_T_1*c_1(x)
    r[2] = 0
    r[3] = 0

    r = r.subs(alpha_grv, 1)

    return r


def generate_tidal_s(B):

    s = -(B[:,2] + l_e*B[:,3])*y_T_1*c_1(x)
    s[1] = 0

    s = s.subs(alpha_grv, 1)

    return s


def generate_tidal_t(t):

    t = sp.zeros(6, 1)

    return t


# Declare equation matrices

# Differential Jacobian matrix and inhomogeneous vector

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
     3 - U(x) - l_i,
     1,
     0,
     0],
    [U(x)*As,
     U(x)*V_g,
     lamda,
     -U(x) - l_i + 2,
     -U(x)*ups_T,
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

r = generate_tidal_r(A)

# Differential Jacobian matrix with turbulent drag term. This is
# formed from the non-adiabatic matrix with a turbulent correction to
# the radial momentum equation:
#
#  x dy/dx = A y + G y_trb - x d/dx(y_trb) + r_nad
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
#
# Introducing the intermediate matrix H so that y_trb = H y, we rewrite
# the differential equations as
#
#  x dz/dx = A_trb z + r
#
# where
#
#   z = y + y_trb
#   A_trb = (A_nad + G H) Q
#   Q = (I + H)^-1
#
# After z is found, y can be reconstructed as y = Q z

G = As + V_g - U(x) + 1 - l_i

H = F_trb * sp.Matrix([
    [0, 0, 0, 0, 0, 0],
    [A[0,0] + l_i - 1, A[0,1], A[0,2], A[0,3], A[0,4], A[0,5]],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
    ])

Q = sp.Inverse(sp.eye(6) + H)

A_trb = (A + G*H) @ Q

# Inner boundary condition matrices and inhomogeneous vectors

IB_regular = sp.Matrix([
    [c_1(x)*alpha_omg*omega_c**2, -l_i, -alpha_grv*l_i, 0, 0, 0],
    [0, 0, l_i, -1, 0, 0],
    [0, 0, 0, 0, 1, 0]
])

Is_regular = generate_tidal_s(IB_regular)

IB_zero_r = sp.Matrix([
    [1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0]
])

Is_zero_r = generate_tidal_s(IB_zero_r)

IB_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0]
])

Is_zero_h = generate_tidal_s(IB_zero_h)

# Outer boundary condition matrices and inhomogeneous vectors

OB_vacuum = sp.Matrix([
    [1, -1, 0, 0, 0, 0],
    [alpha_gbc*U(x), 0, l_e + 1, 1, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

Os_vacuum = generate_tidal_s(OB_vacuum)

OB_zero_r = sp.Matrix([
    [1, 0, 0, 0, 0, 0],
    [alpha_gbc*U(x), 0, l_e + 1, 1, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

Os_zero_r = generate_tidal_s(OB_zero_r)

OB_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0, 0, 0],
    [alpha_gbc*U(x), 0, l_e + 1, 1, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

Os_zero_h = generate_tidal_s(OB_zero_h)

OB_dziem = sp.Matrix([
    [1 + (lamda/(c_1(x)*alpha_omg*omega_c**2) - 4 - c_1(x)*alpha_omg*omega_c**2)/V,
     -1,
     alpha_grv*(lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)/V,
     0,
     0,
     0],
    [alpha_gbc*U(x), 0, l_e + 1, 1, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

Os_dziem = generate_tidal_s(OB_dziem)

OB_decomp = sp.Matrix([
    [-(chi-a_11), a_12, -alpha_grv*G_1, alpha_grv*G_2, 0, 0],
    [alpha_gbc*U(x), 0, l_e + 1, 1, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

Os_decomp = generate_tidal_s(OB_decomp)

OB_jcd = sp.Matrix([
    [chi-b_11,
     -b_12,
     alpha_grv*((lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)*b_12/(V_g + As)),
     0,
     0,
     0],
    [alpha_gbc*U(x), 0, l_e + 1, 1, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

Os_jcd = generate_tidal_s(OB_jcd)

# Match condition matrix and inhomogeneous vector

C = sp.Matrix([
    [1, 0, 0, 0, 0, 0],
    [-U(x), U(x), 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [U(x), 0, 0, 1, 0, 0],
    [-V*nabla_ad, V*nabla_ad, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1]
])

t = generate_tidal_t(C)

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
        [0, lamda/(c_1(x)*alpha_omg*omega_c**2), alpha_grv*lamda/(c_1(x)*alpha_omg*omega_c**2), 0, 0, 0],
        [0, 0, -1, 0, 0, 0],
        [0, 0, U(x)-1, -1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1]
    ])

    T_rjcd = sp.Matrix([
        [1, 0, 0, 0, 0, 0],
        [0, 1/(c_1(x)*alpha_omg*omega_c**2), alpha_grv/(c_1(x)*alpha_omg*omega_c**2), 0, 0, 0],
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

        with open(f'{vars}/A.inc', 'w') as file:
            file.write(generate_A(A, r, T)+'\n')

        with open(f'{vars}/A_trb.inc', 'w') as file:
            file.write(generate_A(A_trb, r, T)+'\n')

        with open(f'{vars}/IB_regular.inc', 'w') as file:
            file.write(generate_IB(IB_regular, Is_regular, T)+'\n')

        with open(f'{vars}/IB_zero_r.inc', 'w') as file:
            file.write(generate_IB(IB_zero_r, Is_zero_r, T)+'\n')

        with open(f'{vars}/IB_zero_h.inc', 'w') as file:
            file.write(generate_IB(IB_zero_h, Is_zero_h, T)+'\n')

        with open(f'{vars}/OB_vacuum.inc', 'w') as file:
            file.write(generate_OB(OB_vacuum, Os_vacuum, T)+'\n')

        with open(f'{vars}/OB_zero_r.inc', 'w') as file:
            file.write(generate_OB(OB_zero_r, Os_zero_r, T)+'\n')

        with open(f'{vars}/OB_zero_h.inc', 'w') as file:
            file.write(generate_OB(OB_zero_h, Os_zero_h, T)+'\n')

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

        with open(f'{vars}/Q.inc', 'w') as f:
            f.write(generate(Q, 'Q')+'\n')
