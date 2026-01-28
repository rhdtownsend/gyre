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

# Declare structure coefficients

x = sp.Symbol('x')
invx = sp.Symbol('invx')

class V_2(sp.Function):
    @classmethod
    def _fcode(self, printer):
        return 'V_2'
    def fdiff(self, argindex=1):
        return V_2(x)*dV_2/self.args[0]

class c_1(sp.Function):
    @classmethod
    def _fcode(self, printer):
        return 'c_1'
    def fdiff(self, argindex=1):
        return c_1(self.args[0])*(3 - U(self.args[0]))/self.args[0]

class U(sp.Function):
    @classmethod
    def _fcode(self, printer):
        return 'U'
    def fdiff(self, argindex=1):
        return U(self.args[0])*(-V_g - As - U(self.args[0]) + 3)/self.args[0]

V = sp.Symbol('V')
V_g = sp.Symbol('V_g')
dV_2 = sp.Symbol('dV_2')
As = sp.Symbol('As')
As_iso = sp.Symbol('As_iso')

nabla = sp.Symbol('nabla')
nabla_ad = sp.Symbol('nabla_ad')
ups_T = sp.Symbol('ups_T')

c_kap_ad = sp.Symbol('c_kap_ad')
c_kap_S = sp.Symbol('c_kap_S')
c_eps_ad = sp.Symbol('c_eps_ad')
c_eps_S = sp.Symbol('c_eps_S')
c_dif = sp.Symbol('c_dif')
c_rad = sp.Symbol('c_rad')
c_thk = sp.Symbol('c_thk')
c_egv = sp.Symbol('c_egv')

f_rht = sp.Symbol('f_rht')
df_rht = sp.Symbol('df_rht')

conv_term = sp.Symbol('conv_term')

lamda = sp.Symbol('lambda')
l_i = sp.Symbol('l_i')
l_e = sp.Symbol('l_e')
omega_c = sp.Symbol('omega_c')
i_omega_c = sp.Symbol('i_omega_c')

eta = sp.Symbol('eta')
om2 = sp.Symbol('om2')

alpha_omg = sp.Symbol('alpha_omg')
alpha_grv = sp.Symbol('alpha_grv')
alpha_gbc = sp.Symbol('alpha_gbc')
alpha_gam = sp.Symbol('alpha_gam')
alpha_hfl = sp.Symbol('alpha_hfl')
alpha_egv = sp.Symbol('alpha_egv')
alpha_thm = sp.Symbol('alpha_thm')

chi = sp.Symbol('chi')
a_11 = sp.Symbol('a_11')
a_12 = sp.Symbol('a_12')
b_11 = sp.Symbol('b_11')
b_12 = sp.Symbol('b_12')
G_1 = sp.Symbol('G_1')
G_2 = sp.Symbol('G_2')

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
     V*c_eps_ad - lamda*c_rad*alpha_hfl*nabla_ad/nabla + conv_term + alpha_egv*c_egv*nabla_ad*V,
     alpha_grv*conv_term,
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

B_i_regular = sp.Matrix([
    [c_1(x)*alpha_omg*omega_c**2, -l_i, -alpha_grv*l_i, 0, 0, 0],
    [0, 0, alpha_grv*l_i, 1-2*alpha_grv, 0, 0],
    [0, 0, 0, 0, 1, 0]
])

B_i_zero_r = sp.Matrix([
    [1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0]
])

B_i_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0]
])

# Outer boundary condition matrices

B_o_vacuum = sp.Matrix([
    [1, -1, 0, 0, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

B_o_zero_r = sp.Matrix([
    [1, 0, 0, 0, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

B_o_zero_h = sp.Matrix([
    [0, 1, alpha_grv, 0, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

B_o_dziem = sp.Matrix([
    [1 + (lamda/(c_1(x)*alpha_omg*omega_c**2) - 4 - c_1(x)*alpha_omg*omega_c**2)/V,
     -1,
     alpha_grv*(lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)/V,
     0,
     0,
     0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

B_o_decomp = sp.Matrix([
    [-(chi-a_11), a_12, -alpha_grv*G_1, alpha_grv*G_2, 0, 0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

B_o_jcd = sp.Matrix([
    [chi-b_11,
     -b_12,
     alpha_grv*((lamda/(c_1(x)*alpha_omg*omega_c**2) - l_e - 1)*b_12/(V_g + As)),
     0,
     0,
     0],
    [alpha_grv*alpha_gbc*U(x), 0, alpha_grv*l_e + 1, alpha_grv, 0, 0],
    [2 - 4*nabla_ad*V, 4*nabla_ad*V, 0, 0, 4*f_rht, -1]
])

# Code generation routines

def generate_A(A, T):

    # Transform the Jacobian matrix

    A = T*(A*T.inv() - x*T.inv().diff(x))

    # Apply substitutions

    A = A.subs(((c_1(x)*alpha_omg*omega_c**2, om2), (lamda/om2, eta)))

    # Convert to Fortran, with the invx scaling

    return spf.fcode(invx*A, assign_to='A', standard=2008, source_format='free')

def generate_B_i(B_i, T):

    # Transform the inner boundary condition matrix

    B_i = B_i*T.inv()

    # Apply substitutions

    B_i = B_i.subs(((V, 0), (V_g, 0), (U(x), 3), (As, 0), (As_iso, 0)))

    # Convert to Fortran

    return spf.fcode(B_i, assign_to='B', standard=2008, source_format='free')

def generate_B_o(B_o, T):

    # Transform the outer boundary condition matrix

    B_o = B_o*T.inv()

    # Convert to Fortran

    return spf.fcode(B_o, assign_to='B', standard=2008, source_format='free')

def generate_C(C, T):

    # Transform the match condition matrix

    C = C*T.inv()

    # Convert to Fortran

    return spf.fcode(C, assign_to='C', standard=2008, source_format='free')

def generate_R(T):

    # Generate the reverse transformation matrix

    R = T.inv()

    # Convert to Fortran

    return spf.fcode(R, assign_to='R', standard=2008, source_format='free')

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

        with open(f'{vars}/B_i_regular.inc', 'w') as f:
            f.write(generate_B_i(B_i_regular, T)+'\n')

        with open(f'{vars}/B_i_zero_r.inc', 'w') as f:
            f.write(generate_B_i(B_i_zero_r, T)+'\n')

        with open(f'{vars}/B_i_zero_h.inc', 'w') as f:
            f.write(generate_B_i(B_i_zero_h, T)+'\n')

        with open(f'{vars}/B_o_vacuum.inc', 'w') as f:
            f.write(generate_B_o(B_o_vacuum, T)+'\n')

        with open(f'{vars}/B_o_zero_r.inc', 'w') as f:
            f.write(generate_B_o(B_o_zero_r, T)+'\n')

        with open(f'{vars}/B_o_zero_h.inc', 'w') as f:
            f.write(generate_B_o(B_o_zero_h, T)+'\n')

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
