#!/usr/bin/env python3
# Program : regen_ad_eqns.fypp
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

omega_c = sp.Symbol('omega_c')

om2 = sp.Symbol('om2')

alpha_omg = sp.Symbol('alpha_omg')

chi = sp.Symbol('chi')
a_11 = sp.Symbol('a_11')
a_12 = sp.Symbol('a_12')
b_11 = sp.Symbol('b_11')
b_12 = sp.Symbol('b_12')

# Declare equation matrices

# Differential Jacobian matrix

A = sp.Matrix([
    [V_g - 1, -V_g],
    [c_1(x)*alpha_omg*omega_c**2 + U(x) - As, As - U(x) + 3]
])

# Match condition matrix

C = sp.Matrix([
    [1, 0],
    [-U(x), U(x)]
])

# Inner boundary condition matrices

B_i_regular = sp.Matrix([
    [c_1(x)*alpha_omg*omega_c**2, 0]
])

B_i_zero_r = sp.Matrix([
    [1, 0]
])

# Outer boundary condition matrices

B_o_vacuum = sp.Matrix([
    [1, -1]
])

B_o_zero_r = sp.Matrix([
    [1, 0]
])

B_o_dziem = sp.Matrix([
    [1 - (4 + c_1(x)*alpha_omg*omega_c**2)/V, -1]
])

B_o_decomp = sp.Matrix([
    [-(chi - a_11), a_12]
])

B_o_jcd = sp.Matrix([
    [chi - b_11, -b_12]
])

# Code generation routines

def generate_A(A, T):

    # Transform the Jacobian matrix

    A = T*(A*T.inv() - x*T.inv().diff(x))

    # Apply substitutions

    A = A.subs(c_1(x)*alpha_omg*omega_c**2, om2)

    # Convert to Fortran, with the invx scaling

    return spf.fcode(invx*A, assign_to='A', standard=2008, source_format='free')

def generate_B_i(B_i, T):

    # Transform the inner boundary condition matrix

    B_i = B_i*T.inv()

    # Apply substitutions

    B_i = B_i.subs(((V, 0), (V_g, 0), (U(x), 3), (As, 0)))

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

    T_gyre = sp.eye(2)

    T_dziem = sp.eye(2)

    T_jcd = sp.Matrix([
        [1, 0],
        [0, 1/(c_1(x)*alpha_omg*omega_c**2)]
    ])

    T_mix = sp.eye(2)

    T_lagp = sp.Matrix([
        [1, 0],
        [-V_2(x), V_2(x)]
    ])

    # Regenerate equation include files

    for vars, T in zip(('gyre', 'dziem', 'jcd', 'mix', 'lagp'), (T_gyre, T_dziem, T_jcd, T_mix, T_lagp)):

        with open(f'{vars}/A.inc', 'w') as f:
            f.write(generate_A(A, T)+'\n')

        with open(f'{vars}/B_i_regular.inc', 'w') as f:
            f.write(generate_B_i(B_i_regular, T)+'\n')

        with open(f'{vars}/B_i_zero_r.inc', 'w') as f:
            f.write(generate_B_i(B_i_zero_r, T)+'\n')

        with open(f'{vars}/B_o_vacuum.inc', 'w') as f:
            f.write(generate_B_o(B_o_vacuum, T)+'\n')

        with open(f'{vars}/B_o_zero_r.inc', 'w') as f:
            f.write(generate_B_o(B_o_zero_r, T)+'\n')

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
