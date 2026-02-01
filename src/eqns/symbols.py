# Module  : symbols.py
# Purpose : common symbols and routines for equation regeneration
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

# Declare symbols & functions

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

f_conv = sp.Symbol('f_conv')

lamda = sp.Symbol('lambda')
l_i = sp.Symbol('l_i')
l_e = sp.Symbol('l_e')
omega_c = sp.Symbol('omega_c')
i_omega_c = sp.Symbol('i_omega_c')

F_trb = sp.Symbol('F_trb')

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

# Code generation routines

def generate_A(A, T):

    # Transform the Jacobian matrix

    A = T*(A*T.inv() - x*T.inv().diff(x))

    # Apply substitutions

    A = A.subs(((c_1(x)*alpha_omg*omega_c**2, om2), (lamda/om2, eta)))

    # Convert to Fortran

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

    # Invert the transformation matrix

    R = T.inv()

    # Convert to Fortran

    return spf.fcode(R, assign_to='R', standard=2008, source_format='free')

def generate(expr, assign_to=None):

    # Convert a generic expression to Fortran

    return spf.fcode(expr, assign_to=assign_to, standard=2008, source_format='free')
