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

# Declare symbols

x = sp.Symbol('x')
invx = sp.Symbol('invx')

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

f_conv = sp.Symbol('f_conv')

f_rht = sp.Symbol('f_rht')
df_rht = sp.Symbol('df_rht')

lamda = sp.Symbol('lambda')
l_i = sp.Symbol('l_i')
l_e = sp.Symbol('l_e')
omega_c = sp.Symbol('omega_c')
i_omega_c = sp.Symbol('i_omega_c')

F_trb = sp.Symbol('F_trb')

eta = sp.Symbol('eta')
om2 = sp.Symbol('om2')

l = sp.Symbol('l')
m = sp.Symbol('m')

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

# Define functions w/ known derivatives

class V_2(sp.Function):
    def fdiff(self, argindex=1):
        return V_2(x)*dV_2/self.args[0]

class c_1(sp.Function):
    def fdiff(self, argindex=1):
        return c_1(self.args[0])*(3 - U(self.args[0]))/self.args[0]

class U(sp.Function):
    def fdiff(self, argindex=1):
        return U(self.args[0])*(-V_g - As - U(self.args[0]) + 3)/self.args[0]

# Define function substitutions (to replace functions with symbols
# after differentiation)

fn_subs = [
        (V_2(x), sp.Symbol('V_2')),
        (U(x), sp.Symbol('U')),
        (c_1(x), sp.Symbol('c_1'))
]

# Declare fcode_symbol, for overriding Fortran output

class fcode_symbol(sp.Symbol):
    def __new__(cls, name, *args):
        return super().__new__(cls, name)
    def __init__(self, name, code=None):
        self.code = code if code is not None else name
    def _fcode(self, printer):
        return self.code

# Code generation routines

def generate_A(A, T, subs=[]):

    # Transform the Jacobian matrix

    A = T*(A*T.inv() - x*T.inv().diff(x))

    # Apply substitutions

    A = A.subs(fn_subs)

    subs += [
        (sp.Symbol('c_1')*alpha_omg*omega_c**2, om2),
        (lamda/om2, eta)
    ]

    A = A.subs(subs)

    # Convert to Fortran

    printer = spf.FCodePrinter({'standard': 2008, 'source_format': 'free'})

    return printer.doprint(invx*A, assign_to='A')

def generate_IB(IB, T, subs=[]):

    # Transform the inner boundary condition matrix

    IB = IB*T.inv()

    # Apply substitutions

    IB = IB.subs(fn_subs)

    subs += [
        (V, 0),
        (As, 0),
        (As_iso, 0),
        (sp.Symbol('U'), 3)
    ]

    IB = IB.subs(subs)

    # Convert to Fortran

    return spf.fcode(IB, assign_to='B', standard=2008, source_format='free')

def generate_OB(OB, T, subs=[]):

    # Transform the outer boundary condition matrix

    OB = OB*T.inv()

    # Apply substitutions

    OB = OB.subs(fn_subs)

    OB = OB.subs(subs)

    # Convert to Fortran

    return spf.fcode(OB, assign_to='B', standard=2008, source_format='free')

def generate_C(C, T, subs=[]):

    # Transform the match condition matrix

    C = C*T.inv()

    # Apply substitutions

    C = C.subs(fn_subs)

    C = C.subs(subs)

    # Convert to Fortran

    return spf.fcode(C, assign_to='C', standard=2008, source_format='free')

def generate_R(T, subs=[]):

    # Invert the transformation matrix

    R = T.inv()

    # Apply substitutions

    R = R.subs(fn_subs)

    R = R.subs(subs)

    # Convert to Fortran

    return spf.fcode(R, assign_to='R', standard=2008, source_format='free')

def generate(expr, assign_to=None, subs=[]):

    # Apply substitutions

    expr = expr.subs(fn_subs)

    expr = expr.subs(subs)

    # Convert a generic expression to Fortran

    return spf.fcode(expr, assign_to=assign_to, standard=2008, source_format='free')
