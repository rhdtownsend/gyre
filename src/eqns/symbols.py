# Module  : symbols.py
# Purpose : common symbols and routines for equation regeneration
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
import textwrap as tw

# Define a printer class that handles (n,1) materices as vectors

class FCodePrinterExt(spf.FCodePrinter):

    def _print_MatrixElement(self, expr):
        # SymPy matrices normally pass (name, row, col)
        # We flatten it to a single index based on your preference
        name = self._print(expr.parent)

        if expr.parent.shape[1] == 1:
            return f"{name}({expr.i+1})"
        else:
            return super()._print_MatrixElement(expr)

# Declare symbols

x = sp.Symbol('x')
del_x = sp.Symbol('del_x')

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

y_T_1 = sp.Symbol('y_T_1')

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

# Define the code printer

printer = FCodePrinterExt({'standard': 2008, 'source_format': 'free'})

# Code generation routines

def generate_E(A, f, T, transpose=False, subs=[]):

    # Transform the Jacobian matrix and inhomogeneous vector

    A = T*(A*T.inv() - x*T.inv().diff(x))
    f = T*f

    # Apply substitutions

    A = A.subs(fn_subs)
    f = f.subs(fn_subs)

    subs += [
        (sp.Symbol('c_1')*alpha_omg*omega_c**2, om2),
        (lamda/om2, eta)
    ]

    A = A.subs(subs)
    f = f.subs(subs)

    # Convert to Fortran

    if transpose:
        code = printer.doprint(del_x*A.T, assign_to='A_t')
    else:
        code = printer.doprint(del_x*A, assign_to='A')

    code += f"""
if (PRESENT(f)) then
{printer.doprint(del_x*f, assign_to='f')}
endif"""

    return code

def generate_IB(IB, g, T, transpose=False, subs=[]):

    # Transform the inner boundary condition matrix

    IB = IB*T.inv()

    # Apply substitutions

    IB = IB.subs(fn_subs)
    g = g.subs(fn_subs)

    subs += [
        (V, 0),
        (As, 0),
        (As_iso, 0),
        (sp.Symbol('U'), 3)
    ]

    IB = IB.subs(subs)
    g = g.subs(subs)

    # Convert to Fortran

    if transpose:
        code = printer.doprint(IB.T, assign_to='B_t')
    else:
        code = printer.doprint(IB, assign_to='B')

    code += f"""
if (PRESENT(f)) then
{printer.doprint(g, assign_to='f')}
endif"""

    return code

def generate_OB(OB, g, T, transpose=False, subs=[]):

    # Transform the outer boundary condition matrix

    OB = OB*T.inv()

    # Apply substitutions

    OB = OB.subs(fn_subs)
    g = g.subs(fn_subs)

    OB = OB.subs(subs)
    g = g.subs(subs)

    # Convert to Fortran

    if transpose:
        code = printer.doprint(OB.T, assign_to='B_t')
    else:
        code = printer.doprint(OB, assign_to='B')

    code += f"""
if (PRESENT(f)) then
{printer.doprint(g, assign_to='f')}
endif"""

    return code

def generate_C(C, T, transpose=False, subs=[]):

    # Transform the match condition matrix

    C = C*T.inv()

    # Apply substitutions

    C = C.subs(fn_subs)

    C = C.subs(subs)

    # Convert to Fortran

    if transpose:
        return printer.doprint(C.T, assign_to='C_t')
    else:
        return printer.doprint(C, assign_to='C')

def generate_R(T, transpose=False, subs=[]):

    # Invert the transformation matrix

    R = T.inv()

    # Apply substitutions

    R = R.subs(fn_subs)

    R = R.subs(subs)

    # Convert to Fortran

    if transpose:
        return printer.doprint(R.T, assign_to='R_t')
    else:
        return printer.doprint(R, assign_to='R')

def generate(expr, assign_to=None, subs=[]):

    # Apply substitutions

    expr = expr.subs(fn_subs)

    expr = expr.subs(subs)

    # Convert a generic expression to Fortran

    return printer.doprint(expr, assign_to=assign_to)
