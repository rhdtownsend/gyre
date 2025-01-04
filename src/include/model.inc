#! Include : model
#! Purpose : model-related macros
#!
#! Copyright 2025 Rich Townsend
#!
#! This file is part of GYRE. GYRE is free software: you can
#! redistribute it and/or modify it under the terms of the GNU General
#! Public License as published by the Free Software Foundation, version 3.
#!
#! GYRE is distributed in the hope that it will be useful, but WITHOUT
#! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
#! License for more details.
#!
#! You should have received a copy of the GNU General Public License
#! along with this program.  If not, see <http://www.gnu.org/licenses/>.

#!!! Canonical coefficient names

#:set COEFF_NAMES = ('V_2', 'As', 'U', 'c_1', 'Gamma_1', 'ups_T', 'nabla_ad', &
   'nabla', 'beta_rad', 'c_lum', 'c_rad', 'c_thn', 'c_thk', 'c_eps', 'c_egv', &
    'eps_rho', 'eps_T', 'kap_rho', 'kap_T', 'Omega_rot')

#!!! Canonical coefficient index names

#:set COEFF_IDX_NAMES = ['I_{:s}'.format(name.upper()) for name in COEFF_NAMES]
