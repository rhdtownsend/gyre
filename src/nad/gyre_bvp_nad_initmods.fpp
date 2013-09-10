! Module   : gyre_bvp_nad_initmods
! Purpose  : modules used by init routine in gyre_bvp_nad (to work around gfortran fixup errors)
!
! Copyright 2013 Rich Townsend
!
! This file is part of GYRE. GYRE is free software: you can
! redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, version 3.
!
! GYRE is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module gyre_bvp_nad_initmods

  use gyre_jacobian_nad_dziem
  use gyre_jacobian_nad_jcd
  use gyre_bound_nad_zero
  use gyre_bound_nad_dziem

end module gyre_bvp_nad_initmods

