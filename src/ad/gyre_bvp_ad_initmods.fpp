! Module   : gyre_bvp_ad_initmods
! Purpose  : modules used by init routine in gyre_bvp_ad (to work around gfortran fixup errors)
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

module gyre_bvp_ad_initmods

  use gyre_jacobian_ad_dziem
  use gyre_jacobian_ad_jcd
  use gyre_bound_ad_zero
  use gyre_bound_ad_dziem
  use gyre_bound_ad_unno
  use gyre_bound_ad_jcd

end module gyre_bvp_ad_initmods
