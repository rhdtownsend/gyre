! Program  : gyre_force_util
! Purpose  : forced oscillation-related utility functions
!
! Copyright 2018-2021 Rich Townsend & The GYRE Team
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

$include 'core.inc'

module gyre_force_util

  ! Uses

  use core_kinds

  use gyre_force_par
  use gyre_mode_par
  use gyre_orbit_par
  use gyre_tide_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: Phi_force

  ! Procedures

contains

  function Phi_force (md_p, fr_p, or_p)

    type(mode_par_t), intent(in)  :: md_p
    type(force_par_t), intent(in) :: fr_p
    type(orbit_par_t), intent(in) :: or_p
    real(WP)                      :: Phi_force

    real(WP) :: R_a
    real(WP) :: eps_tide

    ! Evaluate the surface forcing potential (at x=1), in units of
    ! G*M/R

    select case (fr_p%force_type)

    case ('FIXED')

       Phi_force = -fr_p%Phi

    case ('BINARY')

       R_a = tidal_R_a(or_p%Omega_orb, or_p%q)

       eps_tide = R_a**3*or_p%q

       Phi_force = -eps_tide*tidal_c(R_a, or_p%e, md_p%l, md_p%m, fr_p%k)

    case default

       $ABORT(Invalid force_type)

    end select

    ! Finish

    return

  end function Phi_force

end module gyre_force_util
