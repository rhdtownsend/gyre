! Module   : gyre_ad_discfunc
! Purpose  : adiabatic discriminant root finding
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

$include 'core.inc'

module gyre_ad_discfunc

  ! Uses

  use core_kinds
  use core_func

  use gyre_ad_bvp
  use gyre_ext_arith

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(func_t) :: ad_discfunc_t
     private
     type(ad_bvp_t), pointer :: bp
     integer                 :: e_norm
   contains 
     private
     procedure, public :: init
     procedure, public :: eval_c
  end type ad_discfunc_t

  ! Access specifiers

  private

  public :: ad_discfunc_t

  ! Procedures

contains

  subroutine init (this, bp, omega_norm)

    class(ad_discfunc_t), intent(out)     :: this
    type(ad_bvp_t), intent(inout), target :: bp
    complex(WP), intent(in)               :: omega_norm

    type(ext_complex_t) :: discrim

    ! Initialize the ad_discfunc

    this%bp => bp

    discrim = this%bp%discrim(omega_norm)

    this%e_norm = discrim%e

    ! Finish

    return

  end subroutine init

!****

  function eval_c (this, z) result (f_z)

    class(ad_discfunc_t), intent(inout) :: this
    complex(WP), intent(in)             :: z
    complex(WP)                         :: f_z

    type(ext_complex_t) :: discrim

    ! Evaluate the normalized discriminant

    discrim = this%bp%discrim(z)

    discrim%e = discrim%e - this%e_norm

    f_z = cmplx(discrim)

    ! Finish

    return

  end function eval_c

end module gyre_ad_discfunc
