! Module   : gyre_discfunc
! Purpose  : discriminant root finding
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

module gyre_discfunc

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_ext_arith
  use gyre_ext_func

  ! This should not be needed, but it solves unresolved symbol issues
  ! with gfortran
  use gyre_mode

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(ext_func_t) :: discfunc_t
     private
     class(bvp_t), pointer            :: bp
     complex(WP), allocatable, public :: omega_def(:)
   contains 
     private
     procedure :: eval_er_
     procedure :: eval_ec_
  end type discfunc_t

  ! Interfaces

  interface discfunc_t
     module procedure discfunc_t_
  end interface discfunc_t

  ! Access specifiers

  private

  public :: discfunc_t

  ! Procedures

contains

  function discfunc_t_ (bp) result (df)

    class(bvp_t), pointer, intent(in) :: bp
    type(discfunc_t)                  :: df

    ! Construct the discfunc_t

    df%bp => bp

    ! Finish

    return

  end function discfunc_t_

!****

  function eval_er_ (this, ex) result (f_ex)

    class(discfunc_t), intent(inout) :: this
    type(ext_real_t), intent(in)     :: ex
    type(ext_real_t)                 :: f_ex

    ! Evaluate the discriminant for real frequencies

    associate (omega => REAL(ex))

      f_ex = this%bp%discrim(omega)

      if (ALLOCATED(this%omega_def)) then
         f_ex = f_ex/PRODUCT(omega - REAL(this%omega_def))
      endif

    end associate

    ! Finish

    return

  end function eval_er_

!****

  function eval_ec_ (this, ez) result (f_ez)

    class(discfunc_t), intent(inout) :: this
    type(ext_complex_t), intent(in)  :: ez
    type(ext_complex_t)              :: f_ez

    ! Evaluate the discriminant for complex frequencies

    associate (omega => CMPLX(ez))

      f_ez = this%bp%discrim(omega)

      if(ALLOCATED(this%omega_def)) then
         f_ez = f_ez/PRODUCT(omega - this%omega_def)
      endif

    end associate

    ! Finish

    return

  end function eval_ec_

end module gyre_discfunc
