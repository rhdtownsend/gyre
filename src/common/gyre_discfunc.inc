! Incfile  : gyre_discfunc
! Purpose  : discriminant root finding (template)
!
! Copyright 2013-2014 Rich Townsend
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

$local $VAR ${T}x

module gyre_${T}_discfunc

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_ext
  use gyre_extfunc

  ! This should not be needed, but it solves unresolved symbol issues
  ! with gfortran
  use gyre_mode

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (${T}_extfunc_t) :: ${T}_discfunc_t
     private
     class(${T}_bvp_t), pointer     :: bp
     $TYPE(WP), allocatable, public :: omega_def(:)
   contains 
     private
     procedure, public :: eval => eval_
  end type ${T}_discfunc_t

  ! Interfaces

  interface ${T}_discfunc_t
     module procedure ${T}_discfunc_t_
  end interface ${T}_discfunc_t

  ! Access specifiers

  private

  public :: ${T}_discfunc_t

  ! Procedures

contains

  function ${T}_discfunc_t_ (bp) result (df)

    class(${T}_bvp_t), pointer, intent(in) :: bp
    type(${T}_discfunc_t)                  :: df

    ! Construct the discfunc_t

    df%bp => bp

    ! Finish

    return

  end function ${T}_discfunc_t_

!****

  function eval_ (this, $VAR) result (f_$VAR)

    class(${T}_discfunc_t), intent(inout) :: this
    type(${T}_ext_t), intent(in)          :: $VAR
    type(${T}_ext_t)                      :: f_$VAR

    $TYPE(WP) :: omega

    ! Evaluate the discriminant

    $if ($T eq 'r')
    omega = real($VAR)
    $elsif ($T eq 'c')
    omega = cmplx($VAR)
    $else
    $stop 'Invalid $T'
    $endif

    f_$VAR = this%bp%discrim(omega)

    if (ALLOCATED(this%omega_def)) then
       f_$VAR = f_$VAR/PRODUCT(omega - this%omega_def)
    endif
    
    ! Finish

    return

  end function eval_

end module gyre_${T}_discfunc

$endsub
