! Incfile  : gyre_discrim_func
! Purpose  : discriminant function (template)
!
! Copyright 2013-2016 Rich Townsend
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

module gyre_${T}_discrim_func

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_ext
  use gyre_ext_func
  use gyre_status

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (${T}_ext_func_t) :: ${T}_discrim_func_t
     private
     class(${T}_bvp_t), pointer     :: bp
     real(WP)                       :: omega_min
     real(WP)                       :: omega_max
     $TYPE(WP), allocatable, public :: omega_def(:)
   contains 
     private
     procedure, public :: eval => eval_
  end type ${T}_discrim_func_t

  ! Interfaces

  interface ${T}_discrim_func_t
     module procedure ${T}_discrim_func_t_
  end interface ${T}_discrim_func_t

  ! Access specifiers

  private

  public :: ${T}_discrim_func_t

  ! Procedures

contains

  function ${T}_discrim_func_t_ (bp, omega_min, omega_max) result (df)

    class(${T}_bvp_t), pointer, intent(in) :: bp
    real(WP), intent(in)                   :: omega_min
    real(WP), intent(in)                   :: omega_max
    type(${T}_discrim_func_t)              :: df

    ! Construct the discrim_func_t

    df%bp => bp

    df%omega_min = omega_min
    df%omega_max = omega_max

    ! Finish

    return

  end function ${T}_discrim_func_t_

!****

  subroutine eval_ (this, $VAR, f_$VAR, status)

    class(${T}_discrim_func_t), intent(inout) :: this
    type(${T}_ext_t), intent(in)              :: $VAR
    type(${T}_ext_t), intent(out)             :: f_$VAR
    integer, intent(out)                      :: status

    $TYPE(WP) :: omega

    ! Evaluate the discriminant

    $if ($T eq 'r')
    omega = real($VAR)
    $elsif ($T eq 'c')
    omega = cmplx($VAR)
    $else
    $stop 'Invalid $T'
    $endif

    if (REAL(omega, WP) >= this%omega_min .AND. REAL(omega, WP) <= this%omega_max) then

       call this%bp%build(omega)

       f_$VAR = this%bp%det()

       status = STATUS_OK

       if (ALLOCATED(this%omega_def)) then
          f_$VAR = f_$VAR*PRODUCT(omega/(omega - this%omega_def))
       endif

    else

       status = STATUS_OMEGA_DOMAIN

    endif
    
    ! Finish

    return

  end subroutine eval_

end module gyre_${T}_discrim_func

$endsub
