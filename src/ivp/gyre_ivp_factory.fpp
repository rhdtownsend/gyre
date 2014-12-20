! Module   : gyre_ivp_factory
! Purpose  : factory procedures for r_ivp_t and c_ivp_t types
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

module gyre_ivp_factory

  ! Uses

  use core_kinds

  use gyre_colloc_ivp
  use gyre_findiff_ivp
  use gyre_ivp
  use gyre_jacob
  use gyre_numpar
  use gyre_magnus_ivp

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface r_ivp_t
     module procedure r_ivp_t_
  end interface r_ivp_t

  interface c_ivp_t
     module procedure c_ivp_t_
  end interface c_ivp_t

  ! Access specifiers

  private

  public :: r_ivp_t
  public :: c_ivp_t

  ! Procedures

contains

  $define $IVP_T $sub

  $local $T $1

  function ${T}_ivp_t_ (jc, np) result (iv)

    class(${T}_jacob_t), intent(in) :: jc
    type(numpar_t), intent(in)      :: np
    class(${T}_ivp_t), allocatable  :: iv
    
    ! Create a ${T}_ivp_t

    select case (np%ivp_solver)
    case ('MAGNUS_GL2')
       allocate(iv, SOURCE=${T}_magnus_ivp_t(jc, 'GL2'))
    case ('MAGNUS_GL4')
       allocate(iv, SOURCE=${T}_magnus_ivp_t(jc, 'GL4'))
    case ('MAGNUS_GL6')
       allocate(iv, SOURCE=${T}_magnus_ivp_t(jc, 'GL6'))
    case ('COLLOC_GL2')
       allocate(iv, SOURCE=${T}_colloc_ivp_t(jc, 'GL2'))
    case ('COLLOC_GL4')
       allocate(iv, SOURCE=${T}_colloc_ivp_t(jc, 'GL4'))
    case ('FINDIFF')
       allocate(iv, SOURCE=${T}_findiff_ivp_t(jc))
    case default
       $ABORT(Invalid ivp_solver)
    end select

    ! Finish

    return

  end function ${T}_ivp_t_

  $endsub

  $IVP_T(r)
  $IVP_T(c)

end module gyre_ivp_factory
