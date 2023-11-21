! Module  : sysmtx_factory_m
! Purpose : factory procedures for r_sysmtx_t and c_sysmtx_t types
!
! Copyright 2013-2015 Rich Townsend
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

module sysmtx_factory_m

  ! Uses

  use kinds_m

  use sysmtx_m
  use band_sysmtx_m
  use block_sysmtx_m
  use num_par_m

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface r_sysmtx_t
     module procedure r_sysmtx_t_
  end interface r_sysmtx_t

  interface c_sysmtx_t
     module procedure c_sysmtx_t_
  end interface c_sysmtx_t

  ! Access specifiers

  private

  public :: r_sysmtx_t
  public :: c_sysmtx_t

  ! Procedures

contains

  $define $SYSMTX_T $sub

  $local $T $1

  function ${T}_sysmtx_t_ (n, n_e, n_i, n_o, np) result (sm)

    integer, intent(in)               :: n
    integer, intent(in)               :: n_e
    integer, intent(in)               :: n_i
    integer, intent(in)               :: n_o
    type(num_par_t), intent(in)       :: np
    class(${T}_sysmtx_t), allocatable :: sm
    
    ! Create a ${T}_sysmtx_t

    select case (np%matrix_type)
    case ('BAND')
       allocate(sm, SOURCE=${T}_band_sysmtx_t(n, n_e, n_i, n_o))
    case ('BLOCK')
       allocate(sm, SOURCE=${T}_block_sysmtx_t(n, n_e, n_i, n_o))
    case default
       $ABORT(Invalid matrix_type)
    end select

    ! Finish

    return

  end function ${T}_sysmtx_t_

  $endsub

  $SYSMTX_T(r)
  $SYSMTX_T(c)

end module sysmtx_factory_m
