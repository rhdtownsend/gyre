! Incfile  : gyre_bvp_seg
! Purpose  : boundary-value solver segment (template)
!
! Copyright 2015 Rich Townsend
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

module gyre_${T}_bvp_seg

  ! Uses

  use core_kinds

  use gyre_ext
  use gyre_ivp

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ${T}_bvp_seg_t
     private
     class(${T}_ivp_t), allocatable :: iv
     real(WP), allocatable          :: x(:)
     $TYPE(WP), allocatable         :: y(:,:)
     integer, public                :: n
     integer, public                :: n_e
   contains 
     private
     procedure, public :: build => build_
  end type ${T}_bvp_seg_t

  ! Interfaces

  interface ${T}_bvp_seg_t
     module procedure ${T}_bvp_seg_t_
  end interface ${T}_bvp_seg_t

  ! Access specifiers

  private

  public :: ${T}_bvp_seg_t

contains

  function ${T}_bvp_seg_t_ (x, iv) result (sg)

    real(WP), intent(in)          :: x(:)
    class(${T}_ivp_t), intent(in) :: iv
    type(${T}_bvp_seg_t)          :: sg

    ! Construct the bvp_seg_t

    sg%n = SIZE(x)
    sg%n_e = iv%n_e

    allocate(sg%iv, SOURCE=iv)

    sg%x = x
    allocate(sg%y(sg%n_e,sg%n))

    ! Finish
    
    return

  end function ${T}_bvp_seg_t_

!****

  subroutine build_ (this, omega, k, E_l, E_r, S)

    class(${T}_bvp_t), target, intent(inout) :: this
    $TYPE(WP), intent(in)                    :: omega
    integer, intent(in)                      :: k
    $TYPE(WP), intent(out)                   :: E_l(:,:)
    $TYPE(WP), intent(out)                   :: E_r(:,:)
    type(${T}_ext_t), intent(out)            :: S

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)
    
    ! Build the k'th block within the segment

    call this%iv%shoot(omega, this%x(k), this%x(k+1), E_l, E_r, S)

    ! Finish

    return

  end subroutine build_

end module gyre_${T}_bvp_seg

$endsub
