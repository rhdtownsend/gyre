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

  use gyre_eqns
  use gyre_ext
  use gyre_ivp
  use gyre_ivp_factory
  use gyre_num_par
  use gyre_seg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (seg_t) :: ${T}_bvp_seg_t
     private
     class(${T}_ivp_t), allocatable :: iv
     $TYPE(WP), allocatable, public :: y(:,:)
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

  function ${T}_bvp_seg_t_ (x, eq, nm_p) result (bs)

    real(WP), intent(in)           :: x(:)
    class(${T}_eqns_t), intent(in) :: eq
    type(num_par_t), intent(in)    :: nm_p
    type(${T}_bvp_seg_t)           :: bs

    ! Construct the bvp_seg_t

    bs%x = x
    bs%n = SIZE(x)

    allocate(sg%iv, SOURCE=${T}_ivp_t(eq, nm_p))

    bs%n_e = eq%n_e

    allocate(bs%y(bs%n_e,bs%n))

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
