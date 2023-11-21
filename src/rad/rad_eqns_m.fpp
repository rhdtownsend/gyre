! Module  : rad_eqns_m
! Purpose : radial adiabatic differential equations
!
! Copyright 2013-2020 Rich Townsend & The GYRE Team
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

module rad_eqns_m

  ! Uses

  use kinds_m

  use context_m
  use eqns_m
  use model_m
  use model_util_m
  use mode_par_m
  use osc_par_m
  use point_m
  use rad_trans_m
  use state_m
  
  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: J_V = 1
  integer, parameter :: J_AS = 2
  integer, parameter :: J_U = 3
  integer, parameter :: J_C_1 = 4
  integer, parameter :: J_GAMMA_1 = 5

  integer, parameter :: J_LAST = J_GAMMA_1

  ! Derived-type definitions

  type, extends (r_eqns_t) :: rad_eqns_t
     private
     type(context_t), pointer   :: cx => null()
     type(point_t), allocatable :: pt(:)
     type(rad_trans_t)          :: tr
     real(WP), allocatable      :: coeff(:,:)
     real(WP)                   :: alpha_omg
   contains
     private
     procedure, public :: stencil
     procedure, public :: A
     procedure, public :: xA
  end type rad_eqns_t

  ! Interfaces

  interface rad_eqns_t
     module procedure rad_eqns_t_
  end interface rad_eqns_t

  ! Access specifiers

  private

  public :: rad_eqns_t

  ! Procedures

contains

  function rad_eqns_t_ (cx, md_p, os_p) result (eq)

    class(context_t), pointer, intent(in) :: cx
    type(mode_par_t), intent(in)          :: md_p
    type(osc_par_t), intent(in)           :: os_p
    type(rad_eqns_t)                      :: eq

    ! Construct the rad_eqns_t

    eq%cx => cx

    eq%tr = rad_trans_t(cx, md_p, os_p)

    select case (os_p%time_factor)
    case ('OSC')
       eq%alpha_omg = 1._WP
    case ('EXP')
       eq%alpha_omg = -1._WP
    case default
       $ABORT(Invalid time_factor)
    end select

    eq%n_e = 2

    ! Finish

    return

  end function rad_eqns_t_

  !****

  subroutine stencil (this, pt)

    class(rad_eqns_t), intent(inout) :: this
    type(point_t), intent(in)        :: pt(:)

    class(model_t), pointer :: ml
    integer :: n_s
    integer :: i

    ! Calculate coefficients at the stencil points

    ml => this%cx%model()

    call check_model(ml, [I_V_2,I_AS,I_U,I_C_1,I_GAMMA_1])

    n_s = SIZE(pt)
      
    if (ALLOCATED(this%coeff)) deallocate(this%coeff)
    allocate(this%coeff(n_s,J_LAST))

    do i = 1, n_s

       $ASSERT(.NOT. ml%is_vacuum(pt(i)),Attempt to stencil at vacuum point)
       
       this%coeff(i,J_V) = ml%coeff(I_V_2, pt(i))*pt(i)%x**2
       this%coeff(i,J_AS) = ml%coeff(I_AS, pt(i))
       this%coeff(i,J_U) = ml%coeff(I_U, pt(i))
       this%coeff(i,J_C_1) = ml%coeff(I_C_1, pt(i))
       this%coeff(i,J_GAMMA_1) = ml%coeff(I_GAMMA_1, pt(i))

    end do

    ! Set up stencil for the tr component

    call this%tr%stencil(pt)

    ! Store the stencil points for on-the-fly evaluations

    this%pt = pt

    ! Finish

    return

  end subroutine stencil

  !****

  function A (this, i, st)

    class(rad_eqns_t), intent(in) :: this
    integer, intent(in)           :: i
    class(r_state_t), intent(in)  :: st
    real(WP)                      :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A = this%xA(i, st)/this%pt(i)%x

    ! Finish

    return

  end function A

  !****

  function xA (this, i, st)

    class(rad_eqns_t), intent(in) :: this
    integer, intent(in)           :: i
    class(r_state_t), intent(in)  :: st
    real(WP)                      :: xA(this%n_e,this%n_e)

    real(WP) :: omega_c
    
    ! Evaluate the log(x)-space RHS matrix

    associate ( &
         omega => st%omega, &
         V => this%coeff(i,J_V), &
         As => this%coeff(i,J_AS), &
         U => this%coeff(i,J_U), &
         c_1 => this%coeff(i,J_C_1), &
         Gamma_1 => this%coeff(i,J_GAMMA_1), &
         alpha_omg => this%alpha_omg)

      omega_c = omega

      ! Set up the matrix

      xA(1,1) = V/Gamma_1 - 1._WP
      xA(1,2) = -V/Gamma_1
      
      xA(2,1) = c_1*alpha_omg*omega_c**2 + U - As
      xA(2,2) = As - U + 3._WP

    end associate

    ! Apply the variables transformation

    call this%tr%trans_eqns(xA, i, st)

    ! Finish

    return

  end function xA

end module rad_eqns_m
