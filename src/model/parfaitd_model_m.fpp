! Module  : parfaitd_model
! Purpose : dimensioned PARFAIT model
!
! Copyright 2023 Rich Townsend & The GYRE Team
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

module parfaitd_model_m

  ! Uses

  use kinds_m

  use constants_m
  use math_m
  use model_m
  use parfait_model_m
  use point_m

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (parfait_model_t) :: parfaitd_model_t
     private
     real(WP), public :: M_star
     real(WP), public :: R_star
   contains
     private
     procedure, public :: M_r
     procedure, public :: P
     procedure, public :: rho
  end type parfaitd_model_t

  ! Interfaces

  interface parfaitd_model_t
     module procedure parfaitd_model_t_
  end interface parfaitd_model_t

  ! Access specifiers

  private

  public :: parfaitd_model_t

  ! Procedures

contains

  function parfaitd_model_t_ (ml_in, M_star, R_star) result (ml)

    type(parfait_model_t), intent(in) :: ml_in
    real(WP), intent(in)              :: M_star
    real(WP), intent(in)              :: R_star
    type(parfaitd_model_t)            :: ml

    ! Construct the parfaitd_model_t

    ml%parfait_model_t = ml_in

    ml%M_star = M_star
    ml%R_star = R_star

    ! Finish

    return

  end function parfaitd_model_t_

  !****

  function M_r (this, pt)

    class(parfaitd_model_t), intent(in) :: this
    type(point_t), intent(in)           :: pt
    real(WP)                            :: M_r

    ! Evaluate the fractional mass coordinate

    M_r = this%M_star*(pt%x**3/this%coeff(I_C_1, pt))

    ! Finish

    return

  end function M_r
    
  !****

  function P (this, pt)

    class(parfaitd_model_t), intent(in) :: this
    type(point_t), intent(in)           :: pt
    real(WP)                            :: P

    ! Evaluate the total pressure

    P = (G_GRAVITY*this%M_star**2/(4._WP*PI*this%R_star**4))* &
        (this%coeff(I_U, pt)/(this%coeff(I_C_1, pt)**2*this%coeff(I_V_2, pt)))

    ! Finish

    return

  end function P
    
  !****

  function rho (this, pt)

    class(parfaitd_model_t), intent(in) :: this
    type(point_t), intent(in)           :: pt
    real(WP)                            :: rho

    ! Evaluate the density

    rho = (this%M_star/(4._WP*PI*this%R_star**3))*(this%coeff(I_U, pt)/this%coeff(I_C_1, pt))

    ! Finish

    return

  end function rho

end module parfaitd_model_m
