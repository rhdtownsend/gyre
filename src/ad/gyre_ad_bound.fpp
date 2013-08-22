! Module   : gyre_ad_bound
! Purpose  : adiabatic boundary conditions
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

module gyre_ad_bound

  ! Uses

  use core_kinds

  use gyre_coeffs
  use gyre_oscpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ad_bound_t
     private
     class(coeffs_t), pointer :: cf => null()
     type(oscpar_t), pointer  :: op => null()
     integer, public          :: n_e
     integer, public          :: n_i
     integer, public          :: n_o
   contains 
     private
     procedure, public :: init
     procedure, public :: inner_bound
     procedure, public :: outer_bound
     procedure, public :: outer_bound_zero
     procedure, public :: outer_bound_dziem
     procedure, public :: outer_bound_unno
     procedure, public :: outer_bound_jcd
  end type ad_bound_t

  ! Access specifiers

  private

  public :: ad_bound_t
  public :: eval_outer_coeffs_unno
  public :: eval_outer_coeffs_jcd
  public :: outer_wavenumber
  public :: eval_cutoffs

  ! Procedures

contains

  subroutine init (this, cf, op)

    class(ad_bound_t), intent(out)      :: this
    class(coeffs_t), intent(in), target :: cf
    type(oscpar_t), intent(in), target  :: op

    ! Initialize the ad_bound

    this%cf => cf
    this%op => op

    this%n_i = 2
    this%n_o = 2
    this%n_e = this%n_i + this%n_o

    ! Finish

    return
    
  end subroutine init

!****

  function inner_bound (this, x_i, omega) result (B_i)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: x_i
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: B_i(this%n_i,this%n_e)

    $ASSERT(x_i == 0._WP,Boundary condition invalid for x_i /= 0)

    ! Set the inner boundary conditions to enforce non-diverging modes

    associate(c_1 => this%cf%c_1(x_i), l => this%op%l, &
              omega_c => this%cf%omega_c(x_i, this%op%m, omega))
                 
      B_i(1,1) = c_1*omega_c**2
      B_i(1,2) = -l
      B_i(1,3) = 0._WP
      B_i(1,4) = 0._WP
        
      B_i(2,1) = 0._WP
      B_i(2,2) = 0._WP
      B_i(2,3) = l
      B_i(2,4) = -1._WP
      
    end associate

    ! Finish

    return

  end function inner_bound

!****

  function outer_bound (this, x_o, omega) result (B_o)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: x_o
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: B_o(this%n_o,this%n_e)

    ! Set the outer boundary conditions

    select case (this%op%outer_bound_type)
    case ('ZERO')
       B_o = this%outer_bound_zero(x_o, omega)
    case ('DZIEM')
       B_o = this%outer_bound_dziem(x_o, omega)
    case ('UNNO')
       B_o = this%outer_bound_unno(x_o, omega)
    case ('JCD')
       B_o = this%outer_bound_jcd(x_o, omega)
    case default
       $ABORT(Invalid outer_bound_type)
    end select

    ! Finish

    return

  end function outer_bound

!****

  function outer_bound_zero (this, x_o, omega) result (B_o)

    class(ad_bound_t), intent(in) :: this
    real(WP)                      :: x_o
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: B_o(this%n_o,this%n_e)

    ! Set the outer boundary conditions, assuming delta p -> 0

    associate(U => this%cf%U(x_o), l => this%op%l)

      B_o(1,1) = 1._WP
      B_o(1,2) = -1._WP
      B_o(1,3) = 1._WP
      B_o(1,4) = 0._WP
      
      B_o(2,1) = U
      B_o(2,2) = 0._WP
      B_o(2,3) = l + 1._WP
      B_o(2,4) = 1._WP

    end associate

    ! Finish

    return

  end function outer_bound_zero

!****

  function outer_bound_dziem (this, x_o, omega) result (B_o)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: x_o
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: B_o(this%n_o,this%n_e)

    ! Set the outer boundary conditions, assuming Dziembowski's (1971)
    ! condition: d(delta p)/dr -> 0 for an isothermal atmosphere

    associate(V => this%cf%V(x_o), c_1 => this%cf%V(x_o), &
              l => this%op%l, omega_c => this%cf%omega_c(x_o, this%op%m, omega))
        
      B_o(1,1) = 1 + (l*(l+1)/(c_1*omega_c**2) - 4._WP - c_1*omega_c**2)/V
      B_o(1,2) = -1._WP
      B_o(1,3) = 1 + (l*(l+1)/(c_1*omega_c**2) - l - 1._WP)/V
      B_o(1,4) = 0._WP
      
      B_o(2,1) = 0._WP
      B_o(2,2) = 0._WP
      B_o(2,3) = l + 1._WP
      B_o(2,4) = 1._WP

    end associate

    ! Finish

    return

  end function outer_bound_dziem

!****

  function outer_bound_unno (this, x_o, omega) result (B_o)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: x_o
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: B_o(this%n_o,this%n_e)

    real(WP)    :: V_g
    real(WP)    :: As
    real(WP)    :: c_1
    complex(WP) :: lambda
    complex(WP) :: b_11
    complex(WP) :: b_12
    complex(WP) :: b_13
    complex(WP) :: b_21
    complex(WP) :: b_22
    complex(WP) :: b_23
    complex(WP) :: alpha_1
    complex(WP) :: alpha_2

    ! Set the outer boundary conditions, assuming Unno et al.'s (1989,
    ! S18.1) formulation.

    call eval_outer_coeffs_unno(this%cf, x_o, V_g, As, c_1)

    associate(l => this%op%l, omega_c => this%cf%omega_c(x_o, this%op%m, omega))

      lambda = outer_wavenumber(V_g, As, c_1, omega_c, l)
      
      b_11 = V_g - 3._WP
      b_12 = l*(l+1)/(c_1*omega_c**2) - V_g
      b_13 = V_g

      b_21 = c_1*omega_c**2 - As
      b_22 = 1._WP + As
      b_23 = -As
    
      alpha_1 = (b_12*b_23 - b_13*(b_22+l))/((b_11+l)*(b_22+l) - b_12*b_21)
      alpha_2 = (b_21*b_13 - b_23*(b_11+l))/((b_11+l)*(b_22+l) - b_12*b_21)

      B_o(1,1) = (lambda - b_11)/b_12
      B_o(1,2) = -1._WP
      B_o(1,3) = -(alpha_1*(lambda - b_11)/b_12 - alpha_2)
      B_o(1,4) = 0._WP

      B_o(2,1) = 0._WP
      B_o(2,2) = 0._WP
      B_o(2,3) = l + 1._WP
      B_o(2,4) = 1._WP

    end associate

    ! Finish

    return

  end function outer_bound_unno

!****

  function outer_bound_jcd (this, x_o, omega) result (B_o)

    class(ad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: x_o
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: B_o(this%n_o,this%n_e)

    real(WP)    :: V_g
    real(WP)    :: As
    real(WP)    :: c_1
    complex(WP) :: lambda
    complex(WP) :: b_11
    complex(WP) :: b_12

    ! Set the outer boundary conditions, assuming
    ! Christensen-Dalsgaard's formulation (see ADIPLS documentation)

    call eval_outer_coeffs_jcd(this%cf, x_o, V_g, As, c_1)

    associate(l => this%op%l, omega_c => this%cf%omega_c(x_o, this%op%m, omega))

      lambda = outer_wavenumber(V_g, As, c_1, omega_c, l)

      b_11 = V_g - 3._WP
      b_12 = l*(l+1)/(c_1*omega_c**2) - V_g

      if(l /= 0) then
         B_o(1,1) = (lambda - b_11)/b_12
         B_o(1,2) = -1._WP
         B_o(1,3) = 1._WP + (l*(l+1)/(c_1*omega_c**2) - l - 1._WP)/(V_g + As)
         B_o(1,4) = 0._WP
      else
         B_o(1,1) = (lambda - b_11)/b_12
         B_o(1,2) = -1._WP
         B_o(1,3) = 1._WP
         B_o(1,4) = 0._WP
      endif

      B_o(2,1) = 0._WP
      B_o(2,2) = 0._WP
      B_o(2,3) = l + 1._WP
      B_o(2,4) = 1._WP

    end associate

    ! Finish

    return

  end function outer_bound_jcd

!****

  subroutine eval_outer_coeffs_unno (cf, x_o, V_g, As, c_1)

    class(coeffs_t), intent(in) :: cf
    real(WP), intent(in)        :: x_o
    real(WP), intent(out)       :: V_g
    real(WP), intent(out)       :: As
    real(WP), intent(out)       :: c_1

    ! Calculate coefficients at the outer boundary, for use in the
    ! Unno boundary prescription

    V_g = cf%V(x_o)/cf%Gamma_1(x_o)
    As = cf%As(x_o)
    c_1 = cf%c_1(x_o)

    ! Finish

    return

  end subroutine eval_outer_coeffs_unno
    
!****

  subroutine eval_outer_coeffs_jcd (cf, x_o, V_g, As, c_1)

    class(coeffs_t), intent(in) :: cf
    real(WP), intent(in)        :: x_o
    real(WP), intent(out)       :: V_g
    real(WP), intent(out)       :: As
    real(WP), intent(out)       :: c_1

    ! Calculate coefficients at the outer boundary, for use in the
    ! JCD boundary prescription

    V_g = cf%V(x_o)/cf%Gamma_1(x_o)
    As = cf%V(x_o)*(1._WP-1._WP/cf%Gamma_1(x_o))
    c_1 = cf%c_1(x_o)

    ! Finish

    return

  end subroutine eval_outer_coeffs_jcd
  
!****

  function outer_wavenumber (V_g, As, c_1, omega_c, l) result (lambda)

    real(WP)                :: V_g
    real(WP), intent(in)    :: As
    real(WP), intent(in)    :: c_1
    complex(WP), intent(in) :: omega_c
    integer, intent(in)     :: l
    complex(WP)             :: lambda

    real(WP)    :: omega_c_cutoff_lo
    real(WP)    :: omega_c_cutoff_hi
    complex(WP) :: gamma
    complex(WP) :: sgamma

    ! Calculate the wavenumber at the outer boundary

    if(AIMAG(omega_c) == 0._WP) then

       ! Calculate cutoff frequencies

       call eval_cutoffs_from_coeffs(V_g, As, c_1, l, omega_c_cutoff_lo, omega_c_cutoff_hi)

       ! Evaluate the wavenumber

       gamma = -4._WP*V_g*c_1*(omega_c**2 - omega_c_cutoff_lo**2)*(omega_c**2 - omega_c_cutoff_hi**2)/omega_c**2

       if(ABS(REAL(omega_c)) > omega_c_cutoff_hi) then

          ! Acoustic waves

          lambda = 0.5_WP*((V_g + As - 2._WP) - SQRT(gamma))

       elseif(ABS(REAL(omega_c)) < omega_c_cutoff_lo) then

          ! Gravity waves

          lambda = 0.5_WP*((V_g + As - 2._WP) + SQRT(gamma))

       else

          ! Evanescent

          lambda = 0.5_WP*((V_g + As - 2._WP) - SQRT(gamma))

       endif

    else

       ! Evaluate the wavenumber

       gamma = (As - V_g + 4._WP)**2 + 4*(l*(l+1)/(c_1*omega_c**2) - V_g)*(c_1*omega_c**2 - As)
       sgamma = SQRT(gamma)

       if(AIMAG(omega_c) > 0._WP) then

          ! Decaying oscillations; choose the wave with diverging
          ! energy density (see Townsend 2000b)

          if(REAL(sgamma) > 0._WP) then
             lambda = 0.5_WP*((V_g + As - 2._WP) + sgamma)
          else
             lambda = 0.5_WP*((V_g + As - 2._WP) - sgamma)
          endif

       else

          ! Growing oscillations; choose the wave with non-diverging
          ! energy density (see Townsend 2000b)

          if(REAL(sgamma) > 0._WP) then
             lambda = 0.5_WP*((V_g + As - 2._WP) - sgamma)
          else
             lambda = 0.5_WP*((V_g + As - 2._WP) + sgamma)
          endif

       endif

    end if

    ! Finish

    return

  end function outer_wavenumber

!****

  subroutine eval_cutoffs (cf, op, x_o, omega_cutoff_lo, omega_cutoff_hi)

    class(coeffs_t), intent(in) :: cf
    type(oscpar_t), intent(in)  :: op
    real(WP), intent(in)        :: x_o
    real(WP), intent(out)       :: omega_cutoff_lo
    real(WP), intent(out)       :: omega_cutoff_hi

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: c_1
    real(WP) :: omega_c_cutoff_lo
    real(WP) :: omega_c_cutoff_hi

    ! Calculate coefficients at the outer boundary

    select case (op%outer_bound_type)
    case ('ZERO')
       $ABORT(Cutoff frequencies are undefined for ZERO outer_bound_type)
    case ('DZIEM')
       $ABORT(Cutoff frequencies are undefined for DZIEM outer_bound_type)
    case ('UNNO')
       call eval_outer_coeffs_unno(cf, x_o, V_g, As, c_1)
    case('JCD')
       call eval_outer_coeffs_jcd(cf, x_o, V_g, As, c_1)
    case default
       $ABORT(Invalid outer_bound_type)
    end select

    ! Evaluate the cutoff freqs

    call eval_cutoffs_from_coeffs(V_g, As, c_1, op%l, omega_c_cutoff_lo, omega_c_cutoff_hi)

    omega_cutoff_lo = REAL(cf%omega(x_o, op%m, CMPLX(omega_c_cutoff_lo, KIND=WP)))
    omega_cutoff_hi = REAL(cf%omega(x_o, op%m, CMPLX(omega_c_cutoff_hi, KIND=WP)))

    ! Finish

    return

  end subroutine eval_cutoffs

!****

  subroutine eval_cutoffs_from_coeffs (V_g, As, c_1, l, omega_c_cutoff_lo, omega_c_cutoff_hi)

    real(WP), intent(in)  :: V_g
    real(WP), intent(in)  :: As
    real(WP), intent(in)  :: c_1
    integer, intent(in)   :: l
    real(WP), intent(out) :: omega_c_cutoff_lo
    real(WP), intent(out) :: omega_c_cutoff_hi

    real(WP) :: a
    real(WP) :: b
    real(WP) :: c

    ! Evaluate the cutoff frequencies from the coefficients at the
    ! outer boundary

    a = -4._WP*V_g*c_1**2
    b = ((As - V_g + 4._WP)**2 + 4._WP*V_g*As + 4._WP*l*(l+1))*c_1
    c = -4._WP*l*(l+1)*As

    omega_c_cutoff_lo = SQRT((-b + SQRT(b**2 - 4._WP*a*c))/(2._WP*a))
    omega_c_cutoff_hi = SQRT((-b - SQRT(b**2 - 4._WP*a*c))/(2._WP*a))
    
    $ASSERT(omega_c_cutoff_hi >= omega_c_cutoff_lo,Incorrect cutoff frequency ordering)

    ! Finish

    return

  end subroutine eval_cutoffs_from_coeffs

end module gyre_ad_bound
