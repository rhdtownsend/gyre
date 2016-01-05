! Incfile  : gyre_nad_bound
! Purpose  : nonadiabatic boundary conditions
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

module gyre_nad_bound

  ! Uses

  use core_kinds

  use gyre_atmos
  use gyre_bound
  use gyre_ext
  use gyre_mode_par
  use gyre_model
  use gyre_nad_vars
  use gyre_osc_par
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: REGULAR_TYPE = 1
  integer, parameter :: ZERO_TYPE = 2
  integer, parameter :: DZIEM_TYPE = 3
  integer, parameter :: UNNO_TYPE = 4
  integer, parameter :: JCD_TYPE = 5

  ! Derived-type definitions

  type, extends (c_bound_t) :: nad_bound_t
     private
     class(model_t), pointer     :: ml => null()
     class(c_rot_t), allocatable :: rt
     type(nad_vars_t)            :: vr
     integer                     :: type_i
     integer                     :: type_o
     logical                     :: cowling_approx
   contains 
     private
     procedure, public :: build_i
     procedure         :: build_regular_i_
     procedure         :: build_zero_i_
     procedure, public :: build_o
     procedure         :: build_zero_o_
     procedure         :: build_dziem_o_
     procedure         :: build_unno_o_
     procedure         :: build_jcd_o_
  end type nad_bound_t

  ! Interfaces

  interface nad_bound_t
     module procedure nad_bound_t_
  end interface nad_bound_t

  ! Access specifiers

  private

  public :: nad_bound_t

  ! Procedures

contains

  function nad_bound_t_ (ml, md_p, os_p) result (bd)

    class(model_t), pointer, intent(in) :: ml
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(nad_bound_t)                   :: bd

    ! Construct the nad_bound_t

    bd%ml => ml
    
    allocate(bd%rt, SOURCE=c_rot_t(ml, md_p, os_p))
    bd%vr = nad_vars_t(ml, md_p, os_p)

    select case (os_p%inner_bound)
    case ('REGULAR')
       bd%type_i = REGULAR_TYPE
    case ('ZERO')
       bd%type_i = ZERO_TYPE
    case default
       $ABORT(Invalid inner_bound)
    end select

    select case (os_p%outer_bound)
    case ('ZERO')
       bd%type_o = ZERO_TYPE
    case ('DZIEM')
       bd%type_o = DZIEM_TYPE
    case ('UNNO')
       bd%type_o = UNNO_TYPE
    case ('JCD')
       bd%type_o = JCD_TYPE
    case default
       $ABORT(Invalid outer_bound)
    end select

    bd%cowling_approx = os_p%cowling_approx

    bd%n_i = 3
    bd%n_o = 3

    bd%n_e = 6

    ! Finish

    return
    
  end function nad_bound_t_

  !****

  subroutine build_i (this, omega, B_i, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B_i(:,:)
    type(c_ext_t), intent(out)     :: scl

    $CHECK_BOUNDS(SIZE(B_i, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B_i, 2),this%n_e)
    
    ! Evaluate the inner boundary conditions

    select case (this%type_i)
    case (REGULAR_TYPE)
       call this%build_regular_i_(omega, B_i, scl)
    case (ZERO_TYPE)
       call this%build_zero_i_(omega, B_i, scl)
    case default
       $ABORT(Invalid type_i)
    end select

    ! Finish

    return

  end subroutine build_i

  !****

  subroutine build_regular_i_ (this, omega, B_i, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B_i(:,:)
    type(c_ext_t), intent(out)     :: scl

    real(WP)    :: c_1
    complex(WP) :: l_i
    complex(WP) :: omega_c
    real(WP)    :: alpha_gr

    $CHECK_BOUNDS(SIZE(B_i, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B_i, 2),this%n_e)
    
    $ASSERT(this%ml%x_i(1) == 0._WP,Boundary condition invalid for x /= 0)

    ! Evaluate the inner boundary conditions (regular-enforcing)

    associate (s => 1, &
               x => this%ml%x_i(1))

      ! Calculate coefficients

      c_1 = this%ml%c_1(s, x)

      l_i = this%rt%l_i(omega)

      omega_c = this%rt%omega_c(s, x, omega)

      if (this%cowling_approx) then
         alpha_gr = 0._WP
      else
         alpha_gr = 1._WP
      endif

      ! Set up the boundary conditions

      B_i(1,1) = c_1*omega_c**2
      B_i(1,2) = -l_i
      B_i(1,3) = alpha_gr*(0._WP)
      B_i(1,4) = alpha_gr*(0._WP)
      B_i(1,5) = 0._WP
      B_i(1,6) = 0._WP

      B_i(2,1) = alpha_gr*(0._WP)
      B_i(2,2) = alpha_gr*(0._WP)
      B_i(2,3) = alpha_gr*(l_i)
      B_i(2,4) = alpha_gr*(-1._WP) + (1._WP - alpha_gr)
      B_i(2,5) = alpha_gr*(0._WP)
      B_i(2,6) = alpha_gr*(0._WP)

      B_i(3,1) = 0._WP
      B_i(3,2) = 0._WP
      B_i(3,3) = alpha_gr*(0._WP)
      B_i(3,4) = alpha_gr*(0._WP)
      B_i(3,5) = 1._WP
      B_i(3,6) = 0._WP

      scl = c_ext_t(1._WP)

      ! Apply the variables transformation

      B_i = MATMUL(B_i, this%vr%H(s, x, omega))

    end associate

    ! Finish

    return

  end subroutine build_regular_i_

  !****

  subroutine build_zero_i_ (this, omega, B_i, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B_i(:,:)
    type(c_ext_t), intent(out)     :: scl

    real(WP) :: alpha_gr

    $CHECK_BOUNDS(SIZE(B_i, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B_i, 2),this%n_e)

    $ASSERT(this%ml%x_i(1) /= 0._WP,Boundary condition invalid for x == 0)

    ! Evaluate the inner boundary conditions (zero
    ! displacement/gravity)

    associate (s => 1, &
               x => this%ml%x_i(1))

      ! Calculate coefficients

      if (this%cowling_approx) then
         alpha_gr = 0._WP
      else
         alpha_gr = 1._WP
      endif

      ! Set up the boundary conditions

      B_i(1,1) = 1._WP
      B_i(1,2) = 0._WP
      B_i(1,3) = alpha_gr*(0._WP)
      B_i(1,4) = alpha_gr*(0._WP)
      B_i(1,5) = 0._WP
      B_i(1,6) = 0._WP
        
      B_i(2,1) = alpha_gr*(0._WP)
      B_i(2,2) = alpha_gr*(0._WP)
      B_i(2,3) = alpha_gr*(0._WP)
      B_i(2,4) = alpha_gr*(1._WP) + (1._WP - alpha_gr)
      B_i(2,5) = alpha_gr*(0._WP)
      B_i(2,6) = alpha_gr*(0._WP)
      
      B_i(3,1) = 0._WP
      B_i(3,2) = 0._WP
      B_i(3,3) = alpha_gr*(0._WP)
      B_i(3,4) = alpha_gr*(0._WP)
      B_i(3,5) = 1._WP
      B_i(3,6) = 0._WP

      scl = c_ext_t(1._WP)
      
      ! Apply the variables transformation

      B_i = MATMUL(B_i, this%vr%H(s, x, omega))

    end associate

    ! Finish

    return

  end subroutine build_zero_i_

  !****

  subroutine build_o (this, omega, B_o, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B_o(:,:)
    type(c_ext_t), intent(out)     :: scl

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)
    
    ! Evaluate the outer boundary conditions

    select case (this%type_o)
    case (ZERO_TYPE)
       call this%build_zero_o_(omega, B_o, scl)
    case (DZIEM_TYPE)
       call this%build_dziem_o_(omega, B_o, scl)
    case (UNNO_TYPE)
       call this%build_unno_o_(omega, B_o, scl)
    case (JCD_TYPE)
       call this%build_jcd_o_(omega, B_o, scl)
    case default
       $ABORT(Invalid type_o)
    end select

    ! Finish

    return

  end subroutine build_o
  
  !****

  subroutine build_zero_o_ (this, omega, B_o, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B_o(:,:)
    type(c_ext_t), intent(out)     :: scl

    real(WP)    :: V
    real(WP)    :: U
    real(WP)    :: nabla_ad
    complex(WP) :: l_e
    real(WP)    :: alpha_gr

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)

    ! Evaluate the outer boundary conditions (zero-pressure)

    associate (s => this%ml%n_s, &
               x => this%ml%x_o(this%ml%n_s))

      ! Calculate coefficients

      V = this%ml%V_2(s, x)*x**2
      U = this%ml%U(s, x)
      nabla_ad = this%ml%nabla_ad(s, x)

      l_e = this%rt%l_e(s, x, omega)

      if (this%cowling_approx) then
         alpha_gr = 0._WP
      else
         alpha_gr = 1._WP
      endif

      ! Set up the boundary conditions

      B_o(1,1) = 1._WP
      B_o(1,2) = -1._WP
      B_o(1,3) = alpha_gr*(1._WP)
      B_o(1,4) = alpha_gr*(0._WP)
      B_o(1,5) = 0._WP
      B_o(1,6) = 0._WP
      
      B_o(2,1) = alpha_gr*(U)
      B_o(2,2) = alpha_gr*(0._WP)
      B_o(2,3) = alpha_gr*(l_e + 1._WP) + (1._WP - alpha_gr)
      B_o(2,4) = alpha_gr*(1._WP)
      B_o(2,5) = alpha_gr*(0._WP)
      B_o(2,6) = alpha_gr*(0._WP)

      B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
      B_o(3,2) = 4._WP*nabla_ad*V
      B_o(3,3) = alpha_gr*(-4._WP*nabla_ad*V)
      B_o(3,4) = alpha_gr*(0._WP)
      B_o(3,5) = 4._WP
      B_o(3,6) = -1._WP

      scl = c_ext_t(1._WP)

      ! Apply the variables transformation

      B_o = MATMUL(B_o, this%vr%H(s, x, omega))

    end associate

    ! Finish

    return

  end subroutine build_zero_o_

  !****

  subroutine build_dziem_o_ (this, omega, B_o, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B_o(:,:)
    type(c_ext_t), intent(out)     :: scl

    real(WP)    :: V
    real(WP)    :: c_1
    real(WP)    :: nabla_ad
    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: omega_c
    real(WP)    :: alpha_gr

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)

    ! Evaluate the outer boundary conditions ([Dzi1971] formulation)

    associate (s => this%ml%n_s, &
               x => this%ml%x_o(this%ml%n_s))

      ! Calculate coefficients

      V = this%ml%V_2(s, x)*x**2
      c_1 = this%ml%c_1(s, x)
      nabla_ad = this%ml%nabla_ad(s, x)

      lambda = this%rt%lambda(s, x, omega)
      l_e = this%rt%l_e(s, x, omega)
    
      omega_c = this%rt%omega_c(s, x, omega)

      if (this%cowling_approx) then
         alpha_gr = 0._WP
      else
         alpha_gr = 1._WP
      endif

      ! Set up the boundary conditions

      B_o(1,1) = 1._WP + (lambda/(c_1*omega_c**2) - 4._WP - c_1*omega_c**2)/V
      B_o(1,2) = -1._WP
      B_o(1,3) = alpha_gr*(1._WP + (lambda/(c_1*omega_c**2) - l_e - 1._WP)/V)
      B_o(1,4) = alpha_gr*(0._WP)
      B_o(1,5) = 0._WP
      B_o(1,6) = 0._WP
     
      B_o(2,1) = alpha_gr*(0._WP)
      B_o(2,2) = alpha_gr*(0._WP)
      B_o(2,3) = alpha_gr*(l_e + 1._WP) + (1._WP - alpha_gr)
      B_o(2,4) = alpha_gr*(1._WP)
      B_o(2,5) = alpha_gr*(0._WP)
      B_o(2,6) = alpha_gr*(0._WP)

      B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
      B_o(3,2) = 4._WP*nabla_ad*V
      B_o(3,3) = alpha_gr*(-4._WP*nabla_ad*V)
      B_o(3,4) = alpha_gr*(0._WP)
      B_o(3,5) = 4._WP
      B_o(3,6) = -1._WP

      scl = c_ext_t(1._WP)

      ! Apply the variables transformation

      B_o = MATMUL(B_o, this%vr%H(s, x, omega))

    end associate

    ! Finish

    return

  end subroutine build_dziem_o_

  !****

  subroutine build_unno_o_ (this, omega, B_o, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B_o(:,:)
    type(c_ext_t), intent(out)     :: scl

    real(WP)    :: V_g
    real(WP)    :: V
    real(WP)    :: As
    real(WP)    :: c_1
    real(WP)    :: nabla_ad
    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: omega_c
    complex(WP) :: beta
    complex(WP) :: alpha_gr
    complex(WP) :: b_11
    complex(WP) :: b_12
    complex(WP) :: b_13
    complex(WP) :: b_21
    complex(WP) :: b_22
    complex(WP) :: b_23
    complex(WP) :: alpha_1
    complex(WP) :: alpha_2

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)

    ! Evaluate the outer boundary conditions ([Unn1989] formulation)

    associate (s => this%ml%n_s, &
               x => this%ml%x_o(this%ml%n_s))

      ! Calculate coefficients

      call eval_atmos_coeffs_unno(this%ml, V_g, As, c_1)

      V = this%ml%V_2(s, x)*x**2
      nabla_ad = this%ml%nabla_ad(s, x)

      lambda = this%rt%lambda(s, x, omega)
      l_e = this%rt%l_e(s, x, omega)

      omega_c = this%rt%omega_c(s, x, omega)

      beta = atmos_beta(V_g, As, c_1, omega_c, lambda)

      if (this%cowling_approx) then
         alpha_gr = 0._WP
      else
         alpha_gr = 1._WP
      endif

      b_11 = V_g - 3._WP
      b_12 = lambda/(c_1*omega_c**2) - V_g
      b_13 = alpha_gr*(V_g)
      
      b_21 = c_1*omega_c**2 - As
      b_22 = 1._WP + As
      b_23 = alpha_gr*(-As)
    
      alpha_1 = (b_12*b_23 - b_13*(b_22+l_e))/((b_11+l_e)*(b_22+l_e) - b_12*b_21)
      alpha_2 = (b_21*b_13 - b_23*(b_11+l_e))/((b_11+l_e)*(b_22+l_e) - b_12*b_21)

      ! Set up the boundary conditions

      B_o(1,1) = beta - b_11
      B_o(1,2) = -b_12
      B_o(1,3) = -(alpha_1*(beta - b_11) - alpha_2*b_12)
      B_o(1,4) = 0._WP
      B_o(1,5) = 0._WP
      B_o(1,6) = 0._WP

      B_o(2,1) = alpha_gr*(0._WP)
      B_o(2,2) = alpha_gr*(0._WP)
      B_o(2,3) = alpha_gr*(l_e + 1._WP) + (1._WP - alpha_gr)
      B_o(2,4) = alpha_gr*(1._WP)
      B_o(2,5) = alpha_gr*(0._WP)
      B_o(2,6) = alpha_gr*(0._WP)
    
      B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
      B_o(3,2) = 4._WP*nabla_ad*V
      B_o(3,3) = alpha_gr*(-4._WP*nabla_ad*V)
      B_o(3,4) = alpha_gr*(0._WP)
      B_o(3,5) = 4._WP
      B_o(3,6) = -1._WP
 
      scl = c_ext_t(1._WP)

      ! Apply the variables transformation

      B_o = MATMUL(B_o, this%vr%H(s, x, omega))

    end associate

    ! Finish

    return

  end subroutine build_unno_o_

  !****

  subroutine build_jcd_o_ (this, omega, B_o, scl)

    class(nad_bound_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP), intent(out)       :: B_o(:,:)
    type(c_ext_t), intent(out)     :: scl

    real(WP)    :: V_g
    real(WP)    :: V
    real(WP)    :: As
    real(WP)    :: c_1
    real(WP)    :: nabla_ad
    complex(WP) :: lambda
    complex(WP) :: l_e
    complex(WP) :: omega_c
    complex(WP) :: beta
    real(WP)    :: alpha_gr
    complex(WP) :: b_11
    complex(WP) :: b_12

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)

    ! Evaluate the outer boundary conditions ([Chr2008] formulation)

    associate (s => this%ml%n_s, &
               x => this%ml%x_o(this%ml%n_s))

      ! Calculate coefficients

      call eval_atmos_coeffs_jcd(this%ml, V_g, As, c_1)

      V = this%ml%V_2(s, x)*x**2
      nabla_ad = this%ml%nabla_ad(s, x)
    
      lambda = this%rt%lambda(s, x, omega)
      l_e = this%rt%l_e(s, x, omega)

      omega_c = this%rt%omega_c(s, x, omega)

      beta = atmos_beta(V_g, As, c_1, omega_c, lambda)

      if (this%cowling_approx) then
         alpha_gr = 0._WP
      else
         alpha_gr = 1._WP
      endif

      b_11 = V_g - 3._WP
      b_12 = lambda/(c_1*omega_c**2) - V_g

      ! Set up the boundary conditions

      B_o(1,1) = beta - b_11
      B_o(1,2) = -b_12
      B_o(1,3) = alpha_gr*(b_12 + (lambda/(c_1*omega_c**2) - l_e - 1._WP)*b_12/(V_g + As))
      B_o(1,4) = alpha_gr*(0._WP)
      B_o(1,5) = 0._WP
      B_o(1,6) = 0._WP
    
      B_o(2,1) = alpha_gr*(0._WP)
      B_o(2,2) = alpha_gr*(0._WP)
      B_o(2,3) = alpha_gr*(l_e + 1._WP) + (1._WP - alpha_gr)
      B_o(2,4) = alpha_gr*(1._WP)
      B_o(2,5) = alpha_gr*(0._WP)
      B_o(2,6) = alpha_gr*(0._WP)
    
      B_o(3,1) = 2._WP - 4._WP*nabla_ad*V
      B_o(3,2) = 4._WP*nabla_ad*V
      B_o(3,3) = alpha_gr*(-4._WP*nabla_ad*V)
      B_o(3,4) = alpha_gr*(0._WP)
      B_o(3,5) = 4._WP
      B_o(3,6) = -1._WP

      scl = c_ext_t(1._WP)

      ! Apply the variables transformation

      B_o = MATMUL(B_o, this%vr%H(s, x, omega))

    end associate

    ! Finish

    return

  end subroutine build_jcd_o_

end module gyre_nad_bound
