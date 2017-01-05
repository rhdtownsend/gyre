! Incfile  : gyre_rad_bound
! Purpose  : adiabatic radial boundary conditions
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

module gyre_rad_bound

  ! Uses

  use core_kinds

  use gyre_atmos
  use gyre_bound
  use gyre_ext
  use gyre_grid
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rad_vars
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

  type, extends (r_bound_t) :: rad_bound_t
     private
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     type(rad_vars_t)            :: vr
     type(point_t)               :: pt_i
     type(point_t)               :: pt_o
     real(WP)                    :: alpha_om
     integer                     :: type_i
     integer                     :: type_o
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
  end type rad_bound_t

  ! Interfaces

  interface rad_bound_t
     module procedure rad_bound_t_
  end interface rad_bound_t

  ! Access specifiers

  private

  public :: rad_bound_t

  ! Procedures

contains

  function rad_bound_t_ (ml, gr, md_p, os_p) result (bd)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(rad_bound_t)                   :: bd

    ! Construct the ad_bound_t

    bd%ml => ml
    
    allocate(bd%rt, SOURCE=r_rot_t(ml, gr, md_p, os_p))
    bd%vr = rad_vars_t(ml, gr, md_p, os_p)

    bd%pt_i = gr%pt(1)
    bd%pt_o = gr%pt(gr%n_k)

    select case (os_p%inner_bound)
    case ('REGULAR')
       $ASSERT(bd%pt_i%x == 0._WP,Boundary condition invalid for x /= 0)
       bd%type_i = REGULAR_TYPE
    case ('ZERO')
       $ASSERT(bd%pt_i%x /= 0._WP,Boundary condition invalid for x == 0)
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

    select case (os_p%time_factor)
    case ('OSC')
       bd%alpha_om = 1._WP
    case ('EXP')
       bd%alpha_om = -1._WP
    case default
       $ABORT(Invalid time_factor)
    end select

    bd%n_i = 1
    bd%n_o = 1

    bd%n_e = 2

    ! Finish

    return
    
  end function rad_bound_t_

  !****

  subroutine build_i (this, omega, B_i, scl_i)

    class(rad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: omega
    real(WP), intent(out)          :: B_i(:,:)
    type(r_ext_t), intent(out)     :: scl_i(:)

    $CHECK_BOUNDS(SIZE(B_i, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B_i, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl_i),this%n_i)
    
    ! Evaluate the inner boundary conditions

    select case (this%type_i)
    case (REGULAR_TYPE)
       call this%build_regular_i_(omega, B_i, scl_i)
    case (ZERO_TYPE)
       call this%build_zero_i_(omega, B_i, scl_i)
    case default
       $ABORT(Invalid type_i)
    end select

    ! Finish

    return

  end subroutine build_i

  !****

  subroutine build_regular_i_ (this, omega, B_i, scl_i)

    class(rad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: omega
    real(WP), intent(out)          :: B_i(:,:)
    type(r_ext_t), intent(out)     :: scl_i(:)

    real(WP) :: c_1
    real(WP) :: omega_c
    real(WP) :: alpha_om

    $CHECK_BOUNDS(SIZE(B_i, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B_i, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(scl_i),this%n_i)

    ! Evaluate the inner boundary conditions (regular-enforcing)

    associate (pt => this%pt_i)

      ! Calculate coefficients

      c_1 = this%ml%c_1(pt)

      omega_c = this%rt%omega_c(pt, omega)

      alpha_om = this%alpha_om

      ! Set up the boundary conditions

      B_i(1,1) = c_1*alpha_om*omega_c**2
      B_i(1,2) = 0._WP

      scl_i = r_ext_t(1._WP)

      ! Apply the variables transformation

      B_i = MATMUL(B_i, this%vr%H(pt, omega))

    end associate

    ! Finish

    return

  end subroutine build_regular_i_

  !****

  subroutine build_zero_i_ (this, omega, B_i, scl_i)

    class(rad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: omega
    real(WP), intent(out)          :: B_i(:,:)
    type(r_ext_t), intent(out)     :: scl_i(:)

    $CHECK_BOUNDS(SIZE(B_i, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B_i, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl_i),this%n_i)

    ! Evaluate the inner boundary conditions (zero displacement)

    associate (pt => this%pt_i)

      ! Set up the boundary conditions

      B_i(1,1) = 1._WP
      B_i(1,2) = 0._WP

      scl_i = r_ext_t(1._WP)
      
      ! Apply the variables transformation

      B_i = MATMUL(B_i, this%vr%H(pt, omega))

    end associate

    ! Finish

    return

  end subroutine build_zero_i_

  !****

  subroutine build_o (this, omega, B_o, scl_o)

    class(rad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: omega
    real(WP), intent(out)          :: B_o(:,:)
    type(r_ext_t), intent(out)     :: scl_o(:)

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(scl_o),this%n_o)

    ! Evaluate the outer boundary conditions

    select case (this%type_o)
    case (ZERO_TYPE)
       call this%build_zero_o_(omega, B_o, scl_o)
    case (DZIEM_TYPE)
       call this%build_dziem_o_(omega, B_o, scl_o)
    case (UNNO_TYPE)
       call this%build_unno_o_(omega, B_o, scl_o)
    case (JCD_TYPE)
       call this%build_jcd_o_(omega, B_o, scl_o)
    case default
       $ABORT(Invalid type_o)
    end select

    ! Finish

    return

  end subroutine build_o
  
  !****

  subroutine build_zero_o_ (this, omega, B_o, scl_o)

    class(rad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: omega
    real(WP), intent(out)          :: B_o(:,:)
    type(r_ext_t), intent(out)     :: scl_o(:)

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl_o),this%n_o)

    ! Evaluate the outer boundary conditions (zero-pressure)

    associate (pt => this%pt_o)
    
      ! Set up the boundary conditions

      B_o(1,1) = 1._WP
      B_o(1,2) = -1._WP
      
      scl_o = r_ext_t(1._WP)
    
      ! Apply the variables transformation

      B_o = MATMUL(B_o, this%vr%H(pt, omega))

    end associate

    ! Finish

    return

  end subroutine build_zero_o_

  !****

  subroutine build_dziem_o_ (this, omega, B_o, scl_o)

    class(rad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: omega
    real(WP), intent(out)          :: B_o(:,:)
    type(r_ext_t), intent(out)     :: scl_o(:)

    real(WP) :: V
    real(WP) :: c_1
    real(WP) :: omega_c
    real(WP) :: alpha_om

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl_o),this%n_o)

    ! Evaluate the outer boundary conditions ([Dzi1971] formulation)

    associate (pt => this%pt_o)

      if (this%ml%vacuum(pt)) then

         ! For a vacuum, the boundary condition reduces to the zero
         ! condition

         call this%build_zero_o_(omega, B_o, scl_o)

      else

         ! Calculate coefficients

         V = this%ml%V_2(pt)*pt%x**2
         c_1 = this%ml%c_1(pt)

         omega_c = this%rt%omega_c(pt, omega)

         alpha_om = this%alpha_om

         ! Set up the boundary conditions
        
         B_o(1,1) = 1 - (4._WP + c_1*alpha_om*omega_c**2)/V
         B_o(1,2) = -1._WP

         scl_o = r_ext_t(1._WP)

         ! Apply the variables transformation

         B_o = MATMUL(B_o, this%vr%H(pt, omega))

      endif

    end associate

    ! Finish

    return

  end subroutine build_dziem_o_

  !****
  
  subroutine build_unno_o_ (this, omega, B_o, scl_o)

    class(rad_bound_t), intent(in) :: this
    real(WP), intent(in)           :: omega
    real(WP), intent(out)          :: B_o(:,:)
    type(r_ext_t), intent(out)     :: scl_o(:)

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: c_1
    real(WP) :: omega_c
    real(WP) :: beta
    real(WP) :: b_11
    real(WP) :: b_12

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl_o),this%n_o)

    ! Evaluate the outer boundary conditions ([Unn1989] formulation)

    associate (pt => this%pt_o)

      if (this%ml%vacuum(pt)) then

         ! For a vacuum, the boundary condition reduces to the zero
         ! condition

         call this%build_zero_o_(omega, B_o, scl_o)

      else

         ! Calculate coefficients

         call eval_atmos_coeffs_unno(this%ml, pt, V_g, As, c_1)

         omega_c = this%rt%omega_c(pt, omega)

         beta = atmos_beta(V_g, As, c_1, omega_c, 0._WP)
      
         b_11 = V_g - 3._WP
         b_12 = -V_g
    
         ! Set up the boundary conditions
      
         B_o(1,1) = beta - b_11
         B_o(1,2) = -b_12

         scl_o = r_ext_t(1._WP)

         ! Apply the variables transformation

         B_o = MATMUL(B_o, this%vr%H(pt, omega))

      endif

    end associate

    ! Finish

    return

  end subroutine build_unno_o_

  !****

  subroutine build_jcd_o_ (this, omega, B_o, scl_o)

    class(rad_bound_t), intent(in) :: this
    real(WP), intent(in)          :: omega
    real(WP), intent(out)         :: B_o(:,:)
    type(r_ext_t), intent(out)    :: scl_o(:)

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: c_1
    real(WP) :: omega_c
    real(WP) :: beta
    real(WP) :: b_11
    real(WP) :: b_12

    $CHECK_BOUNDS(SIZE(B_o, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B_o, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl_o),this%n_o)

    ! Evaluate the outer boundary conditions ([Chr2008] formulation)

    ! Calculate coefficients

    associate (pt => this%pt_o)

      if (this%ml%vacuum(pt)) then

         ! For a vacuum, the boundary condition reduces to the zero
         ! condition

         call this%build_zero_o_(omega, B_o, scl_o)

      else

         call eval_atmos_coeffs_jcd(this%ml, pt, V_g, As, c_1)

         omega_c = this%rt%omega_c(pt, omega)

         beta = atmos_beta(V_g, As, c_1, omega_c, 0._WP)

         b_11 = V_g - 3._WP
         b_12 = -V_g

         ! Set up the boundary conditions

         B_o(1,1) = beta - b_11
         B_o(1,2) = -b_12

         scl_o = r_ext_t(1._WP)

         ! Apply the variables transformation

         B_o = MATMUL(B_o, this%vr%H(pt, omega))

      endif

    end associate

    ! Finish

    return

  end subroutine build_jcd_o_

end module gyre_rad_bound
