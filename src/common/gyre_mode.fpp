! Module   : gyre_mode
! Purpose  : mode data
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
$include 'core_memory.inc'

module gyre_mode

  ! Uses

  use core_kinds
  use gyre_constants
  use core_parallel

  use gyre_model
  $if ($MPI)
  use gyre_model_mpi
  $endif
  use gyre_ext_arith
  use gyre_modepar
  use gyre_oscpar
  use gyre_grid
  use gyre_mode_funcs
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: mode_t
     class(model_t), pointer  :: ml => null()
     type(modepar_t)          :: mp
     type(oscpar_t)           :: op
     type(ext_complex_t)      :: discrim
     type(ext_real_t)         :: chi
     real(WP), allocatable    :: x(:)
     complex(WP), allocatable :: y(:,:)
     real(WP)                 :: x_ref
     complex(WP)              :: y_ref(6)
     complex(WP)              :: omega
     integer                  :: n_pg
     integer                  :: n_p
     integer                  :: n_g
     integer                  :: n
     integer                  :: n_iter
     logical                  :: pruned
   contains
     private
     $if($GFORTRAN_PR57922)
     procedure, public :: final => final_
     $endif
     procedure, public :: prune => prune_
     procedure, public :: freq => freq_
     procedure, public :: xi_r => xi_r_
     procedure, public :: xi_h => xi_h_
     procedure, public :: eul_phi => eul_phi_
     procedure, public :: deul_phi => deul_phi_
     procedure, public :: lag_S => lag_S_
     procedure, public :: lag_L => lag_L_
     procedure, public :: eul_P => eul_P_
     procedure, public :: lag_P => lag_P_
     procedure, public :: eul_rho => eul_rho_
     procedure, public :: lag_rho => lag_rho_
     procedure, public :: eul_T => eul_T_
     procedure, public :: lag_T => lag_T_
     procedure, public :: dE_dx => dE_dx_
     procedure, public :: dW_dx => dW_dx_
     procedure, public :: F_j_rey => F_j_rey_
     procedure, public :: Yt_1 => Yt_1_
     procedure, public :: Yt_2 => Yt_2_
     procedure, public :: I_0 => I_0_
     procedure, public :: I_1 => I_1_
     procedure, public :: xi_r_ref => xi_r_ref_
     procedure, public :: xi_h_ref => xi_h_ref_
     procedure, public :: eul_phi_ref => eul_phi_ref_
     procedure, public :: deul_phi_ref => deul_phi_ref_
     procedure, public :: lag_T_eff => lag_T_eff_
     procedure, public :: lag_g_eff => lag_g_eff_
     procedure, public :: E => E_
     procedure, public :: E_norm => E_norm_
     procedure, public :: W => W_
     procedure, public :: K => K_
     procedure, public :: beta => beta_
     procedure, public :: omega_int => omega_int_
     procedure, public :: eta => eta_
     procedure, public :: prop_type => prop_type_
     procedure         :: classify_
     procedure, public :: lag_S_en => lag_S_en_
     procedure, public :: lag_L_rd => lag_L_rd_
  end type mode_t

  ! Interfaces

  interface mode_t
     module procedure mode_t_
  end interface mode_t

  interface reallocate
     module procedure reallocate_1_
  end interface reallocate

  $if ($MPI)
  interface bcast
     module procedure bcast_
  end interface bcast
  $endif

  ! Access specifiers

  private

  public :: mode_t
  public :: reallocate
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  function mode_t_ (ml, mp, op, omega, discrim, x, y, x_ref, y_ref) result (md)

    class(model_t), pointer, intent(in) :: ml
    type(modepar_t), intent(in)         :: mp
    type(oscpar_t), intent(in)          :: op
    complex(WP), intent(in)             :: omega
    type(ext_complex_t), intent(in)     :: discrim
    real(WP), intent(in)                :: x(:)
    complex(WP), intent(in)             :: y(:,:)
    real(WP), intent(in)                :: x_ref
    complex(WP), intent(in)             :: y_ref(:)
    type(mode_t)                        :: md

    real(WP)    :: phase
    complex(WP) :: norm_fac
    integer     :: n_p
    integer     :: n_g
    integer     :: n_pg

    $CHECK_BOUNDS(SIZE(y, 1),6)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    ! Construct the mode_t

    md%ml => ml

    md%mp = mp
    md%op = op

    md%x = x
    md%y = y

    md%x_ref = x_ref
    md%y_ref = y_ref

    md%omega = omega
    md%discrim = discrim

    md%n = SIZE(md%x)

    md%pruned = .FALSE.

    ! Normalize by the mode inertia, and so that y_ref(1) is real

    phase = ATAN2(AIMAG(md%y_ref(1)), REAL(md%y_ref(1)))

    norm_fac = 1._WP/SQRT(integrate(md%x, md%dE_dx()))*EXP(CMPLX(0._WP, -phase, KIND=WP))

    md%y = norm_fac*md%y
    md%y_ref = norm_fac*md%y_ref

    ! Classify the mode

    call classify_(md, n_p, n_g, n_pg)

    md%n_p = n_p
    md%n_g = n_g
    md%n_pg = n_pg

    ! Finish

    return

  end function mode_t_

!****

  $REALLOCATE(type(mode_t),1)

!****

  $if($GFORTRAN_PR57922)

  subroutine final_ (this)

    class(mode_t), intent(inout) :: this

    ! Finalize the mode_t

    if (.NOT. this%pruned) then

       deallocate(this%x)
       deallocate(this%y)

    endif

    ! Finish

  end subroutine final_

  $endif

!****

  subroutine prune_ (this)

    class(mode_t), intent(inout) :: this

    ! Prune the mode_t

    if (.NOT. this%pruned) then

       deallocate(this%x)
       deallocate(this%y)

       this%pruned = .TRUE.

    endif

    ! Finish

    return

  end subroutine prune_

!****

  function freq_ (this, freq_units) result (freq)

    class(mode_t), intent(in)    :: this
    character(LEN=*), intent(in) :: freq_units
    complex(WP)                  :: freq

    ! Calculate the frequency

    freq = this%omega*freq_scale(this%ml, this%mp, this%op, this%x(this%n), freq_units)

    ! Finish
    
    return

  end function freq_

!****

  $define $CALC_GRID $sub

  $local $NAME $1
  $local $TYPE $2

  function ${NAME}_ (this)

    class(mode_t), intent(in) :: this
    $TYPE(WP)                 :: ${NAME}_(this%n)

    integer :: k
    
    $ASSERT(.NOT. this%pruned,Cannot calculate ${NAME} from pruned solution)
    
    ! Calculate $NAME on the x grid

    !$OMP PARALLEL DO
    do k = 1, this%n
       ${NAME}_(k) = ${NAME}(this%ml, this%mp, this%omega, this%x(k), this%y(:,k))
    end do
    
    ! Finish

    return

  end function ${NAME}_

  $endsub

  $CALC_GRID(xi_r,complex)
  $CALC_GRID(xi_h,complex)
  $CALC_GRID(eul_phi,complex)
  $CALC_GRID(deul_phi,complex)
  $CALC_GRID(lag_S,complex)
  $CALC_GRID(lag_L,complex)
  $CALC_GRID(lag_P,complex)
  $CALC_GRID(eul_P,complex)
  $CALC_GRID(eul_rho,complex)
  $CALC_GRID(lag_rho,complex)
  $CALC_GRID(eul_T,complex)
  $CALC_GRID(lag_T,complex)
  $CALC_GRID(dE_dx,real)
  $CALC_GRID(dW_dx,real)
  $CALC_GRID(F_j_rey,real)
  $CALC_GRID(Yt_1,complex)
  $CALC_GRID(Yt_2,complex)
  $CALC_GRID(I_0,complex)
  $CALC_GRID(I_1,complex)

!****

  function f_T_ (this) result (f_T)

    class(mode_t), intent(in) :: this
    real(WP)                  :: f_T

    complex(WP) :: C_T

    ! Calculate the non-adiabatic f_T parameter. This is expression is
    ! based on eqn. 5 of [Dup2003]

    C_T = this%lag_T_eff()/this%xi_r_ref()

    f_T = ABS(C_T)

    ! Finish

    return

  end function f_T_

!****

  function f_g_ (this) result (f_g)

    class(mode_t), intent(in) :: this
    real(WP)                  :: f_g

    complex(WP) :: C_g

    ! Calculate the non-adiabatic f_g parameter. This is expression is
    ! based on eqn. 6 of [Dup2003]

    C_g = this%lag_g_eff()/this%xi_r_ref()

    f_g = -ABS(C_g)

    ! Finish

    return

  end function f_g_

!****

  function psi_T_ (this) result (psi_T)

    class(mode_t), intent(in) :: this
    real(WP)                  :: psi_T

    complex(WP) :: C_T

    ! Calculate the non-adiabatic psi_T parameter, in radians. This is
    ! expression is based on eqn. 5 of [Dup2003]

    C_T = this%lag_T_eff()/this%xi_r_ref()

    psi_T = ATAN2(AIMAG(C_T), REAL(C_T))

    ! Finish

    return

  end function psi_T_

!****

  $define $CALC_REF $sub

  $local $NAME $1
  $local $TYPE $2

  function ${NAME}_ref_ (this)

    class(mode_t), intent(in) :: this
    $TYPE(WP)                 :: ${NAME}_ref_

    integer :: k
    
    ! Calculate $NAME at x_ref

    ${NAME}_ref_ = ${NAME}(this%ml, this%mp, this%omega, this%x_ref, this%y_ref)
    
    ! Finish

    return

  end function ${NAME}_ref_

  $endsub

  $CALC_REF(xi_r,complex)
  $CALC_REF(xi_h,complex)
  $CALC_REF(eul_phi,complex)
  $CALC_REF(deul_phi,complex)

!****

  function lag_T_eff_ (this) result (lag_T_eff)

    class(mode_t), intent(in) :: this
    complex(WP)               :: lag_T_eff

    ! Calculate the effective temperature perturbation at x_ref
    ! (assumed to correspond to the photosphere), in units of
    ! T_eff. This expression is based on the standard definition of
    ! effective temperature

    lag_T_eff = 0.25_WP*(lag_L(this%ml, this%mp, this%omega, this%x_ref, this%y_ref) - &
                  2._WP*xi_r(this%ml, this%mp, this%omega, this%x_ref, this%y_ref))

    ! Finish

    return

  end function lag_T_eff_

!****

  function lag_g_eff_ (this) result (lag_g_eff)

    class(mode_t), intent(in) :: this
    complex(WP)               :: lag_g_eff

    ! Calculate the effective gravity perturbation at x_ref (assumed
    ! to correspond to the photosphere), in units of the gravity. This
    ! expression is based on eqn. 24 of [Dup2002]

    associate (c_1 => this%ml%c_1(this%x_ref), U => this%ml%U(this%x_ref), &
               x => this%x_ref)

      lag_g_eff = (c_1/x)*this%deul_phi_ref() + (U - (2._WP + c_1*this%omega**2))*this%xi_r_ref()/x

    end associate

    ! Finish

    return

  end function lag_g_eff_

!****

  function E_ (this) result (E)

    class(mode_t), intent(in) :: this
    real(WP)                  :: E
    
    ! Calculate the total mode inertia, in units of M_star
    ! R_star**2. Due to the mode normalization enforced in mode_t_,
    ! this is always unity

    E = 1._WP

    ! Finish

    return

  end function E_

!****

  function E_norm_ (this) result (E_norm)
 
    class(mode_t), intent(in) :: this
    real(WP)                  :: E_norm

    real(WP)    :: E
    complex(WP) :: xi_r(this%n)
    complex(WP) :: xi_h(this%n)
    complex(WP) :: xi_r_1
    complex(WP) :: xi_h_1
    real(WP)    :: A2

    ! Calculate the normalized mode inertia. This expression is based
    ! on eqn. 3.140 of [Aer2010]

    E = this%E()

    associate(l => this%mp%l)

      select case (this%op%inertia_norm_type)
      case ('RADIAL')
         A2 = ABS(this%xi_r_ref())**2
      case ('HORIZ')
         A2 = l*(l+1)*ABS(this%xi_h_ref())**2
      case ('BOTH')
         A2 = ABS(this%xi_r_ref())**2 + l*(l+1)*ABS(this%xi_h_ref())**2
      case default
         $ABORT(Invalid inertia_norm_type)
      end select

      if(A2 == 0._WP) then
         $WARN(Amplitude at x_ref is zero; not normalizing inertia)
         E_norm = E
      else
         E_norm = E/A2
      endif

    end associate

    ! Finish

    return

  end function E_norm_

!****

  function W_ (this) result (W)

    class(mode_t), intent(in) :: this
    real(WP)                  :: W
    
    ! Calculate the total work by integrating the differential work

    W = integrate(this%x, this%dW_dx())

    ! Finish

    return

  end function W_

!****

  function K_ (this) result (K)

    class(mode_t), intent(in) :: this
    complex(WP)               :: K(this%n)

    complex(WP) :: xi_r(this%n)
    complex(WP) :: xi_h(this%n)

    ! Calculate the rotation splitting kernel

    xi_r = this%xi_r()
    xi_h = this%xi_h()

    associate(x => this%x, U => this%ml%U(this%x), c_1 => this%ml%c_1(this%x), &
              l => this%mp%l)

      K = (ABS(xi_r)**2 + (l*(l+1)-1)*ABS(xi_h)**2 - 2._WP*xi_r*CONJG(xi_h))*U*x**2/c_1

      K = K/integrate(x, K)

    end associate

    ! Finish

    return

  end function K_

!****

  function beta_ (this) result (beta)

    class(mode_t), intent(in) :: this
    complex(WP)               :: beta

    complex(WP) :: xi_r(this%n)
    complex(WP) :: xi_h(this%n)

    ! Calculate the rotation splitting scale

    xi_r = this%xi_r()
    xi_h = this%xi_h()

    associate(x => this%x, U => this%ml%U(this%x), c_1 => this%ml%c_1(this%x), &
              l => this%mp%l)

      beta = integrate(x, (ABS(xi_r)**2 + (l*(l+1)-1)*ABS(xi_h)**2 - 2._WP*xi_r*CONJG(xi_h))*U*x**2/c_1) / &
             integrate(x, (ABS(xi_r)**2 + l*(l+1)*ABS(xi_h)**2)*U*x**2/c_1)

    end associate

    ! Finish

    return

  end function beta_

!****

  function omega_int_ (this) result (omega_int)

    class(mode_t), intent(in) :: this
    complex(WP)               :: omega_int

    real(WP)    :: x4_V(this%n)
    complex(WP) :: W_th
    complex(WP) :: W_re
    complex(WP) :: W_gr
    complex(WP) :: W_xi

    ! Calculate the dimensionless frequency from the integral
    ! expression in eqn. (1.71) of [Dup2003]

    associate(x => this%x, &
              V => this%ml%V(this%x), V_g => this%ml%V(this%x)/this%ml%Gamma_1(this%x), &
              As => this%ml%As(this%x), U => this%ml%U(this%x), c_1 => this%ml%c_1(this%x), &
              xi_r => this%xi_r(), eul_phi => this%eul_phi(), &
              lag_rho => this%lag_rho(), eul_rho => this%eul_rho(), &
              lag_P => this%lag_P())

      where (x /= 0._WP)
         x4_V = x**4/V
      elsewhere
         x4_V = 0._WP
      endwhere

      ! Work-function-like terms

      W_th = integrate(x, CONJG(lag_rho)*lag_P*(U*x4_V/(c_1**2)))

      W_re = integrate(x, 2._WP*REAL(lag_rho*CONJG(xi_r)*(x/c_1)*(x**2*U/c_1)))

      W_gr = integrate(x, CONJG(eul_rho)*eul_phi*(x**2*U/c_1))

      W_xi = integrate(x, -ABS(xi_r)**2*(x/c_1)*(x*U*(-V_g-As)/c_1))

    end associate
    
    omega_int = SQRT((W_th + W_re + W_gr + W_xi)/this%E())

    ! Finish

    return

  end function omega_int_

!****

  function eta_ (this) result (eta)

    class(mode_t), intent(in) :: this
    real(WP)                  :: eta

    real(WP)    :: dW_dx(this%n)

    ! Calculate the normalized growth rate defined by [Stel1978]

    dW_dx = this%dW_dx()

    eta = integrate(this%x, dW_dx)/integrate(this%x, ABS(dW_dx))

    ! Finish

    return

  end function eta_
    
!****

  function prop_type_ (this) result (prop_type)

    class(mode_t), intent(in) :: this
    integer                   :: prop_type(this%n)

    ! Set up the propagation type (0 -> evanescent, 1 -> p, -1 -> g)

    associate(x => this%x, V_g => this%ml%V(this%x)/this%ml%Gamma_1(this%x), &
              As => this%ml%As(this%x), c_1 => this%ml%c_1(this%x), &
              l => this%mp%l, omega_c => REAL(this%ml%omega_c(this%x, this%mp%m, this%omega)))

      prop_type = MERGE(1, 0, c_1*omega_c**2 > As) + &
                   MERGE(-1, 0, l*(l+1)/(c_1*omega_c**2) > V_g)

    end associate

    ! Finish

    return

  end function prop_type_

!****

  subroutine classify_ (md, n_p, n_g, n_pg)

    class(mode_t), intent(in) :: md
    integer, intent(out)      :: n_p
    integer, intent(out)      :: n_g
    integer, intent(out)      :: n_pg

    real(WP) :: y_1(md%n)
    real(WP) :: y_2(md%n)
    real(WP) :: x_turn
    integer  :: i
    integer  :: n_c
    integer  :: n_a

    ! Classify the mode based on its eigenfunctions

    select case (md%mp%l)

    case (0)

       ! Radial modes
       
       ! Look for the first monotonic segment in y_1 (this is to deal with
       ! noisy near-zero solutions at the origin)

       y_1 = REAL(md%y(1,:))
       y_2 = REAL(md%y(2,:))

       mono_loop : do i = 2,md%n-1
          if((y_1(i) >= y_1(i-1) .AND. y_1(i+1) >= y_1(i)) .OR. &
             (y_1(i) <= y_1(i-1) .AND. y_1(i+1) <= y_1(i))) exit mono_loop
       end do mono_loop

       ! Count winding numbers

       call count_windings_(y_1(i:), y_2(i:), n_c, n_a)

       ! Classify (the additional 1 is for the node at the center)

       n_p = n_a + 1
       n_g = n_c

       n_pg = n_p - n_g

    case (1)

       ! Dipole modes

       ! Set up the Takata Y^a_1 and Y^a_2 functions
       
       y_1 = REAL(md%Yt_1())
       y_2 = REAL(md%Yt_2())

       ! Find the inner turning point (this is to deal with noisy
       ! near-zero solutions at the origin)

       call find_x_turn(md%x, md%ml, md%mp, REAL(md%omega), x_turn)

       x_turn_loop : do i = 1,md%n-1
          if(md%x(i) > x_turn) exit x_turn_loop
       end do x_turn_loop

       ! Count winding numbers

!       call count_windings(y_1(i:), y_2(i:), n_c, n_a, md%x)
       call count_windings_(y_1(i:), y_2(i:), n_c, n_a)

       n_p = n_a
       n_g = n_c

       if(n_p >= n_g) then
          n_pg = n_p - n_g + 1
       else
          n_pg = n_p - n_g
       endif

    case default

       ! Other modes

       y_1 = REAL(md%y(1,:))
       y_2 = REAL(md%y(2,:))

       ! Count winding numbers

       call count_windings_(y_1, y_2, n_c, n_a)

       ! Classify

       n_p = n_a
       n_g = n_c

       n_pg = n_p - n_g

    end select

    ! Finish

    return

  contains

    subroutine count_windings_ (y_1, y_2, n_c, n_a, x)

      real(WP), intent(in)           :: y_1(:)
      real(WP), intent(in)           :: y_2(:)
      integer, intent(out)           :: n_c
      integer, intent(out)           :: n_a
      real(WP), optional, intent(in) :: x(:)

      integer  :: i
      real(WP) :: y_2_cross

      $CHECK_BOUNDS(SIZE(y_2),SIZE(y_1))

      if(PRESENT(x)) then
         $CHECK_BOUNDS(SIZE(x),SIZE(y_1))
      endif

      ! Count clockwise (n_c) and anticlockwise (n_a) windings in the (y_1,y_2) plane

      n_c = 0
      n_a = 0

      do i = 1,SIZE(y_1)-1

         ! Look for a node in y_1

         if(y_1(i) >= 0._WP .AND. y_1(i+1) < 0._WP) then

            ! Solve for the crossing ordinate

            y_2_cross = y_2(i) - y_1(i)*(y_2(i+1) - y_2(i))/(y_1(i+1) - y_1(i))

            if(y_2_cross >= 0._WP) then
               n_a = n_a + 1
               if(PRESENT(x)) print *,'A node:',x(i),x(i+1)
            else
               n_c = n_c + 1
               if(PRESENT(x)) print *,'C node:',x(i),x(i+1)
            endif

         elseif(y_1(i) <= 0._WP .AND. y_1(i+1) > 0._WP) then

            ! Solve for the crossing ordinate

            y_2_cross = y_2(i) - y_1(i)*(y_2(i+1) - y_2(i))/(y_1(i+1) - y_1(i))

            if(y_2_cross <= 0._WP) then
               n_a = n_a + 1
               if(PRESENT(x)) print *,'A node:',x(i),x(i+1)
            else
               n_c = n_c + 1
               if(PRESENT(x)) print *,'C node:',x(i),x(i+1)
            endif

         endif

      end do

      ! Finish

      return

    end subroutine count_windings_

  end subroutine classify_

!****

  function lag_S_en_ (this) result (lag_S_en)

    class(mode_t), intent(in) :: this
    complex(WP)               :: lag_S_en(this%n)

    complex(WP) :: A_6(6,this%n)
    complex(WP) :: dy_6(this%n)
    complex(WP) :: y_5(this%n)

    ! Calculate the Lagrangian specific entropy perturbation in units
    ! of c_p, from the energy equation

    associate(V => this%ml%V(this%x), c_1 => this%ml%c_1(this%x), &
              nabla_ad => this%ml%nabla_ad(this%x), nabla => this%ml%nabla(this%x), &
              c_rad => this%ml%c_rad(this%x), dc_rad => this%ml%dc_rad(this%x), c_thm => this%ml%c_thm(this%x), &
              c_eps_ad => this%ml%c_eps_ad(this%x), c_eps_S => this%ml%c_eps_S(this%x), &              
              l => this%mp%l, omega_c => this%ml%omega_c(this%x, this%mp%m, this%omega))

      where(this%x /= 0)
         A_6(1,:) = l*(l+1)*(nabla_ad/nabla - 1._WP)*c_rad - V*c_eps_ad
         A_6(2,:) = V*c_eps_ad - l*(l+1)*c_rad*(nabla_ad/nabla - (3._WP + dc_rad)/(c_1*omega_c**2))
         A_6(3,:) = l*(l+1)*nabla_ad/nabla*c_rad - V*c_eps_ad
         A_6(4,:) = 0._WP
         A_6(5,:) = c_eps_S - l*(l+1)*c_rad/(nabla*V) - (0._WP,1._WP)*omega_c*c_thm
         A_6(6,:) = -1._WP - l
      elsewhere
         A_6(1,:) = l*(l+1)*(nabla_ad/nabla - 1._WP)*c_rad - V*c_eps_ad
         A_6(2,:) = V*c_eps_ad - l*(l+1)*c_rad*(nabla_ad/nabla - (3._WP + dc_rad)/(c_1*omega_c**2))
         A_6(3,:) = l*(l+1)*nabla_ad/nabla*c_rad - V*c_eps_ad
         A_6(4,:) = 0._WP
         A_6(5,:) = -HUGE(0._WP)
         A_6(6,:) = -1._WP - l
      endwhere

      dy_6 = this%x*deriv_(this%x, this%y(6,:))

      y_5 = (dy_6 - (A_6(1,:)*this%y(1,:) + &
                     A_6(2,:)*this%y(2,:) + &
                     A_6(3,:)*this%y(3,:) + &
                     A_6(4,:)*this%y(4,:) + &
                     A_6(6,:)*this%y(6,:)))/A_6(5,:)

      where(this%x /= 0._WP)
         lag_S_en = y_5*this%x**(l-2)
      elsewhere
         lag_S_en = 0._WP
      end where

    end associate
         
    ! Finish

    return

  contains

    function deriv_ (x, y) result (dy_dx)

      real(WP), intent(in)    :: x(:)
      complex(WP), intent(in) :: y(:)
      complex(WP)             :: dy_dx(SIZE(x))
      
      integer :: n
      integer :: i

      $CHECK_BOUNDS(SIZE(y),SIZE(x))

      ! Differentiate y(x) using centered finite differences

      n = SIZE(x)

      dy_dx(1) = (y(2) - y(1))/(x(2) - x(1))

      do i = 2,n-1
         dy_dx(i) = 0.5_WP*((y(i) - y(i-1))/(x(i) - x(i-1)) + &
                            (y(i+1) - y(i))/(x(i+1) - x(i)))
      end do

      dy_dx(n) = (y(n) - y(n-1))/(x(n) - x(n-1))

      ! Finish

      return

    end function deriv_

  end function lag_S_en_

!****

  function lag_L_rd_ (this) result (lag_L_rd)

    class(mode_t), intent(in) :: this
    complex(WP)               :: lag_L_rd(this%n)

    complex(WP) :: A_5(6,this%n)
    complex(WP) :: dy_5(this%n)
    complex(WP) :: y_6(this%n)

    ! Calculate the Lagrangian luminosity perturbation in units of L_star, 
    ! from the radiative diffusion equation

    associate(V => this%ml%V(this%x), U => this%ml%U(this%x), c_1 => this%ml%c_1(this%x), &
              nabla_ad => this%ml%nabla_ad(this%x), nabla => this%ml%nabla(this%x), &
              c_dif => this%ml%c_dif(this%x), c_rad => this%ml%c_rad(this%x), &
              kappa_S => this%ml%kappa_S(this%x), &
              l => this%mp%l, omega_c => this%ml%omega_c(this%x, this%mp%m, this%omega))

      A_5(1,:) = V*(nabla_ad*(U - c_1*omega_c**2) - 4._WP*(nabla_ad - nabla) + c_dif)
      A_5(2,:) = V*(l*(l+1)/(c_1*omega_c**2)*(nabla_ad - nabla) - c_dif)
      A_5(3,:) = V*c_dif
      A_5(4,:) = V*nabla_ad
      A_5(5,:) = V*nabla*(4._WP - kappa_S) - (l - 2._WP)
      A_5(6,:) = -V*nabla/c_rad

      dy_5 = this%x*deriv_(this%x, this%y(5,:))

      where(this%x /= 0._WP)

         y_6 = (dy_5 - (A_5(1,:)*this%y(1,:) + &
                        A_5(2,:)*this%y(2,:) + &
                        A_5(3,:)*this%y(3,:) + &
                        A_5(4,:)*this%y(4,:) + &
                        A_5(5,:)*this%y(5,:)))/A_5(6,:)

      elsewhere

         y_6 = 0._WP

      endwhere

      lag_L_rd = y_6*this%x**(l+1)

    end associate
         
    ! Finish

    return

  contains

    function deriv_ (x, y) result (dy_dx)

      real(WP), intent(in)    :: x(:)
      complex(WP), intent(in) :: y(:)
      complex(WP)             :: dy_dx(SIZE(x))
      
      integer :: n
      integer :: i

      $CHECK_BOUNDS(SIZE(y),SIZE(x))

      ! Differentiate y(x) using centered finite differences

      n = SIZE(x)

      dy_dx(1) = (y(2) - y(1))/(x(2) - x(1))

      do i = 2,n-1
         dy_dx(i) = 0.5_WP*((y(i) - y(i-1))/(x(i) - x(i-1)) + &
                            (y(i+1) - y(i))/(x(i+1) - x(i)))
      end do

      dy_dx(n) = (y(n) - y(n-1))/(x(n) - x(n-1))

      ! Finish

      return

    end function deriv_

  end function lag_L_rd_

!****

  $if($MPI)

  subroutine bcast_ (md, root_rank, ml)

    class(mode_t), intent(inout)       :: md
    integer, intent(in)                :: root_rank
    class(model_t), intent(in), target :: ml

    ! Broadcast the mode_t

    if(MPI_RANK /= root_rank) then
       md%ml => ml
    endif

    call bcast(md%mp, root_rank)
    call bcast(md%op, root_rank)

    call bcast_alloc(md%x, root_rank)
    call bcast_alloc(md%y, root_rank)
    
    call bcast(md%x_ref, root_rank)
    call bcast(md%y_ref, root_rank)

    call bcast(md%omega, root_rank)

    call bcast(md%n_p, root_rank)
    call bcast(md%n_g, root_rank)
    call bcast(md%n_pg, root_rank)

    call bcast(md%n, root_rank)
    call bcast(md%n_iter, root_rank)

    call bcast(md%pruned, root_rank)

    ! Finish

    return

  end subroutine bcast_

  $endif

end module gyre_mode
