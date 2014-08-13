! Module   : gyre_cons_model
! Purpose  : self-consistent model
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

module gyre_scons_model

  ! Uses

  use core_kinds
  use core_parallel
  use core_order

  use gyre_constants
  use gyre_model
  use gyre_cocache

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure :: ${NAME}_1_
    procedure :: ${NAME}_v_
  $endsub

  type, extends(model_t) :: scons_model_t
     private
     real(WP), allocatable :: dt_x(:)
     real(WP), allocatable :: dt_V(:)
     real(WP), allocatable :: dt_U(:)
     real(WP), allocatable :: dt_c_1(:)
     real(WP), allocatable :: dt_Gamma_1(:)
     real(WP)              :: dt_Omega_rot
     real(WP), public      :: M_star
     real(WP), public      :: R_star
     real(WP), public      :: L_star
     real(WP)              :: p_c
     real(WP)              :: rho_c
   contains
     private
     $PROC_DECL(V)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(c_1)
     $PROC_DECL(Gamma_1)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(delta)
     $PROC_DECL(c_rad)
     $PROC_DECL(dc_rad)
     $PROC_DECL(c_thm)
     $PROC_DECL(c_dif)
     $PROC_DECL(c_eps_ad)
     $PROC_DECL(c_eps_S)
     $PROC_DECL(nabla)
     $PROC_DECL(kappa_ad)
     $PROC_DECL(kappa_S)
     $PROC_DECL(tau_thm)
     $PROC_DECL(Omega_rot)
     procedure, public :: pi_c => pi_c_
     procedure, public :: is_zero => is_zero_
     procedure, public :: attach_cache => attach_cache_
     procedure, public :: detach_cache => detach_cache_
     procedure, public :: fill_cache => fill_cache_
  end type scons_model_t

  ! Interfaces

  interface scons_model_t
     module procedure scons_model_t_mech_coeffs_
     module procedure scons_model_t_mech_
  end interface scons_model_t

  $if ($MPI)
  interface bcast
     module procedure bcast_
  end interface bcast
  $endif

  ! Access specifiers

  private

  public :: scons_model_t
  $if ($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  function scons_model_t_mech_coeffs_ (x, V, U, c_1, Gamma_1, Omega_rot) result (ml)

    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: V(:)
    real(WP), intent(in) :: U(:)
    real(WP), intent(in) :: c_1(:)
    real(WP), intent(in) :: Gamma_1(:)
    real(WP), intent(in) :: Omega_rot
    type(scons_model_t)    :: ml

    integer  :: n
    real(WP) :: B
    integer  :: i

    $CHECK_BOUNDS(SIZE(V),SIZE(x))
    $CHECK_BOUNDS(SIZE(U),SIZE(x))
    $CHECK_BOUNDS(SIZE(Gamma_1),SIZE(x))

    ! Construct the scons_model_t

    n = SIZE(x)

    ! Set up U

    allocate(ml%dt_U(n))

    ml%dt_U(1:2) = 3._WP
    ml%dt_U(3:) = U(3:)

    ! Set up c_1 (integration from outside in)

    allocate(ml%dt_c_1(n))

    ml%dt_c_1(n) = c_1(n)

    do i = n-1, 2, -1
       B = LOG(U(i)/U(i+1))/LOG(x(i)/x(i+1))
       ml%dt_c_1(i) = ml%dt_c_1(i+1)*EXP(-(U(i)-U(i+1))/B)*(x(i)/x(i+1))**3
    end do

    ml%dt_c_1(1) = ml%dt_c_1(2)

    ! Set up V and Gamma for linear interpolation

    ml%dt_V = V
    ml%dt_Gamma_1 = Gamma_1

    ! Store the abscissa

    ml%dt_x = x

    ! Finish

    return

  end function scons_model_t_mech_coeffs_

!****

  recursive function scons_model_t_mech_ (M_star, R_star, L_star, r, m, p, rho, T, N2, &
                                          Gamma_1, nabla_ad, delta, Omega_rot, &
                                          add_center) result (ml)

    real(WP), intent(in)          :: M_star
    real(WP), intent(in)          :: R_star
    real(WP), intent(in)          :: L_star
    real(WP), intent(in)          :: r(:)
    real(WP), intent(in)          :: m(:)
    real(WP), intent(in)          :: p(:)
    real(WP), intent(in)          :: rho(:)
    real(WP), intent(in)          :: T(:)
    real(WP), intent(in)          :: N2(:)
    real(WP), intent(in)          :: Gamma_1(:)
    real(WP), intent(in)          :: nabla_ad(:)
    real(WP), intent(in)          :: delta(:)
    real(WP), intent(in)          :: Omega_rot(:)
    logical, optional, intent(in) :: add_center
    type(scons_model_t)           :: ml

    logical  :: add_center_
    integer  :: n
    real(WP) :: V(SIZE(r))
    real(WP) :: U(SIZE(r))
    real(WP) :: c_1(SIZE(r))
    real(WP) :: x(SIZE(r))
    real(WP) :: Omega_rot_

    $CHECK_BOUNDS(SIZE(m),SIZE(r))
    $CHECK_BOUNDS(SIZE(p),SIZE(r))
    $CHECK_BOUNDS(SIZE(rho),SIZE(r))
    $CHECK_BOUNDS(SIZE(T),SIZE(r))
    $CHECK_BOUNDS(SIZE(N2),SIZE(r))
    $CHECK_BOUNDS(SIZE(Gamma_1),SIZE(r))
    $CHECK_BOUNDS(SIZE(nabla_ad),SIZE(r))
    $CHECK_BOUNDS(SIZE(delta),SIZE(r))
    $CHECK_BOUNDS(SIZE(Omega_rot),SIZE(r))

    if(PRESENT(add_center)) then
       add_center_ = add_center
    else
       add_center_ = .FALSE.
    endif

    ! Construct the evol_model_t using the mechanical structure data

    ! See if we need a central point

    if (add_center_) then

       ! Add a central point and initialize using recursion

       ml = scons_model_t(M_star, R_star, L_star, [0._WP,r], [0._WP,m], &
                          prep_center_(r, p), prep_center_(r, rho), prep_center_(r, T), &
                          [0._WP,N2], prep_center_(r, Gamma_1), prep_center_(r, nabla_ad), prep_center_(r, delta), &
                          prep_center_(r, Omega_rot), .FALSE.)

    else
       
       ! Perform basic validations
       
       n = SIZE(r)

       $ASSERT(r(1) == 0._WP,First grid point not at center)
       $ASSERT(m(1) == 0._WP,First grid point not at center)

       $ASSERT(ALL(r(2:) > r(:n-1)),Non-monotonic radius data)
       $ASSERT(ALL(m(2:) >= m(:n-1)),Non-monotonic mass data)

       ! Calculate coefficients

       where(r /= 0._WP)
          V = G_GRAVITY*m*rho/(p*r)
          U = 4._WP*PI*rho*r**3/m
          c_1 = (r/R_star)**3/(m/M_star)
       elsewhere
          V = 0._WP
          U = 3._WP
          c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
       end where

       x = r/R_star

       Omega_rot_ = SQRT(R_star**3/(G_GRAVITY*M_star))*Omega_rot(n)

       ! Initialize the model

       ml = scons_model_t(x, V, U, c_1, Gamma_1, Omega_rot_)

       ml%M_star = M_star
       ml%R_star = R_star
       ml%L_star = L_star

       ml%p_c = p(1)
       ml%rho_c = rho(1)

    endif

    ! Finish

    return

  end function scons_model_t_mech_

!****
  
  function prep_center_ (x, y) result (y_prep)
      
    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: y(:)
    real(WP)             :: y_prep(SIZE(y)+1)

    real(WP) :: y_0

    $CHECK_BOUNDS(SIZE(x),SIZE(y))

    $ASSERT(SIZE(y) >= 2,Insufficient grid points)

    ! Use parabola fitting to interpolate y at the center
      
    y_0 = (x(2)**2*y(1) - x(1)**2*y(2))/(x(2)**2 - x(1)**2)

    ! Preprend this to the array

    y_prep = [y_0,y]

    ! Finish

    return

  end function prep_center_

!****

  function V_1_ (this, x) result (V)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x
    real(WP)                       :: V

    integer  :: i
    real(WP) :: w

    ! Calculate V

    call locate(this%dt_x, x, i)

    $ASSERT(i > 0 .AND. i < SIZE(this%dt_x),Out-of-bounds interpolation)

    w = (x - this%dt_x(i))/(this%dt_x(i+1) - this%dt_x(i))

    V = (1._WP-w)*this%dt_V(i) + w*this%dt_V(i+1)

    ! Finish

    return

  end function V_1_

!****
  
  function V_v_ (this, x) result (V)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x(:)
    real(WP)                       :: V(SIZE(x))

    integer :: i

    ! Calculate V

    x_loop : do i = 1,SIZE(x)
       V(i) = this%V(x(i))
    end do x_loop

    ! Finish

    return

  end function V_v_

!****

  function As_1_ (this, x) result (As)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x
    real(WP)                       :: As

    integer  :: i
    real(WP) :: B

    ! Calculate As

    call locate(this%dt_x, x, i)

    $ASSERT(i > 0 .AND. i < SIZE(this%dt_x),Out-of-bounds interpolation)

    if (i > 1) then
       B = LOG(this%dt_U(i)/this%dt_U(i+1))/LOG(this%dt_x(i)/this%dt_x(i+1))
    else
       B = 0._WP
    endif

    As = -this%V(x)/this%Gamma_1(x) - this%U(x) - B + 3._WP

    ! Finish

    return

  end function As_1_

!****
  
  function As_v_ (this, x) result (As)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x(:)
    real(WP)                       :: As(SIZE(x))

    integer :: i

    ! Calculate As

    x_loop : do i = 1,SIZE(x)
       As(i) = this%As(x(i))
    end do x_loop

    ! Finish

    return

  end function As_v_

!****

  function U_1_ (this, x) result (U)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x
    real(WP)                       :: U

    integer  :: i
    real(WP) :: B

    ! Calculate U

    call locate(this%dt_x, x, i)

    $ASSERT(i > 0 .AND. i < SIZE(this%dt_x),Out-of-bounds interpolation)

    if (i > 1) then
       B = LOG(this%dt_U(i)/this%dt_U(i+1))/LOG(this%dt_x(i)/this%dt_x(i+1))
       U = this%dt_U(i)*(x/this%dt_x(i))**B
    else
       U = 3._WP
    endif

    ! Finish

    return

  end function U_1_

!****
  
  function U_v_ (this, x) result (U)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x(:)
    real(WP)                       :: U(SIZE(x))

    integer :: i

    ! Calculate U

    x_loop : do i = 1,SIZE(x)
       U(i) = this%U(x(i))
    end do x_loop

    ! Finish

    return

  end function U_v_

!****

  function c_1_1_ (this, x) result (c_1)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x
    real(WP)                       :: c_1

    integer  :: i
    real(WP) :: B

    ! Calculate c_1

    call locate(this%dt_x, x, i)

    $ASSERT(i > 0 .AND. i < SIZE(this%dt_x),Out-of-bounds interpolation)

    if (i > 1) then
       B = LOG(this%dt_U(i)/this%dt_U(i+1))/LOG(this%dt_x(i)/this%dt_x(i+1))
       c_1 = this%dt_c_1(i)*EXP(-(this%U(x)-this%dt_U(i))/B)*(x/this%dt_x(i))**3
    else
       c_1 = this%dt_c_1(1)
    endif

    ! Finish

    return

  end function c_1_1_

!****
  
  function c_1_v_ (this, x) result (c_1)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x(:)
    real(WP)                       :: c_1(SIZE(x))

    integer :: i

    ! Calculate c_1

    x_loop : do i = 1,SIZE(x)
       c_1(i) = this%c_1(x(i))
    end do x_loop

    ! Finish

    return

  end function c_1_v_

!****

  function Gamma_1_1_ (this, x) result (Gamma_1)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x
    real(WP)                       :: Gamma_1

    integer  :: i
    real(WP) :: w

    ! Calculate Gamma_1

    call locate(this%dt_x, x, i)

    $ASSERT(i > 0 .AND. i < SIZE(this%dt_x),Out-of-bounds interpolation)

    w = (x - this%dt_x(i))/(this%dt_x(i+1) - this%dt_x(i))

    Gamma_1 = (1._WP-w)*this%dt_Gamma_1(i) + w*this%dt_Gamma_1(i+1)

    ! Finish

    return

  end function Gamma_1_1_

!****
  
  function Gamma_1_v_ (this, x) result (Gamma_1)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x(:)
    real(WP)                       :: Gamma_1(SIZE(x))

    integer :: i

    ! Calculate Gamma_1
    
    x_loop : do i = 1,SIZE(x)
       Gamma_1(i) = this%Gamma_1(x(i))
    end do x_loop

    ! Finish

    return

  end function Gamma_1_v_

!****

  $define $PROC $sub

  $local $NAME $1

  function ${NAME}_1_ (this, x) result ($NAME)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x
    real(WP)                       :: $NAME

    ! Abort with $NAME undefined

    $ABORT($NAME is undefined)

    ! Finish

    return

  end function ${NAME}_1_

!****

  function ${NAME}_v_ (this, x) result ($NAME)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x(:)
    real(WP)                       :: $NAME(SIZE(x))

    ! Abort with $NAME undefined

    $ABORT($NAME is undefined)

    ! Finish

    return

  end function ${NAME}_v_

  $endsub

  $PROC(nabla_ad)
  $PROC(delta)
  $PROC(c_rad)
  $PROC(dc_rad)
  $PROC(c_thm)
  $PROC(c_dif)
  $PROC(c_eps_ad)
  $PROC(c_eps_S)
  $PROC(nabla)
  $PROC(kappa_S)
  $PROC(kappa_ad)
  $PROC(tau_thm)

!****

  function Omega_rot_1_ (this, x) result (Omega_rot)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x
    real(WP)                       :: Omega_rot

    ! Calculate Omega_rot

    Omega_rot = this%dt_Omega_rot

    ! Finish

    return

  end function Omega_rot_1_

!****
  
  function Omega_rot_v_ (this, x) result (Omega_rot)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x(:)
    real(WP)                       :: Omega_rot(SIZE(x))

    integer :: i

    ! Calculate Omega_rot
    
    x_loop : do i = 1,SIZE(x)
       Omega_rot(i) = this%Omega_rot(x(i))
    end do x_loop

    ! Finish

    return

  end function Omega_rot_v_

!****

  function pi_c_ (this) result (pi_c)

    class(scons_model_t), intent(in) :: this
    real(WP)                       :: pi_c

    ! Calculate pi_c = V/x^2 as x -> 0

    pi_c = 2._WP

    ! Finish

    return

  end function pi_c_

!****

  function is_zero_ (this, x) result (is_zero)

    class(scons_model_t), intent(in) :: this
    real(WP), intent(in)           :: x
    logical                        :: is_zero

    ! Determine whether the point at x has a vanishing pressure and/or density

    is_zero = x >= 1._WP

    ! Finish

    return

  end function is_zero_

!****

  subroutine attach_cache_ (this, cc)

    class(scons_model_t), intent(inout)     :: this
    class(cocache_t), pointer, intent(in) :: cc

    ! Attach a coefficient cache (no-op, since we don't cache)

    ! Finish

    return

  end subroutine attach_cache_

!****

  subroutine detach_cache_ (this)

    class(scons_model_t), intent(inout) :: this

    ! Detach the coefficient cache (no-op, since we don't cache)

    ! Finish

    return

  end subroutine detach_cache_

!****

  subroutine fill_cache_ (this, x)

    class(scons_model_t), intent(inout) :: this
    real(WP), intent(in)              :: x(:)

    ! Fill the coefficient cache (no-op, since we don't cache)

    ! Finish

    return

  end subroutine fill_cache_

!****

  $if ($MPI)

  subroutine bcast_ (ml, root_rank)

    class(scons_model_t), intent(inout) :: ml
    integer, intent(in)               :: root_rank

    ! Broadcast the scons_model_t

    call bcast(ml%dt_Gamma_1, root_rank)

    ! Finish

    return

  end subroutine bcast_

  $endif

end module gyre_scons_model
