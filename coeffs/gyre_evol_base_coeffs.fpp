! Module   : gyre_evol_base_coeffs
! Purpose  : base structure coefficients for evolutionary models
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

module gyre_evol_base_coeffs

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel
  use core_spline

  use gyre_base_coeffs
  use gyre_cocache

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $VAR_DECL $sub
    $local $NAME $1
    type(spline_t) :: sp_$NAME
  $endsub
  
  $define $PROC_DECL $sub
    $local $NAME $1
    procedure :: ${NAME}_1
    procedure :: ${NAME}_v
  $endsub

  $define $PROC_DECL_GEN $sub
    $local $NAME $1
    procedure       :: ${NAME}_1
    procedure       :: ${NAME}_v
    generic, public :: ${NAME} => ${NAME}_1, ${NAME}_v
  $endsub

  type, extends(base_coeffs_t) :: evol_base_coeffs_t
     private
     type(cocache_t) :: cc
     $VAR_DECL(m)
     $VAR_DECL(p)
     $VAR_DECL(rho)
     $VAR_DECL(T)
     $VAR_DECL(V)
     $VAR_DECL(As)
     $VAR_DECL(U)
     $VAR_DECL(c_1)
     $VAR_DECL(Gamma_1)
     $VAR_DECL(nabla_ad)
     $VAR_DECL(delta)
     $VAR_DECL(Omega_rot)
     real(WP), public :: M_star
     real(WP), public :: R_star
     real(WP), public :: L_star
     real(WP), public :: G
     real(WP)         :: p_c
     real(WP)         :: rho_c
     logical          :: cc_enabled
   contains
     private
     procedure, public :: init
     $if($GFORTRAN_PR57922)
     procedure, public :: final
     $endif
     $PROC_DECL_GEN(m)
     $PROC_DECL_GEN(p)
     $PROC_DECL_GEN(rho)
     $PROC_DECL_GEN(T)
     $PROC_DECL(V)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(c_1)
     $PROC_DECL(Gamma_1)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(delta)
     $PROC_DECL(Omega_rot)
     procedure, public :: pi_c
     procedure, public :: enable_cache
     procedure, public :: disable_cache
     procedure, public :: fill_cache
  end type evol_base_coeffs_t
 
  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_bc
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: evol_base_coeffs_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  recursive subroutine init (this, G, M_star, R_star, L_star, r, m, p, rho, T, N2, &
                             Gamma_1, nabla_ad, delta, Omega_rot, deriv_type, add_center)

    class(evol_base_coeffs_t), intent(out) :: this
    real(WP), intent(in)                   :: G
    real(WP), intent(in)                   :: M_star
    real(WP), intent(in)                   :: R_star
    real(WP), intent(in)                   :: L_star
    real(WP), intent(in)                   :: r(:)
    real(WP), intent(in)                   :: m(:)
    real(WP), intent(in)                   :: p(:)
    real(WP), intent(in)                   :: rho(:)
    real(WP), intent(in)                   :: T(:)
    real(WP), intent(in)                   :: N2(:)
    real(WP), intent(in)                   :: Gamma_1(:)
    real(WP), intent(in)                   :: nabla_ad(:)
    real(WP), intent(in)                   :: delta(:)
    real(WP), intent(in)                   :: Omega_rot(:)
    character(LEN=*), intent(in)           :: deriv_type
    logical, intent(in), optional          :: add_center

    logical  :: add_center_
    real(WP) :: V(SIZE(r))
    real(WP) :: As(SIZE(r))
    real(WP) :: U(SIZE(r))
    real(WP) :: c_1(SIZE(r))
    real(WP) :: Omega_rot_(SIZE(r))
    real(WP) :: x(SIZE(r))

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

    ! See if we need a central point

    if(add_center_) then

       ! Add a central point and initialize using recursion

       call this%init(G, M_star, R_star, L_star, [0._WP,r], [0._WP,m], &
                      y_centered(r, p), y_centered(r, rho), y_centered(r, T), &
                      [0._WP,N2], y_centered(r, Gamma_1), y_centered(r, nabla_ad), y_centered(r, delta), &
                      y_centered(r, Omega_rot), deriv_type, .FALSE.)

    else

       ! Perform basic validations
       
       $ASSERT(r(1) == 0._WP,First grid point not at center)
       $ASSERT(m(1) == 0._WP,First grid point not at center)

       $ASSERT(ALL(r(2:) >= r(:SIZE(r)-1)),Non-monotonic radius data)
       $ASSERT(ALL(m(2:) >= m(:SIZE(m)-1)),Non-monotonic mass data)

       ! Calculate coefficients

       where(r /= 0._WP)
          V = G*m*rho/(p*r)
          As = r**3*N2/(G*m)
          U = 4._WP*PI*rho*r**3/m
          c_1 = (r/R_star)**3/(m/M_star)
       elsewhere
          V = 0._WP
          As = 0._WP
          U = 3._WP
          c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
       end where

       Omega_rot_ = SQRT(R_star**3/(G*M_star))*Omega_rot

       x = r/R_star

       ! Initialize the base_coeffs

       !$OMP PARALLEL SECTIONS
       !$OMP SECTION
       call this%sp_m%init(x, m, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_p%init(x, p, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_rho%init(x, rho, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_T%init(x, T, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_V%init(x, V, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_As%init(x, As, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_U%init(x, U, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_c_1%init(x, c_1, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_Gamma_1%init(x, Gamma_1, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_nabla_ad%init(x, nabla_ad, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_delta%init(x, delta, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_Omega_rot%init(x, Omega_rot_, deriv_type, dy_dx_a=0._WP)
       !$OMP END PARALLEL SECTIONS

       this%M_star = M_star
       this%R_star = R_star
       this%L_star = L_star

       this%G = G

       this%cc_enabled = .FALSE.

    endif

    ! Finish

    return

  contains

    function y_centered (x, y)
      
      real(WP), intent(in) :: x(:)
      real(WP), intent(in) :: y(:)
      real(WP)             :: y_centered(SIZE(y)+1)

      real(WP) :: y_center

      $CHECK_BOUNDS(SIZE(x),SIZE(y))

      $ASSERT(SIZE(y) >= 2,Insufficient grid points)

      ! Use parabola fitting to interpolate y at the center
      
      y_center = (x(2)**2*y(1) - x(1)**2*y(2))/(x(2)**2 - x(1)**2)

      ! Create the centered array

      y_centered = [y_center,y]

      ! Finish

      return

    end function y_centered

  end subroutine init

!****

  $if($GFORTRAN_PR57922)

  subroutine final (this)

    class(evol_base_coeffs_t), intent(inout) :: this

    ! Finalize the base_coeffs

     call this%sp_m%final()
     call this%sp_p%final()
     call this%sp_rho%final()
     call this%sp_T%final()
     call this%sp_V%final()
     call this%sp_As%final()
     call this%sp_U%final()
     call this%sp_c_1%final()
     call this%sp_Gamma_1%final()
     call this%sp_nabla_ad%final()
     call this%sp_delta%final()
     call this%sp_Omega_rot%final()

     call this%cc%final()

     ! Finish

     return

   end subroutine final

   $endif

!****

  $if($MPI)

  subroutine bcast_bc (bc, root_rank)

    class(evol_base_coeffs_t), intent(inout) :: bc
    integer, intent(in)                      :: root_rank

    ! Broadcast the base_coeffs

    call bcast(bc%cc, root_rank)

    call bcast(bc%sp_m, root_rank)
    call bcast(bc%sp_p, root_rank)
    call bcast(bc%sp_rho, root_rank)
    call bcast(bc%sp_T, root_rank)

    call bcast(bc%sp_V, root_rank)
    call bcast(bc%sp_As, root_rank)
    call bcast(bc%sp_U, root_rank)
    call bcast(bc%sp_c_1, root_rank)
    call bcast(bc%sp_Gamma_1, root_rank)
    call bcast(bc%sp_nabla_ad, root_rank)
    call bcast(bc%sp_delta, root_rank)
    call bcast(bc%sp_Omega_rot, root_rank)

    call bcast(bc%M_star, root_rank)
    call bcast(bc%R_star, root_rank)
    call bcast(bc%L_star, root_rank)

    call bcast(bc%G, root_rank)

    call bcast(bc%p_c, root_rank)
    call bcast(bc%rho_c, root_rank)

    call bcast(bc%cc_enabled, root_rank)

    ! Finish

    return

  end subroutine bcast_bc

  $endif

!****

  $define $PROC $sub

  $local $NAME $1
  $local $I_CC $2

  function ${NAME}_1 (this, x) result ($NAME)

    class(evol_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    real(WP)                              :: $NAME

    ! Interpolate $NAME

    if(this%cc_enabled) then
       $NAME = this%cc%lookup($I_CC, x)
    else
       $NAME = this%sp_$NAME%interp(x)
    endif
       
    ! Finish

    return

  end function ${NAME}_1

!****

  function ${NAME}_v (this, x) result ($NAME)

    class(evol_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x(:)
    real(WP)                              :: $NAME(SIZE(x))

    ! Interpolate $NAME

    $NAME = this%sp_$NAME%interp(x)

    ! Finish

    return

  end function ${NAME}_v

  $endsub

  $PROC(m,1)
  $PROC(p,2)
  $PROC(rho,3)
  $PROC(T,4)
  $PROC(V,5)
  $PROC(As,6)
  $PROC(U,7)
  $PROC(c_1,8)
  $PROC(Gamma_1,9)
  $PROC(nabla_ad,10)
  $PROC(delta,11)
  $PROC(Omega_rot,12)

!****

  function pi_c (this)

    class(evol_base_coeffs_t), intent(in) :: this
    real(WP)                              :: pi_c

    ! Calculate pi_c = V/x^2 as x -> 0

    pi_c = 4._WP*PI*this%G*this%rho(0._WP)**2*this%R_star**2/(3._WP*this%p(0._WP))

    ! Finish

    return

  end function pi_c

!****

  subroutine enable_cache (this)

    class(evol_base_coeffs_t), intent(inout) :: this

    ! Enable the coefficient cache

    this%cc_enabled = .TRUE.

    ! Finish

    return

  end subroutine enable_cache

!****

  subroutine disable_cache (this)

    class(evol_base_coeffs_t), intent(inout) :: this

    ! Disable the coefficient cache

    this%cc_enabled = .FALSE.

    ! Finish

    return

  end subroutine disable_cache

!****

  subroutine fill_cache (this, x)

    class(evol_base_coeffs_t), intent(inout) :: this
    real(WP), intent(in)                     :: x(:)

    real(WP) :: c(12,SIZE(x))

    ! Fill the coefficient cache

    !$OMP PARALLEL SECTIONS
    !$OMP SECTION
    c(1,:) = this%m(x)
    !$OMP SECTION
    c(2,:) = this%p(x)
    !$OMP SECTION
    c(3,:) = this%rho(x)
    !$OMP SECTION
    c(4,:) = this%T(x)
    !$OMP SECTION
    c(5,:) = this%V(x)
    !$OMP SECTION
    c(6,:) = this%As(x)
    !$OMP SECTION
    c(7,:) = this%U(x)
    !$OMP SECTION
    c(8,:) = this%c_1(x)
    !$OMP SECTION
    c(9,:) = this%Gamma_1(x)
    !$OMP SECTION
    c(10,:) = this%nabla_ad(x)
    !$OMP SECTION
    c(11,:) = this%delta(x)
    !$OMP SECTION
    c(12,:) = this%Omega_rot(x)
    !$OMP END PARALLEL SECTIONS

    call this%cc%init(x, c)

    ! Finish

    return

  end subroutine fill_cache

end module gyre_evol_base_coeffs
