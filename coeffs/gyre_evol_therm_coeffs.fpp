! Module   : gyre_evol_therm_coeffs
! Purpose  : thermal structure coefficients for evolutionary models
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

module gyre_evol_therm_coeffs

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel
  use core_spline

  use gyre_therm_coeffs
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

  type, extends(therm_coeffs_t) :: evol_therm_coeffs_t
     type(cocache_t) :: cc
     $VAR_DECL(c_rad)
     $VAR_DECL(c_thm)
     $VAR_DECL(c_dif)
     $VAR_DECL(c_eps_ad)
     $VAR_DECL(c_eps_S)
     $VAR_DECL(nabla)
     $VAR_DECL(kappa_ad)
     $VAR_DECL(kappa_S)
     $VAR_DECL(epsilon_ad)
     $VAR_DECL(epsilon_S)
     $VAR_DECL(tau_thm)
     logical :: cc_enabled
   contains
     private
     procedure, public :: init
     $if($GFORTRAN_PR57922)
     procedure, public :: final
     $PROC_DECL(c_rad)
     $PROC_DECL(dc_rad)
     $PROC_DECL(c_thm)
     $PROC_DECL(c_dif)
     $PROC_DECL(c_eps_ad)
     $PROC_DECL(c_eps_S)
     $PROC_DECL(nabla)
     $PROC_DECL(kappa_ad)
     $PROC_DECL(kappa_S)
     $PROC_DECL(epsilon_ad)
     $PROC_DECL(epsilon_S)
     $PROC_DECL(tau_thm)
     procedure, public :: enable_cache
     procedure, public :: disable_cache
     procedure, public :: fill_cache
  end type evol_therm_coeffs_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_tc
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: evol_therm_coeffs_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains 

  recursive subroutine init (this, G, M_star, R_star, L_star, r, m, p, rho, T, &
                             Gamma_1, nabla_ad, delta, nabla, &
                             kappa, kappa_rho, kappa_T, &
                             epsilon, epsilon_rho, epsilon_T, deriv_type, add_center)

    class(evol_therm_coeffs_t), intent(out) :: this
    real(WP), intent(in)                    :: G
    real(WP), intent(in)                    :: M_star
    real(WP), intent(in)                    :: R_star
    real(WP), intent(in)                    :: L_star
    real(WP), intent(in)                    :: r(:)
    real(WP), intent(in)                    :: m(:)
    real(WP), intent(in)                    :: p(:)
    real(WP), intent(in)                    :: rho(:)
    real(WP), intent(in)                    :: T(:)
    real(WP), intent(in)                    :: Gamma_1(:)
    real(WP), intent(in)                    :: nabla_ad(:)
    real(WP), intent(in)                    :: delta(:)
    real(WP), intent(in)                    :: nabla(:)
    real(WP), intent(in)                    :: kappa(:)
    real(WP), intent(in)                    :: kappa_rho(:)
    real(WP), intent(in)                    :: kappa_T(:)
    real(WP), intent(in)                    :: epsilon(:)
    real(WP), intent(in)                    :: epsilon_rho(:)
    real(WP), intent(in)                    :: epsilon_T(:)
    character(LEN=*), intent(in)            :: deriv_type
    logical, intent(in), optional           :: add_center

    logical  :: add_center_
    integer  :: n
    real(WP) :: V_x2(SIZE(r))
    real(WP) :: V(SIZE(r))
    real(WP) :: c_p(SIZE(r))
    real(WP) :: c_rad(SIZE(r))
    real(WP) :: c_thm(SIZE(r))
    real(WP) :: c_dif(SIZE(r))
    real(WP) :: kappa_ad(SIZE(r))
    real(WP) :: kappa_S(SIZE(r))
    real(WP) :: epsilon_ad(SIZE(r))
    real(WP) :: epsilon_S(SIZE(r))
    real(WP) :: c_eps_ad(SIZE(r))
    real(WP) :: c_eps_S(SIZE(r))
    real(WP) :: x(SIZE(r))
    real(WP) :: dtau_thm(SIZE(r))
    real(WP) :: tau_thm(SIZE(r))
    integer  :: i

    $CHECK_BOUNDS(SIZE(m),SIZE(r))
    $CHECK_BOUNDS(SIZE(p),SIZE(r))
    $CHECK_BOUNDS(SIZE(rho),SIZE(r))
    $CHECK_BOUNDS(SIZE(T),SIZE(r))
    $CHECK_BOUNDS(SIZE(Gamma_1),SIZE(r))
    $CHECK_BOUNDS(SIZE(nabla_ad),SIZE(r))
    $CHECK_BOUNDS(SIZE(delta),SIZE(r))
    $CHECK_BOUNDS(SIZE(nabla),SIZE(r))
    $CHECK_BOUNDS(SIZE(kappa),SIZE(r))
    $CHECK_BOUNDS(SIZE(kappa_T),SIZE(r))
    $CHECK_BOUNDS(SIZE(kappa_rho),SIZE(r))
    $CHECK_BOUNDS(SIZE(epsilon),SIZE(r))
    $CHECK_BOUNDS(SIZE(epsilon_T),SIZE(r))
    $CHECK_BOUNDS(SIZE(epsilon_rho),SIZE(r))

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
                      y_centered(r, Gamma_1), y_centered(r, nabla_ad), y_centered(r, delta), y_centered(r, nabla), &
                      y_centered(r, kappa), y_centered(r, kappa_rho), y_centered(r, kappa_T), &
                      y_centered(r, epsilon), y_centered(r, epsilon_rho), y_centered(r, epsilon_T), &
                      deriv_type, .FALSE.)

    else

       ! Perform basic validations

       $ASSERT(r(1) == 0._WP,First grid point not at center)
       $ASSERT(m(1) == 0._WP,First grid point not at center)

       $ASSERT(ALL(r(2:) >= r(:SIZE(r)-1)),Non-monotonic radius data)
       $ASSERT(ALL(m(2:) >= m(:SIZE(m)-1)),Non-monotonic mass data)

       ! Calculate coefficients

       n = SIZE(r)

       where(r /= 0._WP)
          V_x2 = G*m*rho/(p*r*(r/R_star)**2)
       elsewhere
          V_x2 = 4._WP*PI*G*rho**2*R_star**2/(3._WP*p)
       end where

       x = r/R_star

       V = V_x2*x**2

       c_p = p*delta/(rho*T*nabla_ad)

       c_rad = 16._WP*PI*A_RADIATION*C_LIGHT*T**4*R_star*nabla*V_x2/(3._WP*kappa*rho*L_star)
       c_thm = 4._WP*PI*rho*T*c_p*SQRT(G*M_star/R_star**3)*R_star**3/L_star

       kappa_ad = nabla_ad*kappa_T + kappa_rho/Gamma_1
       kappa_S = kappa_T - delta*kappa_rho

       c_dif = (kappa_ad-4._WP*nabla_ad)*V*nabla + nabla_ad*(dlny_dlnx(x, nabla_ad)+V)

       epsilon_ad = nabla_ad*epsilon_T + epsilon_rho/Gamma_1
       epsilon_S = epsilon_T - delta*epsilon_rho

       c_eps_ad = 4._WP*PI*rho*epsilon_ad*R_star**3/L_star
       c_eps_S = 4._WP*PI*rho*epsilon_S*R_star**3/L_star

       dtau_thm = 4._WP*PI*rho*r**2*T*c_p*SQRT(G*M_star/R_star**3)/L_star

       tau_thm(n) = 0._WP

       do i = n-1,1,-1
          tau_thm(i) = tau_thm(i+1) + &
               0.5_WP*(dtau_thm(i+1) + dtau_thm(i))*(r(i+1) - r(i))
       end do

       ! Initialize the therm_coeffs

       !$OMP PARALLEL SECTIONS
       !$OMP SECTION
       call this%sp_c_rad%init(x, c_rad, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_c_thm%init(x, c_thm, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_c_dif%init(x, c_dif, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_c_eps_ad%init(x, c_eps_ad, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_c_eps_S%init(x, c_eps_S, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_nabla%init(x, nabla, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_kappa_S%init(x, kappa_S, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_kappa_ad%init(x, kappa_ad, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_epsilon_S%init(x, epsilon_S, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_epsilon_ad%init(x, epsilon_ad, deriv_type, dy_dx_a=0._WP)
       !$OMP SECTION
       call this%sp_tau_thm%init(x, tau_thm, deriv_type, dy_dx_a=0._WP)
       !$OMP END PARALLEL SECTIONS

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

    function dlny_dlnx (x, y)

      real(WP), intent(in) :: x(:)
      real(WP), intent(in) :: y(:)
      real(WP)             :: dlny_dlnx(SIZE(x))

      integer :: n
      integer :: i

      $CHECK_BOUNDS(SIZE(y),SIZE(x))

      ! Calculate the logarithmic derivative of y wrt x

      n = SIZE(x)

      dlny_dlnx(1) = 0._WP

      do i = 2,n-1
         dlny_dlnx(i) = x(i)*0.5_WP*((y(i)-y(i-1))/(x(i)-x(i-1)) + (y(i+1)-y(i))/(x(i+1)-x(i)))/y(i)
      end do

      dlny_dlnx(n) = x(n)*(y(n)-y(n-1))/(x(n)-x(n-1))/y(n)

      ! Finish

    end function dlny_dlnx

  end subroutine init

!****

  $if($GFORTRAN_PR57922)
  
  subroutine final (this)

    class(evol_therm_coeffs_t), intent(inout) :: this

    ! Finalize the therm_coeffs

    call this%sp_c_rad%final()
    call this%sp_c_thm%final()
    call this%sp_c_dif%final()
    call this%sp_c_eps_ad%final()
    call this%sp_c_eps_S%final()
    call this%sp_nabla%final()
    call this%sp_kappa_ad%final()
    call this%sp_kappa_S%final()
    call this%sp_epsilon_ad%final()
    call this%sp_epsilon_S%final()
    call this%sp_tau_thm%final()

    if(this%cc_enabled) call this%cc%final()

    ! Finish

    return

  end subroutine final

  $endif

!****

  $if($MPI)

  subroutine bcast_tc (tc, root_rank)

    class(evol_therm_coeffs_t), intent(inout) :: tc
    integer, intent(in)                       :: root_rank

    ! Broadcast the therm_coeffs

    call bcast(tc%cc, root_rank)

    call bcast(tc%sp_c_rad, root_rank)
    call bcast(tc%sp_c_thm, root_rank)
    call bcast(tc%sp_c_dif, root_rank)
    call bcast(tc%sp_c_eps_ad, root_rank)
    call bcast(tc%sp_c_eps_S, root_rank)
    call bcast(tc%sp_nabla, root_rank)
    call bcast(tc%sp_kappa_S, root_rank)
    call bcast(tc%sp_kappa_ad, root_rank)
    call bcast(tc%sp_epsilon_S, root_rank)
    call bcast(tc%sp_epsilon_ad, root_rank)
    call bcast(tc%sp_tau_thm, root_rank)

    call bcast(tc%cc_enabled, root_rank)

    ! Finish

    return

  end subroutine bcast_tc

  $endif

!****

  $define $PROC $sub

  $local $NAME $1
  $local $I_CC $2

  function ${NAME}_1 (this, x) result ($NAME)

    class(evol_therm_coeffs_t), intent(in) :: this
    real(WP), intent(in)                   :: x
    real(WP)                               :: $NAME

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

    class(evol_therm_coeffs_t), intent(in) :: this
    real(WP), intent(in)                   :: x(:)
    real(WP)                               :: $NAME(SIZE(x))

    ! Interpolate $NAME

    $NAME = this%sp_$NAME%interp(x)

    ! Finish

    return

  end function ${NAME}_v

  $endsub

  $PROC(c_rad,1)
  $PROC(c_thm,2)
  $PROC(c_dif,3)
  $PROC(c_eps_ad,4)
  $PROC(c_eps_S,5)
  $PROC(nabla,6)
  $PROC(kappa_S,7)
  $PROC(kappa_ad,8)
  $PROC(epsilon_S,9)
  $PROC(epsilon_ad,10)
  $PROC(tau_thm,11)

!****

  $define $DPROC $sub

  $local $NAME $1
  $local $I_CC $2

  function d${NAME}_1 (this, x) result (d$NAME)

    class(evol_therm_coeffs_t), intent(in) :: this
    real(WP), intent(in)                   :: x
    real(WP)                               :: d$NAME

    ! Evaluate dln$NAME/dlnx

    if(this%cc_enabled) then
       d$NAME = this%cc%lookup($I_CC, x)
    else
       if(x > 0._WP) then
          d$NAME = x*this%sp_$NAME%deriv(x)/this%sp_$NAME%interp(x)
       else
          d$NAME = 0._WP
       endif
    endif

    ! Finish

    return

  end function d${NAME}_1

!****

  function d${NAME}_v (this, x) result (d$NAME)

    class(evol_therm_coeffs_t), intent(in) :: this
    real(WP), intent(in)                   :: x(:)
    real(WP)                               :: d$NAME(SIZE(x))

    ! Interpolate dln$NAME/dlnx

    where(x > 0._WP)
       d$NAME = x*this%sp_$NAME%deriv(x)/this%sp_$NAME%interp(x)
    elsewhere
       d$NAME = 0._WP
    endwhere

    ! Finish

    return

  end function d${NAME}_v

  $endsub

  $DPROC(c_rad,12)

!****

  subroutine enable_cache (this)

    class(evol_therm_coeffs_t), intent(inout) :: this

    ! Enable the coefficient cache

    this%cc_enabled = .TRUE.

    ! Finish

    return

  end subroutine enable_cache

!****

  subroutine disable_cache (this)

    class(evol_therm_coeffs_t), intent(inout) :: this

    ! Disable the coefficient cache

    this%cc_enabled = .FALSE.

    ! Finish

    return

  end subroutine disable_cache

!****

  subroutine fill_cache (this, x)

    class(evol_therm_coeffs_t), intent(inout) :: this
    real(WP), intent(in)                      :: x(:)

    real(WP) :: c(12,SIZE(x))

    ! Fill the coefficient cache

    !$OMP PARALLEL SECTIONS
    !$OMP SECTION
    c(1,:) = this%c_rad(x)
    !$OMP SECTION
    c(2,:) = this%c_thm(x)
    !$OMP SECTION
    c(3,:) = this%c_dif(x)
    !$OMP SECTION
    c(4,:) = this%c_eps_ad(x)
    !$OMP SECTION
    c(5,:) = this%c_eps_S(x)
    !$OMP SECTION
    c(6,:) = this%nabla(x)
    !$OMP SECTION
    c(7,:) = this%kappa_S(x)
    !$OMP SECTION
    c(8,:) = this%kappa_ad(x)
    !$OMP SECTION
    c(9,:) = this%epsilon_S(x)
    !$OMP SECTION
    c(10,:) = this%epsilon_ad(x)
    !$OMP SECTION
    c(11,:) = this%tau_thm(x)
    !$OMP SECTION
    c(12,:) = this%dc_rad(x)
    !$OMP END PARALLEL SECTIONS

    call this%cc%init(x, c)

    ! Finish

    return
  end subroutine fill_cache

end module gyre_evol_therm_coeffs
