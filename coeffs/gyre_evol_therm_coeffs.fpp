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
  use core_hgroup
  use core_spline

  use gyre_therm_coeffs

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
    procedure :: get_${NAME}_1
    procedure :: get_${NAME}_v
  $endsub

  type, extends(therm_coeffs_t) :: evol_therm_coeffs_t
     $VAR_DECL(c_rad)
     $VAR_DECL(c_gen)
     $VAR_DECL(c_thm)
     $VAR_DECL(c_dif)
     $VAR_DECL(nabla)
     $VAR_DECL(kappa_ad)
     $VAR_DECL(kappa_S)
     $VAR_DECL(epsilon_ad)
     $VAR_DECL(epsilon_S)
   contains
     procedure :: init
     $PROC_DECL(c_rad)
     $PROC_DECL(dc_rad)
     $PROC_DECL(c_gen)
     $PROC_DECL(c_thm)
     $PROC_DECL(c_dif)
     $PROC_DECL(nabla)
     $PROC_DECL(kappa_ad)
     $PROC_DECL(kappa_S)
     $PROC_DECL(epsilon_ad)
     $PROC_DECL(epsilon_S)
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

  subroutine init (this, G, M_star, R_star, L_star, r, m, p, T, rho, &
                   Gamma_1, nabla_ad, delta, nabla, &
                   kappa, kappa_T, kappa_rho, &
                   epsilon, epsilon_T, epsilon_rho, deriv_type)

    class(evol_therm_coeffs_t), intent(out) :: this
    real(WP), intent(in)                    :: G
    real(WP), intent(in)                    :: M_star
    real(WP), intent(in)                    :: R_star
    real(WP), intent(in)                    :: L_star
    real(WP), intent(in)                    :: r(:)
    real(WP), intent(in)                    :: m(:)
    real(WP), intent(in)                    :: p(:)
    real(WP), intent(in)                    :: T(:)
    real(WP), intent(in)                    :: rho(:)
    real(WP), intent(in)                    :: Gamma_1(:)
    real(WP), intent(in)                    :: nabla_ad(:)
    real(WP), intent(in)                    :: delta(:)
    real(WP), intent(in)                    :: nabla(:)
    real(WP), intent(in)                    :: kappa(:)
    real(WP), intent(in)                    :: kappa_T(:)
    real(WP), intent(in)                    :: kappa_rho(:)
    real(WP), intent(in)                    :: epsilon(:)
    real(WP), intent(in)                    :: epsilon_T(:)
    real(WP), intent(in)                    :: epsilon_rho(:)
    character(LEN=*), intent(in)            :: deriv_type

    integer  :: n
    real(WP) :: V_x2(SIZE(r))
    real(WP) :: V(SIZE(r))
    real(WP) :: c_p(SIZE(r))
    real(WP) :: c_rad(SIZE(r))
    real(WP) :: c_gen(SIZE(r))
    real(WP) :: c_thm(SIZE(r))
    real(WP) :: c_dif(SIZE(r))
    real(WP) :: kappa_ad(SIZE(r))
    real(WP) :: kappa_S(SIZE(r))
    real(WP) :: epsilon_ad(SIZE(r))
    real(WP) :: epsilon_S(SIZE(r))
    real(WP) :: x(SIZE(r))

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

    ! Perform basic validations

    n = SIZE(r)

    $ASSERT(r(1) <= 0._WP,Invalid radius range)
    $ASSERT(r(n) >= R_star,Invalid radius range)

    $ASSERT(m(1) <= 0._WP,Invalid mass range)
    $ASSERT(m(n) >= M_star,Invalid mass range)

    ! Calculate coefficients

    where(r /= 0._WP)
       V_x2 = G*m*rho/(p*r*(r/R_star)**2)
    elsewhere
       V_x2 = 4._WP*PI*G*rho**2*R_star**2/(3._WP*p)
    end where

    x = r/R_star

    V = V_x2*x**2

    c_p = p*delta/(rho*T*nabla_ad)

    c_rad = 16._WP*PI*A_RADIATION*C_LIGHT*T**4*R_star*nabla*V_x2/(3._WP*kappa*rho*L_star)
    c_gen = 4._WP*PI*rho*epsilon*R_star**3/L_star
    c_thm = 4._WP*PI*rho*T*c_p*SQRT(G*M_star/R_star**3)*R_star**3/L_star

    kappa_ad = nabla_ad*kappa_T + kappa_rho/Gamma_1
    kappa_S = kappa_T - delta*kappa_rho

    epsilon_ad = nabla_ad*epsilon_T + epsilon_rho/Gamma_1
    epsilon_S = epsilon_T - delta*epsilon_rho

    c_dif = (kappa_ad-4._WP*nabla_ad)*V*nabla + nabla_ad*(dlny_dlnx(x, nabla_ad)+V)

    ! Initialize the therm_coeffs

    call this%sp_c_rad%init(x, c_rad, deriv_type, dy_dx_a=0._WP)
    call this%sp_c_gen%init(x, c_gen, deriv_type, dy_dx_a=0._WP)
    call this%sp_c_thm%init(x, c_thm, deriv_type, dy_dx_a=0._WP)
    call this%sp_c_dif%init(x, c_dif, deriv_type, dy_dx_a=0._WP)
    call this%sp_nabla%init(x, nabla, deriv_type, dy_dx_a=0._WP)
    call this%sp_kappa_S%init(x, kappa_S, deriv_type, dy_dx_a=0._WP)
    call this%sp_kappa_ad%init(x, kappa_ad, deriv_type, dy_dx_a=0._WP)
    call this%sp_epsilon_S%init(x, epsilon_S, deriv_type, dy_dx_a=0._WP)
    call this%sp_epsilon_ad%init(x, epsilon_ad, deriv_type, dy_dx_a=0._WP)

    ! Finish

    return

  contains

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

  $if($MPI)

  subroutine bcast_tc (tc, root_rank)

    class(evol_therm_coeffs_t), intent(inout) :: tc
    integer, intent(in)                       :: root_rank

    ! Broadcast the therm_coeffs

    call bcast(tc%sp_c_rad, root_rank)
    call bcast(tc%sp_c_gen, root_rank)
    call bcast(tc%sp_c_thm, root_rank)
    call bcast(tc%sp_c_dif, root_rank)
    call bcast(tc%sp_nabla, root_rank)
    call bcast(tc%sp_kappa_S, root_rank)
    call bcast(tc%sp_kappa_ad, root_rank)
    call bcast(tc%sp_epsilon_S, root_rank)
    call bcast(tc%sp_epsilon_ad, root_rank)

    ! Finish

    return

  end subroutine bcast_tc

  $endif

!****

  $define $PROC $sub

  $local $NAME $1

  function get_${NAME}_1 (this, x) result ($NAME)

    class(evol_therm_coeffs_t), intent(in) :: this
    real(WP), intent(in)                   :: x
    real(WP)                               :: $NAME

    ! Interpolate $NAME

    $NAME = this%sp_$NAME%interp(x)

    ! Finish

    return

  end function get_${NAME}_1

!****

  function get_${NAME}_v (this, x) result ($NAME)

    class(evol_therm_coeffs_t), intent(in) :: this
    real(WP), intent(in)                   :: x(:)
    real(WP)                               :: $NAME(SIZE(x))

    ! Interpolate $NAME

    $NAME = this%sp_$NAME%interp(x)

    ! Finish

    return

  end function get_${NAME}_v

  $endsub

  $PROC(c_rad)
  $PROC(c_gen)
  $PROC(c_thm)
  $PROC(c_dif)
  $PROC(nabla)
  $PROC(kappa_S)
  $PROC(kappa_ad)
  $PROC(epsilon_S)
  $PROC(epsilon_ad)

!****

  $define $DPROC $sub

  $local $NAME $1

  function get_d${NAME}_1 (this, x) result (d$NAME)

    class(evol_therm_coeffs_t), intent(in) :: this
    real(WP), intent(in)                   :: x
    real(WP)                               :: d$NAME

    ! Evaluate dln$NAME/dlnx

    if(x > 0._WP) then
       d$NAME = x*this%sp_$NAME%deriv(x)/this%sp_$NAME%interp(x)
    else
       d$NAME = 0._WP
    endif

    ! Finish

    return

  end function get_d${NAME}_1

!****

  function get_d${NAME}_v (this, x) result (d$NAME)

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

  end function get_d${NAME}_v

  $endsub

  $DPROC(c_rad)

end module gyre_evol_therm_coeffs
