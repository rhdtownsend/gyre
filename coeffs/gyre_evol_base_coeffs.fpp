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
     real(WP), public :: M_star
     real(WP), public :: R_star
     real(WP), public :: L_star
     real(WP)         :: p_c
     real(WP)         :: rho_c
     real(WP)         :: G
   contains
     private
     procedure, public :: init
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
     procedure, public :: pi_c
     procedure, public :: conv_freq
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

  subroutine init (this, G, M_star, R_star, L_star, r, m, p, rho, T, N2, Gamma_1, nabla_ad, delta, deriv_type)

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
    character(LEN=*), intent(in)           :: deriv_type

    integer  :: n
    real(WP) :: V(SIZE(r))
    real(WP) :: As(SIZE(r))
    real(WP) :: U(SIZE(r))
    real(WP) :: c_1(SIZE(r))
    real(WP) :: x(SIZE(r))

    $CHECK_BOUNDS(SIZE(m),SIZE(r))
    $CHECK_BOUNDS(SIZE(p),SIZE(r))
    $CHECK_BOUNDS(SIZE(rho),SIZE(r))
    $CHECK_BOUNDS(SIZE(T),SIZE(r))
    $CHECK_BOUNDS(SIZE(N2),SIZE(r))
    $CHECK_BOUNDS(SIZE(Gamma_1),SIZE(r))
    $CHECK_BOUNDS(SIZE(nabla_ad),SIZE(r))
    $CHECK_BOUNDS(SIZE(delta),SIZE(r))

    ! Perform basic validations

    n = SIZE(r)

    $ASSERT(r(1) <= 0._WP,Invalid radius range)
    $ASSERT(r(n) >= R_star,Invalid radius range)

    $ASSERT(m(1) <= 0._WP,Invalid mass range)
    $ASSERT(m(n) >= M_star,Invalid mass range)

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

    x = r/R_star

    ! Initialize the base_coeffs

    call this%sp_m%init(x, m, deriv_type, dy_dx_a=0._WP)
    call this%sp_p%init(x, p, deriv_type, dy_dx_a=0._WP)
    call this%sp_rho%init(x, rho, deriv_type, dy_dx_a=0._WP)
    call this%sp_T%init(x, T, deriv_type, dy_dx_a=0._WP)

    call this%sp_V%init(x, V, deriv_type, dy_dx_a=0._WP)
    call this%sp_As%init(x, As, deriv_type, dy_dx_a=0._WP)
    call this%sp_U%init(x, U, deriv_type, dy_dx_a=0._WP)
    call this%sp_c_1%init(x, c_1, deriv_type, dy_dx_a=0._WP)
    call this%sp_Gamma_1%init(x, Gamma_1, deriv_type, dy_dx_a=0._WP)
    call this%sp_nabla_ad%init(x, nabla_ad, deriv_type, dy_dx_a=0._WP)
    call this%sp_delta%init(x, delta, deriv_type, dy_dx_a=0._WP)

    this%M_star = M_star
    this%R_star = R_star
    this%L_star = L_star

    this%p_c = p(1)
    this%rho_c = rho(1)

    this%G = G

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_bc (bc, root_rank)

    class(evol_base_coeffs_t), intent(inout) :: bc
    integer, intent(in)                      :: root_rank

    ! Broadcast the base_coeffs

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

    call bcast(bc%M_star, root_rank)
    call bcast(bc%R_star, root_rank)
    call bcast(bc%L_star, root_rank)

    call bcast(bc%p_c, root_rank)
    call bcast(bc%rho_c, root_rank)

    call bcast(bc%G, root_rank)

    ! Finish

    return

  end subroutine bcast_bc

  $endif

!****

  $define $PROC $sub

  $local $NAME $1

  function ${NAME}_1 (this, x) result ($NAME)

    class(evol_base_coeffs_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    real(WP)                              :: $NAME

    ! Interpolate $NAME

    $NAME = this%sp_$NAME%interp(x)

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

  $PROC(m)
  $PROC(p)
  $PROC(rho)
  $PROC(T)
  $PROC(V)
  $PROC(As)
  $PROC(U)
  $PROC(c_1)
  $PROC(Gamma_1)
  $PROC(nabla_ad)
  $PROC(delta)

!****

  function pi_c (this)

    class(evol_base_coeffs_t), intent(in) :: this
    real(WP)                              :: pi_c

    ! Calculate pi_c = V/x^2 as x -> 0

    pi_c = 4._WP*PI*this%G*this%rho_c**2*this%R_star**2/(3._WP*this%p_c)

    ! Finish

    return

  end function pi_c

!****

  function conv_freq (this, freq, from_units, to_units)

    class(evol_base_coeffs_t), intent(in) :: this
    complex(WP), intent(in)               :: freq
    character(LEN=*), intent(in)          :: from_units
    character(LEN=*), intent(in)          :: to_units
    complex(WP)                           :: conv_freq

    ! Convert the frequency

    conv_freq = freq/freq_scale(from_units)*freq_scale(to_units)

    ! Finish

    return

  contains

    function freq_scale (units)

      character(LEN=*), intent(in) :: units
      real(WP)                     :: freq_scale

      ! Calculate the scale factor to convert a dimensionless angular
      ! frequency to a dimensioned frequency

      select case(units)
      case('NONE')
         freq_scale = 1._WP
      case('HZ')
         freq_scale = 1._WP/(TWOPI*SQRT(this%R_star**3/(this%G*this%M_star)))
      case('UHZ')
         freq_scale = 1.E6_WP/(TWOPI*SQRT(this%R_star**3/(this%G*this%M_star)))
      case default
         $ABORT(Invalid units)
      end select

      ! Finish

      return

    end function freq_scale

  end function conv_freq

end module gyre_evol_base_coeffs
