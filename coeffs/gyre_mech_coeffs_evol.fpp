! Module   : gyre_mech_coeffs_evol
! Purpose  : mechanical structure coefficients for evolutionary models

$include 'core.inc'

module gyre_mech_coeffs_evol

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel
  use core_spline

  use gyre_mech_coeffs

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

  type, extends(mech_coeffs_t) :: mech_coeffs_evol_t
     private
     $VAR_DECL(V)
     $VAR_DECL(As)
     $VAR_DECL(U)
     $VAR_DECL(c_1)
     $VAR_DECL(Gamma_1)
     real(WP) :: t_dyn
   contains
     private
     procedure, public :: init
     $if($MPI)
     procedure, public :: bcast => bcast_mc
     $endif
     $PROC_DECL(V)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(c_1)
     $PROC_DECL(Gamma_1)
     procedure, public :: conv_freq
  end type mech_coeffs_evol_t

  ! Access specifiers

  private

  public :: mech_coeffs_evol_t

  ! Procedures

contains

  subroutine init (this, G, R_star, M_star, r, m, p, rho, N2, Gamma_1)

    class(mech_coeffs_evol_t), intent(out) :: this
    real(WP), intent(in), optional         :: G
    real(WP), intent(in)                   :: R_star
    real(WP), intent(in)                   :: M_star
    real(WP), intent(in)                   :: r(:)
    real(WP), intent(in)                   :: m(:)
    real(WP), intent(in)                   :: p(:)
    real(WP), intent(in)                   :: rho(:)
    real(WP), intent(in)                   :: N2(:)
    real(WP), intent(in)                   :: Gamma_1(:)

    real(WP) :: G_
    integer  :: n
    real(WP) :: V(SIZE(r))
    real(WP) :: As(SIZE(r))
    real(WP) :: U(SIZE(r))
    real(WP) :: c_1(SIZE(r))
    real(WP) :: x(SIZE(r))

    $CHECK_BOUNDS(SIZE(m),SIZE(r))
    $CHECK_BOUNDS(SIZE(p),SIZE(r))
    $CHECK_BOUNDS(SIZE(rho),SIZE(r))
    $CHECK_BOUNDS(SIZE(N2),SIZE(r))
    $CHECK_BOUNDS(SIZE(Gamma_1),SIZE(r))

    if(PRESENT(G)) then
       G_ = G
    else
       G_ = G_GRAVITY
    endif

    ! Perform basic validations

    n = SIZE(r)

    $ASSERT(r(1) == 0._WP,Invalid radius range)
    $ASSERT(r(n) == R_star,Invalid radius range)

    $ASSERT(m(1) == 0._WP,Invalid mass range)
    $ASSERT(m(n) == M_star,Invalid mass range)

    ! Calculate coefficients

    where(r /= 0._WP)
       V = G_*m*rho/(p*r)
       As = r**3*N2/(G_*m)
       U = 4._WP*PI*rho*r**3/m
       c_1 = (r/R_star)**3/(m/M_star)
    elsewhere
       V = 0._WP
       As = 0._WP
       U = 3._WP
       c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
    end where

    x = r/R_star

    ! Initialize the mech_coeffs

    call this%sp_V%init(x, V, dy_dx_a=0._WP)
    call this%sp_As%init(x, As, dy_dx_a=0._WP, linear=.TRUE.)
    call this%sp_U%init(x, U, dy_dx_a=0._WP)
    call this%sp_c_1%init(x, c_1, dy_dx_a=0._WP)
    call this%sp_Gamma_1%init(x, Gamma_1, dy_dx_a=0._WP)

    this%t_dyn = SQRT(R_star**3/(G*M_star))

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_mc (this, root_rank)

    class(mech_coeffs_evol_t), intent(inout) :: this
    integer, intent(in)                      :: root_rank

    ! Broadcast the mech_coeffs

    call this%sp_V%bcast(root_rank)
    call this%sp_As%bcast(root_rank)
    call this%sp_U%bcast(root_rank)
    call this%sp_c_1%bcast(root_rank)
    call this%sp_Gamma_1%bcast(root_rank)

    call bcast(this%t_dyn, root_rank)

    ! Finish

    return

  end subroutine bcast_mc

  $endif

!****

  $define $PROC $sub

  $local $NAME $1

  function get_${NAME}_1 (this, x) result ($NAME)

    class(mech_coeffs_evol_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    real(WP)                              :: $NAME

    ! Interpolate $NAME

    $NAME = this%sp_$NAME%interp(x)

    ! Finish

    return

  end function get_${NAME}_1

!****

  function get_${NAME}_v (this, x) result ($NAME)

    class(mech_coeffs_evol_t), intent(in) :: this
    real(WP), intent(in)                  :: x(:)
    real(WP)                              :: $NAME(SIZE(x))

    ! Interpolate $NAME

    $NAME = this%sp_$NAME%interp(x)

    ! Finish

    return

  end function get_${NAME}_v

  $endsub

  $PROC(V)
  $PROC(As)
  $PROC(U)
  $PROC(c_1)
  $PROC(Gamma_1)

!****

  function conv_freq (this, freq, from_units, to_units)

    class(mech_coeffs_evol_t), intent(in) :: this
    real(WP), intent(in)                  :: freq
    character(LEN=*), intent(in)          :: from_units
    character(LEN=*), intent(in)          :: to_units
    real(WP)                              :: conv_freq

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
         freq_scale = 1._WP/(TWOPI*this%t_dyn)
      case('UHZ')
         freq_scale = 1.E6_WP/(TWOPI*this%t_dyn)
      case default
         $ABORT(Invalid units)
      end select

      ! Finish

      return

    end function freq_scale

  end function conv_freq

end module gyre_mech_coeffs_evol
