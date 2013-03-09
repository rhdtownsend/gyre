! Module   : gyre_therm_coeffs_evol
! Purpose  : thermal structure coefficients for evolutionary models

$include 'core.inc'

module gyre_therm_coeffs_evol

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

  type, extends(therm_coeffs_t) :: therm_coeffs_evol_t
     $VAR_DECL(V_x2)
     $VAR_DECL(c_rad)
     $VAR_DECL(c_gen)
     $VAR_DECL(c_thm)
     $VAR_DECL(nabla)
     $VAR_DECL(nabla_ad)
     $VAR_DECL(alpha_T)
     $VAR_DECL(kappa_ad)
     $VAR_DECL(kappa_S)
     $VAR_DECL(epsilon_ad)
     $VAR_DECL(epsilon_S)
     real(WP) :: t_thm
   contains
     procedure :: init
     $if($MPI)
     procedure :: bcast => bcast_tc
     $endif
     $PROC_DECL(V_x2)
     $PROC_DECL(c_rad)
     $PROC_DECL(dc_rad)
     $PROC_DECL(c_gen)
     $PROC_DECL(c_thm)
     $PROC_DECL(nabla)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(dnabla_ad)
     $PROC_DECL(alpha_T)
     $PROC_DECL(kappa_ad)
     $PROC_DECL(kappa_S)
     $PROC_DECL(epsilon_ad)
     $PROC_DECL(epsilon_S)
  end type therm_coeffs_evol_t

  ! Access specifiers

  private

  public :: therm_coeffs_evol_t

  ! Procedures

contains 

  subroutine init (this, G, R_star, M_star, L_star, r, m, p, T, rho, &
                   nabla, Gamma_1, alpha_T, c_p, &
                   kappa, kappa_T, kappa_rho, epsilon, epsilon_T, epsilon_rho)

    class(therm_coeffs_evol_t), intent(out) :: this
    real(WP), intent(in), optional          :: G
    real(WP), intent(in)                    :: R_star
    real(WP), intent(in)                    :: M_star
    real(WP), intent(in)                    :: L_star
    real(WP), intent(in)                    :: r(:)
    real(WP), intent(in)                    :: m(:)
    real(WP), intent(in)                    :: p(:)
    real(WP), intent(in)                    :: T(:)
    real(WP), intent(in)                    :: rho(:)
    real(WP), intent(in)                    :: nabla(:)
    real(WP), intent(in)                    :: Gamma_1(:)
    real(WP), intent(in)                    :: alpha_T(:)
    real(WP), intent(in)                    :: c_p(:)
    real(WP), intent(in)                    :: kappa(:)
    real(WP), intent(in)                    :: kappa_T(:)
    real(WP), intent(in)                    :: kappa_rho(:)
    real(WP), intent(in)                    :: epsilon(:)
    real(WP), intent(in)                    :: epsilon_T(:)
    real(WP), intent(in)                    :: epsilon_rho(:)

    real(WP) :: G_
    integer  :: n
    real(WP) :: V_x2(SIZE(r))
    real(WP) :: nabla_ad(SIZE(r))
    real(WP) :: c_rad(SIZE(r))
    real(WP) :: c_gen(SIZE(r))
    real(WP) :: c_thm(SIZE(r))
    real(WP) :: kappa_ad(SIZE(r))
    real(WP) :: kappa_S(SIZE(r))
    real(WP) :: epsilon_ad(SIZE(r))
    real(WP) :: epsilon_S(SIZE(r))
    real(WP) :: x(SIZE(r))

    $CHECK_BOUNDS(SIZE(m),SIZE(r))
    $CHECK_BOUNDS(SIZE(p),SIZE(r))
    $CHECK_BOUNDS(SIZE(rho),SIZE(r))
    $CHECK_BOUNDS(SIZE(T),SIZE(r))
    $CHECK_BOUNDS(SIZE(nabla),SIZE(r))
    $CHECK_BOUNDS(SIZE(Gamma_1),SIZE(r))
    $CHECK_BOUNDS(SIZE(alpha_T),SIZE(r))
    $CHECK_BOUNDS(SIZE(c_p),SIZE(r))
    $CHECK_BOUNDS(SIZE(kappa),SIZE(r))
    $CHECK_BOUNDS(SIZE(kappa_T),SIZE(r))
    $CHECK_BOUNDS(SIZE(kappa_rho),SIZE(r))
    $CHECK_BOUNDS(SIZE(epsilon),SIZE(r))
    $CHECK_BOUNDS(SIZE(epsilon_T),SIZE(r))
    $CHECK_BOUNDS(SIZE(epsilon_rho),SIZE(r))

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
       V_x2 = G_*m*rho/(p*r*(r/R_star)**2)
    elsewhere
       V_x2 = 4._WP*PI*G_*rho**2*R_star**2/(3._WP*p)
    end where

    nabla_ad = p*alpha_T/(rho*T*c_p)

    c_rad = 16._WP*PI*A_RADIATION*C_LIGHT*T**4*R_star*nabla*V_x2/(3._WP*kappa*rho*L_star)
    c_gen = 4._WP*PI*rho*epsilon*R_star**3/L_star
    c_thm = 4._WP*PI*rho*T*c_P*SQRT(G*M_star/R_star**3)*R_star**3/L_star

    kappa_ad = nabla_ad*kappa_T + kappa_rho/Gamma_1
    kappa_S = kappa_T - alpha_T*kappa_rho

    epsilon_ad = nabla_ad*epsilon_T + epsilon_rho/Gamma_1
    epsilon_S = epsilon_T - alpha_T*epsilon_rho

    x = r/R_star

    ! Initialize the therm_coeffs

    call this%sp_V_x2%init(x, V_x2, dy_dx_a=0._WP)
    call this%sp_c_rad%init(x, c_rad, dy_dx_a=0._WP)
    call this%sp_c_gen%init(x, c_gen, dy_dx_a=0._WP)
    call this%sp_c_thm%init(x, c_thm, dy_dx_a=0._WP)
    call this%sp_nabla%init(x, nabla, dy_dx_a=0._WP)
    call this%sp_nabla_ad%init(x, nabla_ad, dy_dx_a=0._WP)
    call this%sp_alpha_T%init(x, alpha_T, dy_dx_a=0._WP)
    call this%sp_kappa_S%init(x, kappa_S, dy_dx_a=0._WP)
    call this%sp_kappa_ad%init(x, kappa_ad, dy_dx_a=0._WP)
    call this%sp_epsilon_S%init(x, epsilon_S, dy_dx_a=0._WP)
    call this%sp_epsilon_ad%init(x, epsilon_ad, dy_dx_a=0._WP)

    this%t_thm = SQRT(G*M_star/R_star**3)

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_tc (this, root_rank)

    class(therm_coeffs_evol_t), intent(inout) :: this
    integer, intent(in)                       :: root_rank

    ! Broadcast the therm_coeffs

    call this%sp_V_x2%bcast(root_rank)
    call this%sp_c_rad%bcast(root_rank)
    call this%sp_c_gen%bcast(root_rank)
    call this%sp_c_thm%bcast(root_rank)
    call this%sp_nabla%bcast(root_rank)
    call this%sp_nabla_ad%bcast(root_rank)
    call this%sp_alpha_T%bcast(root_rank)
    call this%sp_kappa_S%bcast(root_rank)
    call this%sp_kappa_ad%bcast(root_rank)
    call this%sp_epsilon_S%bcast(root_rank)
    call this%sp_epsilon_ad%bcast(root_rank)

    ! Finish

    return

  end subroutine bcast_tc

  $endif

!****

  $define $PROC $sub

  $local $NAME $1

  function get_${NAME}_1 (this, x) result ($NAME)

    class(therm_coeffs_evol_t), intent(in) :: this
    real(WP), intent(in)                   :: x
    real(WP)                               :: $NAME

    ! Interpolate $NAME

    $NAME = this%sp_$NAME%interp(x)

    ! Finish

    return

  end function get_${NAME}_1

!****

  function get_${NAME}_v (this, x) result ($NAME)

    class(therm_coeffs_evol_t), intent(in) :: this
    real(WP), intent(in)                   :: x(:)
    real(WP)                               :: $NAME(SIZE(x))

    ! Interpolate $NAME

    $NAME = this%sp_$NAME%interp(x)

    ! Finish

    return

  end function get_${NAME}_v

  $endsub

  $PROC(V_x2)
  $PROC(c_rad)
  $PROC(c_gen)
  $PROC(c_thm)
  $PROC(nabla)
  $PROC(nabla_ad)
  $PROC(alpha_T)
  $PROC(kappa_S)
  $PROC(kappa_ad)
  $PROC(epsilon_S)
  $PROC(epsilon_ad)

!****

  $define $DPROC $sub

  $local $NAME $1

  function get_d${NAME}_1 (this, x) result (d$NAME)

    class(therm_coeffs_evol_t), intent(in) :: this
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

    class(therm_coeffs_evol_t), intent(in) :: this
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
  $DPROC(nabla_ad)

end module gyre_therm_coeffs_evol
