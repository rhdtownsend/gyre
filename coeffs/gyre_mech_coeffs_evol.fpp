! Module   : gyre_mech_coeffs_evol
! Purpose  : mechanical structure coefficients for evolutionary models

$include 'core.inc'

module gyre_mech_coeffs_evol

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel
  use core_hgroup
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
     procedure, public :: read
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

  subroutine read (this, file_type, file, G, x)

    class(mech_coeffs_evol_t), intent(out)       :: this
    character(LEN=*), intent(in)                 :: file_type
    character(LEN=*), intent(in)                 :: file
    real(WP), intent(in), optional               :: G
    real(WP), allocatable, intent(out), optional :: x(:)

    ! Read the mech_coeffs from the file

    select case(file_type)
    case('FGONG')
       call read_fgong(this, file, G, x)
    case('MESA')
       call read_mesa(this, file, G, x)
    case('B3')
       call read_b3(this, file, G, x)
    case default
       $ABORT(Invalid file_type)
    end select

    ! Finish

    return

  end subroutine read

!****

  subroutine read_fgong (mc, file, G, x)

    class(mech_coeffs_evol_t), intent(out)       :: mc
    character(LEN=*), intent(in)                 :: file
    real(WP), intent(in), optional               :: G
    real(WP), allocatable, intent(out), optional :: x(:)

    real(WP)              :: G_
    integer               :: unit
    integer               :: n
    integer               :: iconst
    integer               :: ivar
    integer               :: ivers
    real(WP), allocatable :: glob(:)
    real(WP), allocatable :: var(:,:)
    real(WP)              :: R_star
    real(WP)              :: M_star
    real(WP), allocatable :: r(:)
    real(WP), allocatable :: m(:)
    real(WP), allocatable :: p(:)
    real(WP), allocatable :: rho(:) 
    real(WP), allocatable :: N2(:)
    real(WP), allocatable :: Gamma_1(:)

    if(PRESENT(G)) then
       G_ = G
    else
       G_ = G_GRAVITY
    endif

    ! Read the model from the FGONG-format file

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    ! Read the header

    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)

    read(unit, *) n, iconst, ivar, ivers

    ! Read the data

    allocate(glob(iconst))
    allocate(var(ivar,n))

    read(unit, 100) glob
    read(unit, 100) var

100 format(1P5E16.9)

    close(unit)

    var = var(:,n:1:-1)

    R_star = glob(1)
    M_star = glob(2)

    r = var(1,:)
    m = EXP(var(2,:))*M_star
    p = var(4,:)
    rho = var(5,:)
    N2 = G_*m*var(15,:)/r**3
    Gamma_1 = var(10,:)

    ! If necessary, add central data

    if(r(1) > 0._WP) then

       m = [0._WP,m]
       call add_center(r, p)
       call add_center(r, rho)
       N2 = [0._WP,N2]
       call add_center(r, Gamma_1)

       r = [0._WP,r]

    endif

    ! Initialize the mech_coeffs

    call mc%init(G, M_star, R_star, r, m, p, rho, N2, Gamma_1)

    ! If necessary, return the grid

    if(PRESENT(x)) then
       x = r/R_star
    end if

    ! Finish

    return

  end subroutine read_fgong

!****

  subroutine read_mesa (mc, file, G, x)

    class(mech_coeffs_evol_t), intent(out)       :: mc
    character(LEN=*), intent(in)                 :: file
    real(WP), intent(in), optional               :: G
    real(WP), allocatable, intent(out), optional :: x(:)

    integer               :: unit
    integer               :: n
    real(WP)              :: M_star
    real(WP)              :: R_star
    real(WP)              :: L_star
    real(WP), allocatable :: var(:,:)
    integer               :: k
    integer               :: k_chk
    real(WP), allocatable :: r(:)
    real(WP), allocatable :: m(:)
    real(WP), allocatable :: p(:)
    real(WP), allocatable :: rho(:)
    real(WP), allocatable :: N2(:)
    real(WP), allocatable :: Gamma_1(:)

    ! Read the model from the MESA-format file

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    ! Read the header

    read(unit, *) n, M_star, R_star, L_star

    ! Read the data

    allocate(var(18,n))

    read_loop : do k = 1,n
       read(unit, *) k_chk, var(:,k)
       $ASSERT(k == k_chk,Index mismatch)
    end do read_loop

    close(unit)

    r = var(1,:)
    m = var(2,:)/(1._WP+var(2,:))*M_star
    p = var(4,:)
    rho = var(6,:)
    N2 = var(8,:)
    Gamma_1 = var(12,:)*var(10,:)/var(9,:)

    ! If necessary, add central data

    if(r(1) > 0._WP) then

       m = [0._WP,m]
       call add_center(r, p)
       call add_center(r, rho)
       N2 = [0._WP,N2]
       call add_center(r, Gamma_1)

       r = [0._WP,r]

    endif

    ! Initialize the mech_coeffs

    call mc%init(G, R_star, M_star, r, m, p, rho, N2, Gamma_1)

    ! If necessary, return the grid

    if(PRESENT(x)) then
       x = r/R_star
    end if

    ! Finish

    return

  end subroutine read_mesa

!****

  subroutine read_b3 (mc, file, G, x)

    class(mech_coeffs_evol_t), intent(out)       :: mc
    character(LEN=*), intent(in)                 :: file
    real(WP), intent(in), optional               :: G
    real(WP), allocatable, intent(out), optional :: x(:)

    type(hgroup_t)        :: hg
    integer               :: n
    real(WP)              :: M_star
    real(WP)              :: R_star
    real(WP), allocatable :: r(:)
    real(WP), allocatable :: w(:)
    real(WP), allocatable :: p(:)
    real(WP), allocatable :: rho(:)
    real(WP), allocatable :: N2(:)
    real(WP), allocatable :: c_V(:)
    real(WP), allocatable :: c_p(:)
    real(WP), allocatable :: chi_rho(:)
    real(WP), allocatable :: m(:)
    real(WP), allocatable :: Gamma_1(:)

    ! Read the model from the B3-format file

    call hg%init(file, OPEN_FILE)

    ! Read the header

    call read_attr(hg, 'n_shells', n)

    call read_attr(hg, 'R_star', R_star)
    call read_attr(hg, 'M_star', M_star)

    ! Read the data

    call read_dset(hg, 'r', r, alloc=.TRUE.)
    call read_dset(hg, 'w', w, alloc=.TRUE.)
    call read_dset(hg, 'p', p, alloc=.TRUE.)
    call read_dset(hg, 'rho', rho, alloc=.TRUE.)
    call read_dset(hg, 'N2', N2, alloc=.TRUE.)
    call read_dset(hg, 'c_V', c_V, alloc=.TRUE.)
    call read_dset(hg, 'c_p', c_p, alloc=.TRUE.)
    call read_dset(hg, 'chi_rho', chi_rho, alloc=.TRUE.)

    call hg%final()

    M_star = M_star*1.E3_WP
    R_star = R_star*1.E2_WP

    r = r*1.E2_WP
    p = p*1.E1_WP
    rho = rho*1.E-3_WP
    c_V = c_V*1.E4_WP
    c_p = c_p*1.E4_WP
    
    m = [w(:n-1)/(1._WP+w(:n-1))*M_star,M_star]
    Gamma_1 = chi_rho*c_p/c_V

    ! If necessary, add central data

    if(r(1) > 0._WP) then

       m = [0._WP,m]
       call add_center(r, p)
       call add_center(r, rho)
       N2 = [0._WP,N2]
       call add_center(r, Gamma_1)

       r = [0._WP,r]

    endif

    ! Initialize the mech_coeffs

    call mc%init(G, R_star, M_star, r, m, p, rho, N2, Gamma_1)

    ! If necessary, return the grid

    if(PRESENT(x)) then
       x = r/R_star
    end if

    ! Finish

    return

  end subroutine read_b3

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

!****

  subroutine add_center (x, y)

    real(WP), intent(in)                 :: x(:)
    real(WP), intent(inout), allocatable :: y(:)

    real(WP) :: y_0

    ! Add center (x=0) data to the array y(x), incrementing the
    ! dimension of y by 1. x is not altered.

    y_0 = (x(2)**2*y(1) - x(1)**2*y(2))/(x(2)**2 - x(1)**2)

    y = [y_0,y]

    ! Finish

    return

  end subroutine add_center

end module gyre_mech_coeffs_evol
