! Module   : gyre_fgong_file
! Purpose  : read FGONG files

$include 'core.inc'

module gyre_fgong_file

  ! Uses

  use core_kinds
  use core_constants

  use gyre_mech_coeffs_evol

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_fgong_file

  ! Procedures

contains

  subroutine read_fgong_file (file, G, mc, x)

    character(LEN=*), intent(in)                     :: file
    real(WP), intent(in), optional                   :: G
    class(mech_coeffs_evol_t), intent(out), optional :: mc
    real(WP), allocatable, intent(out), optional     :: x(:)

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

    write(OUTPUT_UNIT, *) 'Reading from FGONG file ', TRIM(file)

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    ! Read the header

    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)

    read(unit, *) n, iconst, ivar, ivers

    write(OUTPUT_UNIT, *) '  Initial points :', n
    write(OUTPUT_UNIT, *) '  File version   :', ivers

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

       write(OUTPUT_UNIT, *) '  Added central point'

    endif

    ! If necessary, rescale the stellar radius and mass

    n = SIZE(r)

    if(r(n) /= R_star) then
       R_star = r(n)
       write(OUTPUT_UNIT, *) '  Forced R_star == r(n)'
    endif

    if(m(n) /= M_star) then
       M_star = m(n)
       write(OUTPUT_UNIT, *) '  Forced M_star == m(n)'
    endif

    ! Initialize the mech_coeffs

    if(PRESENT(mc)) call mc%init(G, M_star, R_star, r, m, p, rho, N2, Gamma_1)

    ! Set up the grid

    if(PRESENT(x)) x = r/R_star

    ! Finish

    return

  end subroutine read_fgong_file

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

end module gyre_fgong_file

