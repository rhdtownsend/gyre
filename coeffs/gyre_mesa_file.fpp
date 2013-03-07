! Module   : gyre_mesa_file
! Purpose  : read MESA files

$include 'core.inc'

module gyre_mesa_file

  ! Uses

  use core_kinds
  use core_constants

  use gyre_mech_coeffs_evol
  use gyre_therm_coeffs_evol

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_mesa_file

  ! Procedures

contains

  subroutine read_mesa_file (file, G, mc, tc, x)

    character(LEN=*), intent(in)                      :: file
    real(WP), intent(in), optional                    :: G
    class(mech_coeffs_evol_t), intent(out), optional  :: mc
    class(therm_coeffs_evol_t), intent(out), optional :: tc
    real(WP), allocatable, intent(out), optional      :: x(:)

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
    real(WP), allocatable :: T(:)
    real(WP), allocatable :: rho(:)
    real(WP), allocatable :: nabla(:)
    real(WP), allocatable :: N2(:)
    real(WP), allocatable :: Gamma_1(:)
    real(WP), allocatable :: alpha_T(:)
    real(WP), allocatable :: c_p(:)
    real(WP), allocatable :: kappa(:)
    real(WP), allocatable :: kappa_T(:)
    real(WP), allocatable :: kappa_rho(:)
    real(WP), allocatable :: epsilon(:)
    real(WP), allocatable :: epsilon_T(:)
    real(WP), allocatable :: epsilon_rho(:)

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
    T = var(5,:)
    rho = var(6,:)
    nabla = var(7,:)
    N2 = var(8,:)
    Gamma_1 = var(12,:)*var(10,:)/var(9,:)
    alpha_T = var(11,:)/var(12,:)
    c_p = var(10,:)
    kappa = var(13,:)
    kappa_T = var(14,:)
    kappa_rho = var(15,:)
    epsilon = var(16,:)
    epsilon_T = var(17,:)/var(16,:)
    epsilon_rho = var(18,:)/var(16,:)

    ! If necessary, add central data

    if(r(1) > 0._WP) then

       m = [0._WP,m]
       N2 = [0._WP,N2]

       call add_center(r, p)
       call add_center(r, T)
       call add_center(r, rho)
       call add_center(r, nabla)
       call add_center(r, Gamma_1)
       call add_center(r, alpha_T)
       call add_center(r, c_p)
       call add_center(r, kappa)
       call add_center(r, kappa_T)
       call add_center(r, kappa_rho)
       call add_center(r, epsilon)
       call add_center(r, epsilon_T)
       call add_center(r, epsilon_rho)

       r = [0._WP,r]

    endif

    ! Initialize the mech_coeffs

    if(PRESENT(mc)) call mc%init(G, R_star, M_star, r, m, p, rho, N2, Gamma_1)

    ! Initialize the therm_coeffs

    if(PRESENT(tc)) call tc%init(G, R_star, M_star, L_star, r, m, p, T, rho, &
                                 nabla, Gamma_1, alpha_T, c_p, &
                                 kappa, kappa_T, kappa_rho, epsilon, epsilon_T, epsilon_rho)

    ! Set up the grid

    if(PRESENT(x)) x = r/R_star

    ! Finish

    return

  end subroutine read_mesa_file

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

end module gyre_mesa_file
