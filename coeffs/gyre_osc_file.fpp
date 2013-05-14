! Module   : gyre_osc_file
! Purpose  : read OSC files
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

module gyre_osc_file

  ! Uses

  use core_kinds
  use core_constants

  use gyre_base_coeffs
  use gyre_therm_coeffs
  use gyre_evol_base_coeffs
  use gyre_evol_therm_coeffs

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_osc_file

  ! Procedures

contains

  subroutine read_osc_file (file, G, deriv_type, bc, tc, x)

    character(LEN=*), intent(in)                              :: file
    real(WP), intent(in)                                      :: G
    character(LEN=*), intent(in)                              :: deriv_type
    class(base_coeffs_t), allocatable, intent(out)            :: bc
    class(therm_coeffs_t), allocatable, intent(out), optional :: tc
    real(WP), allocatable, intent(out), optional              :: x(:)

    integer               :: unit
    integer               :: n
    integer               :: iconst
    integer               :: ivar
    integer               :: iabund
    integer               :: ivers
    real(WP), allocatable :: glob(:)
    real(WP), allocatable :: var(:,:)
    real(WP)              :: M_star
    real(WP)              :: R_star
    real(WP)              :: L_star
    real(WP), allocatable :: r(:)
    real(WP), allocatable :: m(:)
    real(WP), allocatable :: p(:)
    real(WP), allocatable :: rho(:) 
    real(WP), allocatable :: T(:) 
    real(WP), allocatable :: N2(:)
    real(WP), allocatable :: Gamma_1(:)
    real(WP), allocatable :: nabla_ad(:)
    real(WP), allocatable :: delta(:)
    real(WP), allocatable :: nabla(:)
    real(WP), allocatable :: kappa(:)
    real(WP), allocatable :: kappa_rho(:)
    real(WP), allocatable :: kappa_T(:)
    real(WP), allocatable :: epsilon_(:)
    real(WP), allocatable :: epsilon_rho(:)
    real(WP), allocatable :: epsilon_T(:)

    ! Read the model from the OSC-format file

    write(OUTPUT_UNIT, *) 'Reading from OSC file ', TRIM(file)

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    ! Read the header

    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)

    read(unit, *) n, iconst, ivar, iabund, ivers

    write(OUTPUT_UNIT, *) '  Initial points :', n
    write(OUTPUT_UNIT, *) '  File version   :', ivers

    ! Read the data

    allocate(glob(iconst))
    allocate(var(ivar,n))

    read(unit, 100) glob
    read(unit, 100) var

100 format(1P5E19.12)

    close(unit)

    var = var(:,n:1:-1)

    M_star = glob(1)
    R_star = glob(2)
    L_star = glob(3)

    r = var(1,:)
    m = EXP(var(2,:))*M_star
    T = var(3,:)
    p = var(4,:)
    rho = var(5,:)
    N2 = G*m*var(15,:)/r**3
    Gamma_1 = var(10,:)
    nabla_ad = var(11,:)
    delta = var(12,:)
    nabla = var(6,:)
    kappa = var(8,:)
    kappa_T = var(17,:)
    kappa_rho = var(18,:)
    epsilon_ = var(9,:)

    allocate(epsilon_T(n))
    allocate(epsilon_rho(n))

    where(var(9,:) /= 0._WP)
       epsilon_T = var(19,:)/var(9,:)
       epsilon_rho = var(20,:)/var(9,:)
    elsewhere
       epsilon_T = 0._WP
       epsilon_rho = 0._WP
    endwhere

    ! If necessary, add central data

    if(r(1)/R_star < EPSILON(0._WP)) r(1) = 0._WP
    if(m(1)/M_star < EPSILON(0._WP)) m(1) = 0._WP

    if(r(1) /= 0._WP) then

       m = [0._WP,m]
       N2 = [0._WP,N2]

       call add_center(r, p)
       call add_center(r, rho)
       call add_center(r, T)
       call add_center(r, Gamma_1)
       call add_center(r, nabla_ad)
       call add_center(r, delta)
       call add_center(r, nabla)
       call add_center(r, kappa)
       call add_center(r, kappa_rho)
       call add_center(r, kappa_T)
       call add_center(r, epsilon_)
       call add_center(r, epsilon_rho)
       call add_center(r, epsilon_T)

       r = [0._WP,r]

       write(OUTPUT_UNIT, *) '  Added central point'

    endif

    ! Initialize the base_coeffs

    allocate(evol_base_coeffs_t::bc)

    select type (bc)
    type is (evol_base_coeffs_t)
       call bc%init(G, M_star, R_star, L_star, r, m, p, rho, T, &
                    N2, Gamma_1, nabla_ad, delta, deriv_type)
    class default
       $ABORT(Invalid bc type)
    end select

    ! Initialize the therm_coeffs

    if(PRESENT(tc)) then

       allocate(evol_therm_coeffs_t::tc)

       select type (tc)
       type is (evol_therm_coeffs_t)
          call tc%init(G, M_star, R_star, L_star, r, m, p, rho, T, &
                       Gamma_1, nabla_ad, delta, nabla,  &
                       kappa, kappa_rho, kappa_T, &
                       epsilon_, epsilon_rho, epsilon_T, deriv_type)
       class default
          $ABORT(Invalid tc type)
       end select

    endif

    ! Set up the grid

    if(PRESENT(x)) x = r/R_star

    ! Finish

    return

  end subroutine read_osc_file

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

end module gyre_osc_file
