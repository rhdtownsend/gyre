! Module   : gyre_b3_file
! Purpose  : read B3 files
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

module gyre_b3_file

  ! Uses

  use core_kinds
  use core_constants
  use core_hgroup

  use gyre_mech_coeffs_evol
  use gyre_therm_coeffs_evol

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_b3_file

  ! Procedures

contains

  subroutine read_b3_file (file, G, mc, tc, x)

    character(LEN=*), intent(in)                      :: file
    real(WP), intent(in), optional                    :: G
    class(mech_coeffs_evol_t), intent(out)            :: mc
    class(therm_coeffs_evol_t), intent(out), optional :: tc
    real(WP), allocatable, intent(out), optional      :: x(:)

    type(hgroup_t)        :: hg
    integer               :: n
    real(WP)              :: M_star
    real(WP)              :: R_star
    real(WP)              :: L_star
    real(WP), allocatable :: r(:)
    real(WP), allocatable :: w(:)
    real(WP), allocatable :: p(:)
    real(WP), allocatable :: T(:)
    real(WP), allocatable :: rho(:)
    real(WP), allocatable :: nabla(:)
    real(WP), allocatable :: N2(:)
    real(WP), allocatable :: c_V(:)
    real(WP), allocatable :: c_p(:)
    real(WP), allocatable :: chi_T(:)
    real(WP), allocatable :: chi_rho(:)
    real(WP), allocatable :: kappa(:)
    real(WP), allocatable :: kappa_T(:)
    real(WP), allocatable :: kappa_rho(:)
    real(WP), allocatable :: epsilon(:)
    real(WP), allocatable :: epsilon_T(:)
    real(WP), allocatable :: epsilon_rho(:)
    real(WP), allocatable :: m(:)
    real(WP), allocatable :: Gamma_1(:)
    real(WP), allocatable :: alpha_T(:)

    ! Read the model from the B3-format file

    call hg%init(file, OPEN_FILE)

    ! Read the header

    call read_attr(hg, 'n_shells', n)

    call read_attr(hg, 'R_star', R_star)
    call read_attr(hg, 'M_star', M_star)
    call read_attr(hg, 'L_star', L_star)

    ! Read the data

    call read_dset(hg, 'r', r, alloc=.TRUE.)
    call read_dset(hg, 'w', w, alloc=.TRUE.)
    call read_dset(hg, 'p', p, alloc=.TRUE.)
    call read_dset(hg, 'T', T, alloc=.TRUE.)
    call read_dset(hg, 'rho', rho, alloc=.TRUE.)
    call read_dset(hg, 'nabla', nabla, alloc=.TRUE.)
    call read_dset(hg, 'N2', N2, alloc=.TRUE.)
    call read_dset(hg, 'c_V', c_V, alloc=.TRUE.)
    call read_dset(hg, 'c_p', c_p, alloc=.TRUE.)
    call read_dset(hg, 'chi_T', chi_T, alloc=.TRUE.)
    call read_dset(hg, 'chi_rho', chi_rho, alloc=.TRUE.)
    call read_dset(hg, 'epsilon', epsilon, alloc=.TRUE.)
    call read_dset(hg, 'epsilon_T', epsilon_T, alloc=.TRUE.)
    call read_dset(hg, 'epsilon_rho', epsilon_rho, alloc=.TRUE.)
    call read_dset(hg, 'kappa', kappa, alloc=.TRUE.)
    call read_dset(hg, 'kappa_T', kappa_T, alloc=.TRUE.)
    call read_dset(hg, 'kappa_rho', kappa_rho, alloc=.TRUE.)

    call hg%final()

    R_star = R_star*1.E2_WP
    M_star = M_star*1.E3_WP
    L_star = L_star*1.E7_WP

    r = r*1.E2_WP
    p = p*1.E1_WP
    rho = rho*1.E-3_WP
    c_V = c_V*1.E4_WP
    c_p = c_p*1.E4_WP
    kappa = kappa*1.E1_WP
    epsilon = epsilon*1.E4_WP
    
    m = [w(:n-1)/(1._WP+w(:n-1))*M_star,M_star]

    Gamma_1 = chi_rho*c_p/c_V
    alpha_T = chi_T/chi_rho

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

    call mc%init(G, R_star, M_star, r, m, p, rho, N2, Gamma_1)

    ! Initialize the therm_coeffs

    if(PRESENT(tc)) call tc%init(G, R_star, M_star, L_star, r, m, p, T, rho, &
                                 nabla, Gamma_1, alpha_T, c_p, &
                                 kappa, kappa_T, kappa_rho, epsilon, epsilon_T, epsilon_rho)

    ! Set up the grid

    if(PRESENT(x)) x = r/R_star

    ! Finish

    return

  end subroutine read_b3_file

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

end module gyre_b3_file
