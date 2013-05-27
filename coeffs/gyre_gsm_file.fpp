! Module   : gyre_gsm_file
! Purpose  : read GSM (GYRE Stellar Model) files
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

module gyre_gsm_file

  ! Uses

  use core_kinds
  use core_constants
  use core_hgroup

  use gyre_base_coeffs
  use gyre_therm_coeffs
  use gyre_evol_base_coeffs
  use gyre_evol_therm_coeffs
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_gsm_file

  ! Procedures

contains

  subroutine read_gsm_file (file, G, deriv_type, bc, tc, x)

    character(LEN=*), intent(in)                              :: file
    real(WP), intent(in)                                      :: G
    character(LEN=*), intent(in)                              :: deriv_type
    class(base_coeffs_t), allocatable, intent(out)            :: bc
    class(therm_coeffs_t), allocatable, intent(out), optional :: tc
    real(WP), allocatable, intent(out), optional              :: x(:)

    type(hgroup_t)        :: hg
    integer               :: n
    real(WP)              :: M_star
    real(WP)              :: R_star
    real(WP)              :: L_star
    real(WP), allocatable :: r(:)
    real(WP), allocatable :: w(:)
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
    real(WP), allocatable :: epsilon(:)
    real(WP), allocatable :: epsilon_rho(:)
    real(WP), allocatable :: epsilon_T(:)
    real(WP), allocatable :: Omega_rot(:)
    real(WP), allocatable :: m(:)
    logical               :: add_center

    ! Read the model from the GSM-format file

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from GSM file', TRIM(file)
100    format(A,1X,A)
    endif

    call hg%init(file, OPEN_FILE)

    ! Read the header

    call read_attr(hg, 'n', n)

    call read_attr(hg, 'M_star', M_star)
    call read_attr(hg, 'R_star', R_star)
    call read_attr(hg, 'L_star', L_star)

    ! Read the data

    call read_dset_alloc(hg, 'r', r)
    call read_dset_alloc(hg, 'w', w)
    call read_dset_alloc(hg, 'p', p)
    call read_dset_alloc(hg, 'rho', rho)
    call read_dset_alloc(hg, 'T', T)
    call read_dset_alloc(hg, 'N2', N2)
    call read_dset_alloc(hg, 'Gamma_1', Gamma_1)
    call read_dset_alloc(hg, 'nabla_ad', nabla_ad)
    call read_dset_alloc(hg, 'delta', delta)
    call read_dset_alloc(hg, 'nabla', nabla)
    call read_dset_alloc(hg, 'epsilon', epsilon)
    call read_dset_alloc(hg, 'epsilon_rho', epsilon_rho)
    call read_dset_alloc(hg, 'epsilon_T', epsilon_T)
    call read_dset_alloc(hg, 'kappa', kappa)
    call read_dset_alloc(hg, 'kappa_rho', kappa_rho)
    call read_dset_alloc(hg, 'kappa_T', kappa_T)
    call read_dset_alloc(hg, 'Omega_rot', Omega_rot)

    call hg%final()

    m = [w(:n-1)/(1._WP+w(:n-1))*M_star,M_star]

    add_center = r(1) /= 0._WP .OR. m(1) /= 0._WP

    if(add_center .AND. check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Adding central point'
110    format(2X,A)
    endif

    ! Initialize the base_coeffs

    allocate(evol_base_coeffs_t::bc)

    select type (bc)
    type is (evol_base_coeffs_t)
       call bc%init(G, M_star, R_star, L_star, r, m, p, rho, T, &
                    N2, Gamma_1, nabla_ad, delta, deriv_type, add_center)
    class default
       $ABORT(Invalid bc type)
    end select

    ! Initialize the therm_coeffs

    if(PRESENT(tc)) then

       allocate(evol_therm_coeffs_t::tc)

       select type (tc)
       type is (evol_therm_coeffs_t)
          call tc%init(G, M_star, R_star, L_star, r, m, p, rho, T, &
                       Gamma_1, nabla_ad, delta, nabla, &
                       kappa, kappa_rho, kappa_T, &
                       epsilon, epsilon_rho, epsilon_T, deriv_type, add_center)
       class default
          $ABORT(Invalid tc type)
       end select

    endif

    ! Set up the grid

    if(PRESENT(x)) then
       if(add_center) then
          x = [0._WP,r/R_star]
       else
          x = r/R_star
       endif
    endif

    ! Finish

    return

  end subroutine read_gsm_file

end module gyre_gsm_file
