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
  use core_hgroup

  use gyre_constants
  use gyre_model
  use gyre_model_evol
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_b3_file

  ! Procedures

contains

  subroutine read_b3_file (file, deriv_type, ml, x)

    character(LEN=*), intent(in)                 :: file
    character(LEN=*), intent(in)                 :: deriv_type
    type(model_evol_t), intent(out)              :: ml
    real(WP), allocatable, intent(out), optional :: x(:)

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
    real(WP), allocatable :: c_V(:)
    real(WP), allocatable :: c_p(:)
    real(WP), allocatable :: chi_rho(:)
    real(WP), allocatable :: chi_T(:)
    real(WP), allocatable :: nabla(:)
    real(WP), allocatable :: kappa(:)
    real(WP), allocatable :: kappa_rho(:)
    real(WP), allocatable :: kappa_T(:)
    real(WP), allocatable :: epsilon(:)
    real(WP), allocatable :: epsilon_rho(:)
    real(WP), allocatable :: epsilon_T(:)
    real(WP), allocatable :: m(:)
    real(WP), allocatable :: Gamma_1(:)
    real(WP), allocatable :: nabla_ad(:)
    real(WP), allocatable :: delta(:)
    logical               :: add_center

    ! Read the model from the B3-format file

    if(check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from B3 file', TRIM(file)
100    format(A,1X,A)
    endif

    hg = hgroup_t(file, OPEN_FILE)

    ! Read the header

    call read_attr(hg, 'n_shells', n)

    call read_attr(hg, 'R_star', R_star)
    call read_attr(hg, 'M_star', M_star)
    call read_attr(hg, 'L_star', L_star)

    ! Read the data

    call read_dset_alloc(hg, 'r', r)
    call read_dset_alloc(hg, 'w', w)
    call read_dset_alloc(hg, 'p', p)
    call read_dset_alloc(hg, 'rho', rho)
    call read_dset_alloc(hg, 'T', T)
    call read_dset_alloc(hg, 'nabla', nabla)
    call read_dset_alloc(hg, 'N2', N2)
    call read_dset_alloc(hg, 'c_V', c_V)
    call read_dset_alloc(hg, 'c_p', c_p)
    call read_dset_alloc(hg, 'chi_rho', chi_rho)
    call read_dset_alloc(hg, 'chi_T', chi_T)
    call read_dset_alloc(hg, 'epsilon', epsilon)
    call read_dset_alloc(hg, 'epsilon_rho', epsilon_rho)
    call read_dset_alloc(hg, 'epsilon_T', epsilon_T)
    call read_dset_alloc(hg, 'kappa', kappa)
    call read_dset_alloc(hg, 'kappa_rho', kappa_rho)
    call read_dset_alloc(hg, 'kappa_T', kappa_T)

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
    
    epsilon_rho = epsilon_rho*epsilon
    epsilon_T = epsilon_T*epsilon

    m = [w(:n-1)/(1._WP+w(:n-1))*M_star,M_star]

    Gamma_1 = chi_rho*c_p/c_V
    delta = chi_T/chi_rho
    nabla_ad = p*delta/(rho*T*c_p)

    add_center = r(1) /= 0._WP .OR. m(1) /= 0._WP

    if(add_center .AND. check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Adding central point'
110    format(2X,A)
    endif

    ! Initialize the model

    ml = model_evol_t(M_star, R_star, L_star, r, m, p, rho, T, &
                      N2, Gamma_1, nabla_ad, delta, SPREAD(0._WP, DIM=1, NCOPIES=n), &
                      nabla, kappa, kappa_rho, kappa_T, &
                      epsilon, epsilon_rho, epsilon_T, &
                       deriv_type, add_center)

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

  end subroutine read_b3_file

end module gyre_b3_file
