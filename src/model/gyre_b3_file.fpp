! Module   : gyre_b3_file
! Purpose  : read B3 files
!
! Copyright 2013-2017 Rich Townsend
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
  use gyre_evol_model
  use gyre_model
  use gyre_model_par
  use gyre_model_util
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_b3_model

  ! Procedures

contains

  subroutine read_b3_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    type(hgroup_t)              :: hg
    integer                     :: n
    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP)                    :: L_star
    real(WP), allocatable       :: r(:)
    real(WP), allocatable       :: w(:)
    real(WP), allocatable       :: M_r(:)
    real(WP), allocatable       :: L_r(:)
    real(WP), allocatable       :: P(:)
    real(WP), allocatable       :: rho(:)
    real(WP), allocatable       :: T(:)
    real(WP), allocatable       :: N2(:)
    real(WP), allocatable       :: c_V(:)
    real(WP), allocatable       :: c_P(:)
    real(WP), allocatable       :: chi_rho(:)
    real(WP), allocatable       :: chi_T(:)
    real(WP), allocatable       :: nabla(:)
    real(WP), allocatable       :: kap(:)
    real(WP), allocatable       :: kap_rho(:)
    real(WP), allocatable       :: kap_T(:)
    real(WP), allocatable       :: eps(:)
    real(WP), allocatable       :: eps_rho(:)
    real(WP), allocatable       :: eps_T(:)
    real(WP), allocatable       :: Gamma_1(:)
    real(WP), allocatable       :: nabla_ad(:)
    real(WP), allocatable       :: delta(:)
    real(WP), allocatable       :: x(:)
    real(WP), allocatable       :: V_2(:)
    real(WP), allocatable       :: As(:)
    real(WP), allocatable       :: U(:)
    real(WP), allocatable       :: c_1(:)
    real(WP), allocatable       :: beta_rad(:)
    real(WP), allocatable       :: c_lum(:)
    real(WP), allocatable       :: c_rad(:)
    real(WP), allocatable       :: c_thn(:)
    real(WP), allocatable       :: c_thk(:)
    real(WP), allocatable       :: c_eps(:)
    real(WP), allocatable       :: Omega_rot(:)
    type(evol_model_t), pointer :: em

    ! Open the B3-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from B3 file', TRIM(ml_p%file)
100    format(A,1X,A)
    endif

    hg = hgroup_t(ml_p%file, OPEN_FILE)

    ! Read the header

    call read_attr(hg, 'n_shells', n)

    call read_attr(hg, 'R_star', R_star)
    call read_attr(hg, 'M_star', M_star)
    call read_attr(hg, 'L_star', L_star)

     if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Grid points  :', n
110    format(3X,A,1X,I0)
    endif

    ! Read the data

    call read_dset_alloc(hg, 'r', r)
    call read_dset_alloc(hg, 'w', w)
    call read_dset_alloc(hg, 'L_r', L_r)
    call read_dset_alloc(hg, 'p', P)
    call read_dset_alloc(hg, 'rho', rho)
    call read_dset_alloc(hg, 'T', T)
    call read_dset_alloc(hg, 'nabla', nabla)
    call read_dset_alloc(hg, 'N2', N2)
    call read_dset_alloc(hg, 'c_V', c_V)
    call read_dset_alloc(hg, 'c_p', c_P)
    call read_dset_alloc(hg, 'chi_rho', chi_rho)
    call read_dset_alloc(hg, 'chi_T', chi_T)
    call read_dset_alloc(hg, 'epsilon', eps)
    call read_dset_alloc(hg, 'epsilon_rho', eps_rho)
    call read_dset_alloc(hg, 'epsilon_T', eps_T)
    call read_dset_alloc(hg, 'kappa', kap)
    call read_dset_alloc(hg, 'kappa_rho', kap_rho)
    call read_dset_alloc(hg, 'kappa_T', kap_T)

    call hg%final()

    R_star = R_star*1.E2_WP
    M_star = M_star*1.E3_WP
    L_star = L_star*1.E7_WP

    r = r*1.E2_WP
    M_r = w/(1._WP+w)*M_star
    L_r = L_r*1.E7_WP

    P = p*1.E1_WP
    rho = rho*1.E-3_WP
    
    c_V = c_V*1.E4_WP
    c_P = c_P*1.E4_WP
    kap = kap*1.E1_WP
    eps = eps*1.E4_WP
    
    Gamma_1 = chi_rho*c_p/c_V
    delta = chi_T/chi_rho
    nabla_ad = p*delta/(rho*T*c_p)

    ! Snap grid points

    x = r/R_star

    call snap_points(MAX(ml_p%dx_snap, EPSILON(0._WP)), x, M_r)
  
    ! Calculate dimensionless structure data

    allocate(V_2(n))
    allocate(As(n))
    allocate(U(n))
    allocate(c_1(n))
    allocate(c_lum(n))

    where (x /= 0._WP)
       V_2 = G_GRAVITY*M_r*rho/(P*r*x**2)
       As = r**3*N2/(G_GRAVITY*M_r)
       U = 4._WP*PI*rho*r**3/M_r
       c_1 = (r/R_star)**3/(M_r/M_star)
       c_lum = (L_r/L_star)/x**3
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._WP*P(1))
       As = 0._WP
       U = 3._WP
       c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
       c_lum = 4._WP*PI*rho(1)*eps(1)*R_star**3/L_star
    end where

    beta_rad = A_RADIATION*T**4/(3._WP*P)

    c_rad = 16._WP*PI*A_RADIATION*C_LIGHT*T**4*R_star*nabla*V_2/(3._WP*kap*rho*L_star)
    c_thn = c_P*SQRT(G_GRAVITY*M_star/R_star**3)/(A_RADIATION*C_LIGHT*kap*T**3)
    c_thk = 4._WP*PI*rho*T*c_P*SQRT(G_GRAVITY*M_star/R_star**3)*R_star**3/L_star
    c_eps = 4._WP*PI*rho*eps*R_star**3/L_star

    allocate(Omega_rot(n))

    if (ml_p%uniform_rot) then
       Omega_rot = uniform_Omega_rot(ml_p, M_star, R_star)
    else
       Omega_rot = 0._WP
    endif

    ! Initialize the evol_model_t

    allocate(em, SOURCE=evol_model_t(x, M_star, R_star, L_star, ml_p))

    call em%define(I_V_2, V_2)
    call em%define(I_AS, As)
    call em%define(I_U, U)
    call em%define(I_C_1, c_1)

    call em%define(I_GAMMA_1, Gamma_1)
    call em%define(I_DELTA, delta)
    call em%define(I_NABLA_AD, nabla_ad)
    call em%define(I_NABLA, nabla)
    call em%define(I_BETA_RAD, beta_rad)

    call em%define(I_C_LUM, c_lum)
    call em%define(I_C_RAD, c_rad)
    call em%define(I_C_THN, c_thn)
    call em%define(I_C_THK, c_thk)
    call em%define(I_C_EPS, c_eps)

    call em%define(I_EPS_RHO, eps_rho)
    call em%define(I_EPS_T, eps_T)

    call em%define(I_KAP_RHO, kap_rho)
    call em%define(I_KAP_T, kap_T)

    call em%define(I_OMEGA_ROT, Omega_rot)

    ! Return a pointer

    ml => em

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, *)
    endif

    ! Finish

    return

  end subroutine read_b3_model

end module gyre_b3_file
