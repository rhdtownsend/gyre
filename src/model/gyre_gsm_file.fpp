! Module   : gyre_gsm_file
! Purpose  : read GSM (GYRE Stellar Model) files
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

module gyre_gsm_file

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

  public :: read_gsm_model
  public :: read_gsm_data

  ! Procedures

contains

  subroutine read_gsm_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP)                    :: L_star
    real(WP), allocatable       :: r(:)
    real(WP), allocatable       :: M_r(:)
    real(WP), allocatable       :: L_r(:)
    real(WP), allocatable       :: P(:)
    real(WP), allocatable       :: rho(:)
    real(WP), allocatable       :: T(:)
    real(WP), allocatable       :: N2(:)
    real(WP), allocatable       :: Gamma_1(:)
    real(WP), allocatable       :: nabla_ad(:)
    real(WP), allocatable       :: delta(:)
    real(WP), allocatable       :: nabla(:)
    real(WP), allocatable       :: kap(:)
    real(WP), allocatable       :: kap_rho(:)
    real(WP), allocatable       :: kap_T(:)
    real(WP), allocatable       :: eps(:)
    real(WP), allocatable       :: eps_eps_rho(:)
    real(WP), allocatable       :: eps_eps_T(:)
    real(WP), allocatable       :: Omega_rot(:)
    integer                     :: n
    real(WP), allocatable       :: x(:)
    real(WP), allocatable       :: V_2(:)
    real(WP), allocatable       :: As(:)
    real(WP), allocatable       :: U(:)
    real(WP), allocatable       :: c_1(:)
    real(WP), allocatable       :: beta_rad(:)
    real(WP), allocatable       :: c_P(:)
    real(WP), allocatable       :: c_lum(:)
    real(WP), allocatable       :: c_rad(:)
    real(WP), allocatable       :: c_thn(:)
    real(WP), allocatable       :: c_thk(:)
    real(WP), allocatable       :: c_eps(:)
    real(WP), allocatable       :: eps_rho(:)
    real(WP), allocatable       :: eps_T(:)
    type(evol_model_t), pointer :: em

    ! Read data from the GSM-format file

    call read_gsm_data(ml_p%file, M_star, R_star, L_star, r, M_r, L_r, P, rho, T, &
                       N2, Gamma_1, nabla_ad, delta, nabla,  &
                       kap, kap_rho, kap_T, eps, eps_eps_rho, eps_eps_T, &
                       Omega_rot)

    ! Snap grid points

    x = r/R_star

    call snap_points(MAX(ml_p%dx_snap, EPSILON(0._WP)), M_r)
  
    ! Calculate dimensionless structure data

    n = SIZE(x)

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

    c_P = P*delta/(rho*T*nabla_ad)

    c_rad = 16._WP*PI*A_RADIATION*C_LIGHT*T**4*R_star*nabla*V_2/(3._WP*kap*rho*L_star)
    c_thn = c_P*SQRT(G_GRAVITY*M_star/R_star**3)/(A_RADIATION*C_LIGHT*kap*T**3)
    c_thk = 4._WP*PI*rho*T*c_P*SQRT(G_GRAVITY*M_star/R_star**3)*R_star**3/L_star
    c_eps = 4._WP*PI*rho*eps*R_star**3/L_star

    allocate(eps_rho(n))
    allocate(eps_T(n))

    where (eps /= 0._WP)
       eps_rho = eps_eps_rho/eps
       eps_T = eps_eps_T/eps
    elsewhere
       eps_rho = 0._WP
       eps_T = 0._WP
    end where

    if (ml_p%uniform_rot) then
       allocate(Omega_rot(n))
       Omega_rot = uniform_Omega_rot(ml_p, M_star, R_star)
    else
       Omega_rot = Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
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

  end subroutine read_gsm_model

  !****

  subroutine read_gsm_data (file, M_star, R_star, L_star, r, M_r, L_r, P, rho, T, &
                            N2, Gamma_1, nabla_ad, delta, nabla,  &
                            kap, kap_rho, kap_T, eps, eps_eps_rho, eps_eps_T, &
                            Omega_rot)

    character(*), intent(in)           :: file
    real(WP), intent(out)              :: M_star
    real(WP), intent(out)              :: R_star
    real(WP), intent(out)              :: L_star
    real(WP), allocatable, intent(out) :: r(:)
    real(WP), allocatable, intent(out) :: M_r(:)
    real(WP), allocatable, intent(out) :: L_r(:)
    real(WP), allocatable, intent(out) :: P(:)
    real(WP), allocatable, intent(out) :: rho(:)
    real(WP), allocatable, intent(out) :: T(:)
    real(WP), allocatable, intent(out) :: N2(:)
    real(WP), allocatable, intent(out) :: Gamma_1(:)
    real(WP), allocatable, intent(out) :: nabla_ad(:)
    real(WP), allocatable, intent(out) :: delta(:)
    real(WP), allocatable, intent(out) :: nabla(:)
    real(WP), allocatable, intent(out) :: kap(:)
    real(WP), allocatable, intent(out) :: kap_rho(:)
    real(WP), allocatable, intent(out) :: kap_T(:)
    real(WP), allocatable, intent(out) :: eps(:)
    real(WP), allocatable, intent(out) :: eps_eps_rho(:)
    real(WP), allocatable, intent(out) :: eps_eps_T(:)
    real(WP), allocatable, intent(out) :: Omega_rot(:)

    type(hgroup_t) :: hg
    integer        :: n
    integer        :: version

    ! Read data from the GSM-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from GSM file', TRIM(file)
100    format(A)
       write(OUTPUT_UNIT, 110) 'File name', TRIM(file)
110    format(3X,A,1X,A)
    endif

    hg = hgroup_t(file, OPEN_FILE)

    ! Read the header and determine the version

    call read_attr(hg, 'n', n)

    if (attr_exists(hg, 'version')) then
       call read_attr(hg, 'version', version)
    else
       version = 0
    endif

    call read_attr(hg, 'M_star', M_star)
    call read_attr(hg, 'R_star', R_star)
    call read_attr(hg, 'L_star', L_star)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'File version', version/100._WP
120    format(3X,A,1X,F4.2,1X,A)
    endif

    ! Read the data

    select case (version)
    case (0)
       call read_gsm_data_v0_00_()
    case (100)
       call read_gsm_data_v1_00_()
    case (110)
       call read_gsm_data_v1_10_()
    case default
       $ABORT(Unrecognized GSM file version)
    end select

    call hg%final()

    ! Finish

    return

  contains

    subroutine read_gsm_data_v0_00_ ()

      real(WP), allocatable :: w(:)

      ! Read data from the version-0.00 file

      call read_dset_alloc(hg, 'r', r)
      call read_dset_alloc(hg, 'w', w)
      call read_dset_alloc(hg, 'L_r', L_r)
      call read_dset_alloc(hg, 'p', P)
      call read_dset_alloc(hg, 'rho', rho)
      call read_dset_alloc(hg, 'T', T)
      call read_dset_alloc(hg, 'N2', N2)
      call read_dset_alloc(hg, 'Gamma_1', Gamma_1)
      call read_dset_alloc(hg, 'nabla_ad', nabla_ad)
      call read_dset_alloc(hg, 'delta', delta)
      call read_dset_alloc(hg, 'nabla', nabla)
      call read_dset_alloc(hg, 'kappa', kap)
      call read_dset_alloc(hg, 'kappa_rho', kap_rho)
      call read_dset_alloc(hg, 'kappa_T', kap_T)
      call read_dset_alloc(hg, 'epsilon', eps)
      call read_dset_alloc(hg, 'epsilon_rho', eps_eps_rho)
      call read_dset_alloc(hg, 'epsilon_T', eps_eps_T)
      call read_dset_alloc(hg, 'Omega_rot', Omega_rot)

      M_r = w/(1._WP+w)*M_star

      ! Finish

      return

    end subroutine read_gsm_data_v0_00_

    subroutine read_gsm_data_v1_00_ ()

      real(WP), allocatable :: kap_kap_rho(:)
      real(WP), allocatable :: kap_kap_T(:)

      ! Read data from the version-1.00 file

      call read_dset_alloc(hg, 'r', r)
      call read_dset_alloc(hg, 'M_r', M_r)
      call read_dset_alloc(hg, 'L_r', L_r)
      call read_dset_alloc(hg, 'P', P)
      call read_dset_alloc(hg, 'rho', rho)
      call read_dset_alloc(hg, 'T', T)
      call read_dset_alloc(hg, 'N2', N2)
      call read_dset_alloc(hg, 'Gamma_1', Gamma_1)
      call read_dset_alloc(hg, 'nabla_ad', nabla_ad)
      call read_dset_alloc(hg, 'delta', delta)
      call read_dset_alloc(hg, 'nabla', nabla)
      call read_dset_alloc(hg, 'kap', kap)
      call read_dset_alloc(hg, 'kap_rho', kap_kap_rho)
      call read_dset_alloc(hg, 'kap_T', kap_kap_T)
      call read_dset_alloc(hg, 'eps', eps)
      call read_dset_alloc(hg, 'eps_rho', eps_eps_rho)
      call read_dset_alloc(hg, 'eps_T', eps_eps_T)
      call read_dset_alloc(hg, 'Omega_rot', Omega_rot)

      kap_rho = kap_kap_rho/kap
      kap_T = kap_kap_T/kap

      ! Finish

      return

    end subroutine read_gsm_data_v1_00_

    subroutine read_gsm_data_v1_10_ ()

      real(WP), allocatable :: kap_kap_rho(:)
      real(WP), allocatable :: kap_kap_T(:)

      ! Read data from the version-1.10 file

      call read_dset_alloc(hg, 'r', r)
      call read_dset_alloc(hg, 'M_r', M_r)
      call read_dset_alloc(hg, 'L_r', L_r)
      call read_dset_alloc(hg, 'P', P)
      call read_dset_alloc(hg, 'rho', rho)
      call read_dset_alloc(hg, 'T', T)
      call read_dset_alloc(hg, 'N2', N2)
      call read_dset_alloc(hg, 'Gamma_1', Gamma_1)
      call read_dset_alloc(hg, 'nabla_ad', nabla_ad)
      call read_dset_alloc(hg, 'delta', delta)
      call read_dset_alloc(hg, 'nabla', nabla)
      call read_dset_alloc(hg, 'kap', kap)
      call read_dset_alloc(hg, 'kap_kap_rho', kap_kap_rho)
      call read_dset_alloc(hg, 'kap_kap_T', kap_kap_T)
      call read_dset_alloc(hg, 'eps', eps)
      call read_dset_alloc(hg, 'eps_eps_rho', eps_eps_rho)
      call read_dset_alloc(hg, 'eps_eps_T', eps_eps_T)
      call read_dset_alloc(hg, 'Omega_rot', Omega_rot)

      kap_rho = kap_kap_rho/kap
      kap_T = kap_kap_T/kap

      ! Finish

      return

    end subroutine read_gsm_data_v1_10_

  end subroutine read_gsm_data

end module gyre_gsm_file
