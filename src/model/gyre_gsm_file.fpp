! Module   : gyre_gsm_file
! Purpose  : read GSM (GYRE Stellar Model) files
!
! Copyright 2013-2014 Rich Townsend
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
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: read_gsm_model
  public :: read_gsm_data
  public :: write_gsm_data

  ! Procedures

contains

  subroutine read_gsm_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP)                    :: L_star
    real(WP), allocatable       :: r(:)
    real(WP), allocatable       :: w(:)
    real(WP), allocatable       :: p(:)
    real(WP), allocatable       :: rho(:)
    real(WP), allocatable       :: T(:)
    real(WP), allocatable       :: N2(:)
    real(WP), allocatable       :: Gamma_1(:)
    real(WP), allocatable       :: nabla_ad(:)
    real(WP), allocatable       :: delta(:)
    real(WP), allocatable       :: nabla(:)
    real(WP), allocatable       :: kappa(:)
    real(WP), allocatable       :: kappa_rho(:)
    real(WP), allocatable       :: kappa_T(:)
    real(WP), allocatable       :: epsilon(:)
    real(WP), allocatable       :: epsilon_rho(:)
    real(WP), allocatable       :: epsilon_T(:)
    real(WP), allocatable       :: Omega_rot(:)
    integer                     :: n
    real(WP), allocatable       :: x(:)
    real(WP), allocatable       :: m(:)
    real(WP), allocatable       :: V_2(:)
    real(WP), allocatable       :: As(:)
    real(WP), allocatable       :: U(:)
    real(WP), allocatable       :: c_1(:)
    real(WP), allocatable       :: beta_rad(:)
    real(WP), allocatable       :: c_P(:)
    real(WP), allocatable       :: c_rad(:)
    real(WP), allocatable       :: c_thm(:)
    real(WP), allocatable       :: c_dif(:)
    real(WP), allocatable       :: c_eps_ad(:)
    real(WP), allocatable       :: c_eps_S(:)
    real(WP), allocatable       :: kappa_ad(:)
    real(WP), allocatable       :: kappa_S(:)
    type(evol_model_t), pointer :: em

    ! Read data from the GSM-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from GSM file', TRIM(ml_p%file)
100    format(A,1X,A)
    endif

    call read_gsm_data(ml_p%file, M_star, R_star, L_star, r, w, p, rho, T, &
                       N2, Gamma_1, nabla_ad, delta, nabla,  &
                       kappa, kappa_rho, kappa_T, &
                       epsilon, epsilon_rho, epsilon_T, &
                       Omega_rot)

    ! Calculate dimensionless structure data

    x = r/R_star
    m = w/(1._WP+w)

    n = SIZE(x)

    allocate(V_2(n))
    allocate(As(n))
    allocate(U(n))
    allocate(c_1(n))

    where (x /= 0._WP)
       V_2 = G_GRAVITY*m*rho/(P*r*x**2)
       As = r**3*N2/(G_GRAVITY*m)
       U = 4._WP*PI*rho*r**3/m
       c_1 = (r/R_star)**3/(m/M_star)
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._WP*P(1))
       As = 0._WP
       U = 3._WP
       c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
    end where

    beta_rad = A_RADIATION*T**4/(3._WP*P)

    c_P = P*delta/(rho*T*nabla_ad)

    kappa_ad = nabla_ad*kappa_T + kappa_rho/Gamma_1
    kappa_S = kappa_T - delta*kappa_rho

    c_rad = 16._WP*PI*A_RADIATION*C_LIGHT*T**4*R_star*nabla*V_2/(3._WP*kappa*rho*L_star)
    c_thm = 4._WP*PI*rho*T*c_P*SQRT(G_GRAVITY*M_star/R_star**3)*R_star**3/L_star
    c_dif = (kappa_ad-4._WP*nabla_ad)*V_2*x**2*nabla + V_2*x**2*nabla_ad

    c_eps_ad = 4._WP*PI*rho*(nabla_ad*epsilon_T + epsilon_rho/Gamma_1)*R_star**3/L_star
    c_eps_S = 4._WP*PI*rho*(epsilon_T - delta*epsilon_rho)*R_star**3/L_star

    if (ml_p%uniform_rot) then
       Omega_rot = ml_p%Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
    else
       Omega_rot = Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
    endif

    ! Initialize the evol_model_t

    allocate(em, SOURCE=evol_model_t(x, M_star, R_star, L_star, ml_p))

    call em%set_V_2(V_2)
    call em%set_As(As)
    call em%set_U(U)
    call em%set_c_1(c_1)

    call em%set_Gamma_1(Gamma_1)
    call em%set_delta(delta)
    call em%set_nabla_ad(nabla_ad)
    call em%set_nabla(nabla)
    call em%set_beta_rad(beta_rad)

    call em%set_c_rad(c_rad)
    call em%set_c_thm(c_thm)
    call em%set_c_dif(c_dif)
    call em%set_c_eps_ad(c_eps_ad)
    call em%set_c_eps_S(c_eps_S)
    call em%set_kappa_ad(kappa_ad)
    call em%set_kappa_S(kappa_S)

    call em%set_Omega_rot(Omega_rot)

    ! Return a pointer

    ml => em

    ! Finish

    return

  end subroutine read_gsm_model

  !****

  subroutine read_gsm_data (file, M_star, R_star, L_star, r, w, p, rho, T, &
                            N2, Gamma_1, nabla_ad, delta, nabla,  &
                            kappa, kappa_rho, kappa_T, &
                            epsilon, epsilon_rho, epsilon_T, &
                            Omega_rot)

    character(*), intent(in)           :: file
    real(WP), intent(out)              :: M_star
    real(WP), intent(out)              :: R_star
    real(WP), intent(out)              :: L_star
    real(WP), allocatable, intent(out) :: r(:)
    real(WP), allocatable, intent(out) :: w(:)
    real(WP), allocatable, intent(out) :: p(:)
    real(WP), allocatable, intent(out) :: rho(:)
    real(WP), allocatable, intent(out) :: T(:)
    real(WP), allocatable, intent(out) :: N2(:)
    real(WP), allocatable, intent(out) :: Gamma_1(:)
    real(WP), allocatable, intent(out) :: nabla_ad(:)
    real(WP), allocatable, intent(out) :: delta(:)
    real(WP), allocatable, intent(out) :: nabla(:)
    real(WP), allocatable, intent(out) :: kappa(:)
    real(WP), allocatable, intent(out) :: kappa_rho(:)
    real(WP), allocatable, intent(out) :: kappa_T(:)
    real(WP), allocatable, intent(out) :: epsilon(:)
    real(WP), allocatable, intent(out) :: epsilon_rho(:)
    real(WP), allocatable, intent(out) :: epsilon_T(:)
    real(WP), allocatable, intent(out) :: Omega_rot(:)

    type(hgroup_t) :: hg
    integer        :: n

    ! Read data from the GSM-format file

    hg = hgroup_t(file, OPEN_FILE)

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

    $CHECK_BOUNDS(SIZE(r),n)
    $CHECK_BOUNDS(SIZE(w),n)
    $CHECK_BOUNDS(SIZE(p),n)
    $CHECK_BOUNDS(SIZE(rho),n)
    $CHECK_BOUNDS(SIZE(T),n)
    $CHECK_BOUNDS(SIZE(N2),n)
    $CHECK_BOUNDS(SIZE(Gamma_1),n)
    $CHECK_BOUNDS(SIZE(nabla_ad),n)
    $CHECK_BOUNDS(SIZE(delta),n)
    $CHECK_BOUNDS(SIZE(nabla),n)
    $CHECK_BOUNDS(SIZE(epsilon),n)
    $CHECK_BOUNDS(SIZE(epsilon_rho),n)
    $CHECK_BOUNDS(SIZE(epsilon_T),n)
    $CHECK_BOUNDS(SIZE(kappa),n)
    $CHECK_BOUNDS(SIZE(kappa_rho),n)
    $CHECK_BOUNDS(SIZE(kappa_T),n)
    $CHECK_BOUNDS(SIZE(Omega_rot),n)

    ! Finish

    return

  end subroutine read_gsm_data

!****

  subroutine write_gsm_data (file, M_star, R_star, L_star, r, w, p, rho, T, &
                             N2, Gamma_1, nabla_ad, delta, nabla,  &
                             kappa, kappa_rho, kappa_T, &
                             epsilon, epsilon_rho, epsilon_T, &
                             Omega_rot)

    character(*), intent(in) :: file
    real(WP), intent(in)     :: M_star
    real(WP), intent(in)     :: R_star
    real(WP), intent(in)     :: L_star
    real(WP), intent(in)     :: r(:)
    real(WP), intent(in)     :: w(:)
    real(WP), intent(in)     :: p(:)
    real(WP), intent(in)     :: rho(:)
    real(WP), intent(in)     :: T(:)
    real(WP), intent(in)     :: N2(:)
    real(WP), intent(in)     :: Gamma_1(:)
    real(WP), intent(in)     :: nabla_ad(:)
    real(WP), intent(in)     :: delta(:)
    real(WP), intent(in)     :: nabla(:)
    real(WP), intent(in)     :: kappa(:)
    real(WP), intent(in)     :: kappa_rho(:)
    real(WP), intent(in)     :: kappa_T(:)
    real(WP), intent(in)     :: epsilon(:)
    real(WP), intent(in)     :: epsilon_rho(:)
    real(WP), intent(in)     :: epsilon_T(:)
    real(WP), intent(in)     :: Omega_rot(:)

    type(hgroup_t) :: hg

    ! Write data to the GSM-format file

    hg = hgroup_t(file, CREATE_FILE)

    ! Write the header

    call write_attr(hg, 'n', SIZE(r))

    call write_attr(hg, 'M_star', M_star)
    call write_attr(hg, 'R_star', R_star)
    call write_attr(hg, 'L_star', L_star)

    ! Write the data

    call write_dset(hg, 'r', r)
    call write_dset(hg, 'w', w)
    call write_dset(hg, 'p', p)
    call write_dset(hg, 'rho', rho)
    call write_dset(hg, 'T', T)
    call write_dset(hg, 'N2', N2)
    call write_dset(hg, 'Gamma_1', Gamma_1)
    call write_dset(hg, 'nabla_ad', nabla_ad)
    call write_dset(hg, 'delta', delta)
    call write_dset(hg, 'nabla', nabla)
    call write_dset(hg, 'epsilon', epsilon)
    call write_dset(hg, 'epsilon_rho', epsilon_rho)
    call write_dset(hg, 'epsilon_T', epsilon_T)
    call write_dset(hg, 'kappa', kappa)
    call write_dset(hg, 'kappa_rho', kappa_rho)
    call write_dset(hg, 'kappa_T', kappa_T)
    call write_dset(hg, 'Omega_rot', Omega_rot)

    call hg%final()

    ! Finish

    return

  end subroutine write_gsm_data

end module gyre_gsm_file
