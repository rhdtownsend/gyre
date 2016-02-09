! Module   : gyre_osc_file
! Purpose  : read OSC files
!
! Copyright 2013-2016 Rich Townsend
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
  use core_order

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

  public :: read_osc_model

  ! Procedures

contains

  subroutine read_osc_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    character(:), allocatable   :: data_format
    integer                     :: unit
    integer                     :: n
    integer                     :: iconst
    integer                     :: ivar
    integer                     :: iabund
    integer                     :: ivers
    real(WP), allocatable       :: glob(:)
    real(WP), allocatable       :: var(:,:)
    integer                     :: i
    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP)                    :: L_star
    real(WP), allocatable       :: r(:)
    real(WP), allocatable       :: m(:)
    real(WP), allocatable       :: p(:)
    real(WP), allocatable       :: rho(:) 
    real(WP), allocatable       :: T(:) 
    real(WP), allocatable       :: Gamma_1(:)
    real(WP), allocatable       :: nabla_ad(:)
    real(WP), allocatable       :: delta(:)
    real(WP), allocatable       :: As(:)
    real(WP), allocatable       :: nabla(:)
    real(WP), allocatable       :: kappa(:)
    real(WP), allocatable       :: kappa_rho(:)
    real(WP), allocatable       :: kappa_T(:)
    real(WP), allocatable       :: epsilon_(:)
    real(WP), allocatable       :: epsilon_rho(:)
    real(WP), allocatable       :: epsilon_T(:)
    real(WP), allocatable       :: Omega_rot(:)
    real(WP), allocatable       :: x(:)
    real(WP), allocatable       :: V_2(:)
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

    ! Open the OSC-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from OSC file', TRIM(ml_p%file)
100    format(A,1X,A)
    endif
          
    open(NEWUNIT=unit, FILE=ml_p%file, STATUS='OLD')

    ! Read the header

    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *)

    read(unit, *) n, iconst, ivar, iabund, ivers

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Initial points :', n
       write(OUTPUT_UNIT, 110) 'File version   :', ivers
110    format(3X,A,1X,I0)
    endif

    ! Read the data

    if (ml_p%data_format /= '') then
       data_format = ml_p%data_format
    else
       data_format = '(1P5E19.12)'
    endif

    allocate(glob(iconst))
    allocate(var(ivar+iabund,n))

    read(unit, data_format) glob

    read_loop : do i = 1,n
       read(unit, data_format) var(:,i)
    end do read_loop

    close(unit)

    var = var(:,n:1:-1)

    ! Extract structure data

    M_star = glob(1)
    R_star = glob(2)
    L_star = glob(3)

    r = var(1,:)
    m = EXP(var(2,:))

    T = var(3,:)
    P = var(4,:)
    rho = var(5,:)

    Gamma_1 = var(10,:)
    nabla_ad = var(11,:)
    delta = var(12,:)

    As = var(15,:)

    nabla = var(6,:)
    kappa = var(8,:)
    kappa_T = var(17,:)
    kappa_rho = var(18,:)
    epsilon_ = var(9,:)
    epsilon_T = var(19,:)
    epsilon_rho = var(20,:)

    Omega_rot = var(16,:)

    ! Snap grid points

    call snap_points(MAX(ml_p%dx_snap, EPSILON(0._WP)), x, m)

    ! Calculate dimensionless structure data

    x = r/R_star

    allocate(V_2(n))
    allocate(U(n))
    allocate(c_1(n))

    where (x /= 0._WP)
       V_2 = G_GRAVITY*m*rho/(P*r*x**2)
       U = 4._WP*PI*rho*r**3/m
       c_1 = (r/R_star)**3/(m/M_star)
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._WP*P(1))
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

  end subroutine read_osc_model

end module gyre_osc_file
