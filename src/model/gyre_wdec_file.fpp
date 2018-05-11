! Module   : gyre_wdec_file
! Purpose  : read WDEC files
!
! Copyright 2018 Rich Townsend
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

module gyre_wdec_file

  ! Uses

  use core_kinds

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

  public :: read_wdec_model

  ! Procedures

contains

  subroutine read_wdec_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    integer                     :: unit
    integer                     :: n
    real(WP), allocatable       :: point_data(:,:)
    integer                     :: i
    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP)                    :: L_star
    real(WP), allocatable       :: r(:)
    real(WP), allocatable       :: m(:)
    real(WP), allocatable       :: L_r(:)
    real(WP), allocatable       :: p(:)
    real(WP), allocatable       :: rho(:) 
    real(WP), allocatable       :: T(:) 
    real(WP), allocatable       :: Gamma_1(:)
    real(WP), allocatable       :: nabla_ad(:)
    real(WP), allocatable       :: delta(:)
    real(WP), allocatable       :: nabla(:)
    real(WP), allocatable       :: B(:)
    real(WP), allocatable       :: kap(:)
    real(WP), allocatable       :: kap_rho(:)
    real(WP), allocatable       :: kap_T(:)
    real(WP), allocatable       :: eps(:)
    real(WP), allocatable       :: eps_rho(:)
    real(WP), allocatable       :: eps_T(:)
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
    real(WP), allocatable       :: Omega_rot(:)
    type(evol_model_t), pointer :: em

    ! Open the WDEC-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from WDEC file'
100    format(A)
       write(OUTPUT_UNIT, 110) 'File name', TRIM(ml_p%file)
110    format(3X,A,1X,A)
    endif
          
    open(NEWUNIT=unit, FILE=ml_p%file, STATUS='OLD')

    ! Read the header

    read(unit, *)

    read(unit, *) n

    ! Read the data

    allocate(point_data(20,n))

    read_loop : do i = 1, 20
       read(unit, 120) point_data(i,:)
120 format(1P,4E22.15)
    end do read_loop

    close(unit)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 130) 'Read', n, 'points'
130    format(3X,A,1X,I0,1X,A)
    endif
    
    ! Extract structure data

    R_star = point_data(1,n)
    M_star = point_data(2,n)
    L_star = point_data(3,n)

    r = point_data(1,:)
    m = point_data(2,:)
    L_r = point_data(3,:)

    T = point_data(4,:)
    P = point_data(6,:)
    rho = point_data(5,:)

    Gamma_1 = point_data(9,:)/(1._WP - point_data(10,:)*point_data(16,:))
    nabla_ad = point_data(16,:)
    delta = point_data(10,:)/point_data(9,:)

    nabla = point_data(15,:)
    B = point_data(19,:)

    kap = point_data(18,:)
    kap_rho = point_data(13,:)
    kap_T = point_data(14,:)

    eps = point_data(7,:)
    eps_rho = point_data(12,:)
    eps_T = point_data(11,:)

    ! Snap grid points

    x = r/R_star

    call snap_points(MAX(ml_p%dx_snap, EPSILON(0._WP)), x, m)

    ! Calculate dimensionless structure data

    allocate(V_2(n))
    allocate(U(n))
    allocate(c_1(n))
    allocate(c_lum(n))

    where (x /= 0._WP)
       V_2 = G_GRAVITY*m*rho/(P*r*x**2)
       U = 4._WP*PI*rho*r**3/m
       c_1 = (r/R_star)**3/(m/M_star)
       c_lum = (L_r/L_star)/x**3
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._WP*P(1))
       U = 3._WP
       c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
       c_lum = 4._WP*PI*rho(1)*eps(1)*R_star**3/L_star
    end where

    As = V_2*x**2*delta*(nabla_ad - nabla + B)

    beta_rad = A_RADIATION*T**4/(3._WP*P)

    c_P = P*delta/(rho*T*nabla_ad)

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

  end subroutine read_wdec_model

end module gyre_wdec_file
