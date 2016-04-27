! Module   : gyre_mesa_file
! Purpose  : read MESA files
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

module gyre_mesa_file

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

  public :: read_mesa_model
  public :: read_mesa_data

  ! Procedures

contains

  subroutine read_mesa_model (ml_p, ml)

    type(model_par_t), intent(in)        :: ml_p
    class(model_t), pointer, intent(out) :: ml

    real(WP)                    :: M_star
    real(WP)                    :: R_star
    real(WP)                    :: L_star
    real(WP), allocatable       :: r(:)
    real(WP), allocatable       :: M_r(:)
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
    real(WP), allocatable       :: eps_rho(:)
    real(WP), allocatable       :: eps_T(:)
    real(WP), allocatable       :: Omega_rot(:)
    integer                     :: n
    real(WP), allocatable       :: x(:)
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
    real(WP), allocatable       :: kap_ad(:)
    real(WP), allocatable       :: kap_S(:)
    type(evol_model_t), pointer :: em

    ! Read data from the MESA-format file

    call read_mesa_data(ml_p%file, M_star, R_star, L_star, r, M_r, P, rho, T, &
                        N2, Gamma_1, nabla_ad, delta, nabla,  &
                        kap, kap_rho, kap_T, eps, eps_rho, eps_T, &
                        Omega_rot)

    ! Snap grid points

    x = r/R_star

    call snap_points(MAX(ml_p%dx_snap, EPSILON(0._WP)), x, M_r)

    ! Calculate dimensionless structure data

    n = SIZE(x)

    allocate(V_2(n))
    allocate(As(n))
    allocate(U(n))
    allocate(c_1(n))

    where (x /= 0._WP)
       V_2 = G_GRAVITY*M_r*rho/(P*r*x**2)
       As = r**3*N2/(G_GRAVITY*M_r)
       U = 4._WP*PI*rho*r**3/M_r
       c_1 = (r/R_star)**3/(M_r/M_star)
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._WP*P(1))
       As = 0._WP
       U = 3._WP
       c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
    end where

    beta_rad = A_RADIATION*T**4/(3._WP*P)

    c_P = P*delta/(rho*T*nabla_ad)

    kap_ad = nabla_ad*kap_T + kap_rho/Gamma_1
    kap_S = kap_T - delta*kap_rho

    c_rad = 16._WP*PI*A_RADIATION*C_LIGHT*T**4*R_star*nabla*V_2/(3._WP*kap*rho*L_star)
    c_thm = 4._WP*PI*rho*T*c_P*SQRT(G_GRAVITY*M_star/R_star**3)*R_star**3/L_star
    c_dif = (kap_ad-4._WP*nabla_ad)*V_2*x**2*nabla + V_2*x**2*nabla_ad

    c_eps_ad = 4._WP*PI*rho*(nabla_ad*eps_T + eps_rho/Gamma_1)*R_star**3/L_star
    c_eps_S = 4._WP*PI*rho*(eps_T - delta*eps_rho)*R_star**3/L_star

    if (ml_p%uniform_rot) then
       Omega_rot = uniform_Omega_rot(ml_p, M_star, R_star)
    else
       Omega_rot = Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
    endif

    ! Initialize the evol_model_t

    allocate(em, SOURCE=evol_model_t(x, M_star, R_star, L_star, .TRUE., ml_p))

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
    call em%set_kap_ad(kap_ad)
    call em%set_kap_S(kap_S)

    call em%set_Omega_rot(Omega_rot)

    ! Return a pointer

    ml => em

    ! Finish

    return

  end subroutine read_mesa_model

  !****

  subroutine read_mesa_data (file, M_star, R_star, L_star, r, M_r, P, rho, T, &
                             N2, Gamma_1, nabla_ad, delta, nabla,  &
                             kap, kap_rho, kap_T, eps, eps_rho, eps_T, &
                             Omega_rot)

    character(*), intent(in)           :: file
    real(WP), intent(out)              :: M_star
    real(WP), intent(out)              :: R_star
    real(WP), intent(out)              :: L_star
    real(WP), allocatable, intent(out) :: r(:)
    real(WP), allocatable, intent(out) :: M_r(:)
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
    real(WP), allocatable, intent(out) :: eps_rho(:)
    real(WP), allocatable, intent(out) :: eps_T(:)
    real(WP), allocatable, intent(out) :: Omega_rot(:)

    integer  :: unit
    integer  :: n
    integer  :: version
    
    ! Read data from the MESA-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from MESA file', TRIM(file)
100    format(A,1X,A)
    endif

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    ! Read the header and determine the version

    read(unit, *) n, M_star, R_star, L_star, version

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 110) 'Detected version', version, 'file'
110    format(3X,A,1X,F4.2,1X,A)
    endif

    ! Read the data

    select case (version)
    case (1)
       backspace(unit)
       call read_mesa_data_v0_01_()
    case (19)
       call read_mesa_data_v0_19_()
    case (100)
       call read_mesa_data_v1_00_()
    case default
       $ABORT(Unrecognized MESA file version)
    end select

    close(unit)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) 'Read', n, 'points'
120    format(3X,A,1X,I0,1X,A)
    endif
    
    ! Finish

    return

  contains

    subroutine read_mesa_data_v0_01_ ()

      real(WP), allocatable :: point_data(:,:)
      integer               :: k
      integer               :: k_chk

      ! Read data from the version-0.01 file (this is the old variant
      ! with no version number; the 1 comes from reading the first
      ! field of the following record)

      allocate(point_data(18,n))

      read_loop : do k = 1,n
         read(unit, *) k_chk, point_data(:,k)
         $ASSERT(k == k_chk,Index mismatch)
      end do read_loop

      r = point_data(1,:)
      M_r = point_data(2,:)/(1._WP+point_data(2,:))*M_star
      P = point_data(4,:)
      T = point_data(5,:)
      rho = point_data(6,:)
      nabla = point_data(7,:)
      N2 = point_data(8,:)
      Gamma_1 = point_data(12,:)*point_data(10,:)/point_data(9,:)
      delta = point_data(11,:)/point_data(12,:)
      kap = point_data(13,:)
      kap_T = point_data(14,:)
      kap_rho = point_data(15,:)
      eps = point_data(16,:)
      eps_T = point_data(17,:)
      eps_rho = point_data(18,:)

      nabla_ad = p*delta/(rho*T*point_data(10,:))

      allocate(Omega_rot(n))
      Omega_rot = 0._WP

      ! Decide whether eps_T and eps_rho need rescaling

      k = MAXLOC(ABS(eps_T), DIM=1)

      if (ABS(eps_T(k)) < 1E-3*ABS(eps(k))) then

         eps_T = eps_T*eps
         eps_rho = eps_rho*eps

         if(check_log_level('INFO')) then
            write(OUTPUT_UNIT, 120) 'Rescaled epsilon derivatives'
120         format(3X,A)
         endif

      endif

      ! Finish

      return

    end subroutine read_mesa_data_v0_01_

    subroutine read_mesa_data_v0_19_ ()

      real(WP), allocatable :: point_data(:,:)
      integer               :: k
      integer               :: k_chk

      ! Read data from the version-0.01 file (this is the old variant
      ! with no version number; the 19 comes from reading the column count)

      allocate(point_data(18,n))

      read_loop : do k = 1,n
         read(unit, *) k_chk, point_data(:,k)
         $ASSERT(k == k_chk,Index mismatch)
      end do read_loop

      r = point_data(1,:)
      M_r = point_data(2,:)/(1._WP+point_data(2,:))*M_star
      P = point_data(4,:)
      T = point_data(5,:)
      rho = point_data(6,:)
      nabla = point_data(7,:)
      N2 = point_data(8,:)
      Gamma_1 = point_data(9,:)
      nabla_ad = point_data(10,:)
      delta = point_data(11,:)
      kap = point_data(12,:)
      kap_T = point_data(13,:)
      kap_rho = point_data(14,:)
      eps = point_data(15,:)
      eps_T = point_data(16,:)
      eps_rho = point_data(17,:)
      Omega_rot = point_data(18,:)

      ! Finish

      return

    end subroutine read_mesa_data_v0_19_

    subroutine read_mesa_data_v1_00_ ()

      real(WP), allocatable :: point_data(:,:)
      integer               :: k
      integer               :: k_chk

      ! Read data from the version-1.00 file

      allocate(point_data(18,n))

      read_loop : do k = 1,n
         read(unit, *) k_chk, point_data(:,k)
         $ASSERT(k == k_chk,Index mismatch)
      end do read_loop

      r = point_data(1,:)
      M_r = point_data(2,:)
      P = point_data(4,:)
      T = point_data(5,:)
      rho = point_data(6,:)
      nabla = point_data(7,:)
      N2 = point_data(8,:)
      Gamma_1 = point_data(9,:)
      nabla_ad = point_data(10,:)
      delta = point_data(11,:)
      kap = point_data(12,:)
      kap_T = point_data(13,:)/point_data(12,:)
      kap_rho = point_data(14,:)/point_data(12,:)
      eps = point_data(15,:)
      eps_T = point_data(16,:)
      eps_rho = point_data(17,:)
      Omega_rot = point_data(18,:)

      ! Finish

      return

    end subroutine read_mesa_data_v1_00_

  end subroutine read_mesa_data

end module gyre_mesa_file
