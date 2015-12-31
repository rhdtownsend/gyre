! Module   : gyre_mesa_file
! Purpose  : read MESA files
!
! Copyright 2013-2015 Rich Townsend
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
    real(WP), allocatable       :: m(:)
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
    real(WP), allocatable       :: V_2(:)
    real(WP), allocatable       :: As(:)
    real(WP), allocatable       :: U(:)
    real(WP), allocatable       :: c_1(:)
    type(evol_model_t), pointer :: em

    ! Read data from the MESA-format file

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Reading from MESA file', TRIM(ml_p%file)
100    format(A,1X,A)
    endif

    call read_mesa_data(ml_p%file, M_star, R_star, L_star, r, m, p, rho, T, &
                        N2, Gamma_1, nabla_ad, delta, nabla,  &
                        kappa, kappa_rho, kappa_T, &
                        epsilon, epsilon_rho, epsilon_T, &
                        Omega_rot)

    ! Calculate dimensionless structure data

    n = SIZE(r)

    allocate(x(n))

    allocate(V_2(n))
    allocate(As(n))
    allocate(U(n))
    allocate(c_1(n))

    x = r/R_star

    where (x /= 0._WP)
       V_2 = G_GRAVITY*m*rho/(p*r*x**2)
       As = r**3*N2/(G_GRAVITY*m)
       U = 4._WP*PI*rho*r**3/m
       c_1 = (r/R_star)**3/(m/M_star)
    elsewhere
       V_2 = 4._WP*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._WP*p(1))
       As = 0._WP
       U = 3._WP
       c_1 = 3._WP*(M_star/R_star**3)/(4._WP*PI*rho)
    end where

    if (ml_p%uniform_rot) then
       Omega_rot = ml_p%Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
    else
       Omega_rot = Omega_rot*SQRT(R_star**3/(G_GRAVITY*M_star))
    endif

    ! Initialize the model

    allocate(em, SOURCE=evol_model_t(x, M_star, R_star, L_star, ml_p))

    call em%set_V_2(V_2)
    call em%set_As(As)
    call em%set_U(U)
    call em%set_c_1(c_1)

    call em%set_Gamma_1(Gamma_1)
    call em%set_delta(delta)
    call em%set_nabla_ad(nabla_ad)

    call em%set_Omega_rot(Omega_rot)

    ! Return a pointer to the model

    ml => em

    ! Finish

    return

  end subroutine read_mesa_model

!****

  subroutine read_mesa_data (file, M_star, R_star, L_star, r, m, p, rho, T, &
                             N2, Gamma_1, nabla_ad, delta, nabla,  &
                             kappa, kappa_rho, kappa_T, &
                             epsilon, epsilon_rho, epsilon_T, &
                             Omega_rot)

    character(*), intent(in)           :: file
    real(WP), intent(out)              :: M_star
    real(WP), intent(out)              :: R_star
    real(WP), intent(out)              :: L_star
    real(WP), allocatable, intent(out) :: r(:)
    real(WP), allocatable, intent(out) :: m(:)
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

    integer  :: unit
    integer  :: n
    integer  :: n_cols
    
    ! Read data from the MESA-format file

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    ! Read the header

    read(unit, *) n, M_star, R_star, L_star, n_cols

    if(n_cols == 1) then

       ! Old variant (n_cols not specified)

       n_cols = 18

       backspace(unit)

       if(check_log_level('INFO')) then
          write(OUTPUT_UNIT, 100) 'Detected old-variant file'
100       format(3X,A)
       endif

       call read_mesa_data_old_()

    else

       ! New variant (n_cols specified)

       if(check_log_level('INFO')) then
          write(OUTPUT_UNIT, 100) 'Detected new-variant file'
       endif

       call read_mesa_data_new_()

    endif

    close(unit)

    ! Finish

    return

  contains

    subroutine read_mesa_data_old_ ()

      real(WP), allocatable :: var(:,:)
      integer               :: k
      integer               :: k_chk
      integer, allocatable  :: ind(:)

      ! Read data from the old-variant file

      allocate(var(18,n))

      read_loop : do k = 1,n
         read(unit, *) k_chk, var(:,k)
         $ASSERT(k == k_chk,Index mismatch)
      end do read_loop

      close(unit)

      ind = unique_indices(var(1,:))

      if (SIZE(ind) < n) then

         if(check_log_level('WARN')) then
            write(OUTPUT_UNIT, 110) 'WARNING: Duplicate x-point(s) found, using innermost value(s)'
110         format('!!',1X,A)
         endif

         n = SIZE(var, 2)

      endif
       
      var = var(:,ind)

      r = var(1,:)
      m = var(2,:)/(1._WP+var(2,:))*M_star
      p = var(4,:)
      T = var(5,:)
      rho = var(6,:)
      nabla = var(7,:)
      N2 = var(8,:)
      Gamma_1 = var(12,:)*var(10,:)/var(9,:)
      delta = var(11,:)/var(12,:)
      kappa = var(13,:)
      kappa_T = var(14,:)
      kappa_rho = var(15,:)
      epsilon = var(16,:)
      epsilon_T = var(17,:)
      epsilon_rho = var(18,:)

      nabla_ad = p*delta/(rho*T*var(10,:))

      allocate(Omega_rot(n))
      Omega_rot = 0._WP

      ! Decide whether epsilon_T and epsilon_rho need rescaling

      k = MAXLOC(ABS(epsilon_T), DIM=1)

      if(ABS(epsilon_T(k)) < 1E-3*ABS(epsilon(k))) then

         epsilon_T = epsilon_T*epsilon
         epsilon_rho = epsilon_rho*epsilon

         if(check_log_level('INFO')) then
            write(OUTPUT_UNIT, 120) 'Rescaled epsilon derivatives'
120         format(3X,A)
         endif

      endif

      ! Finish

      return

    end subroutine read_mesa_data_old_

    subroutine read_mesa_data_new_ ()

      real(WP), allocatable :: var(:,:)
      integer               :: k
      integer               :: k_chk
      integer, allocatable  :: ind(:)

      $ASSERT(n_cols >= 18,Too few columns)

      ! Read data from the new-variant file

      allocate(var(n_cols-1,n))

      read_loop : do k = 1,n
         read(unit, *) k_chk, var(:,k)
         $ASSERT(k == k_chk,Index mismatch)
      end do read_loop

      close(unit)

      ind = unique_indices(var(1,:))

      if (SIZE(ind) < n) then

         if(check_log_level('WARN')) then
            write(OUTPUT_UNIT, 110) 'WARNING: Duplicate x-point(s) found, using innermost value(s)'
110         format('!!',1X,A)
         endif

         n = SIZE(var, 2)

      endif
       
      var = var(:,ind)

      r = var(1,:)
      m = var(2,:)/(1._WP+var(2,:))*M_star
      p = var(4,:)
      T = var(5,:)
      rho = var(6,:)
      nabla = var(7,:)
      N2 = var(8,:)
      Gamma_1 = var(9,:)
      nabla_ad = var(10,:)
      delta = var(11,:)
      kappa = var(12,:)
      kappa_T = var(13,:)
      kappa_rho = var(14,:)
      epsilon = var(15,:)
      epsilon_T = var(16,:)
      epsilon_rho = var(17,:)
      Omega_rot = var(18,:)

      ! Finish

      return

    end subroutine read_mesa_data_new_

  end subroutine read_mesa_data

end module gyre_mesa_file
