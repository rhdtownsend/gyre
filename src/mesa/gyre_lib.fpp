! Module   : gyre_lib
! Purpose  : library interface for use in MESA
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

module gyre_lib

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_ad_bvp
  use gyre_bvp
  use gyre_constants, gyre_set_constant => set_constant
  use gyre_context
  use gyre_evol_model
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_mesa_file
  use gyre_mode
  use gyre_mode_par
  use gyre_model
  use gyre_model_par
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_osc_par
  use gyre_rad_bvp
  use gyre_scan_par
  use gyre_search
  use gyre_tide
  use gyre_tide_par
  use gyre_util
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  type(model_par_t), save             :: ml_p_m
  type(mode_par_t), allocatable, save :: md_p_m(:)
  type(osc_par_t), allocatable, save  :: os_p_m(:)
  type(num_par_t), allocatable, save  :: nm_p_m(:)
  type(grid_par_t), allocatable, save :: gr_p_m(:)
  type(scan_par_t), allocatable, save :: sc_p_m(:)

  class(evol_model_t), pointer, save :: ml_m => null()

  ! Access specifiers

  private

  public :: WP

  public :: G_GRAVITY
  public :: C_LIGHT
  public :: A_RADIATION
  
  public :: M_SUN
  public :: R_SUN
  public :: L_SUN

  public :: mode_t
  public :: gyre_init
  public :: gyre_final
  public :: gyre_set_constant
  public :: gyre_read_model
  public :: gyre_set_model
  public :: gyre_get_modes

  ! Procedures

contains

  subroutine gyre_init (file)

    character(*), intent(in) :: file

    integer :: unit

    ! Initialize

    call init_parallel()

    call set_log_level('WARN')

    ! Read the namelist file

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    call read_model_par(unit, ml_p_m)
    call read_mode_par(unit, md_p_m)
    call read_osc_par(unit, os_p_m)
    call read_num_par(unit, nm_p_m)
    call read_grid_par(unit, gr_p_m)
    call read_scan_par(unit, sc_p_m)

    ! Finish

    return

  end subroutine gyre_init

  !****

  subroutine gyre_final()

    ! Finalize

    if (ASSOCIATED(ml_m)) deallocate(ml_m)

    if (ALLOCATED(md_p_m)) deallocate(md_p_m)
    if (ALLOCATED(os_p_m)) deallocate(os_p_m)
    if (ALLOCATED(nm_p_m)) deallocate(nm_p_m)
    if (ALLOCATED(gr_p_m)) deallocate(gr_p_m)
    if (ALLOCATED(sc_p_m)) deallocate(sc_p_m)

    call final_parallel()

    ! Finish

    return

  end subroutine gyre_final

  !****

  subroutine gyre_read_model (file)

    character(LEN=*), intent(in) :: file

    class(model_t), pointer :: ml

    ! Read the model

    if (ASSOCIATED(ml_m)) deallocate(ml_m)

    ml_p_m%file = file

    call read_mesa_model(ml_p_m, ml)

    select type (ml)
    class is (evol_model_t)
       ml_m => ml
    class default
       $ABORT(Invalid class)
    end select

    ! Finish

    return

  end subroutine gyre_read_model
  
  !****

  subroutine gyre_set_model (global_data, point_data, version)

    real(WP), intent(in) :: global_data(:)
    real(WP), intent(in) :: point_data(:,:)
    integer, intent(in)  :: version

    ! Initialize the model

    if (ASSOCIATED(ml_m)) deallocate(ml_m)

    call init_mesa_model(ml_p_m, global_data, point_data, version, ml_m)

    ! Finsh

    return

  end subroutine gyre_set_model

  !****

  subroutine gyre_get_modes (l, user_sub, ipar, rpar)

    integer, intent(in)     :: l
    interface
       subroutine user_sub (md, ipar, rpar, retcode)
         import mode_t
         import WP
         type(mode_t), intent(in) :: md
         integer, intent(inout)   :: ipar(:)
         real(WP), intent(inout)  :: rpar(:)
         integer, intent(out)     :: retcode
       end subroutine user_sub
    end interface
    integer, intent(inout)  :: ipar(:)
    real(WP), intent(inout) :: rpar(:)

    type(context_t), pointer      :: cx(:) => null()
    type(mode_t), allocatable     :: md_ad(:)
    integer                       :: n_md_ad
    integer                       :: d_md_ad
    integer                       :: i
    type(osc_par_t)               :: os_p_sel
    type(num_par_t)               :: nm_p_sel
    type(grid_par_t)              :: gr_p_sel
    type(scan_par_t), allocatable :: sc_p_sel(:)
    type(grid_t)                  :: gr
    real(WP), allocatable         :: omega(:)
    real(WP)                      :: omega_min
    real(WP)                      :: omega_max
    class(r_bvp_t), allocatable   :: bp_ad
    class(c_bvp_t), allocatable   :: bp_nad

    $ASSERT(ASSOCIATED(ml_m),No model provided)

    ! Allocate the contexts array (will be initialized later on)

    allocate(cx(SIZE(md_p_m)))

    ! Loop through modepars

    d_md_ad = 128
    n_md_ad = 0

    allocate(md_ad(d_md_ad))

    md_p_loop : do i = 1, SIZE(md_p_m)

       if (md_p_m(i)%l == l) then

          ! Select parameters according to tags

          call select_par(os_p_m, md_p_m(i)%tag, os_p_sel)
          call select_par(nm_p_m, md_p_m(i)%tag, nm_p_sel)
          call select_par(gr_p_m, md_p_m(i)%tag, gr_p_sel)
          call select_par(sc_p_m, md_p_m(i)%tag, sc_p_sel)

          ! Create the scaffold grid (used in setting up the frequency array)

          gr = grid_t(ml_m%grid(), gr_p_sel%x_i, gr_p_sel%x_o)

          ! Set up the frequency array

          call build_scan(ml_m, gr, md_p_m(i), os_p_sel, sc_p_sel, omega)

          call check_scan(ml_m, gr, omega, md_p_m(i), os_p_sel)

          ! Set frequency bounds for solutions

          if (nm_p_sel%restrict_roots) then
             omega_min = MINVAL(omega)
             omega_max = MAXVAL(omega)
          else
             omega_min = -HUGE(0._WP)
             omega_max = HUGE(0._WP)
          endif

          ! Create the full grid

          gr = grid_t(ml_m, omega, gr_p_sel, md_p_m(i), os_p_sel)

          ! Set up the context

          cx(i) = context_t(ml_m, gr%pt_i(), gr%pt_o(), md_p_m(i), os_p_sel)

          ! Set up the bvp's

          if (md_p_m(i)%l == 0 .AND. os_p_sel%reduce_order) then
             allocate(bp_ad, SOURCE=rad_bvp_t(cx(i), gr, md_p_m(i), nm_p_sel, os_p_sel))
          else
             allocate(bp_ad, SOURCE=ad_bvp_t(cx(i), gr, md_p_m(i), nm_p_sel, os_p_sel))
          endif

          if (os_p_sel%nonadiabatic) then
             allocate(bp_nad, SOURCE=nad_bvp_t(cx(i), gr, md_p_m(i), nm_p_sel, os_p_sel))
          endif

          ! Find modes

          if (os_p_sel%nonadiabatic) then
             n_md_ad = 0
             call scan_search(bp_ad, omega, omega_min, omega_max, process_mode_ad, nm_p_sel)
             call prox_search(bp_nad, md_ad(:n_md_ad), omega_min, omega_max, process_mode_nad, md_p_m(i), nm_p_sel, os_p_sel)
          else
             call scan_search(bp_ad, omega, omega_min, omega_max, process_mode_ad, nm_p_sel)
          endif

          ! Clean up

          deallocate(bp_ad)
          if (ALLOCATED(bp_nad)) deallocate(bp_nad)
          
       end if

       ! Loop around

    end do md_p_loop

    ! Finish

    deallocate(cx)

    return

  contains

    subroutine process_mode_ad (md, n_iter, chi)

      type(mode_t), intent(in)  :: md
      integer, intent(in)       :: n_iter
      type(r_ext_t), intent(in) :: chi

      integer :: retcode

      if (md%n_pg < md_p_m(i)%n_pg_min .OR. md%n_pg > md_p_m(i)%n_pg_max) return

      ! Store or process the adiabatic mode

      if (os_p_sel%nonadiabatic) then

         n_md_ad = n_md_ad + 1

         if (n_md_ad > d_md_ad) then
            d_md_ad = 2*d_md_ad
            call reallocate(md_ad, [d_md_ad])
         endif
       
         md_ad(n_md_ad) = md

         call md_ad(n_md_ad)%prune()

      else

         retcode = 0

         call user_sub(md, ipar, rpar, retcode)

      endif

      ! Finish

      return

    end subroutine process_mode_ad

    !****

    subroutine process_mode_nad (md, n_iter, chi)

      type(mode_t), intent(in)  :: md
      integer, intent(in)       :: n_iter
      type(r_ext_t), intent(in) :: chi

      integer :: retcode

      ! Process the non-adiabatic mode

      retcode = 0

      call user_sub(md, ipar, rpar, retcode)

      ! Finish

      return

    end subroutine process_mode_nad

  end subroutine gyre_get_modes

  ! !****

  ! subroutine gyre_get_tide (R_a, eps_T, Omega_orb, l_max, k_max, tau, work)

  !   real(WP), intent(in)               :: R_a
  !   real(WP), intent(in)               :: eps_T
  !   real(WP), intent(in)               :: Omega_orb
  !   integer, intent(in)                :: l_max
  !   integer, intent(in)                :: k_max
  !   real(WP), allocatable, intent(out) :: tau(:)
  !   real(WP), allocatable, intent(out) :: work(:)

  !   ! Create the tide_par_t

  !   td_p = tide_par_t(R_a=R_a, &
  !                     eps_T=eps_T, &
  !                     Omega_orb, &
  !                     l_max=l_max, &
  !                     k_max=k_max)

  !   ! Initialize the net differential torque and work arrays

  !   dwrk_dx_net = 0._WP
  !   dtau_dx_net = 0._WP

  !   ! Evaluate the tide

  !   call eval_tide(XXXXXXXXXXXXXXXX)

  !   ! Finish

  !   return

  ! contains

  !   subroutine process_wave_tide (wv)

  !     type(wave_t), intent(in) : wv

  !     real(WP) :: dwrk_dx(wv%n_k)
  !     real(WP) :: dtau_dx(wv%n_k)

  !     ! Evaluate the torque and rate-of-work functions

  !     !$OMP DO
  !     k_loop : do k = 1, wv%n_k

  !        dwrk_dx(k) = G_GRAVITY*XXXX * wv%dW_dx(k)/(XXXXXX)

  !        dtau_dx(k)= G_GRAVITY*M_star_m**2/R_star_m * wv%dtau_dx_ss(k)

  !     end do k_loop

  !     ! Interpolate these functions onto the original star grid, and
  !     ! add the contributions to the net differential work and torque

  !     XXXXX

  !     !  Finish

  !     return

  !   end subroutine process_wave_tide

  ! end subroutine gyre_get_tide

end module gyre_lib
