! Module   : gyre_lib
! Purpose  : library interface for use in MESA
!
! Copyright 2013-2020 Rich Townsend & The GYRE Team
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
  use core_memory

  use gyre_ad_bvp
  use gyre_bvp
  use gyre_constants, gyre_set_constant => set_constant
  use gyre_context
  use gyre_evol_model
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_math
  use gyre_mesa_file
  use gyre_mode
  use gyre_mode_par
  use gyre_model
  use gyre_model_par
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_osc_par
  use gyre_rad_bvp
  use gyre_rot_par
  use gyre_scan
  use gyre_scan_par
  use gyre_search
  use gyre_util
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  type(model_par_t), save             :: ml_p_m
  type(mode_par_t), allocatable, save :: md_p_m(:)
  type(osc_par_t), allocatable, save  :: os_p_m(:)
  type(rot_par_t), allocatable, save  :: rt_p_m(:)
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
  public :: grid_t
  public :: model_par_t
  public :: mode_par_t
  public :: osc_par_t
  public :: rot_par_t
  public :: num_par_t
  public :: grid_par_t
  public :: scan_par_t
  public :: gyre_init
  public :: gyre_final
  public :: gyre_get_par
  public :: gyre_set_par
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

    call init_math()

    call set_log_level('WARN')

    ! Read the namelist file

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    call read_model_par(unit, ml_p_m)
    call read_mode_par(unit, md_p_m)
    call read_osc_par(unit, os_p_m)
    call read_rot_par(unit, rt_p_m)
    call read_num_par(unit, nm_p_m)
    call read_grid_par(unit, gr_p_m)
    call read_scan_par(unit, sc_p_m)

    close(unit)

    ! Finish

    return

  end subroutine gyre_init

  !****

  subroutine gyre_final()

    ! Finalize

    if (ASSOCIATED(ml_m)) deallocate(ml_m)

    if (ALLOCATED(md_p_m)) deallocate(md_p_m)
    if (ALLOCATED(os_p_m)) deallocate(os_p_m)
    if (ALLOCATED(rt_p_m)) deallocate(rt_p_m)
    if (ALLOCATED(nm_p_m)) deallocate(nm_p_m)
    if (ALLOCATED(gr_p_m)) deallocate(gr_p_m)
    if (ALLOCATED(sc_p_m)) deallocate(sc_p_m)

    call final_parallel()

    ! Finish

    return

  end subroutine gyre_final

  !****

  subroutine gyre_get_par (model_par, mode_par, osc_par, rot_par, &
       num_par, grid_par, scan_par)

     type(model_par_t), intent(out), optional             :: model_par
     type(mode_par_t), allocatable, intent(out), optional :: mode_par(:)
     type(osc_par_t), allocatable, intent(out), optional  :: osc_par(:)
     type(rot_par_t), allocatable, intent(out), optional  :: rot_par(:)
     type(num_par_t), allocatable, intent(out), optional  :: num_par(:)
     type(grid_par_t), allocatable, intent(out), optional :: grid_par(:)
     type(scan_par_t), allocatable, intent(out), optional :: scan_par(:)

     ! Get parameter set(s)

     if (PRESENT(model_par)) model_par = ml_p_m
     if (PRESENT(mode_par)) mode_par = md_p_m
     if (PRESENT(model_par)) osc_par = os_p_m
     if (PRESENT(rot_par)) rot_par = rt_p_m
     if (PRESENT(num_par)) num_par = nm_p_m
     if (PRESENT(grid_par)) grid_par = gr_p_m
     if (PRESENT(scan_par)) scan_par = sc_p_m

     ! Finish

     return

  end subroutine gyre_get_par

  !****

  subroutine gyre_set_par (model_par, mode_par, osc_par, rot_par, &
       num_par, grid_par, scan_par)

     type(model_par_t), intent(in), optional :: model_par
     type(mode_par_t), intent(in), optional  :: mode_par(:)
     type(osc_par_t), intent(in), optional   :: osc_par(:)
     type(rot_par_t), intent(in), optional   :: rot_par(:)
     type(num_par_t), intent(in), optional   :: num_par(:)
     type(grid_par_t), intent(in), optional  :: grid_par(:)
     type(scan_par_t), intent(in), optional  :: scan_par(:)

     ! Get parameter set(s)

     if (PRESENT(model_par)) ml_p_m = model_par
     if (PRESENT(mode_par)) md_p_m = mode_par
     if (PRESENT(model_par)) os_p_m = osc_par
     if (PRESENT(rot_par)) rt_p_m = rot_par
     if (PRESENT(num_par)) nm_p_m = num_par
     if (PRESENT(grid_par)) gr_p_m = grid_par
     if (PRESENT(scan_par)) sc_p_m = scan_par

     ! Finish

     return

  end subroutine gyre_set_par

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

    type(context_t), pointer      :: cx => null()
    integer                       :: i
    type(osc_par_t)               :: os_p_sel
    type(rot_par_t)               :: rt_p_sel
    type(num_par_t)               :: nm_p_sel
    type(grid_par_t)              :: gr_p_sel
    type(scan_par_t), allocatable :: sc_p_sel(:)
    type(grid_t)                  :: gr
    real(WP), allocatable         :: omega_re(:)
    real(WP), allocatable         :: omega_im(:)
    real(WP)                      :: omega_min
    real(WP)                      :: omega_max
    class(r_bvp_t), allocatable   :: bp_ad
    class(c_bvp_t), allocatable   :: bp_nad
    integer                       :: n_ad
    integer                       :: d_ad
    complex(WP), allocatable      :: omega_ad(:)
    integer, allocatable          :: id_ad(:)

    $ASSERT(ASSOCIATED(ml_m),No model provided)

    ! Check that GYRE_DIR is set

    $ASSERT(GYRE_DIR /= '',The GYRE_DIR environment variable is not set)

    ! Allocate the context (will be initialized later on)

    allocate(cx)

    ! Loop through modepars

    md_p_loop : do i = 1, SIZE(md_p_m)

       if (md_p_m(i)%l /= l) cycle md_p_loop

       ! Select parameters according to tags

       call select_par(os_p_m, md_p_m(i)%tag, os_p_sel)
       call select_par(rt_p_m, md_p_m(i)%tag, rt_p_sel)
       call select_par(nm_p_m, md_p_m(i)%tag, nm_p_sel)
       call select_par(gr_p_m, md_p_m(i)%tag, gr_p_sel)
       call select_par(sc_p_m, md_p_m(i)%tag, sc_p_sel)

       ! Set up the context

       cx = context_t(ml_m, gr_p_sel, md_p_m(i), os_p_sel, rt_p_sel)

       ! Set up the frequency array
          
       call build_scan(cx, md_p_m(i), os_p_sel, sc_p_sel, 'REAL', omega_re)
       call build_scan(cx, md_p_m(i), os_p_sel, sc_p_sel, 'IMAG', omega_im)

       ! Create the grid

       gr = grid_t(cx, omega_re, gr_p_sel, os_p_sel)

       ! Set frequency bounds and perform checks
       
       if (nm_p_sel%restrict_roots) then
          omega_min = MINVAL(omega_re)
          omega_max = MAXVAL(omega_re)
       else
          omega_min = -HUGE(0._WP)
          omega_max = HUGE(0._WP)
       endif

       call check_scan(cx, gr, omega_re, md_p_m(i), os_p_sel)

       ! Find adiabatic modes

       if (os_p_sel%adiabatic) then
        
          d_ad = 128
          n_ad = 0

          allocate(omega_ad(d_ad))
          allocate(id_ad(d_ad))

          if (md_p_m(i)%l == 0 .AND. os_p_sel%reduce_order) then
             allocate(bp_ad, SOURCE=rad_bvp_t(cx, gr, md_p_m(i), nm_p_sel, os_p_sel))
          else
             allocate(bp_ad, SOURCE=ad_bvp_t(cx, gr, md_p_m(i), nm_p_sel, os_p_sel))
          endif

          select case (nm_p_sel%ad_search)

          case ('BRACKET')

             if (SIZE(omega_re) > 2) then
                call bracket_search(bp_ad, omega_re, omega_min, omega_max, nm_p_sel, process_mode_ad)
             endif
                
          case default

             $ABORT(Invalid ad_search)

          end select

          deallocate(bp_ad)

       endif

       ! Find non-adiabatic modes

       if (os_p_sel%nonadiabatic) then

          allocate(bp_nad, SOURCE=nad_bvp_t(cx, gr, md_p_m(i), nm_p_sel, os_p_sel))

          select case (nm_p_sel%nad_search)

          case ('AD')

             $ASSERT(os_p_sel%adiabatic,No adiabatic modes to start from)
             call prox_search(bp_nad, omega_ad(:n_ad), id_ad(:n_ad), omega_min, omega_max, nm_p_sel, process_mode_nad)

          case ('MINMOD')

             if (SIZE(omega_re) > 2) then
                call minmod_search(bp_nad, omega_re, omega_min, omega_max, nm_p_sel, process_mode_nad)
             endif

          case ('CONTOUR')

             if (SIZE(omega_re) > 2 .AND. SIZE(omega_im) > 2) then
                call contour_search(bp_nad, omega_re, omega_im, omega_min, omega_max, i, nm_p_sel, process_mode_nad)
             endif

          case default

             $ABORT(Invalid nad_start)

          end select

          deallocate(bp_nad)

       endif

       ! Deallocate adiabatic data

       if (os_p_sel%adiabatic) then
          deallocate(omega_ad)
          deallocate(id_ad)
       endif

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

         n_ad = n_ad + 1

         if (n_ad > d_ad) then
            d_ad = 2*d_ad
            call reallocate(omega_ad, [d_ad])
            call reallocate(id_ad, [d_ad])
         endif

         omega_ad(n_ad) = md%omega
         id_ad(n_ad) = md%id

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

end module gyre_lib
