! Program  : gyre_contour
! Purpose  : discriminant contouring code
!
! Copyright 2014-2018 Rich Townsend
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

program gyre_contour

  ! Uses

  use core_kinds
  use core_hgroup
  use core_order
  use core_parallel
  use core_system

  use gyre_bvp
  use gyre_constants
  use gyre_context
  use gyre_contour_map
  use gyre_contour_seg
  use gyre_discrim_func
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_mode
  use gyre_mode_par
  use gyre_model
  use gyre_model_factory
  use gyre_model_par
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_osc_par
  use gyre_out_par
  use gyre_output
  use gyre_root
  use gyre_scan_par, only : scan_par_t
  use gyre_search
  use gyre_state
  use gyre_status
  use gyre_util
  use gyre_version
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable        :: filename
  integer                          :: unit
  type(model_par_t)                :: ml_p
  type(mode_par_t), allocatable    :: md_p(:)
  type(osc_par_t), allocatable     :: os_p(:)
  type(num_par_t), allocatable     :: nm_p(:)
  type(grid_par_t), allocatable    :: gr_p(:)
  type(scan_par_t), allocatable    :: sc_p_re(:)
  type(scan_par_t), allocatable    :: sc_p_im(:)
  type(out_par_t)                  :: ot_p
  character(FILENAME_LEN)          :: map_filename
  character(FILENAME_LEN)          :: re_seg_filename
  character(FILENAME_LEN)          :: im_seg_filename
  class(model_t), pointer          :: ml => null()
  type(osc_par_t)                  :: os_p_sel
  type(num_par_t)                  :: nm_p_sel
  type(grid_par_t)                 :: gr_p_sel
  type(scan_par_t), allocatable    :: sc_p_re_sel(:)
  type(scan_par_t), allocatable    :: sc_p_im_sel(:)
  type(context_t), pointer         :: cx => null()
  real(WP), allocatable            :: omega_re(:)
  real(WP), allocatable            :: omega_im(:)
  real(WP)                         :: omega_min
  real(WP)                         :: omega_max
  type(grid_t)                     :: gr
  type(nad_bvp_t), target          :: bp
  type(c_discrim_func_t)           :: df
  type(c_state_t)                  :: st
  type(c_ext_t), allocatable       :: discrim_map(:,:)
  type(contour_map_t)              :: cm_re
  type(contour_map_t)              :: cm_im
  type(contour_seg_t), allocatable :: is_re(:)
  type(contour_seg_t), allocatable :: is_im(:)
  type(hgroup_t)                   :: hg
  integer                          :: n_md_nad
  integer                          :: d_md_nad
  type(mode_t), allocatable        :: md_nad(:)

  ! Read command-line arguments

  $ASSERT(n_arg() == 1,Syntax: gyre_map <filename>)

  call get_arg(1, filename)

  ! Initialize

  call init_parallel()

  call set_log_level($str($LOG_LEVEL))

  if (check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_contour ['//VERSION//']', '=')
100  format(A)

     if (check_log_level('DEBUG')) then
        write(OUTPUT_UNIT, 110) 'Compiler         :', COMPILER_VERSION()
        write(OUTPUT_UNIT, 110) 'Compiler options :', COMPILER_OPTIONS()
110     format(A,1X,A)
     endif

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   :', OMP_SIZE_MAX
     write(OUTPUT_UNIT, 120) 'MPI Processors   :', MPI_SIZE
120  format(A,1X,I0)
     
     write(OUTPUT_UNIT, 110) 'Input filename   :', filename
     write(OUTPUT_UNIT, 110) 'GYRE_DIR         :', gyre_dir

     write(OUTPUT_UNIT, 100) form_header('Initialization', '=')

  endif

  ! Process arguments

  open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

  call read_constants(unit)

  call read_model_par(unit, ml_p)
  call read_mode_par(unit, md_p)
  call read_osc_par(unit, os_p)
  call read_num_par(unit, nm_p)
  call read_grid_par(unit, gr_p)
  call read_scan_par(unit, 're', sc_p_re)
  call read_scan_par(unit, 'im', sc_p_im)
  call read_out_par(unit, 'nad', ot_p)
  call read_map_par(unit, map_filename, re_seg_filename, im_seg_filename)

  $ASSERT(SIZE(md_p) == 1,Must be exactly one mode parameter)

  ! Initialize the model

  ml => model_t(ml_p)

  ! Select parameters according to tags

  call select_par(os_p, md_p(1)%tag, os_p_sel)
  call select_par(nm_p, md_p(1)%tag, nm_p_sel)
  call select_par(gr_p, md_p(1)%tag, gr_p_sel)
  call select_par(sc_p_re, md_p(1)%tag, sc_p_re_sel)
  call select_par(sc_p_im, md_p(1)%tag, sc_p_im_sel)
  
  ! Create the context

  allocate(cx, SOURCE=context_t(ml, gr_p_sel, md_p(1), os_p_sel))

  ! Set up the frequency arrays

  call build_scan(cx, md_p(1), os_p_sel, sc_p_re_sel, omega_re)
  call build_scan(cx, md_p(1), os_p_sel, sc_p_im_sel, omega_im)

  ! Create the grid

  gr = grid_t(cx, omega_re, gr_p_sel)

  ! Set frequency bounds

  if (nm_p_sel%restrict_roots) then
     omega_min = MINVAL(omega_re)
     omega_max = MAXVAL(omega_re)
  else
     omega_min = -HUGE(0._WP)
     omega_max = HUGE(0._WP)
  endif

  ! Set up the bvp

  bp = nad_bvp_t(cx, gr, md_p(1), nm_p_sel, os_p_sel)

  ! Set up the discriminant function

  st = c_state_t(omega=0._WP, omega_r=0._WP)
  df = c_discrim_func_t(bp, st, omega_min, omega_max)

  ! Evaluate the discriminant map

  call eval_map(bp, omega_re, omega_im, discrim_map)

  ! (Next steps are root rank only)

  if (MPI_RANK == 0) then

     ! Write out the discriminant map

     if (map_filename /= '') then

        hg = hgroup_t(map_filename, CREATE_FILE)

        call write_dset(hg, 'omega_re', omega_re)
        call write_dset(hg, 'omega_im', omega_im)

        call write_dset(hg, 'discrim_map_f', FRACTION(discrim_map))
        call write_dset(hg, 'discrim_map_e', EXPONENT(discrim_map))

        call hg%final()

     endif

     ! Create the contour maps

     cm_re = contour_map_t(r_ext_t(omega_re), r_ext_t(omega_im), real_part(discrim_map))
     cm_im = contour_map_t(r_ext_t(omega_re), r_ext_t(omega_im), imag_part(discrim_map))
     
     ! Write out the segments

     if (re_seg_filename /= '') then
        call write_segs(cm_re, re_seg_filename)
     endif

     if (im_seg_filename /= '') then
        call write_segs(cm_im, im_seg_filename)
     endif

     ! Search for contour intersections

     call find_isects(cm_re, cm_im, is_re, is_im)

     ! Find roots

     d_md_nad = 128
     n_md_nad = 0

     allocate(md_nad(d_md_nad))

     call find_roots(bp, md_p(1), nm_p_sel, os_p_sel, is_re, is_im, process_mode)

     ! Write the summary file

     call write_summary(md_nad(:n_md_nad), ot_p)

  end if

  ! Finish

  call final_parallel()

contains

  subroutine read_scan_par (unit, axis, sc_p)

    integer, intent(in)                        :: unit
    character(*), intent(in)                   :: axis
    type(scan_par_t), allocatable, intent(out) :: sc_p(:)

    integer                             :: n_sc_p
    integer                             :: i
    real(WP)                            :: freq_min
    real(WP)                            :: freq_max
    integer                             :: n_freq
    character(LEN(sc_p%freq_min_units)) :: freq_min_units
    character(LEN(sc_p%freq_max_units)) :: freq_max_units
    character(LEN(sc_p%freq_frame))     :: freq_frame
    character(LEN(sc_p%grid_type))      :: grid_type
    character(LEN(sc_p%grid_frame))     :: grid_frame
    character(LEN(sc_p%tag_list))       :: tag_list

    namelist /re_scan/ freq_min, freq_max, n_freq, freq_min_units, freq_max_units, &
         freq_frame, grid_type, grid_frame, tag_list

    namelist /im_scan/ freq_min, freq_max, n_freq, freq_min_units, freq_max_units, &
         freq_frame, grid_type, grid_frame, tag_list

    ! Count the number of scan namelists

    rewind(unit)

    n_sc_p = 0

    count_loop : do
       select case (axis)
       case ('re')
          read(unit, NML=re_scan, END=100)
       case ('im')
          read(unit, NML=im_scan, END=100)
       case default
          $ABORT(Invalid axis)
       end select
       n_sc_p = n_sc_p + 1
    end do count_loop

100 continue

    ! Read scan parameters

    rewind(unit)

    allocate(sc_p(n_sc_p))

    read_loop : do i = 1, n_sc_p

       freq_min = 1._WP
       freq_max = 10._WP
       n_freq = 10
          
       freq_max_units = 'NONE'
       freq_min_units = 'NONE'
       freq_frame = 'INERTIAL'

       grid_type = 'LINEAR'
       grid_frame = 'INERTIAL'

       tag_list = ''

       select case (axis)
       case ('re')
          read(unit, NML=re_scan)
       case ('im')
          read(unit, NML=im_scan)
       case default
          $ABORT(Invalid axis)
       end select

       ! Initialize the scan_par

       sc_p(i) = scan_par_t(freq_min=freq_min, &
                            freq_max=freq_max, &
                            n_freq=n_freq, &
                            freq_min_units=freq_min_units, &
                            freq_max_units=freq_max_units, &
                            freq_frame=freq_frame, &
                            grid_type=grid_type, &
                            grid_frame=grid_frame, &
                            tag_list=tag_list)

    end do read_loop

    ! Finish

    return

  end subroutine read_scan_par

  !****

  subroutine read_map_par (unit, map_filename, re_seg_filename, im_seg_filename)

    integer, intent(in)       :: unit
    character(*), intent(out) :: map_filename
    character(*), intent(out) :: re_seg_filename
    character(*), intent(out) :: im_seg_filename

    namelist /map_output/ map_filename, re_seg_filename, im_seg_filename

    ! Read output parameters

    rewind(unit)

    read(unit, NML=map_output)

    ! Finish
    
    return

  end subroutine read_map_par
          
  !****

  subroutine eval_map (bp, omega_re, omega_im, discrim_map)

    class(c_bvp_t), intent(inout)           :: bp
    real(WP), intent(in)                    :: omega_re(:)
    real(WP), intent(in)                    :: omega_im(:)
    type(c_ext_t), allocatable, intent(out) :: discrim_map(:,:)

    integer                  :: n_omega_re
    integer                  :: n_omega_im
    integer, allocatable     :: k_part(:)
    integer                  :: n_percent
    integer                  :: k
    integer                  :: i(2)
    complex(WP)              :: omega
    type(c_ext_t)            :: discrim
    integer                  :: i_percent
    integer                  :: status
    complex(WP), allocatable :: discrim_map_f(:,:)
    integer, allocatable     :: discrim_map_e(:,:)
    $if ($MPI)
    integer                  :: p
    $endif

    ! Map the discriminant

    n_omega_re = SIZE(omega_re)
    n_omega_im = SIZE(omega_im)

    allocate(discrim_map_f(n_omega_re,n_omega_im))
    allocate(discrim_map_e(n_omega_re,n_omega_im))

    allocate(k_part(MPI_SIZE+1))

    call partition_tasks(n_omega_re*n_omega_im, 1, k_part)

    n_percent = 0

    do k = k_part(MPI_RANK+1), k_part(MPI_RANK+2)-1

       i = index_nd(k, [n_omega_re,n_omega_im])

       omega = CMPLX(omega_re(i(1)), omega_im(i(2)), KIND=WP)

       if (MPI_RANK == 0) then
          i_percent = FLOOR(100._WP*REAL(k-k_part(MPI_RANK+1))/REAL(k_part(MPI_RANK+2)-k_part(MPI_RANK+1)-1))
          if (i_percent > n_percent) then
             print *,'Percent complete: ', i_percent
             n_percent = i_percent
          end if
       endif
     
       call df%eval(c_ext_t(omega), discrim, status)
       $ASSERT(status == STATUS_OK,Invalid status)
       
       discrim_map_f(i(1),i(2)) = FRACTION(discrim)
       discrim_map_e(i(1),i(2)) = EXPONENT(discrim)

    end do

    $if ($MPI)

    do p = 1,MPI_SIZE
       call bcast_seq(discrim_map_f, k_part(p), k_part(p+1)-1, p-1)
       call bcast_seq(discrim_map_e, k_part(p), k_part(p+1)-1, p-1)
    end do
    
    $endif

    discrim_map = scale(c_ext_t(discrim_map_f), discrim_map_e)

    ! Finish

    return

  end subroutine eval_map

  !****

  subroutine find_isects (cm_re, cm_im, is_re, is_im)

    type(contour_map_t), intent(in)               :: cm_re
    type(contour_map_t), intent(in)               :: cm_im
    type(contour_seg_t), allocatable, intent(out) :: is_re(:)
    type(contour_seg_t), allocatable, intent(out) :: is_im(:)

    integer                          :: i_x
    integer                          :: i_y
    type(contour_seg_t), allocatable :: cs_re(:)
    type(contour_seg_t), allocatable :: cs_im(:)
    integer                          :: j_re
    integer                          :: j_im
    type(r_ext_t)                    :: M(2,2)
    type(r_ext_t)                    :: R(2)
    type(r_ext_t)                    :: D
    type(r_ext_t)                    :: w_im
    type(r_ext_t)                    :: w_re

    $CHECK_BOUNDS(cm_im%n_x,cm_re%n_x)
    $CHECK_BOUNDS(cm_im%n_y,cm_re%n_y)

    ! Loop through the cells of the real/imaginary contour maps,
    ! looking for intersections between the respective contour
    ! segments

    allocate(is_re(0))
    allocate(is_im(0))

    do i_x = 1, cm_re%n_x-1
       do i_y = 1, cm_im%n_y-1

          ! Get the contour segments for the current cell

          call cm_re%get_segs(i_x, i_y, cs_re)
          call cm_im%get_segs(i_x, i_y, cs_im)

          ! Look for pairwise intersections

          do j_re = 1, SIZE(cs_re)
             do j_im = 1, SIZE(cs_im)

                ! Solve for the intersection between the two line segments

                M(1,1) = cs_re(j_re)%x(2) - cs_re(j_re)%x(1)
                M(1,2) = cs_im(j_im)%x(1) - cs_im(j_im)%x(2)
                
                M(2,1) = cs_re(j_re)%y(2) - cs_re(j_re)%y(1)
                M(2,2) = cs_im(j_im)%y(1) - cs_im(j_im)%y(2)

                R(1) = -cs_re(j_re)%x(1) + cs_im(j_im)%x(1)
                R(2) = -cs_re(j_re)%y(1) + cs_im(j_im)%y(1)

                D = M(1,1)*M(2,2) - M(1,2)*M(2,1)

                if (D /= 0._WP) then

                   ! Lines intersect

                   w_re = (M(2,2)*R(1) - M(1,2)*R(2))/D
                   w_im = (M(1,1)*R(2) - M(2,1)*R(1))/D

                   if (w_re >= 0._WP .AND. w_re <= 1._WP .AND. &
                       w_im >= 0._WP .AND. w_im <= 1._WP) then

                      ! Lines intersect within the cell; save the
                      ! intersecting segments

                      is_re = [is_re,cs_re(j_re)]
                      is_im = [is_im,cs_im(j_im)]

                   end if

                endif
                         
             end do
          end do

       end do
    end do

    ! Finish
    
    return

  end subroutine find_isects

  !****

  subroutine find_roots (bp, md_p, nm_p, os_p, is_re, is_im, process_mode)

    type(nad_bvp_t), target, intent(inout)       :: bp
    type(mode_par_t), intent(in)                 :: md_p
    type(num_par_t), intent(in)                  :: nm_p
    type(osc_par_t), intent(in)                  :: os_p
    type(contour_seg_t), allocatable, intent(in) :: is_re(:)
    type(contour_seg_t), allocatable, intent(in) :: is_im(:)
    interface
       subroutine process_mode (md, n_iter, chi)
         use core_kinds
         use gyre_ext
         use gyre_mode
         type(mode_t), intent(in)  :: md
         integer, intent(in)       :: n_iter
         type(r_ext_t), intent(in) :: chi
       end subroutine process_mode
    end interface

    integer       :: j
    type(c_ext_t) :: omega_a
    type(c_ext_t) :: omega_b
    type(c_ext_t) :: discrim_a
    type(c_ext_t) :: discrim_b
    integer       :: status
    type(c_ext_t) :: omega_root
    integer       :: n_iter
    type(wave_t)  :: wv
    type(mode_t)  :: md
    type(r_ext_t) :: chi
    
    $CHECK_BOUNDS(SIZE(is_re),SIZE(is_im))

    ! Loop over pairs of intersecting segments

    intseg_loop : do j = 1, SIZE(is_re)

       n_iter = 0

       ! Set up initial guesses

       omega_a = omega_initial(is_re(j), .TRUE.)
       omega_b = omega_initial(is_im(j), .FALSE.)

       call df%eval(omega_a, discrim_a, status)
       if (status /= STATUS_OK) then
          call report_status(status, 'initial guess (a)')
          cycle intseg_loop
       endif

       call df%eval(omega_b, discrim_b, status)
       if (status /= STATUS_OK) then
          call report_status(status, 'initial guess (b)')
          cycle intseg_loop
       endif

       ! Find the discriminant root

       call solve(df, nm_p, omega_a, omega_b, r_ext_t(0._WP), omega_root, status, &
                  n_iter=n_iter, n_iter_max=nm_p%n_iter_max, f_cx_a=discrim_a, f_cx_b=discrim_b)
       if (status /= STATUS_OK) then
          call report_status(status, 'solve')
          cycle intseg_loop
       endif

       ! Construct the mode_t

       st = c_state_t(omega=cmplx(omega_root), omega_r=0._WP)

       wv = wave_t(bp, st)
       md = mode_t(wv, j)

       ! Process it

       chi = abs(md%discrim)/max(abs(discrim_a), abs(discrim_b))
       
       call process_mode(md, n_iter, chi)

    end do intseg_loop

    ! Finish

    return

  end subroutine find_roots

  !****

  subroutine report_status (status, stage_str)

    integer, intent(in)      :: status
    character(*), intent(in) :: stage_str

    ! Report the status

    if (check_log_level('WARN')) then

       write(OUTPUT_UNIT, 100) 'Failed during ', stage_str, ' : ', status_str(status)
100    format(4A)

    endif
      
    ! Finish

    return

  end subroutine report_status

  !****

  subroutine process_mode (md, n_iter, chi)

    type(mode_t), intent(in)  :: md
    integer, intent(in)       :: n_iter
    type(r_ext_t), intent(in) :: chi

    ! Process the non-adiabatic mode

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) md%l, md%m, md%n_pg, md%n_p, md%n_g, &
            md%omega, real(chi), n_iter
100    format(1X,I3,1X,I4,1X,I7,1X,I6,1X,I6,1X,E15.8,1X,E15.8,1X,E10.4,1X,I6)
    endif

    ! Store it

    n_md_nad = n_md_nad + 1

    if (n_md_nad > d_md_nad) then
       d_md_nad = 2*d_md_nad
       call reallocate(md_nad, [d_md_nad])
    endif

    md_nad(n_md_nad) = md

    ! Write it

    call write_mode(md_nad(n_md_nad), ot_p)

    ! If necessary, prune it

    if (ot_p%prune_modes) call md_nad(n_md_nad)%prune()

    ! Finish

    return

  end subroutine process_mode

  !****

  function omega_initial (cs, real_seg) result (omega)

    type(contour_seg_t), intent(in) :: cs
    logical, intent(in)             :: real_seg
    type(c_ext_t)                   :: omega

    type(c_ext_t) :: omega_a
    type(c_ext_t) :: omega_b
    type(c_ext_t) :: discrim_a
    type(c_ext_t) :: discrim_b
    integer       :: status
    type(r_ext_t) :: w

    ! Evaluate the discriminant at the segment endpoints

    omega_a = c_ext_t(cs%x(1), cs%y(1))
    omega_b = c_ext_t(cs%x(2), cs%y(2))

    call df%eval(omega_a, discrim_a, status)
    $ASSERT(status == STATUS_OK,Invalid status)

    call df%eval(omega_b, discrim_b, status)
    $ASSERT(status == STATUS_OK,Invalid status)

    ! Look for the point on the segment where the real/imaginary part of the discriminant changes sign

    if (real_seg) then
       w = -imag_part(discrim_a)/(imag_part(discrim_b) - imag_part(discrim_a))
    else
       w = -real_part(discrim_a)/(real_part(discrim_b) - real_part(discrim_a))
    endif
    
    omega = c_ext_t((1._WP-w)*cs%x(1) + w*cs%x(2), (1._WP-w)*cs%y(1) + w*cs%y(2))

    ! Finish
      
    return

  end function omega_initial

  !****

  subroutine write_segs (cm, filename)

    type(contour_map_t), intent(in) :: cm
    character(*), intent(in)        :: filename

    integer                          :: unit
    integer                          :: i_x
    integer                          :: i_y
    type(contour_seg_t), allocatable :: cs(:)
    integer                          :: j

    ! Open the file

    open(NEWUNIT=unit, FILE=filename, STATUS='REPLACE')

    ! Loop through the map, writing out segment info

    do i_x = 1, cm%n_x-1
       do i_y = 1, cm%n_y-1
          call cm%get_segs(i_x, i_y, cs)
          do j = 1, SIZE(cs)
             write(unit, 100) real(cs(j)%x), real(cs(j)%y)
100          format(4(1X,E16.8))
          end do
       end do
    end do

    close(unit)

    ! Finish
    
    return

  end subroutine write_segs

end program gyre_contour
