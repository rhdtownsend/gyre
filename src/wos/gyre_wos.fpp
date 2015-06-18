! Program  : gyre_wos
! Purpose  : wave-on-string demonstration code
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

program gyre_wos

  ! Uses

  use core_kinds, SP_ => SP
  use core_parallel
  use core_hgroup
  use core_system

  use gyre_constants
  use gyre_discrim_func
  use gyre_wos_bvp
  use gyre_ext
  use gyre_num_par
  use gyre_root
  use gyre_status
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable    :: filename
  integer                      :: unit
  real(WP)                     :: c
  real(WP)                     :: x_i
  real(WP)                     :: x_o
  integer                      :: n_x
  real(WP)                     :: omega_min
  real(WP)                     :: omega_max
  integer                      :: n_omega
  character(256)               :: out_prefix
  type(num_par_t), allocatable :: np(:)

  integer                 :: i
  real(WP), allocatable   :: x(:)
  type(wos_bvp_t), target :: bp
  type(r_discrim_func_t)  :: df
  real(WP)                :: omega
  real(WP)                :: omega_prev
  type(r_ext_t)           :: discrim
  type(r_ext_t)           :: discrim_prev
  integer                 :: status
  integer                 :: n_root
  integer                 :: n_iter
  type(r_ext_t)           :: omega_root
  type(r_ext_t)           :: discrim_root
  real(WP), allocatable   :: y(:,:)
  real(WP)                :: y_ref(2)
  character(256)          :: out_file
  type(hgroup_t)          :: hg

  ! Initialize

  call init_parallel()

  ! Get command-line arguments

  $ASSERT(n_arg() == 1,Invalid number of arguments)

  call get_arg(1, filename)

  call set_log_level($str($LOG_LEVEL))

  if (check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_wos', '=')
100  format(A)

     write(OUTPUT_UNIT, 110) 'Compiler         :', COMPILER_VERSION()
     write(OUTPUT_UNIT, 110) 'Compiler options :', COMPILER_OPTIONS()
110  format(A,1X,A)

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   :', OMP_SIZE_MAX
120  format(A,1X,I0)

     write(OUTPUT_UNIT, 110) 'Input filename   :', filename

     write(OUTPUT_UNIT, 100) form_header('Initialization', '=')

  endif

  ! Process arguments

  open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

  call read_string_par(unit, c, x_i, x_o, n_x)
  call read_scan_par(unit, omega_min, omega_max, n_omega)
  call read_out_par(unit, out_prefix)

  call read_num_par(unit, np)

  close(unit)

  $ASSERT(SIZE(np) == 1,Must be exactly one &num_par namelist)

  ! Set up the grid

  allocate(x(n_x))

  grid_loop : do i = 1, n_x
     x(i) = (x_i*(n_x-i) + x_o*(i-1))/(n_x-1)
  end do grid_loop

  allocate(y(2,n_x))

  ! Set up the bvp

  bp = wos_bvp_t(x, c, np(1), omega_min, omega_max)

  ! Set up the discriminant function

  df = r_discrim_func_t(bp)

  ! Evaluate the discriminant between omega_min and omega_max, looking
  ! for roots

  omega = omega_min

  call df%eval(r_ext_t(omega), discrim, status)
  $ASSERT(status == STATUS_OK,Non-OK status returned from df%eval)

  n_root = 0

  omega_loop : do i = 2, n_omega

     ! Update omega

     omega_prev = omega

     omega = (omega_min*(n_omega-i) + omega_max*(i-1))/(n_omega-1)

     ! Update the discriminant

     discrim_prev = discrim

     call df%eval(r_ext_t(omega), discrim, status)
     $ASSERT(status == STATUS_OK,Non-OK status returned from df%eval)

     ! Look for a sign change

     if (discrim*discrim_prev <= 0._WP) then

        ! Solve for the discriminant root

        n_iter = 0

        call solve(df, np(1), r_ext_t(omega_prev), r_ext_t(omega), r_ext_t(0._WP), omega_root, &
                   status, n_iter=n_iter, n_iter_max=np(1)%n_iter_max, f_rx_a=discrim_prev, f_rx_b=discrim)
        $ASSERT(status == STATUS_OK,Non-PK status returned from solve)

       ! Reconstruct the eigenfunction

       call bp%recon(REAL(omega_root), x, x_o, y, y_ref, discrim_root)

       ! Print out the eigenfrequency

       print *,'Found eigenfrequency:', REAL(omega_root)

       ! Write out the eigenfunction

       n_root = n_root + 1

       write(out_file, 130) TRIM(out_prefix), n_root, '.h5'
130    format(A, I3.3, A)

       hg = hgroup_t(out_file, CREATE_FILE)

       call write_attr(hg, 'n_x', n_x)
       call write_attr(hg, 'omega', REAL(omega_root))

       call write_dset(hg, 'x', x)

       call write_dset(hg, 'y_1', y(1,:))
       call write_dset(hg, 'y_2', y(2,:))

       call hg%final()

    endif

  end do omega_loop

  ! Finish

  call final_parallel()

contains

  subroutine read_string_par (unit, c, x_i, x_o, n_x)

    integer, intent(in)   :: unit
    real(WP), intent(out) :: c
    real(WP), intent(out) :: x_i
    real(WP), intent(out) :: x_o
    integer, intent(out)  :: n_x

    namelist /string/ c, x_i, x_o, n_x

    ! Read string parameters

    rewind(unit)

    c = 1._WP

    x_i = 0._WP
    x_o = 1._WP

    n_x = 100

    read(unit, NML=string)

    ! Finish

    return

  end subroutine read_string_par

!****

  subroutine read_scan_par (unit, omega_min, omega_max, n_omega)

    integer, intent(in)   :: unit
    real(WP), intent(out) :: omega_min
    real(WP), intent(out) :: omega_max
    integer, intent(out)  :: n_omega

    namelist /scan/ omega_min, omega_max, n_omega

    ! Read eigenfrequency scan parameters

    rewind(unit)

    omega_min = 0._WP
    omega_max = 10._WP

    n_omega = 100

    read(unit, NML=scan)

    ! Finish

    return

  end subroutine read_scan_par

!****

  subroutine read_out_par (unit, out_prefix)

    integer, intent(in)       :: unit
    character(*), intent(out) :: out_prefix

    namelist /out/ out_prefix

    ! Read output parameters

    rewind(unit)

    out_prefix = 'out-'

    read(unit, NML=out)

    ! Finish

    return

  end subroutine read_out_par

end program gyre_wos
