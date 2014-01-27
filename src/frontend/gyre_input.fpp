! Program  : gyre_input
! Purpose  : input routines
!
! Copyright 2013 Rich Townsend
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

module gyre_input

  ! Uses

  use core_kinds
  use core_constants
  use core_order
  use core_parallel

  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_scanpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: parse_args
  public :: read_coeffs
  public :: read_oscpar
  public :: read_numpar
  public :: read_scanpar
  public :: read_shoot_gridpar
  public :: read_recon_gridpar

contains

  subroutine parse_args (filename)

    character(LEN=:), allocatable, intent(out) :: filename

    integer :: n
    integer :: length

    ! Parse the command-line arguments

    n = COMMAND_ARGUMENT_COUNT()

    $ASSERT(n == 1,Invalid number of arguments)

    call GET_COMMAND_ARGUMENT(1, LENGTH=length)
    allocate(character(LEN=length) :: filename)

    call GET_COMMAND_ARGUMENT(1, VALUE=filename)

    ! Finish

    return

  end subroutine parse_args

!****

  subroutine read_coeffs (unit, x_bc, cf)

    use gyre_coeffs
    use gyre_coeffs_evol
    use gyre_coeffs_poly
    use gyre_coeffs_hom
    use gyre_mesa_file
    use gyre_osc_file
    use gyre_fgong_file
    use gyre_famdl_file
    $if ($HDF5)
    use gyre_b3_file
    use gyre_gsm_file
    use gyre_poly_file
    $endif

    integer, intent(in)                   :: unit
    real(WP), allocatable, intent(out)    :: x_bc(:)
    class(coeffs_t), pointer, intent(out) :: cf

    character(LEN=256)          :: coeffs_type
    character(LEN=256)          :: file_format
    character(LEN=256)          :: data_format
    character(LEN=256)          :: deriv_type
    character(LEN=FILENAME_LEN) :: file
    real(WP)                    :: G
    real(WP)                    :: Gamma_1
    type(coeffs_evol_t)         :: ec
    type(coeffs_poly_t)         :: pc
    type(coeffs_hom_t)          :: hc

    namelist /coeffs/ coeffs_type, file_format, data_format, deriv_type, file, G, Gamma_1

    ! Read structure coefficients parameters

    coeffs_type = ''
    file_format = ''
    data_format = ''
    deriv_type = 'MONO'

    file = ''

    G = G_GRAVITY
    Gamma_1 = 5._WP/3._WP

    rewind(unit)
    read(unit, NML=coeffs, END=900)

    ! Read/initialize the base_coeffs

    select case (coeffs_type)
    case ('EVOL')

       select case (file_format)
       case ('MESA')
          call read_mesa_file(file, G, deriv_type, ec, x=x_bc)
       case('B3')
          $if($HDF5)
          call read_b3_file(file, G, deriv_type, ec, x=x_bc)
          $else
          $ABORT(No HDF5 support, therefore cannot read B3-format files)
          $endif
       case ('GSM')
          $if($HDF5)
          call read_gsm_file(file, G, deriv_type, ec, x=x_bc)
          $else
          $ABORT(No HDF5 support, therefore cannot read GSM-format files)
          $endif
       case ('OSC')
          call read_osc_file(file, G, deriv_type, data_format, ec, x=x_bc)
       case ('FGONG')
          call read_fgong_file(file, G, deriv_type, data_format, ec, x=x_bc) 
       case ('FAMDL')
          call read_famdl_file(file, G, deriv_type, data_format, ec, x=x_bc) 
       case default
          $ABORT(Invalid file_format)
       end select

       allocate(cf, SOURCE=ec)
       
    case ('POLY')

       $if($HDF5)
       call read_poly_file(file, deriv_type, pc, x=x_bc)
       $else
       $ABORT(No HDF5 support, therefore cannot read POLY files)
       $endif

       allocate(cf, SOURCE=pc)

    case ('HOM')

       hc = coeffs_hom_t(Gamma_1)

       allocate(cf, SOURCE=hc)

    case default

       $ABORT(Invalid coeffs_type)

    end select

    ! Finish

    return

    ! Jump-in point for end-of-file

900 continue

    $ABORT(No &coeffs namelist in input file)

  end subroutine read_coeffs

!****

  subroutine read_oscpar (unit, op)

    integer, intent(in)                      :: unit
    type(oscpar_t), allocatable, intent(out) :: op(:)

    integer           :: n_op
    integer           :: i
    integer           :: l
    integer           :: m
    integer           :: X_n_pg_min
    integer           :: X_n_pg_max
    character(LEN=64) :: variables_type
    character(LEN=64) :: outer_bound_type
    character(LEN=64) :: inertia_norm_type
    character(LEN=64) :: tag
    real(WP)          :: x_ref

    namelist /osc/ x_ref, l, m, X_n_pg_min, X_n_pg_max, &
         outer_bound_type, variables_type, &
         inertia_norm_type, tag, x_ref

    ! Count the number of grid namelists

    rewind(unit)

    n_op = 0

    count_loop : do
       read(unit, NML=osc, END=100)
       n_op = n_op + 1
    end do count_loop

100 continue

    ! Read oscillation parameters

    rewind(unit)

    allocate(op(n_op))

    read_loop : do i = 1,n_op

       l = 0
       m = 0

       X_n_pg_min = -HUGE(0)
       X_n_pg_max = HUGE(0)

       variables_type = 'DZIEM'
       outer_bound_type = 'ZERO'
       inertia_norm_type = 'BOTH'
       tag = ''

       x_ref = HUGE(0._WP)

       read(unit, NML=osc)

       ! Initialize the oscpar

       op(i) = oscpar_t(l=l, m=m, X_n_pg_min=X_n_pg_min, X_n_pg_max=X_n_pg_max, &
                        variables_type=variables_type, outer_bound_type=outer_bound_type, &
                        inertia_norm_type=inertia_norm_type, tag=tag, x_ref=x_ref)

    end do read_loop

    ! Finish

    return

  end subroutine read_oscpar

!****

  subroutine read_numpar (unit, np)

    integer, intent(in)                      :: unit
    type(numpar_t), allocatable, intent(out) :: np(:)

    integer             :: n_np
    integer             :: n_iter_max
    real(WP)            :: theta_ad
    logical             :: reduce_order
    logical             :: use_banded
    logical             :: use_trad_approx
    character(LEN=64)   :: ivp_solver_type
    character(LEN=2048) :: tag_list
    integer             :: i

    namelist /num/ n_iter_max, theta_ad, &
         reduce_order, use_banded, use_trad_approx, ivp_solver_type, tag_list

    ! Count the number of num namelists

    rewind(unit)

    n_np = 0

    count_loop : do
       read(unit, NML=num, END=100)
       n_np = n_np + 1
    end do count_loop

100 continue

    ! Read numerical parameters

    rewind(unit)

    allocate(np(n_np))

    read_loop : do i = 1,n_np

       n_iter_max = 50
       theta_ad = 0._WP

       reduce_order = .TRUE.
       use_banded = .FALSE.
       use_trad_approx = .FALSE.

       ivp_solver_type = 'MAGNUS_GL2'
       tag_list = ''

       read(unit, NML=num)

       ! Initialize the numpar

       np(i) = numpar_t(n_iter_max=n_iter_max, theta_ad=theta_ad, &
                        reduce_order=reduce_order, use_banded=use_banded, use_trad_approx=use_trad_approx, &
                       ivp_solver_type=ivp_solver_type, tag_list=tag_list)

    end do read_loop

    ! Finish

    return

  end subroutine read_numpar

!****

  $define $READ_GRIDPAR $sub

  $local $NAME $1

  subroutine read_${NAME}_gridpar (unit, gp)

    integer, intent(in)                       :: unit
    type(gridpar_t), allocatable, intent(out) :: gp(:)

    integer                     :: n_gp
    real(WP)                    :: alpha_osc
    real(WP)                    :: alpha_exp
    real(WP)                    :: alpha_thm
    real(WP)                    :: s
    integer                     :: n
    character(LEN=FILENAME_LEN) :: file
    character(LEN=64)           :: op_type
    character(LEN=2048)         :: tag_list
    integer                     :: i

    namelist /${NAME}_grid/ alpha_osc, alpha_exp, alpha_thm, s, n, file, op_type, tag_list

    ! Count the number of grid namelists

    rewind(unit)

    n_gp = 0

    count_loop : do
       read(unit, NML=${NAME}_grid, END=100)
       n_gp = n_gp + 1
    end do count_loop

100 continue

    ! Read grid parameters

    rewind(unit)

    allocate(gp(n_gp))

    read_loop : do i = 1,n_gp

       alpha_osc = 0._WP
       alpha_exp = 0._WP
       alpha_thm = 0._WP

       s = 0._WP

       n = 0

       file = ''

       op_type = 'CREATE_CLONE'
       tag_list = ''

       read(unit, NML=${NAME}_grid)

       ! Initialize the gridpar

       gp(i) = gridpar_t(alpha_osc=alpha_osc, alpha_exp=alpha_exp, alpha_thm=alpha_thm, &
                         omega_a=0._WP, omega_b=0._WP, &
                         s=s, n=n, file=file, op_type=op_type, tag_list=tag_list)

    end do read_loop

    ! Finish

    return

  end subroutine read_${NAME}_gridpar

  $endsub

  $READ_GRIDPAR(shoot)
  $READ_GRIDPAR(recon)

!****

  subroutine read_scanpar (unit, sp)

    integer, intent(in)                       :: unit
    type(scanpar_t), allocatable, intent(out) :: sp(:)

    integer             :: n_sp
    integer             :: i
    real(WP)            :: freq_min
    real(WP)            :: freq_max
    integer             :: n_freq
    character(LEN=64)   :: freq_units
    character(LEN=64)   :: grid_type
    character(LEN=2048) :: tag_list

    namelist /scan/ freq_min, freq_max, n_freq, freq_units, grid_type, tag_list

    ! Count the number of scan namelists

    rewind(unit)

    n_sp = 0

    count_loop : do
       read(unit, NML=scan, END=100)
       n_sp = n_sp + 1
    end do count_loop

100 continue

    ! Read scan parameters

    rewind(unit)

    allocate(sp(n_sp))

    read_loop : do i = 1,n_sp

       freq_min = 1._WP
       freq_max = 10._WP
       n_freq = 10
          
       freq_units = 'NONE'
       grid_type = 'LINEAR'
       tag_list = ''

       read(unit, NML=scan)

       ! Initialize the scanpar

       sp(i) = scanpar_t(freq_min=freq_min, freq_max=freq_max, n_freq=n_freq, &
                         freq_units=freq_units, grid_type=grid_type, tag_list=tag_list)

    end do read_loop

    ! Finish

    return

  end subroutine read_scanpar

end module gyre_input
