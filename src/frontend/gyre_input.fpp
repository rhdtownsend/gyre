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
  use core_order
  use core_parallel

  use gyre_constants
  use gyre_modepar
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_scanpar
  use gyre_outpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: parse_args
  public :: read_constants
  public :: read_model
  public :: read_modepar
  public :: read_oscpar
  public :: read_numpar
  public :: read_scanpar
  public :: read_shoot_gridpar
  public :: read_recon_gridpar
  public :: read_outpar

contains

  subroutine parse_args (filename)

    character(:), allocatable, intent(out) :: filename

    integer :: n
    integer :: length

    ! Parse the command-line arguments

    n = COMMAND_ARGUMENT_COUNT()

    $ASSERT(n == 1,Invalid number of arguments)

    call GET_COMMAND_ARGUMENT(1, LENGTH=length)
    allocate(character(length) :: filename)

    call GET_COMMAND_ARGUMENT(1, VALUE=filename)

    ! Finish

    return

  end subroutine parse_args

!****

  subroutine read_model (unit, x_bc, ml)

    use gyre_model
    use gyre_evol_model
    use gyre_scons_model
    use gyre_poly_model
    use gyre_hom_model
    use gyre_mesa_file
    use gyre_osc_file
    use gyre_losc_file
    use gyre_fgong_file
    use gyre_famdl_file
    use gyre_amdl_file
    $if ($HDF5)
    use gyre_b3_file
    use gyre_gsm_file
    use gyre_poly_file
    $endif

    integer, intent(in)                  :: unit
    real(WP), allocatable, intent(out)   :: x_bc(:)
    class(model_t), pointer, intent(out) :: ml

    integer                 :: n_ml
    character(256)          :: model_type
    character(256)          :: file_format
    character(256)          :: data_format
    character(256)          :: deriv_type
    character(FILENAME_LEN) :: file
    real(WP)                :: Gamma_1
    real(WP)                :: Omega_rot
    logical                 :: regularize
    type(evol_model_t)      :: ec
    type(scons_model_t)     :: sc
    type(poly_model_t)      :: pc
    type(hom_model_t)       :: hc

    namelist /model/ model_type, file_format, data_format, deriv_type, file, Gamma_1, Omega_rot, regularize

    ! Count the number of model namelists

    rewind(unit)

    n_ml = 0

    count_loop : do
       read(unit, NML=model, END=100)
       n_ml = n_ml + 1
    end do count_loop

100 continue

    $ASSERT(n_ml == 1,Input file should contain exactly one &model namelist)

    ! Read model parameters

    model_type = ''
    file_format = ''
    data_format = ''
    deriv_type = 'MONO'
    regularize = .FALSE.

    file = ''

    Gamma_1 = 5._WP/3._WP
    Omega_rot = 0._WP

    rewind(unit)
    read(unit, NML=model)

    ! Read/initialize the model

    select case (model_type)
    case ('EVOL')

       select case (file_format)
       case ('MESA')
          call read_mesa_model(file, deriv_type, regularize, ec, x=x_bc)
       case('B3')
          $if($HDF5)
          call read_b3_model(file, deriv_type, regularize, ec, x=x_bc)
          $else
          $ABORT(No HDF5 support, therefore cannot read B3-format files)
          $endif
       case ('GSM')
          $if($HDF5)
          call read_gsm_model(file, deriv_type, regularize, ec, x=x_bc)
          $else
          $ABORT(No HDF5 support, therefore cannot read GSM-format files)
          $endif
       case ('OSC')
          call read_osc_model(file, deriv_type, data_format, regularize, ec, x=x_bc)
       case ('LOSC')
          call read_losc_model(file, deriv_type, regularize, ec, x=x_bc)
       case ('FGONG')
          call read_fgong_model(file, deriv_type, data_format, regularize, ec, x=x_bc) 
       case ('FAMDL')
          call read_famdl_model(file, deriv_type, data_format, regularize, ec, x=x_bc)
       case ('AMDL')
          call read_amdl_model(file, deriv_type, regularize, ec, x=x_bc)
       case default
          $ABORT(Invalid file_format)
       end select

       allocate(ml, SOURCE=ec)

    case ('SCONS')

       select case (file_format)
       case ('MESA')
          call read_mesa_model(file, sc, x=x_bc)
       case ('FGONG')
          call read_fgong_model(file, data_format, sc, x=x_bc)
       case default
          $ABORT(Invalid file_format)
       end select

       allocate(ml, SOURCE=sc)

    case ('POLY')

       $if($HDF5)
       call read_poly_model(file, deriv_type, pc, x=x_bc)
       $else
       $ABORT(No HDF5 support, therefore cannot read POLY files)
       $endif

       allocate(ml, SOURCE=pc)

    case ('HOM')

       hc = hom_model_t(Gamma_1, Omega_rot)

       allocate(ml, SOURCE=hc)

    case default

       $ABORT(Invalid model_type)

    end select

    ! Finish

    return

  end subroutine read_model

!****

  subroutine read_constants (unit)

    integer, intent(in) :: unit

    integer :: n_cn

    namelist /constants/ G_GRAVITY, C_LIGHT, A_RADIATION, &
                         M_SUN, R_SUN, L_SUN

    ! Count the number of constants namelists

    rewind(unit)

    n_cn = 0

    count_loop : do
       read(unit, NML=constants, END=100)
       n_cn = n_cn + 1
    end do count_loop

100 continue

    $ASSERT(n_cn == 1,Input file should contain exactly one &constants namelist)

    ! Read constants

    rewind(unit)
    read(unit, NML=constants)

    ! Finish

    return

  end subroutine read_constants

!****

  subroutine read_modepar (unit, mp)

    integer, intent(in)                       :: unit
    type(modepar_t), allocatable, intent(out) :: mp(:)

    integer       :: n_mp
    integer       :: i
    integer       :: l
    integer       :: m
    integer       :: n_pg_min
    integer       :: n_pg_max
    character(64) :: tag

    namelist /mode/ l, m, n_pg_min, n_pg_max, tag

    ! Count the number of mode namelists

    rewind(unit)

    n_mp = 0

    count_loop : do
       read(unit, NML=mode, END=100)
       n_mp = n_mp + 1
    end do count_loop

100 continue

    ! Read mode parameters

    rewind(unit)

    allocate(mp(n_mp))

    read_loop : do i = 1,n_mp

       l = 0
       m = 0

       n_pg_min = -HUGE(0)
       n_pg_max = HUGE(0)

       tag = ''

       read(unit, NML=mode)

       ! Initialize the modepar

       mp(i) = modepar_t(l=l, m=m, n_pg_min=n_pg_min, n_pg_max=n_pg_max, tag=tag)

    end do read_loop

    ! Finish

    return

  end subroutine read_modepar

!****

  subroutine read_oscpar (unit, op)

    integer, intent(in)                      :: unit
    type(oscpar_t), allocatable, intent(out) :: op(:)

    integer         :: n_op
    integer         :: i
    logical         :: reduce_order
    character(64)   :: rot_method
    character(64)   :: variables_type
    character(64)   :: inner_bound_type
    character(64)   :: outer_bound_type
    character(64)   :: inertia_norm_type
    character(2048) :: tag_list
    real(WP)        :: x_ref

    namelist /osc/ x_ref, rot_method, inner_bound_type, outer_bound_type, variables_type, &
         inertia_norm_type, tag_list, reduce_order

    ! Count the number of osc namelists

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

       x_ref = HUGE(0._WP)

       rot_method = 'DOPPLER'
       variables_type = 'DZIEM'
       inner_bound_type = 'REGULAR'
       outer_bound_type = 'ZERO'
       inertia_norm_type = 'BOTH'
       tag_list = ''

       reduce_order = .TRUE.

       read(unit, NML=osc)

       ! Initialize the oscpar

       op(i) = oscpar_t(x_ref=x_ref, rot_method=rot_method, variables_type=variables_type, &
                        inner_bound_type=inner_bound_type, outer_bound_type=outer_bound_type, &
                        inertia_norm_type=inertia_norm_type, tag_list=tag_list, &
                        reduce_order=reduce_order)

    end do read_loop

    ! Finish

    return

  end subroutine read_oscpar

!****

  subroutine read_numpar (unit, np)

    integer, intent(in)                      :: unit
    type(numpar_t), allocatable, intent(out) :: np(:)

    integer         :: n_np
    integer         :: i
    integer         :: n_iter_max
    logical         :: use_banded
    logical         :: use_trad_approx
    logical         :: deflate_roots
    character(64)   :: ivp_solver_type
    character(2048) :: tag_list

    namelist /num/ n_iter_max, &
         use_banded, use_trad_approx, deflate_roots, &
         ivp_solver_type, tag_list

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

       use_banded = .FALSE.
       use_trad_approx = .FALSE.
       deflate_roots = .TRUE.

       ivp_solver_type = 'MAGNUS_GL2'
       tag_list = ''

       read(unit, NML=num)

       ! Initialize the numpar

       np(i) = numpar_t(n_iter_max=n_iter_max, &
                        use_banded=use_banded, &
                        use_trad_approx=use_trad_approx, deflate_roots=deflate_roots, &
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

    integer                 :: n_gp
    integer                 :: i
    real(WP)                :: alpha_osc
    real(WP)                :: alpha_exp
    real(WP)                :: alpha_thm
    real(WP)                :: alpha_str
    real(WP)                :: s
    integer                 :: n
    character(FILENAME_LEN) :: file
    character(64)           :: op_type
    character(2048)         :: tag_list

    namelist /${NAME}_grid/ alpha_osc, alpha_exp, alpha_thm, alpha_str, s, n, file, op_type, tag_list

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

    read_loop : do i = 1, n_gp

       alpha_osc = 0._WP
       alpha_exp = 0._WP
       alpha_thm = 0._WP
       alpha_str = 0._WP

       s = 0._WP

       n = 0

       file = ''

       op_type = 'CREATE_CLONE'
       tag_list = ''

       read(unit, NML=${NAME}_grid)

       ! Initialize the gridpar

       gp(i) = gridpar_t(alpha_osc=alpha_osc, alpha_exp=alpha_exp, alpha_thm=alpha_thm, alpha_str=alpha_str, &
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

    integer         :: n_sp
    integer         :: i
    real(WP)        :: freq_min
    real(WP)        :: freq_max
    integer         :: n_freq
    character(64)   :: freq_units
    character(64)   :: freq_frame
    character(64)   :: grid_type
    character(2048) :: tag_list

    namelist /scan/ freq_min, freq_max, n_freq, freq_units, freq_frame, grid_type, tag_list

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

    read_loop : do i = 1, n_sp

       freq_min = 1._WP
       freq_max = 10._WP
       n_freq = 10
          
       freq_units = 'NONE'
       freq_frame = 'INERTIAL'
       grid_type = 'LINEAR'
       tag_list = ''

       read(unit, NML=scan)

       ! Initialize the scanpar

       sp(i) = scanpar_t(freq_min=freq_min, freq_max=freq_max, n_freq=n_freq, &
                         freq_units=freq_units, freq_frame=freq_frame, &
                         grid_type=grid_type, tag_list=tag_list)

    end do read_loop

    ! Finish

    return

  end subroutine read_scanpar

!****

  subroutine read_outpar (unit, up)

    integer, intent(in)         :: unit
    type(outpar_t), intent(out) :: up

    integer                 :: n_up
    character(256)          :: freq_units
    character(FILENAME_LEN) :: summary_file
    character(256)          :: summary_file_format
    character(2048)         :: summary_item_list
    character(FILENAME_LEN) :: mode_prefix
    character(FILENAME_LEN) :: mode_template
    character(256)          :: mode_file_format
    character(2048)         :: mode_item_list
    logical                 :: prune_modes

    namelist /output/ freq_units, summary_file, summary_file_format, summary_item_list, &
                      mode_prefix, mode_template, mode_file_format, mode_item_list, prune_modes

    ! Count the number of output namelists

    rewind(unit)

    n_up = 0

    count_loop : do
       read(unit, NML=output, END=100)
       n_up = n_up + 1
    end do count_loop

100 continue

    $ASSERT(n_up == 1,Input file should contain exactly one &output namelist)

    ! Read output parameters

    freq_units = 'NONE'

    summary_file = ''
    summary_file_format = 'HDF'
    summary_item_list = 'l,n_pg,omega,freq'
    
    mode_prefix = ''
    mode_template = ''
    mode_file_format = 'HDF'
    mode_item_list = TRIM(summary_item_list)//',x,xi_r,xi_h'

    prune_modes = .FALSE.

    rewind(unit)
    read(unit, NML=output)

    ! Initialize the outpar

    up = outpar_t(freq_units=freq_units, &
                  summary_file=summary_file, summary_file_format=summary_file_format, summary_item_list=summary_item_list, &
                  mode_prefix=mode_prefix, mode_template=mode_template, mode_file_format=mode_file_format, mode_item_list=mode_item_list, &
                  prune_modes=prune_modes)

    ! Finish

    return

  end subroutine read_outpar

end module gyre_input
