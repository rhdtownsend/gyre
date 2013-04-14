! Program  : gyre_frontend
! Purpose  : frontend routines
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

module gyre_frontend

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel
  use core_order
  use core_memory

  use gyre_base_coeffs
  use gyre_therm_coeffs
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_eigfunc
  use gyre_output

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: open_input
  public :: write_header
  public :: init_coeffs
  public :: init_oscpar
  public :: init_numpar
  public :: init_scan
  public :: init_shoot_grid
  public :: init_recon_grid
  public :: write_data

contains

  subroutine open_input (unit)

    integer, intent(out) :: unit

    character(LEN=1024) :: line

    ! Make standard input available through a scratch file

    open(NEWUNIT=unit, STATUS='SCRATCH')

    do
       read(INPUT_UNIT, 100, END=200) line
100    format(A)
       write(unit, *) line
    end do

200 continue

    ! Finish

    return

  end subroutine open_input

!****

  subroutine write_header (header, underchar)

    character(LEN=*), intent(in)           :: header
    character(LEN=*), intent(in), optional :: underchar

    ! Write out the header

    if(MPI_RANK == 0) then

       write(OUTPUT_UNIT, '()')

       write(OUTPUT_UNIT, '(A)') header

       if(PRESENT(underchar)) then
          if(underchar == '') then
             write(OUTPUT_UNIT, '(A)') REPEAT(' ', LEN_TRIM(header))
          else
             write(OUTPUT_UNIT, '(A)') REPEAT(underchar, LEN_TRIM(header)/LEN_TRIM(underchar))
          endif
       endif

       write(OUTPUT_UNIT, '()')

    endif

    ! Finish

    return

  end subroutine write_header

!****

  subroutine init_coeffs (unit, x_bc, bc, tc)

    use gyre_mesa_file
    use gyre_b3_file
    use gyre_gsm_file
    use gyre_fgong_file
    use gyre_osc_file
    use gyre_poly_file
    use gyre_hom_base_coeffs

    integer, intent(in)                                         :: unit
    real(WP), allocatable, intent(out)                          :: x_bc(:)
    $if($GFORTRAN_PR56218)
    class(base_coeffs_t), allocatable, intent(inout)            :: bc
    class(therm_coeffs_t), allocatable, intent(inout), optional :: tc
    $else
    class(base_coeffs_t), allocatable, intent(out)              :: bc
    class(therm_coeffs_t), allocatable, intent(out), optional   :: tc
    $endif

    character(LEN=256)          :: coeffs_type
    character(LEN=256)          :: deriv_type
    character(LEN=FILENAME_LEN) :: file
    real(WP)                    :: G
    real(WP)                    :: Gamma_1

    namelist /coeffs/ coeffs_type, deriv_type, file, G, Gamma_1

    ! Read structure coefficients parameters

    coeffs_type = ''
    deriv_type = 'MONO'

    file = ''

    G = G_GRAVITY
    Gamma_1 = 5._WP/3._WP

    rewind(unit)
    read(unit, NML=coeffs, END=900)

    ! Read/initialize the base_coeffs

    select case(coeffs_type)
    case('MESA')
       call read_mesa_file(file, G, deriv_type, bc, tc, x=x_bc)
    case('B3')
       call read_b3_file(file, G, deriv_type, bc, tc, x=x_bc)
    case('GSM')
       call read_gsm_file(file, G, deriv_type, bc, tc, x=x_bc)
    case('OSC')
       call read_osc_file(file, G, deriv_type, bc, tc, x=x_bc)
    case('FGONG')
       call read_fgong_file(file, G, deriv_type, bc, x=x_bc) 
    case('POLY')
       call read_poly_file(file, deriv_type, bc, x=x_bc)
    case('HOM')
       allocate(hom_base_coeffs_t::bc)
       select type (bc)
       type is (hom_base_coeffs_t)
          call bc%init(Gamma_1)
       end select
    end select

    ! Finish

    return

    ! Jump-in point for end-of-file

900 continue

    $ABORT(No &coeffs namelist in input file)

  end subroutine init_coeffs

!****

  subroutine init_oscpar (unit, op)

    integer, intent(in)         :: unit
    type(oscpar_t), intent(out) :: op

    integer           :: l
    character(LEN=64) :: outer_bound_type

    namelist /osc/ l, outer_bound_type

    ! Read oscillation parameters

    l = 0
    outer_bound_type = 'ZERO'

    rewind(unit)
    read(unit, NML=osc, END=900)

    ! Initialize the oscpar

    op = oscpar_t(l=l, outer_bound_type=outer_bound_type)

    ! Finish

    return

    ! Jump-in point for end-of-file

900 continue

    $ABORT(No &osc namelist in input file)

  end subroutine init_oscpar

!****

  subroutine init_numpar (unit, np)

    integer, intent(in)         :: unit
    type(numpar_t), intent(out) :: np

    integer           :: n_iter_max
    real(WP)          :: theta_ad
    character(LEN=64) :: ivp_solver_type

    namelist /num/ n_iter_max, theta_ad, ivp_solver_type

    ! Read numerical parameters

    n_iter_max = 50
    theta_ad = 0._WP

    ivp_solver_type = 'MAGNUS_GL2'

    rewind(unit)
    read(unit, NML=num, END=900)

    ! Initialize the numpar

    np = numpar_t(n_iter_max=n_iter_max, theta_ad=theta_ad, ivp_solver_type=ivp_solver_type)

    ! Finish

    return

    ! Jump-in point for end-of-file

900 continue

    $ABORT(No &num namelist in input file)

  end subroutine init_numpar

!****

  subroutine init_scan (unit, bc, omega)

    integer, intent(in)                :: unit
    class(base_coeffs_t), intent(in)   :: bc
    real(WP), allocatable, intent(out) :: omega(:)

    character(LEN=256) :: grid_type
    real(WP)           :: freq_min
    real(WP)           :: freq_max
    integer            :: n_freq
    character(LEN=256) :: freq_units
    real(WP)           :: omega_min
    real(WP)           :: omega_max
    integer            :: i

    namelist /scan/ grid_type, freq_min, freq_max, n_freq, freq_units

    ! Read scan parameters

    rewind(unit)

    allocate(omega(0))

    read_loop : do 

       grid_type = 'LINEAR'

       freq_min = 1._WP
       freq_max = 10._WP
       n_freq = 10
          
       freq_units = 'NONE'

       read(unit, NML=scan, END=100)
          
       ! Set up the frequency grid

       omega_min = REAL(bc%conv_freq(CMPLX(freq_min, KIND=WP), freq_units, 'NONE'))
       omega_max = REAL(bc%conv_freq(CMPLX(freq_max, KIND=WP), freq_units, 'NONE'))
       
       select case(grid_type)
       case('LINEAR')
          omega = [omega,(((n_freq-i)*omega_min + (i-1)*omega_max)/(n_freq-1), i=1,n_freq)]
       case('INVERSE')
          omega = [omega,((n_freq-1)/((n_freq-i)/omega_min + (i-1)/omega_max), i=1,n_freq)]
       case default
          $ABORT(Invalid freq_grid)
       end select

    end do read_loop

100 continue

    ! Sort the frequencies

    omega = omega(sort_indices(omega))

    ! Finish

    return

  end subroutine init_scan

!****

  $define $INIT_GRID $sub

  $local $NAME $1

  subroutine init_${NAME}_grid (unit, omega_a, omega_b, gp)

    integer, intent(in)                       :: unit
    real(WP), intent(in)                      :: omega_a
    real(WP), intent(in)                      :: omega_b
    type(gridpar_t), allocatable, intent(out) :: gp(:)

    character(LEN=256) :: op_type
    real(WP)           :: alpha_osc
    real(WP)           :: alpha_exp
    real(WP)           :: s
    integer            :: n
    integer            :: n_grid
    integer            :: i

    namelist /${NAME}_grid/ op_type, alpha_osc, alpha_exp, s, n

    ! Count the number of grid namelists

    rewind(unit)

    n_grid = 0

    count_loop : do
       read(unit, NML=${NAME}_grid, END=100)
       n_grid = n_grid + 1
    end do count_loop

100 continue

    $ASSERT(n_grid >= 1,At least one ${NAME}_grid namelist is required)

    ! Read grid parameters

    rewind(unit)

    allocate(gp(n_grid))

    read_loop : do i = 1,n_grid

       op_type = 'CREATE_CLONE'

       read(unit, NML=${NAME}_grid)

       gp(i) = gridpar_t(op_type=op_type, &
                         alpha_osc=alpha_osc, alpha_exp=alpha_exp, &
                         omega_a=omega_a, omega_b=omega_b, &
                         s=s, n=n)

    end do read_loop

    ! Finish

    return

  end subroutine init_${NAME}_grid

  $endsub

  $INIT_GRID(shoot)
  $INIT_GRID(recon)

!****

  subroutine write_data (unit, ef, bc)

    integer, intent(in)              :: unit
    type(eigfunc_t), intent(in)      :: ef(:)
    class(base_coeffs_t), intent(in) :: bc

    character(LEN=256)          :: freq_units
    character(LEN=FILENAME_LEN) :: summary_file
    character(LEN=2048)         :: summary_item_list
    character(LEN=FILENAME_LEN) :: mode_prefix
    character(LEN=2048)         :: mode_item_list
    character(LEN=FILENAME_LEN) :: mode_file
    integer                     :: j

    namelist /output/ freq_units, summary_file, summary_item_list, mode_prefix, mode_item_list

    ! Read output parameters

    freq_units = 'NONE'

    summary_file = ''
    summary_item_list = 'l,n_p,n_g,omega,freq'

    mode_prefix = ''
    mode_item_list = TRIM(summary_item_list)//',x,xi_r,xi_h'

    rewind(unit)
    read(unit, NML=output, END=900)

    ! Write output files

    if(summary_file /= '') call write_summary(summary_file, ef, bc, split_item_list(summary_item_list), freq_units)

    if(mode_prefix /= '') then

       mode_loop : do j = 1,SIZE(ef)

          write(mode_file, 100) TRIM(mode_prefix), j, '.h5'
100       format(A,I4.4,A)

          call write_mode(mode_file, ef(j), bc, split_item_list(mode_item_list), freq_units, j)

       end do mode_loop
       
    end if

    ! Finish

    return

    ! Jump-in point for end-of-file

900 continue

    $ABORT(No &output namelist in input file)

  end subroutine write_data

!****

  function split_item_list (item_list) result (items)

    character(LEN=*), intent(in)               :: item_list
    character(LEN=LEN(item_list)), allocatable :: items(:)

    character(LEN=LEN(item_list)) :: item_list_
    integer                       :: d
    integer                       :: n
    integer                       :: j
    
    ! Split the comma-sepated list of items into an array of
    ! individual items
 
    d = 16

    allocate(items(d))

    n = 0

    ! Repeatedly split on commas

    item_list_ = item_list

    split_loop : do

       if(item_list_ == ' ') exit split_loop

       j = INDEX(item_list_, ',')

       if(j <= 0) then
          n = n + 1
          items(n) = item_list_
          exit split_loop
       endif

       n = n + 1

       ! Chop out the item name

       items(n) = item_list_(:j-1)
       item_list_ = item_list_(j+1:)

       ! If necessary, expand the list
          
       if(n >= d) then
          d = 2*d
          call reallocate(items, [d])
       end if

    end do split_loop

    ! Reallocate item_list to the correct length

    call reallocate(items, [n])

    ! Finish

    return

  end function split_item_list

end module gyre_frontend
