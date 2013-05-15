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
  use gyre_mode
  use gyre_output

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: parse_args
  public :: write_header
  public :: init_coeffs
  public :: init_oscpar
  public :: init_numpar
  public :: init_scan
  public :: init_shoot_grid
  public :: init_recon_grid
  public :: write_data
  public :: freq_scale

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
    character(LEN=256)          :: file_format
    character(LEN=256)          :: deriv_type
    character(LEN=FILENAME_LEN) :: file
    real(WP)                    :: G
    real(WP)                    :: Gamma_1

    namelist /coeffs/ coeffs_type, file_format, deriv_type, file, G, Gamma_1

    ! Read structure coefficients parameters

    coeffs_type = ''
    file_format = ''
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
          call read_mesa_file(file, G, deriv_type, bc, tc, x=x_bc)
       case('B3')
          call read_b3_file(file, G, deriv_type, bc, tc, x=x_bc)
       case ('GSM')
          call read_gsm_file(file, G, deriv_type, bc, tc, x=x_bc)
       case ('OSC')
          call read_osc_file(file, G, deriv_type, bc, tc, x=x_bc)
       case ('FGONG')
          call read_fgong_file(file, G, deriv_type, bc, x=x_bc) 
       case default
          $ABORT(Invalid file_format)
       end select
    case ('POLY')
       call read_poly_file(file, deriv_type, bc, x=x_bc)
    case ('HOM')
       allocate(hom_base_coeffs_t::bc)
       select type (bc)
       type is (hom_base_coeffs_t)
          call bc%init(Gamma_1)
       end select
    case default
       $ABORT(Invalid coeffs_type)
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

  subroutine init_scan (unit, bc, op, omega)

    integer, intent(in)                :: unit
    class(base_coeffs_t), intent(in)   :: bc
    type(oscpar_t), intent(in)         :: op
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

       omega_min = freq_min/freq_scale(bc, op, freq_units)
       omega_max = freq_max/freq_scale(bc, op, freq_units)
       
       select case(grid_type)
       case('LINEAR')
          omega = [omega,(((n_freq-i)*omega_min + (i-1)*omega_max)/(n_freq-1), i=1,n_freq)]
       case('INVERSE')
          omega = [omega,((n_freq-1)/((n_freq-i)/omega_min + (i-1)/omega_max), i=1,n_freq)]
       case default
          $ABORT(Invalid grid_type)
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

    integer            :: n_gp
    character(LEN=256) :: op_type
    real(WP)           :: alpha_osc
    real(WP)           :: alpha_exp
    real(WP)           :: s
    integer            :: n
    integer            :: i

    namelist /${NAME}_grid/ op_type, alpha_osc, alpha_exp, s, n

    ! Count the number of grid namelists

    rewind(unit)

    n_gp = 0

    count_loop : do
       read(unit, NML=${NAME}_grid, END=100)
       n_gp = n_gp + 1
    end do count_loop

100 continue

    $ASSERT(n_gp >= 1,At least one ${NAME}_grid namelist is required)

    ! Read grid parameters

    rewind(unit)

    allocate(gp(n_gp))

    read_loop : do i = 1,n_gp

       op_type = 'CREATE_CLONE'

       alpha_osc = 0._WP
       alpha_exp = 0._WP

       s = 0._WP

       n = 0

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

  subroutine write_data (unit, md)

    integer, intent(in)         :: unit
    type(mode_t), intent(in)    :: md(:)

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

    if(summary_file /= '') call write_summary(summary_file, md, split_item_list(summary_item_list), &
                                              freq_scale(md(1)%bc, md(1)%op, freq_units))

    if(mode_prefix /= '') then

       mode_loop : do j = 1,SIZE(md)

          write(mode_file, 100) TRIM(mode_prefix), j, '.h5'
100       format(A,I4.4,A)

          call write_mode(mode_file, md(j), split_item_list(mode_item_list), &
                          freq_scale(md(j)%bc, md(j)%op, freq_units), j)

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

!****

  function freq_scale (bc, op, freq_units)

    use gyre_evol_base_coeffs
    use gyre_poly_base_coeffs
    use gyre_hom_base_coeffs

    class(base_coeffs_t), intent(in) :: bc
    type(oscpar_t), intent(in)       :: op
    character(LEN=*), intent(in)     :: freq_units
    real(WP)                         :: freq_scale

    ! Calculate the scale factor to convert a dimensionless angular
    ! frequency to a dimensioned frequency

    select type (bc)
    class is (evol_base_coeffs_t)
       freq_scale = evol_freq_scale(bc, op, freq_units)
    class is (poly_base_coeffs_t)
       freq_scale = poly_freq_scale(freq_units)
    class is (hom_base_coeffs_t)
       freq_scale = hom_freq_scale(freq_units)
    class default
       $ABORT(Invalid bc type)
    end select

    ! Finish

    return

  contains

    function evol_freq_scale (bc, op, freq_units) result (freq_scale)

      use gyre_ad_bound

      class(evol_base_coeffs_t), intent(in) :: bc
      type(oscpar_t), intent(in)            :: op
      character(LEN=*), intent(in)          :: freq_units
      real(WP)                              :: freq_scale

      real(WP) :: omega_cutoff_lo
      real(WP) :: omega_cutoff_hi

      ! Calculate the scale factor to convert a dimensionless angular
      ! frequency to a dimensioned frequency

      select case(freq_units)
      case('NONE')
         freq_scale = 1._WP
      case('HZ')
         freq_scale = 1._WP/(TWOPI*SQRT(bc%R_star**3/(bc%G*bc%M_star)))
      case('UHZ')
         freq_scale = 1.E6_WP/(TWOPI*SQRT(bc%R_star**3/(bc%G*bc%M_star)))
      case('ACOUSTIC_CUTOFF')
         call eval_cutoffs(bc, op, omega_cutoff_lo, omega_cutoff_hi)
         freq_scale = 1._WP/omega_cutoff_hi
      case('GRAVITY_CUTOFF')
         call eval_cutoffs(bc, op, omega_cutoff_lo, omega_cutoff_hi)
         freq_scale = 1._WP/omega_cutoff_lo
      case default
         $ABORT(Invalid freq_units)
      end select

      ! Finish

      return

    end function evol_freq_scale

!****

    function poly_freq_scale (freq_units) result (freq_scale)

      character(LEN=*), intent(in) :: freq_units
      real(WP)                     :: freq_scale

      ! Calculate the scale factor to convert a dimensionless angular
      ! frequency to a dimensioned frequency

      select case (freq_units)
      case ('NONE')
         freq_scale = 1._WP
      case default
         $ABORT(Invalid freq_units)
      end select

      ! Finish

      return

    end function poly_freq_scale

!****

    function hom_freq_scale (freq_units) result (freq_scale)

      character(LEN=*), intent(in) :: freq_units
      real(WP)                     :: freq_scale

      ! Calculate the scale factor to convert a dimensionless angular
      ! frequency to a dimensioned frequency

      select case (freq_units)
      case ('NONE')
         freq_scale = 1._WP
      case default
         $ABORT(Invalid freq_units)
      end select

      ! Finish

      return

    end function hom_freq_scale

  end function freq_scale
    
end module gyre_frontend
