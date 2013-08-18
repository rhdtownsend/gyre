! Program  : gyre_util
! Purpose  : miscellaneous utility routines
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

module gyre_util

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel
  use core_memory

  use gyre_base_coeffs
  use gyre_evol_base_coeffs
  use gyre_poly_base_coeffs
  use gyre_hom_base_coeffs
  use gyre_oscpar
  use gyre_ad_bound

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  character(LEN=64), save :: log_level_m

  ! Interfaces

  interface sprint
     module procedure sprint_i
  end interface sprint

  ! Access specifiers

  private

  public :: form_header
  public :: set_log_level
  public :: check_log_level
  public :: freq_scale
  public :: split_item_list
  public :: join_fmts
  public :: sprint
  public :: rjust

contains

  function form_header (header, underchar)

    character(LEN=*), intent(in)           :: header
    character(LEN=*), intent(in), optional :: underchar
    character(LEN=:), allocatable          :: form_header

    ! Format the header string

    if(PRESENT(underchar)) then

       if(underchar == '') then

          form_header = NEW_LINE('') // &
                        TRIM(header) // NEW_LINE('') // &
                        REPEAT(' ', LEN(header)) // NEW_LINE('')

       else

          form_header = NEW_LINE('') // &
                        TRIM(header) // NEW_LINE('') // &
                        REPEAT(underchar, LEN(header)/LEN(underchar)) // NEW_LINE('')

       endif

    else
       
       form_header = NEW_LINE('') // &
                     TRIM(header) // NEW_LINE('')
       
    endif

    ! Finish

    return

  end function form_header

!****

  subroutine set_log_level (log_level)

    character(LEN=*), intent(in) :: log_level
    
    ! Set the log level

    select case (log_level)
    case ('DEBUG')
    case ('INFO')
    case ('WARN')
    case default
       $ABORT(Invalid log_level)
    end select

    log_level_m = log_level

    ! Finish

    return

  end subroutine set_log_level

!****

  function check_log_level (log_level, rank)

    character(LEN=*), intent(in)  :: log_level
    integer, intent(in), optional :: rank
    logical                       :: check_log_level

    integer :: rank_

    if(PRESENT(rank)) then
       rank_ = rank
    else
       rank_ = 0
    endif

    ! Check whether we should write log output

    if(MPI_RANK == rank_) then
       
       select case (log_level)
       case ('DEBUG')
          check_log_level = log_level_m == 'DEBUG'
       case ('INFO')
          check_log_level = log_level_m == 'INFO' .OR. &
                            log_level_m == 'DEBUG'
       case ('WARN')
          check_log_level = log_level_m == 'WARN' .OR. &
                            log_level_m == 'INFO' .OR. &
                            log_level_m == 'DEBUG'
       case default
          $ABORT(Invalid log_level)
       end select

    else

       check_log_level = .FALSE.

    endif

    ! Finish

    return

  end function check_log_level

!****

  function freq_scale (bc, op, x_o, freq_units)

    class(base_coeffs_t), intent(in) :: bc
    type(oscpar_t), intent(in)       :: op
    real(WP), intent(in)             :: x_o
    character(LEN=*), intent(in)     :: freq_units
    real(WP)                         :: freq_scale

    ! Calculate the scale factor to convert a dimensionless angular
    ! frequency to a dimensioned frequency

    select type (bc)
    class is (evol_base_coeffs_t)
       freq_scale = evol_freq_scale(bc, op, x_o, freq_units)
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

    function evol_freq_scale (bc, op, x_o, freq_units) result (freq_scale)

      class(evol_base_coeffs_t), intent(in) :: bc
      type(oscpar_t), intent(in)            :: op
      real(WP), intent(in)                  :: x_o
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
      case('PER_DAY')
         freq_scale = 86400._WP/(TWOPI*SQRT(bc%R_star**3/(bc%G*bc%M_star)))
      case('ACOUSTIC_CUTOFF')
         call eval_cutoffs(bc, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
         freq_scale = 1._WP/omega_cutoff_hi
      case('GRAVITY_CUTOFF')
         call eval_cutoffs(bc, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
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

  function join_fmts (fmts, n) result (fmt)
    
    character(LEN=*), intent(in)  :: fmts(:)
    integer, intent(in)           :: n(:)
    character(LEN=:), allocatable :: fmt

    integer :: i

    $CHECK_BOUNDS(SIZE(n),SIZE(fmts))

    ! Join format strings with the appropriate repeat counts

    if(SUM(n) > 0) then

       do i = 1, SIZE(fmts)

          if(ALLOCATED(fmt)) then
             fmt = fmt//','//sprint(n(i))//fmts(i)
          else
             fmt = sprint(n(i))//fmts(i)
          endif

       end do

    else

       fmt = ''

    endif

    ! Add wrap-around parens

    fmt = '('//fmt//')'

    ! Finish

    return

  end function join_fmts

!****

  function sprint_i (i) result (a)

    integer, intent(in)           :: i
    character(LEN=:), allocatable :: a

    integer :: n

    ! Print an integer into a character

    ! First, determine the length

    if(i > 0) then
       n = FLOOR(LOG10(REAL(i))) + 1
    elseif(i < 0) then
       n = FLOOR(LOG10(REAL(ABS(i)))) + 2
    else
       n = 1
    endif

    allocate(character(LEN=n)::a)

    ! Do the conversion

    write(a, 100) i
100 format(I0)

    ! Finish

    return

  end function sprint_i

!****

  function rjust (a, n) result (a_just)

    character(LEN=*), intent(in) :: a
    integer, intent(in)          :: n
    character(LEN=n)             :: a_just

    ! Right-justify a in a field width of n

    a_just = REPEAT(' ', MAX(n-LEN_TRIM(a), 0))//a

    ! Finish

    return

  end function rjust

end module gyre_util
