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

  ! Access specifiers

  private

  public :: write_header
  public :: freq_scale
  public :: split_item_list

contains

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

end module gyre_util
