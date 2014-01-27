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

  use gyre_coeffs
  use gyre_coeffs_evol
  use gyre_coeffs_poly
  use gyre_coeffs_hom
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_scanpar
  use gyre_atmos

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  character(LEN=64), save :: log_level_m

  ! Interfaces

  interface select_par
     module procedure select_par_np
     module procedure select_par_gp
     module procedure select_par_sp
  end interface select_par

  interface sprint
     module procedure sprint_i
  end interface sprint

  interface integrate
     module procedure integrate_r
     module procedure Integrate_c
  end interface integrate

  ! Access specifiers

  private

  public :: form_header
  public :: set_log_level
  public :: check_log_level
  public :: freq_scale
  public :: eval_cutoff_freqs
  public :: select_par
  public :: split_list
  public :: join_fmts
  public :: sprint
  public :: rjust
  public :: integrate

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

  function freq_scale (cf, op, x_o, freq_units)

    class(coeffs_t), intent(in)  :: cf
    type(oscpar_t), intent(in)   :: op
    real(WP), intent(in)         :: x_o
    character(LEN=*), intent(in) :: freq_units
    real(WP)                     :: freq_scale

    ! Calculate the scale factor to convert a dimensionless angular
    ! frequency to a dimensioned frequency

    select type (cf)
    class is (coeffs_evol_t)
       freq_scale = evol_freq_scale(cf, op, x_o, freq_units)
    class is (coeffs_poly_t)
       freq_scale = poly_freq_scale(freq_units)
    class is (coeffs_hom_t)
       freq_scale = hom_freq_scale(freq_units)
    class default
       $ABORT(Invalid cf type)
    end select

    ! Finish

    return

  contains

    function evol_freq_scale (cf, op, x_o, freq_units) result (freq_scale)

      class(coeffs_evol_t), intent(in) :: cf
      type(oscpar_t), intent(in)       :: op
      real(WP), intent(in)             :: x_o
      character(LEN=*), intent(in)     :: freq_units
      real(WP)                         :: freq_scale

      real(WP) :: omega_c_cutoff_lo
      real(WP) :: omega_c_cutoff_hi

      ! Calculate the scale factor to convert a dimensionless angular
      ! frequency to a dimensioned frequency

      select case(freq_units)
      case('NONE')
         freq_scale = 1._WP
      case('HZ')
         freq_scale = 1._WP/(TWOPI*SQRT(cf%R_star**3/(cf%G*cf%M_star)))
      case('UHZ')
         freq_scale = 1.E6_WP/(TWOPI*SQRT(cf%R_star**3/(cf%G*cf%M_star)))
      case('PER_DAY')
         freq_scale = 86400._WP/(TWOPI*SQRT(cf%R_star**3/(cf%G*cf%M_star)))
      case('ACOUSTIC_CUTOFF')
         call eval_cutoff_freqs(cf, op, x_o, omega_c_cutoff_lo, omega_c_cutoff_hi)
         freq_scale = 1._WP/omega_c_cutoff_hi
      case('GRAVITY_CUTOFF')
         call eval_cutoff_freqs(cf, op, x_o, omega_c_cutoff_lo, omega_c_cutoff_hi)
         freq_scale = 1._WP/omega_c_cutoff_lo
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

  subroutine eval_cutoff_freqs (cf, op, x_o, omega_c_cutoff_lo, omega_c_cutoff_hi)

    class(coeffs_t), intent(in) :: cf
    type(oscpar_t), intent(in)  :: op
    real(WP), intent(in)        :: x_o
    real(WP), intent(out)       :: omega_c_cutoff_lo
    real(WP), intent(out)       :: omega_c_cutoff_hi

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: c_1

     ! Evaluate the cutoff frequencies

     select case (op%outer_bound_type)
     case ('ZERO')
        omega_c_cutoff_lo = 0._WP
        omega_c_cutoff_hi = HUGE(0._WP)
     case ('DZIEM')
        omega_c_cutoff_lo = 0._WP
        omega_c_cutoff_hi = HUGE(0._WP)
     case ('UNNO')
        call eval_atmos_coeffs_unno(cf, x_o, V_g, As, c_1)
        call eval_atmos_cutoff_freqs(V_g, As, c_1, op%l, omega_c_cutoff_lo, omega_c_cutoff_hi)
     case('JCD')
        call eval_atmos_coeffs_jcd(cf, x_o, V_g, As, c_1)
        call eval_atmos_cutoff_freqs(V_g, As, c_1, op%l, omega_c_cutoff_lo, omega_c_cutoff_hi)
     case default
        $ABORT(Invalid outer_bound_type)
     end select

     ! Finish

     return

   end subroutine eval_cutoff_freqs

!****
   
   $define $SELECT_PAR $sub

   $local $SUFFIX $1
   $local $PAR_TYPE $2

   subroutine select_par_$SUFFIX (par, tag, par_sel, last)

     type($PAR_TYPE), intent(in)               :: par(:)
     character(LEN=*), intent(in)              :: tag
     type($PAR_TYPE), allocatable, intent(out) :: par_sel(:)
     logical, intent(in), optional             :: last

     logical :: last_
     integer :: i
     logical :: mask(SIZE(par))
     integer :: n_par_sel
     integer :: j

     if(PRESENT(last)) then
        last_ = last
     else
        last_ = .FALSE.
     endif

     ! Select all parameters whose tag_list matches tag

     mask_loop : do i = 1,SIZE(par)
        mask(i) = (par(i)%tag_list == '') .OR. &
                  (tag /= '' .AND. ANY(split_list(par(i)%tag_list, ',') == tag))
    end do mask_loop

    n_par_sel = COUNT(mask)

    allocate(par_sel(n_par_sel))

    j = 0

    select_loop : do i = 1,SIZE(par)
       if(mask(i)) then
          j = j + 1
          par_sel(j) = par(i)
       endif
    end do select_loop

    ! If necessary, shrink par_sel to contain only the last element

    if(last_ .AND. n_par_sel > 1) then
       par_sel = par_sel(1:1)
    endif

    ! Finish

    return

  end subroutine select_par_$SUFFIX

  $endsub

  $SELECT_PAR(np,numpar_t)
  $SELECT_PAR(gp,gridpar_t)
  $SELECT_PAR(sp,scanpar_t)

!****

  function split_list (list, delim) result (elems)

    character(LEN=*), intent(in)          :: list
    character(LEN=1), intent(in)          :: delim
    character(LEN=LEN(list)), allocatable :: elems(:)

    character(LEN=LEN(list)) :: list_
    integer                  :: d
    integer                  :: n
    integer                  :: j
    
    ! Split the delimited list into an array of elements
 
    d = 16

    allocate(elems(d))

    n = 0

    ! Repeatedly split on delimiters

    list_ = list

    split_loop : do

       if(list_ == '') exit split_loop

       j = INDEX(list_, delim)

       if(j <= 0) then
          n = n + 1
          elems(n) = list_
          exit split_loop
       endif

       n = n + 1

       ! Chop out the element

       elems(n) = list_(:j-1)
       list_ = list_(j+1:)

       ! If necessary, expand the array
          
       if(n >= d) then
          d = 2*d
          call reallocate(elems, [d])
       end if

    end do split_loop

    ! Reallocate elems to the correct length

    call reallocate(elems, [n])

    ! Finish

    return

  end function split_list

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

!****

  $define $INTEGRATE $sub

  $local $SUFFIX $1
  $local $TYPE $2

  function integrate_$SUFFIX (x, y) result (int_y)

    real(WP), intent(in)  :: x(:)
    $TYPE(WP), intent(in) :: y(:)
    $TYPE(WP)             :: int_y

    integer :: n

    $CHECK_BOUNDS(SIZE(y),SIZE(x))

    ! Integrate y(x) using trapezoidal quadrature

    n = SIZE(x)

    int_y = SUM(0.5_WP*(y(2:) + y(:n-1))*(x(2:) - x(:n-1)))

    ! Finish

    return

  end function integrate_$SUFFIX

  $endsub

  $INTEGRATE(r,real)
  $INTEGRATE(c,complex)

end module gyre_util
