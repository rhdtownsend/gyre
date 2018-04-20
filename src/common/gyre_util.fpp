! Module   : gyre_util
! Purpose  : miscellaneous utility routines
!
! Copyright 2013-2017 Rich Townsend
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
  use core_parallel
  use core_memory

  use gyre_constants
  use gyre_grid_par
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_out_par
  use gyre_scan_par

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  character(64), save :: log_level_m

  ! Interfaces

  interface select_par
     module procedure select_par_gr_1_
     module procedure select_par_nm_1_
     module procedure select_par_os_1_
     module procedure select_par_sc_1_
     module procedure select_par_gr_v_
     module procedure select_par_nm_v_
     module procedure select_par_os_v_
     module procedure select_par_sc_v_
  end interface select_par

  interface sprint
     module procedure sprint_
  end interface sprint

  interface integrate
     module procedure integrate_r_
     module procedure Integrate_c_
  end interface integrate

  interface integral
     module procedure integral_r_
     module procedure Integral_c_
  end interface integral

  ! Access specifiers

  private

  public :: form_header
  public :: set_log_level
  public :: check_log_level
  public :: select_par
  public :: split_list
  public :: join_fmts
  public :: sprint
  public :: rjust
  public :: phase
  public :: integrate
  public :: integral

  ! Procedures

contains

  function form_header (header, underchar)

    character(*), intent(in)           :: header
    character(*), optional, intent(in) :: underchar
    character(:), allocatable          :: form_header

    ! Format the header string

    if(PRESENT(underchar)) then

       if(underchar == '') then

          form_header = TRIM(header) // NEW_LINE('') // &
                        REPEAT(' ', LEN(header)) // NEW_LINE('')

       else

          form_header = TRIM(header) // NEW_LINE('') // &
                        REPEAT(underchar, LEN(header)/LEN(underchar)) // NEW_LINE('')

       endif

    else
       
       form_header = TRIM(header) // NEW_LINE('')
       
    endif

    ! Finish

    return

  end function form_header

  !****

  subroutine set_log_level (log_level)

    character(*), intent(in) :: log_level
    
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

    character(*), intent(in)      :: log_level
    integer, optional, intent(in) :: rank
    logical                       :: check_log_level

    integer :: rank_

    if(PRESENT(rank)) then
       rank_ = rank
    else
       rank_ = 0
    endif

    ! Check whether we should write log output

    if (MPI_RANK == rank_) then
       
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
   
  $define $SELECT_PAR_1 $sub

  $local $INFIX $1
  $local $PAR_TYPE $2

  subroutine select_par_${INFIX}_1_ (par, tag, par_sel)

    type($PAR_TYPE), intent(in)   :: par(:)
    character(*), intent(in)      :: tag
    type($PAR_TYPE), intent(out)  :: par_sel

    type($PAR_TYPE), allocatable :: par_sel_(:)

    ! Select the last parameter whose tag_list matches tag

    call select_par(par, tag, par_sel_)

    par_sel = par_sel_(SIZE(par_sel_))

    ! Finish
     
    return

  end subroutine select_par_${INFIX}_1_

  $endsub

  $SELECT_PAR_1(gr,grid_par_t)
  $SELECT_PAR_1(nm,num_par_t)
  $SELECT_PAR_1(os,osc_par_t)
  $SELECT_PAR_1(sc,scan_par_t)

  !****
   
  $define $SELECT_PAR_V $sub

  $local $INFIX $1
  $local $VAR $2
  $local $PAR_TYPE $3
  $local $PAR_NAME $4

  subroutine select_par_${INFIX}_v_ (par, $VAR, par_sel)

    type($PAR_TYPE), intent(in)               :: par(:)
    character(*), intent(in)                  :: $VAR
    type($PAR_TYPE), allocatable, intent(out) :: par_sel(:)

    integer :: i
    logical :: mask(SIZE(par))
    integer :: n_par_sel
    integer :: j

    ! Select all parameters whose $VAR_list matches $VAR

    mask_loop : do i = 1,SIZE(par)
       mask(i) = (par(i)%${VAR}_list == '') .OR. &
                 ($VAR /= '' .AND. ANY(split_list(par(i)%${VAR}_list, ',') == $VAR))
    end do mask_loop

    n_par_sel = COUNT(mask)

    $ASSERT(n_par_sel >= 1,No matching $PAR_NAME namelists)
    
    allocate(par_sel(n_par_sel))

    j = 0

    select_loop : do i = 1,SIZE(par)
       if (mask(i)) then
          j = j + 1
          par_sel(j) = par(i)
       endif
    end do select_loop

    ! Finish
     
    return

  end subroutine select_par_${INFIX}_v_

  $endsub

  $SELECT_PAR_V(gr,tag,grid_par_t,&grid)
  $SELECT_PAR_V(nm,tag,num_par_t,&num)
  $SELECT_PAR_V(os,tag,osc_par_t,&osc)
  $SELECT_PAR_V(sc,tag,scan_par_t,&scan)

  !****

  function split_list (list, delim) result (elems)

    character(*), intent(in)          :: list
    character(1), intent(in)          :: delim
    character(LEN(list)), allocatable :: elems(:)

    character(LEN(list)) :: list_
    integer              :: d
    integer              :: n
    integer              :: j
    
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
          elems(n) = ADJUSTL(list_)
          exit split_loop
       endif

       n = n + 1

       ! Chop out the element

       elems(n) = ADJUSTL(list_(:j-1))
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
    
    character(*), intent(in)  :: fmts(:)
    integer, intent(in)       :: n(:)
    character(:), allocatable :: fmt

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

  function sprint_ (i) result (a)

    integer, intent(in)       :: i
    character(:), allocatable :: a

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

    allocate(character(n)::a)

    ! Do the conversion

    write(a, 100) i
100 format(I0)

    ! Finish

    return

  end function sprint_

!****

  function rjust (a, n) result (a_just)

    character(*), intent(in) :: a
    integer, intent(in)      :: n
    character(n)             :: a_just

    ! Right-justify a in a field width of n

    a_just = REPEAT(' ', MAX(n-LEN_TRIM(a), 0))//a

    ! Finish

    return

  end function rjust

!****

  function phase (z)

    complex(WP), intent(in) :: z
    real(WP)                :: phase

    ! Calculate the phase (in radians) of the complex number z

    phase = ATAN2(AIMAG(z), REAL(z))

    ! Finish

    return

  end function phase

!****

  $define $INTEGRATE $sub

  $local $INFIX $1
  $local $TYPE $2

  function integrate_${INFIX}_ (x, y, mask) result (int_y)

    real(WP), intent(in)          :: x(:)
    $TYPE(WP), intent(in)         :: y(:)
    logical, optional, intent(in) :: mask(:)
    $TYPE(WP)                     :: int_y

    integer :: n

    $CHECK_BOUNDS(SIZE(y),SIZE(x))

    if (PRESENT(mask)) then
       $CHECK_BOUNDS(SIZE(mask),SIZE(x))
    endif

    ! Integrate y(x) using trapezoidal quadrature, applying the
    ! optional mask

    n = SIZE(x)

    if (PRESENT(mask)) then
       int_y = SUM(0.5_WP*(y(2:) + y(:n-1))*(x(2:) - x(:n-1)), MASK=mask(2:) .AND. mask(:n-1))
    else
       int_y = SUM(0.5_WP*(y(2:) + y(:n-1))*(x(2:) - x(:n-1)))
    endif

    ! Finish

    return

  end function integrate_${INFIX}_

  $endsub

  $INTEGRATE(r,real)
  $INTEGRATE(c,complex)

!****

  $define $INTEGRAL $sub

  $local $INFIX $1
  $local $TYPE $2

  function integral_${INFIX}_ (x, y) result (int_y)

    real(WP), intent(in)  :: x(:)
    $TYPE(WP), intent(in) :: y(:)
    $TYPE(WP)             :: int_y(SIZE(x))

    integer :: n
    integer :: i

    $CHECK_BOUNDS(SIZE(y),SIZE(x))

    ! Calculate the integral of y(x) using trapezoidal quadrature

    n = SIZE(x)

    int_y(1) = 0._WP

    int_loop : do i = 2, n
       int_y(i) = int_y(i-1) + 0.5_WP*(y(i) + y(i-1))*(x(i) - x(i-1))
    end do int_loop

    ! Finish

    return

  end function integral_${INFIX}_

  $endsub

  $INTEGRAL(r,real)
  $INTEGRAL(c,complex)

end module gyre_util
