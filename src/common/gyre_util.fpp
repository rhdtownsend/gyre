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
  use core_parallel
  use core_memory

  use gyre_atmos
  use gyre_constants
  use gyre_evol_model
  use gyre_gridpar
  use gyre_hom_model
  use gyre_model
  use gyre_modepar
  use gyre_numpar
  use gyre_oscpar
  use gyre_poly_model
  use gyre_scanpar
  use gyre_scons_model
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  character(64), save :: log_level_m

  ! Interfaces

  interface omega_from_freq
     module procedure omega_from_freq_r_
     module procedure omega_from_freq_c_
  end interface omega_from_freq

interface freq_from_omega
     module procedure freq_from_omega_r_
     module procedure freq_from_omega_c_
  end interface freq_from_omega

  interface select_par
     module procedure select_par_op_
     module procedure select_par_np_
     module procedure select_par_gp_
     module procedure select_par_sp_
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
  public :: omega_from_freq
  public :: freq_from_omega
  public :: eval_cutoff_freqs
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

  $define $OMEGA_FROM_FREQ $sub

  $local $T $1
  $local $TYPE $2

  function omega_from_freq_${T}_ (freq, ml, mp, op, x_i, x_o, freq_units, freq_frame) result (omega)

    $TYPE(WP), intent(in)               :: freq
    class(model_t), pointer, intent(in) :: ml
    type(modepar_t), intent(in)         :: mp
    type(oscpar_t), intent(in)          :: op
    real(WP), intent(in)                :: x_i
    real(WP), intent(in)                :: x_o
    character(*), intent(in)            :: freq_units
    character(*), intent(in)            :: freq_frame
    $TYPE(WP)                           :: omega

    $TYPE(WP)                      :: omega_l
    class(${T}_rot_t), allocatable :: rt
    real(WP)                       :: omega_cutoff_lo
    real(WP)                       :: omega_cutoff_hi

    ! Calculate the dimensionless inertial-frame frequency omega from
    ! the dimensioned local-frame frequency freq

    ! First calculate the dimensionless frequency in the local frame

    select type (ml)

    class is (evol_model_t)

       select case(freq_units)
       case('NONE')
          omega_l = freq
       case('HZ')
          omega_l = freq*(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))
       case('UHZ')
          omega_l = freq*(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))/1E6_WP
       case('PER_DAY')
          omega_l = freq*(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))/86400._WP
       case('ACOUSTIC_CUTOFF')
          call eval_cutoff_freqs(ml, mp, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
          omega_l = freq*omega_cutoff_hi
       case('GRAVITY_CUTOFF')
          call eval_cutoff_freqs(ml, mp, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
          omega_l = freq*omega_cutoff_lo
       case default
          $ABORT(Invalid freq_units)
       end select

    class is (scons_model_t)

       select case(freq_units)
       case('NONE')
          omega_l = freq
       case('HZ')
          omega_l = freq*(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))
       case('UHZ')
          omega_l = freq*(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))/1E6_WP
       case('PER_DAY')
          omega_l = freq*(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))/86400._WP
       case('ACOUSTIC_CUTOFF')
          call eval_cutoff_freqs(ml, mp, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
          omega_l = freq*omega_cutoff_hi
       case('GRAVITY_CUTOFF')
          call eval_cutoff_freqs(ml, mp, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
          omega_l = freq*omega_cutoff_lo
       case default
          $ABORT(Invalid freq_units)
       end select

    class is (poly_model_t)

       select case (freq_units)
       case ('NONE')
          omega_l = freq
       case default
         $ABORT(Invalid freq_units)
      end select

   class is (hom_model_t)

       select case (freq_units)
       case ('NONE')
          omega_l = freq
       case default
         $ABORT(Invalid freq_units)
      end select

    class default

       $ABORT(Invalid ml type)

    end select

    ! Now convert to the inertial frame

    allocate(rt, SOURCE=${T}_rot_t(ml, mp, op))

    select case (freq_frame)
    case ('INERTIAL')
       omega = omega_l
    case ('COROT_I')
       omega = rt%omega(x_i, omega_l)
    case ('COROT_O')
       omega = rt%omega(x_o, omega_l)
    case default
       $ABORT(Invalid freq_frame)
    end select

    ! Finish

    return

  end function omega_from_freq_${T}_

  $endsub

  $OMEGA_FROM_FREQ(r,real)
  $OMEGA_FROM_FREQ(c,complex)

!****

  $define $FREQ_FROM_OMEGA $sub

  $local $T $1
  $local $TYPE $2

  function freq_from_omega_${T}_ (omega, ml, mp, op, x_i, x_o, freq_units, freq_frame) result (freq)

    $TYPE(WP), intent(in)               :: omega
    class(model_t), pointer, intent(in) :: ml
    type(modepar_t), intent(in)         :: mp
    type(oscpar_t), intent(in)          :: op
    real(WP), intent(in)                :: x_i
    real(WP), intent(in)                :: x_o
    character(*), intent(in)            :: freq_units
    character(*), intent(in)            :: freq_frame
    $TYPE(WP)                           :: freq

    $TYPE(WP)                      :: omega_l
    class(${T}_rot_t), allocatable :: rt
    real(WP)                       :: omega_cutoff_lo
    real(WP)                       :: omega_cutoff_hi

    ! Calculate the dimensioned local-frame frequency freq from the
    ! dimensionless inertial-frame frequency omega

    ! First convert from the inertial frame

    allocate(rt, SOURCE=${T}_rot_t(ml, mp, op))

    select case (freq_frame)
    case ('INERTIAL')
       omega_l = omega
    case ('COROT_I')
       omega_l = rt%omega_c(x_i, omega)
    case ('COROT_O')
       omega_l = rt%omega_c(x_o, omega)
    case default
       $ABORT(Invalid freq_frame)
    end select

    ! Now calculate the dimensionless frequency in the local frame

    select type (ml)

    class is (evol_model_t)

       select case(freq_units)
       case('NONE')
          freq = omega_l
       case('HZ')
          freq = omega_l/(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))
       case('UHZ')
          freq = omega_l/(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))*1E6_WP
       case('PER_DAY')
          freq = omega_l/(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))*86400._WP
       case('ACOUSTIC_CUTOFF')
          call eval_cutoff_freqs(ml, mp, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
          freq = omega_l/omega_cutoff_hi
       case('GRAVITY_CUTOFF')
          call eval_cutoff_freqs(ml, mp, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
          freq = omega_l/omega_cutoff_lo
       case default
          $ABORT(Invalid freq_units)
       end select

    class is (scons_model_t)

       select case(freq_units)
       case('NONE')
          freq = omega_l
       case('HZ')
          freq = omega_l/(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))
       case('UHZ')
          freq = omega_l/(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))*1E6_WP
       case('PER_DAY')
          freq = omega_l/(TWOPI*SQRT(ml%R_star**3/(G_GRAVITY*ml%M_star)))*86400._WP
       case('ACOUSTIC_CUTOFF')
          call eval_cutoff_freqs(ml, mp, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
          freq = omega_l/omega_cutoff_hi
       case('GRAVITY_CUTOFF')
          call eval_cutoff_freqs(ml, mp, op, x_o, omega_cutoff_lo, omega_cutoff_hi)
          freq = omega_l/omega_cutoff_lo
       case default
          $ABORT(Invalid freq_units)
       end select

    class is (poly_model_t)

       select case (freq_units)
       case ('NONE')
          freq = omega_l
       case default
         $ABORT(Invalid freq_units)
      end select

    class is (hom_model_t)

       select case (freq_units)
       case ('NONE')
          freq = omega_l
       case default
         $ABORT(Invalid freq_units)
      end select

    class default

       $ABORT(Invalid ml type)

    end select

    ! Finish

    return

  end function freq_from_omega_${T}_

  $endsub

  $FREQ_FROM_OMEGA(r,real)
  $FREQ_FROM_OMEGA(c,complex)

!****

  subroutine eval_cutoff_freqs (ml, mp, op, x_o, omega_cutoff_lo, omega_cutoff_hi)

    class(model_t), intent(in)  :: ml
    type(modepar_t), intent(in) :: mp
    type(oscpar_t), intent(in)  :: op
    real(WP), intent(in)        :: x_o
    real(WP), intent(out)       :: omega_cutoff_lo
    real(WP), intent(out)       :: omega_cutoff_hi

    real(WP)      :: V_g
    real(WP)      :: As
    real(WP)      :: c_1
    logical, save :: warned = .FALSE.

     ! Evaluate the cutoff frequencies

     select case (op%outer_bound_type)

     case ('ZERO')

        omega_cutoff_lo = 0._WP
        omega_cutoff_hi = HUGE(0._WP)

     case ('DZIEM')

        omega_cutoff_lo = 0._WP
        omega_cutoff_hi = HUGE(0._WP)

     case ('UNNO')

        call eval_atmos_coeffs_unno(ml, x_o, V_g, As, c_1)
        call eval_atmos_cutoff_freqs(V_g, As, c_1, mp%l*(mp%l+1._WP), omega_cutoff_lo, omega_cutoff_hi)

     case('JCD')

        call eval_atmos_coeffs_jcd(ml, x_o, V_g, As, c_1)
        call eval_atmos_cutoff_freqs(V_g, As, c_1, mp%l*(mp%l+1._WP), omega_cutoff_lo, omega_cutoff_hi)

     case default

        $ABORT(Invalid outer_bound_type)

     end select

     if (.not. warned) then
        $WARN(WARNING: Cutoff frequencies do not account for rotation effects)
        warned = .TRUE.
     endif

     ! Finish

     return

   end subroutine eval_cutoff_freqs

!****
   
   $define $SELECT_PAR $sub

   $local $INFIX $1
   $local $PAR_TYPE $2

   subroutine select_par_${INFIX}_ (par, tag, par_sel, last)

     type($PAR_TYPE), intent(in)               :: par(:)
     character(*), intent(in)                  :: tag
     type($PAR_TYPE), allocatable, intent(out) :: par_sel(:)
     logical, optional, intent(in)             :: last

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

  end subroutine select_par_${INFIX}_

  $endsub

  $SELECT_PAR(op,oscpar_t)
  $SELECT_PAR(np,numpar_t)
  $SELECT_PAR(gp,gridpar_t)
  $SELECT_PAR(sp,scanpar_t)

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

  function integrate_${INFIX}_ (x, y) result (int_y)

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
