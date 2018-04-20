! Module   : gyre_context 
! Purpose  : context (model, data, functions) for solving pulsation equations
!
! Copyright 2017 Rich Townsend
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

module gyre_context

  ! Uses

  use core_kinds
  use core_order

  use gyre_freq
  use gyre_grid
  use gyre_grid_par
  use gyre_interp
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_rot_factory
  use gyre_state
  
  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: context_t
     class(model_t), pointer, public :: ml => null()
     class(c_rot_t), allocatable     :: rt
     type(point_t), public           :: pt_i
     type(point_t), public           :: pt_o
     integer                         :: m
     logical                         :: complex_lambda
     type(c_interp_t)                :: in_eps_rho
     type(c_interp_t)                :: in_eps_T
   contains
     private
     procedure       :: omega_c_r_
     procedure       :: omega_c_c_
     generic, public :: omega_c => omega_c_r_, omega_c_c_
     procedure       :: lambda_r_
     procedure       :: lambda_c_
     generic, public :: lambda => lambda_r_, lambda_c_
     procedure       :: l_e_r_
     procedure       :: l_e_c_
     generic, public :: l_e => l_e_r_, l_e_c_
     procedure       :: eps_rho_r_
     procedure       :: eps_rho_c_
     generic, public :: eps_rho => eps_rho_r_, eps_rho_c_
     procedure       :: eps_T_r_
     procedure       :: eps_T_c_
     generic, public :: eps_T => eps_T_r_, eps_T_c_
  end type context_t

  ! Interfaces

  interface context_t
     module procedure context_t_
  end interface context_t

  ! Access specifiers

  private

  public :: context_t

  ! Procedures

contains

  function context_t_ (ml, gr_p, md_p, os_p) result (cx)

    class(model_t), pointer, intent(in) :: ml
    type(grid_par_t), intent(in)        :: gr_p
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(context_t)                     :: cx

    type(grid_t) :: gr

    ! Construct the context_t

    cx%ml => ml

    allocate(cx%rt, SOURCE=c_rot_t(md_p, os_p))

    gr = grid_t(ml%grid(), gr_p%x_i, gr_p%x_o)
    cx%pt_i = gr%pt_i()
    cx%pt_o = gr%pt_o()
 
    cx%m = md_p%m

    cx%complex_lambda = os_p%complex_lambda

    ! If necessary, initialize the nuclear burning partials

    select case (os_p%deps_scheme)
    case ('FILE')
       call read_deps_(ml, md_p, os_p, cx%pt_i, cx%pt_o, cx%in_eps_rho, cx%in_eps_T)
    end select

    ! Finish

    return

  end function context_t_

  !****

  function omega_c_r_ (this, Omega_rot, st) result (omega_c)

    class(context_t), intent(in) :: this
    real(WP), intent(in)         :: Omega_rot
    class(r_state_t), intent(in) :: st
    real(WP)                     :: omega_c

    ! Evaluate the co-rotating frequency (real)

    omega_c = omega_corot(st%omega, Omega_rot, this%m)

    ! Finish

    return

  end function omega_c_r_

  !****

  function omega_c_c_ (this, Omega_rot, st) result (omega_c)

    class(context_t), intent(in) :: this
    real(WP), intent(in)         :: Omega_rot
    class(c_state_t), intent(in) :: st
    complex(WP)                  :: omega_c

    ! Evaluate the co-rotating frequency (complex)

    omega_c = omega_corot(st%omega, Omega_rot, this%m)

    ! Finish

    return

  end function omega_c_c_

  !****

  function lambda_r_ (this, Omega_rot, st) result (lambda)

    class(context_t), intent(in) :: this
    real(WP), intent(in)         :: Omega_rot
    class(r_state_t), intent(in) :: st
    real(WP)                     :: lambda

    ! Evaluate the angular eigenvalue (real)

    lambda = REAL(this%rt%lambda(Omega_rot, CMPLX(st%omega, KIND=WP)))

    ! Finish

    return

  end function lambda_r_

  !****

  function lambda_c_ (this, Omega_rot, st) result (lambda)

    class(context_t), intent(in) :: this
    real(WP), intent(in)         :: Omega_rot
    class(c_state_t), intent(in) :: st
    complex(WP)                  :: lambda

    ! Evaluate the angular eigenvalue (complex)

    if (this%complex_lambda) then

       lambda = this%rt%lambda(Omega_rot, st%omega)

    else

       lambda = this%rt%lambda(Omega_rot, CMPLX(st%omega_r, KIND=WP))

    endif

    ! Finish

    return

  end function lambda_c_

  !****

  function l_e_r_ (this, Omega_rot, st) result (l_e)

    class(context_t), intent(in) :: this
    real(WP), intent(in)         :: Omega_rot
    class(r_state_t), intent(in) :: st
    real(WP)                     :: l_e

    ! Evaluate the effective harmonic degree (real)

    l_e = REAL(this%rt%l_e(Omega_rot, CMPLX(st%omega, KIND=WP)))

    ! Finish

    return

  end function l_e_r_

  !****

  function l_e_c_ (this, Omega_rot, st) result (l_e)

    class(context_t), intent(in) :: this
    real(WP), intent(in)         :: Omega_rot
    class(c_state_t), intent(in) :: st
    complex(WP)                  :: l_e

    ! Evaluate the effective harmonic degree (complex)

    if (this%complex_lambda) then

       l_e = this%rt%l_e(Omega_rot, st%omega)

    else

       l_e = this%rt%l_e(Omega_rot, CMPLX(st%omega_r, KIND=WP))

    endif

    ! Finish

    return

  end function l_e_c_
  
  !****

  subroutine read_deps_ (ml, md_p, os_p, pt_i, pt_o, in_eps_rho, in_eps_T)
    
    class(model_t), pointer, intent(in) :: ml
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(point_t), intent(in)           :: pt_i
    type(point_t), intent(in)           :: pt_o
    type(c_interp_t), intent(out)       :: in_eps_rho
    type(c_interp_t), intent(out)       :: in_eps_T

    ! Read epsilon partials data from a file

    select case (os_p%deps_file_format)
    case ('WOLF')
       call read_deps_wolf_(ml, md_p, os_p, pt_i, pt_o, in_eps_rho, in_eps_T)
    case default
       $ABORT(Invalid deps_file_format)
    end select

    ! Finish

    return

  end subroutine read_deps_

  !****

  subroutine read_deps_wolf_ (ml, md_p, os_p, pt_i, pt_o, in_eps_rho, in_eps_T)

    class(model_t), pointer, intent(in) :: ml
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(point_t), intent(in)           :: pt_i
    type(point_t), intent(in)           :: pt_o
    type(c_interp_t), intent(out)       :: in_eps_rho
    type(c_interp_t), intent(out)       :: in_eps_T

    integer                  :: unit
    integer                  :: n
    real(WP), allocatable    :: omega(:)
    real(WP), allocatable    :: A_norm(:)
    real(WP), allocatable    :: B_norm(:)
    real(WP), allocatable    :: A_phase(:)
    real(WP), allocatable    :: B_phase(:)
    integer                  :: i
    real(WP)                 :: period
    complex(WP), allocatable :: eps_rho(:)
    complex(WP), allocatable :: eps_T(:)
    integer, allocatable     :: j(:)
    
    ! Read epsilon partials data from a file in the format used in the
    ! preparation of Wolf, Townsend & Bildsten (2017)

    ! Open the file and skip the header line

    open(NEWUNIT=unit, FILE=os_p%deps_file, STATUS='OLD')

    read(unit, *)

    ! Count lines

    n = 0

    count_loop : do
       read(unit, *, END=100)
       n = n + 1
    end do count_loop

100 continue

    rewind(unit)

    read(unit, *)

    ! Read data

    allocate(omega(n))

    allocate(A_norm(n))
    allocate(B_norm(n))

    allocate(A_phase(n))
    allocate(B_phase(n))

    read_loop : do i = 1, n
       read(unit, *) period, A_norm(i), B_norm(i), A_phase(i), B_phase(i)
       omega(i) = omega_from_freq(1._WP/period, ml, pt_i, pt_o, 'HZ', 'INERTIAL', md_p, os_p)
    end do read_loop

    close(unit)

    ! Reverse the phase sign because Wolf et al assume a time
    ! depenence of exp(i omega t), rather than GYRE's exp(-i omega t)

    A_phase = -A_phase
    B_phase = -B_phase

    ! Calculate complex eps_* terms

    eps_rho = A_norm*EXP((0._WP,1._WP)*A_phase)
    eps_T = B_norm*EXP((0._WP,1._WP)*A_phase)

    ! Sort in order of increasing omega

    j = sort_indices(omega)

    omega = omega(j)
    eps_rho = eps_rho(j)
    eps_T = eps_T(j)

    ! Set up interpolating splines

    in_eps_rho = c_interp_t(omega, eps_rho, 'SPLINE')
    in_eps_T = c_interp_t(omega, eps_T, 'SPLINE')

    ! Finish

    return

  end subroutine read_deps_wolf_

  !****

  function eps_rho_r_ (this, st) result (eps_rho)

    class(context_t), intent(in) :: this
    class(r_state_t), intent(in) :: st
    complex(WP)                  :: eps_rho

    real(WP) :: omega_min
    real(WP) :: omega_max
    real(WP) :: omega

    ! Evaluate the eps_rho derivative (real)

    omega_min = this%in_eps_rho%x_min()
    omega_max = this%in_eps_rho%x_min()

    omega = MIN(MAX(st%omega, omega_min), omega_max)

    eps_rho = this%in_eps_rho%f(omega)

    ! Finish

    return

  end function eps_rho_r_

  !****

  function eps_rho_c_ (this, st) result (eps_rho)

    class(context_t), intent(in) :: this
    class(c_state_t), intent(in) :: st
    complex(WP)                  :: eps_rho

    real(WP) :: omega_min
    real(WP) :: omega_max
    real(WP) :: omega

    ! Evaluate the eps_rho derivative (complex)

    omega_min = this%in_eps_rho%x_min()
    omega_max = this%in_eps_rho%x_min()

    omega = MIN(MAX(st%omega_r, omega_min), omega_max)

    eps_rho = this%in_eps_rho%f(omega)

    ! Finish

    return

  end function eps_rho_c_

  !****

  function eps_T_r_ (this, st) result (eps_T)

    class(context_t), intent(in) :: this
    class(r_state_t), intent(in) :: st
    complex(WP)                  :: eps_T

    real(WP) :: omega_min
    real(WP) :: omega_max
    real(WP) :: omega

    ! Evaluate the eps_T derivative (real)

    omega_min = this%in_eps_rho%x_min()
    omega_max = this%in_eps_rho%x_min()

    omega = MIN(MAX(st%omega, omega_min), omega_max)

    eps_T = this%in_eps_T%f(omega)

    ! Finish

    return

  end function eps_T_r_

  !****

  function eps_T_c_ (this, st) result (eps_T)

    class(context_t), intent(in) :: this
    class(c_state_t), intent(in) :: st
    complex(WP)                  :: eps_T

    real(WP) :: omega_min
    real(WP) :: omega_max
    real(WP) :: omega

    ! Evaluate the eps_T derivative (real)

    omega_min = this%in_eps_rho%x_min()
    omega_max = this%in_eps_rho%x_min()

    omega = MIN(MAX(st%omega_r, omega_min), omega_max)

    eps_T = this%in_eps_T%f(omega)

    ! Finish

    return

  end function eps_T_c_

end module gyre_context
