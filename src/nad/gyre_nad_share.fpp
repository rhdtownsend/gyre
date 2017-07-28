! Module   : gyre_nad_share
! Purpose  : nonadiabatic shared data
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

module gyre_nad_share

  ! Uses

  use core_kinds
  use core_order

  use gyre_freq
  use gyre_interp
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: nad_share_t
     private
     class(model_t), pointer, public :: ml
     class(c_rot_t), allocatable     :: rt
     real(WP), allocatable           :: omega_r
     real(WP)                        :: Omega_rot_i
     logical                         :: complex_rot
     type(c_interp_t)                :: in_eps_rho
     type(c_interp_t)                :: in_eps_T
   contains
     private
     procedure, public :: set_omega_r
     procedure, public :: omega_c
     procedure, public :: lambda
     procedure, public :: l_e
     procedure, public :: l_i
     procedure, public :: eps_rho
     procedure, public :: eps_T
  end type nad_share_t

  ! Interfaces

  interface nad_share_t
     module procedure nad_share_t_
  end interface nad_share_t

  ! Access specifiers

  private

  public :: nad_share_t

  ! Procedures

contains

  function nad_share_t_ (ml, pt_i, pt_o, md_p, os_p) result (sh)

    class(model_t), pointer, intent(in) :: ml
    type(point_t), intent(in)           :: pt_i
    type(point_t), intent(in)           :: pt_o
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(nad_share_t)                   :: sh

    ! Construct the nad_share_t

    sh%ml => ml
 
    allocate(sh%rt, SOURCE=c_rot_t(md_p, os_p))

    sh%Omega_rot_i = ml%coeff(I_OMEGA_ROT, pt_i)

    sh%complex_rot = os_p%complex_rot

    ! If necessary, initialize the nuclear burning partials

    select case (os_p%deps_scheme)
    case ('FILE')
       call read_deps_(ml, md_p, os_p, pt_i, pt_o, sh%in_eps_rho, sh%in_eps_T)
    case ('MODEL')
    case default
       $ABORT(Invalid deps_scheme)
    end select

    ! Finish

    return

  end function nad_share_t_

  !****

  subroutine set_omega_r (this, omega_r)

    class(nad_share_t), intent(inout) :: this
    real(WP), intent(in)              :: omega_r

    ! Set the real frequency to be used in rotation, eps_*, etc evaluations

    this%omega_r = omega_r

    ! Finish

    return

  end subroutine set_omega_r

  !****

  function omega_c (this, Omega_rot, omega)

    class(nad_share_t), intent(in) :: this
    real(WP), intent(in)           :: Omega_rot
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: omega_c

    ! Evaluate the corotating-frame frequency

    omega_c = this%rt%omega_c(Omega_rot, omega)

    ! Finish

    return

  end function omega_c

  !****

  function lambda (this, Omega_rot, omega)

    class(nad_share_t), intent(in) :: this
    real(WP), intent(in)           :: Omega_rot
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: lambda

    ! Evaluate the angular eigenvalue

    if (this%complex_rot) then

       lambda = this%rt%lambda(Omega_rot, omega)

    else

       $ASSERT(ALLOCATED(this%omega_r),omega_r has not been set)

       lambda = this%rt%lambda(Omega_rot, CMPLX(this%omega_r, KIND=WP))

    endif

    ! Finish

    return

  end function lambda

  !****

  function l_e (this, Omega_rot, omega)

    class(nad_share_t), intent(in) :: this
    real(WP), intent(in)           :: Omega_rot
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: l_e

    ! Evaluate the effective harmonic degree

    if (this%complex_rot) then

       l_e = this%rt%l_e(Omega_rot, omega)

    else

       $ASSERT(ALLOCATED(this%omega_r),omega_r has not been set)

       l_e = this%rt%l_e(Omega_rot, CMPLX(this%omega_r, KIND=WP))

    endif

    ! Finish

    return

  end function l_e

  !****

  function l_i (this, omega)

    class(nad_share_t), intent(in) :: this
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: l_i

    ! Evaluate the effective harmonic degree at the inner boundary

    l_i = this%l_e(this%Omega_rot_i, omega)

    ! Finish

    return

  end function l_i

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
    ! preparation of Wolf et al. (2017)

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
    allocate(A_norm(n))

    allocate(A_phase(n))
    allocate(A_phase(n))

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

  function eps_rho (this)

    class(nad_share_t), intent(in) :: this
    complex(WP)                    :: eps_rho

    real(WP) :: omega_min
    real(WP) :: omega_max
    real(WP) :: omega

    ! Evaluate the eps_rho derivative

    omega_min = this%in_eps_rho%x_min()
    omega_max = this%in_eps_rho%x_min()

    omega = MIN(MAX(this%omega_r, omega_min), omega_max)

    eps_rho = this%in_eps_rho%f(omega)

    ! Finish

    return

  end function eps_rho

  !****

  function eps_T (this)

    class(nad_share_t), intent(in) :: this
    complex(WP)                    :: eps_T

    real(WP) :: omega_min
    real(WP) :: omega_max
    real(WP) :: omega

    ! Evaluate the eps_T derivative

    omega_min = this%in_eps_rho%x_min()
    omega_max = this%in_eps_rho%x_min()

    omega = MIN(MAX(this%omega_r, omega_min), omega_max)

    eps_T = this%in_eps_T%f(omega)

    ! Finish

    return

  end function eps_T

end module gyre_nad_share
