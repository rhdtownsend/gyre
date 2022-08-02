! Module   : gyre_resp
! Purpose  : response data for a single tidal component
!
! Copyright 2013-2022 Rich Townsend & The GYRE Team
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
$include 'core_memory.inc'

module gyre_resp

  ! Uses

  use core_kinds
  use core_memory

  use gyre_context
  use gyre_freq_context
  use gyre_model
  use gyre_orbit_par
  use gyre_tidal_coeff
  use gyre_tide_par
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (wave_t) :: resp_t
     private
     type(orbit_par_t), public :: or_p
     type(tide_par_t), public  :: td_p
     integer, public           :: k
   contains
     private
     procedure, public :: eul_Phi
     procedure, public :: eul_Psi
     procedure, public :: Phi_T
     procedure, public :: deul_Phi
     procedure, public :: deul_Psi
     procedure, public :: dPhi_T
     procedure, public :: Omega_orb
     procedure, public :: R_a
     procedure, public :: c
     procedure, public :: G_1
     procedure, public :: G_2
     procedure, public :: G_3
     procedure, public :: G_4
  end type resp_t

  ! Interfaces

  interface resp_t
     module procedure resp_t_
  end interface resp_t

  interface reallocate
     module procedure reallocate_1_
  end interface reallocate

  ! Access specifiers

  private

  public :: resp_t
  public :: reallocate

  ! Procedures

contains

  function resp_t_ (wv, or_p, td_p, k) result (rs)

    type(wave_t), intent(in)      :: wv
    type(orbit_par_t), intent(in) :: or_p
    type(tide_par_t), intent(in)  :: td_p
    integer, intent(in)           :: k
    type(resp_t)                  :: rs

    ! Construct the resp_t

    rs%wave_t = wv

    rs%or_p = or_p
    rs%td_p = td_p

    rs%k = k

    ! Finish

    return

  end function resp_t_

  !****

  function eul_Phi (this, j)

    class(resp_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: eul_Phi

    ! Evaluate the Eulerian self gravitational potential
    ! perturbation, in units of G M_star / R_star

    eul_Phi = this%eul_Psi(j) - this%Phi_T(j)

    ! Finish

    return

  end function eul_Phi

  !****

  function eul_Psi (this, j)

    class(resp_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: eul_Psi

    ! Evaluate the Eulerian total gravitational potential
    ! perturbation, in units of G M_star / R_star

    eul_Psi = this%wave_t%eul_Phi(j)

    ! Finish

    return

  end function eul_Psi

  !****

  function Phi_T (this, j)

    class(resp_t), intent(in) :: this
    integer, intent(in)       :: j
    real(WP)                  :: Phi_T

    type(context_t) :: cx

    ! Evaluate the tidal potential, in units of G M_star / R_star

    cx = this%context()

    associate (ml => cx%model(),&
               x => this%x(j))
      Phi_T = tidal_Phi_T(ml, this%or_p, x, this%l, this%m, this%k)
    end associate

    ! Finish

    return

  end function Phi_T

  !****

  function deul_Phi (this, j)

    class(resp_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: deul_Phi

    ! Evaluate the Eulerian potential gradient (gravity) perturbation,
    ! in units of G M_star / R_star**2

    deul_Phi = this%deul_Psi(j) - this%dPhi_T(j)

    ! Finish

    return

  end function deul_Phi

  !****

  function deul_Psi (this, j)

    class(resp_t), intent(in) :: this
    integer, intent(in)       :: j
    complex(WP)               :: deul_Psi

    ! Evaluate the Eulerian total potential gradient (gravity)
    ! potential perturbation, in units of G M_star / R_star**2

    deul_Psi = this%wave_t%deul_Phi(j)

    ! Finish

    return

  end function deul_Psi

  !****

  function dPhi_T (this, j)

    class(resp_t), intent(in) :: this
    integer, intent(in)       :: j
    real(WP)                  :: dPhi_T

    type(context_t) :: cx

    ! Evaluate the tidal potential gradient (gravity), in units of G
    ! M_star / R_star**2

    cx = this%context()

    associate (ml => cx%model(),&
               x => this%x(j))
      dPhi_T = tidal_dPhi_T(ml, this%or_p, x, this%l, this%m, this%k)
    end associate

    ! Finish

    return

  end function dPhi_T

  !****

  function Omega_orb (this, freq_units, freq_frame)

    class(resp_t), intent(in)          :: this
    character(*), intent(in)           :: freq_units
    character(*), optional, intent(in) :: freq_frame
    real(WP)                           :: Omega_orb

    type(context_t) :: cx

    ! Evaluate the orbital frequency

    cx = this%context()

    associate (ml => cx%model())

      if (PRESENT(freq_frame)) then
         Omega_orb = (tidal_Omega_orb(ml, this%or_p) + freq_shift(freq_frame, cx, this%md_p))*freq_scale(freq_units, cx, this%md_p, this%os_p)
      else
         Omega_orb = tidal_Omega_orb(ml, this%or_p)*freq_scale(freq_units, cx, this%md_p, this%os_p)
      end if

    end associate

    ! Finish

    return

  end function Omega_orb

  !****

  function R_a (this)

    class(resp_t), intent(in) :: this
    real(WP)                  :: R_a

    type(context_t) :: cx

    ! Evaluate the ratio of the stellar radius to the semi-major axis

    cx = this%context()

    associate (ml => cx%model())
      R_a = tidal_R_a(ml, this%or_p)
    end associate

    ! Finish

    return

  end function R_a

  !****

  function c (this)

    class(resp_t), intent(in) :: this
    real(WP)                  :: c

    type(context_t) :: cx

    ! Evaluate the tidal potential coefficient

    cx = this%context()

    associate (ml => cx%model())
      c = tidal_c(ml, this%or_p, this%l, this%m, this%k)
    end associate

    ! Finish

    return

  end function c

  !****

  function G_1 (this)

    class(resp_t), intent(in) :: this
    real(WP)                  :: G_1

    type(context_t) :: cx

    ! Evaluate the secular evolution coefficient

    cx = this%context()

    associate (ml => cx%model())
      G_1 = secular_G_1(ml, this%or_p, this%l, this%m, this%k)
    end associate

    ! Finish

    return

  end function G_1

  !****

  function G_2 (this)

    class(resp_t), intent(in) :: this
    real(WP)                  :: G_2

    type(context_t) :: cx

    ! Evaluate the secular evolution coefficient

    cx = this%context()

    associate (ml => cx%model())
      G_2 = secular_G_2(ml, this%or_p, this%l, this%m, this%k)
    end associate

    ! Finish

    return

  end function G_2

  !****

  function G_3 (this)

    class(resp_t), intent(in) :: this
    real(WP)                  :: G_3

    type(context_t) :: cx

    ! Evaluate the secular evolution coefficient

    cx = this%context()

    associate (ml => cx%model())
      G_3 = secular_G_3(ml, this%or_p, this%l, this%m, this%k)
    end associate

    ! Finish

    return

  end function G_3

  !****

  function G_4 (this)

    class(resp_t), intent(in) :: this
    real(WP)                  :: G_4

    type(context_t) :: cx

    ! Evaluate the secular evolution coefficient

    cx = this%context()

    associate (ml => cx%model())
      G_4 = secular_G_4(ml, this%or_p, this%l, this%m, this%k)
    end associate

    ! Finish

    return

  end function G_4

  !****

  $REALLOCATE(type(resp_t),1)

end module gyre_resp
