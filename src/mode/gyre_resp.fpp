! Module   : gyre_resp
! Purpose  : force-oscillation response data
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

  use gyre_force_par
  use gyre_force_util
  use gyre_orbit_par
  use gyre_point
  use gyre_tide_util
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (wave_t) :: resp_t
     private
     type(force_par_t), public :: fr_p
     type(orbit_par_t), public :: or_p
   contains
     private
     procedure, public :: eul_phi
     procedure, public :: eul_psi
     procedure, public :: phi_2
     procedure, public :: F
     procedure, public :: J_dot
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

  function resp_t_ (wv, fr_p, or_p) result (rs)

     type(wave_t), intent(in)      :: wv
     type(force_par_t), intent(in) :: fr_p
     type(orbit_par_t), intent(in) :: or_p
     type(resp_t)                  :: rs

     ! Construct the resp_t

     rs%wave_t = wv

     rs%fr_p = fr_p
     rs%or_p = or_p
     
     ! Finish

     return

  end function resp_t_

  !****

  function eul_phi (this, p)

     class(resp_t), intent(in) :: this
     integer, intent(in)       :: p
     complex(WP)               :: eul_phi

     ! Evaluate the Eulerian self gravitational potential
     ! perturbation, in units of G M_star / R_star

     eul_phi = this%eul_psi(p) - this%phi_2(p)

     ! Finish

     return

  end function eul_phi

  !****

  function eul_psi (this, p)

     class(resp_t), intent(in) :: this
     integer, intent(in)       :: p
     complex(WP)               :: eul_psi

     ! Evaluate the Eulerian total gravitational potential
     ! perturbation, in units of G M_star / R_star

     eul_psi = this%wave_t%eul_phi(p)

     ! Finish

     return

  end function eul_psi

  !****

  function phi_2 (this, p)

     class(resp_t), intent(in) :: this
     integer, intent(in)       :: p
     complex(WP)               :: phi_2

     ! Evaluate the Eulerian secondary gravitational potential
     ! perturbation, in units of G M_star / R_star

     associate ( &
          md_p => this%md_p, &
          fr_p => this%fr_p, &
          or_p => this%or_p)

       phi_2 = Phi_force(md_p, fr_p, or_p)

     end associate

     ! Finish

     return

  end function phi_2

  !****

  function F (this)

     class(resp_t), intent(in) :: this
     complex(WP)               :: F

     ! Evaluate the surface response function

     F = 0.5*this%eul_phi(this%n_p)/Phi_force(this%md_p, this%fr_p, this%or_p)

     ! Finish

     return

  end function F

  !****

  function J_dot (this)

     class(resp_t), intent(in) :: this
     complex(WP)               :: J_dot

     real(WP) :: R_a
     real(WP) :: G_4
     real(WP) :: x

     ! Evaluate the external torque, in units of G*M_star**2/R_star

     associate ( &
          Omega_orb => this%or_p%Omega_orb, &
          q => this%or_p%q, &
          e => this%or_p%e, &
          l => this%md_p%l, &
          m => this%md_p%m, &
          k => this%fr_p%k)

       R_a = tidal_R_a(Omega_orb, q)

       G_4 = secular_G_4(R_a, e, l, m, k)

       x = this%x(this%n_p)

       J_dot = -2._WP*CMPLX(0._WP, 1._WP, KIND=WP)*q**2*(R_a)**(l+4)*x**(l+1)*G_4*this%F()

     end associate

     ! Finish

     return

  end function J_dot

  !****

  $REALLOCATE(type(resp_t),1)

end module gyre_resp
