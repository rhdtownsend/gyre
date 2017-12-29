! Incfile  : gyre_tar_rot
! Purpose  : rotational effects (traditional approximation of rotation, template)
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

module gyre_${T}_tar_rot

  ! Uses

  use core_kinds
  use core_hgroup

  use gyre_constants
  use gyre_freq
  use gyre_mode_par
  use gyre_tar_fit
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (${T}_rot_t) :: ${T}_tar_rot_t
     private
     type(tar_fit_t) :: tf
     integer         :: l
     integer         :: m
     logical         :: rossby
   contains
     private
     procedure, public :: l_e
     procedure, public :: lambda
  end type ${T}_tar_rot_t

  ! Interfaces

  interface ${T}_tar_rot_t
     module procedure ${T}_tar_rot_t_
  end interface ${T}_tar_rot_t

  ! Access specifiers

  private

  public :: ${T}_tar_rot_t

  ! Procedures

contains

  function ${T}_tar_rot_t_ (md_p) result (rt)

    type(mode_par_t), intent(in) :: md_p
    type(${T}_tar_rot_t)         :: rt

    integer                 :: k
    character(FILENAME_LEN) :: filename
    type(hgroup_t)          :: hg

    ! Construct the tar_rot_t

    if (md_p%rossby) then
       k = -(md_p%l - ABS(md_p%m) + 1)
    else
       k = md_p%l - ABS(md_p%m)
    endif

    !$OMP CRITICAL

    write(filename, 100) md_p%m, k
100 format(SP,'tar_fit.m',I0,'.k',I0,'.h5')

    hg = hgroup_t(TRIM(GYRE_DIR)//'/data/tar/'//TRIM(filename), OPEN_FILE)
    call read(hg, rt%tf)
    call hg%final()

    !$OMP END CRITICAL

    rt%l = md_p%l
    rt%m = md_p%m

    rt%rossby = md_p%rossby

    ! Finish

    return

  end function ${T}_tar_rot_t_

  !****

  function l_e (this, Omega_rot, omega)

    class(${T}_tar_rot_t), intent(in) :: this
    real(WP), intent(in)              :: Omega_rot
    $TYPE(WP), intent(in)             :: omega
    $TYPE(WP)                         :: l_e

    ! Evaluate the effective harmonic degree

    l_e = 0.5_WP*(-1._WP + SQRT(1._WP + 4._WP*this%lambda(Omega_rot, omega)))

    ! Finish

    return

  end function l_e

  !****

  function lambda (this, Omega_rot, omega)

    class(${T}_tar_rot_t), intent(in) :: this
    real(WP), intent(in)              :: Omega_rot
    $TYPE(WP), intent(in)             :: omega
    $TYPE(WP)                         :: lambda

    $TYPE(WP) :: nu

    ! Evaluate the angular eigenvalue

    nu = 2._WP*Omega_rot/omega_corot(omega, Omega_rot, this%m)

    lambda = this%tf%lambda(nu)

    ! Finish

    return

  end function lambda

end module gyre_${T}_tar_rot

$endsub
