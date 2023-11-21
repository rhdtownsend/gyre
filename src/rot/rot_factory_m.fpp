! Module  : rot_factory_m
! Purpose : factory procedures for r_rot_t and c_rot_t types
!
! Copyright 2013-2020 Rich Townsend & The GYRE Team
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

module rot_factory_m

  ! Uses

  use kinds_m

  use mode_par_m
  use rot_m
  use rot_par_m

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface r_rot_t
     module procedure r_rot_t_
  end interface r_rot_t

  interface c_rot_t
     module procedure c_rot_t_
  end interface c_rot_t

  ! Access specifiers

  private

  public :: r_rot_t
  public :: c_rot_t

  ! Procedures

contains

  $define $ROT_T $sub

  $local $T $1

  function ${T}_rot_t_ (md_p, rt_p) result (rt)

    use null_rot_m
    $if ($HDF5)
    use tar_rot_m
    $endif

    type(mode_par_t), intent(in)   :: md_p
    type(rot_par_t), intent(in)    :: rt_p
    class(${T}_rot_t), allocatable :: rt
    
    ! Create a ${T}_rot_t

    select case (rt_p%coriolis_method)
    case ('NULL')
       allocate(rt, SOURCE=${T}_null_rot_t(md_p))
    case ('TAR')
       $if ($HDF5)
       allocate(rt, SOURCE=${T}_tar_rot_t(md_p, rt_p))
       $else
       $ABORT(TAR rotation method requires HDF support be enabled)
       $endif
    case default
       $ABORT(Invalid coriolis_method)
    end select

    ! Finish

    return

  end function ${T}_rot_t_

  $endsub

  $ROT_T(r)
  $ROT_T(c)

end module rot_factory_m
