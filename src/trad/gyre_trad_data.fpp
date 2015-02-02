! Program  : gyre_trad_table
! Purpose  : traditional approximation tables
!
! Copyright 2013-2015 Rich Townsend
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

module gyre_trad_data

  ! Uses

  use core_kinds
  $if ($HDF5)
  use core_hgroup
  $endif

  use gyre_trad

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: trad_table_t
     private
     type(trad_t), allocatable :: tt(:,:)
     integer                   :: l_max
   contains
  end type trad_table_t

  ! Interfaces

  interface trad_table_t
     module procedure trad_table_t_
  end interface trad_table_t

  $if ($HDF5)
  interface read
     module procedure read_
  end interface read
  interface write
     module procedure write_
  end interface write
  $endif

  $if ($MPI)
  interface bcast
     module procedure bcast_
  end interface bcast
  $endif

  ! Access specifiers

  private

  public :: trad_table_t
  $if ($HDF5)
  public :: read
  public :: write
  $endif

  ! Procedures

contains

  function trad_table_l_max_ (l_max) result (tt)

    integer, intent(in) :: l_max 
    type(trad_t)        :: tt

    ! Construct the trad_table_t

    allocate(tt(0:l_max,-l_max:l_max))

    tt%l_max = l_max

    ! Finish

    return

  end function trad_table_l_max_

!****

  $if ($HDF5)

  subroutine read_ (hg, tt)

    type(hgroup_t), intent(inout)   :: hg
    type(trad_table_t), intent(out) :: tt

    integer        :: l_max
    integer
    type(hgroup_t) :: hg_comp

    ! Read the trad_table_t

    call read_attr(hg, 'm', tr%m)
    call read_attr(hg, 'k', tr%k)

    hg_comp = hgroup_t(hg, 'cb_neg')
    call read(hg_comp, tr%cb_neg)
    call hg_comp%final()

    hg_comp = hgroup_t(hg, 'cb_pos')
    call read(hg_comp, tr%cb_pos)
    call hg_comp%final()

    hg_comp = hgroup_t(hg, 'cb_ctr')
    call read(hg_comp, tr%cb_ctr)
    call hg_comp%final()

    ! Finish

    return

  end subroutine read_

!****

  subroutine write_ (hg, tr)

    type(hgroup_t), intent(inout) :: hg
    type(trad_t), intent(in)      :: tr

    type(hgroup_t) :: hg_comp

    ! Write the trad_t

    call write_attr(hg, 'm', tr%m)
    call write_attr(hg, 'k', tr%k)

    hg_comp = hgroup_t(hg, 'cb_neg')
    call write(hg_comp, tr%cb_neg)
    call hg_comp%final()

    hg_comp = hgroup_t(hg, 'cb_pos')
    call write(hg_comp, tr%cb_pos)
    call hg_comp%final()

    hg_comp = hgroup_t(hg, 'cb_ctr')
    call write(hg_comp, tr%cb_ctr)
    call hg_comp%final()

    ! Finish

    return

  end subroutine write_

  $endif

!****

  $if ($MPI)

  subroutine bcast_ (tr, root_rank)

    class(trad_t), intent(inout) :: tr
    integer, intent(in)          :: root_rank

    ! Broadcast the trad

    call bcast(tr%m, root_rank)
    call bcast(tr%k, root_rank)

    call bcast(tr%cb_neg, root_rank)
    call bcast(tr%cb_pos, root_rank)
    call bcast(tr%cb_ctr, root_rank)

    ! Finish

    return

  end subroutine bcast_

  $endif

!****

  function lambda_r_ (this, nu) result (lambda)

    class(trad_t), intent(in) :: this
    real(WP), intent(in)      :: nu
    real(WP)                  :: lambda

    ! Evaluate the eigenvalue

    if (nu < -1._WP) then
       lambda = this%cb_neg%eval(1._WP/nu)*lambda_norm_(nu, this%m, this%k)
    elseif (nu > 1._WP) then
       lambda = this%cb_pos%eval(1._WP/nu)*lambda_norm_(nu, this%m, this%m)
    else
       lambda = this%cb_ctr%eval(nu)*lambda_norm_(nu, this%m, this%k)
    endif

    ! Finish

    return

  end function lambda_r_

!****

  function lambda_c_ (this, nu) result (lambda)

    class(trad_t), intent(in) :: this
    complex(WP), intent(in)   :: nu
    complex(WP)               :: lambda

    ! Evaluate the eigenvalue

    if (REAL(nu) < -1._WP) then
       lambda = this%cb_neg%eval(1._WP/nu)*lambda_norm_(REAL(nu), this%m, this%k)
    elseif (REAL(nu) > 1._WP) then
       lambda = this%cb_pos%eval(1._WP/nu)*lambda_norm_(REAL(nu), this%m, this%k)
    else
       lambda = this%cb_ctr%eval(nu)*lambda_norm_(REAL(nu), this%m, this%k)
    endif

    ! Finish

    return

  end function lambda_c_

!****

  function lambda_norm_ (nu, m, k) result (lambda_norm)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: lambda_norm

    ! Evaluate the eigenvalue normalization function (used to scale
    ! the eigenvalue in the Chebychev fits)

    if (nu < -1._WP) then
       lambda_norm = lambda_asymp(nu, m, k)
    elseif (nu > 1._WP) then
       lambda_norm = lambda_asymp(nu, m, k)
    else
       associate (l => ABS(m) + k)
         lambda_norm = l*(l+1)
       end associate
    end if

    ! Finish

    return

  end function lambda_norm_
    
end module gyre_trad
