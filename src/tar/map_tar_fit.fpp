! Program  : map_tar_fit
! Purpose  : maps of tar_fits in the complex plane
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

program map_tar_fit

  ! Uses

  use core_kinds
  use core_hgroup
  use core_system

  use gyre_constants
  use gyre_tar_fit

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Variables

  real(WP)                  :: nu_re_min
  real(WP)                  :: nu_re_max
  integer                   :: n_nu_re
  real(WP)                  :: nu_im_min
  real(WP)                  :: nu_im_max
  integer                   :: n_nu_im
  character(:), allocatable :: in_file
  character(:), allocatable :: out_file

  type(hgroup_t)           :: hg
  type(tar_fit_t)          :: tf
  integer                  :: i
  integer                  :: j
  real(WP), allocatable    :: nu_re(:)
  real(WP), allocatable    :: nu_im(:)
  complex(WP)              :: nu
  complex(WP), allocatable :: lambda(:,:)

  ! Read parameters

  $ASSERT(n_arg() == 8,Syntax: map_tar_fit nu_re_min nu_re_max n_nu_re nu_im_min nu_im_max n_nu_im in_file out_file)

  call get_arg(1, nu_re_min)
  call get_arg(2, nu_re_max)
  call get_arg(3, n_nu_re)
  call get_arg(4, nu_im_min)
  call get_arg(5, nu_im_max)
  call get_arg(6, n_nu_im)
  call get_arg(7, in_file)
  call get_arg(8, out_file)

  ! Read the tar_fit_t

  hg = hgroup_t(in_file, OPEN_FILE)
  call read(hg, tf)
  call hg%final()

  ! Set up axes

  nu_re = [((nu_re_min*(n_nu_re-i) + nu_re_max*(i-1))/(n_nu_re-1),i=1,n_nu_re)]
  nu_im = [((nu_im_min*(n_nu_im-j) + nu_im_max*(j-1))/(n_nu_im-1),j=1,n_nu_im)]

  ! Loop over the complex plane, calculating lambda

  allocate(lambda(n_nu_re,n_nu_im))

  re_loop : do i = 1, n_nu_re
     im_loop : do j = 1, n_nu_im
        nu = CMPLX(nu_re(i), nu_im(j), WP)
        lambda(i,j) = tf%lambda(nu)
     end do im_loop
  end do re_loop

  ! Write out data

  hg = hgroup_t(out_file, CREATE_FILE)
  call write_dset(hg, 'nu_re', nu_re)
  call write_dset(hg, 'nu_im', nu_im)
  call write_dset(hg, 'lambda', lambda)
  call hg%final()

  ! Finish

end program map_tar_fit
  
