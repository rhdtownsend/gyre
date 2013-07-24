! Program  : gyre_output
! Purpose  : output routines
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

module gyre_output

  ! Uses

  use core_kinds
  use core_constants
  use core_hgroup

  use gyre_base_coeffs
  use gyre_evol_base_coeffs
  use gyre_poly_base_coeffs
  use gyre_mode
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: write_data

contains

  subroutine write_data (unit, md)

    integer, intent(in)         :: unit
    type(mode_t), intent(in)    :: md(:)

    character(LEN=256)          :: freq_units
    character(LEN=FILENAME_LEN) :: summary_file
    character(LEN=2048)         :: summary_item_list
    character(LEN=FILENAME_LEN) :: mode_prefix
    character(LEN=2048)         :: mode_item_list
    character(LEN=FILENAME_LEN) :: mode_file
    integer                     :: j

    namelist /output/ freq_units, summary_file, summary_item_list, mode_prefix, mode_item_list

    ! Read output parameters

    freq_units = 'NONE'

    summary_file = ''
    summary_item_list = 'l,n_pg,omega,freq'

    mode_prefix = ''
    mode_item_list = TRIM(summary_item_list)//',x,xi_r,xi_h'

    rewind(unit)
    read(unit, NML=output, END=900)

    ! Write output files

    if(summary_file /= '') call write_summary(summary_file, md, split_item_list(summary_item_list), freq_units)

    if(mode_prefix /= '') then

       mode_loop : do j = 1,SIZE(md)

          write(mode_file, 100) TRIM(mode_prefix), j, '.h5'
100       format(A,I4.4,A)

          call write_mode(mode_file, md(j), split_item_list(mode_item_list), freq_units, j)

       end do mode_loop
       
    end if

    ! Finish

    return

    ! Jump-in point for end-of-file

900 continue

    $ABORT(No &output namelist in input file)

  end subroutine write_data

!****

  subroutine write_summary (file, md, items, freq_units)

    character(LEN=*), intent(in) :: file
    type(mode_t), intent(in)     :: md(:)
    character(LEN=*), intent(in) :: items(:)
    character(LEN=*), intent(in) :: freq_units

    integer        :: n_md
    integer        :: i
    integer        :: n_p(SIZE(md))
    integer        :: n_g(SIZE(md))
    integer        :: n_pg(SIZE(md))
    integer        :: j
    type(hgroup_t) :: hg

    ! Calculate summary data

    n_md = SIZE(md)

    mode_loop : do i = 1,n_md
       call md(i)%classify(n_p(i), n_g(i), n_pg(i))
    end do mode_loop

    ! Open the file

    call hg%init(file, CREATE_FILE)

    ! Write items

    item_loop : do j = 1,SIZE(items)

       select case (items(j))

       ! Attributes

       ! Datasets

       case('l')
          call write_dset(hg, 'l', md%op%l)
       case('n_p')
          call write_dset(hg, 'n_p', n_p)
       case('n_g')
          call write_dset(hg, 'n_g', n_g)
       case('n_pg')
          call write_dset(hg, 'n_pg', n_pg)
       case('omega')
          call write_dset(hg, 'omega', md%omega)
       case('freq')
          call write_dset(hg, 'freq', [(md(i)%freq(freq_units), i=1,n_md)])
       case('beta')
          call write_dset(hg, 'beta', [(md(i)%beta(), i=1,n_md)])
       case('E')
          call write_dset(hg, 'E', [(md(i)%E(), i=1,n_md)])
       case('E_norm')
          call write_dset(hg, 'E_norm', [(md(i)%E_norm(), i=1,n_md)])
       case('W')
          call write_dset(hg, 'W', [(md(i)%W(), i=1,n_md)])
       case('omega_im')
          call write_dset(hg, 'omega_im', [(md(i)%omega_im(), i=1,n_md)])
       case default
          select type (bc => md(1)%bc)
          type is (evol_base_coeffs_t)
             call write_summary_evol(hg, bc, items(j))
          type is (poly_base_coeffs_t)
             call write_summary_poly(hg, bc, items(j))
          class default
             write(ERROR_UNIT, *) 'item:', TRIM(items(j))
             $ABORT(Invalid item)
          end select
       end select

    end do item_loop

    ! Close the file

    call hg%final()

    ! Finish

    return

  contains

    subroutine write_summary_evol (hg, bc, item)

      type(hgroup_t), intent(inout)        :: hg
      type(evol_base_coeffs_t), intent(in) :: bc
      character(LEN=*), intent(in)         :: item

      ! Write the item

      select case (item)
      case ('M_star')
         call write_attr(hg, 'M_star', bc%M_star)
      case ('R_star')
         call write_attr(hg, 'R_star', bc%R_star)
      case ('L_star')
         call write_attr(hg, 'L_star', bc%L_star)
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_summary_evol

    subroutine write_summary_poly (hg, bc, item)

      type(hgroup_t), intent(inout)        :: hg
      type(poly_base_coeffs_t), intent(in) :: bc
      character(LEN=*), intent(in)         :: item

      ! Write the item

      select case (item)
      case ('n_poly')
         call write_attr(hg, 'n_poly', bc%n_poly)
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_summary_poly

  end subroutine write_summary

!****

  subroutine write_mode (file, md, items, freq_units, i)

    character(LEN=*), intent(in) :: file
    type(mode_t), intent(in)     :: md
    character(LEN=*), intent(in) :: items(:)
    character(LEN=*), intent(in) :: freq_units
    integer, intent(in)          :: i

    integer        :: n_p
    integer        :: n_g
    integer        :: n_pg
    type(hgroup_t) :: hg
    integer        :: j

    ! Calculate mode data

    call md%classify(n_p, n_g, n_pg)
    
    ! Open the file

    call hg%init(file, CREATE_FILE)

    ! Write items

    call write_attr(hg, 'i', i)

    item_loop : do j = 1,SIZE(items)

       select case (items(j))
       case ('n')
          call write_attr(hg, 'n', md%n)
       case ('l')
          call write_attr(hg, 'l', md%op%l)
       case ('n_p')
          call write_attr(hg, 'n_p', n_p)
       case ('n_g')
          call write_attr(hg, 'n_g', n_g)
       case ('n_pg')
          call write_attr(hg, 'n_pg', n_pg)
       case ('omega')
          call write_attr(hg, 'omega', md%omega)
       case ('freq')
          call write_attr(hg, 'freq', md%freq(freq_units))
       case ('beta')
          call write_attr(hg, 'beta', md%beta())
       case ('E')
          call write_attr(hg, 'E', md%E())
       case ('E_norm')
          call write_attr(hg, 'E_norm', md%E_norm())
       case ('W')
          call write_attr(hg, 'W', md%W())
       case('omega_im')
          call write_attr(hg, 'omega_im', md%omega_im())
       case ('x')
          call write_dset(hg, 'x', md%x)
       case('V')
          call write_dset(hg, 'V', md%bc%V(md%x))
       case('As')
          call write_dset(hg, 'As', md%bc%As(md%x))
       case('U')
          call write_dset(hg, 'U', md%bc%U(md%x))
       case('c_1')
          call write_dset(hg, 'c_1', md%bc%c_1(md%x))
       case ('Gamma_1')
          call write_dset(hg, 'Gamma_1', md%bc%Gamma_1(md%x))
       case ('nabla_ad')
          call write_dset(hg, 'nabla_ad', md%bc%nabla_ad(md%x))
       case ('delta')
          call write_dset(hg, 'delta', md%bc%delta(md%x))
       case ('xi_r')
          call write_dset(hg, 'xi_r', md%xi_r())
       case ('xi_h')
          call write_dset(hg, 'xi_h', md%xi_h())
       case ('Yb_1')
          call write_dset(hg, 'Yb_1', md%Yb_1())
       case ('Yb_2')
          call write_dset(hg, 'Yb_2', md%Yb_2())
       case ('phip')
          call write_dset(hg, 'phip', md%phip())
       case ('dphip_dx')
          call write_dset(hg, 'dphip_dx', md%dphip_dx())
       case ('delS')
          call write_dset(hg, 'delS', md%delS())
       case ('delS_en')
          call write_dset(hg, 'delS_en', md%delS_en())
       case ('delL')
          call write_dset(hg, 'delL', md%delL())
       case ('delL_rd')
          call write_dset(hg, 'delL_rd', md%delL_rd())
       case ('delp')
          call write_dset(hg, 'delp', md%delp())
       case ('delrho')
          call write_dset(hg, 'delrho', md%delrho())
       case ('delT')
          call write_dset(hg, 'delT', md%delT())
       case ('dE_dx')
          call write_dset(hg, 'dE_dx', md%dE_dx())
       case ('dW_dx')
          call write_dset(hg, 'dW_dx', md%dW_dx())
       case ('prop_type')
          call write_dset(hg, 'prop_type', md%prop_type())
       case ('K')
          call write_dset(hg, 'K', md%K())
       case default
          select type (bc => md%bc)
          type is (evol_base_coeffs_t)
             call write_mode_evol(hg, bc, items(j))
          type is (poly_base_coeffs_t)
             call write_mode_poly(hg, bc, items(j))
          class default
             write(ERROR_UNIT, *) 'item:', TRIM(items(j))
             $ABORT(Invalid item)
          end select
       end select

    end do item_loop

    ! Close the file

    call hg%final()

    ! Finish

    return

  contains

    subroutine write_mode_evol (hg, bc, item)

      type(hgroup_t), intent(inout)        :: hg
      type(evol_base_coeffs_t), intent(in) :: bc
      character(LEN=*), intent(in)         :: item

      ! Write the item

      select case (item)
      case ('M_star')
         call write_attr(hg, 'M_star', bc%M_star)
      case ('R_star')
         call write_attr(hg, 'R_star', bc%R_star)
      case ('L_star')
         call write_attr(hg, 'L_star', bc%L_star)
      case ('m')
         call write_dset(hg, 'm', bc%m(md%x))
      case ('p')
         call write_dset(hg, 'p', bc%p(md%x))
      case ('rho')
         call write_dset(hg, 'rho', bc%rho(md%x))
      case ('T')
         call write_dset(hg, 'T', bc%T(md%x))
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_mode_evol

    subroutine write_mode_poly (hg, bc, item)

      type(hgroup_t), intent(inout)        :: hg
      type(poly_base_coeffs_t), intent(in) :: bc
      character(LEN=*), intent(in)         :: item

      ! Write the item

      select case (item)
      case ('n_poly')
         call write_attr(hg, 'n_poly', bc%n_poly)
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_mode_poly

  end subroutine write_mode

end module gyre_output
