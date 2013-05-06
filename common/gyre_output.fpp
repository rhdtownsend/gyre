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
  use core_hgroup

  use gyre_base_coeffs
  use gyre_evol_base_coeffs
  use gyre_poly_base_coeffs
  use gyre_eigfunc

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: write_summary
  public :: write_mode

contains

  subroutine write_summary (file, ef, items, freq_scale)

    character(LEN=*), intent(in) :: file
    type(eigfunc_t), intent(in)  :: ef(:)
    character(LEN=*), intent(in) :: items(:)
    real(WP), intent(in)         :: freq_scale

    integer        :: n_ef
    integer        :: i
    complex(WP)    :: freq(SIZE(ef))
    integer        :: n_p(SIZE(ef))
    integer        :: n_g(SIZE(ef))
    integer        :: j
    type(hgroup_t) :: hg

    ! Calculate summary data

    n_ef = SIZE(ef)

    ef_loop : do i = 1,n_ef

       freq(i) = ef(i)%omega*freq_scale

       call ef(i)%classify(n_p(i), n_g(i))

    end do ef_loop

    ! Open the file

    call hg%init(file, CREATE_FILE)

    ! Write items

    item_loop : do j = 1,SIZE(items)

       select case (items(j))

       ! Attributes

       ! Datasets

       case('l')
          call write_dset(hg, 'l', ef%op%l)
       case('n_p')
          call write_dset(hg, 'n_p', n_p)
       case('n_g')
          call write_dset(hg, 'n_g', n_g)
       case('omega')
          call write_dset(hg, 'omega', ef%omega)
       case ('freq')
          call write_dset(hg, 'freq', freq)
       case('E')
          call write_dset(hg, 'E', [(ef(i)%E(), i=1,n_ef)])
       case('E_norm')
          call write_dset(hg, 'E_norm', [(ef(i)%E_norm(), i=1,n_ef)])
       case('W')
          call write_dset(hg, 'W', [(ef(i)%W(), i=1,n_ef)])
       case('omega_im')
          call write_dset(hg, 'omega_im', [(ef(i)%omega_im(), i=1,n_ef)])
       case default
          select type (bc => ef(1)%bc)
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

  subroutine write_mode (file, ef, items, freq_scale, i)

    character(LEN=*), intent(in) :: file
    type(eigfunc_t), intent(in)  :: ef
    character(LEN=*), intent(in) :: items(:)
    real(WP), intent(in)         :: freq_scale
    integer, intent(in)          :: i

    integer        :: n_p
    integer        :: n_g
    type(hgroup_t) :: hg
    integer        :: j

    ! Calculate mode data

    call ef%classify(n_p, n_g)
    
    ! Open the file

    call hg%init(file, CREATE_FILE)

    ! Write items

    call write_attr(hg, 'i', i)

    item_loop : do j = 1,SIZE(items)

       select case (items(j))
       case ('n')
          call write_attr(hg, 'n', ef%n)
       case ('l')
          call write_attr(hg, 'l', ef%op%l)
       case ('n_p')
          call write_attr(hg, 'n_p', n_p)
       case ('n_g')
          call write_attr(hg, 'n_g', n_g)
       case ('omega')
          call write_attr(hg, 'omega', ef%omega)
       case ('freq')
          call write_attr(hg, 'freq', ef%omega*freq_scale)
       case ('E')
          call write_attr(hg, 'E', ef%E())
       case ('E_norm')
          call write_attr(hg, 'E_norm', ef%E_norm())
       case ('W')
          call write_attr(hg, 'W', ef%W())
       case('omega_im')
          call write_attr(hg, 'omega_im', ef%omega_im())
       case ('x')
          call write_dset(hg, 'x', ef%x)
       case('V')
          call write_dset(hg, 'V', ef%bc%V(ef%x))
       case('As')
          call write_dset(hg, 'As', ef%bc%As(ef%x))
       case('U')
          call write_dset(hg, 'U', ef%bc%U(ef%x))
       case('c_1')
          call write_dset(hg, 'c_1', ef%bc%c_1(ef%x))
       case ('Gamma_1')
          call write_dset(hg, 'Gamma_1', ef%bc%Gamma_1(ef%x))
       case ('nabla_ad')
          call write_dset(hg, 'nabla_ad', ef%bc%nabla_ad(ef%x))
       case ('delta')
          call write_dset(hg, 'delta', ef%bc%delta(ef%x))
       case ('xi_r')
          call write_dset(hg, 'xi_r', ef%xi_r())
       case ('xi_h')
          call write_dset(hg, 'xi_h', ef%xi_h())
       case ('phip')
          call write_dset(hg, 'phip', ef%phip())
       case ('dphip_dx')
          call write_dset(hg, 'dphip_dx', ef%dphip_dx())
       case ('delS')
          call write_dset(hg, 'delS', ef%delS())
       case ('delL')
          call write_dset(hg, 'delL', ef%delL())
       case ('delp')
          call write_dset(hg, 'delp', ef%delp())
       case ('delrho')
          call write_dset(hg, 'delrho', ef%delrho())
       case ('delT')
          call write_dset(hg, 'delT', ef%delT())
       case ('dE_dx')
          call write_dset(hg, 'dE_dx', ef%dE_dx())
       case ('dW_dx')
          call write_dset(hg, 'dW_dx', ef%dW_dx())
       case default
          select type (bc => ef%bc)
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
         call write_dset(hg, 'm', bc%m(ef%x))
      case ('p')
         call write_dset(hg, 'p', bc%p(ef%x))
      case ('rho')
         call write_dset(hg, 'rho', bc%rho(ef%x))
      case ('T')
         call write_dset(hg, 'T', bc%T(ef%x))
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
