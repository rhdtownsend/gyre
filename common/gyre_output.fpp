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

  use gyre_mech_coeffs
  use gyre_evol_mech_coeffs
  use gyre_poly_mech_coeffs
  use gyre_eigfunc

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: write_summary
  public :: write_mode

contains

  subroutine write_summary (file, ef, mc, items, freq_units)

    character(LEN=*), intent(in)     :: file
    type(eigfunc_t), intent(in)      :: ef(:)
    class(mech_coeffs_t), intent(in) :: mc
    character(LEN=*), intent(in)     :: items(:)
    character(LEN=*), intent(in)     :: freq_units

    integer        :: i
    complex(WP)    :: freq(SIZE(ef))
    integer        :: n_p(SIZE(ef))
    integer        :: n_g(SIZE(ef))
    real(WP)       :: E(SIZE(ef))
    real(WP)       :: K(SIZE(ef))
    integer        :: j
    type(hgroup_t) :: hg

    ! Calculate summary data

    ef_loop : do i = 1,SIZE(ef)

       freq(i) = mc%conv_freq(ef(i)%omega, 'NONE', freq_units)

       call ef(i)%classify(n_p(i), n_g(i))

       E(i) = ef(i)%E(mc)
       K(i) = ef(i)%K(mc)

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
          call write_attr(hg, 'freq_units', freq_units)
       case('E')
          call write_dset(hg, 'E', E)
       case('K')
          call write_dset(hg, 'K', K)
       case default
          select type (mc)
          type is (evol_mech_coeffs_t)
             call write_summary_evol(hg, mc, items(j))
          type is (poly_mech_coeffs_t)
             call write_summary_poly(hg, mc, items(j))
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

    subroutine write_summary_evol (hg, mc, item)

      type(hgroup_t), intent(inout)        :: hg
      type(evol_mech_coeffs_t), intent(in) :: mc
      character(LEN=*), intent(in)         :: item

      ! Write the item

      select case (item)
      case ('M_star')
         call write_attr(hg, 'M_star', mc%M_star)
      case ('R_star')
         call write_attr(hg, 'R_star', mc%R_star)
      case ('L_star')
         call write_attr(hg, 'L_star', mc%L_star)
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_summary_evol

    subroutine write_summary_poly (hg, mc, item)

      type(hgroup_t), intent(inout)        :: hg
      type(poly_mech_coeffs_t), intent(in) :: mc
      character(LEN=*), intent(in)         :: item

      ! Write the item

      select case (item)
      case ('n_poly')
         call write_attr(hg, 'n_poly', mc%n_poly)
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_summary_poly

  end subroutine write_summary

!****

  subroutine write_mode (file, ef, mc, items, freq_units)

    character(LEN=*), intent(in)     :: file
    type(eigfunc_t), intent(in)      :: ef
    class(mech_coeffs_t), intent(in) :: mc
    character(LEN=*), intent(in)     :: items(:)
    character(LEN=*), intent(in)     :: freq_units

    integer        :: n_p
    integer        :: n_g
    real(WP)       :: E
    real(WP)       :: K
    type(hgroup_t) :: hg
    integer        :: j

    ! Calculate mode data

    call ef%classify(n_p, n_g)
    
    E = ef%E(mc)
    K = ef%K(mc)

    ! Open the file

    call hg%init(file, CREATE_FILE)

    ! Write items

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
          call write_attr(hg, 'freq', mc%conv_freq(ef%omega, 'NONE', freq_units))
          call write_attr(hg, 'freq_units', freq_units)
       case ('E')
          call write_attr(hg, 'E', E)
       case ('K')
          call write_attr(hg, 'K', K)
       case ('x')
          call write_dset(hg, 'x', ef%x)
       case('V')
          call write_dset(hg, 'V', mc%V(ef%x))
       case('As')
          call write_dset(hg, 'As', mc%As(ef%x))
       case('U')
          call write_dset(hg, 'U', mc%U(ef%x))
       case('c_1')
          call write_dset(hg, 'c_1', mc%c_1(ef%x))
       case ('Gamma_1')
          call write_dset(hg, 'Gamma_1', mc%Gamma_1(ef%x))
       case ('xi_r')
          call write_dset(hg, 'xi_r', ef%xi_r)
       case ('xi_h')
          call write_dset(hg, 'xi_h', ef%xi_h)
       case ('phip')
          call write_dset(hg, 'phip', ef%phip)
       case ('dphip_dx')
          call write_dset(hg, 'dphip_dx', ef%dphip_dx)
       case ('delS')
          call write_dset(hg, 'delS', ef%delS)
       case ('delL')
          call write_dset(hg, 'delL', ef%delL)
       case ('dK_dx')
          call write_dset(hg, 'dK_dx', ef%dK_dx(mc))
       case default
          select type (mc)
          type is (evol_mech_coeffs_t)
             call write_mode_evol(hg, mc, items(j))
          type is (poly_mech_coeffs_t)
             call write_mode_poly(hg, mc, items(j))
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

    subroutine write_mode_evol (hg, mc, item)

      type(hgroup_t), intent(inout)        :: hg
      type(evol_mech_coeffs_t), intent(in) :: mc
      character(LEN=*), intent(in)         :: item

      ! Write the item

      select case (item)
      case ('M_star')
         call write_attr(hg, 'M_star', mc%M_star)
      case ('R_star')
         call write_attr(hg, 'R_star', mc%R_star)
      case ('L_star')
         call write_attr(hg, 'L_star', mc%L_star)
      case ('m')
         call write_dset(hg, 'm', mc%m(ef%x))
      case ('p')
         call write_dset(hg, 'p', mc%p(ef%x))
      case ('rho')
         call write_dset(hg, 'rho', mc%rho(ef%x))
      case ('T')
         call write_dset(hg, 'T', mc%T(ef%x))
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_mode_evol

    subroutine write_mode_poly (hg, mc, item)

      type(hgroup_t), intent(inout)        :: hg
      type(poly_mech_coeffs_t), intent(in) :: mc
      character(LEN=*), intent(in)         :: item

      ! Write the item

      select case (item)
      case ('n_poly')
         call write_attr(hg, 'n_poly', mc%n_poly)
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_mode_poly

  end subroutine write_mode

end module gyre_output
