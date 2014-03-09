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
  use core_string

  use gyre_constants
  use gyre_model
  use gyre_model_evol
  use gyre_model_poly
  use gyre_mode
  use gyre_util
  use gyre_writer
  use gyre_hdf_writer
  use gyre_txt_writer

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

    character(256)               :: freq_units
    character(FILENAME_LEN)      :: summary_file
    character(256)               :: summary_file_format
    character(2048)              :: summary_item_list
    character(FILENAME_LEN)      :: mode_prefix
    character(FILENAME_LEN)      :: mode_template
    character(256)               :: mode_file_format
    character(2048)              :: mode_item_list
    character(:), allocatable    :: mode_file
    character(64)                :: infix
    class(writer_t), allocatable :: wr
    integer                      :: j

    namelist /output/ freq_units, &
         summary_file, summary_file_format, summary_item_list, &
         mode_prefix, mode_template, mode_file_format, mode_item_list

    ! Read output parameters

    freq_units = 'NONE'

    summary_file = ''
    summary_file_format = 'HDF'
    summary_item_list = 'l,n_pg,omega,freq'

    mode_prefix = ''
    mode_template = ''
    mode_file_format = 'HDF'
    mode_item_list = TRIM(summary_item_list)//',x,xi_r,xi_h'

    rewind(unit)
    read(unit, NML=output, END=900)

    ! Write output files

    if(summary_file /= '') then

       select case (summary_file_format)
       case ('HDF')
          allocate(wr, SOURCE=hdf_writer_t(summary_file))
       case ('TXT')
          allocate(wr, SOURCE=txt_writer_t(summary_file))
       case default
          $ABORT(Invalid summary_file_format)
       end select

       call write_summary_(wr, md, split_list(summary_item_list, ','), freq_units)
       call wr%final()

       deallocate(wr)

    endif

    if(mode_prefix /= '' .OR. mode_template /= '') then

       mode_loop : do j = 1,SIZE(md)

          ! Set up the mode filename

          if (mode_template /= '') then

             mode_file = mode_template

             ! Fixed-width fields

             mode_file = subst_(mode_file, '%J', j, '(I5.5)')
             mode_file = subst_(mode_file, '%L', md(j)%op%l, '(I3.3)')
             mode_file = subst_(mode_file, '%N', md(j)%n_pg, '(SP,I6.5)')

             ! Variable-width fields

             mode_file = subst_(mode_file, '%j', j, '(I0)')
             mode_file = subst_(mode_file, '%l', md(j)%op%l, '(I0)')
             mode_file = subst_(mode_file, '%n', md(j)%n_pg, '(SP,I0)')

          else

             write(infix, 100) j
100          format(I5.5)

             select case (mode_file_format)
             case ('HDF')
                mode_file = TRIM(mode_prefix)//infix//'.h5'
             case ('TXT')
                mode_file = TRIM(mode_prefix)//infix//'.txt'
             case default
                $ABORT(Invalid mode_file_format)
             end select

          endif

          select case (mode_file_format)
          case ('HDF')
             allocate(wr, SOURCE=hdf_writer_t(mode_file))
          case ('TXT')
             allocate(wr, SOURCE=txt_writer_t(mode_file))
          case default
             $ABORT(Invalid mode_file_format)
          end select

          call write_mode_(wr, md(j), split_list(mode_item_list, ','), freq_units, j)
          call wr%final()

          deallocate(wr)

       end do mode_loop

    end if

    ! Finish

    return

    ! Jump-in point for end-of-file

900 continue

    $ABORT(No &output namelist in input file)

  end subroutine write_data

!****

  subroutine write_summary_ (wr, md, items, freq_units)

    class(writer_t), intent(inout) :: wr
    type(mode_t), intent(in)       :: md(:)
    character(*), intent(in)       :: items(:)
    character(*), intent(in)       :: freq_units

    integer            :: n_md
    integer            :: i
    integer            :: j

    ! Write items

    n_md = SIZE(md)

    item_loop : do j = 1, SIZE(items)

       select case (items(j))

       case('l')
          call wr%write('l', md%op%l)
       case('n_p')
          call wr%write('n_p', md%n_p)
       case('n_g')
          call wr%write('n_g', md%n_g)
       case('n_pg')
          call wr%write('n_pg', md%n_pg)
       case('omega')
          call wr%write('omega', md%omega)
       case('freq')
          call wr%write('freq', [(md(i)%freq(freq_units), i=1,n_md)])
       case ('f_T')
          call wr%write('f_T', [(ABS(md(i)%delT_eff()/md(i)%xi_r_ref()), i=1,n_md)])
       case ('f_g')
          call wr%write('f_g', [(ABS(md(i)%delg_eff()/md(i)%xi_r_ref()), i=1,n_md)])
       case ('psi_T')
          call wr%write('psi_T', [(phase(md(i)%delT_eff()/md(i)%xi_r_ref()), i=1,n_md)])
       case ('psi_g')
          call wr%write('psi_g', [(phase(md(i)%delg_eff()/md(i)%xi_r_ref()), i=1,n_md)])
       case('beta')
          call wr%write('beta', [(md(i)%beta(), i=1,n_md)])
       case('E')
          call wr%write('E', [(md(i)%E(), i=1,n_md)])
       case('E_norm')
          call wr%write('E_norm', [(md(i)%E_norm(), i=1,n_md)])
       case('W')
          call wr%write('W', [(md(i)%W(), i=1,n_md)])
       case('omega_im')
          call wr%write('omega_im', [(md(i)%omega_im(), i=1,n_md)])
       case ('xi_r_ref')
          call wr%write('xi_r_ref', md(i)%xi_r_ref())
       case ('xi_h_ref')
          call wr%write('xi_h_ref', md(i)%xi_h_ref())
       case('freq_units')
          call wr%write('freq_units', freq_units)
       case default
          if(n_md >= 1) then
             select type (ml => md(1)%ml)
             type is (model_evol_t)
                call write_summary_evol_(wr, ml, items(j))
             type is (model_poly_t)
                call write_summary_poly_(wr, ml, items(j))
             class default
                write(ERROR_UNIT, *) 'item:', TRIM(items(j))
                $ABORT(Invalid item)
             end select
          endif
       end select

    end do item_loop

    ! Finish

    return

  contains

    subroutine write_summary_evol_ (wr, ml, item)

      class(writer_t), intent(inout) :: wr
      type(model_evol_t), intent(in) :: ml
      character(*), intent(in)       :: item

      ! Write the item

      select case (item)
      case ('M_star')
         call wr%write('M_star', ml%M_star)
      case ('R_star')
         call wr%write('R_star', ml%R_star)
      case ('L_star')
         call wr%write('L_star', ml%L_star)
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_summary_evol_

    subroutine write_summary_poly_ (wr, ml, item)

      class(writer_t), intent(inout) :: wr
      type(model_poly_t), intent(in) :: ml
      character(*), intent(in)       :: item

      ! Write the item

      select case (item)
      case ('n_poly')
         call wr%write('n_poly', ml%n_poly)
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_summary_poly_

  end subroutine write_summary_

!****

  subroutine write_mode_ (wr, md, items, freq_units, i)

    class(writer_t), intent(inout) :: wr
    type(mode_t), intent(in)       :: md
    character(*), intent(in)       :: items(:)
    character(*), intent(in)       :: freq_units
    integer, intent(in)            :: i

    integer :: j

    ! Write items

    call wr%write('i', i)

    item_loop : do j = 1, SIZE(items)

       select case (items(j))
       case ('n')
          call wr%write('n', md%n)
       case ('l')
          call wr%write('l', md%op%l)
       case ('n_p')
          call wr%write('n_p', md%n_p)
       case ('n_g')
          call wr%write('n_g', md%n_g)
       case ('n_pg')
          call wr%write('n_pg', md%n_pg)
       case ('omega')
          call wr%write('omega', md%omega)
       case ('freq')
          call wr%write('freq', md%freq(freq_units))
       case ('f_T')
          call wr%write('f_T', ABS(md%delT_eff()/md%xi_r_ref()))
       case ('f_g')
          call wr%write('f_g', ABS(md%delg_eff()/md%xi_r_ref()))
       case ('psi_T')
          call wr%write('psi_T', phase(md%delT_eff()/md%xi_r_ref()))
       case ('psi_g')
          call wr%write('psi_g', phase(md%delg_eff()/md%xi_r_ref()))
       case ('beta')
          call wr%write('beta', md%beta())
       case ('E')
          call wr%write('E', md%E())
       case ('E_norm')
          call wr%write('E_norm', md%E_norm())
       case ('W')
          call wr%write('W', md%W())
       case('omega_im')
          call wr%write('omega_im', md%omega_im())
       case ('x')
          call wr%write('x', md%x)
       case('V')
          call wr%write('V', md%ml%V(md%x))
       case('As')
          call wr%write('As', md%ml%As(md%x))
       case('U')
          call wr%write('U', md%ml%U(md%x))
       case('c_1')
          call wr%write('c_1', md%ml%c_1(md%x))
       case ('Gamma_1')
          call wr%write('Gamma_1', md%ml%Gamma_1(md%x))
       case ('nabla_ad')
          call wr%write('nabla_ad', md%ml%nabla_ad(md%x))
       case ('delta')
          call wr%write('delta', md%ml%delta(md%x))
       case ('xi_r')
          call wr%write('xi_r', md%xi_r())
       case ('xi_h')
          call wr%write('xi_h', md%xi_h())
       case ('xi_r_ref')
          call wr%write('xi_r_ref', md%xi_r_ref())
       case ('xi_h_ref')
          call wr%write('xi_h_ref', md%xi_h_ref())
       case ('phip')
          call wr%write('phip', md%phip())
       case ('dphip_dx')
          call wr%write('dphip_dx', md%dphip_dx())
       case ('delS')
          call wr%write('delS', md%delS())
       case ('delS_en')
          call wr%write('delS_en', md%delS_en())
       case ('delL')
          call wr%write('delL', md%delL())
       case ('delL_rd')
          call wr%write('delL_rd', md%delL_rd())
       case ('delp')
          call wr%write('delp', md%delp())
       case ('delrho')
          call wr%write('delrho', md%delrho())
       case ('delT')
          call wr%write('delT', md%delT())
       case ('dE_dx')
          call wr%write('dE_dx', md%dE_dx())
       case ('dW_dx')
          call wr%write('dW_dx', md%dW_dx())
       case ('Yt_1')
          call wr%write('Yt_1', md%Yt_1())
       case ('Yt_2')
          call wr%write('Yt_2', md%Yt_2())
       case ('I_0')
          call wr%write('I_0', md%I_0())
       case ('I_1')
          call wr%write('I_1', md%I_1())
       case ('prop_type')
          call wr%write('prop_type', md%prop_type())
       case ('K')
          call wr%write('K', md%K())
       case('freq_units')
          call wr%write('freq_units', freq_units)
       case default
          select type (ml => md%ml)
          type is (model_evol_t)
             call write_mode_evol_(wr, ml, items(j))
          type is (model_poly_t)
             call write_mode_poly_(wr, ml, items(j))
          class default
             write(ERROR_UNIT, *) 'item:', TRIM(items(j))
             $ABORT(Invalid item)
          end select
       end select

    end do item_loop

    ! Finish

    return

  contains

    subroutine write_mode_evol_ (wr, ml, item)

      class(writer_t), intent(inout) :: wr
      type(model_evol_t), intent(in) :: ml
      character(*), intent(in)       :: item

      ! Write the item

      select case (item)
      case ('M_star')
         call wr%write('M_star', ml%M_star)
      case ('R_star')
         call wr%write('R_star', ml%R_star)
      case ('L_star')
         call wr%write('L_star', ml%L_star)
      case ('m')
         call wr%write('m', ml%m(md%x))
      case ('p')
         call wr%write('p', ml%p(md%x))
      case ('rho')
         call wr%write('rho', ml%rho(md%x))
      case ('T')
         call wr%write('T', ml%T(md%x))
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_mode_evol_

    subroutine write_mode_poly_ (wr, ml, item)

      class(writer_t), intent(inout) :: wr
      type(model_poly_t), intent(in) :: ml
      character(*), intent(in)       :: item

      ! Write the item

      select case (item)
      case ('n_poly')
         call wr%write('n_poly', ml%n_poly)
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(item)
         $ABORT(Invalid item)
      end select

      ! Finish

      return

    end subroutine write_mode_poly_
    
  end subroutine write_mode_

!****

  function subst_ (string, pattern, i, format) result (new_string)

    character(*), intent(in)  :: string
    character(*), intent(in)  :: pattern
    integer, intent(in)       :: i
    character(*), intent(in)  :: format
    character(:), allocatable :: new_string

    character(64) :: substring

    ! Write i into the substring buffer

    write(substring, format) i

    ! Do the replacement

    new_string = replace(string, pattern, TRIM(substring), every=.TRUE.)

    ! Finish

    return

  end function subst_

end module gyre_output
