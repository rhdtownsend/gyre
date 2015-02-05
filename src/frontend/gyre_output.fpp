! Program  : gyre_output
! Purpose  : output routines
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

module gyre_output

  ! Uses

  use core_kinds
  use core_string

  use gyre_constants
  use gyre_model
  use gyre_evol_model
  use gyre_poly_model
  use gyre_out_par
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

  public :: write_summary
  public :: write_mode

contains

  subroutine write_summary (up, md)

    type(out_par_t), intent(in) :: up
    type(mode_t), intent(in)    :: md(:)

    class(writer_t), allocatable                      :: wr
    character(LEN(up%summary_item_list)), allocatable :: items(:)
    integer                                           :: n_md
    integer                                           :: i
    integer                                           :: j
    
    if (up%summary_file == '') return

    ! Write the summary file

    ! Open the file

    select case (up%summary_file_format)
    case ('HDF')
       allocate(wr, SOURCE=hdf_writer_t(up%summary_file))
    case ('TXT')
       allocate(wr, SOURCE=txt_writer_t(up%summary_file))
    case default
       $ABORT(Invalid summary_file_format)
    end select

    ! Split the item list

    items = split_list(up%summary_item_list, ',')

    ! Write the items

    n_md = SIZE(md)

    item_loop : do i = 1, SIZE(items)

       select case (items(i))

       case('l')
          call wr%write('l', md%mp%l)
       case('n_p')
          call wr%write('n_p', md%n_p)
       case('n_g')
          call wr%write('n_g', md%n_g)
       case('n_pg')
          call wr%write('n_pg', md%n_pg)
       case('omega')
          call wr%write('omega', md%omega)
       case('freq')
          call wr%write('freq', [(md(j)%freq(up%freq_units, up%freq_frame), j=1,n_md)])
       case ('f_T')
          call wr%write('f_T', [(ABS(md(j)%lag_T_eff()/md(j)%xi_r_ref()), j=1,n_md)])
       case ('f_g')
          call wr%write('f_g', [(ABS(md(j)%lag_g_eff()/md(j)%xi_r_ref()), j=1,n_md)])
       case ('psi_T')
          call wr%write('psi_T', [(phase(md(j)%lag_T_eff()/md(j)%xi_r_ref()), j=1,n_md)])
       case ('psi_g')
          call wr%write('psi_g', [(phase(md(j)%lag_g_eff()/md(j)%xi_r_ref()), j=1,n_md)])
       case('beta')
          call wr%write('beta', [(md(j)%beta(), j=1,n_md)])
       case('E')
          call wr%write('E', [(md(j)%E(), j=1,n_md)])
       case('E_norm')
          call wr%write('E_norm', [(md(j)%E_norm(), j=1,n_md)])
       case('E_ratio')
          call wr%write('E_ratio', [(md(j)%E_ratio(), j=1,n_md)])
       case('W')
          call wr%write('W', [(md(j)%W(), j=1,n_md)])
       case('omega_int')
          call wr%write('omega_int', [(md(j)%omega_int(), j=1,n_md)])
       case('eta')
          call wr%write('eta', [(md(j)%eta(), j=1,n_md)])
       case ('xi_r_ref')
          call wr%write('xi_r_ref', md(j)%xi_r_ref())
       case ('xi_h_ref')
          call wr%write('xi_h_ref', md(j)%xi_h_ref())
       case('freq_units')
          call wr%write('freq_units', up%freq_units)
       case('freq_frame')
          call wr%write('freq_frame', up%freq_frame)
       case default
          if(n_md >= 1) then
             select type (ml => md(1)%ml)
             type is (evol_model_t)
                call write_summary_evol_(items(i), ml, wr)
             type is (poly_model_t)
                call write_summary_poly_(items(i), ml, wr)
             class default
                write(ERROR_UNIT, *) 'item:', TRIM(items(i))
                $ABORT(Invalid item)
             end select
          endif
       end select

    end do item_loop

    ! Close the file

    call wr%final()

    ! Finish

    return

  contains

    subroutine write_summary_evol_ (item, ml, wr)

      character(*), intent(in)       :: item
      type(evol_model_t), intent(in) :: ml
      class(writer_t), intent(inout) :: wr

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

    subroutine write_summary_poly_ (item, ml, wr)

      character(*), intent(in)       :: item
      type(poly_model_t), intent(in) :: ml
      class(writer_t), intent(inout) :: wr

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

  end subroutine write_summary

!****

  subroutine write_mode (up, md, j)

    type(out_par_t), intent(in) :: up
    type(mode_t), intent(in)    :: md
    integer, intent(in)         :: j

    character(:), allocatable                      :: mode_file
    character(64)                                  :: infix
    class(writer_t), allocatable                   :: wr
    character(LEN(up%mode_item_list)), allocatable :: items(:)
    integer                                        :: i

    if (up%mode_template == '' .AND. up%mode_prefix == '') return

    ! Write the mode file

    ! Set up the filename

    if (up%mode_template /= '') then

       mode_file = up%mode_template

       ! Substitute fixed-width fields

       mode_file = subst_(mode_file, '%J', j, '(I5.5)')
       mode_file = subst_(mode_file, '%L', md%mp%l, '(I3.3)')
       mode_file = subst_(mode_file, '%N', md%n_pg, '(SP,I6.5)')

       ! Substitute variable-width fields

       mode_file = subst_(mode_file, '%j', j, '(I0)')
       mode_file = subst_(mode_file, '%l', md%mp%l, '(I0)')
       mode_file = subst_(mode_file, '%n', md%n_pg, '(SP,I0)')

    else

       write(infix, 100) j
100    format(I5.5)

       select case (up%mode_file_format)
       case ('HDF')
          mode_file = TRIM(up%mode_prefix)//TRIM(infix)//'.h5'
       case ('TXT')
          mode_file = TRIM(up%mode_prefix)//TRIM(infix)//'.txt'
       case default
          $ABORT(Invalid mode_file_format)
       end select

    endif

    ! Open the file

    select case (up%mode_file_format)
    case ('HDF')
       allocate(wr, SOURCE=hdf_writer_t(mode_file))
    case ('TXT')
       allocate(wr, SOURCE=txt_writer_t(mode_file))
    case default
       $ABORT(Invalid mode_file_format)
    end select

    ! Split the item list

    items = split_list(up%mode_item_list, ',')

    ! Write the items

    item_loop : do i = 1, SIZE(items)

       select case (items(i))
       case ('n')
          call wr%write('n', md%n)
       case ('l')
          call wr%write('l', md%mp%l)
       case ('n_p')
          call wr%write('n_p', md%n_p)
       case ('n_g')
          call wr%write('n_g', md%n_g)
       case ('n_pg')
          call wr%write('n_pg', md%n_pg)
       case ('omega')
          call wr%write('omega', md%omega)
       case ('freq')
          call wr%write('freq', md%freq(up%freq_units, up%freq_frame))
       case ('f_T')
          call wr%write('f_T', ABS(md%lag_T_eff()/md%xi_r_ref()))
       case ('f_g')
          call wr%write('f_g', ABS(md%lag_g_eff()/md%xi_r_ref()))
       case ('psi_T')
          call wr%write('psi_T', phase(md%lag_T_eff()/md%xi_r_ref()))
       case ('psi_g')
          call wr%write('psi_g', phase(md%lag_g_eff()/md%xi_r_ref()))
       case ('beta')
          call wr%write('beta', md%beta())
       case ('E')
          call wr%write('E', md%E())
       case ('E_norm')
          call wr%write('E_norm', md%E_norm())
       case('E_ratio')
          call wr%write('E_ratio',md%E_ratio())
       case ('W')
          call wr%write('W', md%W())
       case('omega_int')
          call wr%write('omega_int', md%omega_int())
       case('eta')
          call wr%write('eta', md%eta())
       case ('x')
          call wr%write('x', md%x)
       case('V_2')
          call wr%write('V_2', md%ml%V_2(md%x))
       case('As')
          call wr%write('As', md%ml%As(md%x))
       case('U')
          call wr%write('U', md%ml%U(md%x))
       case('c_1')
          call wr%write('c_1', md%ml%c_1(md%x))
       case ('Gamma_1')
          call wr%write('Gamma_1', md%ml%Gamma_1(md%x))
       case ('nabla')
          call wr%write('nabla', md%ml%nabla(md%x))
       case ('nabla_ad')
          call wr%write('nabla_ad', md%ml%nabla_ad(md%x))
       case ('delta')
          call wr%write('delta', md%ml%delta(md%x))
       case ('Omega_rot')
          call wr%write('Omega_rot', md%ml%Omega_rot(md%x))
       case ('xi_r')
          call wr%write('xi_r', md%xi_r())
       case ('xi_h')
          call wr%write('xi_h', md%xi_h())
       case ('xi_r_ref')
          call wr%write('xi_r_ref', md%xi_r_ref())
       case ('xi_h_ref')
          call wr%write('xi_h_ref', md%xi_h_ref())
       case ('eul_phi')
          call wr%write('eul_phi', md%eul_phi())
       case ('deul_phi')
          call wr%write('deul_phi', md%deul_phi())
       case ('lag_S')
          call wr%write('lag_S', md%lag_S())
       case ('lag_S_en')
          call wr%write('lag_S_en', md%lag_S_en())
       case ('lag_L')
          call wr%write('lag_L', md%lag_L())
       case ('lag_L_rd')
          call wr%write('lag_L_rd', md%lag_L_rd())
       case ('eul_P')
          call wr%write('eul_P', md%eul_P())
       case ('lag_P')
          call wr%write('lag_P', md%lag_P())
       case ('eul_rho')
          call wr%write('eul_rho', md%eul_rho())
       case ('lag_rho')
          call wr%write('lag_rho', md%lag_rho())
       case ('eul_T')
          call wr%write('eul_T', md%eul_T())
       case ('lag_T')
          call wr%write('lag_T', md%lag_T())
       case ('dE_dx')
          call wr%write('dE_dx', md%dE_dx())
       case ('dW_dx')
          call wr%write('dW_dx', md%dW_dx())
       case ('F_j')
          call wr%write('F_j', md%F_j())
       case ('div_F_j')
          call wr%write('div_F_j', md%div_F_j())
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
          call wr%write('freq_units', up%freq_units)
       case('freq_frame')
          call wr%write('freq_frame', up%freq_frame)
       case default
          select type (ml => md%ml)
          type is (evol_model_t)
             call write_mode_evol_(items(i), ml, wr)
          type is (poly_model_t)
             call write_mode_poly_(items(i), ml, wr)
          class default
             write(ERROR_UNIT, *) 'item:', TRIM(items(i))
             $ABORT(Invalid item)
          end select
       end select

    end do item_loop

    ! Close the file

    call wr%final()

    ! Finish

    return

  contains

    subroutine write_mode_evol_ (item, ml, wr)

      character(*), intent(in)       :: item
      type(evol_model_t), intent(in) :: ml
      class(writer_t), intent(inout) :: wr

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

    subroutine write_mode_poly_ (item, ml, wr)

      character(*), intent(in)       :: item
      type(poly_model_t), intent(in) :: ml
      class(writer_t), intent(inout) :: wr

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
    
  end subroutine write_mode

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
