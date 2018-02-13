! Module   : gyre_output
! Purpose  : output routines
!
! Copyright 2013-2016 Rich Townsend
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
$include 'gyre_output.inc'
  
module gyre_output

  ! Uses

  use core_kinds
  use core_string

  use gyre_constants
  use gyre_evol_model
  use gyre_freq
  use gyre_hdf_writer
  use gyre_mode
  use gyre_model
  use gyre_out_par
  use gyre_point
  use gyre_poly_model
  use gyre_txt_writer
  use gyre_util
  use gyre_writer

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: write_summary
  public :: write_mode

contains

  subroutine write_summary (md, ot_p)

    type(mode_t), intent(in)    :: md(:)
    type(out_par_t), intent(in) :: ot_p

    class(writer_t), allocatable                        :: wr
    character(LEN(ot_p%summary_item_list)), allocatable :: items(:)
    logical                                             :: invalid_items
    integer                                             :: n_md
    real(WP)                                            :: data_r(SIZE(md))
    complex(WP)                                         :: data_c(SIZE(md))
    integer                                             :: i
    integer                                             :: i_md
    
    ! Write the summary file

    if (SIZE(md) == 0 .OR. ot_p%summary_file == '') return

    ! Open the file

    select case (ot_p%summary_file_format)
    case ('HDF')
       allocate(wr, SOURCE=hdf_writer_t(ot_p%summary_file, ot_p%label))
    case ('TXT')
       allocate(wr, SOURCE=txt_writer_t(ot_p%summary_file, ot_p%label))
    case default
       $ABORT(Invalid summary_file_format)
    end select

    ! Split the item list

    items = split_list(ot_p%summary_item_list, ',')

    ! Write the items

    invalid_items = .FALSE.

    n_md = SIZE(md)

    item_loop : do i = 1, SIZE(items)

       select case (items(i))

       case ('j')

          call wr%write('j', md%j)

       case ('l')

          call wr%write('l', md%l)

       case ('l_i')

          call wr%write('l_i', md%l_i)

       case ('m')

          call wr%write('m', md%m)

       case ('n_p')

          call wr%write('n_p', md%n_p)

       case ('n_g')

          call wr%write('n_g', md%n_g)

       case ('n_pg')

          call wr%write('n_pg', md%n_pg)

       case ('omega')

          call wr%write('omega', md%omega)

       case ('omega_int')

          do i_md = 1, n_md
             data_c(i_md) = md(i_md)%omega_int()
          end do

          call wr%write('omega_int', data_c)

       case ('freq')

          do i_md = 1, n_md
             data_c(i_md) = md(i_md)%freq(ot_p%freq_units, ot_p%freq_frame)
          end do

          call wr%write('freq', data_c)

       case ('freq_units')

          call wr%write('freq_units', ot_p%freq_units)

       case ('freq_frame')

          call wr%write('freq_frame', ot_p%freq_frame)

       case ('Delta_p')

          do i_md = 1, n_md
             data_r(i_md) = md(i_md)%cx%ml%Delta_p(md(i_md)%gr%x_i(), md(i_md)%gr%x_o())
          end do

          call wr%write('Delta_p', data_r)

       case ('Delta_g')

          do i_md = 1, n_md
             data_r(i_md) = md(i_md)%cx%ml%Delta_g(md(i_md)%gr%x_i(), md(i_md)%gr%x_o(), &
                                                   md(i_md)%l*(md(i_md)%l+1._WP))
          end do
          
          call wr%write('Delta_g', data_r)

       $OUTPUT_MODES(r,eta,eta())
       $OUTPUT_MODES(r,f_T,f_T())
       $OUTPUT_MODES(r,f_g,f_g())
       $OUTPUT_MODES(r,psi_T,psi_T())
       $OUTPUT_MODES(r,psi_g,psi_g())
       $OUTPUT_MODES(r,E,E())
       $OUTPUT_MODES(r,E_p,E_p())
       $OUTPUT_MODES(r,E_g,E_g())
       $OUTPUT_MODES(r,E_norm,E_norm())
       $OUTPUT_MODES(r,E_ratio,E_ratio())
       $OUTPUT_MODES(r,H,H())
       $OUTPUT_MODES(r,W,W())
       $OUTPUT_MODES(r,W_eps,W_eps())
       $OUTPUT_MODES(r,tau_ss,tau_ss())
       $OUTPUT_MODES(r,tau_tr,tau_tr())
       $OUTPUT_MODES(r,beta,beta())
       $OUTPUT_MODES(r,x_ref,gr%pt(md(i_md)%k_ref)%x)
       $OUTPUT_MODES(c,xi_r_ref,xi_r(md(i_md)%k_ref))
       $OUTPUT_MODES(c,xi_h_ref,xi_h(md(i_md)%k_ref))
       $OUTPUT_MODES(c,eul_phi_ref,eul_phi(md(i_md)%k_ref))
       $OUTPUT_MODES(c,deul_phi_ref,deul_phi(md(i_md)%k_ref))
       $OUTPUT_MODES(c,lag_S_ref,lag_S(md(i_md)%k_ref))
       $OUTPUT_MODES(c,lag_L_ref,lag_L(md(i_md)%k_ref))

       case default

          if (n_md >= 1) then
             select type (ml => md(1)%cx%ml)
             type is (evol_model_t)
                call write_summary_evol_(items(i), ml, wr)
             class default
                write(ERROR_UNIT, *) 'item:', TRIM(items(i))
                invalid_items = .TRUE.
             end select
          endif

       end select

    end do item_loop

    ! Write restart metadata to HDF5 files

    select case (ot_p%summary_file_format)
    case ('HDF')
       call wr%write('i', [(md(i_md)%md_p%i, i_md=1,n_md)])
    end select

    ! Close the file

    call wr%final()

    ! Check whether any invalid items were found

    if (invalid_items) then
       $ABORT(Invalid item(s) in summary_item_list)
    end if

    ! Finish

    return

  contains

    subroutine write_summary_evol_ (item, ml, wr)

      character(*), intent(in)       :: item
      type(evol_model_t), intent(in) :: ml
      class(writer_t), intent(inout) :: wr

      ! Write the item

      select case (items(i))
      case ('M_star')
         call wr%write('M_star', ml%M_star)
      case ('R_star')
         call wr%write('R_star', ml%R_star)
      case ('L_star')
         call wr%write('L_star', ml%L_star)
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(items(i))
         invalid_items = .TRUE.
      end select

      ! Finish

      return

    end subroutine write_summary_evol_

  end subroutine write_summary

  !****

  subroutine write_mode (md, ot_p)

    type(mode_t), intent(in)    :: md
    type(out_par_t), intent(in) :: ot_p

    character(:), allocatable                        :: mode_file
    class(writer_t), allocatable                     :: wr
    character(LEN(ot_p%mode_item_list)), allocatable :: items(:)
    logical                                          :: invalid_items
    integer                                          :: i
    integer                                          :: k
    real(WP)                                         :: data_r(md%n_k)
    complex(WP)                                      :: data_c(md%n_k)

    ! Write the mode file

    if (ot_p%mode_template == '') return

    if (filter_mode_(md, ot_p)) return

    ! Set up the filename

    mode_file = ot_p%mode_template

    mode_file = subst_(mode_file, '%J', md%j, '(I5.5)')
    mode_file = subst_(mode_file, '%L', md%l, '(I3.3)')
    mode_file = subst_(mode_file, '%M', md%m, '(SP,I3.2)')
    mode_file = subst_(mode_file, '%N', md%n_pg, '(SP,I6.5)')

    mode_file = subst_(mode_file, '%j', md%j, '(I0)')
    mode_file = subst_(mode_file, '%l', md%l, '(I0)')
    mode_file = subst_(mode_file, '%m', md%m, '(SP,I0)')
    mode_file = subst_(mode_file, '%n', md%n_pg, '(SP,I0)')

    ! Open the file

    select case (ot_p%mode_file_format)
    case ('HDF')
       allocate(wr, SOURCE=hdf_writer_t(mode_file, ot_p%label))
    case ('TXT')
       allocate(wr, SOURCE=txt_writer_t(mode_file, ot_p%label))
    case default
       $ABORT(Invalid mode_file_format)
    end select
    
    ! Split the item list

    items = split_list(ot_p%mode_item_list, ',')

    ! Write the items

    invalid_items = .FALSE.

    associate (pt => md%gr%pt)

      item_loop : do i = 1, SIZE(items)

         select case (items(i))

         case ('n')

            call wr%write('n', md%n_k)

         case ('j')

            call wr%write('j', md%j)

         case ('l')

            call wr%write('l', md%l)

         case ('l_i')

            call wr%write('l_i', md%l_i)

         case ('m')

            call wr%write('m', md%m)

         case ('lambda')

            call wr%write('lambda', [(md%lambda(k), k=1,md%n_k)])

         case ('n_p')

            call wr%write('n_p', md%n_p)

         case ('n_g')

            call wr%write('n_g', md%n_g)

         case ('n_pg')

            call wr%write('n_pg', md%n_pg)

         case ('omega')

            call wr%write('omega', md%omega)

         case ('omega_int')

            call wr%write('omega_int', md%omega_int())

         case ('freq')

            call wr%write('freq', md%freq(ot_p%freq_units, ot_p%freq_frame))

         case ('freq_units')

            call wr%write('freq_units', ot_p%freq_units)

         case ('freq_frame')

            call wr%write('freq_frame', ot_p%freq_frame)

         case ('Delta_p')

            call wr%write('Delta_p', md%cx%ml%Delta_p(md%gr%x_i(), md%gr%x_o()))

         case ('Delta_g')

            call wr%write('Delta_g', md%cx%ml%Delta_g(md%gr%x_i(), md%gr%x_o(), md%l*(md%l+1._WP)))

         case ('eta')

            call wr%write('eta', md%eta())

         case ('f_T')

            call wr%write('f_T', md%f_T())

         case ('f_g')

            call wr%write('f_g', md%f_g())

         case ('psi_T')

            call wr%write('psi_T', md%psi_T())

         case ('psi_g')

            call wr%write('psi_g', md%psi_g())

         case ('E')

            call wr%write('E', md%E())

         case ('E_p')

            call wr%write('E_p', md%E_p())

         case ('E_g')

            call wr%write('E_g', md%E_g())

         case ('E_norm')

            call wr%write('E_norm', md%E_norm())

         case ('E_ratio')

            call wr%write('E_ratio', md%E_ratio())

         case ('H')

            call wr%write('H', md%H())

         case ('W')

            call wr%write('W', md%W())

         case ('W_eps')

            call wr%write('W_eps', md%W_eps())

         case ('tau_ss')

            call wr%write('tau_ss', md%tau_ss())

         case ('tau_tr')

            call wr%write('tau_tr', md%tau_tr())

         case ('beta')

            call wr%write('beta', md%beta())

         case ('x')

            $if ($GFORTRAN_PR_49636)
            call wr%write('x', md%gr%pt%x)
            $else
            call wr%write('x', pt%x)
            $endif

         case ('x_ref')

            call wr%write('x_ref', pt(md%k_ref)%x)

         $OUTPUT_POINTS(r,V_2,cx%ml%coeff(I_V_2, pt(k)))
         $OUTPUT_POINTS(r,As,cx%ml%coeff(I_AS, pt(k)))
         $OUTPUT_POINTS(r,U,cx%ml%coeff(I_U, pt(k)))
         $OUTPUT_POINTS(r,c_1,cx%ml%coeff(I_C_1, pt(k)))
         $OUTPUT_POINTS(r,Gamma_1,cx%ml%coeff(I_GAMMA_1, pt(k)))
         $OUTPUT_POINTS(r,nabla,cx%ml%coeff(I_NABLA, pt(k)))
         $OUTPUT_POINTS(r,nabla_ad,cx%ml%coeff(I_NABLA_AD, pt(k)))
         $OUTPUT_POINTS(r,dnabla_ad,cx%ml%dcoeff(I_NABLA_AD, pt(k)))
         $OUTPUT_POINTS(r,delta,cx%ml%coeff(I_DELTA, pt(k)))
         $OUTPUT_POINTS(r,c_lum,cx%ml%coeff(I_C_LUM, pt(k)))
         $OUTPUT_POINTS(r,c_rad,cx%ml%coeff(I_C_RAD, pt(k)))
         $OUTPUT_POINTS(r,c_thn,cx%ml%coeff(I_C_THN, pt(k)))
         $OUTPUT_POINTS(r,c_thk,cx%ml%coeff(I_C_THK, pt(k)))
         $OUTPUT_POINTS(r,c_eps,cx%ml%coeff(I_C_EPS, pt(k)))
         $OUTPUT_POINTS(r,eps_rho,cx%ml%coeff(I_EPS_RHO, pt(k)))
         $OUTPUT_POINTS(r,eps_T,cx%ml%coeff(I_EPS_T, pt(k)))
         $OUTPUT_POINTS(r,kap_rho,cx%ml%coeff(I_KAP_RHO, pt(k)))
         $OUTPUT_POINTS(r,kap_T,cx%ml%coeff(I_KAP_T, pt(k)))
         $OUTPUT_POINTS(r,Omega_rot,cx%ml%coeff(I_OMEGA_ROT, pt(k)))
         $OUTPUT_POINTS(c,y_1,y_i(1, k))
         $OUTPUT_POINTS(c,y_2,y_i(2, k))
         $OUTPUT_POINTS(c,y_3,y_i(3, k))
         $OUTPUT_POINTS(c,y_4,y_i(4, k))
         $OUTPUT_POINTS(c,y_5,y_i(5, k))
         $OUTPUT_POINTS(c,y_6,y_i(6, k))
         $OUTPUT_POINTS(c,xi_r,xi_r(k))
         $OUTPUT_POINTS(c,xi_h,xi_h(k))
         $OUTPUT_POINTS(c,eul_phi,eul_phi(k))
         $OUTPUT_POINTS(c,deul_phi,deul_phi(k))
         $OUTPUT_POINTS(c,eul_P,eul_P(k))
         $OUTPUT_POINTS(c,eul_rho,eul_rho(k))
         $OUTPUT_POINTS(c,eul_T,eul_T(k))
         $OUTPUT_POINTS(c,lag_P,lag_P(k))
         $OUTPUT_POINTS(c,lag_rho,lag_rho(k))
         $OUTPUT_POINTS(c,lag_T,lag_T(k))
         $OUTPUT_POINTS(c,lag_S,lag_S(k))
         $OUTPUT_POINTS(c,lag_L,lag_L(k))
         $OUTPUT_POINTS(r,dE_dx,dE_dx(k))
         $OUTPUT_POINTS(r,dW_dx,dW_dx(k))
         $OUTPUT_POINTS(r,dW_eps_dx,dW_eps_dx(k))
         $OUTPUT_POINTS(r,dbeta_dx,dbeta_dx(k))
         $OUTPUT_POINTS(r,dtau_dx_ss,dtau_dx_ss(k))
         $OUTPUT_POINTS(r,dtau_dx_tr,dtau_dx_tr(k))
         $OUTPUT_POINTS(c,Yt_1,Yt_1(k))
         $OUTPUT_POINTS(c,Yt_2,Yt_2(k))
         $OUTPUT_POINTS(c,I_0,I_0(k))
         $OUTPUT_POINTS(c,I_1,I_1(k))
         $OUTPUT_POINTS(r,alpha_0,alpha_0(k))
         $OUTPUT_POINTS(r,alpha_1,alpha_1(k))
         $OUTPUT_POINTS(c,prop_type,prop_type(k))

         $OUTPUT_REF(xi_r_ref,xi_r)
         $OUTPUT_REF(xi_h_ref,xi_h)
         $OUTPUT_REF(eul_phi_ref,eul_phi)
         $OUTPUT_REF(deul_phi_ref,deul_phi)
         $OUTPUT_REF(lag_S_ref,lag_S)
         $OUTPUT_REF(lag_L_ref,lag_L)

         case default

            select type (ml => md%cx%ml)
            type is (evol_model_t)
               call write_mode_evol_(ml)
            class default
               write(ERROR_UNIT, *) 'item:', TRIM(items(i))
               invalid_items = .TRUE.
            end select

         end select

      end do item_loop

    end associate

    ! Close the file

    call wr%final()

    ! Check whether any invalid items were found

    if (invalid_items) then
       $ABORT(Invalid item(s) in mode_item_list)
    end if

    ! Finish

    return

  contains

    subroutine write_mode_evol_ (ml)

      type(evol_model_t), intent(in) :: ml

      integer :: k
      
      ! Write the item

      associate (pt => md%gr%pt)
      
        select case (items(i))
        case ('M_star')
           call wr%write('M_star', ml%M_star)
        case ('R_star')
           call wr%write('R_star', ml%R_star)
        case ('L_star')
           call wr%write('L_star', ml%L_star)
        case ('M_r')
           call wr%write('M_r', [(ml%M_r(pt(k)), k=1,SIZE(pt))])
        case ('P')
           call wr%write('P', [(ml%P(pt(k)), k=1,SIZE(pt))])
        case ('rho')
           call wr%write('rho', [(ml%rho(pt(k)), k=1,SIZE(pt))])
        case ('T')
           call wr%write('T', [(ml%T(pt(k)), k=1,SIZE(pt))])
        case default
           write(ERROR_UNIT, *) 'item:', TRIM(items(i))
           invalid_items = .TRUE.
        end select

      end associate

      ! Finish

      return

    end subroutine write_mode_evol_

  end subroutine write_mode

  !****

  function filter_mode_ (md, ot_p) result (filter_mode)

    type(mode_t), intent(in)    :: md
    type(out_par_t), intent(in) :: ot_p
    logical                     :: filter_mode

    character(LEN(ot_p%mode_filter_list)), allocatable :: filters(:)
    integer                                            :: i

    ! Decide whether to filter the mode

    filters = split_list(ot_p%mode_filter_list, ',')

    filter_mode = .FALSE.

    item_loop : do i = 1, SIZE(filters)

       select case (filters(i))
       case ('stable')
          filter_mode = filter_mode .OR. AIMAG(md%omega) <= 0._WP
       case ('unstable')
          filter_mode = filter_mode .OR. AIMAG(md%omega) > 0._WP
       case default
          $ABORT(Unrecognized filter in mode_filter_list)
       end select

    end do item_loop

    ! Finish

    return

  end function filter_mode_

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
