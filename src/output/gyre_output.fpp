! Module   : gyre_output
! Purpose  : output routines
!
! Copyright 2013-2019 Rich Townsend
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
  use gyre_context
  use gyre_model
  use gyre_evol_model
  use gyre_freq
  use gyre_grid
  use gyre_hdf_writer
  use gyre_mode
  use gyre_model
  use gyre_out_par
  use gyre_point
  use gyre_poly_model
  use gyre_txt_writer
  use gyre_util
  use gyre_wave
  use gyre_writer

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: write_summary
  public :: write_details

contains

  subroutine write_summary (wv, ot_p)

    class(wave_t), intent(in)   :: wv(:)
    type(out_par_t), intent(in) :: ot_p

    class(writer_t), allocatable                        :: wr
    character(LEN(ot_p%summary_item_list)), allocatable :: items(:)
    type(context_t), pointer                            :: cx
    class(model_t), pointer                             :: ml
    integer                                             :: i
    integer                                             :: c
    
    ! Write the summary file

    if (SIZE(wv) == 0 .OR. ot_p%summary_file == '') return

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

    ! Get the contect and model (from the first wave)

    cx => wv(1)%context()
    ml => cx%model()

    ! Write the items

    item_loop : do i = 1, SIZE(items)

       ! Cache the current write count

       c = wr%count()

       ! Write the item

       call write_summary_wave_(items(i), wv, wr)

       select type (wv)
       class is (mode_t)
          call write_summary_mode_(items(i), wv, wr)
       end select

       select type (ml)
       class is (evol_model_t)
          call write_summary_evol_model_(items(i), ml, wr)
       end select

       ! Check whether an item was written

       if (wr%count() /= c+1) then
          write(ERROR_UNIT, *) 'Ignoring invalid summary item:', TRIM(items(i))
       endif

    end do item_loop

    ! Close the file

    call wr%final()

    ! Finish

    return

  contains

    subroutine write_summary_wave_ (item, wv, wr)

      character(*), intent(in)       :: item
      class(wave_t), intent(in)      :: wv(:)
      class(writer_t), intent(inout) :: wr

      integer                  :: n_wv
      integer                  :: i_wv
      real(WP)                 :: data_r(SIZE(wv))
      complex(WP)              :: data_c(SIZE(wv))
      type(context_t), pointer :: cx
      class(model_t), pointer  :: ml
      type(grid_t)             :: gr

      ! Write the item from wave_t data

      n_wv = SIZE(wv)

      select case (item)

      case ('j')

         call wr%write('j', wv%j)

      case ('l')

         call wr%write('l', wv%l)

      case ('l_i')

         call wr%write('l_i', wv%l_i)

      case ('m')

         call wr%write('m', wv%m)

      case ('omega')

         call wr%write('omega', wv%omega)

      case ('freq')

         do i_wv = 1, n_wv
            data_c(i_wv) = wv(i_wv)%freq(ot_p%freq_units, ot_p%freq_frame)
         end do

         call wr%write('freq', data_c)

      case ('dfreq_rot')

         do i_wv = 1, n_wv
            data_r(i_wv) = wv(i_wv)%dfreq_rot(ot_p%freq_units)
         end do

         call wr%write('dfreq_rot', data_r)

      case ('freq_units')

         call wr%write('freq_units', ot_p%freq_units)

      case ('freq_frame')

         call wr%write('freq_frame', ot_p%freq_frame)

      case ('Delta_p')

         do i_wv = 1, n_wv
            cx => wv(i_wv)%context()
            ml => cx%model()
            gr = wv(i_wv)%grid()
            data_r(i_wv) = ml%Delta_p(gr%x_i(), gr%x_o())
         end do

         call wr%write('Delta_p', data_r)

      case ('Delta_g')

         do i_wv = 1, n_wv
            cx => wv(i_wv)%context()
            ml => cx%model()
            gr = wv(i_wv)%grid()
            data_r(i_wv) = ml%Delta_g(gr%x_i(), gr%x_o(), wv(i_wv)%l*(wv(i_wv)%l+1._WP))
         end do

         call wr%write('Delta_g', data_r)

      case ('x_ref')

         do i_wv = 1, n_wv
            gr = wv(i_wv)%grid()
            data_r(i_wv) = gr%pt(wv(i_wv)%k_ref)%x
         end do

         call wr%write('x_ref', data_r)

      $OUTPUT_WAVES(c,omega_int,omega_int())
      $OUTPUT_WAVES(r,domega_rot,domega_rot())
      $OUTPUT_WAVES(r,eta,eta())
      $OUTPUT_WAVES(r,f_T,f_T())
      $OUTPUT_WAVES(r,f_g,f_g())
      $OUTPUT_WAVES(r,psi_T,psi_T())
      $OUTPUT_WAVES(r,psi_g,psi_g())
      $OUTPUT_WAVES(r,E,E())
      $OUTPUT_WAVES(r,E_p,E_p())
      $OUTPUT_WAVES(r,E_g,E_g())
      $OUTPUT_WAVES(r,E_norm,E_norm())
      $OUTPUT_WAVES(r,E_ratio,E_ratio())
      $OUTPUT_WAVES(r,H,H())
      $OUTPUT_WAVES(r,W,W())
      $OUTPUT_WAVES(r,W_eps,W_eps())
      $OUTPUT_WAVES(r,tau_ss,tau_ss())
      $OUTPUT_WAVES(r,tau_tr,tau_tr())
      $OUTPUT_WAVES(r,beta,beta())
      $OUTPUT_WAVES(c,xi_r_ref,xi_r(wv(i_wv)%k_ref))
      $OUTPUT_WAVES(c,xi_h_ref,xi_h(wv(i_wv)%k_ref))
      $OUTPUT_WAVES(c,eul_phi_ref,eul_phi(wv(i_wv)%k_ref))
      $OUTPUT_WAVES(c,deul_phi_ref,deul_phi(wv(i_wv)%k_ref))
      $OUTPUT_WAVES(c,lag_S_ref,lag_S(wv(i_wv)%k_ref))
      $OUTPUT_WAVES(c,lag_L_ref,lag_L(wv(i_wv)%k_ref))

      end select

      ! Finish

      return

    end subroutine write_summary_wave_

    !****         

    subroutine write_summary_mode_ (item, md, wr)

      character(*), intent(in)       :: item
      class(mode_t), intent(in)      :: md(:)
      class(writer_t), intent(inout) :: wr

      ! Write the item from mode_t data

      select case (item)

      case ('n_p')

         call wr%write('n_p', md%n_p)

      case ('n_g')

         call wr%write('n_g', md%n_g)

      case ('n_pg')

         call wr%write('n_pg', md%n_pg)

      end select

      ! Finish

      return

    end subroutine write_summary_mode_

    !****

    subroutine write_summary_evol_model_ (item, ml, wr)

      character(*), intent(in)        :: item
      class(evol_model_t), intent(in) :: ml
      class(writer_t), intent(inout)  :: wr

      ! Write the item from evol_model_t data

      select case (item)

      case ('M_star')

         call wr%write('M_star', ml%M_star)

      case ('R_star')

         call wr%write('R_star', ml%R_star)

      case ('L_star')

         call wr%write('L_star', ml%L_star)

      end select

      ! Finish

      return

    end subroutine write_summary_evol_model_

  end subroutine write_summary

  !****
                   
  subroutine write_details (wv, ot_p)

    class(wave_t), intent(in)   :: wv
    type(out_par_t), intent(in) :: ot_p

    character(:), allocatable                           :: details_file
    class(writer_t), allocatable                        :: wr
    character(LEN(ot_p%details_item_list)), allocatable :: items(:)
    type(context_t), pointer                            :: cx
    class(model_t), pointer                             :: ml
    integer                                             :: i
    integer                                             :: c

    ! Write the detail file

    if (ot_p%details_template == '') return

    if (filter_wave_(wv, ot_p)) return

    ! Set up the filename

    details_file = ot_p%details_template

    details_file = subst_(details_file, '%J', wv%j, '(I5.5)')
    details_file = subst_(details_file, '%L', wv%l, '(I3.3)')
    details_file = subst_(details_file, '%M', wv%m, '(SP,I3.2)')
    details_file = subst_(details_file, '%j', wv%j, '(I0)')
    details_file = subst_(details_file, '%l', wv%l, '(I0)')
    details_file = subst_(details_file, '%m', wv%m, '(SP,I0)')

    select type (wv)
    class is (mode_t)
       details_file = subst_(details_file, '%N', wv%n_pg, '(SP,I6.5)')
       details_file = subst_(details_file, '%n', wv%n_pg, '(SP,I0)')
    end select

    ! Open the file

    select case (ot_p%details_file_format)
    case ('HDF')
       allocate(wr, SOURCE=hdf_writer_t(details_file, ot_p%label))
    case ('TXT')
       allocate(wr, SOURCE=txt_writer_t(details_file, ot_p%label))
    case default
       $ABORT(Invalid details_file_format)
    end select
    
    ! Split the item list

    items = split_list(ot_p%details_item_list, ',')

    ! Get the context, model and grid

    cx => wv%context()
    ml => cx%model()

    ! Write the items

    item_loop : do i = 1, SIZE(items)

       ! Cache the current write count

       c = wr%count()

       ! Write the item

       call write_details_wave_(items(i), wv, wr)

       select type (wv)
       class is (mode_t)
          call write_details_mode_(items(i), wv, wr)
       end select

       select type (ml)
       class is (evol_model_t)
          call write_details_evol_model_(items(i), ml, wr)
       end select

       ! Check whether an item was written

       if (wr%count() /= c+1) then
          write(ERROR_UNIT, *) 'Ignoring invalid summary item:', TRIM(items(i))
       endif

    end do item_loop

    ! Finish

    return

  contains

    subroutine write_details_wave_ (item, wv, wr)

      character(*), intent(in)       :: item
      class(wave_t), intent(in)      :: wv
      class(writer_t), intent(inout) :: wr

      type(grid_t) :: gr
      integer      :: k
      real(WP)     :: data_r(wv%n_k)
      complex(WP)  :: data_c(wv%n_k)

      ! Write the item from wave_t data

      gr = ml%grid()

      select case (item)

      case ('n')
          
         call wr%write('n', wv%n_k)

      case ('j')

         call wr%write('j', wv%j)

      case ('l')

         call wr%write('l', wv%l)

      case ('l_i')

         call wr%write('l_i', wv%l_i)

      case ('m')

         call wr%write('m', wv%m)

      case ('lambda')

         call wr%write('lambda', [(wv%lambda(k), k=1,wv%n_k)])

      case ('omega')

         call wr%write('omega', wv%omega)

      case ('omega_int')

         call wr%write('omega_int', wv%omega_int())

      case ('domega_rot')

         call wr%write('domega_rot', wv%domega_rot())

      case ('freq')

         call wr%write('freq', wv%freq(ot_p%freq_units, ot_p%freq_frame))

      case ('dfreq_rot')

         call wr%write('dfreq_rot', wv%dfreq_rot(ot_p%freq_units))

      case ('freq_units')

         call wr%write('freq_units', ot_p%freq_units)

      case ('freq_frame')

         call wr%write('freq_frame', ot_p%freq_frame)

      case ('Delta_p')

         call wr%write('Delta_p', ml%Delta_p(gr%x_i(), gr%x_o()))

      case ('Delta_g')

         call wr%write('Delta_g', ml%Delta_g(gr%x_i(), gr%x_o(), wv%l*(wv%l+1._WP)))

      case ('eta')

         call wr%write('eta', wv%eta())

      case ('f_T')

         call wr%write('f_T', wv%f_T())

      case ('f_g')

         call wr%write('f_g', wv%f_g())

      case ('psi_T')

         call wr%write('psi_T', wv%psi_T())

      case ('psi_g')

         call wr%write('psi_g', wv%psi_g())

      case ('E')

         call wr%write('E', wv%E())

      case ('E_p')

         call wr%write('E_p', wv%E_p())

      case ('E_g')

         call wr%write('E_g', wv%E_g())

      case ('E_norm')

         call wr%write('E_norm', wv%E_norm())

      case ('E_ratio')

         call wr%write('E_ratio', wv%E_ratio())

      case ('H')

         call wr%write('H', wv%H())

      case ('W')

         call wr%write('W', wv%W())

      case ('W_eps')

         call wr%write('W_eps', wv%W_eps())

      case ('tau_ss')

         call wr%write('tau_ss', wv%tau_ss())

      case ('tau_tr')

         call wr%write('tau_tr', wv%tau_tr())

      case ('beta')

         call wr%write('beta', wv%beta())

      case ('x')

         call wr%write('x', gr%pt%x)

      case ('x_ref')

         call wr%write('x_ref', gr%pt(wv%k_ref)%x)

      $OUTPUT_POINTS(r,V_2,ml,coeff(I_V_2, gr%pt(k)))
      $OUTPUT_POINTS(r,As,ml,coeff(I_AS, gr%pt(k)))
      $OUTPUT_POINTS(r,U,ml,coeff(I_U, gr%pt(k)))
      $OUTPUT_POINTS(r,c_1,ml,coeff(I_C_1, gr%pt(k)))
      $OUTPUT_POINTS(r,Gamma_1,ml,coeff(I_GAMMA_1, gr%pt(k)))
      $OUTPUT_POINTS(r,nabla,ml,coeff(I_NABLA, gr%pt(k)))
      $OUTPUT_POINTS(r,nabla_ad,ml,coeff(I_NABLA_AD, gr%pt(k)))
      $OUTPUT_POINTS(r,dnabla_ad,ml,dcoeff(I_NABLA_AD, gr%pt(k)))
      $OUTPUT_POINTS(r,delta,ml,coeff(I_DELTA, gr%pt(k)))
      $OUTPUT_POINTS(r,c_lum,ml,coeff(I_C_LUM, gr%pt(k)))
      $OUTPUT_POINTS(r,c_rad,ml,coeff(I_C_RAD, gr%pt(k)))
      $OUTPUT_POINTS(r,c_thn,ml,coeff(I_C_THN, gr%pt(k)))
      $OUTPUT_POINTS(r,c_thk,ml,coeff(I_C_THK, gr%pt(k)))
      $OUTPUT_POINTS(r,c_eps,ml,coeff(I_C_EPS, gr%pt(k)))
      $OUTPUT_POINTS(r,eps_rho,ml,coeff(I_EPS_RHO, gr%pt(k)))
      $OUTPUT_POINTS(r,eps_T,ml,coeff(I_EPS_T, gr%pt(k)))
      $OUTPUT_POINTS(r,kap_rho,ml,coeff(I_KAP_RHO, gr%pt(k)))
      $OUTPUT_POINTS(r,kap_T,ml,coeff(I_KAP_T, gr%pt(k)))
      $OUTPUT_POINTS(r,Omega_rot,ml,coeff(I_OMEGA_ROT, gr%pt(k)))

      $OUTPUT_POINTS(c,y_1,wv,y_i(1, k))
      $OUTPUT_POINTS(c,y_2,wv,y_i(2, k))
      $OUTPUT_POINTS(c,y_3,wv,y_i(3, k))
      $OUTPUT_POINTS(c,y_4,wv,y_i(4, k))
      $OUTPUT_POINTS(c,y_5,wv,y_i(5, k))
      $OUTPUT_POINTS(c,y_6,wv,y_i(6, k))
      $OUTPUT_POINTS(c,xi_r,wv,xi_r(k))
      $OUTPUT_POINTS(c,xi_h,wv,xi_h(k))
      $OUTPUT_POINTS(c,eul_phi,wv,eul_phi(k))
      $OUTPUT_POINTS(c,deul_phi,wv,deul_phi(k))
      $OUTPUT_POINTS(c,eul_P,wv,eul_P(k))
      $OUTPUT_POINTS(c,eul_rho,wv,eul_rho(k))
      $OUTPUT_POINTS(c,eul_T,wv,eul_T(k))
      $OUTPUT_POINTS(c,lag_P,wv,lag_P(k))
      $OUTPUT_POINTS(c,lag_rho,wv,lag_rho(k))
      $OUTPUT_POINTS(c,lag_T,wv,lag_T(k))
      $OUTPUT_POINTS(c,lag_S,wv,lag_S(k))
      $OUTPUT_POINTS(c,lag_L,wv,lag_L(k))
      $OUTPUT_POINTS(r,dE_dx,wv,dE_dx(k))
      $OUTPUT_POINTS(r,dW_dx,wv,dW_dx(k))
      $OUTPUT_POINTS(r,dW_eps_dx,wv,dW_eps_dx(k))
      $OUTPUT_POINTS(c,dzeta_dx,wv,dzeta_dx(k))
      $OUTPUT_POINTS(c,dzeta_dm,wv,dzeta_dm(k))
      $OUTPUT_POINTS(r,dbeta_dx,wv,dbeta_dx(k))
      $OUTPUT_POINTS(r,dtau_dx_ss,wv,dtau_dx_ss(k))
      $OUTPUT_POINTS(r,dtau_dx_tr,wv,dtau_dx_tr(k))
      $OUTPUT_POINTS(c,Yt_1,wv,Yt_1(k))
      $OUTPUT_POINTS(c,Yt_2,wv,Yt_2(k))
      $OUTPUT_POINTS(c,I_0,wv,I_0(k))
      $OUTPUT_POINTS(c,I_1,wv,I_1(k))
      $OUTPUT_POINTS(r,alpha_0,wv,alpha_0(k))
      $OUTPUT_POINTS(r,alpha_1,wv,alpha_1(k))
      $OUTPUT_POINTS(c,prop_type,wv,prop_type(k))

      $OUTPUT_REF(xi_r_ref,xi_r)
      $OUTPUT_REF(xi_h_ref,xi_h)
      $OUTPUT_REF(eul_phi_ref,eul_phi)
      $OUTPUT_REF(deul_phi_ref,deul_phi)
      $OUTPUT_REF(lag_S_ref,lag_S)
      $OUTPUT_REF(lag_L_ref,lag_L)

      end select

      ! Finish

      return
      
    end subroutine write_details_wave_

    !****
      
    subroutine write_details_mode_ (item, md, wr)

      character(*), intent(in)       :: item
      class(mode_t), intent(in)      :: md
      class(writer_t), intent(inout) :: wr

      ! Write the item from mode_t data

      select case (item)

      case ('n_p')

         call wr%write('n_p', md%n_p)

      case ('n_g')

         call wr%write('n_g', md%n_g)

      case ('n_pg')

         call wr%write('n_pg', md%n_pg)
         
      end select

      ! Finish

      return

    end subroutine write_details_mode_

    !****

    subroutine write_details_evol_model_ (item, ml, wr)

      character(*), intent(in)        :: item
      class(evol_model_t), intent(in) :: ml
      class(writer_t), intent(inout)  :: wr

      type(grid_t) :: gr
      integer      :: k
      
      ! Write the item from evol_model_t data

      gr = ml%grid()

      select case (item)

      case ('M_star')

         call wr%write('M_star', ml%M_star)

      case ('R_star')

         call wr%write('R_star', ml%R_star)

      case ('L_star')

         call wr%write('L_star', ml%L_star)

      case ('M_r')

         call wr%write('M_r', [(ml%M_r(gr%pt(k)), k=1,gr%n_k)])

      case ('P')

         call wr%write('P', [(ml%P(gr%pt(k)), k=1,gr%n_k)])

      case ('rho')

         call wr%write('rho', [(ml%rho(gr%pt(k)), k=1,gr%n_k)])

      case ('T')

         call wr%write('T', [(ml%T(gr%pt(k)), k=1,gr%n_k)])

      end select

      ! Finish

      return

    end subroutine write_details_evol_model_

  end subroutine write_details

  !****

  function filter_wave_ (wv, ot_p) result (filter_wave)

    class(wave_t), intent(in)   :: wv
    type(out_par_t), intent(in) :: ot_p
    logical                     :: filter_wave

    character(LEN(ot_p%details_filter_list)), allocatable :: filters(:)
    integer                                               :: i

    ! Decide whether to filter the wave

    filters = split_list(ot_p%details_filter_list, ',')

    filter_wave = .FALSE.

    item_loop : do i = 1, SIZE(filters)

       select case (filters(i))
       case ('stable')
          filter_wave = filter_wave .OR. AIMAG(wv%omega) <= 0._WP
       case ('unstable')
          filter_wave = filter_wave .OR. AIMAG(wv%omega) > 0._WP
       case default
          $ABORT(Unrecognized filter in details_filter_list)
       end select

    end do item_loop

    ! Finish

    return

  end function filter_wave_

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
