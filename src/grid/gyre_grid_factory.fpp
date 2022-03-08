! Module   : gyre_grid_factory
! Purpose  : factory procedures for grid_t type
!
! Copyright 2013-2022 Rich Townsend & The GYRE Team
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

module gyre_grid_factory

  ! Uses

  use core_kinds

  use gyre_context
  use gyre_grid
  use gyre_grid_par
  use gyre_grid_spec
  use gyre_grid_util
  use gyre_math
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_state
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface grid_t
     module procedure grid_t_context_
     module procedure grid_t_grid_spec_
  end interface grid_t

  ! Access specifiers

  private

  public :: grid_t

  ! Procedures

contains

  function grid_t_context_ (cx, omega, gr_p, os_p) result (gr)

    type(context_t), pointer, intent(in) :: cx
    real(WP), intent(in)                 :: omega(:)
    type(grid_par_t), intent(in)         :: gr_p
    type(osc_par_t), intent(in)          :: os_p
    type(grid_t)                         :: gr

    type(grid_spec_t) :: gs

    ! Construct the grid_t using the supplied context

    gs = grid_spec_t(cx, omega)

    gr = grid_t([gs], gr_p, os_p)

    ! Finish

    return

  end function grid_t_context_

  !****

  function grid_t_grid_spec_ (gs, gr_p, os_p) result (gr)

    type(grid_spec_t), intent(in) :: gs(:)
    type(grid_par_t), intent(in)  :: gr_p
    type(osc_par_t), intent(in)   :: os_p
    type(grid_t)                  :: gr

    class(model_t), pointer :: ml => null()
    integer                 :: i
    integer                 :: s

    $ASSERT_DEBUG(SIZE(gs) >= 1,Too few grid_specs)

    ! Construct the grid_t using the supplied array of
    ! grid_specs. Each grid spec supplies a context and the set of
    ! frequencies to adopt for that context

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) 'Building spatial grid'
100    format(A)
    endif

    ! Check that the grid_specs are all associated with the same model

    $ASSERT_DEBUG(ASSOCIATED(gs(1)%cx),Null pointer)

    ml => gs(1)%cx%model()

    check_loop : do i = 1, SIZE(gs)
       $ASSERT_DEBUG(ASSOCIATED(gs(i)%cx),Null pointer)
       $ASSERT_DEBUG(ASSOCIATED(gs(i)%cx%model(), ml),Contexts are associated with different models)
    end do check_loop

    ! Create the scaffold grid

    gr = grid_t(ml%grid(), gr_p%x_i, gr_p%x_o)

    ! Add points

    call add_points_(gs, gr_p, os_p, gr)

    ! Report 

    if (check_log_level('INFO')) then

       write(OUTPUT_UNIT, 120) 'Final grid has', gr%s_o()-gr%s_i()+1, 'segment(s) and', gr%n_p, 'point(s):'
120    format(3X,A,1X,I0,1X,A,1X,I0,1X,A)
       
       seg_loop : do s = gr%s_i(), gr%s_o()

          associate( &
               p_i => gr%p_s_i(s), &
               p_o => gr%p_s_o(s))

            write(OUTPUT_UNIT, 130) 'Segment', s, ': x range', gr%pt(p_i)%x, '->', gr%pt(p_o)%x, &
                 '(', p_i, '->', p_o, ')'
130         format(6X,A,1X,I0,1X,A,1X,F6.4,1X,A,1X,F6.4,1X,A,I0,1X,A,1X,I0,A)

          end associate

       end do seg_loop

       write(OUTPUT_UNIT, *)
        
    end if

    ! Finish

  end function grid_t_grid_spec_

  !****

  subroutine add_points_ (gs, gr_p, os_p, gr)

    type(grid_spec_t), intent(in) :: gs(:)
    type(grid_par_t), intent(in)  :: gr_p
    type(osc_par_t), intent(in)   :: os_p
    type(grid_t), intent(inout)   :: gr

    integer              :: n_p
    integer              :: i_iter
    logical, allocatable :: refine(:)
    integer              :: p
    integer              :: i
    real(WP)             :: dx
    logical, allocatable :: refine_new(:)
    integer              :: p_new

    ! Add points globally

    ! Initialize the refine array

    n_p = gr%n_p

    refine = gr%pt(2:)%s == gr%pt(:n_p-1)%s

    ! Iterate until no more points need be added
    
    iter_loop : do i_iter = 1, gr_p%n_iter_max

       ! Loop through grid subintervals

       !$OMP PARALLEL DO PRIVATE (dx, i) SCHEDULE (DYNAMIC)
       sub_loop : do p = 1, n_p-1

          ! Check/update if the subinterval should be considered for refinement

          if (refine(p)) then

             associate (pt_a => gr%pt(p), &
                        pt_b => gr%pt(p+1))

               dx = pt_b%x - pt_a%x

               if (dx > gr_p%dx_min) then

                  spec_loop : do i = 1, SIZE(gs)

                     associate (cx => gs(i)%cx, &
                                omega => gs(i)%omega)

                       if (pt_a%x /= 0._WP) then

                          refine(p) = refine_mech_(cx, pt_a, pt_b, omega, gr_p, os_p) .OR. &
                                      refine_therm_(cx, pt_a, pt_b, omega, gr_p) .OR. &
                                      refine_struct_(cx, pt_a, pt_b, gr_p) .OR. &
                                      dx > gr_p%dx_max

                       else

                          refine(p) = refine_center_(cx, pt_a, pt_b, omega, gr_p) .OR. &
                                      dx > gr_p%dx_max

                       end if

                     end associate

                     if (refine(p)) exit spec_loop

                  end do spec_loop

               else

                  refine(p) = .FALSE.

               endif

             end associate

          end if

       end do sub_loop

       if (check_log_level('INFO')) then
          write(OUTPUT_UNIT, 100) 'Refined', COUNT(refine), 'subinterval(s) in iteration', i_iter
100       format(3X,A,1X,I0,1X,A,1X,I0)
       endif

       ! Check for completion

       if (COUNT(refine) == 0) exit iter_loop

       ! Perform refinement

       gr = grid_t(gr, refine)

       ! Update the refine array

       allocate(refine_new(n_p + COUNT(refine) - 1))

       p_new = 1

       do p = 1, n_p-1

          refine_new(p_new) = refine(p)
          p_new = p_new + 1

          if (refine(p)) then

             refine_new(p_new) = refine(p)
             p_new = p_new + 1

          endif

       end do

       n_p = n_p + COUNT(refine)

       call MOVE_ALLOC(refine_new, refine)

    end do iter_loop

    ! Finish

    return

  end subroutine add_points_
  
  !****

  function refine_mech_ (cx, pt_a, pt_b, omega, gr_p, os_p) result (refine)

    type(context_t), intent(in)  :: cx
    type(point_t), intent(in)    :: pt_a
    type(point_t), intent(in)    :: pt_b
    real(WP), intent(in)         :: omega(:)
    type(grid_par_t), intent(in) :: gr_p
    type(osc_par_t), intent(in)  :: os_p
    logical                      :: refine

    type(point_t)   :: pt
    real(WP)        :: dlnx
    real(WP)        :: V
    real(WP)        :: As
    real(WP)        :: U
    real(WP)        :: c_1
    real(WP)        :: Gamma_1
    real(WP)        :: Omega_rot
    real(WP)        :: Omega_rot_i
    integer         :: j
    type(r_state_t) :: st
    real(WP)        :: omega_c
    real(WP)        :: lambda
    real(WP)        :: l_i
    real(WP)        :: c_0
    real(WP)        :: c_2
    real(WP)        :: c_4
    real(WP)        :: psi2
    real(WP)        :: chi_r
    real(WP)        :: chi_r_p
    real(WP)        :: chi_r_m
    real(WP)        :: chi_i

    ! Determine whether the subinterval [pt_a,pt_b] warrants
    ! refinement, via a local analysis of the mechanical parts of the
    ! oscillation equations
    
    if (gr_p%w_osc == 0._WP .AND. gr_p%w_exp == 0._WP) then

       ! No refinement necessary

       refine = .FALSE.

    else

       ! Evaluate coefficients at the midpoint
       
       pt%s = pt_a%s
       pt%x = 0.5_WP*(pt_a%x + pt_b%x)

       dlnx = LOG(pt_b%x) - LOG(pt_a%x)

       associate (ml => cx%model())

         V = ml%coeff(I_V_2, pt)*pt%x**2
         As = ml%coeff(I_AS, pt)
         U = ml%coeff(I_U, pt)
         c_1 = ml%coeff(I_C_1, pt)
         Gamma_1 = ml%coeff(I_GAMMA_1, pt)

       end associate

       Omega_rot = cx%Omega_rot(pt)
       Omega_rot_i = cx%Omega_rot(cx%point_i())

       ! Loop over omega, updating refine

       refine = .FALSE.

       omega_loop : do j = 1, SIZE(omega)

          st = r_state_t(omega(j))

          omega_c = cx%omega_c(Omega_rot, st)

          lambda = cx%lambda(Omega_rot, st)
          l_i = cx%l_e(Omega_rot_i, st)

          ! Calculate the propagation discriminant psi2

          c_4 = -4._WP*V/Gamma_1*c_1*os_p%alpha_gam
          c_2 = (As - V/Gamma_1 - U + 4._WP)**2 + 4._WP*V/Gamma_1*As*os_p%alpha_gam*os_p%alpha_pi + 4._WP*lambda
          c_0 = -4._WP*lambda*As/c_1*os_p%alpha_pi

          if (c_0 /= 0._WP) then
             psi2 = (c_4*omega_c**4 + c_2*omega_c**2 + c_0)/omega_c**2
          else
             psi2 = c_4*omega_c**2 + c_2
          endif

          ! Evaluate real and imaginary parts of the mechanical
          ! eigenvalue

          if (psi2 < 0._WP) then

             ! Propagation zone

             chi_r = 0.5_WP*(As + V/Gamma_1 - U + 2._WP - 2._WP*l_i)

             chi_i = 0.5_WP*sqrt(-psi2)
             
          else

             ! Evanescent zone; pick the steeper solution

             chi_r_p = 0.5_WP*(As + V/Gamma_1 - U + 2._WP - 2._WP*l_i + sqrt(psi2))
             chi_r_m = 0.5_WP*(As + V/Gamma_1 - U + 2._WP - 2._WP*l_i - sqrt(psi2))

             chi_r = MERGE(chi_r_p, chi_r_m, abs(chi_r_p) > abs(chi_r_m))

             chi_i = 0._WP

          endif

          ! Update refine

          refine = refine .OR. &
                   dlnx*gr_p%w_exp*abs(chi_r) > TWOPI .OR. &
                   dlnx*gr_p%w_osc*abs(chi_i) > TWOPI

          if (refine) exit omega_loop
          
       end do omega_loop

    end if

    ! Finish

    return

  end function refine_mech_

  !****

  function refine_therm_ (cx, pt_a, pt_b, omega, gr_p) result (refine)

    type(context_t), intent(in)  :: cx
    type(point_t), intent(in)    :: pt_a
    type(point_t), intent(in)    :: pt_b
    real(WP), intent(in)         :: omega(:)
    type(grid_par_t), intent(in) :: gr_p
    logical                      :: refine

    type(point_t)   :: pt
    real(WP)        :: V
    real(WP)        :: nabla
    real(WP)        :: c_rad
    real(WP)        :: c_thk
    real(WP)        :: Omega_rot
    real(WP)        :: tau
    integer         :: j
    type(r_state_t) :: st
    real(WP)        :: omega_c
    real(WP)        :: dlnx

    ! Determine whether the subinterval [pt_a,pt_b] warrants
    ! refinement, via a local analysis of the thermal parts of the
    ! oscillation equations

    if (gr_p%w_thm == 0._WP) then

       ! No refinement necessary

       refine = .FALSE.

    else

       ! Evaluate coefficients at the midpoint
       
       pt%s = pt_a%s
       pt%x = 0.5_WP*(pt_a%x + pt_b%x)

       dlnx = LOG(pt_b%x) - LOG(pt_a%x)

       associate (ml => cx%model())

         V = ml%coeff(I_V_2, pt)*pt%x**2
         nabla = ml%coeff(I_NABLA, pt)

         c_rad = ml%coeff(I_C_RAD, pt)
         c_thk = ml%coeff(I_C_THK, pt)

       end associate

       Omega_rot = cx%Omega_rot(pt)

       ! Loop over omega, updating refine

       refine = .FALSE.

       omega_loop : do j = 1, SIZE(omega)

          st = r_state_t(omega(j))

          omega_c = cx%omega_c(Omega_rot, st)

          ! Evaluate the thermal eigenvalue
         
          tau = sqrt(abs(V*nabla*omega_c*c_thk/c_rad))
 
          ! Update refine
          
          refine = refine .OR. &
                   dlnx*gr_p%w_thm*abs(tau) > TWOPI

          if (refine) exit omega_loop
          
       end do omega_loop

    end if

    ! Finish

    return

  end function refine_therm_

  !****

  function refine_struct_ (cx, pt_a, pt_b, gr_p) result (refine)

    type(context_t), intent(in)  :: cx
    type(point_t), intent(in)    :: pt_a
    type(point_t), intent(in)    :: pt_b
    type(grid_par_t), intent(in) :: gr_p
    logical                      :: refine

    type(point_t) :: pt
    real(WP)      :: dV_2
    real(WP)      :: dAs
    real(WP)      :: dU
    real(WP)      :: dc_1
    real(WP)      :: dGamma_1
    real(WP)      :: dlnx
 
    ! Determine whether the subinterval [pt_a,pt_b] warrants
    ! refinement, via a local structure gradient analysis
 
    if (gr_p%w_str == 0._WP) then

       ! No refinement necessary

       refine = .FALSE.

    else

       ! Evaluate coefficients at the midpoint
       
       pt%s = pt_a%s
       pt%x = 0.5_WP*(pt_a%x + pt_b%x)

       dlnx = LOG(pt_b%x) - LOG(pt_a%x)

       associate (ml => cx%model())

         dV_2 = ml%dcoeff(I_V_2, pt)
         dAs = ml%dcoeff(I_AS, pt)
         dU = ml%dcoeff(I_U, pt)
         dc_1 = ml%dcoeff(I_C_1, pt)
         dGamma_1 = ml%dcoeff(I_GAMMA_1, pt)

       end associate

       ! Decide on refinement

       if (pt_a%x /= 0._WP) then

          refine = dlnx*gr_p%w_str*abs(dV_2) > 1._WP .OR. &
                   dlnx*gr_p%w_str*abs(dAs) > 1._WP .OR. &
                   dlnx*gr_p%w_str*abs(dU) > 1._WP .OR. &
                   dlnx*gr_p%w_str*abs(dc_1) > 1._WP .OR. &
                   dlnx*gr_p%w_str*abs(dGamma_1) > 1._WP

       endif

    end if

    ! Finish

    return

  end function refine_struct_

  !****

  function refine_center_ (cx, pt_a, pt_b, omega, gr_p) result (refine)

    type(context_t), intent(in)  :: cx
    type(point_t), intent(in)    :: pt_a
    type(point_t), intent(in)    :: pt_b
    real(WP), intent(in)         :: omega(:)
    type(grid_par_t), intent(in) :: gr_p
    logical                      :: refine

    type(point_t)   :: pt
    real(WP)        :: V
    real(WP)        :: As
    real(WP)        :: U
    real(WP)        :: c_1
    real(WP)        :: Gamma_1
    real(WP)        :: Omega_rot
    real(WP)        :: Omega_rot_i
    integer         :: j
    type(r_state_t) :: st
    real(WP)        :: omega_c
    real(WP)        :: lambda
    real(WP)        :: l_i
    real(WP)        :: c_0
    real(WP)        :: c_2
    real(WP)        :: c_4
    real(WP)        :: psi2
    real(WP)        :: chi_r
    real(WP)        :: chi_i

    ! Determine whether the (center) subinterval [pt_a,pt_b] warrants
    ! refinement, via a local analysis of the mechanical parts of the
    ! oscillation equations
    
    if (.NOT. gr_p%resolve_ctr .OR. gr_p%w_ctr == 0._WP) then

       ! No refinement necessary

       refine = .FALSE.

    else

       ! Evaluate coefficients at the cell midpoint
       
       pt%s = pt_a%s
       pt%x = 0.5_WP*(pt_a%x + pt_b%x)

       associate (ml => cx%model())

         V = ml%coeff(I_V_2, pt)*pt%x**2
         As = ml%coeff(I_AS, pt)
         U = ml%coeff(I_U, pt)
         c_1 = ml%coeff(I_C_1, pt)
         Gamma_1 = ml%coeff(I_GAMMA_1, pt)

       end associate

       Omega_rot = cx%Omega_rot(pt)
       Omega_rot_i = cx%Omega_rot(cx%point_i())

       ! Loop over omega, updating refine

       refine = .FALSE.

       omega_loop : do j = 1, SIZE(omega)

          st = r_state_t(omega(j))

          omega_c = cx%omega_c(Omega_rot, st)

          lambda = cx%lambda(Omega_rot, st)
          l_i = cx%l_e(Omega_rot_i, st)

          ! Calculate the propagation discriminant psi2

          c_4 = -4._WP*V/Gamma_1*c_1
          c_2 = (As - V/Gamma_1 - U + 4._WP)**2 + 4._WP*V/Gamma_1*As + 4._WP*lambda
          c_0 = -4._WP*lambda*As/c_1

          if (c_0 /= 0._WP) then
             psi2 = (c_4*omega_c**4 + c_2*omega_c**2 + c_0)/omega_c**2
          else
             psi2 = c_4*omega_c**2 + c_2
          endif

          ! Evaluate real and imaginary parts of the mechanical
          ! eigenvalue

          if (psi2 < 0._WP) then

             ! Propagation zone

             chi_r = 0.5_WP*(As + V/Gamma_1 - U + 2._WP - 2._WP*l_i)

             chi_i = 0.5_WP*sqrt(-psi2)
             
          else

             ! Evanescent zone; pick the non-diverging solution

             chi_r = 0.5_WP*(As + V/Gamma_1 - U + 2._WP - 2._WP*l_i + sqrt(psi2))

             chi_i = 0._WP

          endif

          ! Update refine

          refine = refine .OR. &
                   chi_i /= 0._WP .OR. &
                   gr_p%w_ctr*abs(chi_r) > 1._WP

          if (refine) exit omega_loop
          
       end do omega_loop

    end if

    ! Finish

    return

  end function refine_center_

end module gyre_grid_factory
