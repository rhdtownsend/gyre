! Module   : astro_star
! Purpose  : centrifugally disorted star

$include 'core.inc'

module astro_star

  ! Uses

  use core_kinds
  use core_constants
  use core_parallel
  $if ($HDF5)
  use core_hgroup
  $endif
  use core_func

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type star_t
     real(WP) :: omega ! Omega/Omega_crit
   contains
     procedure :: Phi_eff => Phi_eff_
     procedure :: vg_eff => vg_eff_
     procedure :: vn_eff => vn_eff_
     procedure :: D_surf => D_surf_
     procedure :: Sigma_0 => Sigma_0_
     procedure :: Sigma_1 => Sigma_1_
     procedure :: A_surf => A_surf_
     procedure :: R_surf => R_surf_
     procedure :: s_surf => s_surf_
  end type star_t

  type, extends (func_t) :: discfunc_t
     class(star_t), pointer :: st
     real(WP)               :: vr_0(3)
     real(WP)               :: vdr(3)
   contains
     procedure, public :: eval_c_
  end type discfunc_t

  ! Interfaces

  interface star_t
     module procedure star_t_
  end interface star_t

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

  public :: star_t
  $if ($HDF5)
  public :: read
  public :: write
  $endif
  $if ($MPI)
  public :: bcast
  $endif

contains

  function star_t_ (omega) result (st)

    real(WP), intent(in) :: omega
    type(star_t)         :: st
    
    ! Construct the star_t

    st%omega = omega

    ! Finish

    return

  end function star_t_
  
!****

  function Phi_eff_ (this, vr) result (Phi_eff)

    class(star_t), intent(in) :: this
    real(WP), intent(in)      :: vr(:)
    real(WP)                  :: Phi_eff

    $CHECK_BOUNDS(SIZE(vr),3)

    ! Calculate the effective potential (units: G M_star/R_pole) at vr
    ! (units: R_pole)

    Phi_eff = -1._WP/NORM2(vr) - &
               4._WP*this%omega**2*(vr(1)**2+vr(2)**2)/27._WP

    ! Finish

    return

  end function Phi_eff_

!****

  function vg_eff_ (this, vr) result (vg_eff)

    class(star_t), intent(in) :: this
    real(WP), intent(in)      :: vr(:)
    real(WP)                  :: vg_eff(3)

    $CHECK_BOUNDS(SIZE(vr),3)

    ! Calculate the effective gravity vector (units: G M_star/R_pole**2) at vr
    ! (units: R_pole)

    vg_eff = -vr/NORM2(vr)**3 + &
            8._WP*this%omega**2*[vr(1),vr(2),0._WP]/27._WP

    ! Finish

    return

  end function vg_eff_

!****

  function vn_eff_ (this, vr) result (vn_eff)

    class(star_t), intent(in) :: this
    real(WP), intent(in)      :: vr(:)
    real(WP)                  :: vn_eff(3)

    real(WP) :: vg_eff(3)

    $CHECK_BOUNDS(SIZE(vr), 3)

    ! Calculate the equipotential surface normal (unit vector) at vr
    ! (units: R_pole). 

    vg_eff = this%vg_eff(vr)

    vn_eff = -vg_eff/NORM2(vg_eff)

    ! Finish

  end function vn_eff_

!****

  function D_surf_ (this, vr) result (D_surf)

    class(star_t), intent(in) :: this
    real(WP), intent(in)      :: vr(:)
    real(WP)                  :: D_surf

    real(WP) :: r

    ! Calculate the surface discriminant at vr; D_surf > 0 outside the
    ! star and < 0 inside the star

    r = NORM2(vr)

    if (r == 0._WP) then

       D_surf = -HUGE(0._WP)

    else

       if (r > 1.5_WP) then

          D_surf = HUGE(0._WP)

       else

          D_surf = 1._WP - 1._WP/r - 4._WP*(vr(1)**2 + vr(2)**2)*this%omega**2/27._WP

       endif

    endif

    ! Finish

    return

  end function D_surf_

!****

  function Sigma_0_ (this) result (Sigma_0)

    class(star_t), intent(in) :: this
    real(WP)                  :: Sigma_0

    ! Calculate the Sigma_0 surface integral (see Cranmer 1996, PhD
    ! thesis, Chap. 4)

    associate (omega => this%omega)

      Sigma_0 = 1._WP + 0.19444_WP*omega**2 + 0.28053_WP*omega**4 - &
                1.9014_WP*omega**6 + 6.8298_WP*omega**8 - 9.5002_WP*omega**10 + &
                4.6631_WP*omega**12

    end associate

    ! Finish

    return

  end function Sigma_0_

!****

  function Sigma_1_ (this) result (Sigma_1)

    class(star_t), intent(in) :: this
    real(WP)                  :: Sigma_1

    ! Calculate the Sigma_1 surface integral (see Cranmer 1996, PhD
    ! thesis, Chap. 4)

    associate (omega => this%omega)

      Sigma_1 = 1._WP - 0.19696_WP*omega**2 - 0.094292_WP*omega**4 + &
                0.33812_WP*omega**6 - 1.3066_WP*omega**8 + 1.8286_WP*omega**10 - &
                0.92714_WP*omega**12

    end associate

    ! Finish

    return

  end function Sigma_1_

!****

  function A_surf_ (this) result (A_surf)

    class(star_t), intent(in) :: this
    real(WP)                  :: A_surf

    ! Calculate the stellar surface area (units: R_pole**2)

    A_surf = 4.*WP*PI*this%Sigma_0()

    ! Finish

    return

  end function A_surf_

!****

  function R_surf_ (this, theta) result (R_surf)

    class(star_t), intent(in) :: this
    real(WP), intent(in)      :: theta
    real(WP)                  :: R_surf

    real(WP) :: ws

    ! Calculate the stellar radius (units: R_pole) at colatitude theta

    ws = this%omega*SIN(theta)

    if (ws**2/6._WP > EPSILON(0._WP)) then
       R_surf = 3._WP*COS((PI+ACOS(ws))/3._WP)/ws
    else
       R_surf = 1._WP
    endif

    ! Finish

    return

  end function R_surf_

!****

  function s_surf_ (this, vr_0, vdr) result (s_surf)

    class(star_t), target, intent(in) :: this
    real(WP), intent(in)              :: vr_0(:)
    real(WP), intent(in)              :: vdr(:)
    real(WP)                          :: s_surf

    real(WP)         :: s_c
    real(WP)         :: vr_c(3)
    real(WP)         :: r_c
    real(WP)         :: ds
    real(WP)         :: s_b1
    real(WP)         :: s_b2
    type(discfunc_t) :: df
    real(WP)         :: s_min
    real(WP)         :: vr_min(3)

    $CHECK_BOUNDS(SIZE(vr_0),3)
    $CHECK_BOUNDS(SIZE(vdr),3)

    ! Calculate the parameter s where the ray with parametric equation
    ! vr = vr_0 + vdr*s intersects the stellar surface. If there is no
    ! s > 0 intersection, then set s to -HUGE

    ! Determine the ray's point of closest approach

    s_c = -DOT_PRODUCT(vr_0, vdr)/DOT_PRODUCT(vdr, vdr)

    vr_c = vr_0 + vdr*s_c
    r_c = NORM2(vr_c)

    ! Quick check against the bounding sphere (r=1.5) for no
    ! intersection

    if (r_c >= 1.5) then
       s_surf = -HUGE(0._WP)
       return
    endif

    ! Find where the ray intersects this bounding sphere

    ds = SQRT(1.5_WP**2 - r_c**2)

    s_b1 = s_c - ds
    s_b2 = s_c + ds

    ! Find the point between these intersections where the potential
    ! is at a minimum

    df%st => this
    df%vr_0 = vr_0
    df%vdr = vdr

    s_min = df%minimum(s_b1, 0.5_WP*(s_b1+s_b2), s_b2, 0._WP)
    
    ! Check if this minimum is inside the equipotential surface

    vr_min = vr_0 + vdr*s_min

    if (this%Phi_eff(vr_min) < -1._WP) then

       ! Now hunt for the intersection on the equipotential surface

       if (s_b1 > 0._WP) then

          s_surf = df%root(s_b1, s_min, 0._WP)

       elseif (s_min > 0._WP) then

          s_surf = df%root(s_b1, s_min, 0._WP)

          if (s_surf < 0._WP) then
             s_surf = df%root(s_min, s_b2, 0._WP)
          endif

       elseif (s_b2 > 0._WP) then

          s_surf = df%root(s_min, s_b2, 0._WP)

          if (s_surf < 0._WP) then
             s_surf = -HUGE(0._WP)
          endif

       else

          s_surf = -HUGE(0._WP)
          
       endif

    else

       s_surf = -HUGE(0._WP)

    endif

    ! Finish

    return

  end function s_surf_

!****

  function eval_c_ (this, z) result (f_z)

    class(discfunc_t), intent(inout) :: this
    complex(WP), intent(in)          :: z
    complex(WP)                      :: f_z

    real(WP) :: vr(3)

    ! Evaluate the discriminant function defining the stellar surface

    associate (s => REAL(z))

      vr = this%vr_0 + this%vdr*s

      f_z = this%st%Phi_eff(vr) + 1._WP

    end associate

    ! Finish

    return

  end function eval_c_

!****

  $if ($HDF5)

  subroutine read_ (hg, st)

    type(hgroup_t), intent(inout) :: hg
    type(star_t), intent(out)     :: st

    real(WP) :: omega

    ! Read the star

    call read_attr(hg, 'omega', omega)

    st = star_t(omega)

    ! Finish

    return

  end subroutine read_

!****

  subroutine write_ (hg, st)

    type(hgroup_t), intent(inout) :: hg
    type(star_t), intent(in)      :: st

    ! Write the star

    call write_attr(hg, 'omega', st%omega)

    ! Finish

    return

  end subroutine write_

  $endif

!****

  $if ($MPI)

  subroutine bcast_ (st, root_rank)

    type(star_t), intent(inout) :: st
    integer, intent(in)         :: root_rank

    ! Broadcast the star

    call bcast(st%omega, root_rank)

    ! Finish

    return

  end subroutine bcast_

  $endif

end module astro_star
