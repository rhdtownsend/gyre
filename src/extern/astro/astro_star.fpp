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

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type star_t
     real(WP) :: omega ! Omega/Omega_crit
   contains
     procedure :: Phi_eff => Phi_eff_
     procedure :: g_eff => g_eff_
     procedure :: n_eff => n_eff_
     procedure :: D_surf => D_surf_
     procedure :: Sigma_0 => Sigma_0_
     procedure :: Sigma_1 => Sigma_1_
     procedure :: A_surf => A_surf_
     procedure :: R_surf => R_surf_
  end type star_t

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

  function Phi_eff_ (this, r) result (Phi_eff)

    class(star_t), intent(in) :: this
    real(WP), intent(in)      :: r(:)
    real(WP)                  :: Phi_eff

    $CHECK_BOUNDS(SIZE(r),3)

    ! Calculate the effective potential (units: G M_star/R_pole) at r
    ! (units: R_pole)

    Phi_eff = -1._WP/NORM2(r) - &
               4._WP*this%omega**2*(r(1)**2+r(2)**2)/27._WP

    ! Finish

    return

  end function Phi_eff_

!****

  function g_eff_ (this, r) result (g_eff)

    class(star_t), intent(in) :: this
    real(WP), intent(in)      :: r(:)
    real(WP)                  :: g_eff(3)

    $CHECK_BOUNDS(SIZE(r),3)

    ! Calculate the effective gravity (units: G M_star/R_pole**2) at r
    ! (units: R_pole)

    g_eff = -r/NORM2(r)**3 + &
            8._WP*this%omega**2*[r(1)**2,r(2)**2,0._WP]/27._WP

    ! Finish

    return

  end function g_eff_

!****

  function n_eff_ (this, r) result (n_eff)

    class(star_t), intent(in) :: this
    real(WP), intent(in)      :: r(:)
    real(WP)                  :: n_eff(3)

    real(WP) :: g_eff(3)

    $CHECK_BOUNDS(SIZE(r), 3)

    ! Calculate the equipotential surface normal (unit vector) at r
    ! (units: R_pole). 

    g_eff = this%g_eff(r)

    n_eff = -g_eff/NORM2(g_eff)

    ! Finish

  end function n_eff_

!****

  function D_surf_ (this, r) result (D_surf)

    class(star_t), intent(in) :: this
    real(WP), intent(in)      :: r(:)
    real(WP)                  :: D_surf

    real(WP) :: nr

    ! Calculate the surface discriminant; D_surf > 0 outside the star
    ! and < 0 inside the star

    nr = NORM2(r)

    if(nr == 0._WP) then

       D_surf = -HUGE(0._WP)

    else

       if (nr > 1.5_WP) then

          D_surf = HUGE(0._WP)

       else

          D_surf = 1._WP - 1._WP/nr - 4._WP*(r(1)**2 + r(2)**2)*this%omega**2/27._WP

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

    if(ABS(ws) > EPSILON(0._WP)) then
       R_surf = 3._WP*COS((PI+ACOS(ws))/3._WP)/ws
    else
       R_surf = 1._WP
    endif

    ! Finish

    return

  end function R_surf_

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

  $if($MPI)

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
