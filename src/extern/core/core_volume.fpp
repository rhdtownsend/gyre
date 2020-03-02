! Module   : core_volume
! Purpose  : volumetric data

$include 'core.inc'

module core_volume

  ! Uses

  use core_kinds
  $if ($HDF5)
  use core_hgroup
  $endif
  use core_parallel

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: NODE_TYPE_LEN = 16

  ! Derived-type definitions

  type volume_t
     real(WP)                 :: x_0
     real(WP)                 :: y_0
     real(WP)                 :: z_0
     real(WP)                 :: dx
     real(WP)                 :: dy
     real(WP)                 :: dz
     integer                  :: n_x
     integer                  :: n_y
     integer                  :: n_z
     character(NODE_TYPE_LEN) :: node_type
   contains
     procedure :: nodes_shape => nodes_shape_
     procedure :: r_node => r_node_
  end type volume_t

  ! Interfaces

  interface volume_t
     module procedure init_
  end interface volume_t

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

  public :: NODE_TYPE_LEN
  public :: volume_t
  $if ($HDF5)
  public :: read
  public :: write
  $if ($MPI)
  public :: bcast
  $endif
  $endif

contains

  function init_ (n_x, n_y, n_z, node_type) result (vl)

    integer, intent(in)      :: n_x
    integer, intent(in)      :: n_y
    integer, intent(in)      :: n_z
    character(*), intent(in) :: node_type
    type(volume_t)           :: vl

    ! Construct the volume_t

    vl%n_x = n_x
    vl%n_y = n_y
    vl%n_z = n_z

    vl%node_type = node_type

    ! Finish

    return

  end function init_

!****

  $if ($HDF5)

  subroutine read_ (hg, vl)

    type(hgroup_t), intent(inout) :: hg
    type(volume_t), intent(out)   :: vl

    integer                  :: n_x
    integer                  :: n_y
    integer                  :: n_z
    character(NODE_TYPE_LEN) :: node_type

    ! Read the volume_t

    call read_attr(hg, 'n_x', n_x)
    call read_attr(hg, 'n_y', n_y)
    call read_attr(hg, 'n_z', n_z)

    call read_attr(hg, 'node_type', node_type)

    vl = volume_t(n_x, n_y, n_z, node_type)

    call read_attr(hg, 'x_0', vl%x_0)
    call read_attr(hg, 'y_0', vl%y_0)
    call read_attr(hg, 'z_0', vl%z_0)
    
    call read_attr(hg, 'dx', vl%dx)
    call read_attr(hg, 'dy', vl%dy)
    call read_attr(hg, 'dz', vl%dz)
    
    ! Finish

    return

  end subroutine read_

!****

  subroutine write_ (hg, vl)

    type(hgroup_t), intent(inout) :: hg
    type(volume_t), intent(in)    :: vl

    ! Write the volume_t

    call write_attr(hg, 'n_x', vl%n_x)
    call write_attr(hg, 'n_y', vl%n_y)
    call write_attr(hg, 'n_z', vl%n_z)

    call write_attr(hg, 'x_0', vl%x_0)
    call write_attr(hg, 'y_0', vl%y_0)
    call write_attr(hg, 'z_0', vl%z_0)
    
    call write_attr(hg, 'dx', vl%dx)
    call write_attr(hg, 'dy', vl%dy)
    call write_attr(hg, 'dz', vl%dz)

    call write_attr(hg, 'node_type', vl%node_type)
    
    ! Finish

    return

  end subroutine write_

  $endif

!****

  function nodes_shape_ (this)

    class(volume_t) :: this
    integer         :: nodes_shape_(3)

    ! Return the shape to use for node data arrays

    select case (this%node_type)
    case ('CELL')
       nodes_shape_ = [this%n_x,this%n_y,this%n_z]
    case ('VERT')
       nodes_shape_ = [this%n_x,this%n_y,this%n_z] + 1
    case default
       $ABORT(Invalid node_type)
    end select

    ! Finish

    return

  end function nodes_shape_
 
!****

  function r_node_ (this, i_x, i_y, i_z)

    class(volume_t), intent(in) :: this
    integer, intent(in)         :: i_x
    integer, intent(in)         :: i_y
    integer, intent(in)         :: i_z
    real(WP)                    :: r_node_(3)

    ! Calculate the node coordinates

    select case (this%node_type)
    case ('CELL')
       r_node_(1) = this%x_0 + (i_x-0.5_WP)*this%dx
       r_node_(2) = this%y_0 + (i_y-0.5_WP)*this%dy
       r_node_(3) = this%z_0 + (i_z-0.5_WP)*this%dz
    case ('VERT')
       r_node_(1) = this%x_0 + (i_x-1)*this%dx
       r_node_(2) = this%y_0 + (i_y-1)*this%dy
       r_node_(3) = this%z_0 + (i_z-1)*this%dz
    case default
       $ABORT(Invalid node_type)
    end select

    ! Finish

    return

  end function r_node_

!****

  $if ($MPI)

  subroutine bcast_ (vl, root_rank)

    type(volume_t), intent(inout) :: vl
    integer, intent(in)           :: root_rank

    ! Broadcast the volume_t

    call bcast(vl%n_x, root_rank)
    call bcast(vl%n_y, root_rank)
    call bcast(vl%n_z, root_rank)

    call bcast(vl%x_0, root_rank)
    call bcast(vl%y_0, root_rank)
    call bcast(vl%z_0, root_rank)

    call bcast(vl%dx, root_rank)
    call bcast(vl%dy, root_rank)
    call bcast(vl%dz, root_rank)

    call bcast(vl%node_type, root_rank)

    ! Finish

    return

  end subroutine bcast_
  
  $endif

end module core_volume
