! Module   : core_volume
! Purpose  : volumetric data

$include 'core.inc'

module core_volume

  ! Uses

  use core_kinds
  $if($HDF5)
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
     real(WP)                     :: x_0
     real(WP)                     :: y_0
     real(WP)                     :: z_0
     real(WP)                     :: dx
     real(WP)                     :: dy
     real(WP)                     :: dz
     integer                      :: n_x
     integer                      :: n_y
     integer                      :: n_z
     character(LEN=NODE_TYPE_LEN) :: node_type
   contains
     procedure :: init
     $if($HDF5)
     procedure :: read
     procedure :: write
     $endif
     procedure :: nodes_shape
     procedure :: r_node
  end type volume_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_vl
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: NODE_TYPE_LEN
  public :: volume_t
  $if($MPI)
  public :: bcast
  $endif

contains

  subroutine init (this, n_x, n_y, n_z, node_type)

    class(volume_t), intent(out) :: this
    integer, intent(in)          :: n_x
    integer, intent(in)          :: n_y
    integer, intent(in)          :: n_z
    character(LEN=*), intent(in) :: node_type

    ! Initialize the volume

    this%n_x = n_x
    this%n_y = n_y
    this%n_z = n_z

    this%node_type = node_type

    ! Finish

    return

  end subroutine init

!****

  $if($HDF5)

  subroutine read (this, hg)

    class(volume_t), intent(out)  :: this
    type(hgroup_t), intent(inout) :: hg

    integer                      :: n_x
    integer                      :: n_y
    integer                      :: n_z
    character(LEN=NODE_TYPE_LEN) :: node_type

    ! Read the volume

    call read_attr(hg, 'n_x', n_x)
    call read_attr(hg, 'n_y', n_y)
    call read_attr(hg, 'n_z', n_z)

    call read_attr(hg, 'node_type', node_type)

    call this%init(n_x, n_y, n_z, node_type)

    call read_attr(hg, 'x_0', this%x_0)
    call read_attr(hg, 'y_0', this%y_0)
    call read_attr(hg, 'z_0', this%z_0)
    
    call read_attr(hg, 'dx', this%dx)
    call read_attr(hg, 'dy', this%dy)
    call read_attr(hg, 'dz', this%dz)
    
    ! Finish

    return

  end subroutine read

!****

  subroutine write (this, hg)

    class(volume_t), intent(in)   :: this
    type(hgroup_t), intent(inout) :: hg

    ! Write the volume

    call write_attr(hg, 'n_x', this%n_x)
    call write_attr(hg, 'n_y', this%n_y)
    call write_attr(hg, 'n_z', this%n_z)

    call write_attr(hg, 'x_0', this%x_0)
    call write_attr(hg, 'y_0', this%y_0)
    call write_attr(hg, 'z_0', this%z_0)
    
    call write_attr(hg, 'dx', this%dx)
    call write_attr(hg, 'dy', this%dy)
    call write_attr(hg, 'dz', this%dz)

    call write_attr(hg, 'node_type', this%node_type)
    
    ! Finish

    return

  end subroutine write

  $endif

!****

  $if($MPI)

  subroutine bcast_vl (vl, root_rank)

    type(volume_t), intent(inout) :: vl
    integer, intent(in)           :: root_rank

    ! Broadcast the volume

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

  end subroutine bcast_vl
  
  $endif

!****

  function nodes_shape (this)

    class(volume_t) :: this
    integer         :: nodes_shape(3)

    ! Return the shape to use for node data arrays

    select case (this%node_type)
    case ('CELL')
       nodes_shape = [this%n_x,this%n_y,this%n_z]
    case ('VERT')
       nodes_shape = [this%n_x,this%n_y,this%n_z] + 1
    case default
       $ABORT(Invalid node_type)
    end select

    ! Finish

    return

  end function nodes_shape
 
!****

  function r_node (this, i_x, i_y, i_z)

    class(volume_t), intent(in) :: this
    integer, intent(in)         :: i_x
    integer, intent(in)         :: i_y
    integer, intent(in)         :: i_z
    real(WP)                    :: r_node(3)

    ! Calculate the node coordinates

    select case (this%node_type)
    case ('CELL')
       r_node(1) = this%x_0 + (i_x-0.5_WP)*this%dx
       r_node(2) = this%y_0 + (i_y-0.5_WP)*this%dy
       r_node(3) = this%z_0 + (i_z-0.5_WP)*this%dz
    case ('VERT')
       r_node(1) = this%x_0 + (i_x-1)*this%dx
       r_node(2) = this%y_0 + (i_y-1)*this%dy
       r_node(3) = this%z_0 + (i_z-1)*this%dz
    case default
       $ABORT(Invalid node_type)
    end select

    ! Finish

    return

  end function r_node

end module core_volume
