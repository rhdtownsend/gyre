! Module   : core_hgroup
! Purpose  : HDF5 input/output

$include 'core.inc'

$perl
xmacro 'HDF5_CALL', sub {
  my @a = @_; @a = get_args() unless @a;
  die("Invalid number of arguments:".scalar @a) unless scalar @a >= 1;
  my $proc_name = shift @a;
  @pos_args = grep(!/=/, @a);
  @key_args = grep(/=/, @a);
  @args = (@pos_args, 'hdf_err', @key_args);
  my $arg_list = join(', ', @args);
  return <<EOF;
call $proc_name($arg_list)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line $line_num $context_doc:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to $proc_name'
   stop 'Program aborted'
endif
EOF
};
$endperl

module core_hgroup

  ! Uses

  use core_kinds

  use hdf5

  use ISO_FORTRAN_ENV
  use ISO_C_BINDING

  ! No implicit typing

  implicit none

  ! Parameters

  integer, parameter :: CREATE_FILE = 1
  integer, parameter :: OPEN_FILE = 2

  ! Derived-type definitions

  type hgroup_t
     private
     integer(HID_T)   :: file_id = 0
     integer(HID_T)   :: group_id = 0
     integer, pointer :: ref_count => null()
   contains
     private
     procedure, public :: final => final_
     procedure, public :: exists => exists_
  end type hgroup_t

  ! Module variables

  integer, save        :: ref_count = 0
  integer(HID_T), save :: c_sp_mem_type_id
  integer(HID_T), save :: c_dp_mem_type_id
  integer(HID_T), save :: c_sp_file_type_id
  integer(HID_T), save :: c_dp_file_type_id
  
  ! Interfaces

  interface hgroup_t
     module procedure init_file_
     module procedure init_group_
  end interface hgroup_t

  interface read_attr
     module procedure read_attr_i_i4_0_
     module procedure read_attr_i_i4_1_
     module procedure read_attr_i_i4_2_
     module procedure read_attr_i_i4_3_
     module procedure read_attr_i_i4_4_
     module procedure read_attr_i_i8_0_
     module procedure read_attr_i_i8_1_
     module procedure read_attr_i_i8_2_
     module procedure read_attr_i_i8_3_
     module procedure read_attr_i_i8_4_
     module procedure read_attr_r_sp_0_
     module procedure read_attr_r_sp_1_
     module procedure read_attr_r_sp_2_
     module procedure read_attr_r_sp_3_
     module procedure read_attr_r_sp_4_
     module procedure read_attr_r_dp_0_
     module procedure read_attr_r_dp_1_
     module procedure read_attr_r_dp_2_
     module procedure read_attr_r_dp_3_
     module procedure read_attr_r_dp_4_
     module procedure read_attr_c_sp_0_
     module procedure read_attr_c_sp_1_
     module procedure read_attr_c_sp_2_
     module procedure read_attr_c_sp_3_
     module procedure read_attr_c_sp_4_
     module procedure read_attr_c_dp_0_
     module procedure read_attr_c_dp_1_
     module procedure read_attr_c_dp_2_
     module procedure read_attr_c_dp_3_
     module procedure read_attr_c_dp_4_
     module procedure read_attr_a_0_
     module procedure read_attr_a_1_
     module procedure read_attr_a_2_
     module procedure read_attr_a_3_
     module procedure read_attr_a_4_
     module procedure read_attr_l_0_
     module procedure read_attr_l_1_
     module procedure read_attr_l_2_
     module procedure read_attr_l_3_
     module procedure read_attr_l_4_
  end interface read_attr

  interface read_dset
     module procedure read_dset_i_i4_0_
     module procedure read_dset_i_i4_1_
     module procedure read_dset_i_i4_2_
     module procedure read_dset_i_i4_3_
     module procedure read_dset_i_i4_4_
     module procedure read_dset_i_i8_0_
     module procedure read_dset_i_i8_1_
     module procedure read_dset_i_i8_2_
     module procedure read_dset_i_i8_3_
     module procedure read_dset_i_i8_4_
     module procedure read_dset_r_sp_0_
     module procedure read_dset_r_sp_1_
     module procedure read_dset_r_sp_2_
     module procedure read_dset_r_sp_3_
     module procedure read_dset_r_sp_4_
     module procedure read_dset_r_dp_0_
     module procedure read_dset_r_dp_1_
     module procedure read_dset_r_dp_2_
     module procedure read_dset_r_dp_3_
     module procedure read_dset_r_dp_4_
     module procedure read_dset_c_sp_0_
     module procedure read_dset_c_sp_1_
     module procedure read_dset_c_sp_2_
     module procedure read_dset_c_sp_3_
     module procedure read_dset_c_sp_4_
     module procedure read_dset_c_dp_0_
     module procedure read_dset_c_dp_1_
     module procedure read_dset_c_dp_2_
     module procedure read_dset_c_dp_3_
     module procedure read_dset_c_dp_4_
     module procedure read_dset_a_0_
     module procedure read_dset_a_1_
     module procedure read_dset_a_2_
     module procedure read_dset_a_3_
     module procedure read_dset_a_4_
     module procedure read_dset_l_0_
     module procedure read_dset_l_1_
     module procedure read_dset_l_2_
     module procedure read_dset_l_3_
     module procedure read_dset_l_4_
  end interface read_dset

  interface read_attr_alloc
     module procedure read_attr_alloc_i_i4_0_
     module procedure read_attr_alloc_i_i4_1_
     module procedure read_attr_alloc_i_i4_2_
     module procedure read_attr_alloc_i_i4_3_
     module procedure read_attr_alloc_i_i4_4_
     module procedure read_attr_alloc_i_i8_0_
     module procedure read_attr_alloc_i_i8_1_
     module procedure read_attr_alloc_i_i8_2_
     module procedure read_attr_alloc_i_i8_3_
     module procedure read_attr_alloc_i_i8_4_
     module procedure read_attr_alloc_r_sp_0_
     module procedure read_attr_alloc_r_sp_1_
     module procedure read_attr_alloc_r_sp_2_
     module procedure read_attr_alloc_r_sp_3_
     module procedure read_attr_alloc_r_sp_4_
     module procedure read_attr_alloc_r_dp_0_
     module procedure read_attr_alloc_r_dp_1_
     module procedure read_attr_alloc_r_dp_2_
     module procedure read_attr_alloc_r_dp_3_
     module procedure read_attr_alloc_r_dp_4_
     module procedure read_attr_alloc_c_sp_0_
     module procedure read_attr_alloc_c_sp_1_
     module procedure read_attr_alloc_c_sp_2_
     module procedure read_attr_alloc_c_sp_3_
     module procedure read_attr_alloc_c_sp_4_
     module procedure read_attr_alloc_c_dp_0_
     module procedure read_attr_alloc_c_dp_1_
     module procedure read_attr_alloc_c_dp_2_
     module procedure read_attr_alloc_c_dp_3_
     module procedure read_attr_alloc_c_dp_4_
     module procedure read_attr_alloc_a_0_
     module procedure read_attr_alloc_a_1_
     module procedure read_attr_alloc_a_2_
     module procedure read_attr_alloc_a_3_
     module procedure read_attr_alloc_a_4_
     module procedure read_attr_alloc_l_0_
     module procedure read_attr_alloc_l_1_
     module procedure read_attr_alloc_l_2_
     module procedure read_attr_alloc_l_3_
     module procedure read_attr_alloc_l_4_
  end interface read_attr_alloc

  interface read_dset_alloc
     module procedure read_dset_alloc_i_i4_0_
     module procedure read_dset_alloc_i_i4_1_
     module procedure read_dset_alloc_i_i4_2_
     module procedure read_dset_alloc_i_i4_3_
     module procedure read_dset_alloc_i_i4_4_
     module procedure read_dset_alloc_i_i8_0_
     module procedure read_dset_alloc_i_i8_1_
     module procedure read_dset_alloc_i_i8_2_
     module procedure read_dset_alloc_i_i8_3_
     module procedure read_dset_alloc_i_i8_4_
     module procedure read_dset_alloc_r_sp_0_
     module procedure read_dset_alloc_r_sp_1_
     module procedure read_dset_alloc_r_sp_2_
     module procedure read_dset_alloc_r_sp_3_
     module procedure read_dset_alloc_r_sp_4_
     module procedure read_dset_alloc_r_dp_0_
     module procedure read_dset_alloc_r_dp_1_
     module procedure read_dset_alloc_r_dp_2_
     module procedure read_dset_alloc_r_dp_3_
     module procedure read_dset_alloc_r_dp_4_
     module procedure read_dset_alloc_c_sp_0_
     module procedure read_dset_alloc_c_sp_1_
     module procedure read_dset_alloc_c_sp_2_
     module procedure read_dset_alloc_c_sp_3_
     module procedure read_dset_alloc_c_sp_4_
     module procedure read_dset_alloc_c_dp_0_
     module procedure read_dset_alloc_c_dp_1_
     module procedure read_dset_alloc_c_dp_2_
     module procedure read_dset_alloc_c_dp_3_
     module procedure read_dset_alloc_c_dp_4_
     module procedure read_dset_alloc_a_0_
     module procedure read_dset_alloc_a_1_
     module procedure read_dset_alloc_a_2_
     module procedure read_dset_alloc_a_3_
     module procedure read_dset_alloc_a_4_
     module procedure read_dset_alloc_l_0_
     module procedure read_dset_alloc_l_1_
     module procedure read_dset_alloc_l_2_
     module procedure read_dset_alloc_l_3_
     module procedure read_dset_alloc_l_4_
  end interface read_dset_alloc
  
  interface write_attr
     module procedure write_attr_i_i4_0_
     module procedure write_attr_i_i4_1_
     module procedure write_attr_i_i4_2_
     module procedure write_attr_i_i4_3_
     module procedure write_attr_i_i4_4_
     module procedure write_attr_i_i8_0_
     module procedure write_attr_i_i8_1_
     module procedure write_attr_i_i8_2_
     module procedure write_attr_i_i8_3_
     module procedure write_attr_i_i8_4_
     module procedure write_attr_r_sp_0_
     module procedure write_attr_r_sp_1_
     module procedure write_attr_r_sp_2_
     module procedure write_attr_r_sp_3_
     module procedure write_attr_r_sp_4_
     module procedure write_attr_r_dp_0_
     module procedure write_attr_r_dp_1_
     module procedure write_attr_r_dp_2_
     module procedure write_attr_r_dp_3_
     module procedure write_attr_r_dp_4_
     module procedure write_attr_c_sp_0_
     module procedure write_attr_c_sp_1_
     module procedure write_attr_c_sp_2_
     module procedure write_attr_c_sp_3_
     module procedure write_attr_c_sp_4_
     module procedure write_attr_c_dp_0_
     module procedure write_attr_c_dp_1_
     module procedure write_attr_c_dp_2_
     module procedure write_attr_c_dp_3_
     module procedure write_attr_c_dp_4_
     module procedure write_attr_a_0_
     module procedure write_attr_a_1_
     module procedure write_attr_a_2_
     module procedure write_attr_a_3_
     module procedure write_attr_a_4_
     module procedure write_attr_l_0_
     module procedure write_attr_l_1_
     module procedure write_attr_l_2_
     module procedure write_attr_l_3_
     module procedure write_attr_l_4_
  end interface write_attr
  
  interface write_dset
     module procedure write_dset_i_i4_0_
     module procedure write_dset_i_i4_1_
     module procedure write_dset_i_i4_2_
     module procedure write_dset_i_i4_3_
     module procedure write_dset_i_i4_4_
     module procedure write_dset_i_i8_0_
     module procedure write_dset_i_i8_1_
     module procedure write_dset_i_i8_2_
     module procedure write_dset_i_i8_3_
     module procedure write_dset_i_i8_4_
     module procedure write_dset_r_sp_0_
     module procedure write_dset_r_sp_1_
     module procedure write_dset_r_sp_2_
     module procedure write_dset_r_sp_3_
     module procedure write_dset_r_sp_4_
     module procedure write_dset_r_dp_0_
     module procedure write_dset_r_dp_1_
     module procedure write_dset_r_dp_2_
     module procedure write_dset_r_dp_3_
     module procedure write_dset_r_dp_4_
     module procedure write_dset_c_sp_0_
     module procedure write_dset_c_sp_1_
     module procedure write_dset_c_sp_2_
     module procedure write_dset_c_sp_3_
     module procedure write_dset_c_sp_4_
     module procedure write_dset_c_dp_0_
     module procedure write_dset_c_dp_1_
     module procedure write_dset_c_dp_2_
     module procedure write_dset_c_dp_3_
     module procedure write_dset_c_dp_4_
     module procedure write_dset_a_0_
     module procedure write_dset_a_1_
     module procedure write_dset_a_2_
     module procedure write_dset_a_3_
     module procedure write_dset_a_4_
     module procedure write_dset_l_0_
     module procedure write_dset_l_1_
     module procedure write_dset_l_2_
     module procedure write_dset_l_3_
     module procedure write_dset_l_4_
  end interface write_dset
  
  interface attr_exists
     module procedure attr_exists_
  end interface attr_exists

  interface dset_exists
     module procedure dset_exists_
  end interface dset_exists

  ! Access specifiers

  private

  public :: CREATE_FILE
  public :: OPEN_FILE
  public :: hgroup_t
  public :: read_attr
  public :: read_dset
  public :: read_attr_alloc
  public :: read_dset_alloc
  public :: write_attr
  public :: write_dset
  public :: get_attr_shape
  public :: get_dset_shape
  public :: elem_group_name
  public :: attr_exists
  public :: dset_exists

  ! Procedures

contains

  function init_file_ (file_name, access_type) result (hg)

    character(*), intent(in) :: file_name
    integer, intent(in)      :: access_type
    type(hgroup_t)           :: hg

    integer        :: hdf_err
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id

    ! If necessary, open the HDF5 library

    if(ref_count == 0) then
       call open_library_()
    endif

    ref_count = ref_count + 1

    ! Depending on the access_type, open or create the file

    select case(access_type)
    case(CREATE_FILE)
       $HDF5_CALL(h5fcreate_f, file_name, H5F_ACC_TRUNC_F, file_id)
    case(OPEN_FILE)
       $HDF5_CALL(h5fopen_f, file_name, H5F_ACC_RDWR_F, file_id)
    case default
       $ABORT(Invalid access_type)
    end select

    ! Open the root group

    $HDF5_CALL(h5gopen_f, file_id, '/', group_id)

    ! Construct the hgroup_t

    hg%file_id = file_id
    hg%group_id = group_id

    allocate(hg%ref_count)
    hg%ref_count = 1

    ! Finish

    return

  end function init_file_

!****

  function init_group_ (hg_parent, group_name) result (hg)

    type(hgroup_t), intent(inout) :: hg_parent
    character(*), intent(in)      :: group_name
    type(hgroup_t)                :: hg

    integer        :: hdf_err
    integer(HID_T) :: group_id

    ! Depending on whether the group already exists, open or create it

    if(hg_parent%exists(group_name)) then
       $HDF5_CALL(h5gopen_f, hg_parent%group_id, group_name, group_id)
    else
       $HDF5_CALL(h5gcreate_f, hg_parent%group_id, group_name, group_id)
    endif

    ! Construct the hgroup_t

    hg%file_id = hg_parent%file_id
    hg%group_id = group_id

    hg%ref_count => hg_parent%ref_count
    hg%ref_count = hg%ref_count + 1

    ! Finish

    return

  end function init_group_

!****

  subroutine open_library_ ()

    integer :: hdf_err

    ! Open the HDF5 library

    $HDF5_CALL(h5open_f)

    ! Create complex data types

    call create_complex_type_(real_mem_type_(SP), c_sp_mem_type_id)
    call create_complex_type_(real_mem_type_(DP), c_dp_mem_type_id)

    call create_complex_type_(real_mem_type_(SP), c_sp_file_type_id)
    call create_complex_type_(real_mem_type_(DP), c_dp_file_type_id)

    ! Finish

    return

  contains

    subroutine create_complex_type_ (comp_type_id, type_id)

      integer(HID_T), intent(in)  :: comp_type_id
      integer(HID_T), intent(out) :: type_id

      integer         :: hdf_err
      integer(SIZE_T) :: comp_size

      ! Create a complex data type

      $HDF5_CALL(h5tget_size_f, comp_type_id, comp_size)

      $HDF5_CALL(h5tcreate_f, H5T_COMPOUND_F, INT(2*comp_size, SIZE_T), type_id)
      $HDF5_CALL(h5tinsert_f, type_id, 're', INT(0, SIZE_T), comp_type_id)
      $HDF5_CALL(h5tinsert_f, type_id, 'im', INT(comp_size, SIZE_T), comp_type_id)

      ! Finish

      return

    end subroutine create_complex_type_

  end subroutine open_library_

!****

  subroutine final_ (this)

    class(hgroup_t), intent(inout) :: this

    integer :: hdf_err

    ! Close the group

    $HDF5_CALL(h5gclose_f, this%group_id)

    this%ref_count = this%ref_count - 1
    
    ! If necessary, close the file also

    if(this%ref_count == 0) then
       $HDF5_CALL(h5fclose_f, this%file_id)
       deallocate(this%ref_count)
       ref_count = ref_count - 1
    endif

    ! If necessary, close the HDF5 library

    if(ref_count == 0) then
       call close_library_()
    endif

    ! Finish

    return

  end subroutine final_

!****

  subroutine close_library_ ()

    integer :: hdf_err

    ! Close complex data types

    $HDF5_CALL(h5tclose_f, c_sp_mem_type_id)
    $HDF5_CALL(h5tclose_f, c_dp_mem_type_id)

    $HDF5_CALL(h5tclose_f, c_sp_file_type_id)
    $HDF5_CALL(h5tclose_f, c_dp_file_type_id)

    ! Close the HDF5 library

    $HDF5_CALL(h5close_f)

    ! Finish

    return

  end subroutine close_library_

!****

  function exists_ (this, group_name)

    class(hgroup_t), intent(inout) :: this
    character(*), intent(in)       :: group_name
    logical                        :: exists_

    integer        :: hdf_err
    integer(HID_T) :: group_id

    ! Determine whether the named group already exists

    call h5eset_auto_f(0, hdf_err)

    call h5gopen_f(this%group_id, group_name, group_id, hdf_err)

    if(hdf_err >= 0) then
       exists_ = .TRUE.
       call h5gclose_f(group_id, hdf_err)
    else
       exists_ = .FALSE.
    endif

    call h5eset_auto_f(1, hdf_err)

    ! Finish

    return

  end function exists_

!****

  $define $READ_ITEM $sub

  $local $ITEM_TYPE $1
  $local $INFIX $2
  $local $DATA_TYPE $3
  $local $DATA_KIND $4
  $local $DATA_RANK $5

  subroutine read_${ITEM_TYPE}_${INFIX}_${DATA_RANK}_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    $DATA_TYPE($DATA_KIND), target, intent(out) :: data$ARRAY_SPEC($DATA_RANK)
    $if($DATA_RANK > 0)
    contiguous :: data
    $endif

    $if($DATA_RANK > 0)
    integer(HSIZE_T) :: item_shape($DATA_RANK)
    $endif
    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    $if($ITEM_TYPE eq 'attr')

    $define($OPEN_PROC,h5aopen_name_f)
    $define($READ_PROC,h5aread_f)
    $define($CLOSE_PROC,h5aclose_f)
    $define($SHAPE_PROC,get_attr_shape)

    $elsif($ITEM_TYPE eq 'dset')

    $define($OPEN_PROC,h5dopen_f)
    $define($READ_PROC,h5dread_f)
    $define($CLOSE_PROC,h5dclose_f)
    $define($SHAPE_PROC,get_dset_shape)

    $else

    $error(Invalid ITEM_TYPE)

    $endif

    ! Check the shape of the item agrees with that of data
    
    $if($DATA_RANK > 0)

    call ${SHAPE_PROC}(hg, item_name, item_shape)
    $ASSERT(ALL(item_shape == SHAPE(data)),Array shape mismatch)

    $endif

    ! Read the item

    mem_type_id = ${DATA_TYPE}_mem_type_($DATA_KIND)

    data_ptr = C_LOC(data)

    $HDF5_CALL($OPEN_PROC, hg%group_id, item_name, item_id)
    $HDF5_CALL($READ_PROC, item_id, mem_type_id, data_ptr)
    $HDF5_CALL($CLOSE_PROC, item_id)

    ! Finish

    return

  end subroutine read_${ITEM_TYPE}_${INFIX}_${DATA_RANK}_

  $endsub

  $READ_ITEM(attr,i_i4,integer,I4,0)
  $READ_ITEM(attr,i_i4,integer,I4,1)
  $READ_ITEM(attr,i_i4,integer,I4,2)
  $READ_ITEM(attr,i_i4,integer,I4,3)
  $READ_ITEM(attr,i_i4,integer,I4,4)

  $READ_ITEM(attr,i_i8,integer,I8,0)
  $READ_ITEM(attr,i_i8,integer,I8,1)
  $READ_ITEM(attr,i_i8,integer,I8,2)
  $READ_ITEM(attr,i_i8,integer,I8,3)
  $READ_ITEM(attr,i_i8,integer,I8,4)

  $READ_ITEM(attr,r_sp,real,SP,0)
  $READ_ITEM(attr,r_sp,real,SP,1)
  $READ_ITEM(attr,r_sp,real,SP,2)
  $READ_ITEM(attr,r_sp,real,SP,3)
  $READ_ITEM(attr,r_sp,real,SP,4)

  $READ_ITEM(attr,r_dp,real,DP,0)
  $READ_ITEM(attr,r_dp,real,DP,1)
  $READ_ITEM(attr,r_dp,real,DP,2)
  $READ_ITEM(attr,r_dp,real,DP,3)
  $READ_ITEM(attr,r_dp,real,DP,4)

  $READ_ITEM(attr,c_sp,complex,SP,0)
  $READ_ITEM(attr,c_sp,complex,SP,1)
  $READ_ITEM(attr,c_sp,complex,SP,2)
  $READ_ITEM(attr,c_sp,complex,SP,3)
  $READ_ITEM(attr,c_sp,complex,SP,4)

  $READ_ITEM(attr,c_dp,complex,DP,0)
  $READ_ITEM(attr,c_dp,complex,DP,1)
  $READ_ITEM(attr,c_dp,complex,DP,2)
  $READ_ITEM(attr,c_dp,complex,DP,3)
  $READ_ITEM(attr,c_dp,complex,DP,4)

  $READ_ITEM(dset,i_i4,integer,I4,0)
  $READ_ITEM(dset,i_i4,integer,I4,1)
  $READ_ITEM(dset,i_i4,integer,I4,2)
  $READ_ITEM(dset,i_i4,integer,I4,3)
  $READ_ITEM(dset,i_i4,integer,I4,4)

  $READ_ITEM(dset,i_i8,integer,I8,0)
  $READ_ITEM(dset,i_i8,integer,I8,1)
  $READ_ITEM(dset,i_i8,integer,I8,2)
  $READ_ITEM(dset,i_i8,integer,I8,3)
  $READ_ITEM(dset,i_i8,integer,I8,4)

  $READ_ITEM(dset,r_sp,real,SP,0)
  $READ_ITEM(dset,r_sp,real,SP,1)
  $READ_ITEM(dset,r_sp,real,SP,2)
  $READ_ITEM(dset,r_sp,real,SP,3)
  $READ_ITEM(dset,r_sp,real,SP,4)

  $READ_ITEM(dset,r_dp,real,DP,0)
  $READ_ITEM(dset,r_dp,real,DP,1)
  $READ_ITEM(dset,r_dp,real,DP,2)
  $READ_ITEM(dset,r_dp,real,DP,3)
  $READ_ITEM(dset,r_dp,real,DP,4)

  $READ_ITEM(dset,c_sp,complex,SP,0)
  $READ_ITEM(dset,c_sp,complex,SP,1)
  $READ_ITEM(dset,c_sp,complex,SP,2)
  $READ_ITEM(dset,c_sp,complex,SP,3)
  $READ_ITEM(dset,c_sp,complex,SP,4)

  $READ_ITEM(dset,c_dp,complex,DP,0)
  $READ_ITEM(dset,c_dp,complex,DP,1)
  $READ_ITEM(dset,c_dp,complex,DP,2)
  $READ_ITEM(dset,c_dp,complex,DP,3)
  $READ_ITEM(dset,c_dp,complex,DP,4)

!****

  $define $READ_ITEM_SPECIAL $sub

  $local $ITEM_TYPE $1
  $local $DATA_RANK $2

  subroutine read_${ITEM_TYPE}_a_${DATA_RANK}_ (hg, item_name, data)

    type(hgroup_t), intent(inout)     :: hg
    character(*), intent(in)          :: item_name
    character(*), target, intent(out) :: data$ARRAY_SPEC($DATA_RANK)
    $if($DATA_RANK > 0)
    contiguous :: data
    $endif

    $if($DATA_RANK > 0)
    integer(HSIZE_T) :: item_shape($DATA_RANK)
    $endif
    integer          :: hdf_err
    integer(HID_T)   :: mem_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    $if($ITEM_TYPE eq 'attr')

    $define($OPEN_PROC,h5aopen_name_f)
    $define($READ_PROC,h5aread_f)
    $define($CLOSE_PROC,h5aclose_f)
    $define($SHAPE_PROC,get_attr_shape)

    $elsif($ITEM_TYPE eq 'dset')

    $define($OPEN_PROC,h5dopen_f)
    $define($READ_PROC,h5dread_f)
    $define($CLOSE_PROC,h5dclose_f)
    $define($SHAPE_PROC,get_dset_shape)

    $else

    $error(Invalid ITEM_TYPE)

    $endif

    ! Check the shape of the item agrees with that of data
    
    $if($DATA_RANK > 0)

    call ${SHAPE_PROC}(hg, item_name, item_shape)
    $ASSERT(ALL(item_shape == SHAPE(data)),Array shape mismatch)

    $endif

    ! Read the character item

    $HDF5_CALL(h5tcopy_f, H5T_NATIVE_CHARACTER, mem_type_id)
    $HDF5_CALL(h5tset_size_f, mem_type_id, LEN(data, SIZE_T))

    data_ptr = C_LOC(data)

    $HDF5_CALL($OPEN_PROC, hg%group_id, item_name, item_id)
    $HDF5_CALL($READ_PROC, item_id, mem_type_id, data_ptr)
    $HDF5_CALL($CLOSE_PROC, item_id)

    $HDF5_CALL(h5tclose_f, mem_type_id)

    ! Finish

    return

  end subroutine read_${ITEM_TYPE}_a_${DATA_RANK}_

!****

  subroutine read_${ITEM_TYPE}_l_${DATA_RANK}_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(out)          :: data$ARRAY_SPEC($DATA_RANK)

    $if($DATA_RANK > 0)
    integer(HSIZE_T) :: item_shape($DATA_RANK)
    $endif
    integer          :: data_i$ARRAY_SPEC($DATA_RANK,data)

    $if($ITEM_TYPE eq 'attr')

    $define($SHAPE_PROC,get_attr_shape)

    $elsif($ITEM_TYPE eq 'dset')

    $define($SHAPE_PROC,get_dset_shape)

    $else

    $error(Invalid ITEM_TYPE)

    $endif

    ! Check the shape of the item agrees with that of data
    
    $if($DATA_RANK > 0)

    call ${SHAPE_PROC}(hg, item_name, item_shape)
    $ASSERT(ALL(item_shape == SHAPE(data)),Array shape mismatch)

    $endif

    ! Read the logical item

    call read_${ITEM_TYPE}(hg, item_name, data_i)

    data = data_i /= 0

    ! Finish

    return

  end subroutine read_${ITEM_TYPE}_l_${DATA_RANK}_

  $endsub

  $READ_ITEM_SPECIAL(attr,0)
  $READ_ITEM_SPECIAL(attr,1)
  $READ_ITEM_SPECIAL(attr,2)
  $READ_ITEM_SPECIAL(attr,3)
  $READ_ITEM_SPECIAL(attr,4)

  $READ_ITEM_SPECIAL(dset,0)
  $READ_ITEM_SPECIAL(dset,1)
  $READ_ITEM_SPECIAL(dset,2)
  $READ_ITEM_SPECIAL(dset,3)
  $READ_ITEM_SPECIAL(dset,4)

!****

  $define $READ_ITEM_ALLOC $sub

  $local $ITEM_TYPE $1
  $local $INFIX $2
  $local $DATA_TYPE $3
  $local $DATA_RANK $4

  subroutine read_${ITEM_TYPE}_alloc_${INFIX}_${DATA_RANK}_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    $DATA_TYPE, intent(inout), allocatable :: data$ARRAY_SPEC($DATA_RANK)

    $if($DATA_RANK > 0)
    integer(HSIZE_T) :: item_shape($DATA_RANK)
    $endif

    $if($ITEM_TYPE eq 'attr')

    $define($SHAPE_PROC,get_attr_shape)

    $elsif($ITEM_TYPE eq 'dset')

    $define($SHAPE_PROC,get_dset_shape)

    $else

    $error(Invalid ITEM_TYPE)

    $endif

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    $if($DATA_RANK > 0)

    call ${SHAPE_PROC}(hg, item_name, item_shape)
    allocate(data($ARRAY_EXPAND(item_shape,$DATA_RANK)))

    $else

    allocate(data)

    $endif

    ! Read the item

    call read_${ITEM_TYPE}(hg, item_name, data)

    ! Finish

    return

  end subroutine read_${ITEM_TYPE}_alloc_${INFIX}_${DATA_RANK}_

  $endsub

  $READ_ITEM_ALLOC(attr,i_i4,integer(I4),0)
  $READ_ITEM_ALLOC(attr,i_i4,integer(I4),1)
  $READ_ITEM_ALLOC(attr,i_i4,integer(I4),2)
  $READ_ITEM_ALLOC(attr,i_i4,integer(I4),3)
  $READ_ITEM_ALLOC(attr,i_i4,integer(I4),4)

  $READ_ITEM_ALLOC(attr,i_i8,integer(I8),0)
  $READ_ITEM_ALLOC(attr,i_i8,integer(I8),1)
  $READ_ITEM_ALLOC(attr,i_i8,integer(I8),2)
  $READ_ITEM_ALLOC(attr,i_i8,integer(I8),3)
  $READ_ITEM_ALLOC(attr,i_i8,integer(I8),4)

  $READ_ITEM_ALLOC(attr,r_sp,real(SP),0)
  $READ_ITEM_ALLOC(attr,r_sp,real(SP),1)
  $READ_ITEM_ALLOC(attr,r_sp,real(SP),2)
  $READ_ITEM_ALLOC(attr,r_sp,real(SP),3)
  $READ_ITEM_ALLOC(attr,r_sp,real(SP),4)

  $READ_ITEM_ALLOC(attr,r_dp,real(DP),0)
  $READ_ITEM_ALLOC(attr,r_dp,real(DP),1)
  $READ_ITEM_ALLOC(attr,r_dp,real(DP),2)
  $READ_ITEM_ALLOC(attr,r_dp,real(DP),3)
  $READ_ITEM_ALLOC(attr,r_dp,real(DP),4)

  $READ_ITEM_ALLOC(attr,c_sp,complex(SP),0)
  $READ_ITEM_ALLOC(attr,c_sp,complex(SP),1)
  $READ_ITEM_ALLOC(attr,c_sp,complex(SP),2)
  $READ_ITEM_ALLOC(attr,c_sp,complex(SP),3)
  $READ_ITEM_ALLOC(attr,c_sp,complex(SP),4)

  $READ_ITEM_ALLOC(attr,c_dp,complex(DP),0)
  $READ_ITEM_ALLOC(attr,c_dp,complex(DP),1)
  $READ_ITEM_ALLOC(attr,c_dp,complex(DP),2)
  $READ_ITEM_ALLOC(attr,c_dp,complex(DP),3)
  $READ_ITEM_ALLOC(attr,c_dp,complex(DP),4)

  $READ_ITEM_ALLOC(attr,a,character(*),0)
  $READ_ITEM_ALLOC(attr,a,character(*),1)
  $READ_ITEM_ALLOC(attr,a,character(*),2)
  $READ_ITEM_ALLOC(attr,a,character(*),3)
  $READ_ITEM_ALLOC(attr,a,character(*),4)

  $READ_ITEM_ALLOC(attr,l,logical,0)
  $READ_ITEM_ALLOC(attr,l,logical,1)
  $READ_ITEM_ALLOC(attr,l,logical,2)
  $READ_ITEM_ALLOC(attr,l,logical,3)
  $READ_ITEM_ALLOC(attr,l,logical,4)

  $READ_ITEM_ALLOC(dset,i_i4,integer(I4),0)
  $READ_ITEM_ALLOC(dset,i_i4,integer(I4),1)
  $READ_ITEM_ALLOC(dset,i_i4,integer(I4),2)
  $READ_ITEM_ALLOC(dset,i_i4,integer(I4),3)
  $READ_ITEM_ALLOC(dset,i_i4,integer(I4),4)

  $READ_ITEM_ALLOC(dset,i_i8,integer(I8),0)
  $READ_ITEM_ALLOC(dset,i_i8,integer(I8),1)
  $READ_ITEM_ALLOC(dset,i_i8,integer(I8),2)
  $READ_ITEM_ALLOC(dset,i_i8,integer(I8),3)
  $READ_ITEM_ALLOC(dset,i_i8,integer(I8),4)

  $READ_ITEM_ALLOC(dset,r_sp,real(SP),0)
  $READ_ITEM_ALLOC(dset,r_sp,real(SP),1)
  $READ_ITEM_ALLOC(dset,r_sp,real(SP),2)
  $READ_ITEM_ALLOC(dset,r_sp,real(SP),3)
  $READ_ITEM_ALLOC(dset,r_sp,real(SP),4)

  $READ_ITEM_ALLOC(dset,r_dp,real(DP),0)
  $READ_ITEM_ALLOC(dset,r_dp,real(DP),1)
  $READ_ITEM_ALLOC(dset,r_dp,real(DP),2)
  $READ_ITEM_ALLOC(dset,r_dp,real(DP),3)
  $READ_ITEM_ALLOC(dset,r_dp,real(DP),4)

  $READ_ITEM_ALLOC(dset,c_sp,complex(SP),0)
  $READ_ITEM_ALLOC(dset,c_sp,complex(SP),1)
  $READ_ITEM_ALLOC(dset,c_sp,complex(SP),2)
  $READ_ITEM_ALLOC(dset,c_sp,complex(SP),3)
  $READ_ITEM_ALLOC(dset,c_sp,complex(SP),4)

  $READ_ITEM_ALLOC(dset,c_dp,complex(DP),0)
  $READ_ITEM_ALLOC(dset,c_dp,complex(DP),1)
  $READ_ITEM_ALLOC(dset,c_dp,complex(DP),2)
  $READ_ITEM_ALLOC(dset,c_dp,complex(DP),3)
  $READ_ITEM_ALLOC(dset,c_dp,complex(DP),4)

  $READ_ITEM_ALLOC(dset,a,character(*),0)
  $READ_ITEM_ALLOC(dset,a,character(*),1)
  $READ_ITEM_ALLOC(dset,a,character(*),2)
  $READ_ITEM_ALLOC(dset,a,character(*),3)
  $READ_ITEM_ALLOC(dset,a,character(*),4)

  $READ_ITEM_ALLOC(dset,l,logical,0)
  $READ_ITEM_ALLOC(dset,l,logical,1)
  $READ_ITEM_ALLOC(dset,l,logical,2)
  $READ_ITEM_ALLOC(dset,l,logical,3)
  $READ_ITEM_ALLOC(dset,l,logical,4)

!****

  $define $WRITE_ITEM $sub

  $local $ITEM_TYPE $1
  $local $INFIX $2
  $local $DATA_TYPE $3
  $local $DATA_KIND $4
  $local $DATA_RANK $5

  subroutine write_${ITEM_TYPE}_${INFIX}_${DATA_RANK}_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    $DATA_TYPE($DATA_KIND), target, intent(in) :: data$ARRAY_SPEC($DATA_RANK)
    $if($DATA_RANK > 0)
    contiguous :: data
    $endif
    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape($DATA_RANK)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    $if($ITEM_TYPE eq 'attr')

    $define($OPEN_PROC,h5aopen_f)
    $define($CREATE_PROC,h5acreate_f)
    $define($WRITE_PROC,h5awrite_f)
    $define($CLOSE_PROC,h5aclose_f)
    $define($SHAPE_PROC,get_attr_shape)
    $define($EXISTS_PROC,attr_exists_)

    $elsif($ITEM_TYPE eq 'dset')

    $define($OPEN_PROC,h5dopen_f)
    $define($CREATE_PROC,h5dcreate_f)
    $define($WRITE_PROC,h5dwrite_f)
    $define($CLOSE_PROC,h5dclose_f)
    $define($SHAPE_PROC,get_dset_shape)
    $define($EXISTS_PROC,dset_exists_)

    $else

    $error(Invalid ITEM_TYPE)

    $endif

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = $EXISTS_PROC(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call ${SHAPE_PROC}(hg, item_name, item_shape)
       $ASSERT(ALL(item_shape == SHAPE(data)),Array shape mismatch)

       mem_type_id = ${DATA_TYPE}_mem_type_($DATA_KIND)

       data_ptr = C_LOC(data)

       $HDF5_CALL($OPEN_PROC, hg%group_id, item_name, item_id)
       $HDF5_CALL($WRITE_PROC, item_id, mem_type_id, data_ptr)
       $HDF5_CALL($CLOSE_PROC, item_id)

    else

       $if($DATA_RANK > 0)
       $HDF5_CALL(h5screate_simple_f, $DATA_RANK, INT(SHAPE(data), HSIZE_T), dspace_id)
       $else
       $HDF5_CALL(h5screate_f, H5S_SCALAR_F, dspace_id)
       $endif

       $if($ITEM_TYPE eq 'attr')
       $HDF5_CALL(h5pcreate_f, H5P_ATTRIBUTE_CREATE_F, plist_id)
       $ASSERT(.NOT. PRESENT(comp_level),Attributes cannot be compressed)
       $else
       $HDF5_CALL(h5pcreate_f, H5P_DATASET_CREATE_F, plist_id)
       if (PRESENT(comp_level)) then
          $HDF5_CALL(h5pset_chunk_f, plist_id, $DATA_RANK, INT(SHAPE(data), HSIZE_T))
          $HDF5_CALL(h5pset_deflate_f, plist_id, comp_level)
       endif
       $endif
       
       mem_type_id = ${DATA_TYPE}_mem_type_($DATA_KIND)
       file_type_id = ${DATA_TYPE}_file_type_($DATA_KIND)

       data_ptr = C_LOC(data)

       $if ($ITEM_TYPE eq 'attr')
       $HDF5_CALL($CREATE_PROC, hg%group_id, item_name, file_type_id, dspace_id, item_id, acpl_id=plist_id)
       $else
       $HDF5_CALL($CREATE_PROC, hg%group_id, item_name, file_type_id, dspace_id, item_id, dcpl_id=plist_id)
       $endif
       
       $HDF5_CALL($WRITE_PROC, item_id, mem_type_id, data_ptr)
       $HDF5_CALL(h5pclose_f, plist_id)
       $HDF5_CALL($CLOSE_PROC, item_id)
       $HDF5_CALL(h5sclose_f, dspace_id)

    endif

    ! Finish

    return

  end subroutine write_${ITEM_TYPE}_${INFIX}_${DATA_RANK}_

  $endsub

  $WRITE_ITEM(attr,i_i4,integer,I4,0)
  $WRITE_ITEM(attr,i_i4,integer,I4,1)
  $WRITE_ITEM(attr,i_i4,integer,I4,2)
  $WRITE_ITEM(attr,i_i4,integer,I4,3)
  $WRITE_ITEM(attr,i_i4,integer,I4,4)

  $WRITE_ITEM(attr,i_i8,integer,I8,0)
  $WRITE_ITEM(attr,i_i8,integer,I8,1)
  $WRITE_ITEM(attr,i_i8,integer,I8,2)
  $WRITE_ITEM(attr,i_i8,integer,I8,3)
  $WRITE_ITEM(attr,i_i8,integer,I8,4)

  $WRITE_ITEM(attr,r_sp,real,SP,0)
  $WRITE_ITEM(attr,r_sp,real,SP,1)
  $WRITE_ITEM(attr,r_sp,real,SP,2)
  $WRITE_ITEM(attr,r_sp,real,SP,3)
  $WRITE_ITEM(attr,r_sp,real,SP,4)

  $WRITE_ITEM(attr,r_dp,real,DP,0)
  $WRITE_ITEM(attr,r_dp,real,DP,1)
  $WRITE_ITEM(attr,r_dp,real,DP,2)
  $WRITE_ITEM(attr,r_dp,real,DP,3)
  $WRITE_ITEM(attr,r_dp,real,DP,4)

  $WRITE_ITEM(attr,c_sp,complex,SP,0)
  $WRITE_ITEM(attr,c_sp,complex,SP,1)
  $WRITE_ITEM(attr,c_sp,complex,SP,2)
  $WRITE_ITEM(attr,c_sp,complex,SP,3)
  $WRITE_ITEM(attr,c_sp,complex,SP,4)

  $WRITE_ITEM(attr,c_dp,complex,DP,0)
  $WRITE_ITEM(attr,c_dp,complex,DP,1)
  $WRITE_ITEM(attr,c_dp,complex,DP,2)
  $WRITE_ITEM(attr,c_dp,complex,DP,3)
  $WRITE_ITEM(attr,c_dp,complex,DP,4)

  $WRITE_ITEM(dset,i_i4,integer,I4,0)
  $WRITE_ITEM(dset,i_i4,integer,I4,1)
  $WRITE_ITEM(dset,i_i4,integer,I4,2)
  $WRITE_ITEM(dset,i_i4,integer,I4,3)
  $WRITE_ITEM(dset,i_i4,integer,I4,4)

  $WRITE_ITEM(dset,i_i8,integer,I8,0)
  $WRITE_ITEM(dset,i_i8,integer,I8,1)
  $WRITE_ITEM(dset,i_i8,integer,I8,2)
  $WRITE_ITEM(dset,i_i8,integer,I8,3)
  $WRITE_ITEM(dset,i_i8,integer,I8,4)

  $WRITE_ITEM(dset,r_sp,real,SP,0)
  $WRITE_ITEM(dset,r_sp,real,SP,1)
  $WRITE_ITEM(dset,r_sp,real,SP,2)
  $WRITE_ITEM(dset,r_sp,real,SP,3)
  $WRITE_ITEM(dset,r_sp,real,SP,4)

  $WRITE_ITEM(dset,r_dp,real,DP,0)
  $WRITE_ITEM(dset,r_dp,real,DP,1)
  $WRITE_ITEM(dset,r_dp,real,DP,2)
  $WRITE_ITEM(dset,r_dp,real,DP,3)
  $WRITE_ITEM(dset,r_dp,real,DP,4)

  $WRITE_ITEM(dset,c_sp,complex,SP,0)
  $WRITE_ITEM(dset,c_sp,complex,SP,1)
  $WRITE_ITEM(dset,c_sp,complex,SP,2)
  $WRITE_ITEM(dset,c_sp,complex,SP,3)
  $WRITE_ITEM(dset,c_sp,complex,SP,4)

  $WRITE_ITEM(dset,c_dp,complex,DP,0)
  $WRITE_ITEM(dset,c_dp,complex,DP,1)
  $WRITE_ITEM(dset,c_dp,complex,DP,2)
  $WRITE_ITEM(dset,c_dp,complex,DP,3)
  $WRITE_ITEM(dset,c_dp,complex,DP,4)

!****

  $define $WRITE_ITEM_SPECIAL $sub

  $local $ITEM_TYPE $1
  $local $DATA_RANK $2

  subroutine write_${ITEM_TYPE}_a_${DATA_RANK}_ (hg, item_name, data)

    type(hgroup_t), intent(inout)    :: hg
    character(*), intent(in)         :: item_name
    character(*), target, intent(in) :: data$ARRAY_SPEC($DATA_RANK)
    $if($DATA_RANK > 0)
    contiguous :: data
    $endif

    integer        :: hdf_err
    integer(HID_T) :: mem_type_id
    integer(HID_T) :: file_type_id
    integer(HID_T) :: dspace_id
    type(C_PTR)    :: data_ptr
    integer(HID_T) :: item_id

    $if($ITEM_TYPE eq 'attr')

    $define($CREATE_PROC,h5acreate_f)
    $define($WRITE_PROC,h5awrite_f)
    $define($CLOSE_PROC,h5aclose_f)

    $elsif($ITEM_TYPE eq 'dset')

    $define($CREATE_PROC,h5dcreate_f)
    $define($WRITE_PROC,h5dwrite_f)
    $define($CLOSE_PROC,h5dclose_f)

    $else

    $error(Invalid ITEM_TYPE)

    $endif

    ! Write the character item
    
    $HDF5_CALL(h5tcopy_f, H5T_NATIVE_CHARACTER, mem_type_id)
    $HDF5_CALL(h5tset_size_f, mem_type_id, LEN(data, SIZE_T))

    $HDF5_CALL(h5tcopy_f, H5T_NATIVE_CHARACTER, file_type_id)
    $HDF5_CALL(h5tset_size_f, file_type_id, LEN(data, SIZE_T))

    $if($DATA_RANK > 0)
    $HDF5_CALL(h5screate_simple_f, $DATA_RANK, INT(SHAPE(data), HSIZE_T), dspace_id)
    $else
    $HDF5_CALL(h5screate_f, H5S_SCALAR_F, dspace_id)
    $endif

    data_ptr = C_LOC(data)

    $HDF5_CALL($CREATE_PROC, hg%group_id, item_name, file_type_id, dspace_id, item_id)
    $HDF5_CALL($WRITE_PROC, item_id, mem_type_id, data_ptr)
    $HDF5_CALL($CLOSE_PROC, item_id)

    $HDF5_CALL(h5sclose_f, dspace_id)

    $HDF5_CALL(h5tclose_f, mem_type_id)
    $HDF5_CALL(h5tclose_f, file_type_id)

    ! Finish

    return

  end subroutine write_${ITEM_TYPE}_a_${DATA_RANK}_

!****

  subroutine write_${ITEM_TYPE}_l_${DATA_RANK}_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(in)           :: data$ARRAY_SPEC($DATA_RANK)

    ! Write the logical item

    call write_${ITEM_TYPE}(hg, item_name, MERGE(1, 0, MASK=data))

    ! Finish

    return

  end subroutine write_${ITEM_TYPE}_l_${DATA_RANK}_

  $endsub

  $WRITE_ITEM_SPECIAL(attr,0)
  $WRITE_ITEM_SPECIAL(attr,1)
  $WRITE_ITEM_SPECIAL(attr,2)
  $WRITE_ITEM_SPECIAL(attr,3)
  $WRITE_ITEM_SPECIAL(attr,4)

  $WRITE_ITEM_SPECIAL(dset,0)
  $WRITE_ITEM_SPECIAL(dset,1)
  $WRITE_ITEM_SPECIAL(dset,2)
  $WRITE_ITEM_SPECIAL(dset,3)
  $WRITE_ITEM_SPECIAL(dset,4)

!****

  subroutine get_attr_shape (hg, attr_name, shape)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: attr_name
    integer(HSIZE_T), intent(out) :: shape(:)

    integer          :: hdf_err
    integer(HID_T)   :: attr_id
    integer(HID_T)   :: space_id
    integer          :: rank
    integer(HSIZE_T) :: max_shape(SIZE(shape))

    ! Get the shape of the attribute

    $HDF5_CALL(h5aopen_name_f, hg%group_id, attr_name, attr_id)

    $HDF5_CALL(h5aget_space_f, attr_id, space_id)

    $HDF5_CALL(h5sget_simple_extent_ndims_f, space_id, rank)
    $ASSERT(rank == SIZE(shape),Rank mismatch)

    $HDF5_CALL(h5sget_simple_extent_dims_f, space_id, shape, max_shape)

    $HDF5_CALL(h5sclose_f, space_id)
    $HDF5_CALL(h5aclose_f, attr_id)

    ! Finish

    return

  end subroutine get_attr_shape

!****

  subroutine get_dset_shape (hg, dset_name, shape)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: dset_name
    integer(HSIZE_T), intent(out) :: shape(:)

    integer          :: hdf_err
    integer(HID_T)   :: dset_id
    integer(HID_T)   :: space_id
    integer          :: rank
    integer(HSIZE_T) :: max_shape(SIZE(shape))

    ! Get the shape of the dataset

    $HDF5_CALL(h5dopen_f, hg%group_id, dset_name, dset_id)

    $HDF5_CALL(h5dget_space_f, dset_id, space_id)

    $HDF5_CALL(h5sget_simple_extent_ndims_f, space_id, rank)
    $ASSERT(rank == SIZE(shape),Rank mismatch)

    $HDF5_CALL(h5sget_simple_extent_dims_f, space_id, shape, max_shape)

    $HDF5_CALL(h5sclose_f, space_id)
    $HDF5_CALL(h5dclose_f, dset_id)

    ! Finish

    return

  end subroutine get_dset_shape

!****

  function attr_exists_ (hg, attr_name) result (attr_exists)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: attr_name
    logical                       :: attr_exists

    integer :: hdf_err

    ! Check if the attribute exists

    $HDF5_CALL(h5aexists_f, hg%group_id, attr_name, attr_exists)

    ! Finish

    return

  end function attr_exists_

!****

  function dset_exists_ (hg, dset_name) result (dset_exists)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: dset_name
    logical                       :: dset_exists

    integer          :: hdf_err
    logical          :: link_exists
    type(h5o_info_t) :: obj_info

    ! Check if the dataset exists

    $HDF5_CALL(h5lexists_f, hg%group_id, dset_name, link_exists)

    if(link_exists) then

       $HDF5_CALL(h5oget_info_by_name_f, hg%group_id, dset_name, obj_info)

       dset_exists = obj_info%type == H5O_TYPE_DATASET_F

    else

       dset_exists = .FALSE.

    endif

    ! Finish

    return

  end function dset_exists_

!****

  function integer_mem_type_ (kind) result (mem_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: mem_type_id

    ! Determine the memory type for integers

    mem_type_id = h5kind_to_type(kind, H5_INTEGER_KIND)

    ! Finish

    return

  end function integer_mem_type_

!****

  function real_mem_type_ (kind) result (mem_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: mem_type_id

    ! Determine the memory type for reals

    mem_type_id = h5kind_to_type(kind, H5_REAL_KIND)

    ! Finish

    return

  end function real_mem_type_

!****

  function complex_mem_type_ (kind) result (mem_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: mem_type_id

    ! Determine the memory type for complexes

    select case (kind)
    case (SP)
       mem_type_id = c_sp_mem_type_id
    case (DP)
       mem_type_id = c_dp_mem_type_id
    case default
       $ABORT(Unsupported kind)
    end select

    ! Finish

    return

  end function complex_mem_type_

!****

  function integer_file_type_ (kind) result (file_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: file_type_id

    ! Determine the file type for integers

    select case (kind)
    case (I4)
       file_type_id = H5T_STD_I32LE
    case (I8)
       file_type_id = H5T_STD_I64LE
    case default
       $ABORT(Unsupported kind)
    end select

    ! Finish

    return

  end function integer_file_type_

!****

  function real_file_type_ (kind) result (file_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: file_type_id

    ! Determine the file type for reals

    select case (kind)
    case (SP)
       file_type_id = H5T_IEEE_F32LE
    case (DP)
       file_type_id = H5T_IEEE_F64LE
    case default
       $ABORT(Unsupported kind)
    end select

    ! Finish

    return

  end function real_file_type_

!****

  function complex_file_type_ (kind) result (file_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: file_type_id

    ! Determine the file type for complexes

    select case (kind)
    case (SP)
       file_type_id = c_sp_file_type_id
    case (DP)
       file_type_id = c_dp_file_type_id
    case default
       $ABORT(Unsupported kind)
    end select

    ! Finish

    return

  end function complex_file_type_

!****

  function elem_group_name (prefix, indices) result (group_name)

    character(*), intent(in)  :: prefix
    integer, intent(in)       :: indices(:)
    character(:), allocatable :: group_name

    integer                   :: n_indices
    character(:), allocatable :: name_format
    integer                   :: name_len
    integer                   :: i

    ! Set up an array-element group name

    n_indices = SIZE(indices)

    select case(n_indices)
    case(0)
       name_format = '(A,''()'')'
    case(1)
       name_format = '(A,''('',I0,'')'')'
    case default
       name_format = '(A,''('''//REPEAT('I0,'',''', n_indices-1)//'I0,'')'')'
    end select

    name_len = LEN_TRIM(prefix) + n_indices + 1

    do i = 1,SIZE(indices)
       if(indices(i) < 0) then
          name_len = name_len + FLOOR(LOG10(REAL(ABS(indices(i))))) + 2
       elseif(indices(i) > 0) then
          name_len = name_len + FLOOR(LOG10(REAL(indices(i)))) + 1
       else
          name_len = name_len + 1
       endif
    end do

    allocate(character(name_len) :: group_name)

    write(group_name, name_format) TRIM(prefix), indices

    ! Finish

    return

  end function elem_group_name

end module core_hgroup
