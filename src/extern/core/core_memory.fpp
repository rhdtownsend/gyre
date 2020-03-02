! Module   : core_memory
! Purpose  : memory management

$include 'core.inc'
$include 'core_memory.inc'

module core_memory

  ! Uses

  use core_kinds

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Interfaces

  interface reallocate
     module procedure reallocate_i_i4_1
     module procedure reallocate_i_i4_2
     module procedure reallocate_i_i4_3
     module procedure reallocate_i_i4_4
     module procedure reallocate_i_i8_1
     module procedure reallocate_i_i8_2
     module procedure reallocate_i_i8_3
     module procedure reallocate_i_i8_4
     module procedure reallocate_r_sp_1
     module procedure reallocate_r_sp_2
     module procedure reallocate_r_sp_3
     module procedure reallocate_r_sp_4
     module procedure reallocate_r_dp_1
     module procedure reallocate_r_dp_2
     module procedure reallocate_r_dp_3
     module procedure reallocate_r_dp_4
     module procedure reallocate_c_sp_1
     module procedure reallocate_c_sp_2
     module procedure reallocate_c_sp_3
     module procedure reallocate_c_sp_4
     module procedure reallocate_c_dp_1
     module procedure reallocate_c_dp_2
     module procedure reallocate_c_dp_3
     module procedure reallocate_c_dp_4
     module procedure reallocate_a_1
     module procedure reallocate_a_2
     module procedure reallocate_a_3
     module procedure reallocate_a_4
     module procedure reallocate_l_1
     module procedure reallocate_l_2
     module procedure reallocate_l_3
     module procedure reallocate_l_4
  end interface reallocate

  ! Access specifiers

  private

  public :: reallocate

  ! Procedures

contains

  $REALLOCATE(i_i4,integer(I4),1)
  $REALLOCATE(i_i4,integer(I4),2)
  $REALLOCATE(i_i4,integer(I4),3)
  $REALLOCATE(i_i4,integer(I4),4)

  $REALLOCATE(i_i8,integer(I8),1)
  $REALLOCATE(i_i8,integer(I8),2)
  $REALLOCATE(i_i8,integer(I8),3)
  $REALLOCATE(i_i8,integer(I8),4)

  $REALLOCATE(r_sp,real(SP),1)
  $REALLOCATE(r_sp,real(SP),2)
  $REALLOCATE(r_sp,real(SP),3)
  $REALLOCATE(r_sp,real(SP),4)

  $REALLOCATE(r_dp,real(DP),1)
  $REALLOCATE(r_dp,real(DP),2)
  $REALLOCATE(r_dp,real(DP),3)
  $REALLOCATE(r_dp,real(DP),4)

  $REALLOCATE(c_sp,complex(SP),1)
  $REALLOCATE(c_sp,complex(SP),2)
  $REALLOCATE(c_sp,complex(SP),3)
  $REALLOCATE(c_sp,complex(SP),4)

  $REALLOCATE(c_dp,complex(DP),1)
  $REALLOCATE(c_dp,complex(DP),2)
  $REALLOCATE(c_dp,complex(DP),3)
  $REALLOCATE(c_dp,complex(DP),4)

  $REALLOCATE(a,character(LEN=*),1)
  $REALLOCATE(a,character(LEN=*),2)
  $REALLOCATE(a,character(LEN=*),3)
  $REALLOCATE(a,character(LEN=*),4)

  $REALLOCATE(l,logical,1)
  $REALLOCATE(l,logical,2)
  $REALLOCATE(l,logical,3)
  $REALLOCATE(l,logical,4)

end module core_memory
