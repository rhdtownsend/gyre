#! Include : memory
#! Purpose : memory management macros

#:include 'core.inc'

#!!
#! Reallocate an array
#!!

#:def REALLOCATE(type, suffix, rank)

   #:if type == 'character(*)'
     #:set type_new = 'character(LEN(array))'
   #:else
     #:set type_new = type
   #:endif

   subroutine reallocate_${suffix}$_${rank}$_(array, shape, stat)

      ${type}$, allocatable, intent(inout) :: array${ARRAY_SPEC(rank)}$
      integer, intent(in)                  :: shape(:)
      integer, intent(out), optional       :: stat
      
      ${type_new}$, allocatable :: array_new${ARRAY_SPEC(rank)}$
      integer, allocatable      :: i_a(:)
      integer, allocatable      :: i_b(:)
      integer, allocatable      :: i_c(:)

      @:ASSERT_DEBUG(SIZE(shape) == ${rank}$, 'rank mismatch')

      ! Allocate the new array

      i_a = LBOUND(array)
      i_b = i_a + shape - 1

      if (PRESENT(stat)) then
         allocate(array_new${ARRAY_SPEC(rank, 'i_a({0:d}):i_b({0:d})')}$, STAT=stat)
         if (stat /= 0) return
      else
         allocate(array_new${ARRAY_SPEC(rank, 'i_a({0:d}):i_b({0:d})')}$)
      endif

      ! Copy the data from old to new

      i_c = MIN(i_b, UBOUND(array))

      array_new${ARRAY_SPEC(rank, 'i_a({0:d}):i_c({0:d})')}$ = array${ARRAY_SPEC(rank, 'i_a({0:d}):i_c({0:d})')}$

      ! Move the allocation

      call MOVE_ALLOC(array_new, array)

      ! Finish

      return

   end subroutine reallocate_${suffix}$_${rank}$_

#:enddef
