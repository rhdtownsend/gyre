! Module   : core_list
! Purpose  : doubly-linked list
 
$include 'core.inc'

module core_list

  ! Uses

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: list_t
     private
     type(list_t), pointer :: li_next => NULL()
     type(list_t), pointer :: li_prev => NULL()
   contains
     private
     procedure, public :: next => next_
     procedure, public :: prev => prev_
     procedure, public :: head => head_
     procedure, public :: tail => tail_
     procedure, public :: insert_before => insert_before_
     procedure, public :: insert_after => insert_after_
     procedure, public :: delete => delete_
     procedure, public :: delete_all => delete_all_
  end type list_t

  ! Interfaces

  ! Access specifiers

  private

  public :: func_t

  ! Procedures

contains

  function prev (this, n)

    class(list_t), intent(in)     :: this
    integer, intent(in), optional :: n
    class(list_t), pointer        :: prev

    ! Get the n'th previous item before this

    if(PRESENT(n)) then
       i = MAX(n, 0)
    else
       i = 1
    endif

    prev => this%li_prev

    prev_loop : do

       if(i <= 0 .OR. .NOT. ASSOCIATED(prev)) exit prev_loop

       i = i - 1

       prev => prev%li_prev

    end do prev_loop

    ! Finish

    return

  end function prev

!****

  function next (this, n)

    class(list_t), intent(in)     :: this
    integer, intent(in), optional :: n
    class(list_t), pointer        :: next

    ! Get the n'th next item after this

    if(PRESENT(n)) then
       i = MAX(n, 0)
    else
       i = 1
    endif

    next => this%li_next

    next_loop : do

       if(i <= 0 .OR. .NOT. ASSOCIATED(next)) exit next_loop

       i = i - 1

       next => next%li_next

    end do next_loop

    ! Finish

    return

  end function next

!****

  function head (this)

    class(list_t), intent(in), target :: this
    class(list_t), pointer            :: head

    class(list_t), pointer :: prev

    ! Get the list head

    head => this

    head_loop : do

       prev => head%prev()

       if(.NOT. ASSOCIATED(prev)) exit head_loop

       head => prev

    end do head_loop

    ! Finish

    return

  end function head

!****

  function tail (this)

    class(list_t), intent(in), target :: this
    class(list_t), pointer            :: tail

    class(list_t), pointer :: next

    ! Get the list tail

    tail => this

    tail_loop : do

       next => tail%next()

       if(.NOT. ASSOCIATED(next)) exit tail_loop

       tail => next

    end do tail_loop

    ! Finish

    return

  end function tail

!****
       
  function insert_before (this) result (insert)

    class(list_t), intent(inout) :: this
    class(list_t), pointer       :: insert

    class(list_t), pointer :: prev

    ! Insert a list item before this

    prev => this%prev()

    allocate(insert, MOLD=this)

    this%li_prev => insert
    if(ASSOCIATED(prev)) prev%li_next => insert

    ! Finish

    return

  end function insert_before

!****
    
  function insert_after (this) result (insert)

    class(list_t), intent(inout) :: this
    class(list_t), pointer       :: insert

    class(list_t), pointer :: next

    ! Insert a list item after this

    next => this%next()

    allocate(insert, MOLD=this)

    this%li_next => insert
    if(ASSOCIATED(next)) next%li_prev => insert

    ! Finish

    return

  end function insert_after

!****

  subroutine delete (this)

    class(list_t), pointer :: this

    ! Delete the item

    prev => this%prev()
    next => this%next()

    deallocate(this)

    if(ASSOCIATED(prev)) prev%li_next => next
    if(ASSOCIATED(next)) next%li_prev => prev

    ! Finish

    return

  end subroutine delete

!****

  subroutine delete_all (this)

    class(list_t), pointer :: this

    class(list_t), pointer :: next

    ! Delete the whole list

    this => this%head()

    delete_loop : do

       next => this%next()

       call this%delete()

       this => next

       if(.NOT. ASSOCIATED(this)) exit delete_loop

    end do delete_loop

    ! Finish

    return

  end subroutine delete_all

end module core_list
