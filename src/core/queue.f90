module Queue

  implicit none
!***********************************************************************!
  Type Queue_i

    private

    logical :: full
    integer :: start = 1, end = 0, count = 0, sz
    integer, allocatable, dimension(:) :: data

  contains

    procedure :: init => QUEUE_init
    procedure :: push
    procedure :: pop
    procedure :: isEmpty

  end type Queue_i

!***********************************************************************!

contains

!***********************************************************************!
  subroutine QUEUE_init(this, n)

    implicit none
    class(Queue_i), intent(inout) :: this
    integer, intent(in) :: n

    allocate(this%data(n))

    this%sz = n

  end subroutine QUEUE_init
!***********************************************************************!
  subroutine push(this, id, success)

    implicit none
    class(Queue_i), intent(inout) :: this
    integer, intent(in) :: id
    logical, intent(out) :: success

    success = .false.

    if (this%count .ge. this%sz) then
      write(*,'(A)') 'Cluster Queue is full. Cannot push.'
      return
    end if

    !push id into queue
    this%end = mod(this%end, this%sz)+1
    this%data(this%end) = id
    this%count = this%count + 1

    success = .true.

  end subroutine push
!***********************************************************************!
  subroutine pop(this, result, success)

    implicit none
    class(Queue_i), intent(inout) :: this
    integer, intent(out) :: result
    logical, intent(out) :: success

    success = .false.

    if (this%count .le. 0) then
      write(*,'(A)') 'Cluster Queue is empty. Cannot pop.'
      return
    end if

    result = this%data(this%start)
    this%start = mod(this%start, this%sz)+1
    this%count = this%count - 1

    success = .true.

  end subroutine pop
!***********************************************************************!
  function isEmpty(this)

    implicit none
    class(Queue_i), intent(inout) :: this
    logical :: isEmpty

    isEmpty = this%count .eq. 0

  end function isEmpty

end module Queue