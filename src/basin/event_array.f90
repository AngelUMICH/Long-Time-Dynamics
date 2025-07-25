module EventArray

  use Event, only: EventType
  use TranEvent, only: TransientEventType

  implicit none

  private

  public :: EventArrayType, TransientArrayType

  !*******************************************************************!
  ! Type Declarations
  !*******************************************************************!

  ! Dynamic array for EventType
  Type EventArrayType
    private
    integer :: count = 0          ! Current number of elements
    integer :: capacity = 0       ! Current capacity
    type(EventType), allocatable :: data(:)

  contains
    procedure :: initialize => EventArray_initialize
    procedure :: add => EventArray_add
    procedure :: remove => EventArray_remove
    procedure :: get => EventArray_get
    procedure :: getCount => EventArray_getCount
    procedure :: getAllEvents => EventArray_getAllEvents
    procedure :: copyFrom => EventArray_copyFrom
    procedure :: resize => EventArray_resize
  end Type EventArrayType

  ! Dynamic array for TransientEventType
  Type TransientArrayType
    private
    integer :: count = 0          ! Current number of elements
    integer :: capacity = 0       ! Current capacity
    type(TransientEventType), allocatable :: data(:)

  contains
    procedure :: initialize => TransientArray_initialize
    procedure :: add => TransientArray_add
    procedure :: remove => TransientArray_remove
    procedure :: get => TransientArray_get
    procedure :: getCount => TransientArray_getCount
    procedure :: getAllEvents => TransientArray_getAllEvents
    procedure :: copyFrom => TransientArray_copyFrom
    procedure :: resize => TransientArray_resize
  end Type TransientArrayType

contains

  !*******************************************************************!
  ! EventArrayType Methods
  !*******************************************************************!

  subroutine EventArray_initialize(this, initialCapacity)
    class(EventArrayType), intent(inout) :: this
    integer, intent(in), optional :: initialCapacity

    integer :: initCap

    ! Clean up if already initialized
    if (allocated(this % data)) then
      deallocate (this % data)
    end if

    ! Set initial capacity (default 8)
    initCap = 8
    if (present(initialCapacity)) then
      initCap = max(1, initialCapacity) ! Ensure at least 1
    end if

    this % capacity = initCap
    this % count = 0
    allocate (this % data(initCap))
  end subroutine EventArray_initialize

  !*******************************************************************!
  subroutine EventArray_add(this, newEvent)
    class(EventArrayType), intent(inout) :: this
    class(EventType), intent(in) :: newEvent

    type(EventType), allocatable :: temp(:)
    integer :: newCapacity

    ! Initialize if not already done
    if (.not. allocated(this % data)) then
      call this % initialize()
    end if

    ! Grow array if needed using doubling strategy
    if (this % count >= this % capacity) then
      ! Save current data
      allocate (temp(this % capacity))
      temp = this % data

      ! Double capacity (common growth strategy for dynamic arrays)
      newCapacity = max(8, this % capacity * 2)

      ! Reallocate with new size
      deallocate (this % data)
      allocate (this % data(newCapacity))

      ! Copy old data
      this % data(1:this % count) = temp(1:this % count)
      this % capacity = newCapacity

      ! Clean up
      deallocate (temp)
    end if

    ! Add the new event
    this % count = this % count + 1
    this % data(this % count) = newEvent
  end subroutine EventArray_add

  !*******************************************************************!
  subroutine EventArray_remove(this, index)
    class(EventArrayType), intent(inout) :: this
    integer, intent(in) :: index

    integer :: i

    ! Check for valid index
    if (index < 1 .or. index > this % count) then
      return
    end if

    ! Shift elements to fill the gap
    do i = index, this % count - 1
      this % data(i) = this % data(i + 1)
    end do

    ! Decrement count
    this % count = this % count - 1

    ! Shrink array if it's significantly oversized (optional optimization)
    if (this % count > 0 .and. this % count < this % capacity / 4) then
      call this % resize(this % capacity / 2)
    end if
  end subroutine EventArray_remove

  !*******************************************************************!
  function EventArray_get(this, index) result(singleEvent)
    class(EventArrayType), intent(in) :: this
    integer, intent(in) :: index
    type(EventType) :: singleEvent

    ! Check for valid index
    if (index >= 1 .and. index <= this % count) then
      singleEvent = this % data(index)
    else
      write(*,'(A,I0,A,I0)') "EventArray Error: Event with index ", index, &
      " greater than array count ", this%count
    end if
  end function EventArray_get

  !*******************************************************************!
  function EventArray_getCount(this) result(count)
    class(EventArrayType), intent(in) :: this
    integer :: count

    count = this % count
  end function EventArray_getCount

  !*******************************************************************!
  function EventArray_getAllEvents(this) result(events)
    class(EventArrayType), intent(in) :: this
    type(EventType), allocatable :: events(:)

    if (this % count > 0) then
      allocate (events(this % count))
      events = this % data(1:this % count)
    else
      allocate (events(0))
    end if
  end function EventArray_getAllEvents

  !*******************************************************************!
  subroutine EventArray_copyFrom(this, source)
    class(EventArrayType), intent(inout) :: this
    class(EventArrayType), intent(in) :: source

    ! Deallocate current data if allocated
    if (allocated(this % data)) then
      deallocate (this % data)
    end if

    ! Copy attributes
    this % count = source % count
    this % capacity = source % capacity

    ! Allocate new memory and copy data
    allocate (this % data(this % capacity))
    this % data(1:this % count) = source % data(1:source % count)
    
  end subroutine EventArray_copyFrom

  !*******************************************************************!
  subroutine EventArray_resize(this, newCapacity)
    class(EventArrayType), intent(inout) :: this
    integer, intent(in) :: newCapacity

    type(EventType), allocatable :: temp(:)
    integer :: copyCount

    ! Ignore if capacity is already sufficient
    if (allocated(this % data) .and. this % capacity >= newCapacity) return

    ! Save existing data
    if (allocated(this % data) .and. this % count > 0) then
      allocate (temp(this % count))
      temp(1:this % count) = this % data(1:this % count)
      copyCount = this % count
    else
      copyCount = 0
    end if

    ! Reallocate with new capacity
    if (allocated(this % data)) deallocate (this % data)
    allocate (this % data(newCapacity))

    ! Copy data back
    if (copyCount > 0) then
      this % data(1:copyCount) = temp(1:copyCount)
    end if

    this % capacity = newCapacity

    ! Clean up
    if (allocated(temp)) deallocate (temp)
  end subroutine EventArray_resize

  !*******************************************************************!
  ! TransientArrayType Methods
  !*******************************************************************!

  subroutine TransientArray_initialize(this, initialCapacity)
    class(TransientArrayType), intent(inout) :: this
    integer, intent(in), optional :: initialCapacity

    integer :: initCap

    ! Clean up if already initialized
    if (allocated(this % data)) then
      deallocate (this % data)
    end if

    ! Set initial capacity (default 8)
    initCap = 8
    if (present(initialCapacity)) then
      initCap = max(1, initialCapacity) ! Ensure at least 1
    end if

    this % capacity = initCap
    this % count = 0
    allocate (this % data(initCap))
  end subroutine TransientArray_initialize

  !*******************************************************************!
  subroutine TransientArray_add(this, newEvent)
    class(TransientArrayType), intent(inout) :: this
    class(TransientEventType), intent(in) :: newEvent

    type(TransientEventType), allocatable :: temp(:)
    integer :: newCapacity

    ! Initialize if not already done
    if (.not. allocated(this % data)) then
      call this % initialize()
    end if

    ! Grow array if needed using doubling strategy
    if (this % count >= this % capacity) then
      ! Save current data
      allocate (temp(this % capacity))
      temp = this % data

      ! Double capacity (common growth strategy for dynamic arrays)
      newCapacity = max(8, this % capacity * 2)

      ! Reallocate with new size
      deallocate (this % data)
      allocate (this % data(newCapacity))

      ! Copy old data
      this % data(1:this % count) = temp(1:this % count)
      this % capacity = newCapacity

      ! Clean up
      deallocate (temp)
    end if

    ! Add the new event
    this % count = this % count + 1
    this % data(this % count) = newEvent
  end subroutine TransientArray_add

  !*******************************************************************!
  subroutine TransientArray_remove(this, index)
    class(TransientArrayType), intent(inout) :: this
    integer, intent(in) :: index

    integer :: i

    ! Check for valid index
    if (index < 1 .or. index > this % count) then
      return
    end if

    ! Shift elements to fill the gap
    do i = index, this % count - 1
      this % data(i) = this % data(i + 1)
    end do

    ! Decrement count
    this % count = this % count - 1
  end subroutine TransientArray_remove

  !*******************************************************************!
  function TransientArray_get(this, index) result(singleEvent)
    class(TransientArrayType), intent(in) :: this
    integer, intent(in) :: index
    type(TransientEventType) :: singleEvent

    ! Check for valid index
    if (index >= 1 .and. index <= this % count) then
      singleEvent = this % data(index)
    end if
  end function TransientArray_get

  !*******************************************************************!
  function TransientArray_getCount(this) result(count)
    class(TransientArrayType), intent(in) :: this
    integer :: count

    count = this % count
  end function TransientArray_getCount

  !*******************************************************************!
  function TransientArray_getAllEvents(this) result(events)
    class(TransientArrayType), intent(in) :: this
    type(TransientEventType), allocatable :: events(:)

    if (this % count > 0) then
      allocate (events(this % count))
      events = this % data(1:this % count)
    else
      allocate (events(0))
    end if
  end function TransientArray_getAllEvents

  !*******************************************************************!
  subroutine TransientArray_copyFrom(this, source)
    class(TransientArrayType), intent(inout) :: this
    class(TransientArrayType), intent(in) :: source

    ! Deallocate current data if allocated
    if (allocated(this % data)) then
      deallocate (this % data)
    end if

    ! Copy attributes
    this % count = source % count
    this % capacity = source % capacity

    ! Allocate new memory and copy data
    allocate (this % data(this % capacity))
    this % data(1:this % count) = source % data(1:source % count)

  end subroutine TransientArray_copyFrom

  !*******************************************************************!
  subroutine TransientArray_resize(this, newCapacity)
    class(TransientArrayType), intent(inout) :: this
    integer, intent(in) :: newCapacity

    type(TransientEventType), allocatable :: temp(:)
    integer :: copyCount

    ! Ignore if capacity is already sufficient
    if (allocated(this % data) .and. this % capacity >= newCapacity) return

    ! Save existing data
    if (allocated(this % data) .and. this % count > 0) then
      allocate (temp(this % count))
      temp(1:this % count) = this % data(1:this % count)
      copyCount = this % count
    else
      copyCount = 0
    end if

    ! Reallocate with new capacity
    if (allocated(this % data)) deallocate (this % data)
    allocate (this % data(newCapacity))

    ! Copy data back
    if (copyCount > 0) then
      this % data(1:copyCount) = temp(1:copyCount)
    end if

    this % capacity = newCapacity

    ! Clean up
    if (allocated(temp)) deallocate (temp)
  end subroutine TransientArray_resize

end module EventArray
