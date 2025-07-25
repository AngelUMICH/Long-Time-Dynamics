module Basin

  use Event, only: EventType
  use EventArray, only: EventArrayType, TransientArrayType
  use State, only: StateType
  use TranEvent, only: TransientEventType

  implicit none

  private

  public :: BasinType

  !***********************************************************************!

  Type, extends(StateType) :: BasinType
    ! A basin is a state of the system that is associated with a high rate.
    ! The basin has absorbing and transient events. The absorbing events
    ! are those that may leave the basins. Transient events are those that
    ! lead to other basin states.

    private

    ! Dynamic arrays for events
    type(EventArrayType) :: absorbingEvents
    type(TransientArrayType) :: transientEvents

  contains

    procedure :: addEvent => BasinType_add_event
    procedure :: flagEvent => BasinType_flag_event

    ! Accessor methods
    procedure :: getAbsorbingCount => BasinType_get_absorbing_count
    procedure :: getTransientCount => BasinType_get_transient_count
    procedure :: getAbsorbingEvent => BasinType_get_absorbing_event
    procedure :: getTransientEvent => BasinType_get_transient_event
    procedure :: getAllAbsorbingEvents => BasinType_get_all_absorbing_events
    procedure :: getAllTransientEvents => BasinType_get_all_transient_events
    procedure :: copyFrom => BasinType_copy_from

  end Type BasinType

  !***********************************************************************!

contains

  !***********************************************************************!
  subroutine BasinType_add_event(this, event1, wasAdded)
    ! Add an event to the basin only if it has a unique state
    ! Uniqueness is determined by comparing IDs and positions

    class(BasinType), intent(inout) :: this
    class(EventType), intent(in) :: event1
    logical, intent(out), optional :: wasAdded

    type(EventType), allocatable :: existingEvents(:)
    integer :: i, numEvents
    logical :: isEventAdded, isDuplicate

    ! Default is that event will be added
    isEventAdded = .true.

    ! Check if there are any existing events
    numEvents = this % absorbingEvents % getCount()

    if (numEvents > 0) then
      ! Get all existing events
      existingEvents = this % absorbingEvents % getAllEvents()

      ! Check each existing event for duplicate state
      isDuplicate = .false.
      do i = 1, numEvents
        ! Compare states using the equals method (compares IDs and positions)
        if (event1 % equals(existingEvents(i))) then
          isDuplicate = .true.
          isEventAdded = .false.
          exit
        end if
      end do

      ! Only add the event if it's not a duplicate
      if (.not. isDuplicate) then
        call this % absorbingEvents % add(event1)
      end if
    else
      ! No existing events, so add this one
      call this % absorbingEvents % add(event1)
    end if

    ! Set output parameter if present
    if (present(wasAdded)) then
      wasAdded = isEventAdded
    end if

  end subroutine BasinType_add_event
  !***********************************************************************!
  subroutine BasinType_flag_event(this, eventId, basinId)
    ! Flag an event as transient and associated with another basin

    class(BasinType), intent(inout) :: this
    integer, intent(in) :: eventId
    integer, intent(in) :: basinId

    type(EventType) :: absorbingEvent
    type(TransientEventType) :: transientEvent

    ! Check if there are events
    if (this % absorbingEvents % getCount() == 0) then
      write (*, '(A)') 'No absorbing events to flag as transient'
      return
    end if

    ! Check if event index is valid
    if (eventId < 1 .or. eventId > this % absorbingEvents % getCount()) then
      write (*, '(A,I0)') 'Invalid event index: ', eventId
      return
    end if

    ! Get the absorbing event to flag
    absorbingEvent = this % absorbingEvents % get(eventId)

    ! Create a new transient event from the absorbing event
    call transientEvent % defineTranEvent(absorbingEvent, basinId)

    ! Add the transient event to the transient array
    call this % transientEvents % add(transientEvent)

    ! Remove the event from the absorbing events array
    call this % absorbingEvents % remove(eventId)

  end subroutine BasinType_flag_event
  !***********************************************************************!
  function BasinType_get_absorbing_count(this) result(count)
    ! Return the number of absorbing events

    class(BasinType), intent(in) :: this
    integer :: count

    count = this % absorbingEvents % getCount()

  end function BasinType_get_absorbing_count
  !***********************************************************************!
  function BasinType_get_transient_count(this) result(count)
    ! Return the number of transient events

    class(BasinType), intent(in) :: this
    integer :: count

    count = this % transientEvents % getCount()

  end function BasinType_get_transient_count
  !***********************************************************************!
  function BasinType_get_absorbing_event(this, index) result(aEvent)
    ! Return a specific absorbing event

    class(BasinType), intent(in) :: this
    integer, intent(in) :: index
    type(EventType) :: aEvent

    aEvent = this % absorbingEvents % get(index)

  end function BasinType_get_absorbing_event
  !***********************************************************************!
  function BasinType_get_transient_event(this, index) result(tEvent)
    ! Return a specific transient event

    class(BasinType), intent(in) :: this
    integer, intent(in) :: index
    type(TransientEventType) :: tEvent

    tEvent = this % transientEvents % get(index)

  end function BasinType_get_transient_event
  !***********************************************************************!
  function BasinType_get_all_absorbing_events(this) result(events)
    ! Return all absorbing events

    class(BasinType), intent(in) :: this
    type(EventType), allocatable :: events(:)

    events = this % absorbingEvents % getAllEvents()

  end function BasinType_get_all_absorbing_events
  !***********************************************************************!
  function BasinType_get_all_transient_events(this) result(events)
    ! Return all transient events

    class(BasinType), intent(in) :: this
    type(TransientEventType), allocatable :: events(:)

    events = this % transientEvents % getAllEvents()

  end function BasinType_get_all_transient_events
  !***********************************************************************!
  subroutine BasinType_copy_from(this, other)
    ! Copy data from another BasinType object

    class(BasinType), intent(inout) :: this
    class(BasinType), intent(in) :: other

    ! Call parent's setState method to copy state data
    call this % setState(other % getIds(), other % getState())

    ! Deep copy the events
    call this % absorbingEvents % copyFrom(other % absorbingEvents)
    call this % transientEvents % copyFrom(other % transientEvents)
  end subroutine
!***********************************************************************!

end module Basin
