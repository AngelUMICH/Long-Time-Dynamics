module TranEvent

  use Event, only: EventType

  implicit none

  private

  public :: TransientEventType

  !*******************************************************************!
  ! Type Declerations
  !*******************************************************************!
  Type, extends(EventType) :: TransientEventType
    ! Transient event types are found to lead to a new basin state.
    ! The basin id is recorded for each transient event to find the
    ! mean time spent in each basin state.

    private
    integer :: basinId = 0 ! the id of basin event leads to

  contains

    procedure :: defineTranEvent => TranEvent_define_tran_event
    procedure :: setBasinId => TranEvent_set_basin_id
    procedure :: getBasinId => TranEvent_get_basin_id
    procedure :: copyFromTranEvent => TranEvent_copy_from

  end Type TransientEventType

contains

  !*******************************************************************!
  ! Methods
  !*******************************************************************!
  subroutine TranEvent_define_tran_event(this, ievent, bId)
    ! Define the transient event by coping the data from an event and
    ! providing the basin id to which the transient event lead to

    class(TransientEventType), intent(inout) :: this
    class(EventType), intent(in) :: ievent
    integer, optional, intent(in) :: bId

    call this % setEnergy(ievent % getEnergy())
    call this%setState(ievent%getIds(), ievent%getState())
    
    if (present(bId)) this%basinId = bId

  end subroutine TranEvent_define_tran_event
  !*******************************************************************!
  subroutine TranEvent_set_basin_id(this, bId)

    class(TransientEventType), intent(inout) :: this
    integer, intent(in) :: bId

    this%basinId = bId

  end subroutine TranEvent_set_basin_id
  !*******************************************************************!
  function TranEvent_get_basin_id(this) result(basin_id)

    class(TransientEventType), intent(in) :: this
    integer :: basin_id

    basin_id = this%basinId

  end function TranEvent_get_basin_id
  !*******************************************************************!
  subroutine TranEvent_copy_from(this, other)
    ! Copy data from another TransientEventType object
    
    class(TransientEventType), intent(inout) :: this
    class(TransientEventType), intent(in) :: other

    ! Copy the state data first
    call this%setState(other%getIds(), other%getState())
    
    ! Copy the EventType-specific attributes
    call this%setEnergy(other%getEnergy())
    
    ! Copy TransientEventType specific attributes
    call this%setBasinId(other%getBasinId())

  end subroutine TranEvent_copy_from
  !*******************************************************************!

end module TranEvent