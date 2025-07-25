module Event

  use State, only: StateType

  implicit none

  private

  public :: EventType

  !*******************************************************************!
  ! Type Declerations
  !*******************************************************************!
  Type, extends(StateType) :: EventType

    private
    real*8 :: energy = 0.0d0

  contains

    procedure :: setEnergy => EventType_set_Energy
    procedure :: getEnergy => EventType_get_Energy
    procedure :: copyFrom => EventType_copy_from

  end Type EventType

contains

  !*******************************************************************!
  ! Methods
  !*******************************************************************!
  subroutine EventType_set_Energy(this, energy)
    ! Define the energy barrier of the event

    class(EventType), intent(inout) :: this
    real*8, intent(in) :: energy

    this % energy = energy

  end subroutine EventType_set_Energy
  !*******************************************************************!
  function EventType_get_Energy(this) result(e)

    class(EventType), intent(in) :: this
    real*8 :: e

    e = this % energy

  end function EventType_get_Energy
  !*******************************************************************!
  ! Free Standing Methods
  !*******************************************************************!
  subroutine EventType_copy_from(this, other)
    ! Copy data from another EventType object
    
    class(EventType), intent(inout) :: this
    class(EventType), intent(in) :: other

    ! Copy the state related data using parent's methods
    call this % setState(other % getIds(), other % getState())

    ! Copy event related data
    call this % setEnergy(other % getEnergy())

  end subroutine EventType_copy_from
  !*******************************************************************!
end module Event
