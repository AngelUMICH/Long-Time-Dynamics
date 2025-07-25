module BasinNetwork

  use Basin

  implicit none

  private

  public :: BasinNetworkType

  !***********************************************************************!

  Type BasinNetworkType
    ! A network of basins used to implement Mean-Rate-Method.

    private
    integer :: nBasins = 0
    integer :: activeEventId = 0, activeBasinId = 0

    type(BasinType), allocatable :: basins(:)

  contains

    procedure :: addBasin => BasinNetwork_add_basin
    procedure :: isUniqueBasin => BasinNetwork_is_unique_basin
    procedure :: flagBasinEvent => BasinNetwork_flag_basin_event
    procedure :: setActiveIds => BasinNetwork_set_active_ids
    procedure :: getActiveBasinId => BasinNetwork_get_basin_id
    procedure :: getActiveEventId => BasinNetwork_get_event_id

    ! Accessor methods
    procedure :: getBasinCount => BasinNetwork_get_basin_count
    procedure :: getBasin => BasinNetwork_get_basin
    procedure :: getAllBasins => BasinNetwork_get_all_basins

  end type BasinNetworkType

  !***********************************************************************!

contains

  !***********************************************************************!
  subroutine BasinNetwork_add_basin(this, basin1)
    ! Add a basin to the network
    ! Note: This method does not check if the basin is unique
    ! You should call isUniqueBasin first if uniqueness is required

    class(BasinNetworkType), intent(inout) :: this
    class(BasinType), intent(in) :: basin1

    type(BasinType), allocatable :: temp(:)
    integer :: newCapacity, i

    ! Check if basin is valid to add
    if (basin1%getStateSize() == 0) then
      write(*,'(A)') 'Cannot add basin without state to network.'
      return
    end if

    ! Initialize if not already done
    if (.not. allocated(this % basins)) then
      allocate (this % basins(8))  ! Start with capacity of 8
      this % nBasins = 0
    end if

    ! Check if we need to grow the array
    if (this % nBasins >= size(this % basins)) then
      ! Double the array size (common growth strategy)
      newCapacity = max(8, size(this % basins) * 2)

      ! Save current data
      allocate (temp(this % nBasins))
      ! Copy each basin individually using the copyFrom method
      do i = 1, this % nBasins
        call temp(i) % copyFrom(this % basins(i))
      end do

      ! Reallocate with new size
      deallocate (this % basins)
      allocate (this % basins(newCapacity))

      ! Copy data back using copyFrom for each basin
      do i = 1, this % nBasins
        call this % basins(i) % copyFrom(temp(i))
      end do

      ! Clean up
      deallocate (temp)
    end if

    ! Add the new basin
    this % nBasins = this % nBasins + 1
    call this % basins(this % nBasins) % copyFrom(basin1)

  end subroutine BasinNetwork_add_basin
  !***********************************************************************!
  subroutine BasinNetwork_is_unique_basin(this, basin1, isUnique, basinId)
    ! Given a basin1 it returns whether it is unique state and the
    ! corresponding matching basin id.

    class(BasinNetworkType), intent(in) :: this
    class(BasinType), intent(in) :: basin1
    integer, intent(out) :: basinId
    logical :: isUnique

    integer :: i
    type(BasinType) :: currentBasin

    ! Default: assume it's unique
    isUnique = .true.
    basinId = 0

    ! Check against all existing basins
    do i = 1, this % getBasinCount()
      currentBasin = this % getBasin(i)
      if (currentBasin % equals(basin1)) then
        ! Found a match - not unique
        isUnique = .false.
        basinId = i
        return
      end if
    end do

  end subroutine BasinNetwork_is_unique_basin
  !***********************************************************************!
  subroutine BasinNetwork_flag_basin_event(this, basinId, eventId, id)
    ! Flag event eventId in basin basinId.

    class(BasinNetworkType), intent(inout) :: this
    integer, intent(in) :: basinId, eventId, id

    if (basinId > 0 .and. basinId <= this % nBasins) then
      call this % basins(basinId) % flagEvent(eventId, id)
    else
      write(*, '(A,I0)') "Error: Can't flag invalid basin index: ", basinId
    end if

  end subroutine BasinNetwork_flag_basin_event
  !***********************************************************************!
  function BasinNetwork_get_basin_count(this) result(count)
    ! Return the number of basins in the network

    class(BasinNetworkType), intent(in) :: this
    integer :: count

    count = this % nBasins

  end function BasinNetwork_get_basin_count
  !***********************************************************************!
  function BasinNetwork_get_basin(this, index) result(aBasin)
    ! Return a specific basin from the network

    class(BasinNetworkType), intent(in) :: this
    integer, intent(in) :: index
    type(BasinType) :: aBasin

    ! Check if index is valid
    if (index < 1 .or. index > this % nBasins) then
      write (*, '(A,I0)') "Error: Can't get invalid basin index: ", index
      return
    end if

    call aBasin % copyFrom(this % basins(index))

  end function BasinNetwork_get_basin
  !***********************************************************************!
  function BasinNetwork_get_all_basins(this) result(basins)
    ! Return all basins in the network

    class(BasinNetworkType), intent(in) :: this
    type(BasinType), allocatable :: basins(:)

    integer :: i

    if (this % nBasins > 0) then
      allocate (basins(this % nBasins))
      ! Copy each basin individually using copyFrom
      do i = 1, this % nBasins
        call basins(i) % copyFrom(this % basins(i))
      end do
    else
      ! Return empty array
      allocate (basins(0))
    end if

  end function BasinNetwork_get_all_basins
  !***********************************************************************!
  subroutine BasinNetwork_set_active_ids(this, basinId, eventId)
    ! Set the ids for the chosen basin and event at the KMC stage.

    class(BasinNetworkType), intent(inout) :: this
    integer, intent(in) :: basinId, eventId

    this % activeBasinId = basinId
    this % activeEventId = eventId

  end subroutine BasinNetwork_set_active_ids
  !***********************************************************************!
  function BasinNetwork_get_basin_id(this) result(id)
    ! Return the id of the previously chosen basin
    class(BasinNetworkType), intent(in) :: this
    integer :: id

    id = this % activeBasinId

  end function BasinNetwork_get_basin_id
  !***********************************************************************!
  function BasinNetwork_get_event_id(this) result(id)
    ! Return the id of the previously chosen event
    class(BasinNetworkType), intent(in) :: this
    integer :: id

    id = this % activeEventId

  end function BasinNetwork_get_event_id
  !***********************************************************************!
end module BasinNetwork
