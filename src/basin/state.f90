module State

  implicit none

  private

  real*8, parameter :: tol = 1.0d-2

  public :: StateType, isEqual

  !*******************************************************************!
  ! Type Decleration
  !*******************************************************************!
  Type StateType

    private
    integer :: nAtoms
    integer, allocatable, dimension(:) :: ids
    real*8, allocatable, dimension(:, :) :: state

  contains

    procedure :: getStateSize => StateType_get_state_size
    procedure :: getIds => StateType_get_ids
    procedure :: getState => StateType_get_state
    procedure :: setState => StateType_set_state
    procedure :: equals => StateType_equals
    procedure :: copyState => StateType_copy_state

  end Type StateType

contains

  !*******************************************************************!
  ! Type Methods
  !***********************************************************************!
  function StateType_get_state_size(this) result(n)
    ! Return the size of ids and length of state
    class(StateType), intent(in) :: this
    integer :: n

    n = this%nAtoms
  end function StateType_get_state_size
  !***********************************************************************!
  function StateType_get_ids(this) result(ids_copy)
    ! Return a copy of the IDs array
    class(StateType), intent(in) :: this
    integer, allocatable, dimension(:) :: ids_copy

    if (allocated(this % ids)) then
      allocate (ids_copy(size(this % ids)))
      ids_copy = this % ids
    end if
  end function StateType_get_ids
  !***********************************************************************!
  function StateType_get_state(this) result(state_copy)
    ! Return a copy of the state array
    class(StateType), intent(in) :: this
    real*8, allocatable, dimension(:, :) :: state_copy

    if (allocated(this % state)) then
      allocate (state_copy(size(this % state, 1), size(this % state, 2)))
      state_copy = this % state
    end if
  end function StateType_get_state
  !***********************************************************************!
  subroutine StateType_set_state(this, ids, positions)
    ! Set the state and ids arrays
    class(StateType), intent(inout) :: this
    integer, dimension(:), intent(in) :: ids
    real*8, dimension(:, :), intent(in) :: positions

    ! Deallocate if already allocated
    if (allocated(this % ids)) deallocate (this % ids)
    if (allocated(this % state)) deallocate (this % state)

    ! Allocate and copy data
    allocate (this % ids(size(ids)))
    allocate (this % state(size(positions, 1), size(positions, 2)))

    this % nAtoms = size(ids)
    this % ids = ids
    this % state = positions
  end subroutine StateType_set_state
  !***********************************************************************!
  function StateType_equals(this, other) result(is_equal)
    ! Compare two StateType instances
    class(StateType), intent(in) :: this
    class(StateType), intent(in) :: other
    logical :: is_equal

    ! Check if both have allocated arrays
    if (.not. allocated(this % ids) .or. .not. allocated(other % ids) .or. &
        .not. allocated(this % state) .or. .not. allocated(other % state)) then
      is_equal = .false.
      return
    end if

    ! Use the existing isEqual function to compare
    is_equal = isEqual(this % state, this % ids, this % nAtoms, &
                       other % state, other % ids, other % nAtoms)
  end function StateType_equals
  !***********************************************************************!
  subroutine StateType_copy_state(this, istate)

    class(StateType), intent(inout) :: this
    class(StateType), intent(in) :: istate

    call this % setState (istate%getIds(), istate%getState())

  end subroutine StateType_copy_state
  !***********************************************************************!
  ! Free Standing Functions
  !***********************************************************************!
  function isEqual(pos1, ids1, n, pos2, ids2, m)
    ! Checks if the point configurations of the two sets are the same
    ! up to a tolerance. Since the two basins might have different sets
    ! only those that are common in both sets are compared.

    real*8, dimension(3, n), intent(in) :: pos1
    integer, dimension(n), intent(in) :: ids1
    integer, intent(in) :: n
    real*8, dimension(3, m), intent(in) :: pos2
    integer, dimension(m), intent(in) :: ids2
    integer, intent(in) :: m

    logical :: isEqual

    integer :: i, j, count
    real*8, allocatable, dimension(:, :) :: A, B

    isEqual = .false.

    ! Check which atoms are common in both basins
    count = 0
    do i = 1, n
      do j = 1, m
        if (ids1(i) == ids2(j)) then
          count = count + 1
        end if
      end do
    end do

    ! If no common atoms, then basins are not equal
    if (count == 0) then
      isEqual = .false.
      return
    end if

    allocate (A(3, count))
    allocate (B(3, count))

    count = 0
    do i = 1, n
      do j = 1, m
        if (ids1(i) == ids2(j)) then
          count = count + 1
          A(:, count) = pos1(:, i)
          B(:, count) = pos2(:, j)
        end if
      end do
    end do

    ! Check if the common atoms are equal
    isEqual = all(abs(A - B) < tol)

  end function isEqual
  !***********************************************************************!

end module State
