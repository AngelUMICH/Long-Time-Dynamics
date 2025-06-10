Module Bins

  implicit none

  integer :: maxListSize = 400
  integer :: maxNeighborSize = 2000
!***********************************************************************!
  Type Bin
    ! Atom ids cannot start with zero otherwise errors will incur
    ! Could add check to avoid atom id equal 0

    integer :: nAtoms = 0
    logical :: isPrinted = .false.
    integer, allocatable, dimension(:) :: atomList

  contains

    procedure :: addAtom
    procedure :: resetBin

  end type Bin
!***********************************************************************!
  Type BinGrid

    integer :: nBins, nBinsX, nBinsY, nBinsZ
    logical :: isInitialized = .false.
    real*8 :: xBinSize, yBinSize, zBinSize !lattice constant
    integer :: singleListCount
    integer, allocatable, dimension(:) :: singleNeighList
    type(Bin), allocatable, dimension(:, :, :) :: binList

  contains

    procedure :: initBinGrid
    procedure :: assignAtomsToBins
    procedure :: getNearestNeighbors
    procedure :: resetGrid

  end type BinGrid
!***********************************************************************!
contains
!***********************************************************************!
  subroutine resetBin(this)

    implicit none

    class(Bin), intent(inout) :: this

    this%atomList = 0
    this%nAtoms = 0
    this%isPrinted = .false.

  end subroutine resetBin
!***********************************************************************!
  subroutine addAtom(this, atomID)

    implicit none

    class(Bin), intent(inout) :: this
    integer, intent(in) :: atomID

    if (this%nAtoms .ge. maxListSize) then
      if (.not. this%isPrinted) then
        print *, "Make atomList larger. Current size = ", maxListSize
        this%isPrinted = .true.
      end if
      return
    end if

    this%atomList(this%nAtoms + 1) = atomID
    this%nAtoms = this%nAtoms + 1

  end subroutine addAtom
!***********************************************************************!
  subroutine initBinGrid(this, nBinsX, nBinsY, nBinsZ, xBoxLength, yBoxLength, zBoxLength)

    implicit none

    class(BinGrid), intent(inout) :: this
    integer, intent(in) :: nBinsX, nBinsY, nBinsZ
    real*8, intent(in) :: xBoxLength, yBoxLength, zBoxLength

    integer :: i, j, k

    if (.not. allocated(this%binList)) allocate(this%binList(nBinsX, nBinsY, nBinsZ))
    this%nBins = nBinsX*nBinsY*nBinsZ
    this%nBinsX = nBinsX
    this%nBinsY = nBinsY
    this%nBinsZ = nBinsZ
    this%xBinSize = xBoxLength/nBinsX
    this%yBinSize = yBoxLength/nBinsY
    this%zBinSize = zBoxLength/nBinsZ

    do i = 1, nBinsX
      do j = 1, nBinsY
        do k = 1, nBinsZ
          if(.not. allocated(this%binList(i, j, k)%atomList)) &
          &   allocate(this%binList(i, j, k)%atomList(maxListSize))
        end do
      end do
    end do

    if (.not. allocated(this%singleNeighList)) allocate(this%singleNeighList(maxNeighborSize))
    this%singleNeighList = 0

  end subroutine initBinGrid
!***********************************************************************!
  subroutine assignAtomsToBins(this, atoms, nAtoms)

    implicit none

    class(BinGrid), intent(inout) :: this
    real*8, dimension(:, :) :: atoms !reference atoms from the main program

    integer :: nAtoms
    integer :: atomID, xBinId, yBinId, zBinId

    call this%resetGrid

    do atomID = 1, nAtoms
      xBinId = int(atoms(1, atomID)/this%xBinSize) + 1
      yBinId = int(atoms(2, atomID)/this%yBinSize) + 1
      zBinId = int(atoms(3, atomID)/this%zBinSize) + 1
      if (xBinId .ge. this%nBinsX) xBinId = this%nBinsX
      if (yBinId .ge. this%nBinsY) yBinId = this%nBinsY
      if (zBinId .ge. this%nBinsZ) zBinId = this%nBinsZ

      if (xBinId <= 0 .or. xBinId > this%nBinsX .or. &
        & yBinId <= 0 .or. yBinId > this%nBinsY .or. &
        & zBinId <= 0 .or. zBinId > this%nBinsZ) then
        print *, 'Bad bin ids:', xBinId, yBinId, zBinId
        return
      end if

      call this%binList(xBinId, yBinId, zBinId)%addAtom(atomID)
    end do

  end subroutine assignAtomsToBins
!***********************************************************************!
  subroutine resetGrid(this)

    implicit none

    class(BinGrid) :: this

    integer :: i, j, k

    do i = 1, this%nBinsX
      do j = 1, this%nBinsY
        do k = 1, this%nBinsZ
          call this%binList(i,j,k)%resetBin
        end do
      end do
    end do

    this%singleListCount = 0
    this%singleNeighList = 0

  end subroutine resetGrid
!***********************************************************************!
  subroutine getNearestNeighbors(this, x, y, z, l)
    ! Given a cube length l and position return an atom list of
    ! which bins make up the cube.

    implicit none

    class(BinGrid), intent(inout) :: this
    real*8, intent(in) :: x, y, z, l

    integer :: i, j, k
    integer :: xBinId, yBinId, zBinId, cxBinId, cyBinId, czBinId
    integer :: count, binCount
    integer, dimension(3) :: corner1, corner2

    ! Check l is not larger than box
    if (2*int(l/this%xBinSize)+1 .ge. this%nBinsX .or. &
    &   2*int(l/this%yBinSize)+1 .ge. this%nBinsY .or. &
    &   2*int(l/this%zBinSize)+1 .ge. this%nBinsZ) then

      write(*,'(A,1f8.3,A)') 'getNearestNeighbors: Radius L of ',l,&
                        & ' is too large. Consider reducing.'
      return

    end if

    ! Find the center bin id in which the position lies
    cxBinId = int(x/this%xBinSize) + 1
    cyBinId = int(y/this%yBinSize) + 1
    czBinId = int(z/this%zBinSize) + 1
    if (cxBinId .gt. this%nBinsX) cxBinId = this%nBinsX
    if (cyBinId .gt. this%nBinsY) cyBinId = this%nBinsY
    if (czBinId .gt. this%nBinsZ) czBinId = this%nBinsX

    !Find the two diagonal bin corners which span the cube
    !with the bin found previously in the center
    corner1(1) = cxBinId - (int(l/this%xBinSize))
    corner1(2) = cyBinId - (int(l/this%yBinSize))
    corner1(3) = czBinId - (int(l/this%zBinSize))

    corner2(1) = cxBinId + (int(l/this%xBinSize))
    corner2(2) = cyBinId + (int(l/this%yBinSize))
    corner2(3) = czBinId + (int(l/this%zBinSize))

    ! Loop through bins and add atom list to result
    count  = 0; binCount = 0

    do i = corner1(1), corner2(1)
      xBinId = i
      if (xBinId .gt. this%nBinsX) xBinId = mod(i, this%nBinsX)
      if (xBinId .le. 0) xBinId = this%nBinsX + i

      do j = corner1(2), corner2(2)
        yBinId = j
        if (yBinId .gt. this%nBinsY) yBinId = mod(j, this%nBinsY)
        if (yBinId .le. 0) yBinId = this%nBinsY + j

        do k = corner1(3), corner2(3)
          zBinId = k
          if (zBinId .gt. this%nBinsZ) zBinId = mod(k, this%nBinsZ)
          if (zBinId .le. 0) zBinId = this%nBinsZ + k

          binCount = this%binList(xBinId,yBinId,zBinId)%nAtoms

          if (abs(count + binCount) .gt. maxNeighborSize) then
            print *, xBinId,yBinId,zBinId, count, binCount
            write(*,'(A,I0,A,I0)') 'Neighbor list size of ',count+binCount,&
                                  &' exceeded max size: ', maxNeighborSize
            write(*,'(A)') 'Consider increasing grid dimensions.'
            return
          end if

          this%singleNeighList(count + 1:count + binCount) &
          & = this%binList(xBinId,yBinId,zBinId)%atomList(1:binCount)
          count = count + binCount

        end do
      end do
    end do

    this%singleListCount = count

  end subroutine getNearestNeighbors
!***********************************************************************!
end module Bins
