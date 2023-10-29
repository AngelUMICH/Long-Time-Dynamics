Module Bins

  implicit none

  integer :: tempAtomListSize = 400
  integer :: maxNeighborSize = 20000
!***********************************************************************!
  Type Bin

    integer :: nAtoms = 0
    logical :: isPrinted = .false.
    integer, allocatable, dimension(:) :: tempAtomList
    integer, allocatable, dimension(:) :: atomList

  contains

    procedure :: addAtom
    procedure :: resizeAtomList

  end type Bin
!***********************************************************************!
  Type BinGrid

    integer :: nBins, nBinsX, nBinsY, nBinsZ
    logical :: isInitialized = .false.
    real*8 :: xBinSize, yBinSize, zBinSize !lattice constant
    type(Bin), allocatable, dimension(:, :, :) :: binList

  contains

    procedure :: initBinGrid
    procedure :: assignAtomsToBins
    procedure :: getNearestNeighbors

  end type BinGrid
!***********************************************************************!
contains
!***********************************************************************!
  subroutine addAtom(this, atomID)

    class(Bin), intent(inout) :: this
    integer, intent(in) :: atomID

    if (this%nAtoms .eq. size(this%tempAtomList)) then
      if (.not. this%isPrinted) then
        print *, "Make tempAtomList larger. Current size = ", tempAtomListSize
        this%isPrinted = .true.
        return
      end if
    end if

    this%tempAtomList(this%nAtoms + 1) = atomID
    this%nAtoms = this%nAtoms + 1

  end subroutine addAtom
!***********************************************************************!
  subroutine resizeAtomList(this)

    implicit none

    class(Bin) :: this

    allocate (this%atomList(this%nAtoms))
    this%atomList = this%tempAtomList(1:this%nAtoms)
    deallocate (this%tempAtomList)
  end subroutine resizeAtomList
!***********************************************************************!
  subroutine initBinGrid(this, nBinsX, nBinsY, nBinsZ, xBoxLength, yBoxLength, zBoxLength)

    implicit none

    class(BinGrid), intent(inout) :: this
    integer, intent(in) :: nBinsX, nBinsY, nBinsZ
    real*8, intent(in) :: xBoxLength, yBoxLength, zBoxLength

    integer :: i, j, k

    allocate (this%binList(nBinsX, nBinsY, nBinsZ))
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
          allocate (this%binList(i, j, k)%tempAtomList(tempAtomListSize))
        end do
      end do
    end do

    this%isInitialized = .true.

  end subroutine initBinGrid
!***********************************************************************!
  subroutine assignAtomsToBins(this, atoms, nAtoms)

    implicit none

    class(BinGrid), intent(inout) :: this
    real*8, dimension(:, :) :: atoms !reference atoms from the main program

    integer :: nAtoms
    integer :: atomID, xBinId, yBinId, zBinId

    if (nAtoms/this%nBins .lt. 4) write (*, *) "Warning: Too few atoms per bin"

    do atomID = 1, nAtoms
      xBinId = int(atoms(1, atomID)/this%xBinSize) + 1
      yBinId = int(atoms(2, atomID)/this%yBinSize) + 1
      zBinId = int(atoms(3, atomID)/this%zBinSize) + 1
      call this%binList(xBinId, yBinId, zBinId)%addAtom(atomID)
    end do

    do xBinId = 1, this%nBinsX
      do yBinId = 1, this%nBinsY
        do zBinId = 1, this%nBinsZ
          call this%binList(xBinId, yBinId, zBinId)%resizeAtomList()
        end do
      end do
    end do

  end subroutine assignAtomsToBins
!***********************************************************************!
  function getNearestNeighbors(this, x, y, z) result(neighborAtoms)

    implicit none

    class(BinGrid), intent(inout) :: this
    real*8, intent(in) :: x, y, z
    integer, dimension(maxNeighborSize) :: neighborAtoms

    integer :: i, xBinId, yBinId, zBinId
    integer :: topX, topY, topZ, botX, botY, botZ

    xBinId = int(x/this%xBinSize) + 1
    yBinId = int(y/this%yBinSize) + 1
    zBinId = int(z/this%zBinSize) + 1
    ! Check if the atom is at the edge of the box; put into the last bin if it is
    if (xBinId .gt. this%nBinsX) xBinId = 20
    if (yBinId .gt. this%nBinsY) yBinId = 20
    if (zBinId .gt. this%nBinsZ) zBinId = 20
    topX = xBinId + 1
    topY = yBinId + 1
    topZ = zBinId + 1
    botX = xBinId - 1
    botY = yBinId - 1
    botZ = zBinId - 1

    !Check if the atom is at the edge of the box; apply periodic boundary conditions
    if (topX .gt. this%nBinsX) topX = 1
    if (botX .lt. 1) botX = this%nBinsX
    if (topY .gt. this%nBinsY) topY = 1
    if (botY .lt. 1) botY = this%nBinsY
    if (topZ .gt. this%nBinsZ) topZ = 1
    if (botZ .lt. 1) botZ = this%nBinsZ

    i = 0
    neighborAtoms = 0
    !top z plane
    neighborAtoms(i + 1:i + this%binList(topX, topY, topZ)%nAtoms) = this%binList(topX, topY, topZ)%atomList
    i = i + this%binList(topX, topY, topZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(xBinId, topY, topZ)%nAtoms) = this%binList(xBinId, topY, topZ)%atomList
    i = i + this%binList(xBinId, topY, topZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(botX, topY, topZ)%nAtoms) = this%binList(botX, topY, topZ)%atomList
    i = i + this%binList(botX, topY, topZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(topX, yBinId, topZ)%nAtoms) = this%binList(topX, yBinId, topZ)%atomList
    i = i + this%binList(topX, yBinId, topZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(xBinId, yBinId, topZ)%nAtoms) = this%binList(xBinId, yBinId, topZ)%atomList
    i = i + this%binList(xBinId, yBinId, topZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(botX, yBinId, topZ)%nAtoms) = this%binList(botX, yBinId, topZ)%atomList
    i = i + this%binList(botX, yBinId, topZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(topX, botY, topZ)%nAtoms) = this%binList(topX, botY, topZ)%atomList
    i = i + this%binList(topX, botY, topZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(xBinId, botY, topZ)%nAtoms) = this%binList(xBinId, botY, topZ)%atomList
    i = i + this%binList(xBinId, botY, topZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(botX, botY, topZ)%nAtoms) = this%binList(botX, botY, topZ)%atomList
    i = i + this%binList(botX, botY, topZ)%nAtoms
    !middle z plane
    neighborAtoms(i + 1:i + this%binList(topX, topY, zBinId)%nAtoms) = this%binList(topX, topY, zBinId)%atomList
    i = i + this%binList(topX, topY, zBinId)%nAtoms
    neighborAtoms(i + 1:i + this%binList(xBinId, topY, zBinId)%nAtoms) = this%binList(xBinId, topY, zBinId)%atomList
    i = i + this%binList(xBinId, topY, zBinId)%nAtoms
    neighborAtoms(i + 1:i + this%binList(botX, topY, zBinId)%nAtoms) = this%binList(botX, topY, zBinId)%atomList
    i = i + this%binList(botX, topY, zBinId)%nAtoms
    neighborAtoms(i + 1:i + this%binList(topX, yBinId, zBinId)%nAtoms) = this%binList(topX, yBinId, zBinId)%atomList
    i = i + this%binList(topX, yBinId, zBinId)%nAtoms
    neighborAtoms(i + 1:i + this%binList(xBinId, yBinId, zBinId)%nAtoms) = this%binList(xBinId, yBinId, zBinId)%atomList
    i = i + this%binList(xBinId, yBinId, zBinId)%nAtoms
    neighborAtoms(i + 1:i + this%binList(botX, yBinId, zBinId)%nAtoms) = this%binList(botX, yBinId, zBinId)%atomList
    i = i + this%binList(botX, yBinId, zBinId)%nAtoms
    neighborAtoms(i + 1:i + this%binList(topX, botY, zBinId)%nAtoms) = this%binList(topX, botY, zBinId)%atomList
    i = i + this%binList(topX, botY, zBinId)%nAtoms
    neighborAtoms(i + 1:i + this%binList(xBinId, botY, zBinId)%nAtoms) = this%binList(xBinId, botY, zBinId)%atomList
    i = i + this%binList(xBinId, botY, zBinId)%nAtoms
    neighborAtoms(i + 1:i + this%binList(botX, botY, zBinId)%nAtoms) = this%binList(botX, botY, zBinId)%atomList
    i = i + this%binList(botX, botY, zBinId)%nAtoms
    !bottom z plane
    neighborAtoms(i + 1:i + this%binList(topX, topY, botZ)%nAtoms) = this%binList(topX, topY, botZ)%atomList
    i = i + this%binList(topX, topY, botZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(xBinId, topY, botZ)%nAtoms) = this%binList(xBinId, topY, botZ)%atomList
    i = i + this%binList(xBinId, topY, botZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(botX, topY, botZ)%nAtoms) = this%binList(botX, topY, botZ)%atomList
    i = i + this%binList(botX, topY, botZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(topX, yBinId, botZ)%nAtoms) = this%binList(topX, yBinId, botZ)%atomList
    i = i + this%binList(topX, yBinId, botZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(xBinId, yBinId, botZ)%nAtoms) = this%binList(xBinId, yBinId, botZ)%atomList
    i = i + this%binList(xBinId, yBinId, botZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(botX, yBinId, botZ)%nAtoms) = this%binList(botX, yBinId, botZ)%atomList
    i = i + this%binList(botX, yBinId, botZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(topX, botY, botZ)%nAtoms) = this%binList(topX, botY, botZ)%atomList
    i = i + this%binList(topX, botY, botZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(xBinId, botY, botZ)%nAtoms) = this%binList(xBinId, botY, botZ)%atomList
    i = i + this%binList(xBinId, botY, botZ)%nAtoms
    neighborAtoms(i + 1:i + this%binList(botX, botY, botZ)%nAtoms) = this%binList(botX, botY, botZ)%atomList
    i = i + this%binList(botX, botY, botZ)%nAtoms

  end function getNearestNeighbors
!***********************************************************************!
end module Bins
