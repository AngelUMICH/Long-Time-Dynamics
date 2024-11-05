module ReferenceLattice

  use GrainIdentification, only: Grain
  use HelperFuncs, only: createPerfectLattice
  use Bins, only: BinGrid

  implicit none
!***********************************************************************!

  Type RefLattice

    integer :: numPoints
    
    real*8 :: bounds(6)
    real*8, allocatable, dimension(:,:) :: lattice

    logical :: initiated = .false.

  contains

    procedure :: createReferenceLattice
    procedure :: generateReference
    procedure :: getPerfectLattice

  end type RefLattice

!***********************************************************************!

contains


!***********************************************************************!
  subroutine createReferenceLattice(this, grainData, bounds, pos, count)
  ! This will use the initial input system to create the ref lattice.
  ! All that needs to be done initially is to save the atoms without
  ! grain boundary atoms.

    implicit none

    class(RefLattice), intent(inout) :: this
    class(Grain), intent(in) :: grainData
    integer, intent(in) :: count
    real*8, intent(in) :: bounds(6)
    real*8, dimension(:,:), intent(in) :: pos ! Atom positions dim(3,N)

    integer :: i, j, k, grainId

    ! Save the atom positions in the lattice array.

    if (.not. allocated(this%lattice)) allocate(this%lattice(3,count))

    this%numPoints = 0
    this%lattice = 0
    this%bounds = bounds

    do i =1 , count
      grainId = grainData%grainOccupy(i)
      if (.not. grainId == 0) then
        this%numPoints = this%numPoints + 1
        this%lattice(:,this%numPoints) = pos(:, i)
      end if
    end do

  end subroutine createReferenceLattice
!***********************************************************************!
  subroutine generateReference(this, grainData, bounds, pos, atomCount)
  ! Will use the grain data in the Grain object to create reference lat.
  ! For each grain its perf lattice bounded by the grain boundaries
  ! are combined to create a multilattice reference system.

  class(RefLattice), intent(inout) :: this
  class(Grain), intent(in) :: grainData

  integer, intent(in) :: atomCount
  real*8, intent(in) :: bounds(6)
  real*8, dimension(:,:), intent(in) :: pos ! Atom positions dim(3,N)

  integer :: i, j, k, grainId, numGrains, count, neighId
  integer :: maxRefPointCount
  integer, allocatable, dimension(:) :: refCount
  real*8, allocatable, dimension(:,:) :: ref ! Holds one lattice
  real*8, allocatable, dimension(:,:,:) :: references
  logical, allocatable, dimension(:) :: addFlags

  type(BinGrid) :: bg

  real*8 :: rx, ry,rz, rs, PBC

  numGrains = grainData%numGrains

  if (numGrains <= 0) then
    write(*,'(A)') 'Cannot build reference lattice with zero grains'
    return
  end if

  maxRefPointCount = int(1.5*atomCount)

  allocate(refCount(numGrains))
  allocate(ref(3,maxRefPointCount))
  allocate(references(3,maxRefPointCount, numGrains))

  do i = 1, numGrains
    call this%getPerfectLattice(grainData, i, bounds, pos, atomCount, ref, count)
    refCount(i) = count
    references(:,:count,i) = ref(:,:count)
  end do

  if (.not. allocated(this%lattice))allocate(this%lattice(3, 1.5*maxval(refCount,dim=1)))

  ! Only use reference points from each perfect ref lattice that are located
  ! near atoms with correspoding grain id to build multigrain ref lattice.

  allocate(addFlags(maxval(refCount, dim=1)))

  call bg%initBinGrid(40,40,40,bounds(2),bounds(4),bounds(6))

  do i = 1, numGrains
    print *, 'Working on grain #', i
    grainId = i
    count = refCount(grainId)
    ref = references(:,:,grainId)
    call bg%assignAtomsToBins(ref, count)
    
    addFlags = .false.

    do j = 1, atomCount

      if (grainData%grainOccupy(j) .eq. grainId) then

        call bg%getNearestNeighbors(pos(1,j),pos(2,j),pos(3,j),1.0D0)

        do k = 1, bg%singleListCount
          neighId = bg%singleNeighList(k)

          if (.not. addFlags(neighId)) then
            rx = PBC(pos(1,j)-ref(1,neighId),1)
            ry = PBC(pos(2,j)-ref(2,neighId),2)
            rz = PBC(pos(3,j)-ref(3,neighId),3)
            rs = rx*rx + ry*ry + rz*rz

            if (rs <= 0.25)then
              addFlags(neighId) = .true.
              this%numPoints = this%numPoints + 1
              this%lattice(:, this%numPoints) = ref(:, neighId)
            end if

          end if

        end do

      end if

    end do
    print *, 'Reference point count', this%numPoints
  end do

  end subroutine generateReference

!***********************************************************************!

  subroutine getPerfectLattice(this, grainData, grainId, bounds, pos, posCount,ref, refCount)

    implicit none

    class(RefLattice), intent(inout) :: this
    class(Grain), intent(in) :: grainData

    integer, intent(in) :: posCount
    integer, intent(out) :: refCount
    integer, intent(in) :: grainId
    real*8, intent(in) :: bounds(6)
    real*8, dimension(:,:), intent(in) :: pos
    real*8, dimension(:,:), intent(out) :: ref

    integer :: id
    real*8 :: latticeQuat(5), center(3)

    latticeQuat = grainData%grainInfo(:,grainId)
    id = findloc(grainData%grainOccupy, grainId, dim=1)
    call calcCOM(grainData, grainId, pos, posCount, center)
    call createPerfectLattice(latticeQuat, center, bounds, ref, refCount)

  end subroutine getPerfectLattice
!***********************************************************************!
  subroutine calcCOM(grainData, grainId, pos, count, COM)

  implicit none

  class(Grain), intent(in) :: grainData
  integer, intent(in) :: grainId, count
  real*8, dimension(:,:), intent(in) :: pos

  real*8, intent(out) :: COM(3)

  integer :: n, i
  real*8 :: cumCOM(3)

  do i = 1, count
    if (grainData%grainOccupy(i) == grainId) then
      cumCom = cumCOM + pos(:,i)
      n = n + 1
    end if
  end do

  COM = cumCOM/n

  end subroutine calcCOM
!***********************************************************************!
end module ReferenceLattice