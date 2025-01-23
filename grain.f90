module GrainIdentification

  use PolyhedralTemplateMatching, only: PTM

  implicit none
!***********************************************************************!

  Type Grain

    type(PTM) :: PTMCalculator

    real*8 :: cDistThresh = 1.0 !lattice parameter
    integer :: minGrainSize
    integer:: numGrains = 0
    integer:: maxNumGrains = 50

    character(len=50) :: structuresArg = 'default'

    integer, dimension(:), allocatable :: grainOccupy ! grain atom assignments
    real*8, dimension(:,:), allocatable :: grainInfo ! info on grain structures
                                                     ! type, qw, qx, qy, qz

  contains

    procedure :: init => GRAIN_init
    procedure :: identifyGrains
    procedure :: filterSmallGrains
    procedure :: allocateClusterArray
    procedure :: setGrainIdentity
    procedure :: adoptAtoms
    ! procedure :: identifyBoundaries

  end type Grain

!***********************************************************************!

contains

!***********************************************************************!  

  subroutine GRAIN_init(this, lammpsContainer, minGrainSize)

    use LMP_DRIVER, only: LAMMPS_Container

    implicit none

    class(Grain), intent(inout) :: this
    class(LAMMPS_Container), intent(in) :: lammpsContainer
    integer, intent(in) :: minGrainSize

    this%minGrainSize = minGrainSize

    call this%PTMCalculator%init(lammpsContainer, this%structuresArg)

  end subroutine GRAIN_init

!***********************************************************************!

  subroutine identifyGrains(this)
  ! Cluster the atoms with similar crystal structure to create grains.
  !Subroutine will identify grains with the following strucutre types:
  !1) FCC, 2) HCP, 3) BCC, 4) ICOs

    use SYSTEM, only: RA, NATOMC
    use HelperFuncs

    implicit none

    class(Grain), intent(inout) :: this

    integer :: i
    integer :: nClusters_1, nClusters_2
    logical, dimension(:), allocatable :: flags
    integer, dimension(:), allocatable :: clusterData

    this%numGrains = 0
    if (allocated(this%grainOccupy)) this%grainOccupy = 0
    if (allocated(this%grainInfo)) this%grainInfo = 0

    ! Use polyhedral template matching find all grains and their atoms
    call this%PTMCalculator%doPTMCalculation()

    ! For every structure type find clusters
    allocate(flags(NATOMC))
    allocate(clusterData(NATOMC)); clusterData = 0

    nClusters_1 = 0; nClusters_2 = 0

    do i = 1, 4
      flags = (this%PTMCalculator%ptmData(1,:) == i)
      call clusterAtoms(RA, this%PTMCalculator%ptmData, flags, this%cDistThresh, clusterData, nClusters_1)
    end do

    call this%filterSmallGrains(clusterData, nClusters_1, nClusters_2)

    call this%allocateClusterArray(natomc)

    call this%setGrainIdentity(clusterData, nClusters_2)

    ! call this%adoptAtoms(1.0D0)

  end subroutine identifyGrains
!***********************************************************************!
  subroutine filterSmallGrains(this, clusterIds, numGrains, newNumGrains)
  ! For each cluster/grain identified filter out grains which do not
  ! meet the minimum grain size requirement. Also, modify the grains ids
  ! to have ids start at 1.

    implicit none

    class(Grain), intent(inout) :: this

    integer, dimension(:), intent(inout) :: clusterIds
    integer, intent(in) :: numGrains
    integer, intent(out) :: newNumGrains

    integer :: i, j, n, currentId
    integer :: grainSize(numGrains)

    n = size(clusterIds); grainSize = 0

    do i = 1, n
      if (clusterIds(i) .gt. 0) grainSize(clusterIds(i)) = grainSize(clusterIds(i)) + 1
    end do

    newNumGrains = numGrains; currentId = 1
    do i = 1, numGrains
      if (grainSize(currentId) < this%minGrainSize) then

        do j = 1, n
          if (clusterIds(j) .eq. currentId) clusterIds(j) = 0
          if (clusterIds(j) .gt. currentId) clusterIds(j) = clusterIds(j) - 1
        end do

        do j = currentId+1, numGrains
          grainSize(j-1) = grainSize(j)
        end do

        newNumGrains = newNumGrains - 1
      else
        currentId = currentId + 1
      end if
    end do

  end subroutine filterSmallGrains
!***********************************************************************!
  subroutine allocateClusterArray(this, n)
  ! Create the arrays for containing grain ids and grain information.

    implicit none

    class(Grain), intent(inout) :: this
    integer, intent(in) :: n

    if (.not. allocated(this%grainOccupy)) then
        allocate(this%grainOccupy(n))
        this%grainOccupy = 0
    end if
    if (.not. allocated(this%grainInfo)) then
      allocate(this%grainInfo(5, this%maxNumGrains))
      this%grainInfo = 0
    end if

  end subroutine allocateClusterArray
!***********************************************************************!
  subroutine setGrainIdentity(this, clusterIds, numGrains)
  ! Set the official grain ids for each atom. Also, get an average of
  ! quarternion for each grain. The average will be used to create the
  ! reference system later.

    implicit none

    class(Grain), intent(inout) :: this
    integer, dimension(:), intent(inout) :: clusterIds
    integer, intent(in) :: numGrains

    integer :: count, id
    real*8 :: qSum(5), qAvg(5)

    integer :: i, n

    n = size(clusterIds)

    this%numGrains =  numGrains

    this%grainOccupy = clusterIds

    do id = 1, numGrains
      count = 0; qSum = 0
      do i = 1, n
        if (clusterIds(i) .eq. id) then
          qSum(1) = this%PTMCalculator%ptmData(1,i)
          qSum(2:5) = this%PTMCalculator%ptmData(4:7,i) + qSum(2:5)
          count = count + 1
        end if
      end do

      qAvg(1) = qSum(1)
      qAvg(2:5) = qSum(2:5)/count
      this%grainInfo(:,id) = qAvg

    end do


  end subroutine setGrainIdentity
!***********************************************************************!
  subroutine adoptAtoms(this, thresh)
    ! Atoms not assigned to a grain/cluster will be adopted by the nearest
    ! grain. Grain assignment is done by determining the neighbors and
    ! selecting the most popular grain type. This method is not purely
    ! deterministic. 

    use SYSTEM, only: NATOMC, RA, RBOX, EPSX, EPSY, EPSZ
    use Bins

    implicit none

    class(Grain), intent(inout) :: this
    real*8, intent(in) :: thresh

    real*8 :: PBC

    integer :: i,j,k,m,orphanCount,tempCount
    integer, parameter :: neighLimit = 3
    integer, dimension(:), allocatable :: neighAssign, orphanAtoms
    real*8 :: rs,xs,ys,zs,delta

    type(BinGrid) :: bg

    call bg%initBinGrid(10,10,10, RBOX*EPSX, RBOX*EPSY, RBOX*EPSZ)
    call bg%assignAtomsToBins(RA, natomc)

    delta = thresh*thresh
    allocate(neighAssign(this%numGrains))
    allocate(orphanAtoms(natomc))

    orphanAtoms = 0
    orphanCount = 0

    ! Make list of all orphan atoms
    do i = 1, natomc
      if (this%grainOccupy(i) .eq. 0) then
        orphanCount = orphanCount + 1
        orphanAtoms(orphanCount) = i
      end if
    end do


    ! Adopting all orphan atoms.
    do while (orphanCount > 0)

      tempCount = orphanCount

      do i = 1, orphanCount

          j = orphanAtoms(i)

          neighAssign = 0

          call bg%getNearestNeighbors(RA(1,j),RA(2,j),RA(3,j), thresh)

          do k = 1, bg%singleListCount

            m = bg%singleNeighList(k)

            if (this%grainOccupy(m) .eq. 0) cycle ! do not compare to unassigned

            xs = PBC(RA(1,j) - RA(1,m),1)
            ys = PBC(RA(2,j) - RA(2,m),2)
            zs = PBC(RA(3,j) - RA(3,m),3)
            rs = xs*xs + ys*ys + zs*zs

            ! if (rs < delta) then
              neighAssign(this%grainOccupy(m)) = neighAssign(this%grainOccupy(m)) + 1
            ! end if

          end do

        if (maxval(neighAssign) > 3) then
          this%grainOccupy(j) = maxloc(neighAssign,1)
          tempCount = tempCount - 1
        end if

      end do

      orphanCount = tempCount

    end do

  end subroutine adoptAtoms
!***********************************************************************!
end module GrainIdentification