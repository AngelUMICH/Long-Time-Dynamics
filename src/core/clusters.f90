Module Clusters
  ! Clusters module will be responsible for storing all the
  ! cluster information needed for any operation in OPLD.
  ! For example, this class will hold whether a defect is
  ! allowed to be chosen for dimer searches.

  !TODO given the correct member visibility attributes for both derived types

  implicit none

  integer, parameter :: MAXCLUSTERSIZE = 1000

  !***********************************************************************!

  Type SingleCluster

    integer, public :: atomCount = 0, refCount = 0, type = 0
    real*8, public :: radiusAV = 0, comx = 0, comy = 0, comz = 0
    integer, dimension(MAXCLUSTERSIZE) :: atomIds, refIds

    logical :: banned = .False.

  contains

    !Add atom or ref id to set
    procedure :: expandCluster
    !Replace old cluster data with new data
    procedure :: replaceData
    !Reset data
    procedure :: resetData
    !Return whether given cluster is the same
    procedure :: isSameCluster
    !Return if the cluster is banned
    procedure :: isBanned
    !Calculate geometric center
    procedure :: calcCenter
    !Determine in which directions defects should be centered
    procedure :: doCentering
    !Actually center defects for ease of calculation
    procedure :: centerDefect

  end type SingleCluster
  !***********************************************************************!

  Type DefectClusters
    
    integer, public :: numClusters, numBanned
    logical :: initiated = .false.

    class(SingleCluster), allocatable :: selectedCluster
    class(SingleCluster), allocatable, dimension(:) :: clusters, bannedClusters

  contains

    !Initiates the DefectClusters class
    procedure :: init
    !Adds the cluster to the banned list
    procedure :: banCluster
    !Checks to see if the class was already initiated
    procedure :: isInitiated
    !Creates the clusters from Wigner Sietz analysis
    procedure :: clusterPointDefects
    !Used to identify the new set of clusters at each iteration
    procedure :: identifyClusters
    !Checks to see if the banned clusters still exist
    procedure :: updateBannedClusters
    !Save cluster selection for saddle point searching
    procedure :: saveClusterSelection
    !Clear list of banned clusters
    procedure :: clearBannedList
    !Print count of banned clusters
    procedure :: printBannedCount
    !Use MPI to sync cluster data across all ranks
    procedure :: syncClusterData

  end type DefectClusters

  !***********************************************************************!

  type(DefectClusters) :: defects

contains

  !***********************************************************************!
  subroutine expandCluster(this, id, type)
    !Will add the species id to the correct data array

    class(SingleCluster), intent(inout) :: this
    integer, intent(in) ::id, type

    if (this%atomCount+1 > MAXCLUSTERSIZE) then
      write(*,'(A)') 'SingleCluster Error: Atom set bounds exceeded'
      return
    else if (this%refCount+1 > MAXCLUSTERSIZE) then
      write(*,'(A)') 'SingleCluster Error: Ref set bounds exceeded'
      return
    end if

    !If type == 1 then it is an interstitial. If == 2 then it is vacancy.
    if (type == 1) then
      this%atomCount = this%atomCount + 1
      this%atomIds(this%atomCount) = id
    else if (type == -1) then
      this%refCount = this%refCount + 1
      this%refIds(this%refCount) = id
    end if

    this%type = this%atomCount - this%refCount

    this%radiusAV = 3.334 + 0.866*abs(this%type)**(1./3.)

  end subroutine expandCluster
  !***********************************************************************!
  subroutine replaceData(this, cluster2)
    !Gets data from cluster #2 to replace its own data

    class(SingleCluster), intent(inout) :: this
    class(SingleCluster), intent(in) :: cluster2

    call this%resetData()

    this%atomCount = cluster2%atomCount
    this%refCount = cluster2%refCount
    this%type = cluster2%type

    this%radiusAV = cluster2%radiusAV
    this%comx = cluster2%comx
    this%comy = cluster2%comy
    this%comz = cluster2%comz

    this%atomIds = cluster2%atomIds
    this%refIds = cluster2%refIds

  end subroutine replaceData
  !***********************************************************************!
  function isSameCluster (this, cluster2)
    !Returns boolean whether clusters are the same. This is done by checking
    !all the constituents atoms and refs ids are the same in both clusters.

    class(SingleCluster), intent(inout) :: this
    class(SingleCluster), intent(in) :: cluster2
    logical :: isSameCluster

    integer :: i
    
    if (this%atomCount /= cluster2%atomCount .or. &
      this%atomCount /= cluster2%atomCount) then
      isSameCluster = .false.
      return
    end if

    isSameCluster = .True.

    do i = 1, this%atomCount
      isSameCluster = any(cluster2%atomIds==this%atomIds(i)) .and. isSameCluster
    end do

    do i = 1, this%refCount
      isSameCluster = any(cluster2%refIds==this%refIds(i)) .and. isSameCluster
    end do

  end function isSameCluster
  !***********************************************************************!
  subroutine resetData(this)
    !Clears all data for cluster

    class(SingleCluster), intent(inout) :: this

    this%atomCount = 0
    this%refCount = 0
    this%type = 0
    this%radiusAV = 0
    this%comx = 0
    this%comy = 0
    this%comz = 0
    this%atomIds = 0
    this%refIds = 0
    this%banned = .False.

  end subroutine resetData
  !***********************************************************************!
  function isBanned(this)

    class(SingleCluster), intent(inout) :: this
    logical :: isBanned

    isBanned = this%banned

  end function isBanned
  !***********************************************************************!
  subroutine calcCenter(this, defects, numDefects)
    !Used after acquiring all members of the cluster. The positions are
    !passed in to calculate the geometric center of the cluster.

    class(SingleCluster), intent(inout) :: this
    real*8, intent(in), dimension(:,:) :: defects
    integer, intent(in) :: numDefects

    integer :: i, index, id, count
    real*8 :: x, y, z
    logical :: iflagx, iflagy, iflagz
    integer, allocatable, dimension(:) :: intIds, vacIds
    real*8, allocatable, dimension(:,:) :: pos
    real*8, allocatable, dimension(:) :: defectIds, defectType
    logical, allocatable, dimension(:) :: mask

    allocate(pos(3,numDefects))
    pos = defects(:,1:numDefects)
    allocate(defectIds(numDefects))
    defectIds = defects(7,1:numDefects)
    allocate(defectType(numDefects))
    defectType = defects(6,1:numDefects)
    allocate(mask(numDefects))

    allocate(intIds(this%atomCount), vacIds(this%refCount))

    !TODO change intids and vacIds to intLocs and vacLocs

    !Check if any cluster member is on the boundary. If so, then move all
    !member to the middle of the box. This corrects for the inconsistencies
    !of calculating the center around periodic boundary conditions.

    iflagx = .False.; iflagy = .False.; iflagz = .False.

    !Checking if defects need to be centered if near boundary

    !Check ints
    mask = defectType == 1.0
    do i = 1, this%atomCount

      id = this%atomIds(i)
      index = findloc(defectIds, id, dim=1, mask=mask)
      intIds(i) = index
      x = pos(1, index); y = pos(2, index); z = pos(3, index)
      call this%doCentering(x, y, z, iflagx, iflagy, iflagz)
    end do

    !Check vacs
    mask = defectType == -1.0
    do i = 1, this%refCount

      id = this%refIds(i)
      index = findloc(defectIds, id, dim=1, mask=mask)
      vacIds(i) = index
      x = pos(1, index); y = pos(2, index); z = pos(3, index)
      call this%doCentering(x, y, z, iflagx, iflagy, iflagz)
    end do

    !Centering defects if needed

    if (any([iflagx, iflagy, iflagz])) then
      !Do int
      do i = 1, this%atomCount

        id = intIds(i)
        x = pos(1, id); y = pos(2, id); z = pos(3, id)
        call this%centerDefect(pos(1,id),pos(2,id),pos(3,id), iflagx, iflagy, iflagz)
      end do

      !Do vacs
      do i = 1, this%refCount

        id = vacIds(i)
        x = pos(1, id); y = pos(2, id); z = pos(3, id)
        call this%centerDefect(pos(1,id),pos(2,id),pos(3,id), iflagx, iflagy, iflagz)
      end do
    end if
    
    ! Ready to calculate geometric center

    count = 0
    !Sum ints
    do i  = 1, this%atomCount

      count = count + 1
      id = intIds(i)
      x = pos(1, id); y = pos(2, id); z = pos(3, id)
      this%comx = this%comx + x
      this%comy = this%comy + y
      this%comz = this%comz + z
    end do

    !Sum vacs
    do i  = 1, this%refCount

      count = count + 1
      id = vacIds(i)
      x = pos(1, id); y = pos(2, id); z = pos(3, id)
      this%comx = this%comx + x
      this%comy = this%comy + y
      this%comz = this%comz + z
    end do

    this%comx = this%comx / count
    this%comy = this%comy / count
    this%comz = this%comz / count

    !Do centering algorithm to revert positions and geometric center
    !if needed.

    if (any([iflagx, iflagy, iflagz])) then
      !Revert geometric center
      call this%centerDefect(this%comx, this%comy, this%comz, iflagx, iflagy, iflagz)

      !Revert ints
      do i = 1, this%atomCount

        id = intIds(i)
        x = pos(1, id); y = pos(2, id); z = pos(3, id)
        call this%centerDefect(x, y, z, iflagx, iflagy, iflagz)
      end do

      !Revert vacs
      do i = 1, this%refCount

        id = vacIds(i)
        x = pos(1, id); y = pos(2, id); z = pos(3, id)
        call this%centerDefect(x, y, z, iflagx, iflagy, iflagz)
      end do
    end if


  end subroutine calcCenter
  !***********************************************************************!
  subroutine doCentering(this, x, y, z, flagx, flagy, flagz)
    !Checks whether boundary conditions should be respected to calculate
    !geometric center for the 3 coordinates

    use system, only: boxx, boxy, boxz

    class(SingleCluster), intent(inout) :: this
    real*8, intent(in) :: x, y, z
    logical, intent(inout) :: flagx, flagy, flagz


    if (x < 1.0 .or. x > boxx - 1.0) then
      flagx = .True.
    end if

    if (y < 1.0 .or. y > boxy - 1.0) then
      flagy = .True.
    end if

    if (z < 1.0 .or. z > boxz - 1.0) then
      flagz = .True.
    end if

  end subroutine doCentering
  !***********************************************************************!
  subroutine centerDefect(this, x, y, z, flagx, flagy, flagz)
    !Adds or subtracts half the box distance to place the defect away from
    !box boundaries.

    use system, only: boxx, boxy, boxz

    class(SingleCluster), intent(inout) :: this
    real*8, intent(inout) :: x, y, z
    logical, intent(in) :: flagx, flagy, flagz

    real*8 :: half

    if (flagx) then
      half = boxx/2
      if (x < half) then
        x = x + half
      else
        x = x - half
      end if
    end if

    if (flagy) then
      half = boxy/2
      if (y < half) then
        y = y + half
      else
        y = y - half
      end if
    end if

    if (flagz) then
      half = boxz/2
      if (z < half) then
        z = z + half
      else
        z = z - half
      end if
    end if

  end subroutine centerDefect
  !***********************************************************************!
  !***********************************************************************!
  subroutine init(this, maxNumberDefects)
    ! Must be used before any function or subroutine of the derived class

    class(DefectClusters), intent(inout) :: this
    integer, intent(in) :: maxNumberDefects

    if (maxNumberDefects < 0) then
      write(*,'(A)') "Argument must be larger than 0."
      return
    end if

    if (.not. this%isInitiated()) then

      this%numClusters = 0
      this%numBanned = 0
      allocate(this%selectedCluster)
      allocate(this%clusters(maxNumberDefects))
      allocate(this%bannedClusters(maxNumberDefects))

      this%initiated = .True.

    end if

  end subroutine init
  !***********************************************************************!
  subroutine banCluster(this, cluster)
    ! Put the given cluster on the banned list

    class(DefectClusters), intent(inout) :: this
    type(SingleCluster), intent(in) :: cluster

    if (.not. this%isInitiated()) return

    this%numBanned = this%numBanned + 1
    call this%bannedClusters(this%numBanned)%replaceData(cluster)

  end subroutine banCluster
  !***********************************************************************!
  function isInitiated(this)
    ! Return whether the object has been initiated or not.

    class(DefectClusters), intent(inout) :: this
    logical :: isInitiated

    isInitiated = this%initiated

  end function isInitiated
  !***********************************************************************!
    subroutine clusterPointDefects(this, defects, numDefects, thresh, clusters, count)
    !Will identify defects by clustering point defects. The threshold
    !is used to limit which point defects are allowed in the same cluster.
    !Clusters and their constituents will be stored in the class variables.

    use Queue
    use Bins
    use SYSTEM, only: boxx, boxy, boxz

    implicit none

    class(DefectClusters), intent(inout) :: this
    real*8, intent(in), dimension(:,:) :: defects
    real*8, intent(in) :: thresh
    integer, intent(in) :: numDefects
    integer, intent(out) :: count
    type(SingleCluster), intent(inout), dimension(:) :: clusters

    real*8 :: PBC
    integer :: n, clusterId, id, type
    integer :: i,j,k,m
    real*8 :: rs,xs,ys,zs,delta,x,y,z
    logical :: success
    real*8, dimension(3,numDefects) :: pos
    real*8, dimension(numDefects) :: defectType, defectId
    logical, dimension(numDefects) :: checked

    type(Queue_i) :: q
    type(BinGrid) :: bg

    pos = defects(1:3,1:numDefects)
    defectType = defects(6,1:numDefects)
    defectId = defects(7,1:numDefects)

    call bg%initBinGrid(70,70,70, boxx, boxy, boxz)
    call bg%assignAtomsToBins(pos, numDefects)

    !Use breadth first search to create clusters
    
    n = numDefects; clusterId = 0; delta = thresh*thresh; checked = .False.
    call q%init(n)

    do i=1, n

      if (.not. checked(i)) then

        call q%push(i, success)

        clusterId = clusterId + 1

        do while(.not. q%isEmpty())

          call q%pop(k, success)

          if (success .eqv. .true.) then

            x = pos(1,k); y = pos(2,k); z = pos(3,k)
            id = int(defectId(k)); type = int(defectType(k))
            checked(k) = .True.

            call clusters(clusterId)%expandCluster(id, type)

            call bg%getNearestNeighbors(x, y, z, thresh*1.5)!remove 1.5

            do j = 1, bg%singleListCount

              m = bg%singleNeighList(j)

              if (k == m .or. checked(m)) cycle

              xs = PBC(x - pos(1,m),1)
              ys = PBC(y - pos(2,m),2)
              zs = PBC(z - pos(3,m),3)
              rs = xs*xs + ys*ys + zs*zs

              ! If atoms are within delta distance they are clustered together
              if (rs < delta .and. .not. checked(m)) then
                call q%push(m, success)
                checked(m) = .True.
              end if

            end do

          end if

        end do

        !After finding every member of cluster Calculate center of mass (COM)
        ! print *, 'Calculating center for id = ', clusterId
        call clusters(clusterId)%calcCenter(defects, numDefects)
        ! print *, 'Finished calculating center'

        !If cluster is not a true defect cluster then remove from list
        if (clusters(clusterId)%type == 0) then
          call clusters(clusterId)%resetData()
          clusterId = clusterId - 1
        end if

      end if

    end do

    count = clusterId

  end subroutine clusterPointDefects
  !***********************************************************************!
  subroutine identifyClusters(this, pointDefects, numPointDefects, clusteringThresh)
    !Used to identify the new set clusters at each iteration. 

    class(DefectClusters), intent(inout) :: this
    integer, intent(in) :: numPointDefects
    real*8, intent(in) :: clusteringThresh
    real*8, dimension(:,:), intent(in) :: pointDefects

    integer :: i, maxNumClusters, newClusterCount
    type(SingleCluster), allocatable, dimension(:) :: newClusters

    maxNumClusters = size(this%clusters)
    allocate(newClusters(maxNumClusters))

    !Cluster point defects to find new clusters in the system
    call this%clusterPointDefects(pointDefects, numPointDefects, &
    clusteringThresh, newClusters, newClusterCount)

    !Replace cluster data with new data
    do i = 1, newClusterCount
      call this%clusters(i)%replaceData(newClusters(i))
    end do

    !Reset any data not in use anymore
    do i = newClusterCount + 1, this%numClusters 
      call this%clusters(i)%resetData()
    end do

    this%numClusters = newClusterCount

    call this%updateBannedClusters()

  end subroutine identifyClusters
  !***********************************************************************!
  subroutine updateBannedClusters(this)
    !Update the list of banned clusters. Check if previous banned clusters
    !still exist and process new ones.

    class(DefectClusters), intent(inout) :: this

    integer :: i, j, newBannedCount
    logical :: foundMatch
    type(SingleCluster) :: bannedCluster, cluster
    type(SingleCluster), dimension(this%numBanned) :: copyBannedClusters

    if (this%numBanned <= 0) return

    !Compare the banned clusters with the new identified clusters. If a
    !match is made make sure to save the that cluster. Otherwise, if no
    !match is made for a previously banned cluster then it will be removed
    !from the banned list.

    copyBannedClusters = this%bannedClusters(1:this%numBanned)
    newBannedCount = 0

    do i = 1, this%numBanned
      bannedCluster = copyBannedClusters(i)
      foundMatch = .False.

      do j = 1, this%numClusters
        cluster = this%clusters(j)
        foundMatch = cluster%isSameCluster(bannedCluster)

        if (foundMatch) then
          newBannedCount = newBannedCount + 1
          call this%bannedClusters(newBannedCount)%replaceData(cluster)
          this%clusters(j)%banned = .True.
          exit
        end if
      end do
    end do

    !Clear data not in use anymore
    if (newBannedCount < this%numBanned) then
      do i = newBannedCount + 1, this%numBanned
        call this%bannedClusters(i)%resetData()
      end do
    end if

    this%numBanned = newBannedCount

  end subroutine updateBannedClusters
  !***********************************************************************!
  subroutine saveClusterSelection(this, cluster)
    !Save the data of the given cluster to be used for saddle point searches

    class(DefectClusters), intent(inout) :: this
    type(SingleCluster), intent(in) :: cluster

    call this%selectedCluster%replaceData(cluster)
    
  end subroutine saveClusterSelection
  !***********************************************************************!
  subroutine clearBannedList(this)
    !Clear the list of banned clusters completely. This should be done when
    !system-wide changes occurs such as irradiation. The data of all clusters
    !are reset.

    class(DefectClusters), intent(inout) :: this

    integer :: i

    do i = 1, this%numBanned
      call this%bannedClusters(i)%resetData()
    end do

    this%numBanned = 0

  end subroutine clearBannedList
  !***********************************************************************!
  subroutine printBannedCount(this)
    !Print the number of banned clusters in the system

    class(DefectClusters), intent(inout) :: this

    write(*, '(a,i0)') 'Banned Cluster Count: ', this%numBanned
    write(*,*)

  end subroutine printBannedCount
  !***********************************************************************!
  subroutine syncClusterData(this, comm)
    !Sync the data of the selected cluster across all ranks. This is because
    !selection is only done by one rank. Not necessary to sync atom ids or
    !ref ids since this info is not needed for dimer searches.

    use mpi

    class(DefectClusters), intent(inout) :: this
    integer, intent(in) :: comm

    integer :: ierr

    call MPI_BCAST(this%selectedCluster%comx, 1, MPI_REAL8, 0, comm, ierr)
    call MPI_BCAST(this%selectedCluster%comy, 1, MPI_REAL8, 0, comm, ierr)
    call MPI_BCAST(this%selectedCluster%comz, 1, MPI_REAL8, 0, comm, ierr)
    call MPI_BCAST(this%selectedCluster%radiusAV, 1, MPI_REAL8, 0, comm, ierr)
    call MPI_BCAST(this%selectedCluster%type, 1, MPI_INT, 0, comm, ierr)
    call MPI_BCAST(this%selectedCluster%atomCount, 1, MPI_INT, 0, comm, ierr)
    call MPI_BCAST(this%selectedCluster%refCount, 1, MPI_INT, 0, comm, ierr)
    
  end subroutine syncClusterData
  !***********************************************************************!

end module Clusters