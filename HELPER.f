module HelperFuncs

contains
!***********************************************************************!
  subroutine clusterAtoms(pos, ptm, flags, thresh, clusters, numClusters)
    !Given positions of points, return an array that gives cluster ID each
    !point belongs to. Clusters will have ID >= 1. IDs will only be valid 
    !in the context they were assigned in. Points are assigned to the same
    !cluster if they are seperated < thresh.

    use Queue
    use Bins
    use SYSTEM, only: RBOX, EPSX, EPSY, EPSZ

    implicit none

    real*8, intent(in), dimension(:,:) :: pos, ptm
    logical, intent(in), dimension(:) :: flags
    real*8, intent(in) :: thresh

    integer, intent(out) :: numClusters
    integer, intent(out), dimension(size(flags)) :: clusters

    real*8 :: PBC

    integer :: n, clusterId
    integer :: i,j,k,m
    real*8 :: rs,xs,ys,zs,dDist,vc,dpthresh
    logical :: success
    type(Queue_i) :: q
    type(BinGrid) :: bg

    !ALGO CAN BE IMPROVED WITH BINS AND FASTER ALGO CHOICE
    call bg%initBinGrid(30,30,30, RBOX*EPSX, RBOX*EPSY, RBOX*EPSZ)
    call bg%assignAtomsToBins(pos, size(flags))

    !Use breadth first search to create clusters
    
    n = size(flags); clusterId = 0; dDist = thresh*thresh; dpthresh = 0.9999
    call q%init(n)

    do i=1, n

      if ((flags(i) .eqv. .true.) .and. (clusters(i) .eq. 0)) then

        call q%push(i, success)

        clusterId = clusterId + 1

        do while(.not. q%isEmpty())

          call q%pop(k, success)

          if (success .eqv. .true.) then

            clusters(k) = clusterId + numClusters !id is offset by numClusters

            call bg%getNearestNeighbors(pos(1,k),pos(2,k),pos(3,k), thresh*1.5)

            do j = 1, bg%singleListCount

              m = bg%singleNeighList(j) !contains neighbor list of most recent 'get' call

              if (k == m .or. flags(m) .eqv. .false. .or. clusters(m) .ne. 0) cycle

              xs = PBC(pos(1,k) - pos(1,m),1)
              ys = PBC(pos(2,k) - pos(2,m),2)
              zs = PBC(pos(3,k) - pos(3,m),3)
              rs = xs*xs + ys*ys + zs*zs
              vc = dot_product(ptm(4:7,k),ptm(4:7,m))

              ! If atoms are close enough, and of the same orientation, and 
              ! are not already 'looked at' then they are part of the same cluster.
              if (rs < dDist .and. vc > dpthresh .and. clusters(m) .eq. 0) then
                call q%push(m, success)
                clusters(m) = -2 ! to avoid duplicates in queue
              end if

            end do

          end if

        end do

      end if

    end do

    numClusters = numClusters + clusterId

  end subroutine clusterAtoms
!***********************************************************************!
  subroutine createPerfectLattice(lattice, center, bounds, ref, count)
    ! Create a reference system that has the lattice type, lattice orientation
    ! center, that fill a box with the given bounds.

    implicit none

    integer, intent(out) :: count
    real*8, intent(in) :: lattice(5), center(3), bounds(6)
    real*8, dimension(:,:), intent(out) :: ref

    integer :: i, j, k, m, ilo, ihi, jlo, jhi, klo, khi
    integer :: maxLatticePoints, unitCellSize, atomCount
    real*8 :: nndist, PBC, rx, ry, rz, rs
    real*8 :: minlims(3), maxlims(3), templims(3)
    real*8 :: primeX(3), primeY(3), primeZ(3), primeR(3)
    real*8 :: offset(3), offsetX(3), offsetY(3), pos(3)
    real*8 :: quatRot(4), invQuatRot(4)
    real*8 :: unitbox(4,3), boxCorners(3,8)
    real*8, allocatable, dimension(:,:) :: primitive
    character(len=5) :: latticeStructure
    real*8, allocatable, dimension(:,:) :: temp

    ! Validate input arguments
    if ( abs(norm2(lattice(2:5)) - 1.0D0) .gt. 5E-4) then
      print *, lattice(:)
      write(*,'(A)') 'Quaternion in lattice argument is not normalized.'
      write(*,'(A)') 'Reference lattice will have dialation issues.'
    end if

    if (bounds(1) .gt. bounds(2) .or. bounds(3) .gt. bounds(4) .or. bounds(5) .gt. bounds(6)) then
      write(*,'()') 'Cannot create reference lattice with invalid bounds.'
      return
    end if

    ! Check if center argument is within bounds

    if (.not. isInBounds(center,bounds)) then
      write(*,'(A)') 'Center is not located within bounds.'
      return
    end if
    
    ! Create primitive units lattice depending on structure type
    select case(int(lattice(1)))
      case(1)
        ! FCC
        latticeStructure = 'fcc'
        unitCellSize = 4
        nndist = (sqrt(2.0D0) / 2)**2
      case (3)
        ! BCC
        latticeStructure = 'bcc'
        unitCellSize = 2
        nndist = (sqrt(3.0D0) / 2)**2
    end select

    allocate(primitive(4,unitCellSize))

    if (trim(latticeStructure) == 'bcc') then
      primitive(:,1) = (/0.0, 0.0, 0.0, 0.0/)
      primitive(:,2) = (/0.0, 0.5, 0.5, 0.5/)
    end if

    if (trim(latticeStructure) == 'fcc') then
      primitive(:,1) = (/0.0, 0.0, 0.0, 0.0/)
      primitive(:,2) = (/0.0, 0.0, 0.5, 0.5/)
      primitive(:,3) = (/0.0, 0.5, 0.0, 0.5/)
      primitive(:,4) = (/0.0, 0.5, 0.5, 0.0/)
    end if

    ! Rotate the primitive unit lattice using quaternion math: P' = H(H(R, P), R')
    quatRot = lattice(2:5)
    invQuatRot = (/quatRot(1),-1*quatRot(2),-1*quatRot(3),-1*quatRot(4)/)

    do i = 1, unitCellSize
      primitive(:,i) = Hamilton(Hamilton(quatRot, primitive(:,i)), invQuatRot)
    end do
    
    ! Rotate unitbox
    unitbox(:,1) = (/0.0, 1.0, 0.0, 0.0/)
    unitbox(:,2) = (/0.0, 0.0, 1.0, 0.0/)
    unitbox(:,3) = (/0.0, 0.0, 0.0, 1.0/)
    do i = 1, 3
      unitbox(:,i) = Hamilton(Hamilton(quatRot, unitbox(:,i)), invQuatRot)
    end do

    primeX = unitbox(2:4,1)
    primeY = unitbox(2:4,2)
    primeZ = unitbox(2:4,3)

    boxCorners(:,1) = (/bounds(1),bounds(3),bounds(5)/)
    boxCorners(:,2) = (/bounds(1),bounds(3),bounds(6)/)
    boxCorners(:,3) = (/bounds(1),bounds(4),bounds(6)/)
    boxCorners(:,4) = (/bounds(2),bounds(4),bounds(6)/)
    boxCorners(:,5) = (/bounds(2),bounds(3),bounds(5)/)
    boxCorners(:,6) = (/bounds(2),bounds(4),bounds(5)/)
    boxCorners(:,7) = (/bounds(1),bounds(4),bounds(5)/)
    boxCorners(:,8) = (/bounds(2),bounds(3),bounds(6)/)

    ! Determine needed number of unitcells in x' and y' directions

    minlims = huge(1)
    maxlims = -huge(1)

    do i = 1,8 ! select corner
      primeR = boxCorners(:,i) - center
      templims = (/dot_product(primeR,primeX), dot_product(primeR,primeY), dot_product(primeR,primeZ)/)
      templims = ceiling(templims)
      do j = 1, 3
        if (templims(j) >= 0) maxlims(j) = max(maxlims(j), templims(j))
        if (templims(j) < 0) minlims(j) = min(minlims(j), templims(j))
      end do
    end do

    ihi = maxlims(1); ilo = minlims(1)
    jhi = maxlims(2); jlo = minlims(2)
    khi = maxlims(3); klo = minlims(3)

    maxLatticePoints = (ihi-ilo+1)*(jhi-jlo+1)*(khi-klo+1)*unitCellSize

    allocate(temp(3,maxLatticePoints))
    atomCount = 0

    do i = ilo, ihi
      offsetX = i*primeX
      do j = jlo, jhi
        offsetY = j*primeY
        do k = klo, khi
          offset = offsetX + offsetY + k*primeZ
          do m = 1, unitCellSize
            pos = center + offset + primitive(2:4, m)
            if (isInBounds(pos,bounds)) then
              atomCount = atomCount + 1
              temp(:,atomCount) = pos
            end if
          end do
        end do
      end do
    end do

    ! call removeOverlapping(temp, atomCount, bounds)

    ! Check if reference array has enough space for all the atoms
    if (atomCount > size(ref,dim=2)) then
      write(*,'(A)') 'Array does not have enough space for reference lattice.'
      return
    end if

    count = atomCount

    ref(:,:atomCount) = temp(:,:atomCount)
    deallocate(temp)
    
  end subroutine createPerfectLattice
!***********************************************************************!
  function Hamilton(q1,q2) result(R)
  ! Reference: https://en.wikipedia.org/wiki/Quaternion#Hamilton_product
  
    implicit none

    real*8, intent(in), dimension(4) :: q1,q2
    real*8, dimension(4) :: R

    real*8 :: a1,b1,c1,d1,a2,b2,c2,d2

    a1 = q1(1); b1=q1(2); c1=q1(3); d1=q1(4)
    a2 = q2(1); b2=q2(2); c2=q2(3); d2=q2(4)

    R(1) = a1*a2 - b1*b2 - c1*c2 - d1*d2 
    R(2) = a1*b2 + b1*a2 + c1*d2 - d1*c2
    R(3) = a1*c2 - b1*d2 + c1*a2 + d1*b2
    R(4) = a1*d2 + b1*c2 - c1*b2 + d1*a2

  end function Hamilton
!***********************************************************************!
  function isInBounds(pos, bounds)

    implicit none

    real*8 :: pos(3), bounds(6)
    logical :: isInBounds

    if (pos(1) > bounds(1) .and. pos(1) < bounds(2) .and. &
       &pos(2) > bounds(3) .and. pos(2) < bounds(4) .and. &
       &pos(3) > bounds(5) .and. pos(3) < bounds(6)) then
      isInBounds = .true.
    else
      isInBounds = .false.
    end if

  end function isInBounds
!***********************************************************************!
  subroutine removeOverlapping(positions, count, bounds)

    use Bins

    implicit none

    integer, intent(inout) :: count
    real*8, dimension(:,:), intent(inout) :: positions
    real*8, dimension(:), intent(in) :: bounds

    integer :: i, j, id, sz
    real*8 :: x,y,z,rx,ry,rz,rs,overlap
    logical, allocatable, dimension(:) :: remove
    type(BinGrid) :: bg

    real*8 :: PBC

    overlap = 0.1**2
    allocate(remove(count)); remove = .false.

    call bg%initBinGrid(40,40,40,bounds(2),bounds(4),bounds(6))
    call bg%assignAtomsToBins(positions, count)

    ! Check for points which overlap
    do i = 1 , count
      if (remove(i)) cycle
      x = positions(1,i)
      y = positions(2,i)
      z = positions(3,i)
      call bg%getNearestNeighbors(x,y,z,2.0D0)
      do j = 1, bg%singleListCount
        id = bg%singleNeighList(j)
        if (id == i) cycle !avoid comparing to self
        rx = PBC(x - positions(1,id),1)
        ry = PBC(y - positions(2,id),2)
        rz = PBC(z - positions(3,id),3)
        rs = rx*rx + ry*ry + rz*rz
        if (rs <= overlap) then
          remove(id) = .true.
        end if
      end do
    end do

    sz = count
    count = 0
    do i = 1, sz
      if (remove(i)) cycle
      count = count + 1
      positions(:,count) = positions(:,i)
    end do

  end subroutine removeOverlapping
!***********************************************************************!
end module HelperFuncs