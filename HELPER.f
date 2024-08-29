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
    
    n = size(flags); clusterId = 0; dDist = thresh*thresh; dpthresh = 0.999
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
end module HelperFuncs