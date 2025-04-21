module ReferenceLattice

  implicit none
!***********************************************************************!

  Type RefLattice

    integer :: numPoints
    
    real*8 :: bounds(6)
    real*8, allocatable, dimension(:,:) :: lattice

    logical :: initiated = .false.

  contains

    procedure :: createRefLattice
    procedure :: filterReferenceLattice

  end type RefLattice

!***********************************************************************!

contains

!***********************************************************************!
  subroutine createRefLattice(this, grainData)
  !Create the initial reference system that will be used to identify 
  !point defects and their clusters.

    use SYSTEM
    use LMP_DRIVER
    use MPI
    use GrainIdentification

    implicit none

    class(RefLattice), intent(inout) :: this
    type(Grain), intent(in) :: grainData

    integer :: i
    real*8 :: bounds(6)
    logical :: exists

    !Only use atoms which are not apart of grain boundary in the creation
    !of the reference lattice.
    bounds = (/0.0D0, boxx, 0.0D0, boxy, 0.0D0, boxz/)
    call this%filterReferenceLattice(grainData,bounds,RA,natomc)

    NATOMR = this%numPoints
    if (.not. allocated(RR)) allocate(RR(3,NATOMR))
    RR(:,:NATOMR) = this%lattice(:,:NATOMR)

    ! Create the pos_ref.dat reference file in the current directory.
    inquire (file='pos_ref.dat', exist=exists)
    call MPI_BCAST(exists, 1, MPI_logical, 0, lmpSystem%getCommunicator(), ierr)
    if (.not. exists) then
      if (rank == 0) then
        write (*, '(A)') "Creating 'pos_ref.dat' reference file."
      end if
    else
      if (rank == 0) then
        write (*, '(A)') '**************************************************'
        write (*, '(A)') "Please remove file pos_ref.dat from the directory."
        write (*, '(A)') '**************************************************'
      end if
      call MPI_FINALIZE(ierr)
      stop
    end if

    if (rank == 0) then
      open (unit=21, file='pos_ref.dat', status='new')
      write (21, '(I0)') NATOMR
      write (21, *)
      write (21, '(3F14.5)') (RR(1, i), RR(2, i), RR(3, i), i=1, NATOMR)
      close(21)
    end if

  end subroutine createRefLattice
!***********************************************************************!
  subroutine filterReferenceLattice(this, grainData, bounds, pos, count)
  !Create the reference lattice for wigner sietz analysis. Essentially,
  !the ref lattice is made up atoms that have sound structural properties.

    use system
    use GrainIdentification

    implicit none

    class(RefLattice), intent(inout) :: this
    type(Grain), intent(in) :: grainData
    integer, intent(in) :: count
    real*8, intent(in) :: bounds(6)
    real*8, dimension(:,:), intent(in) :: pos ! Atom positions dim(3,N)

    integer :: i, grainId

    !Create the reference lattice
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

  end subroutine filterReferenceLattice
!***********************************************************************!

end module ReferenceLattice