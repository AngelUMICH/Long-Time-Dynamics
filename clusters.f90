Module Clusters
  ! Clusters module will be responsible for storing all the
  ! cluster information needed for any operation in OPLD.
  ! For example, this class will hold whether a defect is
  ! allowed to be chosen for dimer searches.

  implicit none
  !***********************************************************************!

  Type DefectClusters

    private
    
    logical :: initiated = .false.

    real*8, allocatable, dimension(:,:) :: bannedDefectLocs
    integer, allocatable  :: bannedDefectCount

  contains

    procedure :: init
    procedure :: isBanned
    procedure :: banDefect
    procedure :: clearBanList
    procedure :: isInitiated
    procedure :: printBanned

  end type DefectClusters

  !***********************************************************************!

contains

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

      allocate(this%bannedDefectLocs(3,maxNumberDefects))
      this%bannedDefectCount = 0
      this%initiated = .True.

    end if

  end subroutine init
  !***********************************************************************!
  function isBanned(this, x, y, z)
    ! Given coordinates (x,y,z) function returns if defect found near coor.
    ! is not allowed to be chosen for saddle point searching.

    class(DefectClusters), intent(inout) :: this
    real*8, intent(in) :: x,y,z
    logical :: isBanned

    integer :: i
    real*8 :: rx, ry, rz, rs

    if (.not. this%isInitiated()) return

    do i = 1, this%bannedDefectCount

      rx = x - this%bannedDefectLocs(1,i)
      ry = y - this%bannedDefectLocs(2,i)
      rz = z - this%bannedDefectLocs(3,i)
      rs = rx*rx + ry*ry + rz*rz

      rs = sqrt(rs)

      if (rs < 0.01) then
        isBanned = .True.
        exit
      end if

    end do

    isBanned = .False.

  end function isBanned
  !***********************************************************************!
  subroutine banDefect(this, x, y, z)
    ! Add the coordinate (x,y,z) to the banned list of defects.

    class(DefectClusters), intent(inout) :: this
    real*8, intent(in) :: x,y,z

    if (.not. this%isInitiated()) return

    this%bannedDefectCount = this%bannedDefectCount + 1
    this%bannedDefectLocs(:,this%bannedDefectCount) = (/x, y, z/)

  end subroutine banDefect
  !***********************************************************************!
  subroutine clearBanList(this)
    ! Clear the ban list. Used when list is no longer valid such as after
    ! cascade or defect introduction.

    class(DefectClusters), intent(inout) :: this

    if (.not. this%isInitiated()) return

    this%bannedDefectCount = 0
    this%bannedDefectLocs = 0

  end subroutine clearBanList
  !***********************************************************************!
  function isInitiated(this)
    ! Return whether the object has been initiated or not.

    class(DefectClusters), intent(inout) :: this
    logical :: isInitiated

    isInitiated = this%initiated

  end function isInitiated
  !***********************************************************************!
  subroutine printBanned(this)
    ! Will print summary of banned defect locations

    class(DefectClusters), intent(inout) :: this

    write(*,'(A,I0,A)') 'Banned ', this%bannedDefectCount, ' Defect Positions'
    write(*,*)

  end subroutine printBanned
  !***********************************************************************!

end module Clusters