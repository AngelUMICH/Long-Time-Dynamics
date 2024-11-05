Module PolyhedralTemplateMatching

  use LMP_Driver, only: LAMMPS_Container
  use, intrinsic :: iso_c_binding, only: C_NULL_PTR, c_associated, C_ptr

  implicit none
!***********************************************************************!

  Type PTM

    type(LAMMPS_Container) :: lmpCon
    type(C_ptr) :: lmpPtr = C_NULL_PTR
    character(len=100) :: enabledStructures

    real*8, dimension(:,:), allocatable :: ptmData
    real*8 :: ptmThresh = 0.10

  contains

    procedure :: init => PTM_init
    procedure :: doPTMCalculation

  end type PTM

!***********************************************************************!

contains

!***********************************************************************!
  subroutine PTM_init(this, lammpsContainer, ptmStructures)

    class(PTM), intent(inout) :: this
    class(LAMMPS_Container), intent(in) :: lammpsContainer
    character(len=*), intent(in) :: ptmStructures
    logical :: isError = .false.

    if (c_associated(this%lmpPtr)) return

    if (.not. c_associated(lammpsContainer%getLammpsPtr())) then
      write(6,'(A)') 'Invalid Lammps Container provided.'
      isError = .true.
    end if
    if (ptmStructures == '') then
      write(6,'(A)') 'No structures were enabled.'
      isError = .true.
    end if

    if (.not. isError) then
      this%lmpCon = lammpsContainer
      this%lmpPtr = lammpsContainer%getLammpsPtr()
      this%enabledStructures = ptmStructures
    end if

  end subroutine PTM_init

!***********************************************************************!

  subroutine doPTMCalculation(this)

    use mpi
    use LAMMPS
    use, intrinsic :: iso_c_binding, only: C_double, C_int

    class(PTM), intent(inout) :: this
    character(len=200) :: parameters
    real(C_double), pointer, dimension(:,:) :: array => null()

    integer :: ierr
    integer :: i, natoms, nlocal
    real*8, dimension(:,:), allocatable :: dataBuff
    integer (C_int), dimension(:), pointer :: atomId

    ! Setup the ptm/array compute
    write(parameters, '(A,X,F5.3,X,A)') trim(this%enabledStructures), this%ptmThresh, 'all'

    call lammps_command(this%lmpPtr, 'compute 10 all ptm/atom '//trim(parameters))
    call lammps_command(this%lmpPtr, 'run 0')

    ! Extract the compute data
    call lammps_extract_compute(array, this%lmpPtr, '10', LMP_STYLE_ATOM, LMP_TYPE_ARRAY)

    ! Copy data from c array to fortran array
    natoms = int(lammps_get_natoms(this%lmpPtr))
    nlocal = lammps_extract_setting(this%lmpCon%getLammpsPtr(), 'nlocal')

    if (.not. allocated(this%ptmData)) allocate(this%ptmData(7,natoms))

    do i =1 , nlocal
      this%ptmData(:,i) = array(:,i)
    end do

    ! Rearrage data for mpi reduce
    call lammps_extract_atom(atomId, this%lmpPtr, 'id')

    allocate(dataBuff, mold = this%ptmData)

    dataBuff = 0
    dataBuff = this%ptmData
    this%ptmData = 0

    do i =1 , nlocal
      this%ptmData(:,atomId(i)) = dataBuff(:,i)
    end do

    dataBuff = 0
    
    call MPI_ALLREDUCE(this%ptmData, dataBuff, size(this%ptmData), MPI_REAL8, MPI_SUM, this%lmpCon%getCommunicator(), IERR)
    this%ptmData = dataBuff

    ! Clean up
    call lammps_command(this%lmpPtr, 'uncompute 10')
    deallocate(dataBuff)

  end subroutine doPTMCalculation
!***********************************************************************!

end module PolyhedralTemplateMatching