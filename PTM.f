Module PolyhedralTemplateMatching

  use LMP_Driver
  use, intrinsic :: iso_c_binding, only: C_NULL_PTR, c_associated

  implicit none
!***********************************************************************!

  Type PTM

    private

    type(C_ptr) :: lammpsInstance = C_NULL_PTR
    character(len=100) :: enabledStructures

    real*8, dimension(:,:), allocatable :: ptmData
    real*8 :: ptmThresh = 0.25

  contains

    procedure :: init => PTM_init
    procedure :: doPTMCalculation

  end type PTM

!***********************************************************************!

  type(PTM) :: ptmCalculator


!***********************************************************************!

contains

!***********************************************************************!
  subroutine PTM_init(this, lammpsContainer, ptmStructures)

    class(PTM), intent(inout) :: this
    class(LAMMPS_Container), intent(in) :: lammpsContainer
    character(len=*), intent(in) :: ptmStructures
    logical :: isError = .false.

    if (.not. c_associated(lammpsContainer%getLammpsPtr())) then
      write(6,'(A)') 'Invalid Lammps Container provided.'
      isError = .true.
    end if
    if (ptmStructures == '') then
      write(6,'(A)') 'No structures were enabled.'
      isError = .true.
    end if

    if (.not. isError) then
      this%lammpsInstance = lammpsContainer%getLammpsPtr()
      this%enabledStructures = ptmStructures
    end if

  end subroutine PTM_init

!***********************************************************************!

  subroutine doPTMCalculation(this)

    class(PTM), intent(inout) :: this
    character(len=200) :: parameters
    real(C_double), pointer, dimension(:,:) :: array => null()
    real(C_double), pointer :: numCols => null()

    integer :: i, natoms
    logical :: isError = .false.

    ! Setup the ptm/array compute
    write(parameters, '(A,X,F5.3,X,A)') trim(this%enabledStructures), this%ptmThresh, 'all'

    call lammps_command(this%lammpsInstance, 'compute 10 all ptm/atom '//trim(parameters))
    call lammps_command(this%lammpsInstance, 'run 0')

    ! Extract the compute data
    call lammps_extract_compute(array, this%lammpsInstance, '10', LMP_STYLE_ATOM, LMP_TYPE_ARRAY)

    if (.not. associated(array)) then 
      write(6,'(A)') 'PTM :: data could not be extracted.'
      isError = .true.
    end if
    
    if (.not. isError) then

      natoms = int(lammps_get_natoms(this%lammpsInstance))
      if (.not. allocated(this%ptmData)) allocate(this%ptmData(7,natoms))

      do i =1 , natoms
        this%ptmData(:,i) = array(:,i)
      end do

    end if

    ! Clean up
    call lammps_command(this%lammpsInstance, 'uncompute 10')

  end subroutine

end module PolyhedralTemplateMatching