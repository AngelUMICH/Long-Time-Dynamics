Module LMP_Driver

  use LAMMPS
  use, intrinsic :: iso_c_binding, only: C_NULL_PTR

  implicit none
!***********************************************************************!

  Type LAMMPS_Container

    private

    type(C_ptr) :: lammpsInstance = C_NULL_PTR

    integer :: communicator
    integer :: controlCommunicator
    integer :: rank
    real(C_double), pointer :: energy => null()
    real(C_double), pointer, dimension(:) :: pressure => null()
    logical :: isDimer = .false.
    logical :: isSystem = .false.

  contains

    procedure :: init => LAMMPS_Container_init
    procedure :: setup
    procedure :: readSystem
    procedure :: readForce
    procedure :: readRestartSystem
    procedure :: getBoxDims
    procedure :: getLammpsPtr
    procedure :: getCommunicator
    procedure :: getControlComm
    procedure :: getRank
    procedure :: getPositions
    procedure :: assignPositions
    procedure :: relaxSystemEnergy
    procedure :: relaxSystemPressure
    procedure :: changeBoxDimensions
    procedure :: writeData
    procedure :: readData

    end type LAMMPS_Container

!***********************************************************************!

  type(LAMMPS_Container) :: lmpSystem
  type(LAMMPS_Container) :: lmpDimer

!***********************************************************************!

contains
!***********************************************************************!
  
  subroutine LAMMPS_Container_init(this, commToSplit, groupId, compType)

    class(LAMMPS_Container), intent(inout) :: this
    integer, intent(in) :: commToSplit, groupId
    character, intent(in) :: compType
    integer :: origRank, ierr

    call MPI_Comm_rank(commToSplit, origRank, ierr)

    call MPI_Comm_split(commToSplit, groupId, origRank, this%communicator, ierr)

    call MPI_Comm_rank(this%communicator, this%rank, ierr)

    if (this%isDimer) then
      call mpi_comm_split(this%communicator, this%rank, this%rank, this%controlCommunicator, ierr)
    end if

    call lammps_open('lmp -log none -screen none', this%communicator, this%lammpsInstance)

    if (compType == 's') then
      this%isSystem = .true.
    else if (compType == 'd') then
      this%isDimer = .true.
    end if

    !Necessary commands for DRLD to work
    call lammps_command(this%lammpsInstance, 'units metal')
    call lammps_command(this%lammpsInstance, 'boundary p p p')
    call lammps_command(this%lammpsInstance, 'atom_style atomic')
    call lammps_command(this%lammpsInstance, 'atom_modify map array')

  end subroutine LAMMPS_Container_init
!***********************************************************************!

  subroutine setup(this, normalFile, restartFile, mode)
    !Setup the initial system (i.e atoms pos and box size) from a LAMMPS
    !data file. In both modes, normal or restart, the atoms and box size
    !are read from a data file either from the 'in.lammps' or restart file.

    implicit none

    class(LAMMPS_Container), intent(inout) :: this
    integer, intent(in) :: mode
    character(len=*), intent(in) :: normalFile, restartFile

    !If mode == 1 use normal file, otherwise, mode == 1 use restart file.
    !Same forces must be used. The force is described in the normal file.
    if (mode == 1) then
      call this%readSystem(normalFile)
      call this%readForce(normalFile)
    else if (mode == 2) then
      call this%readRestartSystem(restartFile)
      call this%readForce(normalFile)
    end if

    !Initiating system values
    call lammps_command(this%lammpsInstance, 'run 0')

    !Setting up energy and pressure measurements
    call lammps_extract_compute(this%energy, this%lammpsInstance, &
                              'thermo_pe', LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)
    call lammps_extract_compute(this%pressure, this%lammpsInstance, &
                              'thermo_press', LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR)

    call this%getBoxDims()

    !Dimer LAMMPS instace should be empty because it will be filled with 
    !atoms during saddle point searches.
    if (this%isDimer) call lammps_command(this%lammpsInstance, 'delete_atoms group all')

  end subroutine setup
!***********************************************************************!
  subroutine readSystem(this, inputFile)
    !Read the system (atoms and box) from the input file.
    !The inputFile obtains information of the system in the first block.
    !A block a collection of LAMMPS commands that are preceeded by the
    !size of the block.

    use mpi

    class(LAMMPS_Container), intent(inout) :: this
    character(len=*), intent(in) :: inputFile

    integer :: blockSize, i, ierr
    logical :: exists
    character(len=200) command

    inquire(file=inputFile, exist=exists)
    if (.not. exists) then
      write(*,'(A)') 'File '//trim(inputFile)//' cannot be found.'
      call MPI_FINALIZE(ierr)
      stop
    end if

    open(unit=21, file=inputFile)

    read(21, *)blockSize
    do i = 1, blockSize
      read(21, '(A)') command
      call lammps_command(this%lammpsInstance, trim(command))
    end do

    close(21)

  end subroutine readSystem
!***********************************************************************!
  subroutine readForce(this, inputFile)
    !Initiate the force from the input file. The second block in the file
    !contains the necessary information to setup the force.

    use mpi

    class(LAMMPS_Container), intent(inout) :: this
    character(len=*), intent(in) :: inputFile

    integer :: blockSize, i, ierr
    logical :: exists
    character(len=200) command

    inquire(file=inputFile, exist=exists)
    if (.not. exists) then
      write(*,'(A)') 'File '//trim(inputFile)//' cannot be found.'
      call MPI_FINALIZE(ierr)
      stop
    end if

    open(unit=21, file=inputFile)

    !Skip over the first block
    read(21, *) blockSize
    do i = 1, blockSize
      read(21, '(A)') command
    end do
    read(21, *) !skip spacer line

    !Read the force block
    read(21, *) blockSize
    do i = 1, blockSize
      read(21, '(A)') command
      call lammps_command(this%lammpsInstance, trim(command))
    end do

    close(21)

  end subroutine readForce
!***********************************************************************!
  subroutine readRestartSystem(this, file)
    !Setup the system from a restart file. The path of the data file
    !can be found in the first line of the provided restart file.

    class(LAMMPS_Container), intent(inout) :: this
    character(len=*), intent(in) :: file

    logical :: exists
    character(len=200) path

    inquire(file=file, exist=exists)
    if (.not. exists) then
      write(*,'(A)') 'File '//trim(file)//' cannot be found.'
      stop
    end if

    open(unit=22, file=file)
    read(22, '(A)') path
    call this%readData(path)
    close(22)

  end subroutine readRestartSystem
!***********************************************************************!
  subroutine getBoxDims(this)
    ! Retrieve the system dimension lengths from the LAMMPS instance

    use system

    class(LAMMPS_Container), intent(in) :: this

    real*8 :: xy = 0, yz = 0, xz = 0
    real*8, dimension(3) :: lmpBoxHi = 0, lmpBoxLo = 0
    logical, dimension(3) :: periodicity = 0
    logical :: box_change = 0

    call lammps_extract_box(this%lammpsInstance, lmpBoxLo, lmpBoxHi, &
                              xy, yz, xz, periodicity, box_change)

    BOXX = (lmpBoxHi(1) - lmpBoxLo(1))/(rlat)
    BOXY = (lmpBoxHi(2) - lmpBoxLo(2))/(rlat)
    BOXZ = (lmpBoxHi(3) - lmpBoxLo(3))/(rlat)

  end subroutine getBoxDims
!***********************************************************************!
  function getLammpsPtr(this) result(lammpsPtr)

    use, intrinsic :: iso_c_binding, only: C_NULL_PTR

    class(LAMMPS_Container), intent(in) :: this
    type(C_ptr) :: lammpsPtr

    lammpsPtr = this%lammpsInstance

  end function getLammpsPtr

!***********************************************************************!

  function getCommunicator(this) result(comm)

    class(LAMMPS_Container), intent(in) :: this
    integer :: comm

    comm = this%communicator

  end function getCommunicator

!***********************************************************************!v
  
    function getControlComm(this) result(comm)
  
      class(LAMMPS_Container), intent(in) :: this
      integer :: comm

      if (.not. this%isDimer) write(*,'(A)') 'Control communicator only available for dimer containers'
      comm = this%controlCommunicator
  
    end function getControlComm

!***********************************************************************!

  function getRank(this) result(rank)

    class(LAMMPS_Container), intent(in) :: this
    integer :: rank

    rank = this%rank

  end function getRank

!***********************************************************************!

  subroutine getPositions(this)

    use SYSTEM, only: NATOMC, RA, RLAT

    implicit none

    class(LAMMPS_Container), intent(in) :: this
    integer :: i
    real*8, dimension(:), allocatable :: pos

    call lammps_gather_atoms(this%lammpsInstance, 'x', 3, pos)

    DO I = 1, NATOMC
      RA(:, i) = pos(3*(i - 1) + 1:3*(i - 1) + 3)/RLAT
    end do

  end subroutine getPositions

!***********************************************************************!

  subroutine assignPositions(this)

    use MPI
    use SYSTEM, only: NATOMC, RA, RLAT, IERR

    implicit none

    class(LAMMPS_Container), intent(in) :: this
    integer :: i
    real*8, dimension(3*NATOMC) :: pos

    call mpi_bcast(RA, size(RA), MPI_REAL8, 0, this%communicator, ierr)

    do i = 1, NATOMC
      pos(3*(i - 1) + 1:3*(i - 1) + 3) = RA(:, i)*RLAT
    end do

    call lammps_scatter_atoms(this%lammpsInstance, 'x', pos)

  end subroutine assignPositions

!***********************************************************************!

  subroutine relaxSystemEnergy(this)

    ! This subroutine is valid for non-dimer lammps containers

    use SYSTEM, only: ISTEP, NFREQ, ENCALC, PRESS, RANK

    implicit none

    class(LAMMPS_Container), intent(inout) :: this

    call lammps_command(this%lammpsInstance, 'min_style cg')
    call lammps_command(this%lammpsInstance, 'minimize 1.0e-13 1.0e-13 2000 15000')

    call lammps_extract_compute(this%pressure, this%lammpsInstance, 'thermo_press', LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR)

    ENCALC = this%energy
    PRESS = -this%pressure !reversed sign; positive is tension

  end subroutine relaxSystemEnergy

!***********************************************************************!

  subroutine relaxSystemPressure(this)
    !Pressure relax the initial system by changing box size.

    use SYSTEM
    use MPI

    implicit none

    class(LAMMPS_Container), intent(inout) :: this

    character(len=200) :: str

    !Pressure relax if enabled. Otherwise energy relax.
    write (str, '(A,F10.3,A,F10.3,A,F10.3)') 'x ', SIGX, ' y ', SIGY, ' z ', SIGZ
    if (RANK .eq. 0 .and. SigRelOn) write (*, '(A)') &
    'fix f_boxrelax all box/relax '//trim(str)//' fixedpoint 0 0 0'

    if (SigRelOn) call lammps_command(this%lammpsInstance, &
    'fix f_boxrelax all box/relax '//trim(str)//' fixedpoint 0 0 0')
    call this%relaxSystemEnergy
    if (SigRelOn) call lammps_command(this%lammpsInstance, 'unfix f_boxrelax')

    if (rank .eq. 0 .and. SigRelOn) then
      write(*,'(6A14)') 'pxx','pyy','pzz','pxy','pyz','pxz'
      write(*,'(6E14.5)') PRESS(1), PRESS(2), PRESS(3), PRESS(4), PRESS(5), PRESS(6)
    end if

    ! Get the new system dimensions
    call this%getBoxDims()

    if (rank .eq. 0 .and. SigRelOn) write(*,'(A,3F15.8,A)') &
    'New system box size (', BOXX, BOXY, BOXZ, ' )'


  end subroutine relaxSystemPressure

!***********************************************************************!
  subroutine changeBoxDimensions(this, sizeX, sizeY, sizeZ)
    ! Change the system box dimensions with the following scaling factors
    ! for each dimension. The new dimensions are in units of lattice.

    use system, only: rlat

    class(LAMMPS_Container), intent(inout) :: this
    real*8, intent(in) :: sizeX, sizeY, sizeZ

    character(len=200) :: str, command
    
    write(str, '(3(A,F0.14))')  ' x final 0.0 ', sizeX*rlat, &
                                ' y final 0.0 ', sizeY*rlat, &
                                ' z final 0.0 ', sizeZ*rlat

    write(command, '(3A)') 'change_box all', trim(str), ' units box'

    call lammps_command(this%lammpsInstance, command)

  end subroutine
!***********************************************************************!
  subroutine writeData(this, filepath)
    !Create the data restart files needed to continue a simulation on
    !another run. See https://docs.lammps.org/write_data.html for more
    !information is written.

    use system, only: istep

    class(LAMMPS_Container), intent(inout) :: this
    character(len=*), intent(out) :: filepath

    character(len=200) :: writeCommand

    write(filepath, '(A,I0)') './RESTART/data.', istep
    write(writeCommand, '(A)') 'write_data '//trim(filepath)

    call lammps_command(this%lammpsInstance, writeCommand)

  end subroutine writeData
!***********************************************************************!
  subroutine readData(this, filepath)
    !Read a data file to continue a previous simulation. This preserves
    !box size and atom positions more consistently. See the following
    !https://docs.lammps.org/read_data.html for what data is read.

    use system, only: rank

    class(LAMMPS_Container), intent(inout) :: this
    character(len=*), intent(out) :: filepath

    character(len=200) :: readCommand

    if (rank == 0 .and. this%isSystem) then
      write(*,'(A)') 'Reading data file '//trim(filepath)
    end if
    write(readCommand, '(A)') 'read_data '//trim(filepath)

    call lammps_command(this%lammpsInstance, readCommand)

  end subroutine readData
!***********************************************************************!
end module LMP_Driver