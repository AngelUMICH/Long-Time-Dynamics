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
    procedure :: setupSystem
    procedure :: getLammpsPtr
    procedure :: getCommunicator
    procedure :: getControlComm
    procedure :: getRank
    procedure :: getPositions
    procedure :: assignPositions
    procedure :: relaxSystemEnergy
    procedure :: relaxSystemPressure

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

    call lammps_open('lmp -log none', this%communicator, this%lammpsInstance)

    if (compType == 's') then
      this%isSystem = .true.
    else if (compType == 'd') then
      this%isDimer = .true.
    end if

  end subroutine LAMMPS_Container_init
!***********************************************************************!

  subroutine setupSystem(this)

    use SYSTEM

    implicit none

    class(LAMMPS_Container), intent(inout) :: this
    character(len=200) :: str

    !Setting up system in LAMMPS space

    !Necessary commands for OPLD to work
    call lammps_command(this%lammpsInstance, 'units metal')
    call lammps_command(this%lammpsInstance, 'boundary p p p')
    call lammps_command(this%lammpsInstance, 'atom_style atomic')
    call lammps_command(this%lammpsInstance, 'atom_modify map array')

    if (TLAT .eq. 1) then
      write (str, '(a,x,f7.5)') 'bcc', rlat
    else if (TLAT .eq. 2) then
      write (str, '(a,x,f7.5)') 'fcc', rlat
    end if

    if (RANK .eq. 0 .and. this%isSystem) print *, 'lattice '//trim(str)
    call lammps_command(this%lammpsInstance, 'lattice '//trim(str))

    write (str, '(f3.1,x,f8.5,x,f3.1,x,f8.5,x,f3.1,x,f8.5)') 0.0, rbox, 0.0, rbox, 0.0, rbox
    if (RANK .eq. 0 .and. this%isSystem) print *, 'region box block '//trim(str)//' units lattice'
    call lammps_command(this%lammpsInstance, 'region box block '//trim(str)//' units lattice')
    call lammps_command(this%lammpsInstance, 'create_box 1 box')

    !Read from input file using LAMMPS
    if (RANK .eq. 0 .and. this%isSystem) print *, "Reading in initial positions"
    call lammps_file(this%lammpsInstance, 'in.lammps')
    ! call lammps_command(this%lammpsInstance, 'thermo_style custom step pe lx ly lz pxx pyy pzz')

    !Including necessary computes for OPLD
    call lammps_command(this%lammpsInstance, 'compute TotEn all pe')
    call lammps_command(this%lammpsInstance, 'compute TotPress all pressure NULL virial')

    !Initiating system values
    call lammps_command(this%lammpsInstance, 'run 0')

    !Extracting MD values for cascade
    call lammps_extract_compute(this%energy, this%lammpsInstance, 'thermo_pe', LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)
    call lammps_extract_compute(this%pressure, this%lammpsInstance, 'thermo_press', LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR)

    CALL LAMMPS_COMMAND(this%lammpsInstance, 'fix storePotVector all property/atom d_pot')

    if (this%isDimer) call lammps_command(this%lammpsInstance, 'delete_atoms group all')

  end subroutine setupSystem

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

    if (mod(ISTEP, NFREQ) .eq. 0 .and. RANK .eq. 0) write (*, '(A)') 'Relaxing System'

    call lammps_command(this%lammpsInstance, 'min_style cg')
    call lammps_command(this%lammpsInstance, 'minimize 1.0e-10 1.0e-10 1000 10000')

    ENCALC = this%energy
    PRESS = -this%pressure !reversed sign; positive is tension

    if (mod(ISTEP, NFREQ) .eq. 0 .and. RANK .eq. 0) write (*, '(A)') 'Finished Relaxation'

  end subroutine relaxSystemEnergy

!***********************************************************************!

  subroutine relaxSystemPressure(this)

    ! This subroutine is valid for non-dimer lammps containers

    use SYSTEM
    use MPI

    implicit none

    class(LAMMPS_Container), intent(inout) :: this
    integer :: i, j
    character(len=200) :: str

    real*8, dimension(3) :: lmpBoxLo = 0, lmpBoxHi = 0
    real*8 :: xy = 0, yz = 0, xz = 0
    logical, dimension(3) :: periodicity = 0
    logical :: box_change = 0

    if (ISTEP .eq. 1 .and. SigRelOn) then

      if (RANK .eq. 0) write (*, '(A)') 'Relaxing System'

      write (str, '(A,F10.3,A,F10.3,A,F10.3)') 'x ', SIGX, ' y ', SIGY, ' z ', SIGZ
      if (RANK .eq. 0) write (*, '(A)') 'fix f_boxrelax all box/relax '//trim(str)//' fixedpoint 0 0 0'

      call lammps_command(this%lammpsInstance, 'fix f_boxrelax all box/relax '//trim(str)//' fixedpoint 0 0 0')
      call lammps_command(this%lammpsInstance, 'min_style cg')
      call lammps_command(this%lammpsInstance, 'minimize 1.0e-8 1.0e-8 1000 10000')
      call lammps_command(this%lammpsInstance, 'unfix f_boxrelax')

      ENCALC = this%energy
      PRESS = -this%pressure !reversed sign; positive is tension

      if (RANK .eq. 0) write (*, '(A)') 'Finished Relaxation'

      ! Update the reference lattice
      do i = 1, NATOMR
        RR(1, i) = RR(1, i)/EPSX
        RR(2, i) = RR(2, i)/EPSY
        RR(3, i) = RR(3, i)/EPSZ
      end do

      call lammps_extract_box(this%lammpsInstance, lmpBoxLo, lmpBoxHi, xy, yz, xz, periodicity, box_change)

      EPSX = (lmpBoxHi(1) - lmpBoxLo(1))/(rbox*rlat)
      EPSY = (lmpBoxHi(2) - lmpBoxLo(2))/(rbox*rlat)
      EPSZ = (lmpBoxHi(3) - lmpBoxLo(3))/(rbox*rlat)

      do i = 1, NATOMR
        RR(1, i) = RR(1, i)*EPSX
        RR(2, i) = RR(2, i)*EPSY
        RR(3, i) = RR(3, i)*EPSZ
      end do

      ! Center box corner at origin
      write (str, '(A,F10.6,A,F10.6,A,F10.6)') 'x final 0 ', RBOX*EPSX, ' y final 0 ', RBOX*EPSY, ' z final 0 ', RBOX*EPSZ
      call mpi_bcast(str, len(str), MPI_CHARACTER, 0, this%communicator, ierr)
      call lammps_command(this%lammpsInstance, 'change_box all '//trim(str)//' units lattice')

    end if

  end subroutine relaxSystemPressure

!***********************************************************************!
end module LMP_Driver