subroutine readInput

  use SYSTEM

  !!PARAMETER INPUT FILE,MASTER CONTROL FILE
  open (UNIT=10, FILE='par_init.dat')

  !!readING SYSTEM PARAMTER, SEE PAR_INIT.DAT FOR DETAILS
  read (10, '( 11(I10/),5(F10.3/),2(I10/),8(F14.6/),(I10/),(ES15.5/),(L/),2(F10.6/),(F10.6) )') &
    REST, TDAM, TPBC, TCRY, NTSTEP, NFREQ, NDEMAX, NAVMAX, NDIMER, TLAT, MGSZ, &
    RLAT, RBOX, PCUTR, DCUTR, DECUT, NDSMAX, NDRMAX, RDIMER, FNRMAX, FNRMIN, DSTEPF, &
    DSTEPM, TEMP, PKAE, DRATE, NFREKL, TIME, SIGRELON, SIGX, SIGY, SIGZ

  close (10)

end subroutine readInput
!***********************************************************************!
subroutine lmpInit

  use MPI
  use LMP_Driver
  use SYSTEM

  implicit none

  integer :: EXTRA

  ! CREATE LAMMPS INSTANCES FOR SYSTEM AND DIMER OBJECTS
  call lmpSystem%init(MPI_COMM_WORLD, 0, 'system')

  if (NDIMER .GE. NPROCS) THEN

    DIMERPR = NDIMER/NPROCS
    EXTRA = MOD(NDIMER, NPROCS)

    LDIMERI = RANK*DIMERPR + MIN(RANK, EXTRA) + 1
    if (RANK .LT. EXTRA) THEN
      UDIMERI = LDIMERI + DIMERPR
    else
      UDIMERI = LDIMERI + DIMERPR - 1
    end if
  else if ( NDIMER .LT. NPROCS) THEN

    LDIMERI = MOD(RANK, NDIMER) + 1
    UDIMERI = LDIMERI
  end if

  call lmpDimer%init(MPI_COMM_WORLD, LDIMERI, 'dimer')

end subroutine lmpInit
!***********************************************************************!
subroutine sysInput

  use SYSTEM
  use LMP_Driver

  implicit none

  character(len=9) file1
  character(len=11) file2

  file1 = 'in.lammps'
  file2= 'RESTART.OUT'

  !Setup the box and atom positions
  call lmpSystem%setup(file1, file2, REST)
  call lmpDimer%setup(file1, file2, REST)

  NATOMC = int(lammps_get_natoms(lmpSystem%getLammpsPtr()))
  if (RANK .eq. 0) write(*, '(1I0, A)') NATOMC, " number of atoms"
  allocate (RA(3, NATOMC + NDEMAX))
  call lmpSystem%getPositions

end subroutine sysInput
!***********************************************************************!
!!SYSTEM INITITIALIZATION

subroutine sysInit

  use MPI
  use SYSTEM
  use DEFECT
  use KMCVAR

  implicit none

  integer :: N, I, NEMIG
  integer, allocatable :: SEED(:)
  logical :: EXISTS

  !!CREATE THE NECESSARY DIRECTORIES

  if (RANK .EQ. 0) THEN
    call execute_command_line('mkdir -p DUMP')
    call execute_command_line('mkdir -p DUMP/DXA')
    call execute_command_line('mkdir -p RESTART')
  end if

  !!EITHER RANDOM SEED OR FIXED RANDOM SEED

  call RANDOM_SEED()
  call RANDOM_SEED(SIZE=N)
  allocate (SEED(N))
  call RANDOM_SEED(GET=SEED)

  BSTEP = 0

  FDCONV = 1.0D-5

  !!HALF OF THE BOX SIZE, FOR PBC use

  RBH1 = RBOX/2.0D0
  RBH2 = -RBOX/2.0D0

  !!POTENTIAL CUTOFF:PAIR TERM

  PCUTR = PCUTR/RLAT

  !!POTENTIAL CUTOFF:DENSITY TERM

  DCUTR = DCUTR/RLAT
  PCUTS = PCUTR*PCUTR
  DCUTS = DCUTR*DCUTR

  !!DEFECT CUTOFF

  DECUT = DECUT*DECUT
  NDVAC = 0
  NDINT = 0
  NNEIGH = 10
  UPNEIGH = .TRUE.

  !!IRRADIATION VARIABLES

  TIME_LAST = TIME
  DTIME = 0
  AVDPA = 0
  DPA = 0
  DPAC = 0
  CALC = .FALSE.

  !!INITIALIZING DEFECT VARIABLES

  ONDSTOT = 0
  NDSTOT = 0
  ONDTOT = 0
  NDTOT = 0

  !! KMC ENERGY BARRIER

  EMSKIP = -1*8.617E-5*TEMP*LOG(DRATE/1E12)

  !! SYSTEM PROPERTIES

  PRESS = 0.0D0
  EPSX = 1.0
  EPSY = 1.0
  EPSZ = 1.0

  NATOMR = NATOMC

  call readRestartParameters()

  allocate(FQ(NDIMER)); FQ = 0.0D0
  allocate(EM(NDIMER)); EM = 0.0D0
  allocate(OE(NDIMER)); OE = 0.0D0
  allocate(FF(3, NATOMC)); FF = 0.0D0
  allocate(EE(3, NATOMC)); EE = 0.0D0
  allocate(TV(NATOMR)); TV = 0
  allocate(TI(NATOMC)); TI = 0
  allocate(RNEIGH(NNEIGH, NATOMC)); RNEIGH = 0
  allocate(IV(3, NDEMAX)); IV = 0
  allocate(II(3, NDEMAX)); II = 0
  allocate(OC(10, NATOMR)); OC = 0
  allocate(OCC(NATOMR)); OCC = 0
  allocate(PV(3, NDEMAX)); PV = 0.0D0
  allocate(PI(3, NDEMAX)); PI = 0.0D0
  allocate(PX(7, NDEMAX)); PX = 0.0D0
  allocate(OX(7)); OX = 0.0D0
  allocate(ND(NDIMER)); ND = 0
  allocate(ON(NDIMER)); ON = 0
  allocate(SR(4, NAVMAX, NDIMER)); SR = 0.0D0
  allocate(OR(4, NAVMAX, NDIMER)); OR = 0.0D0

  allocate(TVBUFF(NATOMR)); TVBUFF = 0
  allocate(TIBUFF(NATOMC)); TIBUFF = 0
  allocate(SRBUFF(4, NAVMAX, NDIMER)); SRBUFF = 0.0D0
  allocate(NDBUFF(NDIMER)); NDBUFF = 0
  allocate(EMBUFF(NDIMER)); EMBUFF = 0.0D0

  !!VARIABLES ASSOCIATED WITH DIMER

  RDIMER = RDIMER/RLAT
  DSTEPF = DSTEPF/RLAT
  DSTEPM = DSTEPM/RLAT

  !!Setup for defect selection energies

  SKIP = .FALSE.

  if (.not. allocated(EMIG)) allocate(EMIG(2, 50)); EMIG = 0
  if (.not. allocated(EMIC)) allocate(EMIC(2, 50)); EMIC = 1

  inquire (FILE='emig.dat', EXIST=EXISTS)
  if (.NOT. EXISTS) THEN
    write (*, *) "File 'emig.dat' cannot be found"
    stop
  end if

  open (UNIT=11, FILE='emig.dat')

  !Read vacancy mig. energies
  read (11, *)
  read (11, *) NEMIG

  do I = 1, NEMIG
    read (11, *) EMIG(1, I)
  end do
  do I = NEMIG + 1, 50
    EMIG(1, I) = 3.5
  end do

  !Read interstitial mig. energies
  read (11, *)
  read (11, *) NEMIG

  do I = 1, NEMIG
    read (11, *) EMIG(2, I)
  end do
  do I = NEMIG + 1, 50
    EMIG(2, I) = 3.0
  end do

  close (11)

  RETURN

end subroutine sysInit
!***********************************************************************!
subroutine readRestartParameters

  use SYSTEM
  use DEFECT
  use LMP_DRIVER
  use MicroStrucFeatures

  implicit none

  integer :: I, J
  logical :: EXISTS

  !!UPDATE PARAMETERS FROM RESTART FILE

  if (REST .EQ. 2) THEN
    inquire (FILE='RESTART.OUT', EXIST=EXISTS)
    if (.NOT. EXISTS) THEN
      write (*, *) "File 'RESTART.OUT' cannot be found."
      stop
    end if
    open (UNIT=27, FILE='RESTART.OUT')

    !Skip the first couple lines
    read(27, *)
    read(27, *)

    !Read simulation parameters
    read(27, *) !skip header
    read(27, *) BSTEP
    read(27, *) TIME
    read(27, *) TIME_LAST
    read(27, *) DTIME
    read(27, *) DPA
    read(27, *) DPAC
    read(27, *) AVDPA
    read(27, *) ONDSTOT
    read(27, *) EPSX
    read(27, *) EPSY
    read(27, *) EPSZ
    read(27, *)

    !Read defect migration energies
    if (.not. allocated(EMIG)) allocate(EMIG(2, 50))
    if (.not. allocated(EMIC)) allocate(EMIC(2, 50))
    
    read (27, *) !skip header
    do i = 1, 50
    
      read (27, *) EMIG(1, I), EMIC(1, I)
    end do

    read (27, *) !skip header
    do i = 1, 50

      read (27, *) EMIG(2, I), EMIC(2, I)
    end do
    read (27, *)
    
    !Read grain assignments
    read (27, *) !skip header
    read (27, *) NATOMC
    call features%grains%allocateClusterArray(NATOMC)
    do i = 1, NATOMC

      read (27, *) features%grains%grainOccupy(I)
    end do

    close (27)

    !Read the reference lattice
    open(unit=21, file='pos_ref.dat')
    read(21, *) natomr
    read(21, *)
    if (.not. allocated(RR)) allocate(RR(3,NATOMR))
    do i= 1, natomr

      read(21, *) (RR(j,i), j=1,3)
    end do
    close(21)

    call lmpSystem%getPositions()

  end if

end subroutine readRestartParameters
!***********************************************************************!