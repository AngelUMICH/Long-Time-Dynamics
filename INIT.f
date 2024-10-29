!***********************************************************************!
!***********************************************************************!
!***          SELF EVOLVING ATOMISTIC KINETIC MONTE CARLO CODE       ***!
!***                        AUTHOR:XU,HAIXUAN                        ***!
!***             MATERIALS SCIENCE AND TECHNOLOGY DIVISION           ***!
!***                      OAK RIDGE NATIONAL LAB                     ***!
!***                          XUH1@ORNL.GOV                          ***!
!***                           VERSION 3.0                           ***!
!***                       ALL RIGHTS RESERVED                       ***!
!***********************************************************************!
!***********************************************************************!

!***********************************************************************!

SUBROUTINE READ_INPUT

  USE SYSTEM

  !!PARAMETER INPUT FILE,MASTER CONTROL FILE
  OPEN (UNIT=10, FILE='par_init.dat')

  !!READING SYSTEM PARAMTER, SEE PAR_INIT.DAT FOR DETAILS
  READ (10, '( 12(I10/),5(F10.3/),2(I10/),8(F14.6/),(I10/),(ES15.5/),(L/),2(F10.6/),(F10.6) )') &
    REST, TJOB, TDAM, TPBC, TCRY, NTSTEP, NFREQ, NDEMAX, NAVMAX, NDIMER, TLAT, MGSZ, &
    RLAT, RBOX, PCUTR, DCUTR, DECUT, NDSMAX, NDRMAX, RDIMER, FNRMAX, FNRMIN, DSTEPF, &
    DSTEPM, TEMP, PKAE, DRATE, NFREKL, TIME, SIGRELON, SIGX, SIGY, SIGZ

  CLOSE (10)

END SUBROUTINE READ_INPUT
!***********************************************************************!
SUBROUTINE LMP_INIT

  USE MPI
  USE LMP_Driver
  USE SYSTEM

  IMPLICIT NONE

  INTEGER :: EXTRA

  ! CREATE LAMMPS INSTANCES FOR SYSTEM AND DIMER OBJECTS
  call lmpSystem%init(MPI_COMM_WORLD, 0, 'system')
  call lmpSystem%setupSystem

  IF (NDIMER .GE. NPROCS) THEN

    DIMERPR = NDIMER/NPROCS
    EXTRA = MOD(NDIMER, NPROCS)

    LDIMERI = RANK*DIMERPR + MIN(RANK, EXTRA) + 1
    IF (RANK .LT. EXTRA) THEN
      UDIMERI = LDIMERI + DIMERPR
    ELSE
      UDIMERI = LDIMERI + DIMERPR - 1
    END IF
  ELSE IF ( NDIMER .LT. NPROCS) THEN

    LDIMERI = MOD(RANK, NDIMER) + 1
    UDIMERI = LDIMERI
  END IF

  call lmpDimer%init(MPI_COMM_WORLD, LDIMERI, 'dimer')
  call lmpDimer%setupSystem

END SUBROUTINE LMP_INIT
!***********************************************************************!
!!READING SYSTEM INPUT FROM VARIOUS FILES

SUBROUTINE SYS_INPUT

  USE SYSTEM
  USE LMP_Driver

  IMPLICIT NONE

  !!READING INTITIAL STRUCTURE

  NATOMC = int(lammps_get_natoms(lmpSystem%getLammpsPtr()))
  if (RANK .eq. 0) write(*, '(1I0, A)') NATOMC, " number of atoms"
  allocate (RA(3, NATOMC + NDEMAX))
  call lmpSystem%getPositions
  

END SUBROUTINE SYS_INPUT
!***********************************************************************!
!!SYSTEM INITITIALIZATION

SUBROUTINE SYS_INIT

  USE MPI
  USE SYSTEM
  USE DEFECT
  USE KMCVAR

  IMPLICIT NONE

  INTEGER :: N, I, NEMIG
  INTEGER, ALLOCATABLE :: SEED(:)
  LOGICAL :: EXISTS

  !!CREATE THE NECESSARY DIRECTORIES

  IF (RANK .EQ. 0) THEN
    CALL EXECUTE_COMMAND_LINE('mkdir -p DUMP')
    CALL EXECUTE_COMMAND_LINE('mkdir -p DUMP/DXA')
    CALL EXECUTE_COMMAND_LINE('mkdir -p RESTART')
  END IF

  !!EITHER RANDOM SEED OR FIXED RANDOM SEED

  CALL RANDOM_SEED()
  CALL RANDOM_SEED(SIZE=N)
  ALLOCATE (SEED(N))
  CALL RANDOM_SEED(GET=SEED)

  BSTEP = 0

  FDCONV = 1.0D-5

  !!HALF OF THE BOX SIZE, FOR PBC USE

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

  ! INPLACE TO ALLOCATE MEMORY. WILL BE CHANGED AFTER RELAXATION.
  NATOMR = NATOMC

  !!ALLOCATE DYNAMIC MATRIX
  !!FROM HERE,SEE DEFINE.F FOR THE MEANING OF VARIABLES

  ALLOCATE (FQ(NDIMER)); FQ = 0.0D0
  ALLOCATE (EM(NDIMER)); EM = 0.0D0
  ALLOCATE (OE(NDIMER)); OE = 0.0D0
  ALLOCATE (FF(3, NATOMC)); FF = 0.0D0
  ALLOCATE (EE(3, NATOMC)); EE = 0.0D0
  ALLOCATE (TV(NATOMR)); TV = 0
  ALLOCATE (TI(NATOMC)); TI = 0
  ALLOCATE (RNEIGH(NNEIGH, NATOMC)); RNEIGH = 0
  ALLOCATE (IV(3, NDEMAX)); IV = 0
  ALLOCATE (II(3, NDEMAX)); II = 0
  ALLOCATE (OC(10, NATOMR)); OC = 0
  ALLOCATE (OCC(NATOMR)); OCC = 0
  ALLOCATE (PV(3, NDEMAX)); PV = 0.0D0
  ALLOCATE (PI(3, NDEMAX)); PI = 0.0D0
  ALLOCATE (PX(7, NDEMAX)); PX = 0.0D0
  ALLOCATE (OX(7)); OX = 0.0D0
  ALLOCATE (ND(NDIMER)); ND = 0
  ALLOCATE (ON(NDIMER)); ON = 0
  ALLOCATE (SR(4, NAVMAX, NDIMER)); SR = 0.0D0
  ALLOCATE (OR(4, NAVMAX, NDIMER)); OR = 0.0D0

  ALLOCATE (TVBUFF(NATOMR)); TVBUFF = 0
  ALLOCATE (TIBUFF(NATOMC)); TIBUFF = 0
  ALLOCATE (SRBUFF(4, NAVMAX, NDIMER)); SRBUFF = 0.0D0
  ALLOCATE (NDBUFF(NDIMER)); NDBUFF = 0
  ALLOCATE (EMBUFF(NDIMER)); EMBUFF = 0.0D0

  !!VARIABLES ASSOCIATED WITH DIMER

  RDIMER = RDIMER/RLAT
  DSTEPF = DSTEPF/RLAT
  DSTEPM = DSTEPM/RLAT

  !!SETUP FOR PRE-SELECTION PROCESS

  SKIP = .FALSE.

  ALLOCATE (EMIG(2, 50)); EMIG = 0
  ALLOCATE (EMIC(2, 50)); EMIC = 1

  INQUIRE (FILE='emig.dat', EXIST=EXISTS)
  IF (.NOT. EXISTS) THEN
    WRITE (*, *) "File 'emig.dat' cannot be found"
    STOP
  END IF

  OPEN (UNIT=11, FILE='emig.dat')

  !READ VAC EMIGS
  READ (11, *)
  READ (11, *) NEMIG

  DO I = 1, NEMIG
    READ (11, *) EMIG(1, I)
  END DO
  DO I = NEMIG + 1, 50
    EMIG(1, I) = 3.5
  END DO

  !READ INT EMIGS
  READ (11, *)
  READ (11, *) NEMIG

  DO I = 1, NEMIG
    READ (11, *) EMIG(2, I)
  END DO
  DO I = NEMIG + 1, 50
    EMIG(2, I) = 3.0
  END DO

  CLOSE (11)

  CALL READ_RESTART

  RETURN

END SUBROUTINE SYS_INIT
!***********************************************************************!
SUBROUTINE READ_RESTART

  USE SYSTEM
  USE DEFECT
  USE LMP_DRIVER

  IMPLICIT NONE

  INTEGER :: I, J
  LOGICAL :: EXISTS
  CHARACTER(LEN=200) :: STR

  !!UPDATE PARAMETERS FROM RESTART FILE

  IF (REST .EQ. 2) THEN

    IF (RANK .EQ. 0) WRITE (*, *) 'Running in restart mode.'

    INQUIRE (FILE='RESTART.OUT', EXIST=EXISTS)

    IF (.NOT. EXISTS) THEN

      WRITE (*, *) "File 'RESTART.OUT' cannot be found."

      STOP

    END IF

    OPEN (UNIT=27, FILE='RESTART.OUT')

    !READ POSITIONS

    READ (27, *) NATOMC

    DO I = 1, NATOMC

      READ (27, '(3F14.5)') (RA(J, I), J=1, 3)

    END DO

    !READ PREVIOUS SIMULATION PARAMETERS

    READ (27, *)
    READ (27, *)

    READ (27, *) BSTEP
    READ (27, *) TIME
    READ (27, *) TIME_LAST
    READ (27, *) DTIME
    READ (27, *) DPA
    READ (27, *) DPAC
    READ (27, *) AVDPA
    READ (27, *) ONDSTOT
    READ (27, *) EPSX
    READ (27, *) EPSY
    READ (27, *) EPSZ

    !READ MIGRATION ENERGIES

    READ (27, *)
    READ (27, *)

    DO I = 1, 50
      READ (27, *) EMIG(1, I), EMIC(1, I)
    END DO

    DO I = 1, 50
      READ (27, *) EMIG(2, I), EMIC(2, I)
    END DO

    CLOSE (27)

    ! READ REFERENCE LATTICE

    open(unit=21, file='pos_ref.dat')
    read(21, *) natomr
    read(21, *)
    if (.not. allocated(RR)) allocate(RR(3,NATOMR))
    do i= 1, natomr
      read(21, *) (RR(j,i), j=1,3)
    end do
    close(21)

    IF (EPSX .NE. 1.0 .OR. EPSY .NE. 1.0 .OR. EPSZ .NE. 1.0) THEN
      DO I = 1, NATOMR
        RR(1, I) = RR(1, I)*EPSX
        RR(2, I) = RR(2, I)*EPSY
        RR(3, I) = RR(3, I)*EPSZ
      END DO

      !Change the simualtion dimensions with corresponding strain
      write (str, '(A,F10.8,A,F10.8,A,F10.8)') 'x final 0 ', RBOX*EPSX, ' y final 0 ', RBOX*EPSY, ' z final 0 ', RBOX*EPSZ
      if (rank .eq. 0) print *, 'change_box all '//trim(str)//' units lattice'
      call lammps_command(lmpSystem%getLammpsPtr(), 'change_box all '//trim(str)//' units lattice')
      call lammps_command( lmpDimer%getLammpsPtr(), 'change_box all '//trim(str)//' units lattice')
    END IF

    !update atom count in lammps
    !need to change atom type if creating different atom types
    !at the moment only atom type 1 are made
    if (natomc .lt. int(lammps_get_natoms(lmpSystem%getLammpsPtr()))) then
      write (str, '(I0)') NATOMC
      call lammps_command(lmpSystem%getLammpsPtr(), 'delete_atoms group all')
      call lammps_command(lmpSystem%getLammpsPtr(), 'create_atoms 1 random '//trim(str)//' 384753 NULL')
    end if

    call lmpSystem%assignPositions

  END IF

END SUBROUTINE READ_RESTART
!***********************************************************************!
!BCC REFERENCE LATTICE CREATION ONLY
SUBROUTINE CREATE_REF_LATTICE
! Create the initial reference system that will be used to identify point defects
! and their clusters.

  USE SYSTEM
  use LMP_DRIVER
  use MicroStrucFeatures
  use MPI

  IMPLICIT NONE

  integer :: i,j
  real*8 :: bounds(6)
  logical :: exists

  real*8 x, y, z

  ! Create the reference lattice by first removing the grain
  ! boundaries. Then create the pos_ref.dat file for future use.

  ! First create the reference lattice thorugh features object
  call features%grains%init(lmpSystem, MGSZ)
  call features%grains%identifyGrains

  bounds = (/0.0D0, RBOX*EPSX, 0.0D0, RBOX*EPSY, 0.0D0, RBOX*EPSZ/)
  call features%reference%createReferenceLattice(features%grains,bounds,RA,natomc)

  NATOMR = features%reference%numPoints
  if (.not. allocated(RR)) allocate(RR(3,NATOMR))
  RR(:,:NATOMR) = features%reference%lattice(:,:NATOMR)

  ! Create the pos_ref.dat reference file in the current directory.
  inquire (file='pos_ref.dat', exist=exists)
  call MPI_BCAST(exists, 1, MPI_LOGICAL, 0, lmpSystem%getCommunicator(), ierr)
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
    CALL MPI_FINALIZE(ierr)
    stop
  end if

  if (rank == 0) then
    open (unit=21, file='pos_ref.dat', status='new')
    write (21, '(I0)') NATOMR
    write (21, *)
    write (21, '(3F14.5)') (RR(1, i), RR(2, i), RR(3, i), i=1, NATOMR)
    close(21)
  end if

END SUBROUTINE CREATE_REF_LATTICE

!***********************************************************************!

!***********************************************************************!
!***********************************************************************!
!***          SELF EVOLVING ATOMISTIC KINETIC MONTE CARLO CODE       ***!
!***                        AUTHOR:XU,HAIXUAN                        ***!
!***             MATERIALS SCIENCE AND TECHNOLOGY DIVISION           ***!
!***                      OAK RIDGE NATIONAL LAB                     ***!
!***                          XUH1@ORNL.GOV                          ***!
!***                           VERSION 3.0                           ***!
!***                       ALL RIGHTS RESERVED                       ***!
!***********************************************************************!
!***********************************************************************!
