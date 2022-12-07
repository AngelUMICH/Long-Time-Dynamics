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
!!READING SYSTEM INPUT FROM VARIOUS FILES

      SUBROUTINE SYS_INPUT

      USE SYSTEM
      USE LMP_DRIVER

      IMPLICIT NONE

      INTEGER I,J

!!PARAMETER INPUT FILE,MASTER CONTROL FILE

      OPEN(UNIT=10,FILE='par_init.dat')

!!READING SYSTEM PARAMTER, SEE PAR_INIT.DAT FOR DETAILS

      READ(10,'(8(I10/),5(F10.3/),(I10/),9(F14.6/),(I10/),1(ES15.5))')
     #TJOB,TDAM,TPBC,NTSTEP,NFREQ,NDEMAX,NAVMAX,NDIMER,
     #RLAT,RBOX,PCUTR,DCUTR,DECUT,NDRMAX,RDIMER,FNRMAX,FNRMIN,DSTEPF,
     #DSTEPM,EMSKIP,TEMP,PKAE,DRATE,NFREKL,TIME

!!READING INTITIAL STRUCTURE

      NATOMC = lammps_get_natoms(lmp)

      ALLOCATE(RA(3,NATOMC+NDEMAX))

      CALL RECV_LMP_POS

!!READING REFERENCE STRUCTURE
            
      NATOMR = 2*RBOX**3

      ALLOCATE(RR(3,NATOMR))

      CALL CREATE_REF_LATTICE

      RETURN

      END SUBROUTINE SYS_INPUT
!***********************************************************************!
!!SYSTEM INITITIALIZATION

      SUBROUTINE SYS_INIT

      USE MPI
      USE SYSTEM
      USE DEFECT
      USE KMCVAR

      IMPLICIT NONE

      INTEGER :: N,I
      INTEGER,ALLOCATABLE :: SEED(:)

!!INITIALIZE OUTPUT FILE

!!DEFECT STRUCTURE FILE

      OPEN(UNIT=22,FILE='std_kdef.out')

!!CONTAINS INITIAL, SADDLE, AND FINAL CONFIGURATION ENERGIES

      OPEN(UNIT=25,FILE='std_evo.out')

!!CONTAINS THE OPT STRUCTURE AFTER RELAXATION

      OPEN(UNIT=26,FILE='std_kopt.out')

!!CREATE THE NECESSARY DIRECTORIES

      IF (RANK .EQ. 0) THEN

        CALL EXECUTE_COMMAND_LINE('mkdir DUMP')

        CALL EXECUTE_COMMAND_LINE('mkdir RESTART')

      END IF
!!EITHER RANDOM SEED OR FIXED RANDOM SEED

      CALL RANDOM_SEED()

!      SEED(1) = 214748722

      CALL RANDOM_SEED(SIZE=N)
      ALLOCATE(SEED(N))
      CALL RANDOM_SEED(GET=SEED)


      FDCONV     = 1.0D-5

!!HALF OF THE BOX SIZE, FOR PBC USE

      RBH1       = RBOX/2.0D0

      RBH2       =-RBOX/2.0D0

!!POTENTIAL CUTOFF:PAIR TERM

      PCUTR      = PCUTR/RLAT

!!POTENTIAL CUTOFF:DENSITY TERM

      DCUTR      = DCUTR/RLAT

      PCUTS      = PCUTR*PCUTR

      DCUTS      = DCUTR*DCUTR

!!DEFECT CUTOFF

      DECUT      = DECUT*DECUT

      NDVAC      = 0

      NDINT      = 0

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

!!VARIABLES ASSSOCIATED WITH MPI
!! DIMER LIMIT VARIABLES BEING SET UP FOR J LOOP

      DIMERPR = NDIMER/NPROCS

      LDIMERI = DIMERPR*RANK + 1

      UDIMERI = DIMERPR*(RANK + 1)

!!SETTING UP COMMS FOR DIMER LAMMPS OBJECTS

      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,RANK,RANK,DIMER_COMM,IERR)

!!ALLOCATE DYNAMIC MATRIX
!!FROM HERE,SEE DEFINE.F FOR THE MEANING OF VARIABLES

      ALLOCATE(FQ(NDIMER));FQ=0.0D0

      ALLOCATE(EM(NDIMER));EM=0.0D0

      ALLOCATE(OE(NDIMER));OE=0.0D0

      ALLOCATE(FF(3,NATOMC));FF=0.0D0

      ALLOCATE(EE(3,NATOMC));EE=0.0D0

      ALLOCATE(TV(NATOMR));TV=0

      ALLOCATE(TI(NATOMC));TI=0

      ALLOCATE(IV(3,NDEMAX));IV=0

      ALLOCATE(II(3,NDEMAX));II=0

      ALLOCATE(PV(3,NDEMAX));PV=0.0D0

      ALLOCATE(PI(3,NDEMAX));PI=0.0D0

      ALLOCATE(PX(7,NDEMAX));PX=0.0D0

      ALLOCATE(OX(7));OX=0.0D0

      ALLOCATE(ND(NDIMER));ND=0

      ALLOCATE(ON(NDIMER));ON=0

      ALLOCATE(SR(4,NAVMAX,NDIMER));SR=0.0D0

      ALLOCATE(OR(4,NAVMAX,NDIMER));OR=0.0D0

      IF (RANK .EQ. 0) THEN

            ALLOCATE(SRBUFF(4,NAVMAX,NDIMER));SRBUFF=0.0D0

            ALLOCATE(NDBUFF(NDIMER));NDBUFF=0

            ALLOCATE(EMBUFF(NDIMER));EMBUFF=0.0D0

      END IF

!!VARIABLES ASSOCIATED WITH DIMER

      RDIMER     = RDIMER/RLAT

      DSTEPF     = DSTEPF/RLAT

      DSTEPM     = DSTEPM/RLAT

!!SETUP FOR PRE-SELECTION PROCESS

      SKIP = .FALSE.

      ALLOCATE(EMIG(2,50)); EMIG=0
      ALLOCATE(EMIC(2,50)); EMIC=1

      EMIG(1,1)=0.56
      EMIG(1,2)=0.64
      EMIG(1,3)=0.67
      EMIG(1,4)=0.438

      do i = 5,50
       EMIG(1,i)=3.5
      enddo

      EMIG(2,1)=0.30760
      EMIG(2,2)=0.04050
      EMIG(2,3)=0.03466
      EMIG(2,4)=0.20219
      EMIG(2,5)=0.01495
      EMIG(2,6)=0.08075
      EMIG(2,7)=0.80848
      EMIG(2,8)=0.05914
      EMIG(2,9)=0.13892
      EMIG(2,10)=0.22239
      EMIG(2,11)=0.0535
      EMIG(2,12)=0.09017
      EMIG(2,13)=0.08948
      EMIG(2,14)=0.09223
      EMIG(2,15)=0.031
      EMIG(2,16)=0.031
      EMIG(2,17)=0.031
      EMIG(2,18)=0.031
      EMIG(2,19)=0.031
      EMIG(2,20)=0.031

      do i = 21,50
       EMIG(2,i)=3.0
      enddo

      RETURN

      END SUBROUTINE SYS_INIT
!***********************************************************************!
!BCC REFERENCE LATTICE CREATION ONLY
      SUBROUTINE CREATE_REF_LATTICE

      USE SYSTEM, ONLY: RR, RBOX, NATOMR,RANK

      IMPLICIT NONE

      INTEGER :: I,J,K,COUNT
      REAL*8  :: X1,Y1,Z1,X2,Y2,Z2

      X1 = 0; Y1=0; Z1=0
      X2 = 0.5; Y2=0.5; Z2=0.5
      COUNT=1

      DO I = 1, RBOX
            DO J = 1, RBOX
                  DO K = 1, RBOX

                  RR(1:3,COUNT) = (/X1+(I-1),Y1+(J-1),Z1+(K-1)/)
                  RR(1:3,COUNT+1) = (/X2+(I-1),Y2+(J-1),Z2+(K-1)/)

                  COUNT = COUNT + 2

                  END DO
            END DO
      END DO

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
