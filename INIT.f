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
      OPEN(UNIT=10,FILE='par_init.dat')

      !!READING SYSTEM PARAMTER, SEE PAR_INIT.DAT FOR DETAILS
      READ(10,'( 10(I10/),5(F10.3/),2(I10/),8(F14.6/),(I10/),(ES15.5/),2(F10.6/),(F10.6) )')
     #REST,TJOB,TDAM,TPBC,TCRY,NTSTEP,NFREQ,NDEMAX,NAVMAX,NDIMER,
     #RLAT,RBOX,PCUTR,DCUTR,DECUT,NDSMAX,NDRMAX,RDIMER,FNRMAX,FNRMIN,DSTEPF,
     #DSTEPM,TEMP,PKAE,DRATE,NFREKL,TIME,SIGX,SIGY,SIGZ

      END SUBROUTINE READ_INPUT
!***********************************************************************!
!!READING SYSTEM INPUT FROM VARIOUS FILES

      SUBROUTINE SYS_INPUT

      USE SYSTEM
      USE LMP_DRIVER

      IMPLICIT NONE

!!READING INTITIAL STRUCTURE

      NATOMC =int(lammps_get_natoms(lmp))

      ALLOCATE(RA(3,NATOMC+NDEMAX))

      CALL RECV_LMP_POS

!!CREATING REFERENCE STRUCTURE

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
!      USE LMP_DRIVER

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
        
        CALL EXECUTE_COMMAND_LINE('mkdir DUMP/DXA')

        CALL EXECUTE_COMMAND_LINE('mkdir RESTART')

      END IF

!!EITHER RANDOM SEED OR FIXED RANDOM SEED

      CALL RANDOM_SEED()

      CALL RANDOM_SEED(SIZE=N)
      ALLOCATE(SEED(N))
      CALL RANDOM_SEED(GET=SEED)

      BSTEP = 0

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

!! KMC ENERGY BARRIER

      EMSKIP = -1*8.617E-5*TEMP*LOG(DRATE/1E12)

!! SYSTEM PROPERTIES

      PRESS=0.0D0

      EPSX = 1

      EPSY = 1

      EPSZ = 1

!!VARIABLES ASSSOCIATED WITH MPI
!! DIMER LIMIT VARIABLES BEING SET UP FOR J LOOP

      DIMERPR = NDIMER/NPROCS

      LDIMERI = DIMERPR*RANK + 1

      UDIMERI = DIMERPR*(RANK + 1)

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

      ALLOCATE(OC(10,NATOMR));OC=0

      ALLOCATE(OCC(NATOMR));OCC=0

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

      EMIG(1,1)=0.637
      EMIG(1,2)=0.462
      EMIG(1,3)=0.604
      EMIG(1,4)=0.612

      do i = 5,50
       EMIG(1,i)=3.5
      enddo

      EMIG(2,1)=0.167
      EMIG(2,2)=0.083
      EMIG(2,3)=0.054
      EMIG(2,4)=0.043
      EMIG(2,5)=0.021
      EMIG(2,6)=0.006
      EMIG(2,7)=0.005
      EMIG(2,8)=0.002
      EMIG(2,9)= 0.006
      EMIG(2,10)=0.016
      EMIG(2,11)=0.0535
      EMIG(2,12)=0.031
      EMIG(2,13)=0.031
      EMIG(2,14)=0.031
      EMIG(2,15)=0.031
      EMIG(2,16)=0.031
      EMIG(2,17)=0.031
      EMIG(2,18)=0.031
      EMIG(2,19)=0.031
      EMIG(2,20)=0.031

      do i = 21,50
       EMIG(2,i)=3.0
      enddo

      CALL READ_RESTART

      !UPDATE THE REFERENCE LATTICE
      !IF DUE TO RESTART READ STRAIN WILL NOT BE 1.
      DO I=1, NATOMC
        RR(1,I) = RR(1,I)*EPSX
        RR(2,I) = RR(2,I)*EPSY
        RR(3,I) = RR(3,I)*EPSZ
      END DO

      RETURN

      END SUBROUTINE SYS_INIT
!***********************************************************************!
      SUBROUTINE READ_RESTART

      USE SYSTEM
      USE DEFECT
      USE LMP_DRIVER

      IMPLICIT NONE

      INTEGER :: I,J
      LOGICAL :: EXISTS
      CHARACTER (LEN=200) :: STR

      !!UPDATE PARAMETERS FROM RESTART FILE

      IF (REST . EQ. 2) THEN

      INQUIRE(FILE='RESTART.OUT', EXIST=EXISTS)

      IF (.NOT. EXISTS) THEN

            WRITE(*,*) "File 'RESTART.OUT' cannot be found."
      
            STOP
      
      END IF

      OPEN(UNIT=27,FILE='RESTART.OUT')

      !READ POSITIONS

      READ(27,*) NATOMC

      DO I = 1,NATOMC

            READ(27,'(3F14.5)') (RA(J,I), J=1,3)

      END DO

      !READ PREVIOUS SIMULATION PARAMETERS

      READ(27,*)

      READ(27,*)

      READ(27,*) BSTEP
      READ(27,*) TIME
      READ(27,*) TIME_LAST
      READ(27,*) DTIME
      READ(27,*) DPA
      READ(27,*) DPAC
      READ(27,*) AVDPA
      READ(27,*) ONDSTOT
      READ(27,*) EPSX
      READ(27,*) EPSY
      READ(27,*) EPSZ

      !READ MIGRATION ENERGIES

      READ(27,*)

      READ(27,*)

      DO I = 1, 50
          READ(27,*) EMIG(1,I), EMIC(1,I)
      END DO

      DO I = 1, 50
          READ(27,*) EMIG(2,I), EMIC(2,I)
      END DO

      CLOSE(27)

      !Change the simualtion dimensions with corresponding strain
      WRITE(STR,'(A,F10.6,A,F10.6,A,F10.6)') ' x final 0 ',RBOX*EPSX,' y final 0 ',RBOX*EPSY,' z final 0 ',RBOX*EPSZ
      CALL lammps_command(lmp, 'change_box all'//trim(str)//' remap units lattice')
      CALL lammps_command(lmpDimer, 'change_box all'//trim(str)//' remap units lattice')

      CALL SEND_OPLD_POS

      END IF

      END SUBROUTINE READ_RESTART
!***********************************************************************!
!BCC REFERENCE LATTICE CREATION ONLY
      SUBROUTINE CREATE_REF_LATTICE

      USE SYSTEM, ONLY: RR,RBOX,NATOMR,TCRY

      IMPLICIT NONE

      INTEGER :: I,J,K,COUNT,SIZE
      REAL*8  :: X1,Y1,Z1,X2,Y2,Z2
      LOGICAL :: EXISTS

      X1 = 0; Y1=0; Z1=0
      X2 = 0.5; Y2=0.5; Z2=0.5
      COUNT=1
      SIZE = INT(RBOX)

      IF (TCRY .EQ. 1) THEN

            NATOMR = 2*SIZE**3
            ALLOCATE(RR(3,NATOMR))
            
            DO I = 1, SIZE
                  DO J = 1, SIZE
                        DO K = 1, SIZE

                        RR(1:3,COUNT) = (/X1+(I-1),Y1+(J-1),Z1+(K-1)/)
                        RR(1:3,COUNT+1) = (/X2+(I-1),Y2+(J-1),Z2+(K-1)/)
                        COUNT = COUNT + 2

                        END DO
                  END DO
            END DO

      ELSE IF (TCRY .EQ. 2) THEN

            INQUIRE(FILE='pos_ref.dat', EXIST=EXISTS)
            IF (.NOT. EXISTS) THEN
                  WRITE(*,*) "File 'pos_ref.dat' cannot be found."
                  WRITE(*,*) "A reference file is needed for polycrystal."
                  STOP
            END IF

            OPEN(UNIT=21, FILE='pos_ref.dat')
            READ(21,*) NATOMR
            ALLOCATE(RR(3,NATOMR))
            READ(21,*)(RR(1,I),RR(2,I),RR(3,I),I=1,NATOMR)

      END IF 

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
