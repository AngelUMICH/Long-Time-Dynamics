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
!!APPLY THE PERIODIC BOUNDARY CONDITION TO THE SYSTEM

      REAL*8 FUNCTION PBC(X,DIM)

      USE SYSTEM,ONLY : RBOX,RBH1,RBH2,EPSX,EPSY,EPSZ

      IMPLICIT NONE

      INTEGER DIM
      REAL*8 X, EPS

      IF (DIM .EQ. 1) THEN
        EPS = EPSX
      ELSEIF (DIM .EQ. 2) THEN
        EPS = EPSY
      ELSEIF (DIM .EQ. 3) THEN
        EPS = EPSZ
      END IF
      
      IF((X.LE.RBH1*EPS).AND.(X.GE.RBH2*EPS)) THEN
        PBC = X
      ELSEIF(X.LT.RBH2*EPS) THEN
        PBC = X + RBOX*EPS
      ELSE
        PBC = X - RBOX*EPS
      ENDIF

      RETURN

      END
      
!***********************************************************************!
!!PERIODIC CONDIITONS ONLY USED FOR FINDING DEFECTS

      REAL*8 FUNCTION PBC_DEFECT(X)

      USE SYSTEM,ONLY : RBOX,RBH1,RBH2

      IMPLICIT NONE

      REAL*8 X

!      PBC_DEFECT = X-RBOX*NINT(X/RBOX)

      IF((X.LE.RBH1).AND.(X.GE.RBH2)) THEN
        PBC_DEFECT = X
      ELSEIF(X.LT.RBH2) THEN
        PBC_DEFECT = X + RBOX
      ELSE
        PBC_DEFECT = X - RBOX
      ENDIF

      RETURN

      END FUNCTION
!***********************************************************************!
!!IDENTIFY THE DEFECT/CENTER OF AV IN THE SYSTEM
!!THIS ALGORITHM IS ONLY VALID FOR BCC SYSTEM,E.G. IRON
!!A GENERALIZED VERSION IS ALSO AVALABLE

      SUBROUTINE FIND_DEFECT

      USE SYSTEM

      USE DEFECT

      IMPLICIT NONE

      INTEGER I

      NDINT = 0

      NDVAC = 0

      TI = 0

      TV = 0

!!SCALE DOWN POSITIONS TO CUBE
      DO I=1,NATOMC
        RA(1,I) = RA(1,I)/EPSX
        RA(2,I) = RA(2,I)/EPSY
        RA(3,I) = RA(3,I)/EPSZ
      END DO

      IF (TCRY .EQ. 1) THEN
        CALL SINGLE_CRYSTAL_WS
      ELSE IF (TCRY .EQ. 2) THEN
        CALL POLY_CRYSTAL_WS
      END IF
      
!!TALLY NUMBER OF VACANCIES AND INTERSTITIALS

      DO I=1,NATOMC

        IF(TI(I).NE.1) THEN

          NDINT = NDINT + 1

        ENDIF

      ENDDO

      DO I=1,NATOMR

        IF(TV(I).NE.1) THEN

          NDVAC = NDVAC + 1

        ENDIF

      ENDDO

!!RESCALE ATOM POSITIONS TO ORIGINAL VALUES
      DO I=1,NATOMC
        RA(1,I) = RA(1,I)*EPSX
        RA(2,I) = RA(2,I)*EPSY
        RA(3,I) = RA(3,I)*EPSZ
      END DO

      RETURN
      END SUBROUTINE FIND_DEFECT
!***********************************************************************!
      SUBROUTINE SINGLE_CRYSTAL_WS

      USE SYSTEM

      USE DEFECT

      IMPLICIT NONE

      REAL*8 RX,RY,RZ,RS,PBC_DEFECT

      INTEGER I,J,NBOX

      INTEGER IX,IY,IZ,JX,JY,JZ,RX1,RY1,RZ1,RX2,RY2,RZ2

      NBOX = INT(RBOX)*2

      DO I=1,NATOMC

        IF(RA(1,I).LT.0.0) RA(1,I)=RA(1,I)+RBOX

        IF(RA(2,I).LT.0.0) RA(2,I)=RA(2,I)+RBOX

        IF(RA(3,I).LT.0.0) RA(3,I)=RA(3,I)+RBOX

        IF(RA(1,I).GE.RBOX) RA(1,I)=RA(1,I)-RBOX

        IF(RA(2,I).GE.RBOX) RA(2,I)=RA(2,I)-RBOX

        IF(RA(3,I).GE.RBOX) RA(3,I)=RA(3,I)-RBOX

        IX = FLOOR(RA(1,I)/0.5D0)

        IY = FLOOR(RA(2,I)/0.5D0)

        IZ = FLOOR(RA(3,I)/0.5D0)

        IF(IX.GE.NBOX) IX = IX - NBOX

        IF(IY.GE.NBOX) IY = IY - NBOX

        IF(IZ.GE.NBOX) IZ = IZ - NBOX

        JX = IX + 1

        JY = IY + 1

        JZ = IZ + 1

        IF(JX.GE.NBOX) JX = JX - NBOX

        IF(JY.GE.NBOX) JY = JY - NBOX

        IF(JZ.GE.NBOX) JZ = JZ - NBOX

        IF(MOD(IX,2).EQ.1) THEN

          RX1 = IX

          RX2 = JX

        ELSE

          RX1 = JX

          RX2 = IX

        ENDIF

        IF(MOD(IY,2).EQ.1) THEN

          RY1 = IY

          RY2 = JY

        ELSE

          RY1 = JY

          RY2 = IY

        ENDIF

        IF(MOD(IZ,2).EQ.1) THEN

          RZ1 = IZ

          RZ2 = JZ

        ELSE

          RZ1 = JZ

          RZ2 = IZ

        ENDIF

        RX=PBC_DEFECT(DFLOAT(RX1)*0.5D0-RA(1,I))

        RY=PBC_DEFECT(DFLOAT(RY1)*0.5D0-RA(2,I))

        RZ=PBC_DEFECT(DFLOAT(RZ1)*0.5D0-RA(3,I))

        RS=RX*RX+RY*RY+RZ*RZ

        IF(RS.LE.DECUT) THEN

            J = FLOOR(DFLOAT(RX1)*0.5D0)*2*RBOX**2
     #        + FLOOR(DFLOAT(RY1)*0.5D0)*2*RBOX
     #        + RZ1 + 1

            TI(I) = TI(I) + 1

            TV(J) = TV(J) + 1

        ENDIF

        RX=PBC_DEFECT(DFLOAT(RX2)*0.5D0-RA(1,I))

        RY=PBC_DEFECT(DFLOAT(RY2)*0.5D0-RA(2,I))

        RZ=PBC_DEFECT(DFLOAT(RZ2)*0.5D0-RA(3,I))

        RS=RX*RX+RY*RY+RZ*RZ

        IF(RS.LE.DECUT) THEN

            J = FLOOR(DFLOAT(RX2)*0.5D0)*2*RBOX**2
     #        + FLOOR(DFLOAT(RY2)*0.5D0)*2*RBOX
     #        + RZ2 + 1

            TI(I) = TI(I) + 1

            TV(J) = TV(J) + 1

        ENDIF

      ENDDO

      END SUBROUTINE SINGLE_CRYSTAL_WS
!***********************************************************************!
      SUBROUTINE POLY_CRYSTAL_WS

      USE OMP_LIB
      USE SYSTEM
      USE DEFECT

      IMPLICIT NONE

      REAL*8 RX,RY,RZ,RS,RMIN,PBC_DEFECT
      INTEGER I,J,MIN,N,ATOM

      OC = 0
      OCC= 0

!!DETERMINE THE OCCUPANCY OF EACH SITE

!$OMP PARALLEL PRIVATE(RS,RX,RY,RZ,RMIN,MIN,J,N,ATOM) SHARED(OC,OCC)
!$OMP DO
      DO I = 1, NATOMC

        RS=0;RMIN=HUGE(RMIN);MIN=0
        DO J = 1, NATOMR
          RX = PBC_DEFECT(RA(1,I) - RR(1,J))
          RY = PBC_DEFECT(RA(2,I) - RR(2,J))
          RZ = PBC_DEFECT(RA(3,I) - RR(3,J))
          RS = RX*RX + RY*RY + RZ*RZ
          IF (RS .LT. RMIN) THEN
            MIN = J
            RMIN = RS
          END IF
        END DO

!$OMP CRITICAL
        OCC(MIN) = OCC(MIN) + 1
        OC(OCC(MIN),MIN) = I
!$OMP END CRITICAL

      END DO
!$OMP END DO

!!USING OCCUPANCY TO DETERMINE DEFECTS
!$OMP DO
      DO I=1,NATOMR

        N = OCC(I)
        IF (N .GT. 0) THEN
          TV(I) = N
          DO J=1,N
            ATOM = OC(J,I)
            TI(ATOM) = N
          END DO
        END IF

      END DO
!$OMP END DO
!$OMP END PARALLEL

      END SUBROUTINE POLY_CRYSTAL_WS
!***********************************************************************!
!!OUTPUT DEFECTS FOR OVITO

      SUBROUTINE OVITO

      USE SYSTEM

      USE DEFECT

      IMPLICIT NONE

      INTEGER I, J

      CHARACTER(len = 200) FILENAME

!!CREATING THE CURRENT DEFECT POSITION FILE

      WRITE(FILENAME, '(A8, I0, A4)') './DUMP/D', ISTEP, '.OUT'

      OPEN(UNIT=29,FILE=TRIM(FILENAME))

      WRITE(29,'(A)') 'ITEM: TIMESTEP'
      WRITE(29,*) ISTEP
      WRITE(29,'(A)') 'ITEM: NUMBER OF ATOMS'
      WRITE(29,*) NDVAC+NDINT
      WRITE(29,'(A)') 'ITEM: BOX BOUNDS PP PP PP'
      WRITE(29,'(2F10.4,I5)') 0.00,RBOX*EPSX,3
      WRITE(29,'(2F10.4,I5)') 0.00,RBOX*EPSY,3
      WRITE(29,'(2F10.4,I5)') 0.00,RBOX*EPSZ,3
      WRITE(29,'(A)') 'ITEM: ATOMS x y z type'

      DO I=1,NATOMC

        IF(TI(I).NE.1) THEN

          WRITE(29,'(3F12.5,I5)') RA(:,I),2

        ENDIF

      ENDDO

      DO I=1,NATOMR

        IF(TV(I).NE.1) THEN

          WRITE(29,'(3F12.5,I5)') RR(:,I),1

        ENDIF

      ENDDO

      CLOSE(29)

!!CREATING DISLOCATION ANALYSIS FILE

      WRITE(FILENAME, '(A12, I0, A4)') './DUMP/DXA/A', ISTEP, '.OUT'

      OPEN(UNIT=30,FILE=TRIM(FILENAME))

      WRITE(30,'(A)') 'ITEM: TIMESTEP'
      WRITE(30,*) ISTEP
      WRITE(30,'(A)') 'ITEM: NUMBER OF ATOMS'
      WRITE(30,*) NATOMC
      WRITE(30,'(A)') 'ITEM: BOX BOUNDS PP PP PP'
      WRITE(30,'(2F10.4,I5)') 0.00,RBOX*EPSX,3
      WRITE(30,'(2F10.4,I5)') 0.00,RBOX*EPSY,3
      WRITE(30,'(2F10.4,I5)') 0.00,RBOX*EPSZ,3
      WRITE(30,'(A)') 'ITEM: ATOMS x y z'

      DO I = 1, NATOMC
        WRITE(30,'(3F14.5)') (RA(J,I), J=1,3)
      END DO

      CLOSE(30)

      END SUBROUTINE OVITO
!***********************************************************************!
      SUBROUTINE BASIC_DAMAGE(NFRKL)

      USE SYSTEM, ONLY: RA,NATOMC,RANK,DECUT,TDAM

      IMPLICIT NONE

      INTEGER :: I, J, ID, NFRKL, N
      INTEGER :: MOVEID(NFRKL)
      REAL *8 :: DLIM
      REAL *8 :: RND, NEWPOS(3)
      LOGICAL :: INSERT

      DLIM = DECUT
      
      N = 0

      MOVEID = 0

      CALL RANDOM_SEED()

      !!PICK NFRKL ATOMS TO MOVE (SHOULD CREATE NFKRL # OF FRENKEL PAIRS)

      IF (RANK .EQ. 0) THEN

        N = 1

        DO WHILE (N .LE. NFRKL)

          CALL RANDOM_NUMBER(RND)
          ID = CEILING(RND*NATOMC)
          IF (ID .EQ. 0) ID = 1

          INSERT = .TRUE.
          DO J = 1,N-1
            IF (ID .EQ. MOVEID(J)) INSERT = .FALSE.
          END DO
          IF (INSERT) THEN
            MOVEID(N) = ID
            N = N + 1
          END IF

        END DO

      END IF

      !!MOVE CHOSEN ATOMS TO RANDOM LOCATIONS

      IF (TDAM .EQ. 2) THEN

        IF (RANK .EQ. 0) THEN
        
          WRITE(*,'(A,I0,A)') 'Generating ',NFRKL,' Frenkel pairs'

          DO I = 1, NFRKL

            CALL FRENKEL_POS_UPDATE(NEWPOS,DLIM)
            RA(:,MOVEID(I)) = NEWPOS

          END DO

        END IF

      ELSE IF (TDAM .EQ. 3) THEN

        IF (RANK .EQ.0) WRITE(*,'(A,I0,A)') 'Generating ',NFRKL,' single vacancies'

        CALL CREATE_VACANCIES(NFRKL, MOVEID)
        
      END IF

      END SUBROUTINE BASIC_DAMAGE
!***********************************************************************!
      SUBROUTINE CREATE_VACANCIES(NFRKL,IDS)

      USE SYSTEM
      USE LMP_DRIVER
      USE MPI

      IMPLICIT NONE

      INTEGER :: I,NFRKL
      INTEGER :: IDS(NFRKL)
      CHARACTER(LEN=10) :: IDSTR

      CALL MPI_BCAST(IDS, NFRKL, MPI_INT, 0, MPI_COMM_WORLD, IERR)
      
      !Update lammps on what atoms were deleted
      DO I = 1, NFRKL
        WRITE(IDSTR,'(I0)') IDS(I)
        CALL lammps_command(lmp, 'group deleteList id '//TRIM(IDSTR))
      END DO
      CALL lammps_command(lmp, 'delete_atoms group deleteList')

      !Upadate all other processes on changed parameters
      NATOMC = int(lammps_get_natoms(lmp))
      CALL RECV_LMP_POS

      END SUBROUTINE CREATE_VACANCIES
!***********************************************************************!
      SUBROUTINE FRENKEL_POS_UPDATE(POS,THRESH)

      USE SYSTEM, ONLY: RA, NATOMC, RBOX, EPSX, EPSY, EPSZ

      IMPLICIT NONE

      INTEGER :: I, ID_NEAR
      REAL *8, INTENT(INOUT) :: POS(3)
      REAL *8, INTENT(IN) :: THRESH
      REAL *8 :: RX, RY, RZ, RS, RS_NEAR
      REAL *8 :: UNIT(3)
      LOGICAL :: FLAG

      ID_NEAR = 0; RS_NEAR = HUGE(RS_NEAR); FLAG = .TRUE.

      CALL RANDOM_NUMBER(POS)

      POS(1) = POS(1)*RBOX*EPSX
      POS(2) = POS(2)*RBOX*EPSY
      POS(3) = POS(3)*RBOX*EPSZ 

      !!FIND CLOSEST ATOM TO (X,Y,Z)

      DO I = 1, NATOMC

        RX = RA(1,I) - POS(1)

        RY = RA(2,I) - POS(2)

        RZ = RA(3,I) - POS(3)

        RS = RX*RX + RY*RY + RZ*RZ

        IF (RS .LT. RS_NEAR) THEN

          RS_NEAR = RS

          ID_NEAR = I

        END IF

      END DO

      !!! UPDATE THE POSITION UNTIL IT IS SUITABLE FOR RELAXATION

      DO WHILE (FLAG)

        FLAG = .FALSE.

        UNIT(1) = POS(1) - RA(1,ID_NEAR) 

        UNIT(2) = POS(2) - RA(2,ID_NEAR) 

        UNIT(3) = POS(3) - RA(3,ID_NEAR)

        UNIT = UNIT/NORM2(UNIT)

        POS = POS + SQRT(THRESH)/3*UNIT

        RX = POS(1) - RA(1,ID_NEAR)

        RY = POS(2) - RA(2,ID_NEAR)

        RZ = POS(3) - RA(3,ID_NEAR)

        RS = RX*RX + RY*RY + RZ*RZ

        IF (RS .LT. THRESH) FLAG = .TRUE.

      END DO

      END SUBROUTINE FRENKEL_POS_UPDATE
!***********************************************************************!
      SUBROUTINE MD_DPA

      USE SYSTEM, ONLY: NATOMC, DPA, AVDPA, DPAC, DTIME, DRATE, TDAM,CALC
      USE DEFECT, ONLY: ONDSTOT, NDSTOT

      IMPLICIT NONE

      REAL *8 :: NEWDSP

      !UPDATE AVERAGE DPA VALUE

      NEWDSP = NDSTOT - ONDSTOT

      AVDPA = (AVDPA*DPAC + NEWDSP/NATOMC) / (DPAC + 1)

      DPAC = DPAC + 1

      !UPDATE TIME BETWEEN CASCADES

      DTIME = AVDPA / DRATE
      IF (TDAM .EQ. 3) DTIME = HUGE(DTIME)

      !ADD AVERAGE DPA PER CASCADE TO CURRENT DPA VALUE

      DPA = DPA + NEWDSP/NATOMC

      CALC = .FALSE.

      END SUBROUTINE MD_DPA
!***********************************************************************!
      !!CREATES RESTART FILES AT FIXED INTERVALS. THESE RESTART FILES
      !!INCLUDE THE SYSTEM CONFIGURATION INCASE ONE WANT TO RESTART THE 
      !!SIMULATION AT A CERTAIN POINT FOR ANY REASON.
      SUBROUTINE RESTART
      
      USE SYSTEM
      USE DEFECT

      IMPLICIT NONE

      INTEGER I, J

      CHARACTER(len = 200) FILENAME

      IF (RANK .EQ. 0) THEN

      IF (MOD(ISTEP, 1000) .EQ. 0) THEN

        WRITE(FILENAME, '(A17, I0, A4)') './RESTART/RESTART', 
     &  ISTEP, '.OUT'

        OPEN(UNIT=27,FILE=TRIM(FILENAME))

        WRITE(27,*) NATOMC

        DO I = 1, NATOMC

          WRITE(27,'(3F14.5)') (RA(J,I), J=1,3)

        ENDDO

        !SIMULATION MEASUREMENTS AND PARAMETERS

        WRITE(27,*)

        WRITE(27,*) 
     #    'ISTEP, TIME, TIME_LAST, DTIME, DPA, DPAC, AVDPA, ONDSTOT'

        WRITE(27,*) ISTEP
        WRITE(27,*) TIME
        WRITE(27,*) TIME_LAST
        WRITE(27,*) DTIME
        WRITE(27,*) DPA
        WRITE(27,*) DPAC
        WRITE(27,*) AVDPA
        WRITE(27,*) ONDSTOT
        WRITE(27,*) EPSX
        WRITE(27,*) EPSY
        WRITE(27,*) EPSZ

        !AVERAGE DEFECT MIGRATION ENERGIES

        WRITE(27,*)

        WRITE(27,*) 'MIGRATION ENERGIES (VACANCIES THEN INTERSTITIALS)'

        DO I = 1, 50
          WRITE(27,*) EMIG(1,I), EMIC(1,I)
        END DO

        DO I = 1, 50
          WRITE(27,*) EMIG(2,I), EMIC(2,I)
        END DO

      END IF

      CLOSE(27)

      END IF

      END SUBROUTINE

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
