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

      REAL*8 FUNCTION PBC(X)

      USE SYSTEM,ONLY : RBOX,RBH1,RBH2

      IMPLICIT NONE

      REAL*8 X

!      PBC = X-RBOX*NINT(X/RBOX)

      IF((X.LE.RBH1).AND.(X.GE.RBH2)) THEN

        PBC = X

      ELSEIF(X.LT.RBH2) THEN

        PBC = X + RBOX

      ELSE

        PBC = X - RBOX

      ENDIF

      RETURN

      END
!***********************************************************************!
!!IDENTIFY THE DEFECT/CENTER OF AV IN THE SYSTEM
!!THIS ALGORITHM IS ONLY VALID FOR BCC SYSTEM,E.G. IRON
!!A GENERALIZED VERSION IS ALSO AVALABLE

      SUBROUTINE FIND_DEFECT

      USE SYSTEM

      USE DEFECT

      IMPLICIT NONE

      REAL*8 RX,RY,RZ,RD,RS,PBC

      INTEGER I,J,NBOX

      INTEGER IX,IY,IZ,JX,JY,JZ,RX1,RY1,RZ1,RX2,RY2,RZ2

      NDINT = 0

      NDVAC = 0

      TI = 0

      TV = 0

      NBOX = INT(RBOX)*2

!!! COULD USE OPENMP

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

        RX=PBC(DFLOAT(RX1)*0.5D0-RA(1,I))

        RY=PBC(DFLOAT(RY1)*0.5D0-RA(2,I))

        RZ=PBC(DFLOAT(RZ1)*0.5D0-RA(3,I))

        RS=RX*RX+RY*RY+RZ*RZ

        IF(RS.LE.DECUT) THEN

            J = FLOOR(DFLOAT(RX1)*0.5D0)*2*RBOX**2
     #        + FLOOR(DFLOAT(RY1)*0.5D0)*2*RBOX
     #        + RZ1 + 1

            TI(I) = TI(I) + 1

            TV(J) = TV(J) + 1

        ENDIF

        RX=PBC(DFLOAT(RX2)*0.5D0-RA(1,I))

        RY=PBC(DFLOAT(RY2)*0.5D0-RA(2,I))

        RZ=PBC(DFLOAT(RZ2)*0.5D0-RA(3,I))

        RS=RX*RX+RY*RY+RZ*RZ

        IF(RS.LE.DECUT) THEN

            J = FLOOR(DFLOAT(RX2)*0.5D0)*2*RBOX**2
     #        + FLOOR(DFLOAT(RY2)*0.5D0)*2*RBOX
     #        + RZ2 + 1

            TI(I) = TI(I) + 1

            TV(J) = TV(J) + 1

        ENDIF

      ENDDO
      
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

      RETURN
      END SUBROUTINE FIND_DEFECT
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
      WRITE(29,'(2F10.4,I5)') 0.00,RBOX,3
      WRITE(29,'(2F10.4,I5)') 0.00,RBOX,3
      WRITE(29,'(2F10.4,I5)') 0.00,RBOX,3
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
      WRITE(30,'(2F10.4,I5)') 0.00,RBOX,3
      WRITE(30,'(2F10.4,I5)') 0.00,RBOX,3
      WRITE(30,'(2F10.4,I5)') 0.00,RBOX,3
      WRITE(30,'(A)') 'ITEM: ATOMS x y z'

      DO I = 1, NATOMC

        WRITE(30,'(3F14.5)') (RA(J,I), J=1,3)

      END DO

      CLOSE(30)

      END SUBROUTINE OVITO
!***********************************************************************!
      SUBROUTINE ELECTRON_DAMAGE(NFRKL)

      USE SYSTEM, ONLY: RA,NATOMC,RBOX,RANK,DECUT
      USE DEFECT, ONLY: IV,II,NDINT,NDVAC

      IMPLICIT NONE

      INTEGER :: I, J, ID, NFRKL, N
      INTEGER :: MOVEID(NFRKL)
      REAL *8 :: RX, RY, RZ, RS, PBC, DLIM
      REAL *8 :: RND, NEWPOS(3)
      LOGICAL :: FLAG

      DLIM = DECUT
      
      N = 0

      MOVEID = 0

      CALL RANDOM_SEED()

      WRITE(*,'(A,I0,A)') 'Generating ',NFRKL,' Frenkel pairs'

      !!PICK NFRKL ATOMS TO MOVE (SHOULD CREATE NFKRL # OF FRENKEL PAIRS)

      DO I = 1, NFRKL

        CALL RANDOM_NUMBER(RND)

        ID = CEILING(RND*NATOMC)

        IF (ID .EQ. 0) ID = 1

        MOVEID(I) = ID

      END DO

      !!MOVE CHOSEN ATOMS TO RANDOM LOCATIONS

      DO I = 1, NFRKL

        CALL FRENKEL_POS_UPDATE(NEWPOS,DLIM)

        RA(:,MOVEID(I)) = NEWPOS

      END DO

      END SUBROUTINE ELECTRON_DAMAGE
!***********************************************************************!
      SUBROUTINE FRENKEL_POS_UPDATE(POS,THRESH)

      USE SYSTEM, ONLY: RA, NATOMC, RBOX

      IMPLICIT NONE

      INTEGER :: I, ID_NEAR
      REAL *8, INTENT(INOUT) :: POS(3)
      REAL *8, INTENT(IN) :: THRESH
      REAL *8 :: RX, RY, RZ, RS, RS_NEAR
      REAL *8 :: UNIT(3)
      LOGICAL :: FLAG

      ID_NEAR = 0; RS_NEAR = HUGE(RS_NEAR); FLAG = .TRUE.

      CALL RANDOM_NUMBER(POS)

      POS = POS * RBOX

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

      USE SYSTEM, ONLY: NATOMC, DPA, AVDPA, DPAC, DTIME, DRATE, TDAM,
     #  CALC
      USE DEFECT, ONLY: ONDSTOT, NDSTOT

      IMPLICIT NONE

      REAL *8 :: NEWDSP

      !UPDATE AVERAGE DPA VALUE

      NEWDSP = NDSTOT - ONDSTOT

      AVDPA = (AVDPA*DPAC + NEWDSP/NATOMC) / (DPAC + 1)

      DPAC = DPAC + 1

      !UPDATE TIME BETWEEN CASCADES

      DTIME = AVDPA / DRATE

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
