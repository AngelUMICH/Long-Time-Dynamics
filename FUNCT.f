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

!!COMPARE THE CURRENT LATTICE WITH THE REFERENCE LATTICE

C      DO J=1,NIONSREF
C        DO I = 1,NIONS
C          RX=PBC(POSION(1,I)-POSREF(1,J))
C          IF(ABS(RX).GT.RCUTDEF) GOTO 100
C          RY=PBC(POSION(2,I)-POSREF(2,J))
C          IF(ABS(RY).GT.RCUTDEF) GOTO 100
C          RZ=PBC(POSION(3,I)-POSREF(3,J))
C          IF(ABS(RZ).GT.RCUTDEF) GOTO 100
C          RD=RX*RX+RY*RY+RZ*RZ
C          IF(RD.LE.RCUTDEFS) THEN
C            IFLAG(I)=IFLAG(I)+1
C            JFLAG(J)=JFLAG(J)+1
C          ENDIF
C100       CONTINUE
C        ENDDO
C      ENDDO

      RETURN
      END SUBROUTINE FIND_DEFECT
!***********************************************************************!
!***********************************************************************!
!!OUTPUT THE DEFECT STRUCTURE TO FILE THAT CAN BE SEEN BY ATOMTV

      SUBROUTINE VISUAL

      USE SYSTEM,ONLY : NFRAME

      IMPLICIT NONE

      CALL FIND_DEFECT

      CALL ATOMTV

      NFRAME = NFRAME + 1

      RETURN
      END SUBROUTINE VISUAL
!***********************************************************************!
!!OUTPUT THE DEFECT STUCTURE USING ATOMTV

      SUBROUTINE ATOMTV

      USE SYSTEM

      USE DEFECT

      IMPLICIT NONE

      INTEGER I,NBOND

      NBOND = 0

      NDVAC = 0

      NDINT = 0

!!! COULD USE OPENMP

!!COUNTING THE NUMBER OF INTERSTITIALS

      DO I=1,NATOMC

        IF(TI(I).NE.1) THEN

          NDINT = NDINT + 1

        ENDIF

      ENDDO

!!COUNTING THE NUMBER OF VACANCIES

      DO I=1,NATOMR

        IF(TV(I).NE.1) THEN

          NDVAC = NDVAC + 1

        ENDIF

      ENDDO

      WRITE(22,'(2I8,A5,E14.5,A1)') NDVAC+NDINT+8,NBOND,'"',TIME,'"'

      WRITE(22,'(3F10.4,2I5)') 0.00,0.00,0.00, 3,3
      WRITE(22,'(3F10.4,2I5)') 0.00,0.00,RBOX, 3,3
      WRITE(22,'(3F10.4,2I5)') 0.00,RBOX,0.00, 3,3
      WRITE(22,'(3F10.4,2I5)') RBOX,0.00,0.00, 3,3
      WRITE(22,'(3F10.4,2I5)') 0.00,RBOX,RBOX, 3,3
      WRITE(22,'(3F10.4,2I5)') RBOX,0.00,RBOX, 3,3
      WRITE(22,'(3F10.4,2I5)') RBOX,RBOX,0.00, 3,3
      WRITE(22,'(3F10.4,2I5)') RBOX,RBOX,RBOX, 3,3

      DO I=1,NATOMC

        IF(TI(I).NE.1) THEN

          WRITE(22,'(3F12.5,2I5)') RA(:,I),1,1

        ENDIF

      ENDDO

      DO I=1,NATOMR

        IF(TV(I).NE.1) THEN

          WRITE(22,'(3F12.5,2I5)') RR(:,I),2,2

        ENDIF

      ENDDO

      RETURN
      END SUBROUTINE ATOMTV
!***********************************************************************!
!!OUTPUT DEFECTS FOR OVITO

      SUBROUTINE OVITO

      USE SYSTEM

      USE DEFECT

      IMPLICIT NONE

      INTEGER I

      CHARACTER(len = 200) FILENAME, STR

      IF (ISTEP .EQ. 0) THEN

        CALL EXECUTE_COMMAND_LINE('mkdir DUMP')

      ENDIF

!!CREATING THE CURRENT DEFECT POSITION FILE

      WRITE(FILENAME, '(A8, I5.5, A4)') './DUMP/D', ISTEP, '.OUT'

      OPEN(UNIT=29,FILE=TRIM(FILENAME))

      WRITE(29,'(A)') 'ITEM: TIMESTEP'
      WRITE(29,*) ISTEP+1
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

      END SUBROUTINE OVITO
!***********************************************************************!
      !!CREATES THE SAME AMOUNT OF VACANCIES AND INTERSTITIALS WHICH ARE ADDED TO
      !!SYSTEM THAT MAY OR MAY NOT HAVE DEFECTS INSIDE ALREADY. DUE TO FITTING
      !!ISSUES OF DEFECTS THE DESIRED AMOUNT OF DEFECTS MAY NOT BE ACHIEVED BUT
      !!FRENEKEL PAIRS ARE CREATED.
      SUBROUTINE FRENKELGENERATOR(NDFRKL)

      USE SYSTEM,ONLY: RA,RR,NATOMC,RBOX
      USE DEFECT,ONLY: IV,II,NDINT,NDVAC

      IMPLICIT NONE

      INTEGER NDFRKL, ADD_VACS, ADD_INTS
      INTEGER I, J, TRIES, ICLOSE
      INTEGER :: MOVEI(NDFRKL)
      LOGICAL TOO_CLOSE
      REAL*8 RX, RY, RZ, RS, RS_LOWEST, RND(3)
      REAL*8 ::XVAC(NDFRKL), YVAC(NDFRKL), ZVAC(NDFRKL),
     &         XINT(NDFRKL), YINT(NDFRKL), ZINT(NDFRKL)

      !!NDFRKL::NUMBER OF FRENKEL PAIRS TO CREATE
      !!NDINT:: THE NUMBER OF CURRENT INTERSTITIALS
      !!NDVAC:: NUMBER OF CURRENT VACANCIES
      !!RBOX:: SIMULATION SIZE IN LATTICE UNITS
      !!II::INDEX OF INTERSTITIALS
      !!IV::INDEX OF VACANCIES
      !!MOVEI::INDEXES OF ATOMS TO RELOCATE TO CREATE NEW INTERSITITALS
      !!ADD_VACS::NUMBER OF NEW VACANCIES TO BE CREATED

      ADD_VACS = 0
      ADD_INTS = 0

      CALL RANDOM_SEED()

      !!CREATING THE VACANCIES
      DO I = 1, NDFRKL
        TRIES = 0
        RS = 0
        TOO_CLOSE = .TRUE.
        DO WHILE((TRIES .LT. 100) .AND. (TOO_CLOSE .EQV. .TRUE.))
          CALL RANDOM_NUMBER(RND)
          XVAC(I) = INT(RND(1)*RBOX/0.5D0)*0.5D0
          !CORNER ATOM POSITIONS WHERE COORDINATES CAN BE INTEGERS ONLY
          IF (XVAC(I) - INT(XVAC(I)) .LT. 1D-2) THEN
            YVAC(I) = INT(RND(2)*RBOX)
            ZVAC(I) = INT(RND(3)*RBOX)
          !CENTER ATOM POSITIONS (I.5, J.5, K.5) WHERE I J OR K CAN BE ANY INTEGER
          ELSE IF (XVAC(I) - INT(XVAC(I)) .GE. 0.5D0) THEN
            YVAC(I) = INT(RND(2)*RBOX) + 0.5D0
            ZVAC(I) = INT(RND(3)*RBOX) + 0.5D0
          ENDIF

          TOO_CLOSE = .FALSE.
          !!CHECKING IF SELECTED VACANCY LOCATION IS TOO CLOSE TO CURRENT INTERSITITALS
          DO J = 1, NDINT
            RX = XVAC(I) - RA(1,II(1,J))
            RY = YVAC(I) - RA(2,II(1,J))
            RZ = ZVAC(I) - RA(3,II(1,J))
            RS = RX**2 + RY**2 + RZ**2
            IF (RS .LT. 2.0**2) THEN
              TOO_CLOSE = .TRUE.
            ENDIF
          ENDDO
          !!CHECKING IF SELECTED VACANCY LOCATION IS TOO CLOSE TO CURRENT VACANCIES
          DO J = 1, NDVAC
            RX = XVAC(I) - RR(1,IV(1,J))
            RY = YVAC(I) - RR(2,IV(1,J))
            RZ = ZVAC(I) - RR(3,IV(1,J))
            RS = RX**2 + RY**2 + RZ**2
            IF (RS .LT. 2.0**2) THEN
              TOO_CLOSE = .TRUE.
            ENDIF
          ENDDO
          !!CHECKING IF SELECTED VACANCY IS TOO CLOSE TO PREVIOUS NEW VACANCIES
          IF (I .NE. 1) THEN
            DO J = 1, I-1
              RX = XVAC(I) - XVAC(J)
              RY = YVAC(I) - YVAC(J)
              RZ = ZVAC(I) - ZVAC(J)
              RS = RX**2 + RY**2 + RZ**2
              IF (RS .LT. 2.0**2) THEN
                TOO_CLOSE = .TRUE.
              ENDIF
            ENDDO
          ENDIF
          TRIES = TRIES + 1
        ENDDO
        IF (TRIES .LT. 100) THEN
          ADD_VACS = ADD_VACS + 1
        ENDIF
      ENDDO

      !WRITE(*,'(A11,I5)') '# OF VACS= ',ADD_VACS
      !WRITE(*,'(A7,I8)') 'TRIES= ',TRIES

      !!SHOULD HAVE BY NOW LOCATION OF NEW VANCANIES
      !!DETERMINE THE INDEX IN RA ARRRAY
      !!THE INDEXES FOUND WILL SERVE AS THE INDEXES OF ATOMS TO MOVE
      !!TO CREATE NEW INTERSITITALS

      !!ICLOSE :: INDEX OF THE CLOSEST ATOM

      DO I = 1, ADD_VACS
        ICLOSE = 1
        DO J = 1, NATOMC
          RX = XVAC(I) - RA(1,J)
          RY = YVAC(I) - RA(2,J)
          RZ = ZVAC(I) - RA(3,J)
          RS = RX**2 + RY**2 + RZ**2

          RX = XVAC(I) - RA(1,ICLOSE)
          RY = YVAC(I) - RA(2,ICLOSE)
          RZ = ZVAC(I) - RA(3,ICLOSE)
          RS_LOWEST = RX**2 + RY**2 + RZ**2
          IF (RS .LT. RS_LOWEST) THEN
            ICLOSE = J
          ENDIF
        ENDDO
        MOVEI(I) = ICLOSE
      ENDDO

      !!CREATING INTERSTITIALS ANYWHERE IN THE SYSTEM EXCLUDING CLOSE TO EDGES
      DO I = 1, ADD_VACS  !USING ADD_VACS INSTEAD OF NDFRKL TO AVOID EXTRA WORK
        TRIES = 0
        RS = 0
        TOO_CLOSE = .TRUE.
        DO WHILE((TRIES .LT. 100) .AND. (TOO_CLOSE .EQV. .TRUE.))
          CALL RANDOM_NUMBER(RND)
          XINT(I) = RND(1)*(RBOX-4) + 2
          YINT(I) = RND(2)*(RBOX-4) + 2
          ZINT(I) = RND(3)*(RBOX-4) + 2

          TOO_CLOSE = .FALSE.
          !!CHECKING IF SELECTED INT LOCATION IS TOO CLOSE TO CURRENT INTERSITITALS
          DO J = 1, NDINT
            RX = XINT(I) - RA(1,II(1,J))
            RY = YINT(I) - RA(2,II(1,J))
            RZ = ZINT(I) - RA(3,II(1,J))
            RS = RX**2 + RY**2 + RZ**2
            IF (RS .LT. 2.0**2) THEN
              TOO_CLOSE = .TRUE.
            ENDIF
          ENDDO
          !!CHECKING IF SELECTED INT LOCATION IS TOO CLOSE TO CURRENT VACANCIES
          DO J = 1, NDVAC
            RX = XINT(I) - RR(1,IV(1,J))
            RY = YINT(I) - RR(2,IV(1,J))
            RZ = ZINT(I) - RR(3,IV(1,J))
            RS = RX**2 + RY**2 + RZ**2
            IF (RS .LT. 2.0**2) THEN
              TOO_CLOSE = .TRUE.
            ENDIF
          ENDDO
          !!CHECKING IF SELECTED INT LOCATION IS TOO CLOSE TO NEW VACANCIES
          DO J = 1, ADD_VACS
            RX = XINT(I) - XVAC(J)
            RY = YINT(I) - YVAC(J)
            RZ = ZINT(I) - ZVAC(J)
            RS = RX**2 + RY**2 + RZ**2
            IF (RS .LT. 2.0**2) THEN
              TOO_CLOSE = .TRUE.
            ENDIF
          ENDDO
          !!CHECKING IF SELECTED INT LOCATION IS TOO CLOSE TO PREVIOUS NEW INTS
          IF (I .NE. 1) THEN
            DO J = 1, I-1
              RX = XINT(I) - XINT(J)
              RY = YINT(I) - YINT(J)
              RZ = ZINT(I) - ZINT(J)
              RS = RX**2 + RY**2 + RZ**2
              IF (RS .LT. 2.0**2) THEN
                TOO_CLOSE = .TRUE.
              ENDIF
            ENDDO
          ENDIF
          TRIES = TRIES + 1
        ENDDO
        IF (TRIES .LT. 100) THEN
          ADD_INTS = ADD_INTS + 1
        ENDIF
      ENDDO

      !!NOW I SHOULD HAVE LOCATIONS FOR NEW INTS AND VACS
      !!NOW ALTER INT POSITIONS FOR BETTER RELAXATION

      IF (ADD_INTS .NE. 0) THEN

        !!MOVE THE CLOSEST ATOM TO THE INT AND THE INT A SLIGHT BIT FOR BETTER RELAX
        DO I = 1, ADD_INTS
          ICLOSE = 1
          DO J = 1, NATOMC
            RX = XINT(I) - RA(1,J)
            RY = YINT(I) - RA(2,J)
            RZ = ZINT(I) - RA(3,J)
            RS = RX**2 + RY**2 + RZ**2

            RX = XINT(I) - RA(1,ICLOSE)
            RY = YINT(I) - RA(2,ICLOSE)
            RZ = ZINT(I) - RA(3,ICLOSE)
            RS_LOWEST = RX**2 + RY**2 + RZ**2
            IF (RS .LT. RS_LOWEST) THEN
              ICLOSE = J
            ENDIF
          ENDDO

          !!CALCULATE THE UNIT VECTOR FOR MOVE STEP
          RX = RA(1,ICLOSE) - XINT(I)
          RY = RA(2,ICLOSE) - YINT(I)
          RZ = RA(3,ICLOSE) - ZINT(I)
          RS_LOWEST = RX**2 + RY**2 + RZ**2
          RX = RX / SQRT(RS_LOWEST)
          RY = RY / SQRT(RS_LOWEST)
          RZ = RZ / SQRT(RS_LOWEST)

          !!THIS ACTUALLY CREATES THE FRENKEL DEFECTS (#VACS CREATED = #INTS CREATED)
          RA(1,MOVEI(I)) = RA(1,ICLOSE) - RX*0.25
          RA(2,MOVEI(I)) = RA(2,ICLOSE) - RY*0.25
          RA(3,MOVEI(I)) = RA(3,ICLOSE) - RZ*0.25

          !!MOVE BOTH THE CLOSEST ATOM AND THE INT
          RA(1,ICLOSE) = RA(1,ICLOSE) + RX*0.25
          RA(2,ICLOSE) = RA(2,ICLOSE) + RY*0.25
          RA(3,ICLOSE) = RA(3,ICLOSE) + RZ*0.25

        ENDDO

      ENDIF
      !!VACANCIES SHOULD BE CREATED AND THE SAME AMOUNT OF SLIGHTLY MORE STABLE
      !!INTS WERE ALSO CREATED. 
      !THIS SUB GUARANTEED TO CREATE ADD_INTS AMOUNT OF FRENKEL PAIRS
      WRITE(*,'(I4,X,A33)')ADD_INTS, "FRENKEL PAIR DEFECTS WERE CREATED"

      END SUBROUTINE

!***********************************************************************!
      SUBROUTINE ADDFRENKEL
      !THIS WAS CREATED TO TEST OPT SUBROUTINE FOR CREATED FRENKEL PAIR

      USE SYSTEM,ONLY: RA,RR,NATOMC,RBOX

      IMPLICIT NONE

      INTEGER CLOSE, I
      REAL*8 RX, RY, RZ, RS1, RS2

      WRITE(*,*) 'STARTING ADD FRENKEL SUBROUTINE'

      RA(1,4500) = RA(1,6500) - 0.2
      RA(2,4500) = RA(2,6500) - 0.2
      RA(3,4500) = RA(3,6500) - 0.2

      RA(1,6500) = RA(1,6500) + 0.2
      RA(2,6500) = RA(2,6500) + 0.2
      RA(3,6500) = RA(3,6500) + 0.2

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
