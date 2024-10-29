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
REAL*8 FUNCTION PBC(X, DIM)

  USE SYSTEM, ONLY: RBOX, RBH1, RBH2, EPSX, EPSY, EPSZ

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

  IF ((X .LE. RBH1*EPS) .AND. (X .GE. RBH2*EPS)) THEN
    PBC = X
  ELSEIF (X .LT. RBH2*EPS) THEN
    PBC = X + RBOX*EPS
  ELSE
    PBC = X - RBOX*EPS
  END IF

  RETURN

END FUNCTION PBC
!***********************************************************************!
SUBROUTINE FIND_DEFECT

  USE SYSTEM
  USE DEFECT
  use MicroStrucFeatures, only: features
  use LMP_DRIVER, only: lmpSystem

  IMPLICIT NONE

  INTEGER I

  NDINT = 0
  NDVAC = 0
  TI = 0
  TV = 0

  CALL WignerSeitzAnalysis

  !!TALLY NUMBER OF VACANCIES AND INTERSTITIALS

  DO I = 1, NATOMC
    IF (TI(I) .NE. 1) NDINT = NDINT + 1
  END DO

  DO I = 1, NATOMR
    IF (TV(I) .NE. 1) NDVAC = NDVAC + 1
  END DO

  RETURN
END SUBROUTINE FIND_DEFECT
!***********************************************************************!
SUBROUTINE WignerSeitzAnalysis

  USE OMP_LIB
  USE SYSTEM
  USE DEFECT
  USE MPI
  USE LMP_DRIVER
  use MicroStrucFeatures

  IMPLICIT NONE

  REAL*8 RX, RY, RZ, RS, RMIN, FACTOR, PBC
  INTEGER I, J, MIN, N, ATOM
  INTEGER IDHI, IDLO, RID
  INTEGER OCCBUFF(SIZE(OCC))

  IDLO = NATOMC/NPROCS*(RANK) + 1
  IDHI = NATOMC/NPROCS*(RANK + 1)
  IF (RANK .eq. NPROCS - 1) IDHI = IDHI + MOD(NATOMC, NPROCS)

  CALL MPI_BCAST(RA, SIZE(RA), MPI_REAL8, 0, lmpSystem%getCommunicator(), IERR)

  CALL RFGRID%initBinGrid(30, 30, 30, RBOX*EPSX, RBOX*EPSY, RBOX*EPSZ)
  CALL RFGRID%assignAtomsToBins(RR, NATOMR)

  OC = 0
  OCC = 0
  OCCBUFF = 0
  TVBUFF = 0
  TIBUFF = 0

  FACTOR = maxval((/epsx,epsy,epsz/),dim=1)
  
!!DETERMINE THE OCCUPANCY OF EACH SITE

!$OMP PARALLEL PRIVATE(RS,RX,RY,RZ,RMIN,RID,MIN,J,N,ATOM) SHARED(OC,OCC)
!$OMP DO
  DO I = IDLO, IDHI
    if (features%grains%grainOccupy(I) == 0) cycle
    RS = 0; RMIN = HUGE(RMIN); MIN = 0
    call RFGRID%getNearestNeighbors(RA(1, I), RA(2, I), RA(3, I), 2.0D0*FACTOR)
    if (RFGRID%singleListCount == 0) print *, RA(1, I), RA(2, I), RA(3, I)
    DO J = 1, RFGRID%singleListCount
      RID = RFGRID%singleNeighList(j)
      RX = PBC(RA(1, I) - RR(1, RID), 1)
      RY = PBC(RA(2, I) - RR(2, RID), 2)
      RZ = PBC(RA(3, I) - RR(3, RID), 3)
      RS = RX*RX + RY*RY + RZ*RZ
      IF (RS .LT. RMIN) THEN
        MIN = RID
        RMIN = RS
      END IF
    END DO

!$OMP CRITICAL
    OCC(MIN) = OCC(MIN) + 1
    OC(OCC(MIN), MIN) = I
!$OMP END CRITICAL

  END DO
!$OMP END DO

!$OMP SINGLE
  CALL MPI_ALLREDUCE(OCC, OCCBUFF, SIZE(OCC), MPI_INT, MPI_SUM, lmpSystem%getCommunicator(), IERR)
!$OMP END SINGLE

!!USING OCCUPANCY TO DETERMINE DEFECTS
!$OMP DO
  DO I = 1, NATOMR

    N = OCC(I)
    IF (N .GT. 0) THEN
      TV(I) = N
      DO J = 1, N
        ATOM = OC(J, I)
        TI(ATOM) = OCCBUFF(I)
      END DO
    END IF

  END DO
!$OMP END DO
!$OMP END PARALLEL

  ! Remove the grain boundary ints from defect list
  do i = IDLO, IDHI
    if (features%grains%grainOccupy(i) == 0)then
      TI(i) = 1
    end if
  end do

  CALL MPI_REDUCE(TV, TVBUFF, SIZE(TV), MPI_INT, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)
  CALL MPI_REDUCE(TI, TIBUFF, SIZE(TI), MPI_INT, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)

  TV = TVBUFF
  TI = TIBUFF

  CALL MPI_BCAST(TV, SIZE(TV), MPI_INT, 0, lmpSystem%getCommunicator(), ierr)
  CALL MPI_BCAST(TI, SIZE(TI), MPI_INT, 0, lmpSystem%getCommunicator(), ierr)

END SUBROUTINE WignerSeitzAnalysis
!***********************************************************************!
recursive subroutine quicksort_int(array1, array2, lo, hi)
    implicit none
    integer :: array1(*), array2(*)
    integer :: pivot, right, left, temp, n, lo, hi

    n = size(array1(lo:hi))

    if (n > 1) then
      pivot = array1(lo)
      left = lo
      right = hi
      do while (left .lt. right)
        do while (array1(left) < pivot)
          left = left + 1
        end do
        do while (array1(right) > pivot)
          right = right - 1
        end do
        if (left <= right) then
          temp = array1(left)
          array1(left) = array1(right)
          array1(right) = temp
          temp = array2(left)
          array2(left) = array2(right)
          array2(right) = temp
          left = left + 1
          right = right - 1
        end if
      end do

      call quicksort_int(array1, array2, lo, right)
      call quicksort_int(array1, array2, left, hi)

    end if

end subroutine quicksort_int
!***********************************************************************!
!!OUTPUT DEFECTS FOR OVITO

SUBROUTINE OVITO

  USE SYSTEM

  USE DEFECT

  USE MicroStrucFeatures, only: features

  IMPLICIT NONE

  INTEGER I

  CHARACTER(len=200) FILENAME

!!CREATING THE CURRENT DEFECT POSITION FILE

  WRITE (FILENAME, '(A8, I0, A4)') './DUMP/D', ISTEP, '.OUT'

  OPEN (UNIT=29, FILE=TRIM(FILENAME))

  WRITE (29, '(A)') 'ITEM: TIMESTEP'
  WRITE (29, *) ISTEP
  WRITE (29, '(A)') 'ITEM: NUMBER OF ATOMS'
  WRITE (29, *) NDVAC + NDINT
  WRITE (29, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
  WRITE (29, '(2F10.4,I5)') 0.00, RBOX*EPSX, 3
  WRITE (29, '(2F10.4,I5)') 0.00, RBOX*EPSY, 3
  WRITE (29, '(2F10.4,I5)') 0.00, RBOX*EPSZ, 3
  WRITE (29, '(A)') 'ITEM: ATOMS x y z type'

  DO I = 1, NATOMC

    IF (TI(I) .NE. 1) THEN

      WRITE (29, '(3F12.5,I5)') RA(:, I), 2

    END IF

  END DO

  DO I = 1, NATOMR

    IF (TV(I) .NE. 1) THEN

      WRITE (29, '(3F12.5,I5)') RR(:, I), 1

    END IF

  END DO

  CLOSE (29)

!!CREATING DISLOCATION ANALYSIS FILE

  WRITE (FILENAME, '(A12, I0, A4)') './DUMP/DXA/A', ISTEP, '.OUT'

  OPEN (UNIT=30, FILE=TRIM(FILENAME))

  WRITE (30, '(A)') 'ITEM: TIMESTEP'
  WRITE (30, *) ISTEP
  WRITE (30, '(A)') 'ITEM: NUMBER OF ATOMS'
  WRITE (30, *) NATOMC
  WRITE (30, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
  WRITE (30, '(2F10.4,I5)') 0.00, RBOX*EPSX, 3
  WRITE (30, '(2F10.4,I5)') 0.00, RBOX*EPSY, 3
  WRITE (30, '(2F10.4,I5)') 0.00, RBOX*EPSZ, 3
  WRITE (30, '(A)') 'ITEM: ATOMS x y z grain'

  DO I = 1, NATOMC
    WRITE (30, '(3F14.5,1I14)') RA(1, I),RA(2, I),RA(3, I),features%grains%grainOccupy(i)
  END DO

  CLOSE (30)

END SUBROUTINE OVITO
!***********************************************************************!
SUBROUTINE BASIC_DAMAGE(NFRKL)

  USE SYSTEM, ONLY: RA, NATOMC, RANK, DECUT, TDAM
  USE MicroStrucFeatures

  IMPLICIT NONE

  INTEGER :: I, J, ID, NFRKL, N
  INTEGER :: MOVEID(NFRKL)
  REAL*8 :: DLIM
  REAL*8 :: RND, NEWPOS(3)
  LOGICAL :: INSERT

  DLIM = DECUT

  N = 0

  MOVEID = 0

  IF (NFRKL .LE. 0) RETURN

  CALL RANDOM_SEED()

      !!PICK NFRKL ATOMS TO MOVE (SHOULD CREATE NFKRL # OF FRENKEL PAIRS)

  IF (RANK .EQ. 0) THEN

    N = 1

    DO WHILE (N .LE. NFRKL)

      CALL RANDOM_NUMBER(RND)
      ID = CEILING(RND*NATOMC)
      IF (ID .EQ. 0) ID = 1

      ! Do not insert chosen atom if it already is in the list
      ! or if it grain boundary atom
      INSERT = .TRUE.
      DO J = 1, N - 1
        IF (ID .EQ. MOVEID(J) .or. features%grains%grainOccupy(ID) == 0) INSERT = .FALSE.
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

      WRITE (*, '(A,I0,A)') 'Generating ', NFRKL, ' Frenkel pairs'

      DO I = 1, NFRKL

        CALL FRENKEL_POS_UPDATE(NEWPOS, DLIM)
        RA(:, MOVEID(I)) = NEWPOS

      END DO

    END IF

  ELSE IF (TDAM .EQ. 3) THEN

    IF (RANK .EQ. 0) WRITE (*, '(A,I0,A)') 'Generating ', NFRKL, ' single vacancies'

    CALL CREATE_VACANCIES(NFRKL, MOVEID)

  END IF

END SUBROUTINE BASIC_DAMAGE
!***********************************************************************!
SUBROUTINE CREATE_VACANCIES(NFRKL, IDS)

  USE SYSTEM
  USE LMP_DRIVER
  USE MPI
  use MicroStrucFeatures

  IMPLICIT NONE

  INTEGER :: I, NFRKL
  INTEGER :: IDS(NFRKL)
  CHARACTER(LEN=10) :: IDSTR

  INTEGER :: ONATOMC
  REAL*8, DIMENSION(3,NATOMC) :: ORA

  CALL MPI_BCAST(IDS, NFRKL, MPI_INT, 0, lmpSystem%getCommunicator(), IERR)

  !Update lammps on what atoms were deleted
  DO I = 1, NFRKL
    WRITE (IDSTR, '(I0)') IDS(I)
    CALL lammps_command(lmpSystem%getLammpsPtr(), 'group deleteList id '//TRIM(IDSTR))
  END DO
  CALL lammps_command(lmpSystem%getLammpsPtr(), 'delete_atoms group deleteList')

  ONATOMC = NATOMC
  ORA = RA

  !Upadate all other processes on changed parameters
  NATOMC = int(lammps_get_natoms(lmpSystem%getLammpsPtr()))
  call lmpSystem%getPositions

  !Fix grain assignment since ordering got messed up
  call grainAssignmentFix(ORA,ONATOMC)

END SUBROUTINE CREATE_VACANCIES
!***********************************************************************!
SUBROUTINE grainAssignmentFix(OLDPOS, OLDC)

  USE SYSTEM
  USE MicroStrucFeatures
  USE BINS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: OLDC
  REAL*8, DIMENSION(3,OLDC), INTENT(IN) :: OLDPOS

  INTEGER :: I,J,K
  INTEGER, DIMENSION(OLDC) :: OGA ! OLD GRAIN ASSIGNMENTS
  TYPE(BinGrid) :: BG

  OGA = features%grains%grainOccupy

  CALL BG%initBinGrid(30, 30, 30, RBOX*EPSX, RBOX*EPSY, RBOX*EPSZ)
  CALL BG%assignAtomsToBins(OLDPOS, OLDC)

  DO I=1,NATOMC

    CALL BG%getNearestNeighbors(RA(1,I),RA(2,I),RA(3,I),1.0D0)
    DO J = 1, BG%singleListCount
      K = BG%singleNeighList(J)
      IF (norm2(RA(:,I)-OLDPOS(:,K)) < 0.01) THEN
        features%grains%grainOccupy(i) = OGA(k)
        if (OGA(k) > 3 .and. rank == 0) print *, OGA(k), k
        exit
      END IF
    END DO

  END DO

END SUBROUTINE grainAssignmentFix
!***********************************************************************!
SUBROUTINE FRENKEL_POS_UPDATE(POS, THRESH)

  USE SYSTEM, ONLY: RA, NATOMC, RBOX, EPSX, EPSY, EPSZ

  IMPLICIT NONE

  INTEGER :: I, ID_NEAR
  REAL*8, INTENT(INOUT) :: POS(3)
  REAL*8, INTENT(IN) :: THRESH
  REAL*8 :: RX, RY, RZ, RS, RS_NEAR
  REAL*8 :: UNIT(3)
  LOGICAL :: FLAG

  ID_NEAR = 0; RS_NEAR = HUGE(RS_NEAR); FLAG = .TRUE.

  CALL RANDOM_NUMBER(POS)

  POS(1) = POS(1)*RBOX*EPSX
  POS(2) = POS(2)*RBOX*EPSY
  POS(3) = POS(3)*RBOX*EPSZ

      !!FIND CLOSEST ATOM TO (X,Y,Z)

  DO I = 1, NATOMC

    RX = RA(1, I) - POS(1)

    RY = RA(2, I) - POS(2)

    RZ = RA(3, I) - POS(3)

    RS = RX*RX + RY*RY + RZ*RZ

    IF (RS .LT. RS_NEAR) THEN

      RS_NEAR = RS

      ID_NEAR = I

    END IF

  END DO

      !!! UPDATE THE POSITION UNTIL IT IS SUITABLE FOR RELAXATION

  DO WHILE (FLAG)

    FLAG = .FALSE.

    UNIT(1) = POS(1) - RA(1, ID_NEAR)

    UNIT(2) = POS(2) - RA(2, ID_NEAR)

    UNIT(3) = POS(3) - RA(3, ID_NEAR)

    UNIT = UNIT/NORM2(UNIT)

    POS = POS + SQRT(THRESH)/3*UNIT

    RX = POS(1) - RA(1, ID_NEAR)

    RY = POS(2) - RA(2, ID_NEAR)

    RZ = POS(3) - RA(3, ID_NEAR)

    RS = RX*RX + RY*RY + RZ*RZ

    IF (RS .LT. THRESH) FLAG = .TRUE.

  END DO

END SUBROUTINE FRENKEL_POS_UPDATE
!***********************************************************************!
SUBROUTINE MD_DPA

  USE SYSTEM, ONLY: NATOMC, DPA, AVDPA, DPAC, DTIME, DRATE, TDAM, CALC, NFREKL
  USE DEFECT, ONLY: ONDSTOT, NDSTOT

  IMPLICIT NONE

  REAL*8 :: NEWDSP

  !UPDATE AVERAGE DPA VALUE

  NEWDSP = NDSTOT - ONDSTOT

  AVDPA = (AVDPA*DPAC + NEWDSP/NATOMC)/(DPAC + 1)

  DPAC = DPAC + 1

  !UPDATE TIME BETWEEN CASCADES

  IF (TDAM .EQ. 1) THEN
    DTIME = AVDPA/DRATE
    DPA = DPA + NEWDSP/NATOMC
  ELSE IF (TDAM .EQ. 2 .OR. TDAM .EQ. 3) THEN
    DTIME = (FLOAT(NFREKL)/FLOAT(NATOMC))/DRATE
    DPA = DPA + (FLOAT(NFREKL)/FLOAT(NATOMC))
  END IF

  !ADD AVERAGE DPA PER CASCADE TO CURRENT DPA VALUE

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

  CHARACTER(len=200) FILENAME

  IF (RANK .EQ. 0) THEN

    IF (MOD(ISTEP, 1000) .EQ. 0) THEN

      WRITE (FILENAME, '(A17, I0, A4)') './RESTART/RESTART', &
        ISTEP, '.OUT'

      OPEN (UNIT=27, FILE=TRIM(FILENAME))

      WRITE (27, *) NATOMC

      DO I = 1, NATOMC

        WRITE (27, '(3F14.5)') (RA(J, I), J=1, 3)

      END DO

      !SIMULATION MEASUREMENTS AND PARAMETERS

      WRITE (27, *)

      WRITE (27, *) &
        'ISTEP, TIME, TIME_LAST, DTIME, DPA, DPAC, AVDPA, ONDSTOT, EPSX, EPSY, EPSZ'

      WRITE (27, *) ISTEP
      WRITE (27, *) TIME
      WRITE (27, *) TIME_LAST
      WRITE (27, *) DTIME
      WRITE (27, *) DPA
      WRITE (27, *) DPAC
      WRITE (27, *) AVDPA
      WRITE (27, *) ONDSTOT
      WRITE (27, *) EPSX
      WRITE (27, *) EPSY
      WRITE (27, *) EPSZ

      !AVERAGE DEFECT MIGRATION ENERGIES

      WRITE (27, *)

      WRITE (27, *) 'MIGRATION ENERGIES (VACANCIES THEN INTERSTITIALS)'

      DO I = 1, 50
        WRITE (27, *) EMIG(1, I), EMIC(1, I)
      END DO

      DO I = 1, 50
        WRITE (27, *) EMIG(2, I), EMIC(2, I)
      END DO

    END IF

    CLOSE (27)

  END IF

END SUBROUTINE
!***********************************************************************!
SUBROUTINE DEFECT_PRINT(DEFECTLIST, NDEFECTS)

  USE SYSTEM, ONLY: NDEMAX

  IMPLICIT NONE

  REAL*8 :: DEFECTLIST(7, NDEMAX)
  INTEGER :: NDEFECTS

  INTEGER :: I, J
  INTEGER :: INTSIZES(MAXVAL(INT(DEFECTLIST(6, :)))), VACSIZES(ABS(MINVAL(INT(DEFECTLIST(6, :)))))
  INTEGER :: INTCLUSTERCOUNTS(MAXVAL(INT(DEFECTLIST(6, :)))), VACCLUSTERCOUNTS(ABS(MINVAL(INT(DEFECTLIST(6, :)))))
  INTEGER :: INTCOUNT, VACCOUNT
  CHARACTER(LEN=20) :: PRINTFORMAT

  INTCOUNT = 0; VACCOUNT = 0
  INTCLUSTERCOUNTS = 0; VACCLUSTERCOUNTS = 0

  !DETERMINE WHAT CLUSTER SIZES TO PRINT IN THE TABLE
  DO I = 1, NDEFECTS
    IF (.NOT. ANY(INTSIZES == INT(DEFECTLIST(6, I))) .AND. INT(DEFECTLIST(6, I)) > 0) THEN
      INTSIZES(INTCOUNT + 1) = INT(DEFECTLIST(6, I))
      INTCOUNT = INTCOUNT + 1
    ELSE IF (.NOT. ANY(VACSIZES == INT(DEFECTLIST(6, I))) .AND. INT(DEFECTLIST(6, I)) < 0) THEN
      VACSIZES(VACCOUNT + 1) = INT(DEFECTLIST(6, I))
      VACCOUNT = VACCOUNT + 1
    END IF
  END DO

  ! DETERMINE THE COUNT FOR EACH DEFECT SIZE
  DO I = 1, NDEFECTS
    IF (DEFECTLIST(6, I) > 0) THEN
      DO J = 1, INTCOUNT
        IF (INT(DEFECTLIST(6, I)) == INTSIZES(J)) THEN
          INTCLUSTERCOUNTS(J) = INTCLUSTERCOUNTS(J) + 1
        END IF
      END DO
    ELSE
      DO J = 1, VACCOUNT
        IF (INT(DEFECTLIST(6, I)) == VACSIZES(J)) THEN
          VACCLUSTERCOUNTS(J) = VACCLUSTERCOUNTS(J) + 1
        END IF
      END DO
    END IF
  END DO

  ! SORT THE DEFECT LISTS
  CALL QUICKSORT_INT(INTSIZES, INTCLUSTERCOUNTS, 1, INTCOUNT)
  CALL QUICKSORT_INT(VACSIZES, VACCLUSTERCOUNTS, 1, VACCOUNT)

  PRINTFORMAT = '(X,I4,X)'

  !WRITE THE INT TABLE
  IF (INTCOUNT > 0) THEN

    WRITE (*, '(A9)', ADVANCE="NO") 'Int Size|'
    DO I = 1, INTCOUNT
      WRITE (*, PRINTFORMAT, ADVANCE="NO") INTSIZES(I)
    END DO
    WRITE (*, *) !NEWLINE

    WRITE (*, '(A9)', ADVANCE="NO") 'Count|'
    DO I = 1, INTCOUNT
      WRITE (*, PRINTFORMAT, ADVANCE="NO") INTCLUSTERCOUNTS(I)
    END DO
    WRITE (*, *) !NEWLINE
    WRITE (*, *) !NEWLINE

  END IF

  !WRITE THE VAC TABLE
  IF (VACCOUNT > 0) THEN

    WRITE (*, '(A9)', ADVANCE="NO") 'Vac Size|'
    DO I = VACCOUNT, 1, -1
      WRITE (*, PRINTFORMAT, ADVANCE="NO") VACSIZES(I)
    END DO
    WRITE (*, *) !NEWLINE

    WRITE (*, '(A9)', ADVANCE="NO") 'Count|'
    DO I = VACCOUNT, 1, -1
      WRITE (*, PRINTFORMAT, ADVANCE="NO") VACCLUSTERCOUNTS(I)
    END DO
    WRITE (*, *) !NEWLINE
    WRITE (*, *) !NEWLINE

  END IF
END SUBROUTINE DEFECT_PRINT