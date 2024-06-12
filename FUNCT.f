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
!!PERIODIC CONDIITONS ONLY USED FOR FINDING DEFECTS
REAL*8 FUNCTION PBC_DEFECT(X)

  USE SYSTEM, ONLY: RBOX, RBH1, RBH2

  IMPLICIT NONE

  REAL*8 X

  IF ((X .LE. RBH1) .AND. (X .GE. RBH2)) THEN
    PBC_DEFECT = X
  ELSEIF (X .LT. RBH2) THEN
    PBC_DEFECT = X + RBOX
  ELSE
    PBC_DEFECT = X - RBOX
  END IF

  RETURN

END FUNCTION PBC_DEFECT
!***********************************************************************!
SUBROUTINE FIND_DEFECT

  USE SYSTEM
  USE DEFECT

  IMPLICIT NONE

  INTEGER I

  NDINT = 0
  NDVAC = 0
  TI = 0
  TV = 0

  IF (TCRY .EQ. 1) THEN
    CALL SINGLE_CRYSTAL_WS
  ELSE IF (TCRY .EQ. 2) THEN
    CALL POLY_CRYSTAL_WS
  END IF

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
SUBROUTINE SINGLE_CRYSTAL_WS

  USE SYSTEM

  USE DEFECT

  IMPLICIT NONE

  REAL*8 RX, RY, RZ, RS, PBC_DEFECT

  INTEGER I, J, NBOX

  INTEGER IX, IY, IZ, JX, JY, JZ, RX1, RY1, RZ1, RX2, RY2, RZ2

  NBOX = INT(RBOX)*2

      !!SCALE DOWN POSITIONS TO CUBE
  DO I = 1, NATOMC
    RA(1, I) = RA(1, I)/EPSX
    RA(2, I) = RA(2, I)/EPSY
    RA(3, I) = RA(3, I)/EPSZ
  END DO

  DO I = 1, NATOMC

    IF (RA(1, I) .LT. 0.0) RA(1, I) = RA(1, I) + RBOX

    IF (RA(2, I) .LT. 0.0) RA(2, I) = RA(2, I) + RBOX

    IF (RA(3, I) .LT. 0.0) RA(3, I) = RA(3, I) + RBOX

    IF (RA(1, I) .GE. RBOX) RA(1, I) = RA(1, I) - RBOX

    IF (RA(2, I) .GE. RBOX) RA(2, I) = RA(2, I) - RBOX

    IF (RA(3, I) .GE. RBOX) RA(3, I) = RA(3, I) - RBOX

    IX = FLOOR(RA(1, I)/0.5D0)

    IY = FLOOR(RA(2, I)/0.5D0)

    IZ = FLOOR(RA(3, I)/0.5D0)

    IF (IX .GE. NBOX) IX = IX - NBOX

    IF (IY .GE. NBOX) IY = IY - NBOX

    IF (IZ .GE. NBOX) IZ = IZ - NBOX

    JX = IX + 1

    JY = IY + 1

    JZ = IZ + 1

    IF (JX .GE. NBOX) JX = JX - NBOX

    IF (JY .GE. NBOX) JY = JY - NBOX

    IF (JZ .GE. NBOX) JZ = JZ - NBOX

    IF (MOD(IX, 2) .EQ. 1) THEN

      RX1 = IX

      RX2 = JX

    ELSE

      RX1 = JX

      RX2 = IX

    END IF

    IF (MOD(IY, 2) .EQ. 1) THEN

      RY1 = IY

      RY2 = JY

    ELSE

      RY1 = JY

      RY2 = IY

    END IF

    IF (MOD(IZ, 2) .EQ. 1) THEN

      RZ1 = IZ

      RZ2 = JZ

    ELSE

      RZ1 = JZ

      RZ2 = IZ

    END IF

    RX = PBC_DEFECT(DFLOAT(RX1)*0.5D0 - RA(1, I))

    RY = PBC_DEFECT(DFLOAT(RY1)*0.5D0 - RA(2, I))

    RZ = PBC_DEFECT(DFLOAT(RZ1)*0.5D0 - RA(3, I))

    RS = RX*RX + RY*RY + RZ*RZ

    IF (RS .LE. DECUT) THEN

      J = FLOOR(DFLOAT(RX1)*0.5D0)*2*RBOX**2 &
          + FLOOR(DFLOAT(RY1)*0.5D0)*2*RBOX &
          + RZ1 + 1

      TI(I) = TI(I) + 1

      TV(J) = TV(J) + 1

    END IF

    RX = PBC_DEFECT(DFLOAT(RX2)*0.5D0 - RA(1, I))

    RY = PBC_DEFECT(DFLOAT(RY2)*0.5D0 - RA(2, I))

    RZ = PBC_DEFECT(DFLOAT(RZ2)*0.5D0 - RA(3, I))

    RS = RX*RX + RY*RY + RZ*RZ

    IF (RS .LE. DECUT) THEN

      J = FLOOR(DFLOAT(RX2)*0.5D0)*2*RBOX**2 &
          + FLOOR(DFLOAT(RY2)*0.5D0)*2*RBOX &
          + RZ2 + 1

      TI(I) = TI(I) + 1

      TV(J) = TV(J) + 1

    END IF

  END DO

      !!RESCALE ATOM POSITIONS TO ORIGINAL VALUES
  DO I = 1, NATOMC
    RA(1, I) = RA(1, I)*EPSX
    RA(2, I) = RA(2, I)*EPSY
    RA(3, I) = RA(3, I)*EPSZ
  END DO

END SUBROUTINE SINGLE_CRYSTAL_WS
!***********************************************************************!
SUBROUTINE POLY_CRYSTAL_WS

  USE OMP_LIB
  USE SYSTEM
  USE DEFECT
  USE MPI
  USE LMP_DRIVER

  IMPLICIT NONE

  REAL*8 RX, RY, RZ, RS, RMIN, PBC
  INTEGER I, J, MIN, N, ATOM
  INTEGER IDHI, IDLO, RID
  INTEGER OCCBUFF(SIZE(OCC))

  IDLO = NATOMC/NPROCS*(RANK) + 1
  IDHI = NATOMC/NPROCS*(RANK + 1)
  IF (RANK .eq. NPROCS - 1) IDHI = IDHI + MOD(NATOMC, NPROCS)

  CALL MPI_BCAST(RA, SIZE(RA), MPI_REAL8, 0, lmpSystem%getCommunicator(), IERR)
  CALL MPI_BCAST(UPNEIGH, 1, MPI_LOGICAL, 0, lmpSystem%getCommunicator(), IERR)

  IF (MOD(ISTEP - 1, 50) .EQ. 0 .and. UPNEIGH) CALL DEFECT_NEIGH_UPDATE(IDLO, IDHI)

  OC = 0
  OCC = 0
  OCCBUFF = 0

!!DETERMINE THE OCCUPANCY OF EACH SITE

!$OMP PARALLEL PRIVATE(RS,RX,RY,RZ,RMIN,MIN,J,N,ATOM) SHARED(OC,OCC)
!$OMP DO
  DO I = IDLO, IDHI

    RS = 0; RMIN = HUGE(RMIN); MIN = 0
    DO J = 1, NNEIGH
      RID = RNEIGH(J, I)
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

  CALL MPI_REDUCE(TV, TVBUFF, SIZE(TV), MPI_INT, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)
  CALL MPI_REDUCE(TI, TIBUFF, SIZE(TI), MPI_INT, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)

  IF (RANK .eq. 0) THEN
    TV = TVBUFF
    TI = TIBUFF
  END IF

END SUBROUTINE POLY_CRYSTAL_WS
!***********************************************************************!
SUBROUTINE DEFECT_NEIGH_UPDATE(LO, HI)

  USE SYSTEM
  USE DEFECT
  USE BINS, ONLY: maxNeighborSize

  IMPLICIT NONE

  !LO & HI ARE THE IDS FOR WHICH NEIGH LIST WILL BE UPDATED FOR
  INTEGER :: LO, HI
  INTEGER :: I, J, K
  REAL*8 :: RX, RY, RZ, RS, PBC
  REAL*8 START, FINISH
  REAL*8 :: RSLIST(NNEIGH)
  INTEGER :: RFNEIGH(maxNeighborSize)
  if (rank .eq. 0) WRITE (*, '(a)') 'Updating defect neighbor lists'

  UPNEIGH = .FALSE.

  CALL CPU_TIME(START)

  RX = 0; RY = 0; RZ = 0; RS = 0
  RSLIST = 0; RNEIGH = 0

  IF (.NOT. RFGRID%isInitialized) THEN
    CALL RFGRID%initBinGrid(20, 20, 20, RBOX*EPSX, RBOX*EPSY, RBOX*EPSZ)
    CALL RFGRID%assignAtomsToBins(RR, NATOMR)
  END IF

  !FOR EACH ATOM IN THE SYSTEM FIND N NEAREST REFERENCE ATOMS
  DO I = LO, HI
    RSLIST = HUGE(RS)
    RFNEIGH = RFGRID%getNearestNeighbors(RA(1, I), RA(2, I), RA(3, I))
    J = 1
    DO WHILE (RFNEIGH(J) .NE. 0 .AND. J .LE. SIZE(RFNEIGH))
      RX = PBC(RA(1, I) - RR(1, RFNEIGH(J)), 1)
      RY = PBC(RA(2, I) - RR(2, RFNEIGH(J)), 2)
      RZ = PBC(RA(3, I) - RR(3, RFNEIGH(J)), 3)
      RS = RX*RX + RY*RY + RZ*RZ
      IF (RS .lt. RSLIST(NNEIGH)) THEN
        RSLIST(NNEIGH) = RS
        RNEIGH(NNEIGH, I) = RFNEIGH(J)
        CALL quicksort(RSLIST, RNEIGH(:, I), 1, NNEIGH)
      END IF
      J = J + 1
    END DO
  END DO

  CALL CPU_TIME(FINISH)
  IF (RANK .eq. NPROCS - 1) WRITE (*, '(a,f14.5,a)') 'Neighbor list update time= ', FINISH - START, ' s'

END SUBROUTINE DEFECT_NEIGH_UPDATE
!***********************************************************************!
recursive subroutine quicksort(a, b, first, last)
  ! This quicksort routine sorts both a and b arrays based on the values in a
  implicit none
  real*8 a(*), x, t
  integer b(*)
  integer first, last
  integer i, j

  x = a((first + last)/2)
  i = first
  j = last
  do
    do while (a(i) .lt. x)
      i = i + 1
    end do
    do while (x .lt. a(j))
      j = j - 1
    end do
    if (i .ge. j) exit
    t = a(i); a(i) = a(j); a(j) = t
    t = b(i); b(i) = b(j); b(j) = t
    i = i + 1
    j = j - 1
  end do
  if (first .lt. i - 1) call quicksort(a, b, first, i - 1)
  if (j + 1 .lt. last) call quicksort(a, b, j + 1, last)
end subroutine quicksort
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

  IMPLICIT NONE

  INTEGER I, J

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
  WRITE (30, '(A)') 'ITEM: ATOMS x y z'

  DO I = 1, NATOMC
    WRITE (30, '(3F14.5)') (RA(J, I), J=1, 3)
  END DO

  CLOSE (30)

END SUBROUTINE OVITO
!***********************************************************************!
SUBROUTINE BASIC_DAMAGE(NFRKL)

  USE SYSTEM, ONLY: RA, NATOMC, RANK, DECUT, TDAM

  IMPLICIT NONE

  INTEGER :: I, J, ID, NFRKL, N
  INTEGER :: MOVEID(NFRKL)
  REAL*8 :: DLIM
  REAL*8 :: RND, NEWPOS(3)
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
      DO J = 1, N - 1
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

  IMPLICIT NONE

  INTEGER :: I, NFRKL
  INTEGER :: IDS(NFRKL)
  CHARACTER(LEN=10) :: IDSTR

  CALL MPI_BCAST(IDS, NFRKL, MPI_INT, 0, lmpSystem%getCommunicator(), IERR)

  !Update lammps on what atoms were deleted
  DO I = 1, NFRKL
    WRITE (IDSTR, '(I0)') IDS(I)
    CALL lammps_command(lmpSystem%getLammpsPtr(), 'group deleteList id '//TRIM(IDSTR))
  END DO
  CALL lammps_command(lmpSystem%getLammpsPtr(), 'delete_atoms group deleteList')

  !Upadate all other processes on changed parameters
  NATOMC = int(lammps_get_natoms(lmpSystem%getLammpsPtr()))
  call lmpSystem%getPositions

END SUBROUTINE CREATE_VACANCIES
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
  ELSE IF (TDAM .EQ. 2) THEN
    DTIME = (FLOAT(NFREKL)/FLOAT(NATOMC))/DRATE
    DPA = DPA + (NFREKL/NATOMC)
  ELSE IF (TDAM .EQ. 3) THEN
    DTIME = HUGE(DTIME)
    DPA = 0
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
