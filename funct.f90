!!APPLY THE PERIODIC BOUNDARY CONDITION TO THE SYSTEM
REAL*8 FUNCTION PBC(X, DIM)

  use system, only: BOXX, BOXY, BOXZ

  implicit none

  INTEGER DIM
  REAL*8 X

  real*8 half

  if (DIM .EQ. 1) then
    half = BOXX/2
  elseif (DIM .EQ. 2) then
    half = BOXY/2
  elseif (DIM .EQ. 3) then
    half = BOXZ/2
  end if

  if ((X .LE. half) .AND. (X .GE. -half)) then
    PBC = X
  elseif (X .LT. -half) then
    PBC = X + half*2
  else
    PBC = X - half*2
  end if

  RETURN

end FUNCTION PBC
!***********************************************************************!
subroutine FIND_DEFECT

  use system
  use defect

  implicit none

  INTEGER I

  NDINT = 0
  NDVAC = 0
  TI = 0
  TV = 0

  CALL WignerSeitzAnalysis

  !!TALLY NUMBER OF VACANCIES AND INTERSTITIALS

  do I = 1, NATOMC
    if (TI(I) .NE. 1) NDINT = NDINT + 1
  end do

  do I = 1, NATOMR
    if (TV(I) .NE. 1) NDVAC = NDVAC + 1
  end do

  RETURN
end subroutine FIND_DEFECT
!***********************************************************************!
subroutine WignerSeitzAnalysis

  USE OMP_LIB
  USE SYSTEM
  USE DEFECT
  USE MPI
  USE LMP_DRIVER
  use Clusters

  implicit none

  REAL*8 RX, RY, RZ, RS, RMIN, PBC
  INTEGER I, J, MIN, N, ATOM
  INTEGER IDHI, IDLO, RID
  INTEGER OCCBUFF(SIZE(OCC))

  call defects%init(NDEMAX)

  IDLO = NATOMC/NPROCS*(RANK) + 1
  IDHI = NATOMC/NPROCS*(RANK + 1)
  if (RANK .eq. NPROCS - 1) IDHI = IDHI + MOD(NATOMC, NPROCS)

  CALL MPI_BCAST(RA, SIZE(RA), MPI_REAL8, 0, lmpSystem%getCommunicator(), IERR)

  CALL RFGRID%initBinGrid(30, 30, 30, BOXX, BOXY, BOXZ)
  CALL RFGRID%assignAtomsToBins(RR, NATOMR)

  OC = 0
  OCC = 0
  OCCBUFF = 0
  TVBUFF = 0
  TIBUFF = 0
  
!!DETERMINE THE OCCUPANCY OF EACH SITE

!$OMP PARALLEL PRIVATE(RS,RX,RY,RZ,RMIN,RID,MIN,J,N,ATOM) SHARED(OC,OCC)
!$OMP do
  do I = IDLO, IDHI
    RS = 0; RMIN = HUGE(RMIN); MIN = 0
    call RFGRID%getNearestNeighbors(RA(1, I), RA(2, I), RA(3, I), 2.0D0)
    if (RFGRID%singleListCount == 0) print *, RA(1, I), RA(2, I), RA(3, I)
    do J = 1, RFGRID%singleListCount
      RID = RFGRID%singleNeighList(j)
      RX = PBC(RA(1, I) - RR(1, RID), 1)
      RY = PBC(RA(2, I) - RR(2, RID), 2)
      RZ = PBC(RA(3, I) - RR(3, RID), 3)
      RS = RX*RX + RY*RY + RZ*RZ
      if (RS .LT. RMIN) then
        MIN = RID
        RMIN = RS
      end if
    end do

!$OMP CRITICAL
    OCC(MIN) = OCC(MIN) + 1
    OC(OCC(MIN), MIN) = I
!$OMP end CRITICAL

  end do
!$OMP end do

!$OMP SINGLE
  CALL MPI_ALLREDUCE(OCC, OCCBUFF, SIZE(OCC), MPI_INT, MPI_SUM, lmpSystem%getCommunicator(), IERR)
!$OMP end SINGLE

!!USING OCCUPANCY TO DETERMINE DEFECTS
!$OMP do
  do I = 1, NATOMR

    N = OCC(I)
    if (N .GT. 0) then
      TV(I) = N
      do J = 1, N
        ATOM = OC(J, I)
        TI(ATOM) = OCCBUFF(I)
      end do
    end if

  end do
!$OMP end do
!$OMP end PARALLEL

  CALL MPI_REDUCE(TV, TVBUFF, SIZE(TV), MPI_INT, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)
  CALL MPI_REDUCE(TI, TIBUFF, SIZE(TI), MPI_INT, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)

  TV = TVBUFF
  TI = TIBUFF

  CALL MPI_BCAST(TV, SIZE(TV), MPI_INT, 0, lmpSystem%getCommunicator(), ierr)
  CALL MPI_BCAST(TI, SIZE(TI), MPI_INT, 0, lmpSystem%getCommunicator(), ierr)

end subroutine WignerSeitzAnalysis
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

subroutine OVITO

  USE SYSTEM

  USE DEFECT

  implicit none

  INTEGER I

  CHARACTER(len=200) FILENAME

!!CREATING THE CURRENT DEFECT POSITION FILE

  write (FILENAME, '(A8, I0, A4)') './DUMP/D', ISTEP, '.OUT'

  open (UNIT=29, FILE=TRIM(FILENAME))

  write (29, '(A)') 'ITEM: TIMESTEP'
  write (29, *) ISTEP
  write (29, '(A)') 'ITEM: NUMBER OF ATOMS'
  write (29, *) NDVAC + NDINT
  write (29, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
  write (29, '(2F10.4,I5)') 0.00, BOXX, 3
  write (29, '(2F10.4,I5)') 0.00, BOXY, 3
  write (29, '(2F10.4,I5)') 0.00, BOXZ, 3
  write (29, '(A)') 'ITEM: ATOMS x y z type'

  do I = 1, NATOMC

    if (TI(I) .NE. 1) then

      write (29, '(3F12.5,I5)') RA(:, I), 2

    end if

  end do

  do I = 1, NATOMR

    if (TV(I) .NE. 1) then

      write (29, '(3F12.5,I5)') RR(:, I), 1

    end if

  end do

  close (29)

!!CREATING DISLOCATION ANALYSIS FILE

  write (FILENAME, '(A12, I0, A4)') './DUMP/DXA/A', ISTEP, '.OUT'

  open (UNIT=30, FILE=TRIM(FILENAME))

  write (30, '(A)') 'ITEM: TIMESTEP'
  write (30, *) ISTEP
  write (30, '(A)') 'ITEM: NUMBER OF ATOMS'
  write (30, *) NATOMC
  write (30, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
  write (30, '(2F15.8,I1)') 0.00, BOXX, 3
  write (30, '(2F15.8,I1)') 0.00, BOXY, 3
  write (30, '(2F15.8,I1)') 0.00, BOXZ, 3
  write (30, '(A)') 'ITEM: ATOMS x y z'

  do I = 1, NATOMC
    write (30, '(3F14.5)') RA(1, I),RA(2, I),RA(3, I)
  end do

  close (30)

end subroutine OVITO
!***********************************************************************!
subroutine BASIC_DAMAGE(NFRKL)

  USE SYSTEM, ONLY: RA, NATOMC, RANK, DECUT, TDAM

  implicit none

  INTEGER :: I, J, ID, NFRKL, N
  INTEGER :: MOVEID(NFRKL)
  REAL*8 :: DLIM
  REAL*8 :: RND, NEWPOS(3)
  LOGICAL :: INSERT

  DLIM = DECUT

  N = 0

  MOVEID = 0

  if (NFRKL .LE. 0) RETURN

  CALL RANdoM_SEED()

      !!PICK NFRKL ATOMS TO MOVE (SHOULD CREATE NFKRL # OF FRENKEL PAIRS)

  if (RANK .EQ. 0) then

    N = 1

    do WHILE (N .LE. NFRKL)

      CALL RANdoM_NUMBER(RND)
      ID = CEILING(RND*NATOMC)
      if (ID .EQ. 0) ID = 1

      ! Do not insert chosen atom if it already is in the list
      INSERT = .TRUE.
      do J = 1, N - 1
        if (ID .EQ. MOVEID(J)) INSERT = .FALSE.
      end do
      if (INSERT) then
        MOVEID(N) = ID
        N = N + 1
      end if

    end do

  end if

      !!MOVE CHOSEN ATOMS TO RANdoM LOCATIONS

  if (TDAM .EQ. 2) then

    if (RANK .EQ. 0) then

      write (*, '(A,I0,A)') 'Generating ', NFRKL, ' Frenkel pairs'

      do I = 1, NFRKL

        CALL FRENKEL_POS_UPDATE(NEWPOS, DLIM)
        RA(:, MOVEID(I)) = NEWPOS

      end do

    end if

  ELSE if (TDAM .EQ. 3) then

    if (RANK .EQ. 0) write (*, '(A,I0,A)') 'Generating ', NFRKL, ' single vacancies'

    CALL CREATE_VACANCIES(NFRKL, MOVEID)

  end if

end subroutine BASIC_DAMAGE
!***********************************************************************!
subroutine CREATE_VACANCIES(NFRKL, IDS)

  USE SYSTEM
  USE LMP_DRIVER
  USE MPI

  implicit none

  INTEGER :: I, NFRKL
  INTEGER :: IDS(NFRKL)
  CHARACTER(LEN=10) :: IDSTR

  INTEGER :: ONATOMC
  REAL*8, DIMENSION(3,NATOMC) :: ORA

  CALL MPI_BCAST(IDS, NFRKL, MPI_INT, 0, lmpSystem%getCommunicator(), IERR)

  !Update lammps on what atoms were deleted
  do I = 1, NFRKL
    write (IDSTR, '(I0)') IDS(I)
    CALL lammps_command(lmpSystem%getLammpsPtr(), 'group deleteList id '//TRIM(IDSTR))
  end do
  CALL lammps_command(lmpSystem%getLammpsPtr(), 'delete_atoms group deleteList')

  ONATOMC = NATOMC
  ORA = RA

  !Upadate all other processes on changed parameters
  NATOMC = int(lammps_get_natoms(lmpSystem%getLammpsPtr()))
  call lmpSystem%getPositions

end subroutine CREATE_VACANCIES
!***********************************************************************!
subroutine FRENKEL_POS_UPDATE(POS, THRESH)

  USE SYSTEM, ONLY: RA, NATOMC, BOXX, BOXY, BOXZ

  implicit none

  INTEGER :: I, ID_NEAR
  REAL*8, INTENT(INOUT) :: POS(3)
  REAL*8, INTENT(IN) :: THRESH
  REAL*8 :: RX, RY, RZ, RS, RS_NEAR
  REAL*8 :: UNIT(3)
  LOGICAL :: FLAG

  ID_NEAR = 0; RS_NEAR = HUGE(RS_NEAR); FLAG = .TRUE.

  call random_number(POS)

  POS(1) = POS(1)*BOXX
  POS(2) = POS(2)*BOXY
  POS(3) = POS(3)*BOXZ

      !!FIND closeST ATOM TO (X,Y,Z)

  do I = 1, NATOMC

    RX = RA(1, I) - POS(1)

    RY = RA(2, I) - POS(2)

    RZ = RA(3, I) - POS(3)

    RS = RX*RX + RY*RY + RZ*RZ

    if (RS .LT. RS_NEAR) then

      RS_NEAR = RS

      ID_NEAR = I

    end if

  end do

      !!! UPDATE THE POSITION UNTIL IT IS SUITABLE FOR RELAXATION

  do while (FLAG)

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

    if (RS .LT. THRESH) FLAG = .TRUE.

  end do

end subroutine FRENKEL_POS_UPDATE
!***********************************************************************!
subroutine MD_DPA

  USE SYSTEM, ONLY: NATOMC, DPA, AVDPA, DPAC, DTIME, DRATE, TDAM, CALC, NFREKL
  USE DEFECT, ONLY: ONDSTOT, NDSTOT

  implicit none

  REAL*8 :: NEWDSP

  !UPDATE AVERAGE DPA VALUE

  NEWDSP = NDSTOT - ONDSTOT

  AVDPA = (AVDPA*DPAC + NEWDSP/NATOMC)/(DPAC + 1)

  DPAC = DPAC + 1

  !UPDATE TIME BETWEEN CASCADES

  if (TDAM .EQ. 1) then
    DTIME = AVDPA/DRATE
    DPA = DPA + NEWDSP/NATOMC
  ELSE if (TDAM .EQ. 2 .OR. TDAM .EQ. 3) then
    DTIME = (FLOAT(NFREKL)/FLOAT(NATOMC))/DRATE
    DPA = DPA + (FLOAT(NFREKL)/FLOAT(NATOMC))
  end if

  !ADD AVERAGE DPA PER CASCADE TO CURRENT DPA VALUE

  CALC = .FALSE.

end subroutine MD_DPA
!***********************************************************************!
subroutine RESTART
  !Create all the necessary restart files at a frequency matching
  !'writeFreq'. A LAMMPS data restart file is created first. Then
  !internal DRLD paramters and variables are writen.

  use system
  use defect
  use lmp_driver

  implicit none

  integer :: i, writeFreq

  character(len=200) filename, datapath

  writeFreq = 1000

  if (MOD(ISTEP, writeFreq) .EQ. 0) then
    !Write the atom pos and box size from LAMMPS
    call lmpSystem%writeData(dataPath)

    if (RANK .EQ. 0) then
      !write restart file with rest of simulation info
      write(filename, '(A17, I0, A4)') './RESTART/RESTART', ISTEP, '.OUT'
      open(unit=27, file=trim(filename))

      write(27, *) trim(datapath)
      write(27, *)

      write(27, *) &
        'ISTEP, TIME, TIME_LAST, DTIME, DPA, DPAC, AVDPA, ONDSTOT'
      write(27, *) ISTEP
      write(27, *) TIME
      write(27, *) TIME_LAST
      write(27, *) DTIME
      write(27, *) DPA
      write(27, *) DPAC
      write(27, *) AVDPA
      write(27, *) ONDSTOT
      write(27, *)

      !Write defect migration energies
      write (27, *) 'Vacancy Migration Energies'
      do i = 1, 50

        write (27, *) EMIG(1, I), EMIC(1, I)
      end do

      write (27, *) 'Interstitial Migration Energies'
      do i = 1, 50

        write (27, *) EMIG(2, I), EMIC(2, I)
      end do

      close (27)
    end if

  end if

end subroutine
!***********************************************************************!
subroutine DEFECT_PRINT(DEFECTLIST, NDEFECTS)

  USE SYSTEM, ONLY: NDEMAX

  implicit none

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
  do I = 1, NDEFECTS
    if (.NOT. ANY(INTSIZES == INT(DEFECTLIST(6, I))) .AND. INT(DEFECTLIST(6, I)) > 0) then
      INTSIZES(INTCOUNT + 1) = INT(DEFECTLIST(6, I))
      INTCOUNT = INTCOUNT + 1
    ELSE if (.NOT. ANY(VACSIZES == INT(DEFECTLIST(6, I))) .AND. INT(DEFECTLIST(6, I)) < 0) then
      VACSIZES(VACCOUNT + 1) = INT(DEFECTLIST(6, I))
      VACCOUNT = VACCOUNT + 1
    end if
  end do

  ! DETERMINE THE COUNT FOR EACH DEFECT SIZE
  do I = 1, NDEFECTS
    if (DEFECTLIST(6, I) > 0) then
      do J = 1, INTCOUNT
        if (INT(DEFECTLIST(6, I)) == INTSIZES(J)) then
          INTCLUSTERCOUNTS(J) = INTCLUSTERCOUNTS(J) + 1
        end if
      end do
    ELSE
      do J = 1, VACCOUNT
        if (INT(DEFECTLIST(6, I)) == VACSIZES(J)) then
          VACCLUSTERCOUNTS(J) = VACCLUSTERCOUNTS(J) + 1
        end if
      end do
    end if
  end do

  ! SORT THE DEFECT LISTS
  CALL QUICKSORT_INT(INTSIZES, INTCLUSTERCOUNTS, 1, INTCOUNT)
  CALL QUICKSORT_INT(VACSIZES, VACCLUSTERCOUNTS, 1, VACCOUNT)

  PRINTFORMAT = '(X,I4,X)'

  !write THE INT TABLE
  if (INTCOUNT > 0) then

    write (*, '(A9)', ADVANCE="NO") 'Int Size|'
    do I = 1, INTCOUNT
      write (*, PRINTFORMAT, ADVANCE="NO") INTSIZES(I)
    end do
    write (*, *) !NEWLINE

    write (*, '(A9)', ADVANCE="NO") 'Count|'
    do I = 1, INTCOUNT
      write (*, PRINTFORMAT, ADVANCE="NO") INTCLUSTERCOUNTS(I)
    end do
    write (*, *) !NEWLINE
    write (*, *) !NEWLINE

  end if

  !write THE VAC TABLE
  if (VACCOUNT > 0) then

    write (*, '(A9)', ADVANCE="NO") 'Vac Size|'
    do I = VACCOUNT, 1, -1
      write (*, PRINTFORMAT, ADVANCE="NO") VACSIZES(I)
    end do
    write (*, *) !NEWLINE

    write (*, '(A9)', ADVANCE="NO") 'Count|'
    do I = VACCOUNT, 1, -1
      write (*, PRINTFORMAT, ADVANCE="NO") VACCLUSTERCOUNTS(I)
    end do
    write (*, *) !NEWLINE
    write (*, *) !NEWLINE

  end if
end subroutine DEFECT_PRINT