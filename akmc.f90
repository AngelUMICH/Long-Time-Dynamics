SUBROUTINE OPLD

  USE SYSTEM, ONLY: BSTEP, NTSTEP, ISTEP, RANK, SKIP
  USE DEFECT, ONLY: UPNEIGH

  IMPLICIT NONE

  DO ISTEP = (BSTEP + 1), (NTSTEP + BSTEP)

    CALL DOSE

    CALL FIND_DEFECT

    IF (RANK .EQ. 0) THEN
      CALL DCIDNT
      CALL DFSEL
    END IF

    CALL SPS
    CALL KMC
    IF (.NOT. SKIP) CALL OPT
    CALL RESTART

  END DO

END SUBROUTINE
!***********************************************************************!
!THIS SUBROUTINE IS TO PERFORM SADDLE POINT SEARCHING

SUBROUTINE SPS

  USE MPI
  USE SYSTEM
  USE DEFECT
  USE KMCVAR

  USE CLASS_DIMER
  USE LMP_DRIVER
  use MicroStrucFeatures

  IMPLICIT NONE

  TYPE(DIMER) :: ADIMER

  INTEGER :: J, K

  INTEGER :: RA_LEN, SR_LEN

  !CLEAR ALL DATA

  SR = 0.0D0
  IF (RANK .EQ. 0) SRBUFF = 0.0D0
  EM = 0.0D0
  IF (RANK .EQ. 0) EMBUFF = 0.0D0
  ND = 0.0D0
  IF (RANK .EQ. 0) NDBUFF = 0.0D0

  !MPI VARIABLES

  RA_LEN = SIZE(RA)

  SR_LEN = SIZE(SR)

  !TODO Eliminate obtaining info from features to put in px
  if (rank .eq. 0) then
    px(1,1) = features%defects%selectedCluster%comx
    px(2,1) = features%defects%selectedCluster%comy
    px(3,1) = features%defects%selectedCluster%comz
    px(4,1) = features%defects%selectedCluster%radiusAV
  end if

  CALL MPI_BCAST(RA, RA_LEN, MPI_REAL8, 0, lmpSystem%getCommunicator(), IERR)

  CALL MPI_BCAST(NATOMC, 1, MPI_INT, 0, lmpSystem%getCommunicator(), IERR)

  CALL MPI_BCAST(PX, 7*NDEMAX, MPI_REAL8, 0, lmpSystem%getCommunicator(), IERR)

  !CARRY OUT SADDLE POINT SEARCHING

  DO J = LDIMERI, UDIMERI

    CALL DIMER_INIT(ADIMER, PX(1, 1), PX(2, 1), PX(3, 1), PX(4, 1), lmpDimer)

    CALL DIMER_SEARCH(ADIMER)

    !STORE INFORMATION,EM:MIGRATION ENERGY BARRIER
    !ND:TOTAL NUMBER OF ATOMS IN THE AV

    IF (lmpDimer%getRank() .EQ. 0) THEN
      ND(J) = ADIMER%NATOM1
      EM(J) = ADIMER%EMDIFF
    

      !STORE THE SADDLE POINT CONFIGURATION,SR:POSITION OF SADDLE POINT

      DO K = 1, ND(J)

        SR(1, K, J) = ADIMER%SX(1, K)
        SR(2, K, J) = ADIMER%SX(2, K)
        SR(3, K, J) = ADIMER%SX(3, K)
        SR(4, K, J) = ADIMER%SI(1, K)

      END DO

    END IF

    IF (MOD(ISTEP, NFREQ) .EQ. 0 .AND. lmpDimer%getRank() .EQ. 0) THEN
      WRITE (*, '(A,1I12,F14.5)') 'SPS', J, EM(J)
    END IF

    CALL DIMER_FINISH(ADIMER)

  END DO

  !GATHER RESULTS ON MAIN RANK
  CALL MPI_REDUCE(SR, SRBUFF, SR_LEN, MPI_REAL8, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)
  CALL MPI_REDUCE(EM, EMBUFF, NDIMER, MPI_REAL8, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)
  CALL MPI_REDUCE(ND, NDBUFF, NDIMER, MPI_INT, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)

  IF (RANK .EQ. 0) THEN

    SR = SRBUFF

    EM = EMBUFF

    ND = NDBUFF

  END IF

  RETURN
END SUBROUTINE
!***********************************************************************!

!STATIC RELAXATION OF THE SYSTEM

SUBROUTINE OPT

  USE SYSTEM
  USE KMCVAR
  USE DEFECT
  USE LMP_DRIVER
  USE MPI

  IMPLICIT NONE

  INTEGER I, J, NSTMP

  IF (RANK .EQ. 0) THEN

    !THE TOTAL NUMBER OF THE ATOMS INSIDE THE CHOSEN AV

    NSTMP = ND(EID)

    !PASSING THE POSITIONS OF ATOMS INSIDE AV TO THE SYSTEM

    DO I = 1, NSTMP

      J = INT(SR(4, I, EID))

      RA(1, J) = SR(1, I, EID)
      RA(2, J) = SR(2, I, EID)
      RA(3, J) = SR(3, I, EID)

    END DO

  END IF

  call lmpSystem%assignPositions

  call lmpSystem%relaxSystemEnergy

  call lmpSystem%getPositions

  IF (RANK .EQ. 0) THEN

    !ISTEP:CURRENT NUMBER OF STEPS,NDTOT:TOTAL NUMBER OF DEFECTS
    !JID:DEFECT ID,EID:EVENT ID,EM:MIGRATION ENERGY BARRIER
    !ENCALC:THE SYSTEM ENERGY AFTER OPT,TIME:CURRENT TIME

    IF (MOD(ISTEP, NFREQ) .EQ. 0) THEN
      WRITE (*, '(3X, 2A12, 7A14)') 'STEP', 'DIMER ID', 'BARRIER', &
        'SYS ENERGY', 'TIME', 'Px', 'Py', 'Pz', 'DPA'

      WRITE (*, '(A,2I12,2F14.5,4E14.5,E14.5)') 'ITR', ISTEP, EID, EM(EID), &
        ENCALC, TIME, PRESS(1), PRESS(2), PRESS(3), DPA

      WRITE(*,*)

    END IF

  END IF

  !OUTPUT DEFECTS FOR OVITO

  CALL FIND_DEFECT

  IF (MOD(ISTEP, NFREQ) .EQ. 0 .AND. RANK .EQ. 0) CALL OVITO

  RETURN

END SUBROUTINE OPT

!***********************************************************************!
!THIS SUBROUTINE IS TO IDENTIFY DEFECT CLUSTERS IN THE SYSTEM

SUBROUTINE DCIDNT
  use SYSTEM
  use DEFECT
  use MicroStrucFeatures

  implicit none

  integer i,j
  integer :: sumInt, sumVac

  if (ISTEP .eq. 1) then
    ISTEP = 0
    call OVITO
    ISTEP = 1
  end if

  NDVAC = 0
  NDINT = 0
  II = 0
  IV = 0

!IDENTIFY INTERSTITIAL,PI:POSITION OF INTERSITITAL

  do i = 1, NATOMC

    if (TI(i) .ne. 1) then
      NDINT = NDINT + 1
      II(1, NDINT) = i
      PI(1, NDINT) = RA(1, i)
      PI(2, NDINT) = RA(2, i)
      PI(3, NDINT) = RA(3, i)
    end if
  end do

!IDENTIFY VACANCY,PV:POSITION OF VACANCY

  do j = 1, NATOMR

    if (TV(j) .ne. 1) then
      NDVAC = NDVAC + 1
      IV(1, NDVAC) = J
      PV(1, NDVAC) = RR(1, J)
      PV(2, NDVAC) = RR(2, J)
      PV(3, NDVAC) = RR(3, J)
    end if
  end do

  NDTOT = 0

  do i = 1, NDINT

    NDTOT = NDTOT + 1
    PX(1, NDTOT) = PI(1, i)
    PX(2, NDTOT) = PI(2, i)
    PX(3, NDTOT) = PI(3, i)
    PX(6, NDTOT) = 1.0D0  !interstitial identifier
    PX(7, NDTOT) = II(1, i)
  end do

  do j = 1, NDVAC

    NDTOT = NDTOT + 1
    PX(1, NDTOT) = PV(1, j)
    PX(2, NDTOT) = PV(2, j)
    PX(3, NDTOT) = PV(3, j)
    PX(6, NDTOT) = -1.0D0  !vacancy identifier
    PX(7, NDTOT) = IV(1, j)
  end do
  
  call features%defects%identifyClusters(PX, NDTOT, 1.0D0)

  NDTOT = features%defects%numClusters

  !Put the clusters data in the PX array
  !TODO remove dependancy of PX array
  do i = 1, NDTOT

    px(1,i) = features%defects%clusters(i)%comx
    px(2,i) = features%defects%clusters(i)%comy
    px(3,i) = features%defects%clusters(i)%comz
    px(4,i) = features%defects%clusters(i)%radiusAV
    px(5,i) = 1.0D0
    px(6,i) = features%defects%clusters(i)%type
    px(7,i) = 0
  end do

  NDSTOT = 0
  sumInt = 0; sumVac = 0

  do i = 1, NDTOT

    if (PX(6, i) .GT. 0) sumInt = sumInt + int(abs(PX(6, i)))
    if (PX(6, i) .LT. 0) sumVac = sumVac + int(abs(PX(6, i)))
  end do

  if (mod(ISTEP, NFREQ) .eq. 0) then
    call DEFECT_PRINT(PX, NDTOT)
    call features%defects%printBannedCount()
  end if

  call flush(6)

  NDSTOT = max(sumInt, sumVac)

  if (CALC) call MD_DPA

  ONDTOT = NDTOT

  ONDSTOT = NDSTOT

!IF NO DEFECT IS PRESENT, STOP THE PROGRM

  if (NDTOT .EQ. 0) then
    write (*, *) 'No defects in system. Stopping.'
    stop
  end if

END SUBROUTINE DCIDNT
!***********************************************************************!
!************************************************************
! Idea is to decrease dimer search for defects using KMC    *
! to select a defect for dimer search - migration energy    *
! was pre-calculated one                                    *
!                                                           *
! Gao - 2021                                                *
!************************************************************
SUBROUTINE DFSEL
  USE SYSTEM
  USE DEFECT
  USE KMCVAR
  use MicroStrucFeatures

  IMPLICIT NONE

  REAL*8 PTEMP(7, 1)
  REAL*8 rnd, speed, dt, scalef, sump, prop
  real*8 :: x,y,z,r
  integer :: i, j, k, jmapsz, numClusters
  integer :: clusterType
  REAL*8, allocatable, dimension(:, :) :: pj

  SKIP = .FALSE.

  speed = 0.0d0
  jmapsz = 1000
  numClusters = features%defects%numClusters

  allocate (pj(1:2, 1:ndtot))

  do i = 1, numClusters

    if (features%defects%clusters(i)%isBanned()) cycle

    clusterType = features%defects%clusters(i)%type
    j = abs(clusterType)

    if (clusterType .lt. 0) then
      k = 1 !vacancy type
    else
      k = 2 !interstitial type
    end if

    if (j .lt. 50) then
      pj(k, i) = 1.0D12*exp(-emig(k, j)/kb/temp)
    else
      pj(k, i) = 1.0D12*exp(-emig(k, 50)/kb/temp)
    end if

    speed = speed + pj(k, i)
  end do

  dt = 1.0/speed
  scalef = dt*float(jmapsz)

  CALL RANDOM_NUMBER(RND)

  prop = rnd*float(jmapsz)

  sump = 0.0

  do i = 1, ndtot

    if (features%defects%clusters(i)%isBanned()) cycle

    clusterType = features%defects%clusters(i)%type

    if (clusterType .lt. 0) then
      k = 1 !vacancy type
    else
      k = 2 !interstitial type
    end if

    sump = sump + pj(k, i)*scalef
    if (sump .ge. prop) exit
  end do

  call features%defects%saveClusterSelection(features%defects%clusters(i))

  OBFQSUM = SPEED - pj(k, i)*scalef !rate sum of all other non chosen clusters

  PX = 0.0D0

!PLACE SELECTED DEFECT IN THE FIRST POSITION

  if (mod(ISTEP, NFREQ) .eq. 0) then
    x = features%defects%selectedCluster%comx
    y = features%defects%selectedCluster%comy
    z = features%defects%selectedCluster%comz
    r = features%defects%selectedCluster%radiusAV
    clusterType = features%defects%selectedCluster%type

    write (*, '(15X,5(A10,X))') 'X', 'Y', 'Z', 'RADIUS', 'TYPE'
    write (*, '(A15,4(F10.3,X),I10)') 'Selected Defect', x, y, z, r, clusterType
  end if

  deallocate (pj)

END SUBROUTINE DFSEL

!***********************************************************************!
!THIS SUBROUTINE IS TO CHOSE AV BASED ON THE MIRGRATION ENERGY BARRIERS
!THE CONVENTIONAL KINETIC MONTE CARLO ALGORITHM IS USED

SUBROUTINE KMC

  USE SYSTEM
  USE KMCVAR
  USE DEFECT
  USE MPI
  USE LMP_DRIVER
  USE MicroStrucFeatures

  IMPLICIT NONE

  INTEGER I, K, jmapsz

  REAL*8 RND, dt, scalef, prop, TOTFQSUM

  logical ban

  CHARACTER(LEN=40) STR

  IF (RANK .EQ. 0) THEN

!FQSUM:SUM OF FREQUENCY,POSUM:SUM OF PROBABILITIES(FOR CHECK ONLY)

    FQSUM = 0.0D0
    TOTFQSUM = 0.0D0
    POSUM = 0.0D0
    jmapsz = 1000

!CALCULATE FREQUENCY BASED ON THE MIGRATION ENERGY BARRIER
!THE PREFACTOR HERE IS USED TO BE CONSTANT (NOT TRUE IN GENERAL)
!THE PREFACTOR CALCUATED ON THE FLY IS IMPLEMENTED IN ANOTHER VERSION
! Now considering good dimer search - Gao 2021

    DO K = 1, NDIMER
      if (EM(k) .lt. 10) then
        FQ(K) = 1.0D12*EXP(-EM(K)/KB/TEMP)
        FQSUM = FQSUM + FQ(K)
      end if
    END DO

    CALL RANDOM_NUMBER(RND)
    prop = rnd*float(jmapsz)

    !DETERMINE TIME STEP IN KMC

    if (FQSUM .gt. 0.0) then
      TMEAN = -LOG10(RND)/(FQSUM+OBFQSUM) !With time correction using other objects FQs

      dt = 1.0/fqsum
      scalef = dt*dfloat(jmapsz)
    end if

    DO K = 1, NDIMER
      if (EM(k) .lt. 10.0) then
        POSUM = POSUM + FQ(K)*scalef
      end if
      if (POSUM .ge. prop) exit
    END DO

    !EID:THE ID OF THE EVENT THAT IS CHOSEN

    EID = K

    !UPDATING THE PRE-SELECTION ENERGIES

    IF (PX(6, 1) .LT. 0) THEN
      K = 1
    ELSE
      K = 2
    END IF

    !CHECKING FOR TOO HIGH TRANSITION ENERGY

    IF (EM(EID) .LE. EMSKIP .AND. EID .LE. NDIMER) THEN

      TIME = TIME + TMEAN

      !ONLY POSITIVE SADDLE POINT ENERGIES WILL BE USED IN AVERAGE

      IF (EM(EID) .GT. 0) THEN

        I = ABS(INT(PX(6, 1)))

        EMIG(K, I) = (EMIG(K, I)*EMIC(K, I) + EM(EID))/(EMIC(K, I) + 1)

        EMIC(K, I) = EMIC(K, I) + 1

      END IF

    ELSE

      SKIP = .TRUE.

      WRITE (STR, '(F8.5)') EMSKIP

      WRITE (6, *) 'MIGRATION ENERGY GREATER THAN '//TRIM(STR)// &
        'eV -> SKIPPING'

    END IF

    IF (MOD(ISTEP, NFREQ) .EQ. 0) THEN
      WRITE (6, '(A,I10,A,F10.4)') 'EID: ', EID, '  RND: ', RND
    END IF

  END IF

  ! If defect had no significant energy barriers ban that defect
  ban = .true.
  do i = 1, NDIMER
    if (EM(i) > 1e-3 .and. EM(i) < 10.0D0) ban = .false.
  end do

  call MPI_BCAST(ban, 1, MPI_LOGICAL, 0, lmpSystem%getCommunicator(), ierr)

  if (ban) then
    call features%defects%banCluster(features%defects%selectedCluster)
    if (rank == 0) write(*,'(A,3F14.5)') 'Banned defect at: ', PX(1,1),PX(2,1),PX(3,1)
  end if

  ! Update important data to all procs
  CALL MPI_BCAST(TIME, 1, MPI_REAL8, 0, lmpSystem%getCommunicator(), IERR)

  CALL MPI_BCAST(DTIME, 1, MPI_REAL8, 0, lmpSystem%getCommunicator(), IERR)

  CALL MPI_BCAST(SKIP, 1, MPI_LOGICAL, 0, lmpSystem%getCommunicator(), IERR)

  RETURN

END SUBROUTINE KMC
!***********************************************************************!
!WILL INTRODUCE MORE DEFECTS AFTER A SPECIFIED AMOUNT OF TIME THAT MATCHES
!EXPERIMENTAL DOSE RATE
SUBROUTINE DOSE

  USE SYSTEM
  USE LMP_DRIVER
  use MicroStrucFeatures

  IMPLICIT NONE

  IF (TIME - TIME_LAST .GE. DTIME .OR. ISTEP .EQ. 1) THEN

    IF (RANK .EQ. 0) THEN

      WRITE (*, '(3A14)') 'TIME', ' TIME_LAST', ' DTIME'

      WRITE (*, '(2E14.5,E14.5)') TIME, TIME_LAST, DTIME

    END IF

    call lmpSystem%relaxSystemPressure

    call lmpSystem%getPositions

    call lmpSystem%updateRef

    IF (TDAM .EQ. 2 .OR. TDAM .EQ. 3) THEN ! Frenkel or Vacancy

      CALL BASIC_DAMAGE(NFREKL)

      call lmpSystem%assignPositions

    END IF

    call lmpSystem%relaxSystemEnergy

    call lmpSystem%getPositions

    TIME_LAST = TIME

    IF (RANK == 0) CALC = .TRUE.

    call features%defects%clearBannedList()

  END IF

END SUBROUTINE
!***********************************************************************!
