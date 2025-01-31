SUBROUTINE OPLD

  use system
  use lmp_driver
  use MicroStrucFeatures

  IMPLICIT NONE

  if (BSTEP == 0) then
    call lmpSystem%relaxSystemPressure()
    call lmpSystem%getPositions()
    call features%grains%init(lmpSystem, MGSZ)
    call features%grains%identifyGrains()
    call features%reflat%createRefLattice(features%grains)
  end if

  do ISTEP = (BSTEP + 1), (NTSTEP + BSTEP)

    call DOSE()

    call FIND_DEFECT()

    if (RANK .EQ. 0) THEN
      call DCIDNT()
      call DFSEL()
    end if

    call SPS()
    call KMC()
    if (.NOT. SKIP) call OPT()
    call RESTART()

  end do

end SUBROUTINE
!***********************************************************************!
subroutine SPS

  use MPI
  use SYSTEM
  use DEFECT
  use KMCVAR

  use CLASS_DIMER
  use LMP_DRIVER
  use MicroStrucFeatures

  implicit none

  type(Dimer) :: ADIMER

  integer :: J, K, SR_LEN
  real*8 :: x,y,z,r

  SR = 0.0D0; if (RANK .EQ. 0) SRBUFF = 0.0D0
  EM = 0.0D0; if (RANK .EQ. 0) EMBUFF = 0.0D0
  ND = 0.0D0; if (RANK .EQ. 0) NDBUFF = 0.0D0
  SR_LEN = SIZE(SR)

  call features%defects%syncClusterData(lmpSystem%getCommunicator())
  x = features%defects%selectedCluster%comx
  y = features%defects%selectedCluster%comy
  z = features%defects%selectedCluster%comz
  r = features%defects%selectedCluster%radiusAV

  !CARRY OUT SADDLE POINT SEARCHING
  do J = LDIMERI, UDIMERI

    call DIMER_INIT(ADIMER, x, y, z, r, lmpDimer)
    call DIMER_SEARCH(ADIMER)

    if (lmpDimer%getRank() .EQ. 0) then
      ND(J) = ADIMER%NATOM1
      EM(J) = ADIMER%EMDIFF
    
      !STORE THE SADDLE POINT CONFIGURATION,SR:POSITION OF SADDLE POINT
      do K = 1, ND(J)

        SR(1, K, J) = ADIMER%SX(1, K)
        SR(2, K, J) = ADIMER%SX(2, K)
        SR(3, K, J) = ADIMER%SX(3, K)
        SR(4, K, J) = ADIMER%SI(1, K)
      end do
    end if

    if (mod(ISTEP, NFREQ) .eq. 0 .and. lmpDimer%getRank() .eq. 0) then
      write (*, '(A,1I12,F14.5)') 'SPS', J, EM(J)
    end if

    call DIMER_FINISH(ADIMER)
  end do

  !Gather data on main rank
  call MPI_REDUCE(SR, SRBUFF, SR_LEN, MPI_REAL8, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)
  call MPI_REDUCE(EM, EMBUFF, NDIMER, MPI_REAL8, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)
  call MPI_REDUCE(ND, NDBUFF, NDIMER, MPI_INT, MPI_SUM, 0, lmpSystem%getCommunicator(), IERR)

  if (rank .eq. 0) then
    SR = SRBUFF
    EM = EMBUFF
    ND = NDBUFF
  end if

end subroutine
!***********************************************************************!
SUBROUTINE OPT

  USE SYSTEM
  USE KMCVAR
  USE DEFECT
  USE LMP_DRIVER
  USE MPI

  IMPLICIT NONE

  INTEGER I, J, NSTMP

  if (RANK .EQ. 0) THEN

    !THE TOTAL NUMBER OF THE ATOMS INSIDE THE CHOSEN AV

    NSTMP = ND(EID)

    !PASSING THE POSITIONS OF ATOMS INSIDE AV TO THE SYSTEM

    do I = 1, NSTMP

      J = INT(SR(4, I, EID))

      RA(1, J) = SR(1, I, EID)
      RA(2, J) = SR(2, I, EID)
      RA(3, J) = SR(3, I, EID)

    end do

  end if

  call lmpSystem%assignPositions

  call lmpSystem%relaxSystemEnergy

  call lmpSystem%getPositions

  if (RANK .EQ. 0) THEN

    !ISTEP:CURRENT NUMBER OF STEPS,NDTOT:TOTAL NUMBER OF DEFECTS
    !JID:DEFECT ID,EID:EVENT ID,EM:MIGRATION ENERGY BARRIER
    !ENCALC:THE SYSTEM ENERGY AFTER OPT,TIME:CURRENT TIME

    if (MOD(ISTEP, NFREQ) .EQ. 0) THEN
      write (*, '(3X, 2A12, 7A14)') 'STEP', 'DIMER ID', 'BARRIER', &
        'SYS ENERGY', 'TIME', 'Px', 'Py', 'Pz', 'DPA'

      write (*, '(A,2I12,2F14.5,4ES14.5,E14.5)') 'ITR', ISTEP, EID, EM(EID), &
        ENCALC, TIME, PRESS(1), PRESS(2), PRESS(3), DPA

      write(*,*)

    end if

  end if

  !OUTPUT DEFECTS FOR OVITO

  call FIND_DEFECT

  if (MOD(ISTEP, NFREQ) .EQ. 0 .AND. RANK .EQ. 0) call OVITO

  RETURN

end SUBROUTINE OPT

!***********************************************************************!
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
    call flush(6)
  end if

  NDSTOT = max(sumInt, sumVac)
  if (CALC) call MD_DPA
  ONDTOT = NDTOT
  ONDSTOT = NDSTOT

  if (NDTOT .EQ. 0) then
    write (*, *) 'No defects in system. Stopping.'
    stop
  end if

end SUBROUTINE DCIDNT
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

  call RANDOM_NUMBER(RND)

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

  x = features%defects%selectedCluster%comx
  y = features%defects%selectedCluster%comy
  z = features%defects%selectedCluster%comz
  r = features%defects%selectedCluster%radiusAV
  clusterType = features%defects%selectedCluster%type

  if (mod(ISTEP, NFREQ) .eq. 0) then
    write (*, '(15X,5(A10,X))') 'X', 'Y', 'Z', 'RADIUS', 'TYPE'
    write (*, '(A15,4(F10.3,X),I10)') 'Selected Defect', x, y, z, r, clusterType
  end if

  deallocate (pj)

end SUBROUTINE DFSEL

!***********************************************************************!
SUBROUTINE KMC

  use SYSTEM
  use KMCVAR
  use DEFECT
  use MPI
  use LMP_DRIVER
  use MicroStrucFeatures

  implicit none

  integer :: I, K, jmapsz, type

  real*8 RND, dt, scalef, prop, TOTFQSUM, x, y, z

  logical ban

  character(LEN=40) str

  if (rank .eq. 0) then

!FQSUM:SUM OF FREQUENCY,POSUM:SUM OF PROBABILITIES(FOR CHECK ONLY)
    FQSUM = 0.0D0
    TOTFQSUM = 0.0D0
    POSUM = 0.0D0
    jmapsz = 1000

!CALCULATE FREQUENCY BASED ON THE MIGRATION ENERGY BARRIER
!THE PREFACTOR HERE IS USED TO BE CONSTANT (NOT TRUE IN GENERAL)
!THE PREFACTOR CALCUATED ON THE FLY IS IMPLEMENTED IN ANOTHER VERSION
! Now considering good dimer search - Gao 2021

    do K = 1, NDIMER

      if (EM(k) .lt. 10) then
        FQ(K) = 1.0D12*EXP(-EM(K)/KB/TEMP)
        FQSUM = FQSUM + FQ(K)
      end if
    end do

    call RANDOM_NUMBER(RND)
    prop = rnd*float(jmapsz)

    !DETERMINE TIME STEP IN KMC
    if (FQSUM .gt. 0.0) then
      TMEAN = -LOG10(RND)/(FQSUM+OBFQSUM) !With time correction using other objects FQs
      dt = 1.0/fqsum
      scalef = dt*dfloat(jmapsz)
    end if

    do K = 1, NDIMER

      if (EM(k) .lt. 10.0) then
        POSUM = POSUM + FQ(K)*scalef
      end if
      if (POSUM .ge. prop) exit
    end do

    !EID:THE ID OF THE EVENT THAT IS CHOSEN

    EID = K

    !UPDATING THE PRE-SELECTION ENERGIES
    x = features%defects%selectedCluster%comx
    y = features%defects%selectedCluster%comy
    z = features%defects%selectedCluster%comz
    type = features%defects%selectedCluster%type
    if (type .lt. 0) then
      K = 1
    else
      K = 2
    end if

    !CHECKING FOR TOO HIGH TRANSITION ENERGY

    if (EM(EID) .LE. EMSKIP .AND. EID .LE. NDIMER) THEN
      TIME = TIME + TMEAN

      !ONLY POSITIVE SADDLE POINT ENERGIES WILL BE USED IN AVERAGE
      if (EM(EID) .GT. 0) THEN
        I = abs(type)
        EMIG(K, I) = (EMIG(K, I)*EMIC(K, I) + EM(EID))/(EMIC(K, I) + 1)
        EMIC(K, I) = EMIC(K, I) + 1
      end if
    else
      SKIP = .TRUE.
      write (STR, '(F8.5)') EMSKIP
      write (6, *) 'MIGRATION ENERGY GREATER THAN '//TRIM(STR)// &
        'eV -> SKIPPING'
    end if

    if (MOD(ISTEP, NFREQ) .EQ. 0) THEN
      write (6, '(A,I10,A,F10.4)') 'EID: ', EID, '  RND: ', RND
    end if

  end if

  ! If defect had no significant energy barriers ban that defect
  ban = .true.
  do i = 1, NDIMER
    if (EM(i) > 1e-3 .and. EM(i) < 10.0D0) ban = .false.
  end do

  call MPI_BCAST(ban, 1, MPI_LOGICAL, 0, lmpSystem%getCommunicator(), ierr)

  if (ban) then
    call features%defects%banCluster(features%defects%selectedCluster)
    if (rank == 0) write(*,'(A,3F14.5)') 'Banned defect at: ', x, y, z
  end if

  ! Update important data to all procs
  call MPI_BCAST(TIME, 1, MPI_REAL8, 0, lmpSystem%getCommunicator(), IERR)
  call MPI_BCAST(DTIME, 1, MPI_REAL8, 0, lmpSystem%getCommunicator(), IERR)
  call MPI_BCAST(SKIP, 1, MPI_LOGICAL, 0, lmpSystem%getCommunicator(), IERR)

  RETURN

end SUBROUTINE KMC
!***********************************************************************!
SUBROUTINE DOSE

  USE SYSTEM
  USE LMP_DRIVER
  use MicroStrucFeatures

  IMPLICIT NONE

  if (TIME - TIME_LAST .GE. DTIME .OR. ISTEP .EQ. 1) THEN

    if (RANK .EQ. 0) THEN
      write (*, '(3A14)') 'TIME', ' TIME_LAST', ' DTIME'
      write (*, '(2E14.5,E14.5)') TIME, TIME_LAST, DTIME
    end if

    if (TDAM .EQ. 2 .OR. TDAM .EQ. 3) THEN

      call BASIC_DAMAGE(NFREKL)

      call lmpSystem%assignPositions

    end if

    call lmpSystem%relaxSystemEnergy

    call lmpSystem%getPositions

    TIME_LAST = TIME

    if (RANK == 0) CALC = .TRUE.

    call features%defects%clearBannedList()

  end if

end SUBROUTINE
!***********************************************************************!
