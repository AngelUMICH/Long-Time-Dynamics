subroutine drld

  use system
  use lmp_driver
  
  implicit none
  
  if (BSTEP == 0) then
    call lmpSystem%relaxSystemPressure()
    call lmpSystem%getPositions()
    call lmpDimer%changeBoxDimensions(boxx, boxy, boxz)
  end if
  
  do ISTEP = (BSTEP + 1), (NTSTEP + BSTEP)

    call DOSE()
    call FIND_DEFECT()

    if (RANK .EQ. 0) THEN
      call DCIDNT()
      call DFSEL()
    end if

    call BEM()
    call RESTART()
    call FIND_DEFECT()
    call OVITO()
    
  end do
  
end subroutine drld
!***********************************************************************!
subroutine DOSE

  use SYSTEM
  use LMP_DRIVER
  use Clusters

  implicit none

  if (TIME - TIME_LAST .GE. DTIME .OR. ISTEP .EQ. 1) THEN

    if (RANK .EQ. 0) THEN
      write (*, '(3A14)') 'TIME', ' TIME_LAST', ' DTIME'
      write (*, '(2E14.5,E14.5)') TIME, TIME_LAST, DTIME
    end if

    if (TDAM .EQ. 2 .OR. TDAM .EQ. 3) THEN
      call BASIC_DAMAGE(NFREKL)
      call lmpSystem%assignPositions
    end if

    if (NFREKL > 0) then
      call lmpSystem%relaxSystemEnergy
      call lmpSystem%getPositions
    end if

    TIME_LAST = TIME

    if (RANK == 0) CALC = .TRUE.

  end if

end subroutine DOSE
!***********************************************************************!
subroutine DCIDNT
  use SYSTEM
  use DEFECT
  use Clusters

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
  
  print *, 'Identifying defect clusters'
  call defects%identifyClusters(PX, NDTOT, 1.2D0)

  NDTOT = defects%numClusters

  !Put the clusters data in the PX array
  !TODO remove dependancy of PX array
  do i = 1, NDTOT

    px(1,i) = defects%clusters(i)%comx
    px(2,i) = defects%clusters(i)%comy
    px(3,i) = defects%clusters(i)%comz
    px(4,i) = defects%clusters(i)%radiusAV
    px(5,i) = 1.0D0
    px(6,i) = defects%clusters(i)%type
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

end subroutine DCIDNT
!***********************************************************************!
subroutine DFSEL
  !************************************************************
  ! Idea is to decrease dimer search for defects using KMC    *
  ! to select a defect for dimer search - migration energy    *
  ! was pre-calculated one                                    *
  !                                                           *
  ! Gao - 2021                                                *
  !************************************************************
  use SYSTEM
  use DEFECT
  use KMCVAR
  use Clusters

  implicit none

  REAL*8 rnd, speed, dt, scalef, sump, prop
  real*8 :: x,y,z,r
  integer :: i, j, k, jmapsz, numClusters
  integer :: clusterType
  REAL*8, allocatable, dimension(:, :) :: pj

  speed = 0.0d0
  jmapsz = 1000
  numClusters = defects%numClusters

  allocate (pj(1:2, 1:ndtot))

  do i = 1, numClusters

    clusterType = defects%clusters(i)%type
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

    clusterType = defects%clusters(i)%type

    if (clusterType .lt. 0) then
      k = 1 !vacancy type
    else
      k = 2 !interstitial type
    end if

    sump = sump + pj(k, i)*scalef
    if (sump .ge. prop) exit
  end do

  call defects%saveClusterSelection(defects%clusters(i))

  OBFQSUM = SPEED - pj(k, i) !rate sum of all other non chosen clusters

  x = defects%selectedCluster%comx
  y = defects%selectedCluster%comy
  z = defects%selectedCluster%comz
  r = defects%selectedCluster%radiusAV
  clusterType = defects%selectedCluster%type

  if (mod(ISTEP, NFREQ) .eq. 0) then
    write (*, '(15X,5(A10,X))') 'X', 'Y', 'Z', 'RADIUS', 'TYPE'
    write (*, '(A15,4(F10.3,X),I10)') 'Selected Defect', x, y, z, r, clusterType
  end if

  deallocate (pj)

end subroutine DFSEL
!***********************************************************************!
subroutine BEM
  ! Basin escape method (BEM) uses MRM to quickly leave basins.

  use BasinNetwork, only: BasinNetworkType
  use Event, only: EventType
  use Basin, only: BasinType
  
  type(BasinNetworkType) :: network !resets every call of the subroutine
  type(EventType) :: event
  type(BasinType) :: basin
  logical :: continueSearch

  continueSearch = .True.

  do while (continueSearch)
    call SPS()
    call EBN(network, continueSearch)
    call KMC(network, continueSearch)
    basin = network%getBasin(network%getActiveBasinId())
    event = basin%getAbsorbingEvent(network%getActiveEventId())
    call OPT(event, continueSearch)
  end do

end subroutine BEM
!***********************************************************************!
subroutine SPS()

  use MPI
  use SYSTEM
  use DEFECT
  use KMCVAR

  use CLASS_DIMER
  use LMP_DRIVER
  use Clusters

  implicit none

  type(Dimer) :: ADIMER

  integer :: J, K, SR_LEN
  real*8 :: x,y,z,r

  SR = 0.0D0; SRBUFF = 0.0D0
  EM = 0.0D0; EMBUFF = 0.0D0
  ND = 0.0D0; NDBUFF = 0.0D0
  SR_LEN = SIZE(SR)

  call defects%syncClusterData(lmpSystem%getCommunicator())
  x = defects%selectedCluster%comx
  y = defects%selectedCluster%comy
  z = defects%selectedCluster%comz
  r = defects%selectedCluster%radiusAV

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

  call MPI_ALLREDUCE(SR, SRBUFF, SR_LEN, MPI_REAL8, MPI_SUM, lmpSystem%getCommunicator(), IERR)
  call MPI_ALLREDUCE(EM, EMBUFF, NDIMER, MPI_REAL8, MPI_SUM, lmpSystem%getCommunicator(), IERR)
  call MPI_ALLREDUCE(ND, NDBUFF, NDIMER, MPI_INT,   MPI_SUM, lmpSystem%getCommunicator(), IERR)

  SR = SRBUFF
  EM = EMBUFF
  ND = NDBUFF

end subroutine SPS
!***********************************************************************!
subroutine EBN(net, continueSearch)
  ! Expand basin network (EBN) with new found saddle points.
  ! Put the events into a basin and then into the network.

  use SYSTEM
  use DEFECT
  use KMCVAR
  use BasinNetwork
  use Basin
  use Event

  type(BasinNetworkType), intent(inout) :: net
  logical, intent(in) :: continueSearch

  integer :: i, j
  integer :: atomCount, atomId, numEvents
  integer :: matchingBasinId, activeBasinId, activeEventId, latestBasinId
  logical :: isUnique
  integer, allocatable :: atomIds(:)
  real*8, allocatable :: atomPositions(:,:)

  type(BasinType) :: basin, tempBasin
  type(EventType) :: event

  ! Set basin state using atom IDs and positions
  
  atomCount = ND(1)
  allocate(atomIds(atomCount))
  allocate(atomPositions(3, atomCount))

  do i = 1, atomCount
    atomId = int(SR(4, i, 1))
    atomIds(i) = atomId
    atomPositions(1, i) = RA(1, atomId)
    atomPositions(2, i) = RA(2, atomId)
    atomPositions(3, i) = RA(3, atomId)
  end do
  
  call basin%setState(atomIds, atomPositions)
  
  ! Now define the events and add them to the basin
  numEvents = NDIMER
  
  do i = 1, numEvents
    if (EM(i) < 10.D0) then
      do j = 1, atomCount
        atomPositions(1, j) = SR(1, j, i)
        atomPositions(2, j) = SR(2, j, i)
        atomPositions(3, j) = SR(3, j, i)
      end do
      
      call event%setState(atomIds, atomPositions)
      call event%setEnergy(EM(i))
      call basin%addEvent(event)
    end if
  end do

  ! Before adding the new basin to the network check each if each event
  ! leads to a new unique basin. If not, then flag the event. Reset
  ! after each check.
  
  i = 1

  do while (i <= basin%getAbsorbingCount())
    
    event = basin%getAbsorbingEvent(i)
    call OPT(event, continueSearch)

    do j = 1, atomCount
      atomId = atomIds(j)
      atomPositions(1, j) = RA(1, atomId)
      atomPositions(2, j) = RA(2, atomId)
      atomPositions(3, j) = RA(3, atomId)
    end do

    call tempBasin%setState(atomIds, atomPositions)
    call net%isUniqueBasin(tempBasin, isUnique, matchingBasinId)
    if (.not. isUnique) then
      call basin%flagEvent(i, matchingBasinId)
      i = i - 1
    end if

    if (allocated(atomPositions)) deallocate(atomPositions)
    atomPositions = basin%getState()
    call event%setState(atomIds, atomPositions)
    call OPT(event, continueSearch)

    i = i + 1

  end do

  call net%addBasin(basin)
  if (net%getBasinCount() > 1) then
    activeBasinId = net%getActiveBasinId()
    activeEventId = net%getActiveEventId()
    latestBasinId = net%getBasinCount()
    call net%flagBasinEvent(activeBasinId, activeEventId, latestBasinId)
  end if
  
end subroutine EBN
!***********************************************************************!
subroutine KMC(net, continueSearch)

  use SYSTEM
  use KMCVAR
  use DEFECT
  use MPI
  use LMP_DRIVER
  use Clusters
  use BasinNetwork
  use Basin
  use Event
  use MRM

  implicit none

  type(BasinNetworkType), intent(inout) :: net
  logical, intent(inout) :: continueSearch

  integer :: I, K, jmapsz, type, basinId, eventId
  integer :: numBasins, maxNumEventsPerBasin
  real*8 :: RND, dt, scalef, prop, rate, TOTFQSUM, TAUSUM
  real*8, allocatable, dimension(:,:) :: newRates

  type(BasinType) :: basin
  type(EventType) :: event

  character(LEN=40) str

  ! Use Mean Rate Method to calculate new rates
  call calcMRM(net, temp, newRates, TAUSUM)

  if (rank .eq. 0) then

    FQSUM = 0.0D0 ! frequency sum of current transition states
    TOTFQSUM = 0.0D0 ! approximate frequency sum of entire system
    POSUM = 0.0D0 ! probability sum
    jmapsz = 1000

    ! Calculate total exiting rate
    numBasins = size(newRates,1)
    maxNumEventsPerBasin = size(newRates,2)

    do basinId = 1, numBasins
      do eventId = 1, maxNumEventsPerBasin
        rate = newRates(basinId, eventId)
        if (rate > 0.0D0) then
          FQSUM = FQSUM + rate
        end if
      end do
    end do

    ! Determine mean time in KMC step
    call RANDOM_NUMBER(RND)
    TOTFQSUM = (1/TAUSUM)+OBFQSUM
    if (FQSUM > 0.0D0 .and. TOTFQSUM > 0.0D0) then
      TMEAN = -LOG(RND)/(TOTFQSUM)
      dt = 1.0/FQSUM
      scalef = dt*dfloat(jmapsz)
    end if

    ! Randomly choosing event based on probability
    prop = rnd*float(jmapsz)
    print *, "PROP =", prop
    do basinId = 1, numBasins
      do eventId = 1, maxNumEventsPerBasin
        rate = newRates(basinId, eventId)
        if (rate > 0.0D0) then
          POSUM = POSUM + rate*scalef
        end if
        if (POSUM >= prop) exit
      end do
      if (POSUM >= prop) exit
    end do

  end if

  call MPI_BCAST(basinId, 1, MPI_INT, 0, lmpSystem%getCommunicator(), IERR)
  call MPI_BCAST(eventId, 1, MPI_INT, 0, lmpSystem%getCommunicator(), IERR)

  basin = net % getBasin(basinId)
  event = basin % getAbsorbingEvent(eventId)
  eventEnergy = event % getEnergy()
  type = defects%selectedCluster%type

  EID = eventId

  call net%setActiveIds(basinId, eventId)

  if (eventEnergy < 10.0 .and. eventEnergy >= merge(VACLIM,INTLIM,type<0)) then
      continueSearch = .false.
  else
      continueSearch = .true.
  end if

  ! Only update time and energies when search for KMC step is done.
  if (.not. continueSearch) then
    TIME = TIME + TMEAN
    if (eventEnergy > 0.0D0) then
      I = abs(type)
      K = merge(1, 2, type < 0)
      EMIG(K, I) = (EMIG(K, I)*EMIC(K, I) + eventEnergy)/(EMIC(K, I) + 1)
      EMIC(K, I) = EMIC(K, I) + 1
    end if
  end if

  call MPI_BCAST(continueSearch, 1, MPI_LOGICAL, 0, lmpSystem%getCommunicator(), IERR)
  call MPI_BCAST(TIME, 1, MPI_REAL8, 0, lmpSystem%getCommunicator(), IERR)
  call MPI_BCAST(DTIME, 1, MPI_REAL8, 0, lmpSystem%getCommunicator(), IERR)

end subroutine KMC
!***********************************************************************!
subroutine OPT(event, continueSearch)

  use SYSTEM
  use KMCVAR
  use DEFECT
  use LMP_DRIVER
  use MPI
  use BasinNetwork
  use Basin
  use Event

  implicit none

  type(EventType), intent(in) :: event
  logical, intent(in) :: continueSearch
  
  integer i, atomId, numAtoms
  integer, allocatable, dimension(:) :: ids
  real*8, allocatable, dimension(:,:) :: pos

  ! Update the system with event state

  numAtoms = event % getStateSize()
  ids = event % getIds()
  pos = event % getState()

  do i = 1, numAtoms
    atomId = ids(i)
    RA(1, atomId) = pos(1, i)
    RA(2, atomId) = pos(2, i)
    RA(3, atomId) = pos(3, i)
  end do

  call lmpSystem%assignPositions
  call lmpSystem%relaxSystemEnergy
  call lmpSystem%getPositions

  if (.not. continueSearch) then
    if (RANK == 0 .and. MOD(ISTEP, NFREQ) .EQ. 0) THEN
      write (*, '(3X, 2A12, 7A18)') 'STEP', 'DIMER ID', 'BARRIER', &
        'SYS ENERGY', 'TIME', 'Px', 'Py', 'Pz', 'DPA'
      write (*, '(A,2I12,2F18.5,4ES18.5,E18.5)') 'ITR', ISTEP, EID, &
        eventEnergy, ENCALC, TIME, PRESS(1), PRESS(2), PRESS(3), DPA
      write(*,*)
    end if
  end if

end subroutine OPT
!***********************************************************************!
