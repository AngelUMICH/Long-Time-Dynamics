module MRM

  implicit none

  private

  public :: calcMRM

contains

  !***********************************************************************!
  subroutine calcMRM(net, temperature, R, tauSum)
    ! Update the rates for absorbing events using jump-first MRM method.
    ! Given a basin network instance and temperature the resultant array
    ! R contains the new rates for the absorbing events. The mean escape
    ! time for leaving the basin is stored in tauSum.
    

    use BasinNetwork, only: BasinNetworkType
    use Basin, only: BasinType
    use TranEvent, only: TransientEventType
    use Event, only: EventType
    use System, only: KB, NU
    use linear_solve, only: inverse

    class(BasinNetworkType), intent(in) :: net
    real*8, intent(in) :: temperature
    real*8, allocatable, dimension(:,:), intent(out) :: R
    real*8, intent(out) :: tauSum

    integer :: eventId, i, j
    integer :: nBasinCount, nAbsorbing, maxNAbsorbing, nTransient

    real*8 :: rateSum, rate, energy
    real*8, allocatable, dimension(:) :: tau, omega, omegaTot
    real*8, allocatable, dimension(:, :) :: T, oneMinusT, oneMinusTInv

    type(BasinType) :: tempb
    type(TransientEventType) :: tempte
    type(EventType) :: tempe

    nBasinCount = net % getBasinCount()

    allocate (T(nBasinCount, nBasinCount))
    allocate (oneMinusT(nBasinCount, nBasinCount))
    allocate (oneMinusTInv(nBasinCount, nBasinCount))
    allocate (tau(nBasinCount))
    allocate (omega(nBasinCount))
    allocate (omegaTot(nBasinCount))

    T = 0.0d0
    tau = 0.0d0
    omega = 0.0d0
    omegaTot = 0.0d0
    oneMinusT = 0.0d0
    oneMinusTInv = 0.0d0
    maxNAbsorbing = 0
    rateSum = 0.0d0
    tauSum = 0.0d0
    i = 0
    j = 0
    
    ! Calculate residence time for each transient event
    do i = 1, nBasinCount
      tempb = net % getBasin(i)
      nTransient = tempb % getTransientCount()
      do eventId = 1, nTransient
        tempte = tempb % getTransientEvent(eventId)
        energy = tempte % getEnergy()
        rate = nu*exp(-energy/KB/temperature)
        rateSum = rateSum + rate
      end do

      nAbsorbing = tempb % getAbsorbingCount()
      do eventId = 1, nAbsorbing
        tempe = tempb % getAbsorbingEvent(eventId)
        energy = tempe % getEnergy()
        rate = nu*exp(-energy/KB/temperature)
        rateSum = rateSum + rate
      end do
      if (nAbsorbing > maxNAbsorbing) maxNAbsorbing = nAbsorbing

      if (rateSum == 0.0d0) then
        write (*, '(A)') 'MRM Error: Rate sum equals zero.'
        exit
      end if
      tau(i) = 1.0d0 / rateSum
    end do

    ! Calculate rate leaving basin i and entering basin j (transients only)
    do i = 1, nBasinCount
      tempb = net % getBasin(i)
      nTransient = tempb % getTransientCount()
      do eventId = 1, nTransient
        tempte = tempb % getTransientEvent(eventId)
        j = tempte % getBasinId()
        energy = tempte % getEnergy()
        rate = nu*exp(-energy/KB/temperature)
        T(j, i) = tau(i) * rate
      end do
    end do

    ! Calculate omega total for each transient event
    do i = 1, nBasinCount
      do j = 1, nBasinCount
        if (i == j) then
          oneMinusT(i, j) = 1.0d0 - T(i, j)
        else
          oneMinusT(i, j) = 0.0d0 - T(i, j)
        end if
      end do
    end do

    omega(1) = 1.0d0
    oneMinusTInv = inverse(oneMinusT)
    omegaTot = MATMUL(oneMinusTInv, omega)

    do i = 1, nBasinCount
      tau(i) = tau(i) * omegaTot(i)
      tauSum = tauSum + tau(i)
    end do

    ! Output new rates for absorbing events
    allocate(R(nBasinCount, maxNAbsorbing))
    R = 0.0d0
    do i = 1, nBasinCount
      tempb = net % getBasin(i)
      nAbsorbing = tempb % getAbsorbingCount()
      do j = 1, nAbsorbing
        tempe = tempb % getAbsorbingEvent(j)
        energy = tempe % getEnergy()
        rate = nu*exp(-energy/KB/temperature)
        R(i, j) = tau(i)/tauSum * rate
      end do

    end do
    
  end subroutine calcMRM
  !***********************************************************************!
end module MRM