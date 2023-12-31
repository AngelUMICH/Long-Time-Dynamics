!***********************************************************************!
module LMP_DRIVER

  use LAMMPS

  implicit none

!lmp: interactable LAMMPS object
!lmpDimer: LAMMPS object for dimer searches
  type(C_ptr) :: lmp
  type(C_ptr) :: lmpDimer
!lmpTemp: current measurement of system temperature
  real(C_double), pointer :: lmpTemp => NULL()
!lmpEnergy: current measurement of system energy
  real(C_double), pointer :: lmpEnergy => NULL()
!lmpPressure: current measurement of system pressure tensor
  real(C_double), dimension(:), pointer :: lmpPressure => NULL()
!lmpL*: length of the system in the respective dimensions
  real*8, dimension(3) :: lmpBoxLo = 0, lmpBoxHi = 0
  real*8 :: xy = 0, yz = 0, xz = 0
  logical, dimension(3) :: periodicity = 0
  logical :: box_change = 0

!***********************************************************************!
contains
!***********************************************************************!
  subroutine LMP_INIT

    use MPI
    use SYSTEM

    implicit none
    character(len=200) :: str

    !attaches lammps object to pointer variable lmp
    call lammps_open('lmp -log none -screen none', MPI_COMM_WORLD, lmp)

    !Setting up system in LAMMPS space

    !Necessary commands for OPLD to work
    call lammps_command(lmp, 'units metal')
    call lammps_command(lmp, 'boundary p p p')
    call lammps_command(lmp, 'atom_style atomic')
    call lammps_command(lmp, 'atom_modify map array')

    if (TLAT .eq. 1) then
      write (str, '(a,x,f7.5)') 'bcc', rlat
    else if (TLAT .eq. 2) then
      write (str, '(a,x,f7.5)') 'fcc', rlat
    end if

    if (rank .eq. 0) print *, 'lattice '//trim(str)
    call lammps_command(lmp, 'lattice '//trim(str))

    write (str, '(f3.1,x,f8.5,x,f3.1,x,f8.5,x,f3.1,x,f8.5)') 0.0, rbox, 0.0, rbox, 0.0, rbox
    if (rank .eq. 0) print *, 'region box block '//trim(str)//' units lattice'
    call lammps_command(lmp, 'region box block '//trim(str)//' units lattice')
    call lammps_command(lmp, 'create_box 1 box')

    !Read from input file using LAMMPS
    if (rank .eq. 0) print *, "Reading in initial positions"
    call lammps_file(lmp, 'in.lammps')
    call lammps_command(lmp, 'thermo_style custom step pe lx ly lz xlo ylo zlo pxx pyy pzz')

    !Including necessary computes for OPLD
    call lammps_command(lmp, 'compute TotEn all pe')
    call lammps_command(lmp, 'compute TotPress all pressure NULL virial')

    !Initiating system values
    call lammps_command(lmp, 'run 0')

    !Extracting MD values for cascade
    call lammps_extract_compute(lmpTemp, lmp, 'thermo_temp', LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)
    call lammps_extract_compute(lmpEnergy, lmp, 'thermo_pe', LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)
    call lammps_extract_compute(lmpPressure, lmp, 'thermo_press', LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR)

    !!SETTING UP COMMS FOR DIMER LAMMPS OBJECTS
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, RANK, RANK, DIMER_COMM, IERR)

    CALL LMP_DIMER_INIT

  end subroutine LMP_INIT
!***********************************************************************!
  subroutine LMP_DIMER_INIT

    use MPI
    use SYSTEM

    implicit none

    character(len=200) :: str

    !attaches lammps object to pointer variable lmp
    call lammps_open('lmp -log none -screen none', DIMER_COMM, lmpDimer)

    !Necessary commands for OPLD to work
    call lammps_command(lmpDimer, 'units metal')
    call lammps_command(lmpDimer, 'boundary p p p')
    call lammps_command(lmpDimer, 'atom_style atomic')
    call lammps_command(lmpDimer, 'atom_modify map array')

    if (TLAT .eq. 1) then
      write (str, '(a,x,f7.5)') 'bcc', rlat
    else if (TLAT .eq. 2) then
      write (str, '(a,x,f7.5)') 'fcc', rlat
    end if

    call lammps_command(lmpDimer, 'lattice '//trim(str))

    write (str, '(f3.1,x,f8.5,x,f3.1,x,f8.5,x,f3.1,x,f8.5)') 0.0, rbox, 0.0, rbox, 0.0, rbox
    call lammps_command(lmpDimer, 'region box block '//trim(str)//' units lattice')
    call lammps_command(lmpDimer, 'create_box 1 box')

    !Read from input file using LAMMPS
    call lammps_file(lmpDimer, 'in.lammps')

    !Dimer lammps object initially has no atoms in it
    call lammps_command(lmpDimer, 'delete_atoms group all')

  end subroutine LMP_DIMER_INIT

!***********************************************************************!
  subroutine LMP_RELAX

    use SYSTEM

    implicit none

    integer :: i, j

    IF (MOD(ISTEP, NFREQ) .EQ. 0 .and. rank .eq. 0) write (*, '(A)') 'Relaxing System'

    !Minimize/relax the system with stress profile considerations
    call lammps_command(lmp, 'min_style cg')
    call lammps_command(lmp, 'minimize 1.0e-10 1.0e-10 5000 50000')

    ENCALC = lmpEnergy
    !reversed sign; positive is tension
    PRESS = -lmpPressure

    IF (MOD(ISTEP, NFREQ) .EQ. 0 .and. rank .eq. 0) write (*, '(A)') 'Finished Relaxation'

    !When ONLY relaxing system
    if (TJOB .EQ. 1) then

      !output relaxed system
      if (RANK .eq. 0) then

        write (*, '(A,F16.5)') 'Total System Energy = ', ENCALC

      !!CONTAINS THE OPT STRUCTURE AFTER RELAXATION
        open (unit=26, file='std_kopt.out')
        write (26, '(I0)') NATOMC
        do i = 1, NATOMC
          write (26, '(I0, 3F14.5)') I, (RA(j, i), j=1, 3)
        end do
      end if

    end if

  end subroutine LMP_RELAX
!***********************************************************************!
  subroutine LMP_STRESS_RELAX

    use system

    implicit none

    integer :: i
    character(len=200) :: str

    if (ISTEP .eq. 1 .and. SigRelOn) then

      if (rank .eq. 0) write (*, '(A)') 'Relaxing System'

      write (str, '(A,F10.3,A,F10.3,A,F10.3)') 'x ', SIGX, ' y ', SIGY, ' z ', SIGZ
      if (rank .eq. 0) write (*, '(A)') 'fix f_boxrelax all box/relax '//trim(str)//' fixedpoint 0 0 0'

      call lammps_command(lmp, 'fix f_boxrelax all box/relax '//trim(str)//' fixedpoint 0 0 0')
      call lammps_command(lmp, 'min_style cg')
      call lammps_command(lmp, 'minimize 1.0e-8 1.0e-8 5000 50000')
      call lammps_command(lmp, 'unfix f_boxrelax')

      ENCALC = lmpEnergy
      !reversed sign; positive is tension
      PRESS = -lmpPressure

      if (rank .eq. 0) write (*, '(A)') 'Finished Relaxation'

      ! Update the reference lattice
      do I = 1, NATOMR
        RR(1, I) = RR(1, I)/EPSX
        RR(2, I) = RR(2, I)/EPSY
        RR(3, I) = RR(3, I)/EPSZ
      end do

      call lammps_extract_box(lmp, lmpBoxLo, lmpBoxHi, xy, yz, xz, periodicity, box_change)

      EPSX = (lmpBoxHi(1) - lmpBoxLo(1))/(rbox*rlat)
      EPSY = (lmpBoxHi(2) - lmpBoxLo(2))/(rbox*rlat)
      EPSZ = (lmpBoxHi(3) - lmpBoxLo(3))/(rbox*rlat)

      do I = 1, NATOMR
        RR(1, I) = RR(1, I)*EPSX
        RR(2, I) = RR(2, I)*EPSY
        RR(3, I) = RR(3, I)*EPSZ
      end do

      write (str, '(A,F10.6,A,F10.6,A,F10.6)') 'x final 0 ', RBOX*EPSX, ' y final 0 ', RBOX*EPSY, ' z final 0 ', RBOX*EPSZ
      call lammps_command(lmpDimer, 'change_box all '//trim(str)//' units lattice')

    end if

  end subroutine LMP_STRESS_RELAX
!***********************************************************************!

  subroutine LMP_CASCADE

    use MPI
    use SYSTEM

    implicit none

    !m: mass of system atoms in grams/mole
    !b: border thickness in lattice units
    !ETOL: energy tolerance (eV) for cascade
    !TTOL: temperature tolerance (K) for equilibration
    real*8 :: m, b, ETOL, TTOL
    !v: magnitude of velocity in angs/picosec
    !vx: x component of velocity
    !rnd: random number
    real*8 :: v, vx, vy, vz, rnd
    !ev_J: conversion from 1 eV to Joules
    !amu_kg: conversion of 1 amu to kilograms
    real*8 :: eV_J, amu_kg
    !d_TEMP: temperature difference
    !d_Energy: energy difference
    !lastEnergy: previous energy measurement
    real*8 :: d_TEMP, d_Energy, lastEnergy
    !xDist,yDist,zDist
    real*8 :: xDist, yDist, zDist
    !str: general string
    !s_rnd: random number str
    !s_TEMP: TEMP str
    !s_border: selected border
    character(len=200) :: str, s_rnd, s_TEMP, s_PKA, s_border

    m = 56
    b = 2.5
    ETOL = 30
    TTOL = 5
    eV_J = 1.6E-19 !energy conversion factor
    amu_kg = 1.66E-27 !mass conversion factor

    !Setting up regions and groups

    !PKA selection
    call random_number(rnd)
    rnd = ceiling(rnd*NATOMC)
    write (s_PKA, '(3F10.3)') RA(:, int(rnd))
    call MPI_BCAST(s_PKA, 200, MPI_CHAR, 0, MPI_COMM_WORLD, IERR)

    call lammps_command(lmp, 'region PKACHOOSING sphere '//trim(s_PKA)//' 0.2 units lattice')
    call lammps_command(lmp, 'group  PKA region PKACHOOSING')

    !Border group selection
    xDist = RA(1, int(rnd)) - RBH1*EPSX
    yDist = RA(2, int(rnd)) - RBH1*EPSY
    zDist = RA(3, int(rnd)) - RBH1*EPSZ

    if (abs(xDist) .le. min(abs(yDist), abs(zDist))) then
      if (xDist .ge. 0) then
        !'xlo'
        write (str, '(A,F10.3,A)') '0 ', b, ' INF INF INF INF'
        s_border = 'xlo'
      else
        !'xhi'
        write (str, '(F10.3,A)') RBOX*EPSX - b, ' INF INF INF INF INF'
        s_border = 'xhi'
      end if
    else if (abs(yDist) .le. min(abs(xDist), abs(zDist))) then
      if (yDist .ge. 0) then
        !'ylo'
        write (str, '(A,F10.3,A)') 'INF INF 0 ', b, ' INF INF'
        s_border = 'ylo'
      else
        !'yhi'
        write (str, '(A,F10.3,A)') 'INF INF ', RBOX*EPSY - b, ' INF INF INF'
        s_border = 'yhi'
      end if
    else if (abs(zDist) .le. min(abs(xDist), abs(yDist))) then
      if (zDist .ge. 0) then
        !'zlo'
        write (str, '(A,F10.3)') 'INF INF INF INF 0 ', b
        s_border = 'zlo'
      else
        !'zhi'
        write (str, '(A,F10.3,A)') 'INF INF INF INF ', RBOX*EPSZ - b, ' INF'
        s_border = 'zhi'
      end if
    end if

    call MPI_BCAST(str, 200, MPI_CHAR, 0, MPI_COMM_WORLD, IERR)

    call lammps_command(lmp, 'region border block '//trim(str)//' units lattice')
    call lammps_command(lmp, 'group border region border')
    call lammps_command(lmp, 'group cascadeRegion subtract all border')

    !Starting dynamics

    !Bring system to temperature
    if (rank .eq. 0) write (*, '(A)') 'Starting cascade simulation'
    if (rank .eq. 0) write (*, '(A)') 'First bringing system up to temperature'
    call random_number(rnd)
    rnd = floor(rnd*100000000)
    write (s_rnd, *) int(rnd)
    write (s_TEMP, '(F9.3)') TEMP

    call lammps_command(lmp, 'velocity all create '//trim(s_TEMP)//' ' &
                        //trim(s_rnd)//' dist gaussian mom yes rot yes')
    call lammps_command(lmp, 'fix fixnpt all nvt temp '//trim(s_TEMP)//' ' &
                        //trim(s_TEMP)//' 0.1')

    d_TEMP = 100

    if (rank .eq. 0) write (*, '(A)') 'Running temp equilibration for 2500 step increments'

    do while (abs(d_TEMP) .gt. TTOL)
      call lammps_command(lmp, 'run 2500')
      d_TEMP = TEMP - lmpTemp
      if (rank .eq. 0) write (*, '(A,F10.3)') 'Temperature diff: ', abs(d_TEMP)
    end do

    call lammps_command(lmp, 'unfix fixnpt')

    !Determining velocity vector for PKA
    v = sqrt(2*PKAE/m)
    !Convert to angs/picosec
    v = v*sqrt(eV_J/amu_kg)*1E-2

    !Determine random direction
    call random_number(vx)
    call random_number(vy)
    vx = (vx*1*2) - 1
    vy = (vy*sqrt(1 - vx*vx)*2) - sqrt(1 - vx*vx)
    vz = sqrt(1 - vx*vx - vy*vy)
    vx = vx*v
    vy = vy*v
    vz = vz*v
    write (str, '(F9.3, X, F9.3, X, F9.3)') vx, vy, vz
    call MPI_BCAST(str, 200, MPI_CHAR, 0, MPI_COMM_WORLD, IERR)

    !Set PKA to velocity
    call lammps_command(lmp, 'velocity PKA set '//trim(str)//' units box')

    if (rank .eq. 0) then
      write (*, '(A)') '!**************************Cascade Info**************************!'
      write (*, '(A18, A11, A18, A19)') 'PKA Energy (eV)', 'Velocity', 'System Temp (K)', 'Border Selection'
      write (*, '(F18.5, F11.3, F18.3, A19)') PKAE, v, lmpTemp, trim(s_border)
      write (*, *)
      write (*, '(6A10)') 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz'
      write (*, '(A,3F10.3)') trim(s_PKA), vx, vy, vz
      write (*, *)
      write (*, '(A)') 'Positions Units: (lattice units)    Velocity Units: (ang/ps)'
      write (*, *)
      if (abs(lmpTemp - TEMP) .gt. 5) write (*, '(A,F8.3,A)') 'Warning: System temp has not reached ', TEMP, ' K'
      write (*, '(A)') '!****************************************************************!'
    end if

    call lammps_command(lmp, 'fix cascade cascadeRegion nve')
    call lammps_command(lmp, 'fix EnergyControl border nvt temp '//trim(s_TEMP) &
                        //' '//trim(s_TEMP)//' 0.1')
    call lammps_command(lmp, 'fix adaptivedt all dt/reset 10 1.0E-10 2.0E-3 0.02 units box')

    !Needed to begin dynamics
    d_Energy = 1.0E12
    lastEnergy = 0

    if (rank .eq. 0) write (*, '(A)') 'Running dynamics for 2500 step increments'

    do while (abs(d_Energy) .gt. ETOL)
      call lammps_command(lmp, 'run 2500')
      d_Energy = lmpEnergy - lastEnergy
      if (rank .eq. 0) write (*, '(A,F14.3)') 'Energy diff: ', abs(d_Energy)
      lastEnergy = lmpEnergy
    end do

    call lammps_command(lmp, 'unfix cascade')
    call lammps_command(lmp, 'unfix EnergyControl')
    call lammps_command(lmp, 'unfix adaptivedt')

    if (rank .eq. 0) write (*, '(A)') 'Finished cascade simulation'

    !Ready system for minimization
    if (rank .eq. 0) write (*, '(A)') 'Zeroing atom velocities'
    call lammps_command(lmp, 'velocity all zero linear')

  end subroutine LMP_CASCADE
!***********************************************************************!
  subroutine RECV_LMP_POS

    use SYSTEM, only: NATOMC, RA, RLAT

    implicit none

    integer :: i
    real*8, dimension(:), allocatable :: pos

    call lammps_gather_atoms(lmp, 'x', 3, pos)

    DO I = 1, NATOMC
      RA(:, i) = pos(3*(i - 1) + 1:3*(i - 1) + 3)/RLAT
    end do

  end subroutine RECV_LMP_POS
!***********************************************************************!
  subroutine SEND_OPLD_POS

    use SYSTEM, only: NATOMC, RA, RLAT, IERR
    use MPI

    implicit none

    integer :: i
    real*8, dimension(:) :: pos(3*NATOMC)

    CALL MPI_BCAST(RA, size(RA), MPI_REAL8, 0, MPI_COMM_WORLD, IERR)

    DO I = 1, NATOMC
      pos(3*(i - 1) + 1:3*(i - 1) + 3) = RA(:, i)*RLAT
    END DO

    CALL lammps_scatter_atoms(lmp, 'x', pos)

  end subroutine SEND_OPLD_POS
!***********************************************************************!
end module LMP_DRIVER
