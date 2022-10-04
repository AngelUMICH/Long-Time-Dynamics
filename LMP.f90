!***********************************************************************!
module LMP_DRIVER

  use LAMMPS

  implicit none

  !lmp: interactable LAMMPS object
  type (C_ptr) :: lmp
  !pe:  potential energy
  real (C_double), pointer :: pe => NULL()
  !compute: a lammps compute value
  real (C_double), pointer :: compute => NULL()

!***********************************************************************!
contains
!***********************************************************************!
  subroutine LMP_INIT
    
    use MPI
    use SYSTEM

    implicit none
    integer :: i, j
    real*8, DIMENSION(:), allocatable :: pos


    !attaches lammps object to pointer variable lmp
    call lammps_open('lmp -log none -screen none', MPI_COMM_WORLD, lmp)

    !Setting up system in LAMMPS space
    call lammps_command(lmp, 'units metal')
    call lammps_command(lmp, 'boundary p p p')
    call lammps_command(lmp, 'atom_style atomic')
    call lammps_command(lmp, 'atom_modify map array')
    call lammps_command(lmp, 'lattice bcc 2.8553')

    call lammps_command(lmp, 'region box block 0 30 0 30 0 30')
    call lammps_command(lmp, 'create_box 1 box')
    call lammps_command(lmp, 'create_atoms 1 box')
    call lammps_command(lmp, 'mass 1 56.0')

    call lammps_command(lmp, 'pair_style eam/fs')
    call lammps_command(lmp, 'pair_coeff * * FeP_mm.eam.fs Fe')
    call lammps_command(lmp, 'neighbor 0 bin')

    call lammps_command(lmp, 'velocity all zero linear')
    call lammps_command(lmp, 'timestep 0.001')
    call lammps_command(lmp, 'thermo 200')

    !Including necessary computes
    call lammps_command(lmp, 'compute TotEn all pe')

    CALL lammps_gather_atoms(lmp, 'x', 3, pos)
    
    do i =1, NATOMC
        RA(:,i) = pos(3*(i-1)+1: 3*(i-1)+3) / RLAT
    end do
  
  end subroutine LMP_INIT
!***********************************************************************!
  subroutine LMP_RELAX

  use SYSTEM

  implicit none

  integer :: i,j

  if (rank == 0) write(*,'(A)') 'Starting Relaxation'

  !Minimize/relax the system
  call lammps_command(lmp, 'min_style cg')
  call lammps_command(lmp, 'minimize 1.0e-4 1.0e-6 100 1000')

  call LMP_ENERGY
  if (rank == 0) write(*,'(A)') 'Ending Relaxation'

  !When ONLY relaxing system
  if (TJOB.EQ.1) then

    !output relaxed system to std_kopt.out
    if (RANK == 0) then

      write(*,'(A,F16.5)') 'Total System Energy = ', ENCALC

      write(26,'(I0)') NATOMC
      do i = 1, NATOMC
        write(26,'(I0, 3F14.5)') I, (RA(j,i), j=1,3)
      end do
    end if

  end if

  end subroutine LMP_RELAX
!***********************************************************************!
  subroutine LMP_ENERGY

  use SYSTEM, only: ENCALC, RANK

  implicit none

  !OENCALC: Old system energy
  real*8 :: OENCALC

  OENCALC = ENCALC

  !Might not need to extract pe everytime. It is already bound
  call lammps_extract_compute (pe, lmp, 'TotEn', &
                                LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)
  if ( .not. associated (pe) .and. RANK == 0) write(*,*) 'pe is not associated'
  ENCALC = pe

  if (OENCALC == ENCALC .and. RANK == 0) write(*,*) 'Warning: sys energy might not be updated'
  
  end subroutine LMP_ENERGY
!***********************************************************************!
  subroutine LMP_CASCADE

  use MPI
  use SYSTEM

  implicit none

  !m: mass of system atoms in grams/mole
  !E: energy of pka in eV
  real*8 :: m, E, i
  !v: magnitude of velocity in angs/picosec
  !vx: x component of velocity
  !rnd: random number
  real*8 :: v, vx, vy, vz, rnd
  !ev_J: conversion from 1 eV to Joules
  !amu_kg: conversion of 1 amu to kilograms
  real*8 :: eV_J, amu_kg
  !str: general string
  !s_rnd: random number str
  !S_TEMP: TEMP str
  character (len = 200):: str, s_rnd, s_TEMP
  !pos: communicating atom positions with lammps
  real*8, dimension(:) :: pos(3*NATOMC)

  m = 56
  E = 1000
  eV_J = 1.6E-19
  amu_kg = 1.66E-27

  !Bring system to temperature
  if (rank == 0) write(*,'(A)') 'Starting cascade simulation'
  if (rank == 0) write(*,'(A)') 'First bringing system up to temperature'
  call random_number(rnd)
  rnd = floor(rnd * 100000000)
  write(s_rnd,*) int(rnd)
  write(s_TEMP,'(F9.3)') TEMP

  call lammps_command(lmp, 'velocity all create '//trim(s_TEMP)//' '&
                            //trim(s_rnd)//' dist gaussian mom yes rot yes')
  call lammps_command(lmp, 'fix fixnpt all npt temp '//trim(s_TEMP)//' '&
                            //trim(s_TEMP)//' 0.1 iso 0 0 0.1')
  call lammps_command(lmp, 'run 2500')
  call lammps_command(lmp, 'unfix fixnpt')
  call lammps_extract_compute (compute, lmp, 'thermo_temp', LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)

  !Determining velocity vector for PKA
  v = sqrt(2*E/m)
  !Convert to angs/picosec
  v = v * sqrt(eV_J/amu_kg) * 1E-2

  !Determine random direction
  call random_number(vx)
  call random_number(vy)
  vx = (vx*1*2) - 1
  vy = (vy*sqrt(1 - vx*vx)*2) - sqrt(1 - vx*vx)
  vz = sqrt(1 - vx*vx - vy*vy)
  vx = vx * v
  vy = vy * v
  vz = vz * v
  write(str,'(F9.3, X, F9.3, X, F9.3)') vx, vy, vz
  call MPI_BCAST(str, 200, MPI_CHAR, 0, MPI_COMM_WORLD, IERR)

  !Initializing PKA setup
  call lammps_command(lmp, 'region PKACHOOSING sphere 15 15 15 0.2 units lattice')
  call lammps_command(lmp, 'group  PKA region PKACHOOSING')
  call lammps_command(lmp, 'velocity PKA set '//trim(str)//' units box')

  if (rank == 0) then
    write(*,'(A)') '!*********************************Cascade Info********************************************!'
    write(*,'(A18, A25, 4A12)') 'PKA Energy (eV)', 'Velocity (ang/picosec)', 'Vx', 'Vy', 'Vz', 'Temp (K)'
    write(*,'(F18.5, F25.5, 4F12.3)') E, v, vx, vy, vz, compute
    write(*,*)
    if (abs(compute - TEMP) > 5) write(*,'(A,F8.3,A)') 'Warning: System temp has not reached ', TEMP, ' K'
    write(*,'(A)') '!*****************************************************************************************!'
  end if

  !Running cascade simulation
  call lammps_command(lmp, 'fix fixcascade all nve')
  call lammps_command(lmp, 'run 5000')
  call lammps_command(lmp, 'unfix fixcascade')
  if (rank == 0) write(*,'(A)') 'Finished cascade simulation'

  !Relax system
  if (rank == 0) write(*,'(A)') 'Zeroing atom velocities'
  call lammps_command(lmp, 'velocity all zero linear')

  end subroutine LMP_CASCADE
!***********************************************************************!
end module LMP_DRIVER