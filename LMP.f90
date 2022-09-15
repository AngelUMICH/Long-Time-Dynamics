!***********************************************************************!
module LMP_DRIVER

  use LAMMPS

  implicit none

  !lmp: interactable LAMMPS object
  type (C_ptr) :: lmp
  !x:   pointer to atom positions within lammps mem
  real (C_double), dimension(:,:), pointer :: x => NULL()
  !pe:  potential energy
  real (C_double), pointer :: pe => NULL()

!***********************************************************************!
contains
!***********************************************************************!
  subroutine LMP_INIT
    
    use MPI
    use SYSTEM

    implicit none

    integer :: i
    real (C_double) :: r(3)

    !attaches lammps object to pointer variable lmp
    call lammps_open('', MPI_COMM_WORLD, lmp)

    !Setting up system in LAMMPS space
    call lammps_command(lmp, 'units metal')
    call lammps_command(lmp, 'atom_style atomic')
    call lammps_command(lmp, 'lattice bcc 2.8553')

    call lammps_command(lmp, 'region box block 0 30 0 30 0 30')
    call lammps_command(lmp, 'create_box 1 box')
    call lammps_command(lmp, 'create_atoms 1 box')
    call lammps_command(lmp, 'mass 1 56.0')

    call lammps_command(lmp, 'pair_style eam/fs')
    call lammps_command(lmp, 'pair_coeff * * FeP_mm.eam.fs Fe')
    call lammps_command(lmp, 'neighbor 0 bin')

    call lammps_command(lmp, 'velocity all zero linear')

    !Including necessary computes
    call lammps_command(lmp, 'compute TotEn all pe')

    !accessing the atom positions within lammps
    call lammps_extract_atom (x, lmp, 'x')
    if ( .not. associated (x) ) print *, 'x is not associated'
  
  end subroutine LMP_INIT
!***********************************************************************!
  subroutine LMP_RELAX

  use SYSTEM

  implicit none

  integer :: i,j

  !Minimize/relax the system
  call lammps_command(lmp, 'min_style cg')
  call lammps_command(lmp, 'minimize 1.0e-4 1.0e-6 100 1000')

  !When ONLY relaxing system
  if (TJOB.EQ.1) then

    call LMP_ENERGY
    
    write(*,*) 'Total System Energy = ', ENCALC

    !output relaxed system to std_kopt.out
    if (rank == 0) then
      write(26,*) NATOMC
        do i = 1, NATOMC
          write(26,'(3F14.5)') (RA(j,i), j=1,3)
        end do
    end if

  end if

  end subroutine LMP_RELAX
!***********************************************************************!
  subroutine LMP_ENERGY

  use SYSTEM, only: ENCALC

  implicit none

  !OENCALC: Old system energy
  real*8 :: OENCALC

  OENCALC = ENCALC

  !Using 'run 0' to update computes
  call lammps_command(lmp, 'run 0')

  !Might not need to extract pe everytime. It is already bound
  call lammps_extract_compute (pe, lmp, 'TotEn', &
                                LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)
  if ( .not. associated (pe) ) print *, 'pe is not associated'
  ENCALC = pe

  if (OENCALC == ENCALC) write(*,*) 'Warning: sys energy might not be updated'
  
  end subroutine LMP_ENERGY
!***********************************************************************!
end module LMP_DRIVER