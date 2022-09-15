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

    !accessing the atom positions within lammps
    call lammps_extract_atom (x, lmp, 'x')
    if ( .not. associated (x) ) print *, 'x is not associated'
  
  end subroutine LMP_INIT
end module LMP_DRIVER