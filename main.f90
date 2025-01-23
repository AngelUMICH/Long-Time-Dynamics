program main

  use MPI
  use SYSTEM, ONLY: IERR, RANK, NPROCS

  implicit none

  call MPI_INIT(IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)

  if (rank == 0) call VERSION

  call readInput
  call lmpInit
  call sysInput
  call sysInit
  call opld

  call MPI_FINALIZE(IERR)

end program main
