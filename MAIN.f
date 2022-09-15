!***********************************************************************!
!***********************************************************************!
!***          SELF EVOLVING ATOMISTIC KINETIC MONTE CARLO CODE       ***!
!***                        AUTHOR:XU,HAIXUAN                        ***!
!***             MATERIALS SCIENCE AND TECHNOLOGY DIVISION           ***!
!***                      OAK RIDGE NATIONAL LAB                     ***!
!***                          XUH1@ORNL.GOV                          ***!
!***                           VERSION 3.0                           ***!
!***                       ALL RIGHTS RESERVED                       ***!
!***********************************************************************!
!***********************************************************************!
      PROGRAM MAIN

      USE MPI
      USE LMP_DRIVER
      USE SYSTEM,ONLY : TJOB, IERR, RANK, NPROCS, NDIMER

      IMPLICIT NONE

      !INTEGER IERR, RANK, NPROCS !OMP VARIABLES

!! OMP INITIALIZATION
!! AT THE MOMENT ONLY EVEN NUMBER OF PROCESSES ARE CONSIDERED
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)


!!READING INPUT FROM VARIOUS FILES

      CALL SYS_INPUT

!!SYSTEM INITIALIZATION

      CALL SYS_INIT

!!LAMMPS INITIALIZATION

      CALL LMP_INIT

!!IF A STATIC RELAXATION IS CHOSEN (TAKEN OUT: ONLY WORKS FOR SERIAL)

      IF(TJOB.EQ.1) THEN

        IF (RANK .EQ. 0) THEN

          CALL LMP_RELAX

        ENDIF
!!IF SELF-EVOVLING ATOMISTIC KMC IS CHOSEN

      ELSE IF(TJOB.EQ.2) THEN

        IF (MOD(NDIMER, NPROCS) .NE. 0) THEN

          WRITE(*,*) "# MOD(NDIMER, NPROCS) .NE. 0"

          WRITE(*,*) "PLEASE CHANGE NDIMER OR NUM OF PROCESSES"

          CALL MPI_FINALIZE(IERR)

          WRITE(*,*) "STOPPING PROGRAM"

          STOP

        ELSE

          CALL SEAKMC

        ENDIF

      ENDIF

!!SYSTEM CLEANUP,DEALLOCATE MEMORIES,CLOSE FILES

      !CALL FINISH

      CALL MPI_FINALIZE(IERR)

      END PROGRAM MAIN

!***********************************************************************!
!***********************************************************************!
!***          SELF EVOLVING ATOMISTIC KINETIC MONTE CARLO CODE       ***!
!***                        AUTHOR:XU,HAIXUAN                        ***!
!***             MATERIALS SCIENCE AND TECHNOLOGY DIVISION           ***!
!***                      OAK RIDGE NATIONAL LAB                     ***!
!***                          XUH1@ORNL.GOV                          ***!
!***                           VERSION 3.0                           ***!
!***                       ALL RIGHTS RESERVED                       ***!
!***********************************************************************!
!***********************************************************************!
