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

      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)

      CALL READ_INPUT
      CALL LMP_INIT
      CALL SYS_INPUT
      CALL SYS_INIT
      CALL LMP_DIMER_INIT

      IF(TJOB.EQ.1) THEN

            CALL LMP_RELAX

      ELSE IF(TJOB.EQ.2) THEN

        IF (MOD(NDIMER, NPROCS) .NE. 0) THEN

          WRITE(*,*) "# MOD(NDIMER, NPROCS) .NE. 0"
          WRITE(*,*) "PLEASE CHANGE NDIMER OR NUM OF PROCESSES"
          CALL MPI_FINALIZE(IERR)
          WRITE(*,*) "STOPPING PROGRAM"
          STOP

        ELSE

          CALL OPLD

        ENDIF

      ENDIF

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
