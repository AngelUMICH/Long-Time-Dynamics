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

!***********************************************************************!
!!MAJOR SUBROUTINE OF SELF-EVOLVING ATOMISTIC KINETIC MONTE CARLO CODE**!

      SUBROUTINE SEAKMC

      USE SYSTEM,ONLY : NTSTEP,ISTEP, RANK,NATOMC,RA, SKIP

      IMPLICIT NONE

      INTEGER I,J

      INCLUDE 'mpif.h'

      DO ISTEP = 1,NTSTEP

        IF (RANK .EQ. 0) THEN
  !!CHARACTERIZE ACTIVE VOLUMES

          CALL DCIDNT

          CALL DFSEL

    !!CHECK BOUNDARY CONDITIONS

          CALL CBC

    !!RECYCLING SADDLE POINTS

          CALL RSP

        END IF

  !!SADDLE POINT SEARCHING - ALL PROCESSES

        CALL SPS

        IF (RANK .EQ. 0) THEN

  !!KINETIC MONTE CARLO

          CALL KMC

  !!INTRODUCE NEW FRENKEL DEFECTS DUE TO DOSE RATE  
          IF (SKIP .EQV. .FALSE.) THEN

            CALL DOSE
          
  !!STATIC RELAXATION

            CALL OPT         
          
          ENDIF

        END IF

      ENDDO

      IF (RANK .EQ. 0) THEN

        CALL VISUAL

        CALL EXECUTE_COMMAND_LINE('mv mig*.log std_evo.out ./DUMP')

       WRITE(26,*) NATOMC

       DO I = 1, NATOMC
        WRITE(26,'(3F14.5)') (RA(J,I), J=1,3)
       ENDDO

      ENDIF

      RETURN

      END SUBROUTINE
!***********************************************************************!
!!THIS SUBROUTINE IS TO CHARACTERIZE THE ACTIVE VOLUMES IN THE SYSTEM

      SUBROUTINE CAV

      USE SYSTEM
      USE DEFECT

      IMPLICIT NONE

      INTEGER I,J,NTEMP

      REAL*8 RX,RY,RZ,RS,PBC,PTEMP(7,NDEMAX)

!!OUTPUT THE DEFECT STRUCTURE TO ATOMTV AND OVITO

      CALL VISUAL

!!OUTPUT INITIAL DEFECTS WITH FILE INDEX 000
      IF (ISTEP .EQ. 1)THEN

        ISTEP = 0

        CALL OVITO

        ISTEP = 1

      END IF

!!INITIALIZATION,ND:NUMBER OF DEFECTS
!!VAC:VACANCY,INT:INTERSTITIAL

      NDVAC = 0
      NDINT = 0

!!INDEX FOR INTERSTITIAL AND VACACNY

      II = 0
      IV = 0

!!IDENTIFY INTERSTITIAL,PI:POSITION OF INTERSITITAL

      DO I=1,NATOMC
        IF(TI(I).NE.1) THEN
          NDINT = NDINT + 1
          II(1,NDINT) = I
          PI(1,NDINT) = RA(1,I)
          PI(2,NDINT) = RA(2,I)
          PI(3,NDINT) = RA(3,I)
        ENDIF
      ENDDO

!!IDENTIFY VACANCY,PV:POSITION OF VACANCY

      DO J=1,NATOMR
        IF(TV(J).NE.1) THEN
          NDVAC = NDVAC + 1
          IV(1,NDVAC) = J
          PV(1,NDVAC) = RR(1,J)
          PV(2,NDVAC) = RR(2,J)
          PV(3,NDVAC) = RR(3,J)
        ENDIF
      ENDDO

      NDTOT = 0

!!COUNT TOTAL NUMBER OF DEFECTS AND DETERMINE THEIR POSITION

      DO I=1,NDINT

        NDTOT = NDTOT + 1

        PX(1,NDTOT) = PI(1,I)
        PX(2,NDTOT) = PI(2,I)
        PX(3,NDTOT) = PI(3,I)

        PX(4,NDTOT) = 4.2D0 !radius of interstitial AV
        PX(5,NDTOT) = 1.0D0 !TELLS IF IN NEED OF SADDLE SEARCHES
        PX(6,NDTOT) = 1.0D0 !THIS MEANS IT IS A INTERSTITIAL

        PX(7,NDTOT) = II(1,I)

      ENDDO

      DO J=1,NDVAC

        NDTOT = NDTOT + 1

        PX(1,NDTOT) = PV(1,J)
        PX(2,NDTOT) = PV(2,J)
        PX(3,NDTOT) = PV(3,J)


        PX(4,NDTOT) = 2.7D0 !radius of vacancy AV
        PX(5,NDTOT) = 1.0D0 !TELLS IF IN NEED OF SADDLE SEARCHES
        PX(6,NDTOT) =-1.0D0 !THIS MEANS IT IS A VACANCY

        PX(7,NDTOT) = IV(1,J)

      ENDDO

!!DETERMINE THE REAL DEFECTS

      DO I=1,NDTOT

        DO J=I+1,NDTOT

          RX=PBC(PX(1,I)-PX(1,J))
          RY=PBC(PX(2,I)-PX(2,J))
          RZ=PBC(PX(3,I)-PX(3,J))

          RS=RX*RX+RY*RY+RZ*RZ

!          IF(RS.LT.(DSQRT(3.0D0)/2.0D0)**2) THEN

          IF(RS.LT.3.000D0) THEN

            IF((PX(5,I).EQ.1.0D0).AND.(PX(5,J).EQ.1.0D0)) THEN

            PX(1,I) = PX(1,I)*ABS(PX(6,I))+PX(1,J)*ABS(PX(6,J))
            PX(2,I) = PX(2,I)*ABS(PX(6,I))+PX(2,J)*ABS(PX(6,J))
            PX(3,I) = PX(3,I)*ABS(PX(6,I))+PX(3,J)*ABS(PX(6,J))

            PX(1,I) = PX(1,I)/(ABS(PX(6,I))+ABS(PX(6,J)))
            PX(2,I) = PX(2,I)/(ABS(PX(6,I))+ABS(PX(6,J)))
            PX(3,I) = PX(3,I)/(ABS(PX(6,I))+ABS(PX(6,J)))

            PX(6,I) = PX(6,I)+PX(6,J)
            PX(5,J) = 0.0D0

            ENDIF

          ENDIF

        ENDDO

      ENDDO

      NTEMP = 0

      DO I=1,NDTOT

        IF((PX(5,I).NE.0.0D0).AND.(PX(6,I).NE.0.0D0)) THEN

        NTEMP = NTEMP + 1

        PTEMP(:,NTEMP) = PX(:,I)

        ENDIF

      ENDDO

      NDTOT = NTEMP

      PX = 0.0D0

!!RESORTING,PX STORES THE POSITIONS OF THE ACTIVE VOLUMES

      DO I=1,NDTOT

        PX(:,I) = PTEMP(:,I)

      WRITE(*,'(A,I12,6(F14.5),I12)') 'CAV',I,PX(1:6,I),
     #INT(PX(7,I))

      ENDDO

      ONDTOT = NDTOT

!!IF NO DEFECT IS PRESENT, STOP THE PROGRM

      IF(NDTOT.EQ.0) THEN

        WRITE(*,*) 'NO DEFECT IN SYSTEM,STOP'

        STOP

      ENDIF

      RETURN
      END SUBROUTINE
!***********************************************************************!
!!THIS SUBROUTINE IS TO CHECK THE BOURNDARY CONDITION
!!ONLY FOR ABSORBING BOUNDARY CONDITION, NOT FOR PBC

      SUBROUTINE CBC

      USE SYSTEM
      USE DEFECT

      INTEGER I

      LOGICAL TOUT,TRELAX

      TRELAX = .FALSE.

!!IF ABSORBING BOUNDARY CONDTION IS USED

      IF(TPBC.EQ.2) THEN

        DO I=1,NDTOT

          TOUT = .FALSE.

!!IF THE DEFECT IS CLOSE TO THE BOUNDARY (LESS THAN 1 LATTICE CONSTANT)
!!IT IS CONSIDERED TO BE REMOVED FROM THE SYSTEM

         IF((PX(1,I).GT.(RBOX-1.0D0)).OR.(PX(1,I).LT.1.0D0)) TOUT=.TRUE.
         IF((PX(2,I).GT.(RBOX-1.0D0)).OR.(PX(2,I).LT.1.0D0)) TOUT=.TRUE.
         IF((PX(3,I).GT.(RBOX-1.0D0)).OR.(PX(3,I).LT.1.0D0)) TOUT=.TRUE.

          IF(TOUT) THEN

            TRELAX = .TRUE.

!!IF THE DEFECT IS A VACANCY, ADD ONE ATOM AT VACANCY POSITION
!!THE SYSTEM ENERGY WILL NOT BE CONSERVED

            IF(PX(6,I).EQ.(-1.0D0)) THEN

              NATOMC = NATOMC + 1

              RA(1,NATOMC) = PX(1,I)
              RA(2,NATOMC) = PX(2,I)
              RA(3,NATOMC) = PX(3,I)

              PX(5,I) = 0.0D0

!!WRITE INFORMATION:WHICH VACANCY IS REMOVED

              WRITE(*,*) 'VAC OBC',PX(1,I),PX(2,I),PX(3,I)

!!IF THE DEFECT IS AN INTERSITITAL, REMOVE THE ATOM FROM THE SYSTEM

            ELSE IF(PX(6,I).EQ.(1.0D0)) THEN

              RA(1,DINT(PX(7,I))) = RA(1,NATOMC)
              RA(2,DINT(PX(7,I))) = RA(2,NATOMC)
              RA(3,DINT(PX(7,I))) = RA(3,NATOMC)

              NATOMC = NATOMC -1

              PX(5,I) = 0.0D0

!!WRITE INFORMATION:WHICH INTERSTITIAL IS REMOVED

              WRITE(*,*) 'INT OBC',PX(1,I),PX(2,I),PX(3,I)

            ENDIF

          ENDIF

        ENDDO

!!RE-OPTIMIZE THE STRUCTURE USING CONJUGATE GRADIENT

        IF(TRELAX) THEN

          CALL OPT_CG_INIT(NATOMC)

          CALL OPT_CG

          CALL OPT_CG_FINISH

          CALL DCIDNT
!          CALL CAV

        ENDIF

      ENDIF

      RETURN

      END SUBROUTINE CBC
!***********************************************************************!
!!THIS SUBROUTINE IS TO RECYCLE THE SADDLE POINTS

      SUBROUTINE RSP

      USE SYSTEM
      USE DEFECT
      USE KMCVAR

      INTEGER I,J

      REAL*8 RX,RY,RZ,RS,TX,TY,TZ,TS,PBC

      SR = 0.0D0 !HAVE TO RESET ARRAY BECAUSE OF MPI_SUM BETWEEN RANK0 AND EVERY OTHER RANK

      EM = 0.0D0

      ND = 0.0D0

!!COMPARE THE CURRENT DEFECTS WITH PREVIOUS DEFECTS

      DO I=1,NDTOT

        DO J=1,ONDTOT

!!IF THE DEFECT INDEX IS THE SAME

          IF(PX(7,I).EQ.OX(7,J)) THEN

            RX = PBC(PX(1,I)-OX(1,J))
            RY = PBC(PX(2,I)-OX(2,J))
            RZ = PBC(PX(3,I)-OX(3,J))

            RS = RX*RX + RY*RY+RZ*RZ

!!IF THE DEFECT POSITION IS THE SAME

            IF(RS.LT.1.0D-8) THEN

              TX = PBC(PX(1,I)-OX(1,JID))
              TY = PBC(PX(2,I)-OX(2,JID))
              TZ = PBC(PX(3,I)-OX(3,JID))

              TS = TX*TX+TY*TY+TZ*TZ

!!IF THE DEFECT IS FAR FROM THE AV CHOSEN BY KMC

              IF(TS.GT.(PX(4,I)+OX(4,JID))**2) THEN

!!THE INFORMATION CAN BE RECYCLYED

              WRITE(*,'(A,I12,4F14.5,I12,4F14.5)') 'RSP',
     #        I,PX(1:4,I),JID,OX(1:4,JID)

              EM(:,I) = OE(:,J)

              ND(:,I) = ON(:,J)

              SR(:,:,:,I) = OR(:,:,:,J)

              PX(5,I) = 0.0D0

              ENDIF

            ENDIF

          ENDIF

        ENDDO

      ENDDO

      RETURN

      END SUBROUTINE RSP
!***********************************************************************!
!!THIS SUBROUTINE IS TO PERFORM SADDLE POINT SEARCHING

      SUBROUTINE SPS

      USE SYSTEM
      USE DEFECT
      USE KMCVAR

      USE CLASS_DIMER

      IMPLICIT NONE

      INCLUDE 'mpif.h'

      TYPE(DIMER) :: ADIMER

      INTEGER I,J,K

      INTEGER RA_LEN, SR_LEN, EM_LEN, ND_LEN !!USED FOR OPENMPI

!!MPI VARIABLES
      RA_LEN = 3*(NATOMC+NDEMAX)

      SR_LEN = 4*NAVMAX*NDIMER*NDEMAX

      EM_LEN = NDIMER*NDEMAX

      ND_LEN = NDIMER*NDEMAX

!!BEFORE WE SEARCH WE NEED TO UPDATE SOME IMPORTANT INFO IN ALL PROCESSES

      CALL MPI_BCAST(RA, RA_LEN, MPI_REAL8, 0, MPI_COMM_WORLD, IERR)

      CALL MPI_BCAST(NATOMC, 1, MPI_INT, 0, MPI_COMM_WORLD, IERR)

      CALL MPI_BCAST(PX, 7*NDEMAX, MPI_REAL8, 0, MPI_COMM_WORLD, IERR)

      CALL MPI_BCAST(NDTOT, 1, MPI_INT, 0, MPI_COMM_WORLD, IERR)

!!DO LOOP OVER ALL THE CURRENT DEFECT

      DO I=1, NDTOT

!!IF THE AV REQUIRES SADDLE POINT SEARCHING

        IF((PX(5,I).NE.0.0D0).AND.(PX(6,I).NE.0.0D0)) THEN

!!CARRY OUT SADDLE POINT SEARCHING

!!$OMP PARALLEL DO PRIVATE(ADIMER)

            DO J=LDIMERI,UDIMERI

            CALL DIMER_INIT(ADIMER,PX(1,I),PX(2,I),PX(3,I),PX(4,I))

            CALL DIMER_SEARCH(ADIMER)

!!$OMP CRITICAL

!!STORE INFORMATION,EM:MIGRATION ENERGY BARRIER
!!ND:TOTAL NUMBER OF ATOMS IN THE AV

            EM(J,I) = ADIMER%EMDIFF

            ND(J,I) = ADIMER%NATOM1

!!STORE THE SADDLE POINT CONFIGURATION,SR:POSITION OF SADDLE POINT

            DO K=1,ND(J,I)

              SR(1,K,J,I) = ADIMER%SX(1,K)
              SR(2,K,J,I) = ADIMER%SX(2,K)
              SR(3,K,J,I) = ADIMER%SX(3,K)
              SR(4,K,J,I) = ADIMER%SI(1,K)

            ENDDO

            WRITE(*,'(A,2I12,F14.5)') 'SPS',I,J,EM(J,I)

!!$OMP END CRITICAL

            CALL DIMER_FINISH(ADIMER)

          ENDDO

!!$OMP END PARALLEL DO

        ENDIF

      ENDDO

!!NOW ALL SR, EM, AND ND WILL GATHER TO RANK 0 SO IT CAN DO KMC

      CALL MPI_REDUCE(SR, SRBUFF, SR_LEN, MPI_REAL8, MPI_SUM, 0,
     &                MPI_COMM_WORLD, IERR)

      CALL MPI_REDUCE(EM, EMBUFF, EM_LEN, MPI_REAL8, MPI_SUM, 0,
     &                MPI_COMM_WORLD, IERR)

      CALL MPI_REDUCE(ND, NDBUFF, ND_LEN, MPI_INT, MPI_SUM, 0,
     &                MPI_COMM_WORLD, IERR)

      !SO THAT CODE CAN CONTINUE AS NORMAL

      IF (RANK .EQ. 0)THEN

        SR = SRBUFF

        EM = EMBUFF

        ND = NDBUFF

      ENDIF

!!STORE INFROMATION FOR SADDLE POINT RECYCLYING

      OE = EM

      ON = ND

      OR = SR

      OX = PX

      RETURN
      END SUBROUTINE
!***********************************************************************
!
!!THIS SUBROUTINE IS TO CHOSE AV BASED ON THE MIRGRATION ENERGY BARRIERS
!!THE CONVENTIONAL KINETIC MONTE CARLO ALGORITHM IS USED
!!The old KMC, but instead by a new one - Gao 2021
!! The problem : it is very likely that em(k,j) is 100 (failure in search
!! saddle point), and thus we neeed to remove these saddle points

      SUBROUTINE KMCold

      USE SYSTEM
      USE KMCVAR
      USE DEFECT

      IMPLICIT NONE

      INTEGER J,K,COUNTER

      REAL*8 RND

!!FQSUM:SUM OF FREQUENCY,POSUM:SUM OF PROBABILITIES(FOR CHECK ONLY)

      FQSUM = 0.0D0
      POSUM = 0.0D0

!!DO LOOP OVER ALL THE DEFECTS AND NUMBER OF SPS

      DO J=1,NDTOT
        DO K=1,NDIMER

!!CALCULATE FREQUENCY BASED ON THE MIGRATION ENERGY BARRIER
!!THE PREFACTOR HERE IS USED TO BE CONSTANT (NOT TRUE IN GENERAL)
!!THE PREFACTOR CALCUATED ON THE FLY IS IMPLEMENTED IN ANOTHER VERSION

            FQ(K,J) = 1.0D12*EXP(-EM(K,J)/KB/TEMP)
            FQSUM = FQSUM + FQ(K,J)
        ENDDO
      ENDDO

      CALL RANDOM_NUMBER(RND)

!!DETERMINE TIME STEP IN KMC

      TMEAN = -LOG10(RND)/FQSUM

      TIME = TIME + TMEAN

      DO J=1,NDTOT
        DO K = 1,NDIMER

!!DETERMINE THE PROBABILITIES FOR EACH EVENT

          PO(K,J) = FQ(K,J)/FQSUM
        ENDDO
      ENDDO

      DO J=1,NDTOT
        DO K=1,NDIMER
          POSUM = POSUM + PO(K,J)
          PO(K,J) = POSUM
        ENDDO
      ENDDO

      CALL RANDOM_NUMBER(RND)

      COUNTER = 0

!!DETERMINE WHICH EVENT IS CHOSEN

      DO J=1,NDTOT
        DO K=1,NDIMER
          IF(RND.GT.PO(K,J)) THEN
             COUNTER = COUNTER + 1
          ENDIF
        ENDDO
      ENDDO

!!JID:THE ID OF THE DEFECT THAT IS CHOSEN

      JID = CEILING(REAL(COUNTER+1)/REAL(NDIMER))

!!EID:THE ID OF THE EVENT THAT IS CHOSEN

      EID = COUNTER+1-(JID-1)*NDIMER

!      WRITE(*,*) 'JID:',JID,'EID:',EID,'RND',RND
      DO K = 1,NDIMER
        WRITE(21, *) (EM(K,J), J = 1, NDTOT)
      ENDDO

      RETURN

      END SUBROUTINE
!***********************************************************************!
!!STATIC RELAXATION OF THE SYSTEM

      SUBROUTINE OPT

      USE SYSTEM
      USE KMCVAR
      USE DEFECT

      IMPLICIT NONE

      INTEGER I,J,NSTMP

!!NEED TO WRITE INITIAL ENERGY BEFORE RECALCULATION
      WRITE(25,'(I3.3, F14.5)', ADVANCE = 'NO') ISTEP, ENCALC

            
!!THE TOTAL NUMBER OF THE ATOMS INSIDE THE CHOSEN AV
        NSTMP = ND(EID,JID)

!!PASSING THE POSITIONS OF ATOMS INSIDE AV TO THE SYSTEM

        DO I=1,NSTMP

          J=INT(SR(4,I,EID,JID))

          RA(1,J) = SR(1,I,EID,JID)
          RA(2,J) = SR(2,I,EID,JID)
          RA(3,J) = SR(3,I,EID,JID)

        ENDDO

!!OPTIMIZATION

      CALL OPT_CG_INIT(NATOMC)

      CALL OPT_CG

      CALL OPT_CG_FINISH

!!RECALCUATE THE FORCE AFTER OPTIMIZATION

      CALL EAM_FORCE

!!ISTEP:CURRENT NUMBER OF STEPS,NDTOT:TOTAL NUMBER OF DEFECTS
!!JID:DEFECT ID,EID:EVENT ID,EM:MIGRATION ENERGY BARRIER
!!ENCALC:THE SYSTEM ENERGY AFTER OPT,TIME:CURRENT TIME

        WRITE(*,'(3X, 4A12, 2A14, A14)') 'STEP', '# DEFECTS',
     & 'DEFECT ID', 'DIMER ID', 'BARRIER', 'SYS ENERGY', 'TIME'

        WRITE(*,'(A,4I12,2F14.5,E14.5)') 'ITR',ISTEP,NDTOT,JID,
     # EID,EM(EID,JID),ENCALC,TIME

!!THIS IS FOR WRITING SADDLE AND FINAL ENERGY AFTER INITIAL ENERGY
     
        WRITE(25, '(2F14.5)') EM(EID, JID), ENCALC

!!OUTPUT DEFECTS FOR OVITO

      CALL VISUAL
      IF (MOD(ISTEP,20) .EQ. 0) THEN
        CALL OVITO
      ENDIF

      RETURN

      END SUBROUTINE

!***********************************************************************!
!!THIS SUBROUTINE IS TO IDENTIFY DEFECT CLUSTERS IN THE SYSTEM
! Gao -2021

      SUBROUTINE DCIDNT
      USE SYSTEM
      USE DEFECT

      IMPLICIT NONE

      INTEGER I,J,NTEMP

      REAL*8 RX,RY,RZ,RS,PBC,PTEMP(7,NDEMAX)
      double precision :: dr
      integer, allocatable, dimension(:) :: idnewdef
      integer :: newdef, nvac1, nint1
      integer :: k, l, m, n, itf
      integer, allocatable, dimension(:) :: Neighbor
      integer, allocatable, dimension(:,:) :: NNeighbor
      integer, allocatable, dimension(:) :: iicn, iiticn
      integer :: iticn, icn,ijpcx,ijpcy,ijpcz

      OPEN(13,file='defecttemp.dat',status='replace')

      PTEMP = 0.0D0

!!OUTPUT THE DEFECT STRUCTURE TO ATOMTV AND OVITO

      CALL VISUAL

!!OUTPUT INITIAL DEFECTS WITH FILE INDEX 000
      IF (ISTEP .EQ. 1)THEN

        ISTEP = 0

        CALL OVITO

        ISTEP = 1

      END IF

!!INITIALIZATION,ND:NUMBER OF DEFECTS
!!VAC:VACANCY,INT:INTERSTITIAL

      NDVAC = 0
      NDINT = 0

!!INDEX FOR INTERSTITIAL AND VACACNY

      II = 0
      IV = 0

!!IDENTIFY INTERSTITIAL,PI:POSITION OF INTERSITITAL

      DO I=1,NATOMC
        IF(TI(I).NE.1) THEN
          NDINT = NDINT + 1
          II(1,NDINT) = I
          PI(1,NDINT) = RA(1,I)
          PI(2,NDINT) = RA(2,I)
          PI(3,NDINT) = RA(3,I)
        ENDIF
      ENDDO

!!IDENTIFY VACANCY,PV:POSITION OF VACANCY

      DO J=1,NATOMR
        IF(TV(J).NE.1) THEN
          NDVAC = NDVAC + 1
          IV(1,NDVAC) = J
          PV(1,NDVAC) = RR(1,J)
          PV(2,NDVAC) = RR(2,J)
          PV(3,NDVAC) = RR(3,J)
        ENDIF
      ENDDO

      NDTOT = 0

!!COUNT TOTAL NUMBER OF DEFECTS AND DETERMINE THEIR POSITION

      DO I=1,NDINT

        NDTOT = NDTOT + 1

        PX(1,NDTOT) = PI(1,I)
        PX(2,NDTOT) = PI(2,I)
        PX(3,NDTOT) = PI(3,I)

        PX(4,NDTOT) = 4.2D0 !radius of interstitial AV
        PX(5,NDTOT) = 1.0D0 !TELLS IF IN NEED OF SADDLE SEARCHES
        PX(6,NDTOT) = 1.0D0 !THIS MEANS IT IS A INTERSTITIAL

        PX(7,NDTOT) = II(1,I)

      ENDDO

      DO J=1,NDVAC

        NDTOT = NDTOT + 1

        PX(1,NDTOT) = PV(1,J)
        PX(2,NDTOT) = PV(2,J)
        PX(3,NDTOT) = PV(3,J)


        PX(4,NDTOT) = 2.7D0 !radius of vacancy AV
        PX(5,NDTOT) = 1.0D0 !TELLS IF IN NEED OF SADDLE SEARCHES
        PX(6,NDTOT) =-1.0D0 !THIS MEANS IT IS A VACANCY

        PX(7,NDTOT) = IV(1,J)

      ENDDO

      allocate(Neighbor(1:NDTOT))
      allocate(NNeighbor(1:NDTOT,1:NDTOT))

!!DETERMINE THE REAL DEFECTS

      DO I=1,NDTOT

       Neighbor(i)=0
       Neighbor(i)=Neighbor(i)+1
       NNeighbor(i,Neighbor(i))=i

        DO J=1,NDTOT

         IF(i==j)CYCLE
          RX=PBC(PX(1,I)-PX(1,J))
          RY=PBC(PX(2,I)-PX(2,J))
          RZ=PBC(PX(3,I)-PX(3,J))

          dr=sqrt(RX*RX+RY*RY+RZ*RZ)

          IF(PX(6,I)>0.0) THEN
           if(dr.LT.1.200D0) then
            Neighbor(i)=Neighbor(i)+1
            NNeighbor(i,Neighbor(i))=j
           endif

          ELSEIF(PX(6,I)<0.0)THEN                   !atom i is a vac
           if(PX(6,J)<0.0) then                     !atom j is a vac
            if(dr.lt.0.95)then
             Neighbor(i)=Neighbor(i)+1
             NNeighbor(i,Neighbor(i))=j
            endif
           endif

           if(PX(6,J)>0.0) then                     !atom j is a int
            if(dr.le.1.200d0)then
             Neighbor(i)=Neighbor(i)+1
             NNeighbor(i,Neighbor(i))=j
            endif
           endif
          ENDIF

        ENDDO
       ENDDO

       nvac1=0
       nint1=0

       DO i=1,NDTOT
        IF(Neighbor(i)==1)THEN
         if(PX(6,I)<0.0)then
          nvac1=nvac1+1
         elseif(PX(6,I)>0.0)then
          nint1=nint1+1
         endif
        ENDIF
       ENDDO

       write(6,"(a35,i8)")"cluster :: 1_vac :: numbers=",nvac1

       DO i=1,NDTOT
        IF(Neighbor(i)==1)THEN
         if(PX(6,i)<0.0)then
          write(22,"(f8.1,3(2x,f15.8))") PX(6,i),PX(1,i),PX(2,i),PX(3,i)
         endif
        ENDIF
       ENDDO

       write(6,"(a35,i8)")"cluster :: 1_int :: numbers=",nint1

       DO i=1,NDTOT
        IF(Neighbor(i)==1)THEN
         if(PX(6,i)>0.0)then
          write(22,"(f8.1,3(2x,f15.8))") PX(6,i),PX(1,i),PX(2,i),PX(3,i)
         endif
        ENDIF
       ENDDO

!******** END write single vac and single int to result.dat
!******** left defects after remove single defect
       newdef=NDTOT-nvac1-nint1

       DO I=1,NDTOT

        if(Neighbor(i).gt.1)then
         write(13,"(i8)")i
        endif

       ENDDO

!********  now find clusters
! DO WHILE

      do while( newdef > 0 )

       if(newdef>0) then
        REWIND(13)
        allocate(idnewdef(1:newdef))

        do i=1,newdef
         read(13,*)idnewdef(i)
        enddo

        close(13)
       endif

          allocate(iicn(1:100*newdef))
          allocate(iiticn(1:100*newdef))

          icn=0
          iicn=0

          DO i=1,Neighbor(idnewdef(1))
           icn=icn+1
           iicn(icn)=NNeighbor(idnewdef(1),i)
          ENDDO

!       DO i=1,10  !!add this cycle to make sure find all belong to CL
       DO i=1,4  !!add this cycle to make sure find all belong to CL
        DO j=2,newdef

           n=0

         do k=1,Neighbor(idnewdef(j))  ! J neighbor number
          do l=1,icn
          if(NNeighbor(idnewdef(j),k)==iicn(l))then
!*****************************************************************
! If any of J int or J neighbor ints is the neighbor of current 1*
! int, all the neighbors of J int belong to 1 cluster.           *
!*****************************************************************
           n=1
          endif
         enddo
        enddo

         IF(n==1)THEN
          do m=1,Neighbor(idnewdef(j))
           icn=icn+1
           iicn(icn)=NNeighbor(idnewdef(j),m)
          enddo
         ENDIF

        ENDDO
       ENDDO

       iticn=0
       iiticn=0

!*****************************************************************
! Delete the interstitial that was repeated many times due to    *
! the above loop - DO i=1,10                                     *
! iticn - number of defects in the current 1 cluster             *
! iiticn - remeber defect index in the current 1 cluster         *
!*****************************************************************

       do j=1,icn
          itf=0
          IF(j==1)THEN
           iticn=iticn+1
           iiticn(iticn)=iicn(j)
          ELSE
           do k=1,iticn
              if(iicn(j)==iiticn(k))then
                itf=1
              endif
           enddo
           if(itf==0)then
            iticn=iticn+1
            iiticn(iticn)=iicn(j)
           endif
          ENDIF
       enddo

! we need to make chnage in PX(5,I) and PX(6,J)
       ijpcx = 0
       ijpcy = 0
       ijpcz = 0

       do K=2,iticn

        I = iiticn(1)
        J = iiticn(K)

        if((px(1,i)<1.0.and.px(1,j)>rbox-1.0).or.
     >     (px(1,j)<1.0.and.px(1,i)>rbox-1.0)) then
           ijpcx = 1
        end if

        if((px(2,i)<1.0.and.px(2,j)>rbox-1.0).or.
     >     (px(2,j)<1.0.and.px(2,i)>rbox-1.0)) then
           ijpcy = 1
        end if

        if((px(3,i)<1.0.and.px(3,j)>rbox-1.0).or.
     >     (px(3,j)<1.0.and.px(3,i)>rbox-1.0)) then
           ijpcz = 1
        end if

       enddo

       if(ijpcx==1.or.ijpcy==1.or.ijpcz==1) then

        do K=1,iticn
         J = iiticn(K)
         PTEMP(:,K) = PX(:,J)
        enddo

        do K=1,iticn

         if(ijpcx==1) then
          if(PTEMP(1,K)<RBH1) then
           PTEMP(1,K) = PTEMP(1,K) + RBH1
          else
           PTEMP(1,K) = PTEMP(1,K) - RBH1
          end if
         endif

         if(ijpcy==1) then
          if(PTEMP(2,K)<RBH1) then
           PTEMP(2,K) = PTEMP(2,K) + RBH1
          else
           PTEMP(2,K) = PTEMP(2,K) - RBH1
          end if
         endif

         if(ijpcz==1) then
          if(PTEMP(3,K)<RBH1) then
           PTEMP(3,K) = PTEMP(3,K) + RBH1
          else
           PTEMP(3,K) = PTEMP(3,K) - RBH1
          end if
         endif

        enddo

        do K=2,iticn

            I = 1
            J = K

            PTEMP(1,I) = PTEMP(1,I)*ABS(PTEMP(6,I))+
     >                   PTEMP(1,J)*ABS(PTEMP(6,J))
            PTEMP(2,I) = PTEMP(2,I)*ABS(PTEMP(6,I))+
     >                   PTEMP(2,J)*ABS(PTEMP(6,J))
            PTEMP(3,I) = PTEMP(3,I)*ABS(PTEMP(6,I))+
     >                   PTEMP(3,J)*ABS(PTEMP(6,J))

            PTEMP(1,I) = PTEMP(1,I)/(ABS(PTEMP(6,I))+ABS(PTEMP(6,J)))
            PTEMP(2,I) = PTEMP(2,I)/(ABS(PTEMP(6,I))+ABS(PTEMP(6,J)))
            PTEMP(3,I) = PTEMP(3,I)/(ABS(PTEMP(6,I))+ABS(PTEMP(6,J)))

            PTEMP(6,I) = PTEMP(6,I)+PTEMP(6,J)
            PTEMP(5,J) = 0.0D0
        enddo

        do K=1,iticn

         if(ijpcx==1) then
          if(PTEMP(1,K)<RBH1) then
           PTEMP(1,K) = PTEMP(1,K) + RBH1
          else
           PTEMP(1,K) = PTEMP(1,K) - RBH1
          end if
         endif

         if(ijpcy==1) then
          if(PTEMP(2,K)<RBH1) then
           PTEMP(2,K) = PTEMP(2,K) + RBH1
          else
           PTEMP(2,K) = PTEMP(2,K) - RBH1
          end if
         endif

         if(ijpcz==1) then
          if(PTEMP(3,K)<RBH1) then
           PTEMP(3,K) = PTEMP(3,K) + RBH1
          else
           PTEMP(3,K) = PTEMP(3,K) - RBH1
          end if
         endif

        enddo

        do K=1,iticn
         J = iiticn(K)
         PX(:,J) = PTEMP(:,K)
        enddo

       else

        do K=2,iticn

            I = iiticn(1)
            J = iiticn(K)

            PX(1,I) = PX(1,I)*ABS(PX(6,I))+PX(1,J)*ABS(PX(6,J))
            PX(2,I) = PX(2,I)*ABS(PX(6,I))+PX(2,J)*ABS(PX(6,J))
            PX(3,I) = PX(3,I)*ABS(PX(6,I))+PX(3,J)*ABS(PX(6,J))

            PX(1,I) = PX(1,I)/(ABS(PX(6,I))+ABS(PX(6,J)))
            PX(2,I) = PX(2,I)/(ABS(PX(6,I))+ABS(PX(6,J)))
            PX(3,I) = PX(3,I)/(ABS(PX(6,I))+ABS(PX(6,J)))

            if(px(1,i)<0.0) px(1,i) = 0.0
            if(px(2,i)<0.0) px(2,i) = 0.0
            if(px(3,i)<0.0) px(3,i) = 0.0

            PX(6,I) = PX(6,I)+PX(6,J)
            PX(5,J) = 0.0D0
        enddo
       endif

       OPEN(13,file='defecttemp.dat',status='replace')


       do i=1,newdef
           k=0

           do j=1,iticn
            if(idnewdef(i)==iiticn(j))then
             k=1
            end if
           enddo

           if(k==0)then
            write(13,"(i8)")idnewdef(i)
           endif

       enddo

           newdef=newdef-iticn

           deallocate(idnewdef)
           deallocate(iicn)
           deallocate(iiticn)

      enddo

      close(13)

           deallocate(Neighbor)
           deallocate(NNeighbor)

      NTEMP = 0
      PTEMP = 0.0D0

      DO I=1,NDTOT

        IF((PX(5,I).NE.0.0D0).AND.(PX(6,I).NE.0.0D0)) THEN

        NTEMP = NTEMP + 1

        PTEMP(:,NTEMP) = PX(:,I)

        ENDIF

      ENDDO

      NDTOT = NTEMP

      PX = 0.0D0

!!RESORTING,PX STORES THE POSITIONS OF THE ACTIVE VOLUMES

      DO I=1,NDTOT

        PX(:,I) = PTEMP(:,I)

      WRITE(*,'(A,I12,6(F14.5),I12)') 'DCIDNT',I,PX(1:6,I),
     #INT(PX(7,I))

      ENDDO

      ONDTOT = NDTOT

!!IF NO DEFECT IS PRESENT, STOP THE PROGRM

      IF(NDTOT.EQ.0) THEN

        WRITE(*,*) 'NO DEFECT IN SYSTEM,STOP'

        STOP

      ENDIF

      RETURN
      END SUBROUTINE

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

      IMPLICIT NONE

      REAL*8 RX,RY,RZ,RS,PBC,PTEMP(7,NDEMAX)
      REAL*8 rnd,speed,dt,scalef,sump,prop
      integer :: i,j,k,jmapsz,ijump
      REAL*8, allocatable, dimension(:,:) :: pj

      SKIP = .FALSE.

      speed = 0.0d0
      jmapsz = 1000

      allocate(pj(1:2,1:ndtot))

      do i= 1,ndtot
       j = abs(px(6,i))
       if(px(6,i).lt.0) then
        k=1
       else
        k=2
       endif

       if(j.lt.50) then
        pj(k,j) = 1.0D12*exp(-emig(k,j)/kb/temp)
       else
        pj(k,j) = 1.0D12*exp(-emig(k,50)/kb/temp)
       endif

       speed = speed + pj(k,j)
      enddo

      dt = 1.0/speed
      scalef = dt*float(jmapsz)

      CALL RANDOM_NUMBER(RND)
      write(6,*) rnd
      prop = rnd*float(jmapsz)

      sump = 0.0

      do i = 1,ndtot
       j = abs(px(6,i))
       if(px(6,i).lt.0) then
        k=1
       else
        k=2
       endif

       sump = sump + pj(k,j)*scalef
       if(sump.ge.prop) goto 200
      enddo

 200  continue

      ijump = i
      NDTOT = 1

      PTEMP(1,1) = px(1,ijump)
      PTEMP(2,1) = px(2,ijump)
      PTEMP(3,1) = px(3,ijump)
      PTEMP(4,1) = px(4,ijump)
      PTEMP(5,1) = px(5,ijump)
      PTEMP(6,1) = px(6,ijump)
      PTEMP(7,1) = px(7,ijump)

      PX = 0.0D0

!!RESORTING,PX STORES THE POSITIONS OF THE ACTIVE VOLUMES

      DO I=1,NDTOT

        PX(:,I) = PTEMP(:,I)

      WRITE(*,'(A,I12,6(F14.5),I12)') 'DFSEL',I,PX(1:6,I),
     #INT(PX(7,I))

      ENDDO

      deallocate(pj)

      RETURN

      END SUBROUTINE


!***********************************************************************!
!!THIS SUBROUTINE IS TO CHOSE AV BASED ON THE MIRGRATION ENERGY BARRIERS
!!THE CONVENTIONAL KINETIC MONTE CARLO ALGORITHM IS USED

      SUBROUTINE KMC

      USE SYSTEM
      USE KMCVAR
      USE DEFECT

      IMPLICIT NONE

      INTEGER I,J,K,COUNTER,jmapsz

      REAL*8 RND,dt,scalef,prop

!!FQSUM:SUM OF FREQUENCY,POSUM:SUM OF PROBABILITIES(FOR CHECK ONLY)

      FQSUM = 0.0D0
      POSUM = 0.0D0
      jmapsz = 1000

!!DO LOOP OVER ALL THE DEFECTS AND NUMBER OF SPS

      DO J=1,NDTOT
        DO K=1,NDIMER

!!CALCULATE FREQUENCY BASED ON THE MIGRATION ENERGY BARRIER
!!THE PREFACTOR HERE IS USED TO BE CONSTANT (NOT TRUE IN GENERAL)
!!THE PREFACTOR CALCUATED ON THE FLY IS IMPLEMENTED IN ANOTHER VERSION
!! Now considering good dimer search - Gao 2021
           if(em(k,j).lt.10) then
            FQ(K,J) = 1.0D12*EXP(-EM(K,J)/KB/TEMP)
            FQSUM = FQSUM + FQ(K,J)
           endif
        ENDDO
      ENDDO

      CALL RANDOM_NUMBER(RND)
      prop = rnd*float(jmapsz)

!!DETERMINE TIME STEP IN KMC

      if(FQSUM.gt.0.0) then
      TMEAN = -LOG10(RND)/FQSUM

      TIME = TIME + TMEAN

      dt = 1.0/fqsum
      scalef = dt*dfloat(jmapsz)
      endif

      DO J=1,NDTOT
        DO K=1,NDIMER
         if(em(k,j).lt.10.0) then
          POSUM = POSUM + FQ(K,J)*scalef
         endif
        if(posum.ge.prop) goto 200
        ENDDO
      ENDDO

 200  continue
!!JID:THE ID OF THE DEFECT THAT IS CHOSEN

      JID = J

!!EID:THE ID OF THE EVENT THAT IS CHOSEN

      EID = K

!!UPDATING THE PRE-SELECTION ENERGIES

      IF(PX(6,1).LT.0) THEN
      K=1
      ELSE
      K=2
      ENDIF

      I = ABS(INT(PX(6,1)))

      EMIG(K,I) = (EMIG(K,I)*EMIC(K,I) + EM(EID,JID)) / (EMIC(K,I)+1)

      EMIC(K,I) = EMIC(K,I) + 1

!!CHECKING FOR TOO HIGH TRANSITION ENERGY

      IF (EM(EID,JID) .GT. 0.8) THEN

        SKIP = .TRUE.

        WRITE(6,*) 'MIGRATION ENERGY GREATER THAN 0.8eV -> SKIPPING'

      ENDIF

!!WRITING THE DIMER ENERGIES TO STD_KDIM.OUT

      WRITE(6,*) 'JID:',JID,'EID:',EID,'RND',RND
      DO K = 1,NDIMER
        WRITE(21, *) (EM(K,J), J = 1, NDTOT)
      ENDDO

      RETURN

      END SUBROUTINE
!***********************************************************************!
!!WILL INTRODUCE MORE DEFECTS AFTER A SPECIFIED AMOUNT OF TIME THAT MATCHES
!!EXPERIMENTAL DOSE RATE
      SUBROUTINE DOSE

      USE SYSTEM
      !!TIME - CURRENT ACCUMILATED KMC TIME
      !!TIME_LAST - TIME OF LAST DEFECT INTRODUCTION
      !!DTIME - PRECALCULATED TIME FOR NEXT DEFECT INTRODUCTION (CONSTANT)

      IMPLICIT NONE

      REAL*8 SIM_DOSE,DOSE_RATE, DTIME

      SIM_DOSE = 100.d0/FLOAT(NATOMC) !NUMBER OF DEFECTS CREATED/TOTAL ATOMS
      DOSE_RATE = 10E-1 !dpa/s
      DTIME = SIM_DOSE/DOSE_RATE
      !DTIME = 1E-16 !RANDOM LOW TIME TO TEST BUGS

      IF(TIME - TIME_LAST .GE. DTIME) THEN

        WRITE(*,*) 'INITIATING FRENKEL GENERATOR'
        CALL FRENKELGENERATOR(INT(SIM_DOSE*NATOMC))
        WRITE(*,*) 'TIME=',TIME,'TIME_LAST=',TIME_LAST,'DTIME=',DTIME
        TIME_LAST = TIME

      ENDIF

      END SUBROUTINE
