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

      SUBROUTINE OPLD

      USE MPI
      USE LMP_DRIVER

      USE SYSTEM,ONLY:BSTEP,NTSTEP,ISTEP,RANK,NATOMC,RA,SKIP
      USE DEFECT,ONLY:UPNEIGH
      USE LMP_DRIVER,ONLY:LMP_STRESS_RELAX

      IMPLICIT NONE

      INTEGER I,J

      DO ISTEP = (BSTEP+1), (NTSTEP+BSTEP)

        CALL DOSE

        IF (MOD(ISTEP-1,50).EQ.0) UPNEIGH = .TRUE.
        CALL FIND_DEFECT

        IF (RANK .EQ. 0) THEN
          CALL DCIDNT
          CALL DFSEL
          CALL CBC
        END IF    

        CALL SPS
        CALL KMC
        IF (.NOT. SKIP) CALL OPT
        CALL RESTART

      ENDDO

      END SUBROUTINE
!***********************************************************************!

!THIS SUBROUTINE IS TO CHECK THE BOURNDARY CONDITION
!ONLY FOR ABSORBING BOUNDARY CONDITION, NOT FOR PBC

      SUBROUTINE CBC

      USE SYSTEM
      USE DEFECT
      USE LMP_DRIVER

      INTEGER I

      LOGICAL TOUT,TRELAX

      TRELAX = .FALSE.

!IF ABSORBING BOUNDARY CONDTION IS USED

      IF(TPBC.EQ.2) THEN

        DO I=1,NDTOT

          TOUT = .FALSE.

!IF THE DEFECT IS CLOSE TO THE BOUNDARY (LESS THAN 1 LATTICE CONSTANT) 
!IT IS CONSIDERED TO BE REMOVED FROM THE SYSTEM

         IF((PX(1,I).GT.(RBOX*EPSX-1.0D0)).OR.(PX(1,I).LT.1.0D0)) TOUT=.TRUE.
         IF((PX(2,I).GT.(RBOX*EPSY-1.0D0)).OR.(PX(2,I).LT.1.0D0)) TOUT=.TRUE.
         IF((PX(3,I).GT.(RBOX*EPSZ-1.0D0)).OR.(PX(3,I).LT.1.0D0)) TOUT=.TRUE.

          IF(TOUT) THEN

            TRELAX = .TRUE.

!IF THE DEFECT IS A VACANCY, ADD ONE ATOM AT VACANCY POSITION
!THE SYSTEM ENERGY WILL NOT BE CONSERVED

            IF(PX(6,I).EQ.(-1.0D0)) THEN

              NATOMC = NATOMC + 1

              RA(1,NATOMC) = PX(1,I)
              RA(2,NATOMC) = PX(2,I)
              RA(3,NATOMC) = PX(3,I)

              PX(5,I) = 0.0D0

!WRITE INFORMATION:WHICH VACANCY IS REMOVED

              WRITE(*,*) 'VAC OBC',PX(1,I),PX(2,I),PX(3,I)

!IF THE DEFECT IS AN INTERSITITAL, REMOVE THE ATOM FROM THE SYSTEM

            ELSE IF(PX(6,I).EQ.(1.0D0)) THEN

              RA(1,DINT(PX(7,I))) = RA(1,NATOMC)
              RA(2,DINT(PX(7,I))) = RA(2,NATOMC)
              RA(3,DINT(PX(7,I))) = RA(3,NATOMC)

              NATOMC = NATOMC -1

              PX(5,I) = 0.0D0

!WRITE INFORMATION:WHICH INTERSTITIAL IS REMOVED

              WRITE(*,*) 'INT OBC',PX(1,I),PX(2,I),PX(3,I)

            ENDIF

          ENDIF

        ENDDO

!RE-OPTIMIZE THE STRUCTURE USING CONJUGATE GRADIENT

        IF(TRELAX) THEN

          CALL SEND_OPLD_POS

          CALL LMP_RELAX

          CALL RECV_LMP_POS

          CALL DCIDNT

        ENDIF

      ENDIF

      RETURN

      END SUBROUTINE CBC
!***********************************************************************!
!THIS SUBROUTINE IS TO PERFORM SADDLE POINT SEARCHING

      SUBROUTINE SPS

      USE MPI

      USE SYSTEM
      USE DEFECT
      USE KMCVAR

      USE CLASS_DIMER

      IMPLICIT NONE

      TYPE(DIMER) :: ADIMER

      INTEGER :: J,K

      INTEGER :: RA_LEN, SR_LEN

!CLEAR ALL DATA

      SR = 0.0D0
      IF (RANK .EQ. 0) SRBUFF = 0.0D0
      EM = 0.0D0
      IF (RANK .EQ. 0) EMBUFF = 0.0D0
      ND = 0.0D0
      IF (RANK .EQ. 0) NDBUFF = 0.0D0
!MPI VARIABLES

      RA_LEN = SIZE(RA)

      SR_LEN = SIZE(SR)

!BEFORE WE SEARCH WE NEED TO UPDATE SOME IMPORTANT INFO IN ALL PROCESSES

      CALL MPI_BCAST(RA, RA_LEN, MPI_REAL8, 0, MPI_COMM_WORLD, IERR)

      CALL MPI_BCAST(NATOMC, 1, MPI_INT, 0, MPI_COMM_WORLD, IERR)

      CALL MPI_BCAST(PX, 7*NDEMAX, MPI_REAL8, 0, MPI_COMM_WORLD, IERR)

!IF THE AV REQUIRES SADDLE POINT SEARCHING

      IF((PX(5,1).NE.0.0D0).AND.(PX(6,1).NE.0.0D0)) THEN

!CARRY OUT SADDLE POINT SEARCHING

        DO J=LDIMERI,UDIMERI

          CALL DIMER_INIT(ADIMER,PX(1,1),PX(2,1),PX(3,1),PX(4,1))

          CALL DIMER_SEARCH(ADIMER)

!STORE INFORMATION,EM:MIGRATION ENERGY BARRIER
!ND:TOTAL NUMBER OF ATOMS IN THE AV

          EM(J) = ADIMER%EMDIFF

          ND(J) = ADIMER%NATOM1

!STORE THE SADDLE POINT CONFIGURATION,SR:POSITION OF SADDLE POINT

          DO K=1,ND(J)

            SR(1,K,J) = ADIMER%SX(1,K)
            SR(2,K,J) = ADIMER%SX(2,K)
            SR(3,K,J) = ADIMER%SX(3,K)
            SR(4,K,J) = ADIMER%SI(1,K)

          ENDDO

          IF (MOD(ISTEP,NFREQ) .EQ. 0) THEN
            WRITE(*,'(A,1I12,F14.5)') 'SPS',J,EM(J)
          ENDIF

          CALL DIMER_FINISH(ADIMER)

        ENDDO

      ENDIF

!GATHER RESULTS ON MAIN RANK

      CALL MPI_REDUCE(SR, SRBUFF, SR_LEN, MPI_REAL8, MPI_SUM, 0,MPI_COMM_WORLD, IERR)

      CALL MPI_REDUCE(EM, EMBUFF, NDIMER, MPI_REAL8, MPI_SUM, 0,MPI_COMM_WORLD, IERR)

      CALL MPI_REDUCE(ND, NDBUFF, NDIMER, MPI_INT, MPI_SUM, 0,MPI_COMM_WORLD, IERR)

      !SO THAT CODE CAN CONTINUE AS NORMAL

      IF (RANK .EQ. 0)THEN

        SR = SRBUFF

        EM = EMBUFF

        ND = NDBUFF

      ENDIF
      
      RETURN
      END SUBROUTINE
!***********************************************************************!

!STATIC RELAXATION OF THE SYSTEM

      SUBROUTINE OPT

      USE SYSTEM
      USE KMCVAR
      USE DEFECT
      USE LMP_DRIVER
      USE MPI

      IMPLICIT NONE

      INTEGER I,J,NSTMP

      IF (RANK.EQ.0) THEN
            
        !THE TOTAL NUMBER OF THE ATOMS INSIDE THE CHOSEN AV

        NSTMP = ND(EID)

        !PASSING THE POSITIONS OF ATOMS INSIDE AV TO THE SYSTEM

        DO I=1,NSTMP
  
          J=INT(SR(4,I,EID))
    
          RA(1,J) = SR(1,I,EID)
          RA(2,J) = SR(2,I,EID)
          RA(3,J) = SR(3,I,EID)
  
        ENDDO

      END IF

      CALL SEND_OPLD_POS

      CALL LMP_RELAX

      CALL RECV_LMP_POS

      IF (RANK .EQ. 0) THEN

        !ISTEP:CURRENT NUMBER OF STEPS,NDTOT:TOTAL NUMBER OF DEFECTS
        !JID:DEFECT ID,EID:EVENT ID,EM:MIGRATION ENERGY BARRIER
        !ENCALC:THE SYSTEM ENERGY AFTER OPT,TIME:CURRENT TIME

        IF (MOD(ISTEP,NFREQ) .EQ. 0) THEN
        WRITE(*,'(3X, 2A12, 7A14)') 'STEP','DIMER ID','BARRIER',
     #    'SYS ENERGY', 'TIME', 'Px', 'Py', 'Pz', 'DPA'

        WRITE(*,'(A,2I12,2F14.5,4E14.5,E14.5)') 'ITR',ISTEP,EID,EM(EID),
     #    ENCALC,TIME,PRESS(1),PRESS(2),PRESS(3),DPA
        ENDIF

      END IF
      
        !OUTPUT DEFECTS FOR OVITO

        CALL FIND_DEFECT

        IF (MOD(ISTEP,NFREQ) .EQ. 0 .AND. RANK .EQ. 0) CALL OVITO

      RETURN

      END SUBROUTINE

!***********************************************************************!
!THIS SUBROUTINE IS TO IDENTIFY DEFECT CLUSTERS IN THE SYSTEM
! Gao -2021
! Modified: Angel Chavira 9/1/2022

      SUBROUTINE DCIDNT
      USE SYSTEM
      USE DEFECT

      IMPLICIT NONE

      INTEGER I,J,NTEMP

      REAL*8 RX,RY,RZ,PBC,PTEMP(7,NDEMAX)
      double precision :: dr
      integer, allocatable, dimension(:) :: idnewdef
      integer :: newdef, nvac1, nint1, sumInt, sumVac
      integer :: k, l, m, n, itf
      integer, allocatable, dimension(:) :: Neighbor
      integer, allocatable, dimension(:,:) :: NNeighbor
      integer, allocatable, dimension(:) :: iicn, iiticn
      integer :: iticn, icn,ijpcx,ijpcy,ijpcz, comdc

      OPEN(13,file='defecttemp.dat',status='replace')

      PTEMP = 0.0D0

!OUTPUT INITIAL DEFECTS WITH FILE INDEX 000
      IF (ISTEP .EQ. 1)THEN

        ISTEP = 0
        CALL OVITO
        ISTEP = 1

      END IF

!INITIALIZATION,ND:NUMBER OF DEFECTS
!VAC:VACANCY,INT:INTERSTITIAL

      NDVAC = 0
      NDINT = 0

!INDEX FOR INTERSTITIAL AND VACACNY

      II = 0
      IV = 0

!IDENTIFY INTERSTITIAL,PI:POSITION OF INTERSITITAL

      DO I=1,NATOMC
        IF(TI(I).NE.1) THEN
          NDINT = NDINT + 1
          II(1,NDINT) = I
          PI(1,NDINT) = RA(1,I)
          PI(2,NDINT) = RA(2,I)
          PI(3,NDINT) = RA(3,I)
        ENDIF
      ENDDO

!IDENTIFY VACANCY,PV:POSITION OF VACANCY

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

!COUNT TOTAL NUMBER OF DEFECTS AND DETERMINE THEIR POSITION

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

      !print *, ndtot, NDVAC, NDINT

      allocate(Neighbor(1:NDTOT))
      allocate(NNeighbor(1:NDTOT,1:NDTOT))

!DETERMINE THE REAL DEFECTS

      DO I=1,NDTOT

       Neighbor(i)=0
       Neighbor(i)=Neighbor(i)+1
       NNeighbor(i,Neighbor(i))=i

        DO J=1,NDTOT

         IF(i==j)CYCLE
          RX=PBC(PX(1,I)-PX(1,J),1)
          RY=PBC(PX(2,I)-PX(2,J),2)
          RZ=PBC(PX(3,I)-PX(3,J),3)

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

       DO i=1,NDTOT
        IF(Neighbor(i)==1)THEN
         if(PX(6,i)<0.0)then
          write(22,"(f8.1,3(2x,f15.8))") PX(6,i),PX(1,i),PX(2,i),PX(3,i)
         endif
        ENDIF
       ENDDO

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

!       DO i=1,10  !add this cycle to make sure find all belong to CL
       DO i=1,4  !add this cycle to make sure find all belong to CL
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

        if((px(1,i)<1.0.or.px(1,j)>rbox*EPSX-1.0).or.
     >     (px(1,j)<1.0.or.px(1,i)>rbox*EPSX-1.0)) then
           ijpcx = 1
        end if

        if((px(2,i)<1.0.or.px(2,j)>rbox*EPSY-1.0).or.
     >     (px(2,j)<1.0.or.px(2,i)>rbox*EPSY-1.0)) then
           ijpcy = 1
        end if

        if((px(3,i)<1.0.or.px(3,j)>rbox*EPSZ-1.0).or.
     >     (px(3,j)<1.0.or.px(3,i)>rbox*EPSZ-1.0)) then
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
          if(PTEMP(1,K)<RBH1*EPSX) then
           PTEMP(1,K) = PTEMP(1,K) + RBH1*EPSX
          else
           PTEMP(1,K) = PTEMP(1,K) - RBH1*EPSX
          end if
         endif

         if(ijpcy==1) then
          if(PTEMP(2,K)<RBH1*EPSY) then
           PTEMP(2,K) = PTEMP(2,K) + RBH1*EPSY
          else
           PTEMP(2,K) = PTEMP(2,K) - RBH1*EPSY
          end if
         endif

         if(ijpcz==1) then
          if(PTEMP(3,K)<RBH1*EPSZ) then
           PTEMP(3,K) = PTEMP(3,K) + RBH1*EPSZ
          else
           PTEMP(3,K) = PTEMP(3,K) - RBH1*EPSZ
          end if
         endif

        enddo

        comdc = 1

        do K=2,iticn

            I = 1
            J = K

            PTEMP(1,I) = PTEMP(1,I)*ABS(comdc)+
     >                   PTEMP(1,J)*ABS(PTEMP(6,J))
            PTEMP(2,I) = PTEMP(2,I)*ABS(comdc)+
     >                   PTEMP(2,J)*ABS(PTEMP(6,J))
            PTEMP(3,I) = PTEMP(3,I)*ABS(comdc)+
     >                   PTEMP(3,J)*ABS(PTEMP(6,J))

            PTEMP(1,I) = PTEMP(1,I)/(ABS(comdc)+ABS(PTEMP(6,J)))
            PTEMP(2,I) = PTEMP(2,I)/(ABS(comdc)+ABS(PTEMP(6,J)))
            PTEMP(3,I) = PTEMP(3,I)/(ABS(comdc)+ABS(PTEMP(6,J)))

            PTEMP(6,I) = PTEMP(6,I)+PTEMP(6,J)
            PTEMP(5,J) = 0.0D0

            comdc = comdc + 1
        enddo

        do K=1,iticn

         if(ijpcx==1) then
          if(PTEMP(1,K)<RBH1*EPSX) then
           PTEMP(1,K) = PTEMP(1,K) + RBH1*EPSX
          else
           PTEMP(1,K) = PTEMP(1,K) - RBH1*EPSX
          end if
         endif

         if(ijpcy==1) then
          if(PTEMP(2,K)<RBH1*EPSY) then
           PTEMP(2,K) = PTEMP(2,K) + RBH1*EPSY
          else
           PTEMP(2,K) = PTEMP(2,K) - RBH1*EPSY
          end if
         endif

         if(ijpcz==1) then
          if(PTEMP(3,K)<RBH1*EPSZ) then
           PTEMP(3,K) = PTEMP(3,K) + RBH1*EPSZ
          else
           PTEMP(3,K) = PTEMP(3,K) - RBH1*EPSZ
          end if
         endif

        enddo

        do K=1,iticn
         J = iiticn(K)
         PX(:,J) = PTEMP(:,K)
        enddo

       else

       comdc = 1

        do K=2,iticn

            I = iiticn(1)
            J = iiticn(K)

            PX(1,I) = PX(1,I)*ABS(comdc)+PX(1,J)*ABS(PX(6,J))
            PX(2,I) = PX(2,I)*ABS(comdc)+PX(2,J)*ABS(PX(6,J))
            PX(3,I) = PX(3,I)*ABS(comdc)+PX(3,J)*ABS(PX(6,J))

            PX(1,I) = PX(1,I)/(ABS(comdc)+ABS(PX(6,J)))
            PX(2,I) = PX(2,I)/(ABS(comdc)+ABS(PX(6,J)))
            PX(3,I) = PX(3,I)/(ABS(comdc)+ABS(PX(6,J)))

            if(px(1,i)<0.0) px(1,i) = 0.0
            if(px(2,i)<0.0) px(2,i) = 0.0
            if(px(3,i)<0.0) px(3,i) = 0.0

            PX(6,I) = PX(6,I)+PX(6,J)
            PX(5,J) = 0.0D0

            comdc = comdc + 1
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

      nvac1 = 0
      nint1 = 0

      DO I=1,NDTOT

        IF((PX(5,I).NE.0.0D0).AND.(PX(6,I).NE.0.0D0)) THEN

          NTEMP = NTEMP + 1

          PX(4,I) = 3.334 + 0.866*abs(PX(6,I))**(1./3.)

          PTEMP(:,NTEMP) = PX(:,I)
          
          IF(PX(6,I).EQ.-1) nvac1 = nvac1 + 1
          IF(PX(6,I).EQ.1) nint1 = nint1 + 1
        ENDIF

      ENDDO

      NDTOT = NTEMP

      PX = 0.0D0

      IF (MOD(ISTEP,NFREQ) .EQ. 0) write(6,"(a35,i8)")"cluster :: 1_vac :: numbers=",nvac1
      IF (MOD(ISTEP,NFREQ) .EQ. 0)write(6,"(a35,i8)")"cluster :: 1_int :: numbers=",nint1
!RESORTING,PX STORES THE POSITIONS OF THE ACTIVE VOLUMES

      NDSTOT = 0
      sumInt = 0; sumVac = 0

      DO I=1,NDTOT

        PX(:,I) = PTEMP(:,I)

        IF (MOD(ISTEP,NFREQ) .EQ. 0) THEN
          WRITE(*,'(A,I12,5(F14.5))') 'DCIDNT',I,PX(1:4,I),PX(6,I)
        ENDIF
        IF (PX(6,I) .GT. 0) sumInt = sumInt + int(abs(PX(6,I)))
        IF (PX(6,I) .LT. 0) sumVac = sumVac + int(abs(PX(6,I)))

      ENDDO

      NDSTOT = MAX(sumInt, sumVac)

      IF(CALC) CALL MD_DPA

      ONDTOT = NDTOT

      ONDSTOT = NDSTOT

!IF NO DEFECT IS PRESENT, STOP THE PROGRM

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

      REAL*8 PTEMP(7,1)
      REAL*8 rnd,speed,dt,scalef,sump,prop
      integer :: i,j,k,jmapsz,ijump
      REAL*8, allocatable, dimension(:,:) :: pj

      SKIP = .FALSE.

      speed = 0.0d0
      jmapsz = 1000

      allocate(pj(1:2,1:ndtot))

      do i= 1,ndtot
       j = int(abs(px(6,i)))
       if(px(6,i).lt.0) then
        k=1
       else
        k=2
       endif

       !if only vacancies then skip ints
       if (tdam ==3 .and. k==2) cycle

       if(j.lt.50) then
        pj(k,i) = 1.0D12*exp(-emig(k,j)/kb/temp)
       else
        pj(k,i) = 1.0D12*exp(-emig(k,50)/kb/temp)
       endif

       speed = speed + pj(k,i)
      enddo

      dt = 1.0/speed
      scalef = dt*float(jmapsz)

      CALL RANDOM_NUMBER(RND)

      prop = rnd*float(jmapsz)

      sump = 0.0

      do i = 1,ndtot
       j = int(abs(px(6,i)))
       if(px(6,i).lt.0) then
        k=1
       else
        k=2
       endif

       if (tdam ==3 .and. k==2) cycle

       sump = sump + pj(k,i)*scalef
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

!PLACE SELECTED DEFECT IN THE FIRST POSITION

      PX(:,1) = PTEMP(:,1)

      IF (MOD(ISTEP,NFREQ) .EQ. 0) THEN
        WRITE(*,'(A,I12,5(F14.5))') 'DFSEL',1,PX(1:4,1),PX(6,1)
      ENDIF

      deallocate(pj)

      RETURN

      END SUBROUTINE


!***********************************************************************!
!THIS SUBROUTINE IS TO CHOSE AV BASED ON THE MIRGRATION ENERGY BARRIERS
!THE CONVENTIONAL KINETIC MONTE CARLO ALGORITHM IS USED

      SUBROUTINE KMC

      USE SYSTEM
      USE KMCVAR
      USE DEFECT
      USE MPI

      IMPLICIT NONE

      INTEGER I,K,jmapsz

      REAL*8 RND,dt,scalef,prop

      CHARACTER(LEN=40) STR

      IF (RANK .EQ. 0) THEN

!FQSUM:SUM OF FREQUENCY,POSUM:SUM OF PROBABILITIES(FOR CHECK ONLY)

        FQSUM = 0.0D0
        POSUM = 0.0D0
        jmapsz = 1000

!CALCULATE FREQUENCY BASED ON THE MIGRATION ENERGY BARRIER
!THE PREFACTOR HERE IS USED TO BE CONSTANT (NOT TRUE IN GENERAL)
!THE PREFACTOR CALCUATED ON THE FLY IS IMPLEMENTED IN ANOTHER VERSION
! Now considering good dimer search - Gao 2021

        DO K=1,NDIMER
            if(EM(k).lt.10) then
            FQ(K) = 1.0D12*EXP(-EM(K)/KB/TEMP)
            FQSUM = FQSUM + FQ(K)
            endif
        ENDDO

        CALL RANDOM_NUMBER(RND)
        prop = rnd*float(jmapsz)

        !DETERMINE TIME STEP IN KMC

        if(FQSUM.gt.0.0) then
        TMEAN = -LOG10(RND)/FQSUM

        dt = 1.0/fqsum
        scalef = dt*dfloat(jmapsz)
        endif

        DO K=1,NDIMER
          if(EM(k).lt.10.0) then
          POSUM = POSUM + FQ(K)*scalef
          endif
          if(POSUM.ge.prop) exit
        ENDDO

        !EID:THE ID OF THE EVENT THAT IS CHOSEN

        EID = K

        !UPDATING THE PRE-SELECTION ENERGIES

        IF(PX(6,1).LT.0) THEN
          K=1
        ELSE
          K=2
        ENDIF


        !CHECKING FOR TOO HIGH TRANSITION ENERGY

        IF (EM(EID) .LE. EMSKIP .AND. EID .LE. NDIMER) THEN

          TIME = TIME + TMEAN

          !ONLY POSITIVE SADDLE POINT ENERGIES WILL BE USED IN AVERAGE

          IF (EM(EID) .GT. 0) THEN

            I = ABS(INT(PX(6,1)))

            EMIG(K,I) = (EMIG(K,I)*EMIC(K,I) + EM(EID)) / (EMIC(K,I)+1)

            EMIC(K,I) = EMIC(K,I) + 1

          END IF

        ELSE

          SKIP = .TRUE.

          WRITE(STR,'(F8.5)') EMSKIP

          WRITE(6,*) 'MIGRATION ENERGY GREATER THAN '//TRIM(STR)//
     #                'eV -> SKIPPING'

        ENDIF

        IF (MOD(ISTEP,NFREQ) .EQ. 0) THEN
        WRITE(6,'(A,I10,A,F10.4)') 'EID: ',EID,'  RND: ',RND
        ENDIF

      END IF

      CALL MPI_BCAST(TIME, 1, MPI_REAL8, 0, MPI_COMM_WORLD, IERR)

      CALL MPI_BCAST(DTIME, 1, MPI_REAL8, 0, MPI_COMM_WORLD, IERR)

      CALL MPI_BCAST(SKIP, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERR)

      RETURN

      END SUBROUTINE
!***********************************************************************!
!WILL INTRODUCE MORE DEFECTS AFTER A SPECIFIED AMOUNT OF TIME THAT MATCHES
!EXPERIMENTAL DOSE RATE
      SUBROUTINE DOSE

      USE SYSTEM
      USE LMP_DRIVER

      IMPLICIT NONE

      IF(TIME - TIME_LAST .GE. DTIME .OR. ISTEP .EQ. 1) THEN

        IF (RANK .EQ. 0) THEN

          WRITE(*,'(3A14)') 'TIME',' TIME_LAST',' DTIME'

          WRITE(*,'(2E14.5,E14.5)') TIME, TIME_LAST, DTIME

        END IF

        CALL LMP_STRESS_RELAX

        IF (TDAM .EQ. 1) THEN ! Cascade

          CALL SEND_OPLD_POS

          CALL LMP_CASCADE

        ELSE IF (TDAM .EQ. 2 .OR. TDAM .EQ. 3) THEN ! Electron or Vacancy

          CALL BASIC_DAMAGE(NFREKL)

          CALL SEND_OPLD_POS

        END IF

        CALL LMP_RELAX

        CALL RECV_LMP_POS

        TIME_LAST = TIME

        IF (RANK==0) CALC = .TRUE.

      ENDIF

      END SUBROUTINE
