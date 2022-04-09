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

!!THIS SUBROUTINE IS AN INTERFACE, AS REQIRED BY ZXCGR
!!THE POSITIONS AND FORCES ARE PASSED ACCORDINGLY

      SUBROUTINE EAM_FORCE_CG(CGN,CGX,CG_EN,GR)
      USE SYSTEM
      IMPLICIT NONE

      INTEGER I,NTMP
      INTEGER CGN
      REAL*8,DIMENSION(CGN) :: CGX,GR
      REAL*8 CG_EN 
      
      NTMP = CGN/3

      DO I = 1, NTMP
         RA(1,I) = CGX(I)
         RA(2,I) = CGX(I+NTMP)
         RA(3,I) = CGX(I+NTMP*2)
      ENDDO

      CALL EAM_FORCE

      DO I = 1,NTMP
        GR(I)        = -FF(1,I)
        GR(I+NTMP)   = -FF(2,I)
        GR(I+2*NTMP) = -FF(3,I)
      ENDDO

      CG_EN = ENCALC

      RETURN
      END SUBROUTINE EAM_FORCE_CG
!***********************************************************************!
      SUBROUTINE EAM_FORCE

      USE SYSTEM

      IMPLICIT NONE

!!MORE OPTIONS ARE AVALALABLE,INCLUDING ACKLAND04,AUKLAND97,FS84
!!THE POTENTIAL TABLE IS GIVEN IN THIS OPTIMIZED VERSION

      IF(TPOT.EQ.1) THEN

        CALL EAM_FORCETABLE

      ENDIF

      RETURN
      END SUBROUTINE EAM_FORCE
!***********************************************************************!
!!CONSTRUCTING LINKCELL,BASED ON THE EXISTING ALGORITHM

      SUBROUTINE LINKCELL

      USE SYSTEM
      USE NBLVAR
      IMPLICIT NONE

      INTEGER I,IT
      REAL*8 RCL

      HH(:) = 0

      DO I=1,NATOMC

        IF(RA(1,I).LT.0.0) RA(1,I)=RA(1,I)+RBOX
        IF(RA(2,I).LT.0.0) RA(2,I)=RA(2,I)+RBOX
        IF(RA(3,I).LT.0.0) RA(3,I)=RA(3,I)+RBOX

        IF(RA(1,I).GE.RBOX) RA(1,I)=RA(1,I)-RBOX
        IF(RA(2,I).GE.RBOX) RA(2,I)=RA(2,I)-RBOX
        IF(RA(3,I).GE.RBOX) RA(3,I)=RA(3,I)-RBOX

        IT=1+INT(RA(1,I)/CSIZE)
     #      +INT(RA(2,I)/CSIZE)*MSIZE
     #      +INT(RA(3,I)/CSIZE)*MSISQ

        LL(I)=HH(IT)

        HH(IT) = I

      ENDDO

      RETURN
      END SUBROUTINE LINKCELL
!***********************************************************************!
!!MAP THE LINK CLEL

      SUBROUTINE MAPCELL

      USE SYSTEM
      USE NBLVAR

      IMPLICIT NONE

      INTEGER I,J,K,M,IMAP,IT
      INTEGER ITOP,IBOT,JTOP,JBOT,KTOP,KBOT

      M = RBOX/CSIZE

      DO I=1,M
        DO J=1,M
          DO K=1,M
     
            IC(I,J,K)=1+MOD(I-1+M,M)+MOD(J-1+M,M)*M+MOD(K-1+M,M)*M*M

            IT = 1+MOD(I-1+M,M)+MOD(J-1+M,M)*M+MOD(K-1+M,M)*M*M

            CP(1,IT) = I*CSIZE
            CP(2,IT) = J*CSIZE
            CP(3,IT) = K*CSIZE

          ENDDO
        ENDDO
      ENDDO

      DO I=1,M
        DO J=1,M
          DO K=1,M

            IMAP=(IC(I,J,K)-1)*26

            ITOP = I+1
            IBOT = I-1
            JTOP = J+1
            JBOT = J-1
            KTOP = K+1
            KBOT = K-1

            IF((I+1).GT.M) ITOP = I+1-M
            IF((I-1).LT.1) IBOT = I-1+M
            IF((J+1).GT.M) JTOP = J+1-M
            IF((J-1).LT.1) JBOT = J-1+M
            IF((K+1).GT.M) KTOP = K+1-M
            IF((K-1).LT.1) KBOT = K-1+M

            MP(IMAP+1) =IC(ITOP,J   ,K   )
            MP(IMAP+2) =IC(IBOT,J   ,K   )
            MP(IMAP+3) =IC(I   ,JTOP,K   )
            MP(IMAP+4) =IC(ITOP,JTOP,K   )
            MP(IMAP+5) =IC(IBOT,JTOP,K   )
            MP(IMAP+6) =IC(I   ,JBOT,K   )
            MP(IMAP+7) =IC(ITOP,JBOT,K   )
            MP(IMAP+8) =IC(IBOT,JBOT,K   )
            MP(IMAP+9) =IC(I   ,J   ,KTOP)
            MP(IMAP+10)=IC(ITOP,J   ,KTOP)
            MP(IMAP+11)=IC(IBOT,J   ,KTOP)
            MP(IMAP+12)=IC(I   ,JTOP,KTOP)
            MP(IMAP+13)=IC(ITOP,JTOP,KTOP)
            MP(IMAP+14)=IC(IBOT,JTOP,KTOP)
            MP(IMAP+15)=IC(I   ,JBOT,KTOP)
            MP(IMAP+16)=IC(ITOP,JBOT,KTOP)
            MP(IMAP+17)=IC(IBOT,JBOT,KTOP)
            MP(IMAP+18)=IC(I   ,J   ,KBOT)
            MP(IMAP+19)=IC(ITOP,J   ,KBOT)
            MP(IMAP+20)=IC(IBOT,J   ,KBOT)
            MP(IMAP+21)=IC(I   ,JTOP,KBOT)
            MP(IMAP+22)=IC(ITOP,JTOP,KBOT)
            MP(IMAP+23)=IC(IBOT,JTOP,KBOT)
            MP(IMAP+24)=IC(I  , JBOT,KBOT)
            MP(IMAP+25)=IC(ITOP,JBOT,KBOT)
            MP(IMAP+26)=IC(IBOT,JBOT,KBOT)

          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE MAPCELL
!***********************************************************************!
!!USE LINKCELL TO BUILD THE NEIGHBOR LIST

      SUBROUTINE VERLIST

      USE SYSTEM
      USE NBLVAR

      IMPLICIT NONE

      INTEGER I,J,K,NB,JN
      INTEGER IT,JC,NA,JCELL

      REAL*8 RX,RY,RZ,RS,RT,PBC
      REAL*8 TX,TY,TZ
 
      NI = 0
      NL = 0

      DO IT=1,NCELL

        I=HH(IT)

        DO WHILE(I.GT.0)

          TX = RA(1,I)

          TY = RA(2,I)

          TZ = RA(3,I)

          J=LL(I)

          DO WHILE(J.GT.0)

            RX=PBC(RA(1,J)-TX)

            RY=PBC(RA(2,J)-TY)

            RZ=PBC(RA(3,J)-TZ)

            RS=RX*RX+RY*RY+RZ*RZ

            IF(RS.LE.SLIST) THEN

              NI(I) = NI(I) +1

              NL(NI(I),I) = J

              NI(J) = NI(J) + 1

              NL(NI(J),J) = I

            ENDIF

            J=LL(J)

          ENDDO

          JCELL=26*(IT-1)

          DO NA=1,26

            JC=MP(JCELL+NA)

            K=HH(JC)

            DO WHILE(K.GT.0)

              RX=PBC(RA(1,K)-TX)

              RY=PBC(RA(2,K)-TY)

              RZ=PBC(RA(3,K)-TZ)

              RS=RX*RX+RY*RY+RZ*RZ

              IF(RS.LE.SLIST) THEN

                NI(I) = NI(I) +1

                NL(NI(I),I) = K

              ENDIF

              K=LL(K)

            ENDDO

          ENDDO

          I=LL(I)

        ENDDO

      ENDDO

      RETURN

      END SUBROUTINE VERLIST
!***********************************************************************!
!!CHECK WHETHER NEIGHBOR LIST NEEDS TO BE UPDATED

      SUBROUTINE CHECKNBLIST

      USE SYSTEM
      USE NBLVAR
      IMPLICIT NONE

      INTEGER I
      REAL*8 DISPMX,DISPMY,DISPMZ

      DISPMX = 0.0D0
      DISPMY = 0.0D0
      DISPMZ = 0.0D0

!!FINDIND THE MAXIMUM DISPLACMENT

!!$OMP PARALLEL DO REDUCTION(MAX:DISPMX,DISMPMY,DISPMZ)

      DO I=1,NATOMC

        DISPMX = MAX(DABS(RA(1,I)-RL(1,I)),DISPMX)
        DISPMY = MAX(DABS(RA(2,I)-RL(2,I)),DISPMY)
        DISPMZ = MAX(DABS(RA(3,I)-RL(3,I)),DISPMZ)

      ENDDO

!!$OMP END PARALLEL DO

      DISPMX = DSQRT(DISPMX*DISPMX+DISPMY*DISPMY+DISPMZ*DISPMZ)

!!IF THE MAXIMUM DISPLACEMENT IS LARGER THEN THE SKIN THICKNESS
!!THE NEIGHBOR LIST NEEDS TO BE UPDATED

      NL_NEW = (DISPMX.GT.RSKIN)

      RETURN

      END SUBROUTINE CHECKNBLIST
!***********************************************************************!
!!CALCULATE THE ENERGY AND FORCE BASED ON A POTENTIAL TABLE

      SUBROUTINE EAM_FORCETABLE

      USE SYSTEM
      USE NBLVAR
      USE POTTAB

      IMPLICIT NONE

      INTEGER I,J,K,NTEMP,JN

      INTEGER PN,DN,EN,TN

      REAL*8 RX,RY,RZ,RS,TX,TY,TZ,PBC

      REAL*8 TMP,TME,TFP,TFD,TFE1,TFE2,TFF

      REAL*8 RP,RD,RE,TE


      FF = 0.0

      EE = 0.0

      ENCALC = 0.0

      CALL CHECKNBLIST

      IF(NL_NEW) THEN

        CALL LINKCELL

        CALL VERLIST

        RL = RA

      ENDIF

!$OMP PARALLEL DO 
!$OMP& PRIVATE(NTEMP,TX,TY,TZ,JN,J,RX,RY,RZ,RS,PN,RP,TMP,DN,RD,TME)

      DO I=1,NATOMC

        NTEMP = NI(I)

        TX = RA(1,I)

        TY = RA(2,I)

        TZ = RA(3,I)

        DO JN = 1,NTEMP

          J =NL(JN,I)

          RX=PBC(RA(1,J)-TX)

          RY=PBC(RA(2,J)-TY)

          RZ=PBC(RA(3,J)-TZ)

          RS=RX*RX+RY*RY+RZ*RZ

!!CALCULATE THE PAIR TERM

          IF(RS.LE.PCUTS) THEN

            PN = (RS-RPMIN)*PSTEI + 1

            RP = RS-PR(PN)

            TMP=PA(PN)+RP*(PB(PN)+RP*(PC(PN)+RP*PD(PN)))

            EE(1,I)=EE(1,I) + TMP

!!DENSITY TERM

            IF(RS.LE.DCUTS) THEN

            DN=(RS-RDMIN)*DSTEI + 1

            RD = RS-DR(DN)

            TME=DA(DN)+RD*(DB(DN)+RD*(DC(DN)+RD*DD(DN)))

            EE(2,I)=EE(2,I) + TME

            ENDIF 

          ENDIF

        ENDDO

      ENDDO

!$OMP END PARALLEL DO

!!CALCULATE FORCE

!$OMP PARALLEL DO REDUCTION(+:ENCALC) PRIVATE(EN,RE,TFE1,NTEMP,TX,TY,TZ,
!$OMP& JN,J,RX,RY,RZ,RS,PN,RP,TFP,DN,RD,TFD,TN,TE,TFE2,TFF)

      DO I=1,NATOMC

        EN = EE(2,I)*ESTEI + 1

        RE = EE(2,I)-(EN-1)*ESTEP

        EE(3,I)=EE(1,I)*0.5D0+EA(EN)+RE*(EB(EN)+RE*(EC(EN)+RE*ED(EN)))

        ENCALC=ENCALC+EE(3,I)

        TFE1=FEB(EN)+RE*(FEC(EN)+RE*FED(EN))

          NTEMP = NI(I)

          TX = RA(1,I)

          TY = RA(2,I)

          TZ = RA(3,I)

          DO JN=1,NTEMP

            J = NL(JN,I)

            RX=PBC(RA(1,J)-TX)

            RY=PBC(RA(2,J)-TY)

            RZ=PBC(RA(3,J)-TZ)

            RS=RX*RX+RY*RY+RZ*RZ

            IF(RS.LE.PCUTS) THEN

              PN=(RS-RPMIN)*PSTEI + 1

              RP=RS-PR(PN)

              TFP=FPB(PN)+RP*(FPC(PN)+RP*FPD(PN))

              IF(RS.LE.DCUTS) THEN

              DN=(RS-RDMIN)*DSTEI + 1

              RD=RS-DR(DN)

              TFD=FDB(DN)+RD*(FDC(DN)+RD*FDD(DN))

              ELSE

              TFD=0.0D0

              ENDIF
              
              TN = EE(2,J)*ESTEI + 1

              TE = EE(2,J)-(TN-1)*ESTEP

              TFE2 = FEB(TN)+TE*(FEC(TN)+TE*FED(TN))

              TFF = TFP + TFD*(TFE1+TFE2)

              FF(1,I)=FF(1,I) + TFF*RX

              FF(2,I)=FF(2,I) + TFF*RY

              FF(3,I)=FF(3,I) + TFF*RZ

            ENDIF

          ENDDO

      ENDDO

!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE EAM_FORCETABLE
!***********************************************************************!

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

