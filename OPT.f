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

!!PERFORM STATIC RELAXATION OF THE SYSTEM

      SUBROUTINE RELAX

      USE SYSTEM

      IMPLICIT NONE

      INTEGER I,J

!      INCLUDE 'mkl_rci.fi'

!      INTEGER I,N,RCI_REQUEST,NUMITR,IPAR(128)

!      REAL*8 :: DPAR(128)

!      REAL*8,ALLOCATABLE ::  X(:),B(:),TMP(:,:)

!!IF CONJUGATE GRADIENT METHOD IS CHOSEN

      IF(TOPT.EQ.1) THEN

        CALL OPT_CG_INIT(NATOMC)

        CALL OPT_CG

        CALL OPT_CG_FINISH

!!A POTENTIAL WAY TO CALL INTEL LIBRARIES

!       CALL DCG_INIT(N,X,B,RCI_REQUEST,IPAR,DPAR,TMP)

!       CALL DCG_CHECK(N,X,B,RCI_REQUEST,IPAR,DPAR,TMP)

!       CALL DCG(NATOMC*3,X,B,RCI_REQUEST,IPAR,DPAR,TMP)

!       CALL DCG_GET(NATOMC*3,X,B,RCI_REQUEST,IPAR,DPAR,TMP,NUMITR)

!!IF LBFGS IS CHOSEN

      ELSE IF(TOPT.EQ.2) THEN

        CALL OPT_LBFGS_INIT

        CALL OPT_LBFGS

        CALL OPT_LBFGS_FINISH

      ENDIF

!!CALCULATE ENERGY AND FORCE AFTER RELAXATION

      CALL EAM_FORCE

      WRITE(*,*) 'ENERGY AFTER RELAXATION:',ENCALC

!!OUTPUT DEFECT STRUCTURE

      CALL VISUAL

!!OUTPUT OPTIMIZED STRUCTURE INTO STD_KOPT.OUT

      WRITE(26,*) NATOMC

      DO I = 1, NATOMC
        WRITE(26,'(3F14.5)') (RA(J,I), J=1,3)
      ENDDO

      RETURN

      END SUBROUTINE RELAX
!***********************************************************************!
      SUBROUTINE OPT_CG

      USE SYSTEM

      USE OPTCGVAR

      IMPLICIT NONE

      INTEGER I

      EXTERNAL EAM_FORCE_CG

!!! COULD USE OPNMP
!! PASSING THE POSITION OF THE ATOMS (RA) INTO XX,USED BY ZXCGR

      DO I = 1, NATOMC

        XX(I)          = RA(1,I)

        XX(I+NATOMC)   = RA(2,I)

        XX(I+2*NATOMC) = RA(3,I)

      ENDDO

!!CALL ZXCGR, USING EXISTING CG ALGORITHM,SEE ZXCGR FOR DETAILS

      CALL ZXCGR (EAM_FORCE_CG,(3*NATOMC), GRADTL, MAXFN, DFPRED,
     #       XX, GG, EN, WW, IER)

!!PASSING THE POSTION OF ATOMS BACK TO RA

      DO I = 1, NATOMC

        RA(1,I) = XX(I)

        RA(2,I) = XX(I+NATOMC)

        RA(3,I) = XX(I+NATOMC*2)

      ENDDO

      RETURN

      END SUBROUTINE OPT_CG
!***********************************************************************!
      SUBROUTINE OPT_CG_INIT(NATOMS)

      USE OPTCGVAR

      IMPLICIT NONE

      INTEGER NATOMS

      DFPRED    = 0.01

      GRADTL    = 1.0D-3

      MAXFN     = 500

!!DIMENSION DEFINED AS REQUIRED

      ALLOCATE(XX(3*NATOMS))

      ALLOCATE(GG(3*NATOMS))

      ALLOCATE(WW(6*(3*NATOMS)))

      RETURN

      END SUBROUTINE OPT_CG_INIT
!***********************************************************************!
      SUBROUTINE OPT_CG_FINISH()

      USE OPTCGVAR

      IMPLICIT NONE

!!FREE MEMORY

      DEALLOCATE(XX)

      DEALLOCATE(GG)

      DEALLOCATE(WW)

      RETURN

      END SUBROUTINE OPT_CG_FINISH
!***********************************************************************!
      SUBROUTINE OPT_LBFGS_INIT

      USE SYSTEM

      USE OPTVAR

      IMPLICIT NONE

      INTEGER I

      N = 3*NATOMC

      M = 3

      IF(TJOB.EQ.2) THEN

        IPRINT(1) = -1

      ELSE

        IPRINT(1) =  1

      ENDIF

      IPRINT(2) = 0

      EPS = 1E-5

      XTOL = 1E-8

      DIAGCO = .FALSE.

      IIFLAG = 0

      ALLOCATE(X(N))

      ALLOCATE(G(N))

      ALLOCATE(DIAG(N))

      ALLOCATE(W(N*(2*M+1)+2*M))

      DO I = 1, NATOMC

        X(I) = RA(1,I)

        X(I+NATOMC) = RA(2,I)

        X(I+NATOMC*2) = RA(3,I)

      ENDDO

      RETURN

      END SUBROUTINE OPT_LBFGS_INIT
!***********************************************************************!
!!LBFGS ALGORITHM,SEE BELOW FOR DETAILS

      SUBROUTINE OPT_LBFGS
      USE SYSTEM
      USE OPTVAR
      IMPLICIT NONE

      INTEGER I,NITER

      DO NITER=1,500

        DO I=1,NATOMC

!!! COULD USE OPENMP

          RA(1,I) = X(I)

          RA(2,I) = X(I+NATOMC)

          RA(3,I) = X(I+NATOMC*2)

        ENDDO

        CALL EAM_FORCE

        DO I=1,NATOMC

          G(I)          = -FF(1,I)

          G(I+NATOMC)   = -FF(2,I)

          G(I+NATOMC*2) = -FF(3,I)

        ENDDO

        F = ENCALC

        CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IIFLAG)

        IF(IIFLAG.LE.0) EXIT

      END DO

      DO I = 1, NATOMC

        RA(1,I) = X(I)

        RA(2,I) = X(I+NATOMC)

        RA(3,I) = X(I+NATOMC*2)

      ENDDO

      RETURN

      END SUBROUTINE OPT_LBFGS
!***********************************************************************!
      SUBROUTINE OPT_LBFGS_FINISH

      USE OPTVAR

      IMPLICIT NONE

      DEALLOCATE(X)

      DEALLOCATE(G)

      DEALLOCATE(DIAG)

      DEALLOCATE(W)

      RETURN

      END SUBROUTINE OPT_LBFGS_FINISH
!***********************************************************************!



!***********************************************************************!
!***********************************************************************!
!***********************************************************************!
C-----------------------------------------------------------------------ZXCG0030
C                                                                       ZXCG0040
C   COMPUTER            - HP9000/DOUBLE                                 ZXCG0050
C                                                                       ZXCG0060
C   LATEST REVISION     - NOVEMBER 1, 1979                              ZXCG0070
C                                                                       ZXCG0080
C   PURPOSE             - A CONJUGATE GRADIENT ALGORITHM FOR FINDING    ZXCG0090
C                           THE MINIMUM OF A FUNCTION OF N VARIABLES    ZXCG0100
C                                                                       ZXCG0110
C   USAGE               - CALL ZXCGR (FUNCT,N,ACC,MAXFN,DFPRED,X,G,F,W, ZXCG0120
C                           IER)                                        ZXCG0130
C                                                                       ZXCG0140
C   ARGUMENTS    FUNCT  - A USER SUPPLIED SUBROUTINE WHICH CALCULATES   ZXCG0150
C                           THE OBJECTIVE FUNCTION AND ITS GRADIENT     ZXCG0160
C                           FOR GIVEN PARAMETER VALUES                  ZXCG0170
C                           X(1),X(2),...,X(N).                         ZXCG0180
C                           THE CALLING SEQUENCE HAS THE FOLLOWING FORM ZXCG0190
C                           CALL FUNCT (N,X,F,G)                        ZXCG0200
C                           WHERE X AND G ARE VECTORS OF LENGTH N.      ZXCG0210
C                           THE SCALAR F IS FOR THE OBJECTIVE FUNCTION. ZXCG0220
C                           G(1), G(2), ..., G(N) ARE FOR THE COMPONENTSZXCG0230
C                           OF THE GRADIENT OF F.                       ZXCG0240
C                           FUNCT MUST APPEAR IN AN EXTERNAL STATEMENT  ZXCG0250
C                           IN THE CALLING PROGRAM. FUNCT MUST NOT      ZXCG0260
C                           ALTER THE VALUES OF X(I),I=1,...,N OR N.    ZXCG0270
C                N      - THE NUMBER OF PARAMETERS OF THE OBJECTIVE     ZXCG0280
C                           FUNCTION. (INPUT) (I.E.,THE LENGTH OF X)    ZXCG0290
C                ACC    - CONVERGENCE CRITERION. (INPUT)                ZXCG0300
C                           THE CALCULATION ENDS WHEN THE SUM OF SQUARESZXCG0310
C                           OF THE COMPONENTS OF G IS LESS THAN ACC.    ZXCG0320
C                MAXFN  - MAXIMUM NUMBER OF FUNCTION EVALUATIONS (I.E., ZXCG0330
C                           CALLS TO SUBROUTINE FUNCT) ALLOWED. (INPUT) ZXCG0340
C                           IF MAXFN IS SET TO ZERO, THEN THERE IS      ZXCG0350
C                           NO RESTRICTION ON THE NUMBER OF FUNCTION    ZXCG0360
C                           EVALUATIONS.                                ZXCG0370
C                DFPRED - A ROUGH ESTIMATE OF THE EXPECTED REDUCTION    ZXCG0380
C                           IN F, WHICH IS USED TO DETERMINE THE SIZE   ZXCG0390
C                           OF THE INITIAL CHANGE TO X. (INPUT)         ZXCG0400
C                           NOTE THAT DFPRED IS THE EXPECTED REDUCTION  ZXCG0410
C                           ITSELF, AND DOES NOT DEPEND ON ANY RATIOS.  ZXCG0420
C                           A BAD VALUE OF DFPRED CAUSES AN ERROR       ZXCG0430
C                           MESSAGE, WITH IER=129, AND A RETURN ON THE  ZXCG0440
C                           FIRST ITERATION. (SEE THE DESCRIPTION OF    ZXCG0450
C                           IER BELOW)                                  ZXCG0460
C                X      - VECTOR OF LENGTH N CONTAINING PARAMETER       ZXCG0470
C                           VALUES.                                     ZXCG0480
C                         ON INPUT, X MUST CONTAIN THE INITIAL          ZXCG0490
C                           PARAMETER ESTIMATES.                        ZXCG0500
C                         ON OUTPUT, X CONTAINS THE FINAL PARAMETER     ZXCG0510
C                           ESTIMATES AS DETERMINED BY ZXCGR.           ZXCG0520
C                G      - A VECTOR OF LENGTH N CONTAINING THE           ZXCG0530
C                           COMPONENTS OF THE GRADIENT OF F AT THE      ZXCG0540
C                           FINAL PARAMETER ESTIMATES. (OUTPUT)         ZXCG0550
C                F      - A SCALAR CONTAINING THE VALUE OF THE FUNCTION ZXCG0560
C                           AT THE FINAL PARAMETER ESTIMATES. (OUTPUT)  ZXCG0570
C                W      - WORK VECTOR OF LENGTH 6*N.                    ZXCG0580
C                IER    - ERROR PARAMETER. (OUTPUT)                     ZXCG0590
C                         IER = 0 IMPLIES THAT CONVERGENCE WAS          ZXCG0600
C                           ACHIEVED AND NO ERRORS OCCURRED.            ZXCG0610
C                         TERMINAL ERROR                                ZXCG0620
C                           IER = 129 IMPLIES THAT THE LINE SEARCH OF   ZXCG0630
C                             AN INTEGRATION WAS ABANDONED. THIS        ZXCG0640
C                             ERROR MAY BE CAUSED BY AN ERROR IN THE    ZXCG0650
C                             GRADIENT.                                 ZXCG0660
C                           IER = 130 IMPLIES THAT THE CALCULATION      ZXCG0670
C                             CANNOT CONTINUE BECAUSE THE SEARCH        ZXCG0680
C                             DIRECTION IS UPHILL.                      ZXCG0690
C                           IER = 131 IMPLIES THAT THE ITERATION WAS    ZXCG0700
C                             TERMINATED BECAUSE MAXFN WAS EXCEEDED.    ZXCG0710
C                           IER = 132 IMPLIES THAT THE CALCULATION      ZXCG0720
C                             WAS TERMINATED BECAUSE TWO CONSECUTIVE    ZXCG0730
C                             ITERATIONS FAILED TO REDUCE F.            ZXCG0740
C                                                                       ZXCG0750
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         ZXCG0760
C                       - SINGLE/H36,H48,H60                            ZXCG0770
C                                                                       ZXCG0780
C   REQD. IMSL ROUTINES - UERTST,UGETIO                                 ZXCG0790
C   (THE CALL IS COMMENTED OUT)
C                                                                       ZXCG0800
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           ZXCG0810
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      ZXCG0820
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  ZXCG0830
C                                                                       ZXCG0840
C   REMARKS  1.  THE ROUTINE INCLUDES NO THOROUGH CHECKS ON THE PART    ZXCG0850
C                OF THE USER PROGRAM THAT CALCULATES THE DERIVATIVES    ZXCG0860
C                OF THE OBJECTIVE FUNCTION. THEREFORE, BECAUSE          ZXCG0870
C                DERIVATIVE CALCULATION IS A FREQUENT SOURCE OF         ZXCG0880
C                ERROR, THE USER SHOULD VERIFY INDEPENDENTLY THE        ZXCG0890
C                CORRECTNESS OF THE DERIVATIVES THAT ARE GIVEN TO       ZXCG0900
C                THE ROUTINE.                                           ZXCG0910
C            2.  BECAUSE OF THE CLOSE RELATION BETWEEN THE CONJUGATE    ZXCG0920
C                GRADIENT METHOD AND THE METHOD OF STEEPEST DESCENTS,   ZXCG0930
C                IT IS VERY HELPFUL TO CHOOSE THE SCALE OF THE          ZXCG0940
C                VARIABLES IN A WAY THAT BALANCES THE MAGNITUDES OF     ZXCG0950
C                THE COMPONENTS OF A TYPICAL DERIVATE VECTOR. IT        ZXCG0960
C                CAN BE PARTICULARLY INEFFICIENT IF A FEW COMPONENTS    ZXCG0970
C                OF THE GRADIENT ARE MUCH LARGER THAN THE REST.         ZXCG0980
C            3.  IF THE VALUE OF THE PARAMETER ACC IN THE ARGUMENT      ZXCG0990
C                LIST OF THE ROUTINE IS SET TO ZERO, THEN THE           ZXCG1000
C                SUBROUTINE WILL CONTINUE ITS CALCULATION UNTIL IT      ZXCG1010
C                STOPS REDUCING THE OBJECTIVE FUNCTION. IN THIS CASE    ZXCG1020
C                THE USUAL BEHAVIOUR IS THAT CHANGES IN THE             ZXCG1030
C                OBJECTIVE FUNCTION BECOME DOMINATED BY COMPUTER        ZXCG1040
C                ROUNDING ERRORS BEFORE PRECISION IS LOST IN THE        ZXCG1050
C                GRADIENT VECTOR. THEREFORE, BECAUSE THE POINT OF       ZXCG1060
C                VIEW HAS BEEN TAKEN THAT THE USER REQUIRES THE         ZXCG1070
C                LEAST POSSIBLE VALUE OF THE FUNCTION, A VALUE OF       ZXCG1080
C                THE OBJECTIVE FUNCTION THAT IS SMALL DUE TO            ZXCG1090
C                COMPUTER ROUNDING ERRORS CAN PREVENT FURTHER           ZXCG1100
C                PROGRESS. HENCE THE PRECISION IN THE FINAL VALUES      ZXCG1110
C                OF THE VARIABLES MAY BE ONLY ABOUT HALF THE            ZXCG1120
C                NUMBER OF SIGNIFICANT DIGITS IN THE COMPUTER           ZXCG1130
C                ARITHMETIC, BUT THE LEAST VALUE OF F IS USUALLY        ZXCG1140
C                FOUND TO QUITE HIGH ACCURACY.                          ZXCG1150
C                                                                       ZXCG1160
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       ZXCG1170
C                                                                       ZXCG1180
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN ZXCG1190
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    ZXCG1200
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        ZXCG1210
C                                                                       ZXCG1220
C-----------------------------------------------------------------------ZXCG1230
C                                                                       ZXCG1240
      SUBROUTINE ZXCGR  (FUNCT,N,ACC,MAXFN,DFPRED,X,G,F,W,IER)          ZXCG1250
C                                  SPECIFICATIONS FOR ARGUMENTS         ZXCG1260
      INTEGER            N,MAXFN,IER                                    ZXCG1270
      DOUBLE PRECISION   ACC,DFPRED,X(N),G(N),F,W(6*N)                    ZXCG1280
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   ZXCG1290
      INTEGER            MAXLIN,MXFCON,I,IGINIT,IGOPT,IRETRY,IRSDG,     ZXCG1300
     1                   IRSDX,ITERC,ITERFM,ITERRS,IXOPT,NCALLS,NFBEG,  ZXCG1310
     2                   NFOPT                                          ZXCG1320
      DOUBLE PRECISION   BETA,DDSPLN,DFPR,FCH,FINIT,FMIN,GAMDEN,GAMA,   ZXCG1330
     1                   GINIT,GMIN,GNEW,GSPLN,GSQRD,SBOUND,STEP,STEPCH,ZXCG1340
     2                   STMIN,SUM,WORK                                 ZXCG1350
      DATA               MAXLIN/5/,MXFCON/2/                            ZXCG1360
C                                  FIRST EXECUTABLE STATEMENT           ZXCG1370
        EXTERNAL FUNCT
      IER = 0                                                           ZXCG1380
C                                  THE WORKING SPACE ARRAY IS SPLIT     ZXCG1390
C                                    INTO SIX VECTORS OF LENGTH N. THE  ZXCG1400
C                                    FIRST PART IS USED FOR THE SEARCH  ZXCG1410
C                                    DIRECTION OF AN ITERATION. THE     ZXCG1420
C                                    SECOND AND THIRD PARTS CONTAIN THE ZXCG1430
C                                    INFORMATION THAT IS REQUIRED BY    ZXCG1440
C                                    THE CONJUGACY CONDITIONS OF THE    ZXCG1450
C                                    RESTART PROCEDURE. THE FOURTH PART ZXCG1460
C                                    CONTAINS THE GRADIENT AT THE START ZXCG1470
C                                    OF AN ITERATION. THE FIFTH PART    ZXCG1480
C                                    CONTAINS THE PARAMETERS THAT GIVE  ZXCG1490
C                                    THE LEAST CALCULATED VALUE OF F.   ZXCG1500
C                                    THE SIXTH PART CONTAINS THE        ZXCG1510
C                                    GRADIENT VECTOR WHERE F IS LEAST.  ZXCG1520
      IRSDX = N                                                         ZXCG1530
      IRSDG = IRSDX+N                                                   ZXCG1540
      IGINIT = IRSDG+N                                                  ZXCG1550
      IXOPT = IGINIT+N                                                  ZXCG1560
      IGOPT = IXOPT+N                                                   ZXCG1570
C                                  SET SOME PARAMETERS TO BEGIN THE     ZXCG1580
C                                    CALCULATION. ITERC AND             ZXCG1590
C                                    NCALLS COUNT THE NUMBER OF         ZXCG1600
C                                    ITERATIONS AND CALLS OF FUNCT.     ZXCG1610
C                                    ITERFM IS THE NUMBER OF THE MOST   ZXCG1620
C                                    RECENT ITERATION THAT DECREASES F. ZXCG1630
      ITERC = 0                                                         ZXCG1640
      NCALLS = 0                                                        ZXCG1650
      ITERFM = ITERC                                                    ZXCG1660
C                                  CALL SUBROUTINE FUNCT. LET THE       ZXCG1670
C                                    INITIAL SEARCH DIRECTION BE MINUS  ZXCG1680
C                                    THE GRADIENT VECTOR. USUALLY THE   ZXCG1690
C                                    PARAMETER ITERRS GIVES THE         ZXCG1700
C                                    ITERATION NUMBER OF THE MOST       ZXCG1710
C                                    RECENT RESTART, BUT IT IS SET TO   ZXCG1720
C                                    ZERO WHEN THE STEEPEST DESCENT     ZXCG1730
C                                    DIRECTION IS USED.                 ZXCG1740
    5 NCALLS = NCALLS+1                                                 ZXCG1750
      CALL FUNCT (N,X,F,G)                                              ZXCG1760
      IF (NCALLS.GE.2) GO TO 20                                         ZXCG1770
   10 DO 15 I=1,N                                                       ZXCG1780
   15 W(I) = -G(I)                                                      ZXCG1790
      ITERRS = 0                                                        ZXCG1800
      IF (ITERC.GT.0) GO TO 80                                          ZXCG1810
C                                  SET SUM TO G SQUARED. GMIN AND GNEW  ZXCG1820
C                                    ARE THE OLD AND THE NEW            ZXCG1830
C                                    DIRECTIONAL DERIVATIVES ALONG THE  ZXCG1840
C                                    CURRENT SEARCH DIRECTION. LET FCH  ZXCG1850
C                                    BE THE DIFFERENCE BETWEEN F AND    ZXCG1860
C                                    THE PREVIOUS BEST VALUE OF THE     ZXCG1870
C                                    OBJECTIVE FUNCTION.                ZXCG1880
   20 GNEW = 0.0D0                                                      ZXCG1890
      SUM = 0.0D0                                                       ZXCG1900
      DO 25 I=1,N                                                       ZXCG1910
         GNEW = GNEW+W(I)*G(I)                                          ZXCG1920
   25 SUM = SUM+G(I)**2                                                 ZXCG1930
      IF (NCALLS.EQ.1) GO TO 35                                         ZXCG1940
      FCH = F-FMIN                                                      ZXCG1950
C                                  STORE THE VALUES OF X, F AND G, IF   ZXCG1960
C                                    THEY ARE THE BEST THAT HAVE BEEN   ZXCG1970
C                                    CALCULATED SO FAR, AND NOTE G      ZXCG1980
C                                    SQUARED AND THE VALUE OF NCALLS.   ZXCG1990
C                                    TEST FOR CONVERGENCE.              ZXCG2000
      IF (FCH) 35,30,50                                                 ZXCG2010
   30 IF (GNEW/GMIN.LT.-1.0D0) GO TO 45                                 ZXCG2020
   35 FMIN = F                                                          ZXCG2030
      GSQRD = SUM                                                       ZXCG2040
      NFOPT = NCALLS                                                    ZXCG2050
      DO 40 I=1,N                                                       ZXCG2060
         W(IXOPT+I) = X(I)                                              ZXCG2070
   40 W(IGOPT+I) = G(I)                                                 ZXCG2080
        WRITE(*,'(A,2I12,2F14.5)')'OPT',ITERC,NCALLS, SUM, F
C        WRITE(51,*)ITERC,NCALLS, SUM, F
   45 IF (SUM.LE.ACC) GO TO 9005                                        ZXCG2090
C                                  TEST IF THE VALUE OF MAXFN ALLOWS    ZXCG2100
C                                    ANOTHER CALL OF FUNCT.             ZXCG2110
   50 IF (NCALLS.NE.MAXFN) GO TO 55                                     ZXCG2120
      IER = 131                                                         ZXCG2130
      GO TO 9000                                                        ZXCG2140
   55 IF (NCALLS.GT.1) GO TO 100                                        ZXCG2150
C                                  SET DFPR TO THE ESTIMATE OF THE      ZXCG2160
C                                    REDUCTION IN F GIVEN IN THE        ZXCG2170
C                                    ARGUMENT LIST, IN ORDER THAT THE   ZXCG2180
C                                    INITIAL CHANGE TO THE PARAMETERS   ZXCG2190
C                                    IS OF A SUITABLE SIZE. THE VALUE   ZXCG2200
C                                    OF STMIN IS USUALLY THE            ZXCG2210
C                                    STEP-LENGTH OF THE MOST RECENT     ZXCG2220
C                                    LINE SEARCH THAT GIVES THE LEAST   ZXCG2230
C                                    CALCULATED VALUE OF F.             ZXCG2240
      DFPR = DFPRED                                                     ZXCG2250
      STMIN = DFPRED/GSQRD                                              ZXCG2260
C                                  BEGIN THE ITERATION                  ZXCG2270
   80 ITERC = ITERC+1                                                   ZXCG2280
C                                  STORE THE INITIAL FUNCTION VALUE AND ZXCG2290
C                                    GRADIENT, CALCULATE THE INITIAL    ZXCG2300
C                                    DIRECTIONAL DERIVATIVE, AND BRANCH ZXCG2310
C                                    IF ITS VALUE IS NOT NEGATIVE. SET  ZXCG2320
C                                    SBOUND TO MINUS ONE TO INDICATE    ZXCG2330
C                                    THAT A BOUND ON THE STEP IS NOT    ZXCG2340
C                                    KNOWN YET, AND SET NFBEG TO THE    ZXCG2350
C                                    CURRENT VALUE OF NCALLS. THE       ZXCG2360
C                                    PARAMETER IRETRY SHOWS THE NUMBER  ZXCG2370
C                                    OF ATTEMPTS AT SATISFYING THE BETA ZXCG2380
C                                    CONDITION.                         ZXCG2390
      FINIT = F                                                         ZXCG2400
      GINIT = 0.0D0                                                     ZXCG2410
      DO 85 I=1,N                                                       ZXCG2420
         W(IGINIT+I) = G(I)                                             ZXCG2430
   85 GINIT = GINIT+W(I)*G(I)                                           ZXCG2440
      IF (GINIT.GE.0.0D0) GO TO 165                                     ZXCG2450
      GMIN = GINIT                                                      ZXCG2460
      SBOUND = -1.0D0                                                   ZXCG2470
      NFBEG = NCALLS                                                    ZXCG2480
      IRETRY = -1                                                       ZXCG2490
C                                  SET STEPCH SO THAT THE INITIAL       ZXCG2500
C                                    STEP-LENGTH IS CONSISTENT WITH THE ZXCG2510
C                                    PREDICTED REDUCTION IN F, SUBJECT  ZXCG2520
C                                    TO THE CONDITION THAT IT DOES NOT  ZXCG2530
C                                    EXCEED THE STEP-LENGTH OF THE      ZXCG2540
C                                    PREVIOUS ITERATION. LET STMIN BE   ZXCG2550
C                                    THE STEP TO THE LEAST CALCULATED   ZXCG2560
C                                    VALUE OF F.                        ZXCG2570
      STEPCH = DMIN1(STMIN,DABS(DFPR/GINIT))                            ZXCG2580
      STMIN = 0.0D0                                                     ZXCG2590
C                                  CALL SUBROUTINE FUNCT AT THE VALUE   ZXCG2600
C                                    OF X THAT IS DEFINED BY THE NEW    ZXCG2610
C                                    CHANGE TO THE STEP-LENGTH, AND LET ZXCG2620
C                                    THE NEW STEP-LENGTH BE STEP. THE   ZXCG2630
C                                    VARIABLE WORK IS USED AS WORK      ZXCG2640
C                                    SPACE.                             ZXCG2650
   90 STEP = STMIN+STEPCH                                               ZXCG2660
      WORK = 0.0D0                                                      ZXCG2670
      DO 95 I=1,N                                                       ZXCG2680
         X(I) = W(IXOPT+I)+STEPCH*W(I)                                  ZXCG2690
   95 WORK = DMAX1(WORK,DABS(X(I)-W(IXOPT+I)))                          ZXCG2700
      IF (WORK.GT.0.0D0) GO TO 5                                        ZXCG2710
C                                  TERMINATE THE LINE SEARCH IF STEPCH  ZXCG2720
C                                    IS EFFECTIVELY ZERO.               ZXCG2730
      IF (NCALLS.GT.NFBEG+1) GO TO 115                                  ZXCG2740
      IF (DABS(GMIN/GINIT)-0.2D0) 170,170,115                           ZXCG2750
C                                  LET SPLN BE THE QUADRATIC SPLINE     ZXCG2760
C                                    THAT INTERPOLATES THE CALCULATED   ZXCG2770
C                                    FUNCTION VALUES AND DIRECTIONAL    ZXCG2780
C                                    DERIVATIVES AT THE POINTS STMIN    ZXCG2790
C                                    AND STEP OF THE LINE SEARCH, WHERE ZXCG2800
C                                    THE KNOT OF THE SPLINE IS AT       ZXCG2810
C                                    0.5*(STMIN+STEP). REVISE STMIN,    ZXCG2820
C                                    GMIN AND SBOUND, AND SET DDSPLN TO ZXCG2830
C                                    THE SECOND DERIVATIVE OF SPLN AT   ZXCG2840
C                                    THE NEW STMIN. HOWEVER, IF FCH IS  ZXCG2850
C                                    ZERO, IT IS ASSUMED THAT THE       ZXCG2860
C                                    MAXIMUM ACCURACY IS ALMOST         ZXCG2870
C                                    ACHIEVED, SO DDSPLN IS CALCULATED  ZXCG2880
C                                    USING ONLY THE CHANGE IN THE       ZXCG2890
C                                    GRADIENT.                          ZXCG2900
  100 WORK = (FCH+FCH)/STEPCH-GNEW-GMIN                                 ZXCG2910
      DDSPLN = (GNEW-GMIN)/STEPCH                                       ZXCG2920
      IF (NCALLS.GT.NFOPT) SBOUND = STEP                                ZXCG2930
      IF (NCALLS.GT.NFOPT) GO TO 105                                    ZXCG2940
      IF (GMIN*GNEW.LE.0.0D0) SBOUND = STMIN                            ZXCG2950
      STMIN = STEP                                                      ZXCG2960
      GMIN = GNEW                                                       ZXCG2970
      STEPCH = -STEPCH                                                  ZXCG2980
  105 IF (FCH.NE.0.0D0) DDSPLN = DDSPLN+(WORK+WORK)/STEPCH              ZXCG2990
C                                                                       ZXCG3000
C                                  TEST FOR CONVERGENCE OF THE LINE     ZXCG3010
C                                    SEARCH, BUT FORCE AT LEAST TWO     ZXCG3020
C                                    STEPS TO BE TAKEN IN ORDER NOT TO  ZXCG3030
C                                    LOSE QUADRATIC TERMINATION.        ZXCG3040
      IF (GMIN.EQ.0.0D0) GO TO 170                                      ZXCG3050
      IF (NCALLS.LE.NFBEG+1) GO TO 120                                  ZXCG3060
      IF (DABS(GMIN/GINIT).LE.0.2D0) GO TO 170                          ZXCG3070
C                                  APPLY THE TEST THAT DEPENDS ON THE   ZXCG3080
C                                    PARAMETER MAXLIN.                  ZXCG3090
  110 IF (NCALLS.LT.NFOPT+MAXLIN) GO TO 120                             ZXCG3100
  115 IER = 129                                                         ZXCG3110
      GO TO 170                                                         ZXCG3120
C                                  SET STEPCH TO THE GREATEST CHANGE TO ZXCG3130
C                                    THE CURRENT VALUE OF STMIN THAT IS ZXCG3140
C                                    ALLOWED BY THE BOUND ON THE LINE   ZXCG3150
C                                    SEARCH. SET GSPLN TO THE GRADIENT  ZXCG3160
C                                    OF THE QUADRATIC SPLINE AT         ZXCG3170
C                                    (STMIN+STEPCH). HENCE CALCULATE    ZXCG3180
C                                    THE VALUE OF STEPCH THAT MINIMIZES ZXCG3190
C                                    THE SPLINE FUNCTION, AND THEN      ZXCG3200
C                                    OBTAIN THE NEW FUNCTION AND        ZXCG3210
C                                    GRADIENT VECTOR, FOR THIS VALUE OF ZXCG3220
C                                    THE CHANGE TO THE STEP-LENGTH.     ZXCG3230
  120 STEPCH = 0.5D0*(SBOUND-STMIN)                                     ZXCG3240
      IF (SBOUND.LT.-0.5D0) STEPCH = 9.0D0*STMIN                        ZXCG3250
      GSPLN = GMIN+STEPCH*DDSPLN                                        ZXCG3260
      IF (GMIN*GSPLN.LT.0.0D0) STEPCH = STEPCH*GMIN/(GMIN-GSPLN)        ZXCG3270
      GO TO 90                                                          ZXCG3280
C                                  CALCULATE THE VALUE OF BETA THAT     ZXCG3290
C                                    OCCURS IN THE NEW SEARCH           ZXCG3300
C                                    DIRECTION.                         ZXCG3310
  125 SUM = 0.0D0                                                       ZXCG3320
      DO 130 I=1,N                                                      ZXCG3330
  130 SUM = SUM+G(I)*W(IGINIT+I)                                        ZXCG3340
      BETA = (GSQRD-SUM)/(GMIN-GINIT)                                   ZXCG3350
C                                  TEST THAT THE NEW SEARCH DIRECTION   ZXCG3360
C                                    CAN BE MADE DOWNHILL. IF IT        ZXCG3370
C                                    CANNOT, THEN MAKE ONE ATTEMPT TO   ZXCG3380
C                                    IMPROVE THE ACCURACY OF THE LINE   ZXCG3390
C                                    SEARCH.                            ZXCG3400
      IF (DABS(BETA*GMIN).LE.0.2D0*GSQRD) GO TO 135                     ZXCG3410
      IRETRY = IRETRY+1                                                 ZXCG3420
      IF (IRETRY.LE.0) GO TO 110                                        ZXCG3430
C                                  APPLY THE TEST THAT DEPENDS ON THE   ZXCG3440
C                                    PARAMETER MXFCON.                  ZXCG3450
C                                    SET DFPR TO THE PREDICTED          ZXCG3460
C                                    REDUCTION IN F ON THE NEXT         ZXCG3470
C                                    ITERATION.                         ZXCG3480
  135 IF (F.LT.FINIT) ITERFM = ITERC                                    ZXCG3490
      IF (ITERC.LT.ITERFM+MXFCON) GO TO 140                             ZXCG3500
      IER = 132                                                         ZXCG3510
      GO TO 9000                                                        ZXCG3520
  140 DFPR = STMIN*GINIT                                                ZXCG3530
C                                  BRANCH IF A RESTART PROCEDURE IS     ZXCG3540
C                                    REQUIRED DUE TO THE ITERATION      ZXCG3550
C                                    NUMBER OR DUE TO THE SCALAR        ZXCG3560
C                                    PRODUCT OF CONSECUTIVE GRADIENTS.  ZXCG3570
      IF (IRETRY.GT.0) GO TO 10                                         ZXCG3580
      IF (ITERRS.EQ.0) GO TO 155                                        ZXCG3590
      IF (ITERC-ITERRS.GE.N) GO TO 155                                  ZXCG3600
      IF (DABS(SUM).GE.0.2D0*GSQRD) GO TO 155                           ZXCG3610
C                                  CALCULATE THE VALUE OF GAMA THAT     ZXCG3620
C                                    OCCURS IN THE NEW SEARCH           ZXCG3630
C                                    DIRECTION, AND SET SUM TO A SCALAR ZXCG3640
C                                    PRODUCT FOR THE TEST BELOW. THE    ZXCG3650
C                                    VALUE OF GAMDEN IS SET BY THE      ZXCG3660
C                                    RESTART PROCEDURE.                 ZXCG3670
      GAMA = 0.0D0                                                      ZXCG3680
      SUM = 0.0D0                                                       ZXCG3690
      DO 145 I=1,N                                                      ZXCG3700
         GAMA = GAMA+G(I)*W(IRSDG+I)                                    ZXCG3710
  145 SUM = SUM+G(I)*W(IRSDX+I)                                         ZXCG3720
      GAMA = GAMA/GAMDEN                                                ZXCG3730
C                                  RESTART IF THE NEW SEARCH DIRECTION  ZXCG3740
C                                    IS NOT SUFFICIENTLY DOWNHILL.      ZXCG3750
C                                                                       ZXCG3760
      IF (DABS(BETA*GMIN+GAMA*SUM).GE.0.2D0*GSQRD) GO TO 155            ZXCG3770
C                                                                       ZXCG3780
C                                  CALCULATE THE NEW SEARCH DIRECTION.  ZXCG3790
      DO 150 I=1,N                                                      ZXCG3800
  150 W(I) = -G(I)+BETA*W(I)+GAMA*W(IRSDX+I)                            ZXCG3810
      GO TO 80                                                          ZXCG3820
C                                  APPLY THE RESTART PROCEDURE.         ZXCG3830
  155 GAMDEN = GMIN-GINIT                                               ZXCG3840
      DO 160 I=1,N                                                      ZXCG3850
         W(IRSDX+I) = W(I)                                              ZXCG3860
         W(IRSDG+I) = G(I)-W(IGINIT+I)                                  ZXCG3870
  160 W(I) = -G(I)+BETA*W(I)                                            ZXCG3880
      ITERRS = ITERC                                                    ZXCG3890
      GO TO 80                                                          ZXCG3900
C                                  SET IER TO INDICATE THAT THE SEARCH  ZXCG3910
C                                    DIRECTION IS UPHILL.               ZXCG3920
  165 IER = 130                                                         ZXCG3930
C                                  ENSURE THAT F, X AND G ARE OPTIMAL.  ZXCG3940
  170 IF (NCALLS.EQ.NFOPT) GO TO 180                                    ZXCG3950
      F = FMIN                                                          ZXCG3960
      DO 175 I=1,N                                                      ZXCG3970
         X(I) = W(IXOPT+I)                                              ZXCG3980
  175 G(I) = W(IGOPT+I)                                                 ZXCG3990
  180 IF (IER.EQ.0) GO TO 125                                           ZXCG4000
 9000 CONTINUE                                                          ZXCG4010
C      CALL UERTST (IER,'ZXCGR ')                                        ZXCG4020
9005    WRITE(*,'(1X,A,I8,7X,A,I8)') 'NCALLS =',NCALLS,'IER =',IER
C     WRITE (IDOUT,'(/, 1X, A, T44, 1I, /, 1X, A, T40, 1G14.6)')
C      WRITE (*,'(/, 1X, A, T44, 1I8, /, 1X, A, T40, 1G14.6)')
C     $     'NUMBER OF MINIMIZATION CALLS:', NCALLS,
C     $     'FINAL RESIDUAL GRADIENT:     ', SUM
      RETURN                                                            ZXCG4030
      END                                                               ZXCG4040

!***********************************************************************!
!***********************************************************************!
!***********************************************************************!
C     ----------------------------------------------------------------------
C     THIS FILE CONTAINS THE LBFGS ALGORITHM AND SUPPORTING ROUTINES
C
C     ****************
C     LBFGS SUBROUTINE
C     ****************
C
      SUBROUTINE LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
C
      INTEGER N,M,IPRINT(2),IFLAG
      DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
      DOUBLE PRECISION F,EPS,XTOL
      LOGICAL DIAGCO
C
C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** JULY 1990 ***
C
C
C     THIS SUBROUTINE SOLVES THE UNCONSTRAINED MINIMIZATION PROBLEM
C
C                      MIN F(X),    X= (X1,X2,...,XN),
C
C      USING THE LIMITED MEMORY BFGS METHOD. THE ROUTINE IS ESPECIALLY
C      EFFECTIVE ON PROBLEMS INVOLVING A LARGE NUMBER OF VARIABLES. IN
C      A TYPICAL ITERATION OF THIS METHOD AN APPROXIMATION HK TO THE
C      INVERSE OF THE HESSIAN IS OBTAINED BY APPLYING M BFGS UPDATES TO
C      A DIAGONAL MATRIX HK0, USING INFORMATION FROM THE PREVIOUS M STEPS.
C      THE USER SPECIFIES THE NUMBER M, WHICH DETERMINES THE AMOUNT OF
C      STORAGE REQUIRED BY THE ROUTINE. THE USER MAY ALSO PROVIDE THE
C      DIAGONAL MATRICES HK0 IF NOT SATISFIED WITH THE DEFAULT CHOICE.
C      THE ALGORITHM IS DESCRIBED IN "ON THE LIMITED MEMORY BFGS METHOD
C      FOR LARGE SCALE OPTIMIZATION", BY D. LIU AND J. NOCEDAL,
C      MATHEMATICAL PROGRAMMING B 45 (1989) 503-528.
C
C      THE USER IS REQUIRED TO CALCULATE THE FUNCTION VALUE F AND ITS
C      GRADIENT G. IN ORDER TO ALLOW THE USER COMPLETE CONTROL OVER
C      THESE COMPUTATIONS, REVERSE  COMMUNICATION IS USED. THE ROUTINE
C      MUST BE CALLED REPEATEDLY UNDER THE CONTROL OF THE PARAMETER
C      IFLAG.
C
C      THE STEPLENGTH IS DETERMINED AT EACH ITERATION BY MEANS OF THE
C      LINE SEARCH ROUTINE MCVSRCH, WHICH IS A SLIGHT MODIFICATION OF
C      THE ROUTINE CSRCH WRITTEN BY MORE' AND THUENTE.
C
C      THE CALLING STATEMENT IS
C
C          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
C
C      WHERE
C
C     N       IS AN INTEGER VARIABLE THAT MUST BE SET BY THE USER TO THE
C             NUMBER OF VARIABLES. IT IS NOT ALTERED BY THE ROUTINE.
C             RESTRICTION: N>0.
C
C     M       IS AN INTEGER VARIABLE THAT MUST BE SET BY THE USER TO
C             THE NUMBER OF CORRECTIONS USED IN THE BFGS UPDATE. IT
C             IS NOT ALTERED BY THE ROUTINE. VALUES OF M LESS THAN 3 ARE
C             NOT RECOMMENDED; LARGE VALUES OF M WILL RESULT IN EXCESSIVE
C             COMPUTING TIME. 3<= M <=7 IS RECOMMENDED. RESTRICTION: M>0.
C
C     X       IS A DOUBLE PRECISION ARRAY OF LENGTH N. ON INITIAL ENTRY
C             IT MUST BE SET BY THE USER TO THE VALUES OF THE INITIAL
C             ESTIMATE OF THE SOLUTION VECTOR. ON EXIT WITH IFLAG=0, IT
C             CONTAINS THE VALUES OF THE VARIABLES AT THE BEST POINT
C             FOUND (USUALLY A SOLUTION).
C
C     F       IS A DOUBLE PRECISION VARIABLE. BEFORE INITIAL ENTRY AND ON
C             A RE-ENTRY WITH IFLAG=1, IT MUST BE SET BY THE USER TO
C             CONTAIN THE VALUE OF THE FUNCTION F AT THE POINT X.
C
C     G       IS A DOUBLE PRECISION ARRAY OF LENGTH N. BEFORE INITIAL
C             ENTRY AND ON A RE-ENTRY WITH IFLAG=1, IT MUST BE SET BY
C             THE USER TO CONTAIN THE COMPONENTS OF THE GRADIENT G AT
C             THE POINT X.
C
C     DIAGCO  IS A LOGICAL VARIABLE THAT MUST BE SET TO .TRUE. IF THE
C             USER  WISHES TO PROVIDE THE DIAGONAL MATRIX HK0 AT EACH
C             ITERATION. OTHERWISE IT SHOULD BE SET TO .FALSE., IN WHICH
C             CASE  LBFGS WILL USE A DEFAULT VALUE DESCRIBED BELOW. IF
C             DIAGCO IS SET TO .TRUE. THE ROUTINE WILL RETURN AT EACH
C             ITERATION OF THE ALGORITHM WITH IFLAG=2, AND THE DIAGONAL
C              MATRIX HK0  MUST BE PROVIDED IN THE ARRAY DIAG.
C
C
C     DIAG    IS A DOUBLE PRECISION ARRAY OF LENGTH N. IF DIAGCO=.TRUE.,
C             THEN ON INITIAL ENTRY OR ON RE-ENTRY WITH IFLAG=2, DIAG
C             IT MUST BE SET BY THE USER TO CONTAIN THE VALUES OF THE
C             DIAGONAL MATRIX HK0.  RESTRICTION: ALL ELEMENTS OF DIAG
C             MUST BE POSITIVE.
C
C     IPRINT  IS AN INTEGER ARRAY OF LENGTH TWO WHICH MUST BE SET BY THE
C             USER.
C
C             IPRINT(1) SPECIFIES THE FREQUENCY OF THE OUTPUT:
C                IPRINT(1) < 0 : NO OUTPUT IS GENERATED,
C                IPRINT(1) = 0 : OUTPUT ONLY AT FIRST AND LAST ITERATION,
C                IPRINT(1) > 0 : OUTPUT EVERY IPRINT(1) ITERATIONS.
C
C             IPRINT(2) SPECIFIES THE TYPE OF OUTPUT GENERATED:
C                IPRINT(2) = 0 : ITERATION COUNT, NUMBER OF FUNCTION
C                                EVALUATIONS, FUNCTION VALUE, NORM OF THE
C                                GRADIENT, AND STEPLENGTH,
C                IPRINT(2) = 1 : SAME AS IPRINT(2)=0, PLUS VECTOR OF
C                                VARIABLES AND  GRADIENT VECTOR AT THE
C                                INITIAL POINT,
C                IPRINT(2) = 2 : SAME AS IPRINT(2)=1, PLUS VECTOR OF
C                                VARIABLES,
C                IPRINT(2) = 3 : SAME AS IPRINT(2)=2, PLUS GRADIENT VECTOR.
C
C
C     EPS     IS A POSITIVE DOUBLE PRECISION VARIABLE THAT MUST BE SET BY
C             THE USER, AND DETERMINES THE ACCURACY WITH WHICH THE SOLUTION
C             IS TO BE FOUND. THE SUBROUTINE TERMINATES WHEN
C
C                         ||G|| < EPS MAX(1,||X||),
C
C             WHERE ||.|| DENOTES THE EUCLIDEAN NORM.
C
C     XTOL    IS A  POSITIVE DOUBLE PRECISION VARIABLE THAT MUST BE SET BY
C             THE USER TO AN ESTIMATE OF THE MACHINE PRECISION (E.G.
C             10**(-16) ON A SUN STATION 3/60). THE LINE SEARCH ROUTINE WILL
C             TERMINATE IF THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C             IS LESS THAN XTOL.
C
C     W       IS A DOUBLE PRECISION ARRAY OF LENGTH N(2M+1)+2M USED AS
C             WORKSPACE FOR LBFGS. THIS ARRAY MUST NOT BE ALTERED BY THE
C             USER.
C
C     IFLAG   IS AN INTEGER VARIABLE THAT MUST BE SET TO 0 ON INITIAL ENTRY
C             TO THE SUBROUTINE. A RETURN WITH IFLAG<0 INDICATES AN ERROR,
C             AND IFLAG=0 INDICATES THAT THE ROUTINE HAS TERMINATED WITHOUT
C             DETECTING ERRORS. ON A RETURN WITH IFLAG=1, THE USER MUST
C             EVALUATE THE FUNCTION F AND GRADIENT G. ON A RETURN WITH
C             IFLAG=2, THE USER MUST PROVIDE THE DIAGONAL MATRIX HK0.
C
C             THE FOLLOWING NEGATIVE VALUES OF IFLAG, DETECTING AN ERROR,
C             ARE POSSIBLE:
C
C              IFLAG=-1  THE LINE SEARCH ROUTINE MCSRCH FAILED. THE
C                        PARAMETER INFO PROVIDES MORE DETAILED INFORMATION
C                        (SEE ALSO THE DOCUMENTATION OF MCSRCH):
C
C                       INFO = 0  IMPROPER INPUT PARAMETERS.
C
C                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
C                                 UNCERTAINTY IS AT MOST XTOL.
C
C                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE
C                                 REQUIRED AT THE PRESENT ITERATION.
C
C                       INFO = 4  THE STEP IS TOO SMALL.
C
C                       INFO = 5  THE STEP IS TOO LARGE.
C
C                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
C                                 THERE MAY NOT BE A STEP WHICH SATISFIES
C                                 THE SUFFICIENT DECREASE AND CURVATURE
C                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL.
C
C
C              IFLAG=-2  THE I-TH DIAGONAL ELEMENT OF THE DIAGONAL INVERSE
C                        HESSIAN APPROXIMATION, GIVEN IN DIAG, IS NOT
C                        POSITIVE.
C
C              IFLAG=-3  IMPROPER INPUT PARAMETERS FOR LBFGS (N OR M ARE
C                        NOT POSITIVE).
C
C
C
C    ON THE DRIVER:
C
C    THE PROGRAM THAT CALLS LBFGS MUST CONTAIN THE DECLARATION:
C
C                       EXTERNAL LB2
C
C    LB2 IS A BLOCK DATA THAT DEFINES THE DEFAULT VALUES OF SEVERAL
C    PARAMETERS DESCRIBED IN THE COMMON SECTION.
C
C
C
C    COMMON:
C
C     THE SUBROUTINE CONTAINS ONE COMMON AREA, WHICH THE USER MAY WISH TO
C    REFERENCE:
C
         COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
C
C    MP  IS AN INTEGER VARIABLE WITH DEFAULT VALUE 6. IT IS USED AS THE
C        UNIT NUMBER FOR THE PRINTING OF THE MONITORING INFORMATION
C        CONTROLLED BY IPRINT.
C
C    LP  IS AN INTEGER VARIABLE WITH DEFAULT VALUE 6. IT IS USED AS THE
C        UNIT NUMBER FOR THE PRINTING OF ERROR MESSAGES. THIS PRINTING
C        MAY BE SUPPRESSED BY SETTING LP TO A NON-POSITIVE VALUE.
C
C    GTOL IS A DOUBLE PRECISION VARIABLE WITH DEFAULT VALUE 0.9, WHICH
C        CONTROLS THE ACCURACY OF THE LINE SEARCH ROUTINE MCSRCH. IF THE
C        FUNCTION AND GRADIENT EVALUATIONS ARE INEXPENSIVE WITH RESPECT
C        TO THE COST OF THE ITERATION (WHICH IS SOMETIMES THE CASE WHEN
C        SOLVING VERY LARGE PROBLEMS) IT MAY BE ADVANTAGEOUS TO SET GTOL
C        TO A SMALL VALUE. A TYPICAL SMALL VALUE IS 0.1.  RESTRICTION:
C        GTOL SHOULD BE GREATER THAN 1.D-04.
C
C    STPMIN AND STPMAX ARE NON-NEGATIVE DOUBLE PRECISION VARIABLES WHICH
C        SPECIFY LOWER AND UPER BOUNDS FOR THE STEP IN THE LINE SEARCH.
C        THEIR DEFAULT VALUES ARE 1.D-20 AND 1.D+20, RESPECTIVELY. THESE
C        VALUES NEED NOT BE MODIFIED UNLESS THE EXPONENTS ARE TOO LARGE
C        FOR THE MACHINE BEING USED, OR UNLESS THE PROBLEM IS EXTREMELY
C        BADLY SCALED (IN WHICH CASE THE EXPONENTS SHOULD BE INCREASED).
C
C
C  MACHINE DEPENDENCIES
C
C        THE ONLY VARIABLES THAT ARE MACHINE-DEPENDENT ARE XTOL,
C        STPMIN AND STPMAX.
C
C
C  GENERAL INFORMATION
C
C    OTHER ROUTINES CALLED DIRECTLY:  DAXPY, DDOT, LB1, MCSRCH
C
C    INPUT/OUTPUT  :  NO INPUT; DIAGNOSTIC MESSAGES ON UNIT MP AND
C                     ERROR MESSAGES ON UNIT LP.
C
C
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      DOUBLE PRECISION GTOL,ONE,ZERO,GNORM,DDOT,STP1,FTOL,STPMIN,
     .                 STPMAX,STP,YS,YY,SQ,YR,BETA,XNORM
      INTEGER MP,LP,ITER,NFUN,POINT,ISPT,IYPT,MAXFEV,INFO,
     .        BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN
      LOGICAL FINISH
C
      SAVE
      DATA ONE,ZERO/1.0D+0,0.0D+0/
C
C     INITIALIZE
C     ----------
C
      IF(IFLAG.EQ.0) GO TO 10
      GO TO (172,100) IFLAG
  10  ITER= 0
      IF(N.LE.0.OR.M.LE.0) GO TO 196
      IF(GTOL.LE.1.D-04) THEN
        IF(LP.GT.0) WRITE(LP,245)
        GTOL=9.D-01
      ENDIF
      NFUN= 1
      POINT= 0
      FINISH= .FALSE.
      IF(DIAGCO) THEN
         DO 30 I=1,N
 30      IF (DIAG(I).LE.ZERO) GO TO 195
      ELSE
         DO 40 I=1,N
 40      DIAG(I)= 1.0D0
      ENDIF
C
C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
C     ---------------------------------------
C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
C         OTHER TEMPORARY INFORMATION.
C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
C         IN THE FORMULA THAT COMPUTES H*G.
C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
C         STEPS.
C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
C         GRADIENT DIFFERENCES.
C
C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
C
      ISPT= N+2*M
      IYPT= ISPT+N*M
      DO 50 I=1,N
 50   W(ISPT+I)= -G(I)*DIAG(I)
      GNORM= DSQRT(DDOT(N,G,1,G,1))
      STP1= ONE/GNORM
C
C     PARAMETERS FOR LINE SEARCH ROUTINE
C
      FTOL= 1.0D-4
      MAXFEV= 20
C
      IF(IPRINT(1).GE.0) CALL LB1(IPRINT,ITER,NFUN,
     *                     GNORM,N,M,X,F,G,STP,FINISH)
C
C    --------------------
C     MAIN ITERATION LOOP
C    --------------------
C
 80   ITER= ITER+1
      INFO=0
      BOUND=ITER-1
      IF(ITER.EQ.1) GO TO 165
      IF (ITER .GT. M)BOUND=M
C
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
      IF(.NOT.DIAGCO) THEN
         YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
         DO 90 I=1,N
   90    DIAG(I)= YS/YY
      ELSE
         IFLAG=2
         RETURN
      ENDIF
 100  CONTINUE
      IF(DIAGCO) THEN
        DO 110 I=1,N
 110    IF (DIAG(I).LE.ZERO) GO TO 195
      ENDIF
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: NOCEDAL, J. 1980,
C     "UPDATING QUASI-NEWTON MATRICES WITH LIMITED STORAGE",
C     MATHEMATICS OF COMPUTATION, VOL.24, NO.151, PP. 773-782.
C     ---------------------------------------------------------
C
      CP= POINT
      IF (POINT.EQ.0) CP=M
      W(N+CP)= ONE/YS
      DO 112 I=1,N
 112  W(I)= -G(I)
      CP= POINT
      DO 125 I= 1,BOUND
         CP=CP-1
         IF (CP.EQ. -1)CP=M-1
         SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
         INMC=N+M+CP+1
         IYCN=IYPT+CP*N
         W(INMC)= W(N+CP+1)*SQ
         CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
 125  CONTINUE
C
      DO 130 I=1,N
 130  W(I)=DIAG(I)*W(I)
C
      DO 145 I=1,BOUND
         YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
         BETA= W(N+CP+1)*YR
         INMC=N+M+CP+1
         BETA= W(INMC)-BETA
         ISCN=ISPT+CP*N
         CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
         CP=CP+1
         IF (CP.EQ.M)CP=0
 145  CONTINUE
C
C     STORE THE NEW SEARCH DIRECTION
C     ------------------------------
C
       DO 160 I=1,N
 160   W(ISPT+POINT*N+I)= W(I)
C
C     OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION
C     BY USING THE LINE SEARCH ROUTINE MCSRCH
C     ----------------------------------------------------
 165  NFEV=0
      STP=ONE
      IF (ITER.EQ.1) STP=STP1
      DO 170 I=1,N
 170  W(I)=G(I)
 172  CONTINUE
      CALL MCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP,FTOL,
     *            XTOL,MAXFEV,INFO,NFEV,DIAG)
      IF (INFO .EQ. -1) THEN
        IFLAG=1
        RETURN
      ENDIF
      IF (INFO .NE. 1) GO TO 190
      NFUN= NFUN + NFEV
C
C     COMPUTE THE NEW STEP AND GRADIENT CHANGE
C     -----------------------------------------
C
      NPT=POINT*N
      DO 175 I=1,N
      W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
 175  W(IYPT+NPT+I)= G(I)-W(I)
      POINT=POINT+1
      IF (POINT.EQ.M)POINT=0
C
C     TERMINATION TEST
C     ----------------
C
      GNORM= DSQRT(DDOT(N,G,1,G,1))
      XNORM= DSQRT(DDOT(N,X,1,X,1))
      XNORM= DMAX1(1.0D0,XNORM)
      IF (GNORM/XNORM .LE. EPS) FINISH=.TRUE.
C
      IF(IPRINT(1).GE.0) CALL LB1(IPRINT,ITER,NFUN,
     *               GNORM,N,M,X,F,G,STP,FINISH)
      IF (FINISH) THEN
         IFLAG=0
         RETURN
      ENDIF
      GO TO 80
C
C     ------------------------------------------------------------
C     END OF MAIN ITERATION LOOP. ERROR EXITS.
C     ------------------------------------------------------------
C
 190  IFLAG=-1
      IF(LP.GT.0) WRITE(LP,200) INFO
      RETURN
 195  IFLAG=-2
      IF(LP.GT.0) WRITE(LP,235) I
      RETURN
 196  IFLAG= -3
      IF(LP.GT.0) WRITE(LP,240)
C
C     FORMATS
C     -------
C
 200  FORMAT(/' IFLAG= -1 ',/' LINE SEARCH FAILED. SEE'
     .          ' DOCUMENTATION OF ROUTINE MCSRCH',/' ERROR RETURN'
     .          ' OF LINE SEARCH: INFO= ',I2,/
     .          ' POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT',/,
     .          ' OR INCORRECT TOLERANCES')
 235  FORMAT(/' IFLAG= -2',/' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,
     .       ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
 240  FORMAT(/' IFLAG= -3',/' IMPROPER INPUT PARAMETERS (N OR M',
     .       ' ARE NOT POSITIVE)')
 245  FORMAT(/'  GTOL IS LESS THAN OR EQUAL TO 1.D-04',
     .       / ' IT HAS BEEN RESET TO 9.D-01')
      RETURN
      END
C
C     LAST LINE OF SUBROUTINE LBFGS
C
C
      SUBROUTINE LB1(IPRINT,ITER,NFUN,
     *                     GNORM,N,M,X,F,G,STP,FINISH)
C
C     -------------------------------------------------------------
C     THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND
C     AMOUNT OF OUTPUT ARE CONTROLLED BY IPRINT.
C     -------------------------------------------------------------
C
      INTEGER IPRINT(2),ITER,NFUN,LP,MP,N,M
      DOUBLE PRECISION X(N),G(N),F,GNORM,STP,GTOL,STPMIN,STPMAX
      LOGICAL FINISH
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
C
      IF (ITER.EQ.0)THEN
           WRITE(MP,10)
           WRITE(MP,20) N,M
           WRITE(MP,30)F,GNORM
                 IF (IPRINT(2).GE.1)THEN
                     WRITE(MP,40)
                     WRITE(MP,50) (X(I),I=1,N)
                     WRITE(MP,60)
                     WRITE(MP,50) (G(I),I=1,N)
                  ENDIF
           WRITE(MP,10)
           WRITE(MP,70)
      ELSE
          IF ((IPRINT(1).EQ.0).AND.(ITER.NE.1.AND..NOT.FINISH))RETURN
              IF (IPRINT(1).NE.0)THEN
                   IF(MOD(ITER-1,IPRINT(1)).EQ.0.OR.FINISH)THEN
                         IF(IPRINT(2).GT.1.AND.ITER.GT.1) WRITE(MP,70)
                         WRITE(MP,80)ITER,NFUN,F,GNORM,STP
                   ELSE
                         RETURN
                   ENDIF
              ELSE
                   IF( IPRINT(2).GT.1.AND.FINISH) WRITE(MP,70)
                   WRITE(MP,80)ITER,NFUN,F,GNORM,STP
              ENDIF
              IF (IPRINT(2).EQ.2.OR.IPRINT(2).EQ.3)THEN
                    IF (FINISH)THEN
                        WRITE(MP,90)
                    ELSE
                        WRITE(MP,40)
                    ENDIF
                      WRITE(MP,50)(X(I),I=1,N)
                  IF (IPRINT(2).EQ.3)THEN
                      WRITE(MP,60)
                      WRITE(MP,50)(G(I),I=1,N)
                  ENDIF
              ENDIF
            IF (FINISH) WRITE(MP,100)
      ENDIF
C
 10   FORMAT('*************************************************')
 20   FORMAT('  N=',I5,'   NUMBER OF CORRECTIONS=',I2,
     .       /,  '       INITIAL VALUES')
 30   FORMAT(' F= ',1PD10.3,'   GNORM= ',1PD10.3)
 40   FORMAT(' VECTOR X= ')
 50   FORMAT(6(2X,1PD10.3))
 60   FORMAT(' GRADIENT VECTOR G= ')
 70   FORMAT(/'   I   NFN',4X,'FUNC',8X,'GNORM',7X,'STEPLENGTH'/)
 80   FORMAT(2(I4,1X),3X,3(1PD10.3,2X))
 90   FORMAT(' FINAL POINT X= ')
 100  FORMAT(/' THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.',
     .       /' IFLAG = 0')
C
      RETURN
      END
C     ******
C
C
C   ----------------------------------------------------------
C     DATA
C   ----------------------------------------------------------
C
      BLOCK DATA LB2
      INTEGER LP,MP
      DOUBLE PRECISION GTOL,STPMIN,STPMAX
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
      DATA MP,LP,GTOL,STPMIN,STPMAX/6,6,9.0D-01,1.0D-20,1.0D+20/
      END
C
C
C   ----------------------------------------------------------
C
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
C
C
C   ----------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
C    ------------------------------------------------------------------
C
C     **************************
C     LINE SEARCH ROUTINE MCSRCH
C     **************************
C
      SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL,MAXFEV,INFO,NFEV,WA)
      INTEGER N,MAXFEV,INFO,NFEV
      DOUBLE PRECISION F,STP,FTOL,GTOL,XTOL,STPMIN,STPMAX
      DOUBLE PRECISION X(N),G(N),S(N),WA(N)
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
      SAVE
C
C                     SUBROUTINE MCSRCH
C
C     A SLIGHT MODIFICATION OF THE SUBROUTINE CSRCH OF MORE' AND THUENTE.
C     THE CHANGES ARE TO ALLOW REVERSE COMMUNICATION, AND DO NOT AFFECT
C     THE PERFORMANCE OF THE ROUTINE.
C
C     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
C     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
C
C     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
C     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
C     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
C     MINIMIZER OF THE MODIFIED FUNCTION
C
C          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
C
C     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
C     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
C     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
C     CONTAINS A MINIMIZER OF F(X+STP*S).
C
C     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
C     THE SUFFICIENT DECREASE CONDITION
C
C           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),
C
C     AND THE CURVATURE CONDITION
C
C           ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).
C
C     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
C     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
C     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
C     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
C     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
C     SATISFIES THE SUFFICIENT DECREASE CONDITION.
C
C     THE SUBROUTINE STATEMENT IS
C
C        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF VARIABLES.
C
C       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
C         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
C         X + STP*S.
C
C       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
C         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
C
C       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
C         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
C         OF F AT X + STP*S.
C
C       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
C         SEARCH DIRECTION.
C
C       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
C         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
C         STP CONTAINS THE FINAL ESTIMATE.
C
C       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (IN THIS REVERSE
C         COMMUNICATION IMPLEMENTATION GTOL IS DEFINED IN A COMMON
C         STATEMENT.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
C         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
C         SATISFIED.
C
C       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
C         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C         IS AT MOST XTOL.
C
C       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
C         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (IN THIS REVERSE
C         COMMUNICATION IMPLEMENTATIN THEY ARE DEFINED IN A COMMON
C         STATEMENT).
C
C       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
C         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
C         MAXFEV BY THE END OF AN ITERATION.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
C
C         INFO = 0  IMPROPER INPUT PARAMETERS.
C
C         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
C
C         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
C                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
C
C         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C                   IS AT MOST XTOL.
C
C         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
C
C         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
C
C         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
C
C         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
C                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
C                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
C                   TOLERANCES MAY BE TOO SMALL.
C
C       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
C         CALLS TO FCN.
C
C       WA IS A WORK ARRAY OF LENGTH N.
C
C     SUBPROGRAMS CALLED
C
C       MCSTEP
C
C       FORTRAN-SUPPLIED...ABS,MAX,MIN
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
C     JORGE J. MORE', DAVID J. THUENTE
C
C     **********
      INTEGER INFOC,J
      LOGICAL BRACKT,STAGE1
      DOUBLE PRECISION DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,
     *       FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY,
     *       STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO
      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/
      IF(INFO.EQ.-1) GO TO 45
      INFOC = 1
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (N .LE. 0 .OR. STP .LE. ZERO .OR. FTOL .LT. ZERO .OR.
     *    GTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. STPMIN .LT. ZERO
     *    .OR. STPMAX .LT. STPMIN .OR. MAXFEV .LE. 0) RETURN
C
C     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
C     AND CHECK THAT S IS A DESCENT DIRECTION.
C
      DGINIT = ZERO
      DO 10 J = 1, N
         DGINIT = DGINIT + G(J)*S(J)
   10    CONTINUE
      IF (DGINIT .GE. ZERO) THEN
         WRITE(LP,15)
   15    FORMAT(/'  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION')
         RETURN
         ENDIF
C
C     INITIALIZE LOCAL VARIABLES.
C
      BRACKT = .FALSE.
      STAGE1 = .TRUE.
      NFEV = 0
      FINIT = F
      DGTEST = FTOL*DGINIT
      WIDTH = STPMAX - STPMIN
      WIDTH1 = WIDTH/P5
      DO 20 J = 1, N
         WA(J) = X(J)
   20    CONTINUE
C
C     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
C     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
C     THE INTERVAL OF UNCERTAINTY.
C     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
C
      STX = ZERO
      FX = FINIT
      DGX = DGINIT
      STY = ZERO
      FY = FINIT
      DGY = DGINIT
C
C     START OF ITERATION.
C
   30 CONTINUE
C
C        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
C        TO THE PRESENT INTERVAL OF UNCERTAINTY.
C
         IF (BRACKT) THEN
            STMIN = MIN(STX,STY)
            STMAX = MAX(STX,STY)
         ELSE
            STMIN = STX
            STMAX = STP + XTRAPF*(STP - STX)
            END IF
C
C        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
C
         STP = MAX(STP,STPMIN)
         STP = MIN(STP,STPMAX)
C
C        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
C        STP BE THE LOWEST POINT OBTAINED SO FAR.
C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     *      .OR. NFEV .GE. MAXFEV-1 .OR. INFOC .EQ. 0
     *      .OR. (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)) STP = STX
C
C        EVALUATE THE FUNCTION AND GRADIENT AT STP
C        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
C        WE RETURN TO MAIN PROGRAM TO OBTAIN F AND G.
C
         DO 40 J = 1, N
            X(J) = WA(J) + STP*S(J)
   40       CONTINUE
         INFO=-1
         RETURN
C
   45    INFO=0
         NFEV = NFEV + 1
         DG = ZERO
         DO 50 J = 1, N
            DG = DG + G(J)*S(J)
   50       CONTINUE
         FTEST1 = FINIT + STP*DGTEST
C
C        TEST FOR CONVERGENCE.
C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     *      .OR. INFOC .EQ. 0) INFO = 6
         IF (STP .EQ. STPMAX .AND.
     *       F .LE. FTEST1 .AND. DG .LE. DGTEST) INFO = 5
         IF (STP .EQ. STPMIN .AND.
     *       (F .GT. FTEST1 .OR. DG .GE. DGTEST)) INFO = 4
         IF (NFEV .GE. MAXFEV) INFO = 3
         IF (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX) INFO = 2
         IF (F .LE. FTEST1 .AND. ABS(DG) .LE. GTOL*(-DGINIT)) INFO = 1
C
C        CHECK FOR TERMINATION.
C
         IF (INFO .NE. 0) RETURN
C
C        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
C
         IF (STAGE1 .AND. F .LE. FTEST1 .AND.
     *       DG .GE. MIN(FTOL,GTOL)*DGINIT) STAGE1 = .FALSE.
C
C        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
C        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
C        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
C        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
C
         IF (STAGE1 .AND. F .LE. FX .AND. F .GT. FTEST1) THEN
C
C           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
C
            FM = F - STP*DGTEST
            FXM = FX - STX*DGTEST
            FYM = FY - STY*DGTEST
            DGM = DG - DGTEST
            DGXM = DGX - DGTEST
            DGYM = DGY - DGTEST
C
C           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C
            CALL MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,
     *                 BRACKT,STMIN,STMAX,INFOC)
C
C           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
C
            FX = FXM + STX*DGTEST
            FY = FYM + STY*DGTEST
            DGX = DGXM + DGTEST
            DGY = DGYM + DGTEST
         ELSE
C
C           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C
            CALL MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,
     *                 BRACKT,STMIN,STMAX,INFOC)
            END IF
C
C        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
C        INTERVAL OF UNCERTAINTY.
C
         IF (BRACKT) THEN
            IF (ABS(STY-STX) .GE. P66*WIDTH1)
     *         STP = STX + P5*(STY - STX)
            WIDTH1 = WIDTH
            WIDTH = ABS(STY-STX)
            END IF
C
C        END OF ITERATION.
C
         GO TO 30
C
C     LAST LINE OF SUBROUTINE MCSRCH.
C
      END
      SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
     *                 STPMIN,STPMAX,INFO)
      INTEGER INFO
      DOUBLE PRECISION STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
      LOGICAL BRACKT,BOUND
C
C     SUBROUTINE MCSTEP
C
C     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
C     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
C     A MINIMIZER OF THE FUNCTION.
C
C     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
C     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
C     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
C     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
C     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
C     WITH ENDPOINTS STX AND STY.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
C                        STPMIN,STPMAX,INFO)
C
C     WHERE
C
C       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
C         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
C         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
C         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
C         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
C
C       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
C         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
C         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
C         UPDATED APPROPRIATELY.
C
C       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
C         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
C         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
C         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
C
C       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
C         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
C         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
C         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
C
C       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
C         AND UPPER BOUNDS FOR THE STEP.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
C         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
C         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
C         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
C     JORGE J. MORE', DAVID J. THUENTE
C
      DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
      INFO = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF ((BRACKT .AND. (STP .LE. MIN(STX,STY) .OR.
     *    STP .GE. MAX(STX,STY))) .OR.
     *    DX*(STP-STX) .GE. 0.0 .OR. STPMAX .LT. STPMIN) RETURN
C
C     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
C
      SGND = DP*(DX/ABS(DX))
C
C     FIRST CASE. A HIGHER FUNCTION VALUE.
C     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
C     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
C     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
C
      IF (FP .GT. FX) THEN
         INFO = 1
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .LT. STX) GAMMA = -GAMMA
         P = (GAMMA - DX) + THETA
         Q = ((GAMMA - DX) + GAMMA) + DP
         R = P/Q
         STPC = STX + R*(STP - STX)
         STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)
         IF (ABS(STPC-STX) .LT. ABS(STPQ-STX)) THEN
            STPF = STPC
         ELSE
           STPF = STPC + (STPQ - STPC)/2
           END IF
         BRACKT = .TRUE.
C
C     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
C     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
C     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
C     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
C
      ELSE IF (SGND .LT. 0.0) THEN
         INFO = 2
         BOUND = .FALSE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = ((GAMMA - DP) + GAMMA) + DX
         R = P/Q
         STPC = STP + R*(STX - STP)
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN
            STPF = STPC
         ELSE
            STPF = STPQ
            END IF
         BRACKT = .TRUE.
C
C     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
C     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
C     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
C     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
C     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
C     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
C     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
C
      ELSE IF (ABS(DP) .LT. ABS(DX)) THEN
         INFO = 3
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
C
C        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
C        TO INFINITY IN THE DIRECTION OF THE STEP.
C
         GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = (GAMMA + (DX - DP)) + GAMMA
         R = P/Q
         IF (R .LT. 0.0 .AND. GAMMA .NE. 0.0) THEN
            STPC = STP + R*(STX - STP)
         ELSE IF (STP .GT. STX) THEN
            STPC = STPMAX
         ELSE
            STPC = STPMIN
            END IF
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (BRACKT) THEN
            IF (ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
               END IF
         ELSE
            IF (ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
               END IF
            END IF
C
C     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
C     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
C     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
C
      ELSE
         INFO = 4
         BOUND = .FALSE.
         IF (BRACKT) THEN
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
            IF (STP .GT. STY) GAMMA = -GAMMA
            P = (GAMMA - DP) + THETA
            Q = ((GAMMA - DP) + GAMMA) + DY
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
         ELSE IF (STP .GT. STX) THEN
            STPF = STPMAX
         ELSE
            STPF = STPMIN
            END IF
         END IF
C
C     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
C     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
C
      IF (FP .GT. FX) THEN
         STY = STP
         FY = FP
         DY = DP
      ELSE
         IF (SGND .LT. 0.0) THEN
            STY = STX
            FY = FX
            DY = DX
            END IF
         STX = STP
         FX = FP
         DX = DP
         END IF
C
C     COMPUTE THE NEW STEP AND SAFEGUARD IT.
C
      STPF = MIN(STPMAX,STPF)
      STPF = MAX(STPMIN,STPF)
      STP = STPF
      IF (BRACKT .AND. BOUND) THEN
         IF (STY .GT. STX) THEN
            STP = MIN(STX+0.66*(STY-STX),STP)
         ELSE
            STP = MAX(STX+0.66*(STY-STX),STP)
            END IF
         END IF
      RETURN
C
C     LAST LINE OF SUBROUTINE MCSTEP.
C
      END
!***********************************************************************!
