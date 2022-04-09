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
!!READING SYSTEM INPUT FROM VARIOUS FILES

      SUBROUTINE SYS_INPUT

      USE SYSTEM
      USE POTTAB

      IMPLICIT NONE

      INTEGER I,J

!!PARAMETER INPUT FILE,MASTER CONTROL FILE

      OPEN(UNIT=10,FILE='par_init.dat')

!!FILE CONTAINS INTITIAL STRUCTURE

      OPEN(UNIT=15,FILE='pos_init.dat')

!!REFERENCE STRUCTURE FILE

      OPEN(UNIT=16,FILE='pos_refs.dat')

!!POTENTIAL FILE, IN A TABLE FORMAT

!!PAIR TERM

      OPEN(UNIT=17,FILE='pot_iapp.dat')

!!DENSITY TERM

      OPEN(UNIT=18,FILE='pot_iapn.dat')

!!EMBEDDED ENERGY TERM

      OPEN(UNIT=19,FILE='pot_iape.dat')

!!READING SYSTEM PARAMTER, SEE PAR_INIT.DAT FOR DETAILS

      READ(10,'(9(I10/),5(F10.3/),1(I10/),6(F10.3/),1(ES15.5))')
     #TJOB,TOPT,TPOT,TPBC,NTSTEP,NDEMAX,NAVMAX,NNLMAX,NDIMER,RLAT,RBOX,
     #PCUTR,DCUTR,DECUT,NDRMAX,RDIMER,FNRMAX,FNRMIN,DSTEPF,DSTEPM,TEMP,
     #TIME

!!READING INTITIAL STRUCTURE

      READ(15,*) NATOMC

      ALLOCATE(RA(3,NATOMC+NDEMAX))

      READ(15,*)(RA(1,I),RA(2,I),RA(3,I),I=1,NATOMC)

!!READING REFERENCE STRUCTURE

      READ(16,*) NATOMR

      ALLOCATE(RR(3,NATOMR))

      READ(16,*)(RR(1,I),RR(2,I),RR(3,I),I=1,NATOMR)

!!READING POTENTIAL FILES

      READ(17,*) RPCUT,SPCUT,RPMIN,PSTEP,NPNODE

      ALLOCATE(PA(NPNODE),PB(NPNODE),PC(NPNODE),PD(NPNODE),PR(NPNODE))

      READ(17,*)(PA(J),PB(J),PC(J),PD(J),PR(J),J=1,NPNODE)


      READ(18,*) RDCUT,SDCUT,RDMIN,DSTEP,NDNODE

      ALLOCATE(DA(NDNODE),DB(NDNODE),DC(NDNODE),DD(NDNODE),DR(NDNODE))

      READ(18,*)(DA(J),DB(J),DC(J),DD(J),DR(J),J=1,NDNODE)


      READ(19,*) RECUT,ESTEP,NENODE

      ALLOCATE(EA(NENODE),EB(NENODE),EC(NENODE),ED(NENODE))

      READ(19,*)(EA(J),EB(J),EC(J),ED(J),J=1,NENODE)

      RETURN

      END SUBROUTINE SYS_INPUT
!***********************************************************************!
!!SYSTEM INITITIALIZATION

      SUBROUTINE SYS_INIT

      USE SYSTEM
      USE DEFECT
      USE NBLVAR
      USE POTTAB
      USE KMCVAR

      IMPLICIT NONE

      INTEGER :: N
      INTEGER,ALLOCATABLE :: SEED(:)

!!INITIALIZE OUTPUT FILE

!!DIMER OUTPUT FILE,NOT USED ANYMORE

      OPEN(UNIT=21,FILE='std_kdim.out')

!!DEFECT STRUCTURE FILE

      OPEN(UNIT=22,FILE='std_kdef.out')

!!CONTAINS INITIAL, SADDLE, AND FINAL CONFIGURATION ENERGIES

      OPEN(UNIT=25,FILE='std_evo.out')

!!CONTAINS THE OPT STRUCTURE AFTER RELAXATION

      OPEN(UNIT=26,FILE='std_kopt.out')

!!EITHER RANDOM SEED OR FIXED RANDOM SEED

      CALL RANDOM_SEED()

!      SEED(1) = 214748722

      CALL RANDOM_SEED(SIZE=N)
      ALLOCATE(SEED(N))
      CALL RANDOM_SEED(GET=SEED)


      FDCONV     = 1.0D-5

!!HALF OF THE BOX SIZE, FOR PBC USE

      RBH1       = RBOX/2.0D0

      RBH2       =-RBOX/2.0D0

      NFRAME     = 0

!!POTENTIAL CUTOFF:PAIR TERM

      PCUTR      = PCUTR/RLAT

!!POTENTIAL CUTOFF:DENSITY TERM

      DCUTR      = DCUTR/RLAT

      PCUTS      = PCUTR*PCUTR

      DCUTS      = DCUTR*DCUTR

!!DEFECT CUTOFF

      DECUT      = DECUT*DECUT

!!LENGTH OF NEIBHBOR LIST

      RLIST      = 1.3D0*PCUTR

      SLIST      = RLIST*RLIST

!!THICKNESS OF THE SKIN FOR NL

      RSKIN      = RLIST-PCUTR

!!TIME OF LAST IRRADIATION SET TO INTIAL TIME

      TIME_LAST = TIME

!!ALLOCATE DYNAMIC MATRIX
!!FROM HERE,SEE DEFINE.F FOR THE MEANING OF VARIABLES

      ALLOCATE(RL(3,NATOMC));RL=0.0D0

      ALLOCATE(NI(NATOMC));NI=0

      ALLOCATE(NL(NNLMAX,NATOMC));NL=0

      ALLOCATE(FQ(NDIMER,NDEMAX));FQ=0.0D0

      ALLOCATE(PO(NDIMER,NDEMAX));PO=0.0D0

      ALLOCATE(EM(NDIMER,NDEMAX));EM=0.0D0

      ALLOCATE(OE(NDIMER,NDEMAX));OE=0.0D0

      ALLOCATE(FF(3,NATOMC));FF=0.0D0

      ALLOCATE(EE(3,NATOMC));EE=0.0D0

      ALLOCATE(TV(NATOMR));TV=0

      ALLOCATE(TI(NATOMC));TI=0

      ALLOCATE(IV(3,NDEMAX));IV=0

      ALLOCATE(II(3,NDEMAX));II=0

      ALLOCATE(PV(3,NDEMAX));PV=0.0D0

      ALLOCATE(PI(3,NDEMAX));PI=0.0D0

      ALLOCATE(PX(7,NDEMAX));PX=0.0D0

      ALLOCATE(OX(7,NDEMAX));OX=0.0D0

      ALLOCATE(ND(NDIMER,NDEMAX));ND=0

      ALLOCATE(ON(NDIMER,NDEMAX));ON=0

      ALLOCATE(SR(4,NAVMAX,NDIMER,NDEMAX));SR=0.0D0

      ALLOCATE(OR(4,NAVMAX,NDIMER,NDEMAX));OR=0.0D0

      IF (RANK .EQ. 0)THEN

        ALLOCATE(SRBUFF(4,NAVMAX,NDIMER,NDEMAX));SRBUFF=0.0D0

        ALLOCATE(EMBUFF(NDIMER,NDEMAX));EMBUFF=0.0D0

        ALLOCATE(NDBUFF(NDIMER,NDEMAX));NDBUFF=0

      ENDIF

!!VARIABLES ASSSOCIATED WITH MPI
!! DIMER LIMIT VARIABLES BEING SET UP FOR J LOOP

      DIMERPR = NDIMER/NPROCS

      LDIMERI = DIMERPR*RANK + 1

      UDIMERI = DIMERPR*(RANK + 1)

!!VARIABLES ASSOCIATED WITH DIMER

      RDIMER     = RDIMER/RLAT

      DSTEPF     = DSTEPF/RLAT

      DSTEPM     = DSTEPM/RLAT

!!VARIABLES ASSOCIATED POTENTIAL TABLE
!!UNIT IS IN LATTICE CONSTANT

      ALLOCATE(FPB(NPNODE));FPB=0.0D0

      ALLOCATE(FPC(NPNODE));FPC=0.0D0

      ALLOCATE(FPD(NPNODE));FPD=0.0D0

      RPCUT  = RPCUT/RLAT

      SPCUT  = SPCUT/RLAT**2

      RPMIN  = RPMIN/RLAT**2

      PSTEP  = 1.0D0/PSTEP/RLAT**2

      PSTEI  = 1.0D0/PSTEP

      PR     = PR/RLAT**2

      PB     = PB*RLAT**2

      PC     = PC*RLAT**4

      PD     = PD*RLAT**6

      FPB    = PB*2.0D0

      FPC    = PC*4.0D0

      FPD    = PD*6.0D0

      ALLOCATE(FDB(NDNODE));FDB=0.0D0

      ALLOCATE(FDC(NDNODE));FDC=0.0D0

      ALLOCATE(FDD(NDNODE));FDD=0.0D0

      RDCUT  = RDCUT/RLAT

      SDCUT  = SDCUT/RLAT**2

      RDMIN  = RDMIN/RLAT**2

      DSTEP  = 1.0D0/DSTEP/RLAT**2

      DSTEI  = 1.0D0/DSTEP

      DR     = DR/RLAT**2

      DB     = DB*RLAT**2

      DC     = DC*RLAT**4

      DD     = DD*RLAT**6

      FDB    = DB*2.0D0

      FDC    = DC*4.0D0

      FDD    = DD*6.0D0

      ALLOCATE(FEB(NENODE));FEB=0.0D0

      ALLOCATE(FEC(NENODE));FEC=0.0D0

      ALLOCATE(FED(NENODE));FED=0.0D0

      ESTEI  = 1.0D0/ESTEP

      FEB    = EB*1.0D0

      FEC    = EC*2.0D0

      FED    = ED*3.0D0

!!VARIABLES ASSOCIATED WITH LINKCELL

      CSIZE    = 2.0D0

      MSIZE    = RBOX/CSIZE

      MSISQ    = MSIZE**2

      NCELL    = MSIZE**3

      NMAPS    = NCELL*26

      ALLOCATE(LL(NATOMC))

      ALLOCATE(HH(NCELL))

      ALLOCATE(CP(3,NCELL))

      ALLOCATE(MP(NMAPS))

      ALLOCATE(IC(MSIZE,MSIZE,MSIZE))

      CALL MAPCELL

      RETURN

      END SUBROUTINE SYS_INIT
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
