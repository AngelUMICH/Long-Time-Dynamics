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
      SUBROUTINE FINISH

      USE SYSTEM
      USE DEFECT
      USE NBLVAR
      USE POTTAB
      USE KMCVAR

!!DEALLOCATE ALL THE DYNAMIC MATRIX

      DEALLOCATE(RA)
      DEALLOCATE(RR)
      DEALLOCATE(PA)
      DEALLOCATE(PB)
      DEALLOCATE(PC)
      DEALLOCATE(PD)
      DEALLOCATE(PR)
      DEALLOCATE(DA)
      DEALLOCATE(DB)
      DEALLOCATE(DC)
      DEALLOCATE(DD)
      DEALLOCATE(DR)
      DEALLOCATE(EA)
      DEALLOCATE(EB)
      DEALLOCATE(EC)
      DEALLOCATE(ED)
      DEALLOCATE(RL)
      DEALLOCATE(NI)
      DEALLOCATE(NL)
      DEALLOCATE(FQ)
      DEALLOCATE(PO)
      DEALLOCATE(EM)
      DEALLOCATE(FF)
      DEALLOCATE(EE)
      DEALLOCATE(TV)
      DEALLOCATE(TI)
      DEALLOCATE(IV) 
      DEALLOCATE(II)
      DEALLOCATE(PV)
      DEALLOCATE(PI)
      DEALLOCATE(PX)
      DEALLOCATE(ND)
      DEALLOCATE(SR)
      DEALLOCATE(LL)
      DEALLOCATE(HH)
      DEALLOCATE(CP)
      DEALLOCATE(MP)
      DEALLOCATE(IC)

      DEALLOCATE(FPB)
      DEALLOCATE(FPC)
      DEALLOCATE(FPD)
      DEALLOCATE(FDB)
      DEALLOCATE(FDC)
      DEALLOCATE(FDD)
      DEALLOCATE(FEB)
      DEALLOCATE(FEC)
      DEALLOCATE(FED)

!!CLOSE ALL FILES

      CLOSE(10)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
      CLOSE(19)

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
      CLOSE(21)
      CLOSE(22)
      CLOSE(25)
      CLOSE(26)

      RETURN
      END SUBROUTINE FINISH
!***********************************************************************!

