C
      SUBROUTINE PREPOT
C
C   System:          H3
C   Common name:     LSTH surface b
C   Cross reference: D. G. Truhlar and C. J. Horowitz
C                    J. Chem. Phys. 68, 2466 (1978); 71, 1514(E) (1979)
C
C   Note: This is not the original LSTH potential energy surface described
C         in Liu et al.  This version of the potential uses Tully's global
C         fit to H2 instead of a spline fit.                    
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in block data subprogram PTPACM
C   Coordinates, potential energy, and derivatives are passed 
C   through the common block PT31CM:
C                  COMMON /PT31CM/ RSV(3), ENERGY, DEDR(3)
C   The potential energy in the three asymptotic valleys are 
C   stored in the common block PT35CM:
C                  COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
C   The potential energy in the AB valley, EASYAB, is equal to the potential 
C   energy of the H "infinitely" far from the H2 diatomic, with the 
C   H2 diatomic at its equilibrium configuration.  For this potential energy
C   surface the potential energy in the three asymptotic valleys are equivalent.
C   All the information passed through the common blocks PT31CM and PT35CM
C   is in Hartree atomic units.  
C
C   This potential is written such that:
C                  RSV(1) = R(H-H)
C                  RSV(2) = R(H-H)
C                  RSV(3) = R(H-H)
C   The zero of energy is defined at H "infinitely" far from the H2 diatomic.
C
C   The flags that indicate what calculations should be carried out in 
C   the potential routine are passed through the common block PT32CM:
C                  /PT32CM/ NSURF, NDER, NFLAG(20)
C   where:
C        NSURF - which electronic state should be used.
C                This option is not used for this potential as only the 
C                ground electronic state is available.
C        NDER  = 0 => no derivatives should be calculated
C        NDER  = 1 => calculate first derivatives
C        NFLAG - these integer values can be used to flag options
C                within the potential; in this potential these options 
C                are not used.
C
C   The common block PT34CM contains the FORTRAN unit numbers for the 
C   potential output.  In this potential PT34CM contains one variable, IPRT,
C                      /PT34CM/ IPRT
C
C   Potential parameters' default settings
C                  Variable            Default value
C                  NSURF               0
C                  NDER                1
C                  IPRT                6
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT31CM/ RSV(3), ENERGY, DEDR(3)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /PT34CM/ IPRT
      COMMON /PT35CM/ EASYAB, EASYBC, EASYAC
C
      DIMENSION S1(3),S2(3),S3(3)
      COMMON /POTCOM/ C6,C8
      COMMON/VCOM/C,A,A1,F,FNS,F1,F2,F3,AN1,AN2,AN3,AN4,B1,B2,B3,
     1            W1,W2,W3,D1,D2,D3,D4,XL1,XL2
C
C   Echo the name of the potential to the file linked to the FORTRAN unit IPRT
C
      WRITE (IPRT, 600)
600   FORMAT (/,2X,T5,'PREPOT has been called for the H3 ',
     *                'potential energy surface LSTH surface b')
C
      R1 = 1.40105D0
      IC = 0
    5 CALL VTULL(R1,S1)
      IF(ABS(S1(2)) .LT. 1.D-10) GO TO 10
      DE1 = S1(2)
      IC = IC + 1
      IF(IC.GT.1) GO TO 6
      R2 = R1
      R1 = R1 + 0.0001D0
      DE2 = DE1
      GO TO 5
    6 IF(IC.LT.30) GO TO 7
      WRITE(IPRT,6000)
6000  FORMAT(/,2X,T5,'Error: In PREPOT the minimum of VH2 could ',
     *               'not be found')
      STOP 'PREPOT 1'
    7 R = (R2*DE1 - R1*DE2)/(DE1-DE2)
      R2 = R1
      R1 = R
      DE2 = DE1
      GO TO 5
   10 EZ = S1(1)
C
C   Initialize the energy in the asymptotic valleys
C   The energy in the asymptotic valley is set equal to -EZ because the
C   potential energy is defined as ENERGY  = E - EZ in this routine,
C   but the convention is ENERGY = E + EZ
C
      EASYAB = -EZ
      EASYBC = EASYAB
      EASYAC = EASYAB
C
      RETURN
C
      ENTRY POT
C
C   Check the values of NSURF and NDER for validity.
C
         IF (NSURF .NE. 0) THEN
             WRITE (IPRT, 900) NSURF
             STOP 'POT 1'
         ENDIF
         IF (NDER .GT. 1) THEN
             WRITE (IPRT, 910) NDER
             STOP 'POT 2'
         ENDIF
C
C E LONDON
C
      EF1=EXP(-F*RSV(1))
      EF2=EXP(-F*RSV(2))
      EF3=EXP(-F*RSV(3))
      X21=RSV(1)*RSV(1)
      X22=RSV(2)*RSV(2)
      X23=RSV(3)*RSV(3)
      T1=C*(A+RSV(1)+A1*X21)*EF1
      T2=C*(A+RSV(2)+A1*X22)*EF2
      T3=C*(A+RSV(3)+A1*X23)*EF3
      CALL VH2(RSV,S1,S2,S3)
      XQ1=S1(1)+T1
      XQ2=S2(1)+T2
      XQ3=S3(1)+T3
      XJ1=S1(1)-T1
      XJ2=S2(1)-T2
      XJ3=S3(1)-T3
      XQ=(XQ1+XQ2+XQ3)/2.D0
      XJ=SQRT(((XJ1-XJ2)**2+(XJ2-XJ3)**2+(XJ3-XJ1)**2)/8.D0)
      ELOND=XQ-XJ
C
C ENS
C
      WNT=(RSV(1)-RSV(2))*(RSV(2)-RSV(3))*(RSV(3)-RSV(1))
      WN=ABS(WNT)
      WN2=WN*WN
      WN3=WN2*WN
      WN4=WN3*WN
      WN5=WN4*WN
      R=RSV(1)+RSV(2)+RSV(3)
      R2=R*R
      R3=R2*R
      EXNS=EXP(-FNS*R3)
      ENS=(AN1*WN2+AN2*WN3+AN3*WN4+AN4*WN5)*EXNS
C
C NONLINEAR CORRECTIONS
C
      COS =(X21+X22+X23)/2.D0
      COS1=(X21-COS)/RSV(2)/RSV(3)
      COS2=(X22-COS)/RSV(1)/RSV(3)
      COS3=(X23-COS)/RSV(1)/RSV(2)
      WB=COS1+COS2+COS3+1.D0
      WB2=WB*WB
      WB3=WB2*WB
      WB4=WB3*WB
      EXF1=EXP(-F1*R)
      EXF2=EXP(-F2*R2)
      EXF3=EXP(-F3*R)
      EB1T=(B1+B2*R)*EXF1
      EB3T=(XL1+XL2*R2)*EXF3
      EB1=WB*(EB1T+EB3T)
      EB2=(WB2*W1+WB3*W2+WB4*W3)*EXF2
      EQ=(RSV(1)-RSV(2))**2+(RSV(2)-RSV(3))**2+(RSV(3)-RSV(1))**2
      RI=1.D0/RSV(1)+1.D0/RSV(2)+1.D0/RSV(3)
      EB4A=WB*D1*EXF1+WB2*D2*EXF2
      EB4B=D3*EXF1+D4*EXF2
      EB4=EB4A*RI+EB4B*WB*EQ
C
      ENERGY = ELOND+ENS+EB1+EB2+EB4 - EZ
C
C DERIVATIVES
C
C   If NDER =/= 1, then the derivatives of the energy
C   with respect to the coordinates are not computed.
C
      IF (NDER .NE. 1) GO TO 700
C
C E LONDON DERIVATIVES
C
      IF (XJ.EQ.0.0D0) GO TO 1
      XJS=(XJ1+XJ2+XJ3)/8.D0
      T1P=C*(1.D0+2.D0*A1*RSV(1))*EF1-F*T1
      T2P=C*(1.D0+2.D0*A1*RSV(2))*EF2-F*T2
      T3P=C*(1.D0+2.D0*A1*RSV(3))*EF3-F*T3
      ELON1P=(S1(2)+T1P)/2.D0-(S1(2)-T1P)*(.375D0*XJ1-XJS)/XJ
      ELON2P=(S2(2)+T2P)/2.D0-(S2(2)-T2P)*(.375D0*XJ2-XJS)/XJ
      ELON3P=(S3(2)+T3P)/2.D0-(S3(2)-T3P)*(.375D0*XJ3-XJS)/XJ
C
C ENS DERIVATIVES
C
      ENSPWN=(AN1*WN*2.D0+AN2*3.D0*WN2+AN3*4.D0*WN3+AN4*5.D0*WN4)*EXNS
      ENSPR=(-3.D0)*FNS*R2*ENS
C
C WN DERIVATIVES
C
      DELTA =-1.D0
      IF (WN.EQ.WNT) DELTA =1.D0
      WNP1=(2.D0*RSV(1)*(RSV(3)-RSV(2))+X22-X23)*DELTA
      WNP2=(2.D0*RSV(2)*(RSV(1)-RSV(3))+X23-X21)*DELTA
      WNP3=(2.D0*RSV(3)*(RSV(2)-RSV(1))+X21-X22)*DELTA
C
C DENS/DXI=(DENS/DWN)(DWN/DXI)+(DENS/DR)(DR/DXI)
C
      ENSP1=ENSPWN*WNP1+ENSPR
      ENSP2=ENSPWN*WNP2+ENSPR
      ENSP3=ENSPWN*WNP3+ENSPR
C
C WB DERIVATIVES
C
      WB1P=(RSV(1)/RSV(3)-1.D0)/RSV(2)-1.D0/RSV(3)-(COS2+COS3)/RSV(1)
      W23P=(RSV(3)/RSV(2)-1.D0)/RSV(1)-1.D0/RSV(2)-(COS1+COS2)/RSV(3)
      WB2P=(RSV(2)/RSV(3)-1.D0)/RSV(1)-1.D0/RSV(3)-(COS1+COS3)/RSV(2)
C
      WB3P=(RSV(3)/RSV(2)-1.D0)/RSV(1)-1.D0/RSV(2)-(COS1+COS2)/RSV(3)
C
C DEB1/DXI = (DEB1/DWB)(DWB/DXI)+(DEB1/DR)(DR/DXI)
C
      EB1PR=WB*(F1*EB1T+F3*EB3T-B2*EXF1-2.D0*R*XL2*EXF3)
      EB1PWB=EB1T+EB3T
      EB1P1=EB1PWB*WB1P-EB1PR
      EB1P2=EB1PWB*WB2P-EB1PR
      EB1P3=EB1PWB*WB3P-EB1PR
      EB2PWB=(2.D0*WB*W1+3.D0*WB2*W2+4.D0*WB3*W3)*EXF2
      EB2PR=F2*(-2.D0)*R*EB2
      EB2P1=EB2PWB*WB1P+EB2PR
      EB2P2=EB2PWB*WB2P+EB2PR
      EB2P3=EB2PWB*WB3P+EB2PR
C
C DEB4A/DXI= (DEB4A/DWB)(DWB/DXI)+(DEB4A/DRI)(DRI/DXI)+
C             (DEB4A/DR)(DR/DXI)
C
      EB4APW=(D1*EXF1+2.D0*WB*D2*EXF2)*RI
      EB4APR=RI*(WB2*F2*(-2.D0)*R*D2*EXF2-F1*D1*WB*EXF1)
      EB4AP1=EB4APW*WB1P-EB4A/X21+EB4APR
      EB4AP2=EB4APW*WB2P-EB4A/X22+EB4APR
      EB4AP3=EB4APW*WB3P-EB4A/X23+EB4APR
      EB4BPW=EB4B*EQ
      B4BPEQ=EB4B*WB
      EB4BPR=EQ*WB*((-2.D0)*F2*R*D4*EXF2-F1*D3*EXF1)
      EB4BP1=EB4BPW*WB1P+B4BPEQ*(6.D0*RSV(1)-2.D0*R)+EB4BPR
      EB4BP2=EB4BPW*WB2P+B4BPEQ*(6.D0*RSV(2)-2.D0*R)+EB4BPR
      EB4BP3=EB4BPW*WB3P+B4BPEQ*(6.D0*RSV(3)-2.D0*R)+EB4BPR
C
      DEDR(1)=ELON1P+ENSP1+EB1P1+EB2P1+EB4AP1+EB4BP1
      DEDR(2)=ELON2P+ENSP2+EB1P2+EB2P2+EB4AP2+EB4BP2
      DEDR(3)=ELON3P+ENSP3+EB1P3+EB2P3+EB4AP3+EB4BP3
      RETURN
    1 WRITE (IPRT,2)
2     FORMAT(/,2X,T5,'Warning: The coordinates form an equilateral ',
     *               'triangle and ',
     *       /,2X,T14,'the derivatives will be infinite.',
     *       /,2X,T14,'The derivatives at this geometry have been set ',
     *                'equal to zero.')
      DEDR(1) = 0.0D0
      DEDR(2) = 0.0D0
      DEDR(3) = 0.0D0
700   CONTINUE
C
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not allowed in this ',
     *              'version of the potential.')
C
      RETURN
      END
C
      SUBROUTINE VH2(X,S1,S2,S3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /POTCOM/ C6,C8
      DIMENSION X(3),S1(3),S2(3),S3(3)
      CALL VTULL(X(1),S1)
      CALL VTULL(X(2),S2)
      CALL VTULL(X(3),S3)
      RETURN
      END
C
      SUBROUTINE VTULL(R,V)
C
C   TULLY'S GLOBAL FIT TO KOLOS-WOLNEIWICZ H2
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(3)
      DATA A1,A2,A3,A4,A5,A6/149.480D0,-59.6557D0,-4.13792D0,
     1                       -23.7299D0,3.91747D0,-1.41350D0/
      DATA EEV/27.21161D0/
      EX1 = EXP(A3*R)
      EX2 = EXP(A6*R)
      T1 = A1 + A2*R
      T2 = R*R*(A4+A5*R)
      V(1) = (T1*EX1 + T2*EX2)/EEV
      V(2) = ((A2+A3*T1)*EX1 + (R*(2.0D0*A4+3.0D0*A5*R)+A6*T2)*EX2)/
     1       EEV
      RETURN
      END
C
      BLOCK DATA PTPACM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /PT34CM/ IPRT
      COMMON /POTCOM/ C6,C8
      COMMON/VCOM/C,A,A1,F,FN,FB1,FB2,FB3,XN1,XN2,XN3,XN4,B1,B2,B3
     1             ,G1,G2,G3,D1,D2,D3,D4,XL1,XL2
C
C   Initialize the control flags for the potential
C
      DATA IPRT /6/
      DATA NDER /1/
      DATA NFLAG /20*0/
C
C   Initialize the parameters for the potential
C
      DATA C6,C8/6.89992032D0,219.9997304D0/
      DATA C,A,A1,F/-1.2148730613D0,-1.514663474D0,-1.46D0,2.088442D0/
      DATA FN,XN1,XN2,XN3,XN4/.0035D0,.0012646477D0,-.0001585792D0,
     1 .0000079707D0,-.0000001151D0/
      DATA FB1,B1,B2/.52D0,3.0231771503D0,-1.08935219D0/
      DATA FB2,G1,G2,G3/.052D0,1.7732141742D0,-2.0979468223D0,
     1 -3.978850217D0/
      DATA D1,D2,D3,D4/.4908116374D0,-.8718696387D0,.1612118092D0,
     1 -.1273731045D0/
      DATA FB3,XL2,XL1/.79D0,.9877930913D0,-13.3599568553D0/
      END
