C
      SUBROUTINE PREPOT
C
C   System:          H3
C   Common name:     LSTH surface a
C   Cross reference: D. G Truhlar and C. J. Horowitz
C                    J. Chem. Phys. 68, 2466 (1978); 71, 1514(E) (1979)
C
C   Notes:  This is not the original LSTH potential energy surface described
C           in the Liu et al. reference.  In this routine the spline fit for
C           the H2 diatomic is replaced by a 6th. order polynomial fit to 
C           Kolos and Wolneiwicz.
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in the block data subprogram PTPACM
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
      DIMENSION S1(3),S2(3),S3(3)
      COMMON /POTCOM/ C6,C8,RKW(85),EKW(85),DEKW(85),COEF(6,83),NTABL,
     1                NTABL1
      COMMON/VCOM/C,A,A1,F,FNS,F1,F2,F3,AN1,AN2,AN3,AN4,B1,B2,B3,
     1            W1,W2,W3,D1,D2,D3,D4,XL1,XL2
C
C   Echo the name of the potential routine to the file linked to FORTRAN unit
C 
      WRITE (IPRT, 600)
600   FORMAT (/,2X,T5,'PREPOT has been called for the H3 ',
     *                'potential energy surface LSTH surface a')
C
      CALL VSET
      R1 = 1.40105D0
      IC = 0
    5 CALL VKW(R1,S1)
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
C   If NDER =/= 1, then the derivatives of the energy with respect to the 
C   coordinates are not computed.
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
      COMMON /POTCOM/ C6,C8,RKW(85),EKW(85),DEKW(85),COEF(6,83),NTABL,
     1                NTABL1
      DIMENSION X(3),S1(3),S2(3),S3(3)
1     IF (X(1).GT.10.D0) CALL VBIGR(X(1),S1)
      IF (X(1).GT.10.D0) GO TO 2
      CALL VKW(X(1),S1)
    2 IF (X(2).GT.10.D0) CALL VBIGR(X(2),S2)
      IF (X(2).GT.10.D0) GO TO 3
      CALL VKW(X(2),S2)
    3 IF (X(3).GT.10.D0) CALL VBIGR(X(3),S3)
      IF (X(3).GT.10.D0) RETURN
      CALL VKW(X(3),S3)
      RETURN
      END
C
      SUBROUTINE VKW(R,V)
C
C   6TH ORDER POLYNOMIAL FIT OF KOLOS-WOLNEIWICZ H2 DATA
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(3)
      COMMON /POTCOM/ C6,C8,RKW(85),EKW(85),DEKW(85),COEF(6,83),NTABL,
     1                NTABL1
      K4 = -1
      DO 10 I = 2,NTABL1
      K4 = K4 + 1
      IF(R.LT.RKW(I)) GO TO 20
   10 CONTINUE
      K4 = K4 + 1
      GO TO 30
   20 K4 = MAX0(1,K4)
      IF((R-RKW(K4+1)).GT.(RKW(K4+2)-R)) K4 = K4 + 1
   30 V(1) = COEF(1,K4)
      T = 5.0D0
      V(2) = T*V(1)
      DO 40 J = 2,5
      T = T - 1.0D0
      V(1) = V(1)*R + COEF(J,K4)
   40 V(2) = V(2)*R + T*COEF(J,K4)
      V(1) = V(1)*R + COEF(6,K4)
      RETURN
      END
C
      SUBROUTINE VSET
C
C   SOLVE FOR COEFFIECIENTS OF 6TH ORDER POLYNOMIAL 
C   FIT TO KOLOS-WOLNIEWICZ H2 DATA
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AB(6,7),TEMP(7)
      COMMON /PT34CM/ IPRT
      COMMON /POTCOM/ C6,C8,RKW(85),EKW(85),DEKW(85),COEF(6,83),NTABL,
     1                NTABL1
      IND=NTABL-2
      NTABL1 = NTABL - 1
      ID=6
      N=6
      EPS=1.D-14
      JRANK=N
      M=1
      K5=0
   10 K5=K5+1
      IF(K5.GT.IND) GO TO 999
      K4=K5-1
      DO 80 LD=1,5,2
      K4=K4+1
      AB(LD,6)=1.0D0
      AB(LD+1,6)=0.0D0
      X1=RKW(K4)
      AB(LD,7)=EKW(K4)
      AB(LD+1,7)=DEKW(K4)
      DO 60 LM=1,5
      LM6M=6-LM
   60 AB(LD,LM6M)=X1*AB(LD,LM6M+1)
      AB(LD+1,5)=1.0D0
      DO  80 LM=1,4
      LM5M=5-LM
   80 AB(LD+1,LM5M)=AB(LD,LM5M+1)*DBLE(LM+1)
C
      CALL MXLNEQ(AB,N,ID,DET,JRANK,EPS,TEMP,M)
      IF(JRANK .NE. N) GO TO 1000
      DO 90 J=1,6
   90 COEF(J,K5)=AB(J,7)
      GO TO 10
  999 RETURN
 1000 WRITE(IPRT,4) K4
      STOP 'VSET 1'
4     FORMAT(/,2X,T5,'Error: Unsuccessful matrix operation for ',
     *               'K4 = ',I3)
      END
C
      SUBROUTINE VBIGR(X,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /POTCOM/ C6,C8,RKW(85),EKW(85),DEKW(85),COEF(6,83),NTABL,
     1                NTABL1
      DIMENSION S(3)
      X2=X*X
      X3=X2*X
      X6=X3*X3
      C8A=C8/X2
      S(1) = -(C6+C8A)/X6
      S(2)=(C6*6.D0+C8A*8.D0)/X6/X
      RETURN
      END
C
      BLOCK DATA PTPACM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT32CM/ NSURF, NDER, NFLAG(20)
      COMMON /PT34CM/ IPRT
      COMMON /POTCOM/ C6,C8,RKW(85),EKW(85),DEKW(85),COEF(6,83),NTABL,
     1                NTABL1
      COMMON/VCOM/C,A,A1,F,FN,FB1,FB2,FB3,XN1,XN2,XN3,XN4,B1,B2,B3,
     1            G1,G2,G3,D1,D2,D3,D4,XL1,XL2
C
C    Initialize the control flags for the potential routine
C
         DATA IPRT /6/
         DATA NDER /1/
         DATA NFLAG /20*0/
C
C   Initialize the parameters for the potential routine
C
      DATA C6,C8/6.89992032D0,219.9997304D0/
      DATA C,A,A1,F/-1.2148730613D0,-1.514663474D0,-1.46D0,2.088442D0/
      DATA FN,XN1,XN2,XN3,XN4/.0035D0,.0012646477D0,-.0001585792D0,
     1                        .0000079707D0,-.0000001151D0/
      DATA FB1,B1,B2/.52D0,3.0231771503D0,-1.08935219D0/
      DATA FB2,G1,G2,G3/.052D0,1.7732141742D0,-2.0979468223D0,
     1                -3.978850217D0/
      DATA D1,D2,D3,D4/.4908116374D0,-.8718696387D0,.1612118092D0,
     1                -.1273731045D0/
      DATA FB3,XL2,XL1/.79D0,.9877930913D0,-13.3599568553D0/
C
C     KOLOS-WOLNEIWICZ H2 DATA
C
      DATA NTABL/ 85/
      DATA (RKW(I), I = 1, 42) /       4.00000000D-01, 4.50000000D-01,
     1 5.00000000D-01, 5.50000000D-01, 6.00000000D-01, 6.50000000D-01,
     1 7.00000000D-01, 7.50000000D-01, 8.00000000D-01, 9.00000000D-01,
     1 1.00000000D+00, 1.10000000D+00, 1.20000000D+00, 1.30000000D+00,
     1 1.35000000D+00, 1.39000000D+00, 1.40110000D+00, 1.41000000D+00,
     1 1.45000000D+00, 1.50000000D+00, 1.60000000D+00, 1.70000000D+00,
     1 1.80000000D+00, 1.90000000D+00, 2.00000000D+00, 2.10000000D+00,
     1 2.20000000D+00, 2.30000000D+00, 2.40000000D+00, 2.50000000D+00,
     1 2.60000000D+00, 2.70000000D+00, 2.80000000D+00, 2.90000000D+00,
     1 3.00000000D+00, 3.10000000D+00, 3.20000000D+00, 3.30000000D+00,
     1 3.40000000D+00, 3.50000000D+00, 3.60000000D+00, 3.70000000D+00/
       DATA (RKW(I), I = 43, 85) /
     1 3.80000000D+00, 3.90000000D+00, 4.00000000D+00, 4.10000000D+00,
     1 4.20000000D+00, 4.30000000D+00, 4.40000000D+00, 4.50000000D+00,
     1 4.60000000D+00, 4.70000000D+00, 4.80000000D+00, 4.90000000D+00,
     1 5.00000000D+00, 5.10000000D+00, 5.20000000D+00, 5.30000000D+00,
     1 5.40000000D+00, 5.50000000D+00, 5.60000000D+00, 5.70000000D+00,
     1 5.80000000D+00, 5.90000000D+00, 6.00000000D+00, 6.10000000D+00,
     1 6.20000000D+00, 6.30000000D+00, 6.40000000D+00, 6.50000000D+00,
     1 6.60000000D+00, 6.70000000D+00, 6.80000000D+00, 6.90000000D+00,
     1 7.00000000D+00, 7.20000000D+00, 7.40000000D+00, 7.60000000D+00,
     1 7.80000000D+00, 8.00000000D+00, 8.25000000D+00, 8.50000000D+00,
     1 9.00000000D+00, 9.50000000D+00, 1.00000000D+01/
      DATA (EKW(I), I = 1, 42) /       8.79797200D-01, 6.49071800D-01,
     1 4.73373000D-01, 3.37229300D-01, 2.30365900D-01, 1.45638600D-01,
     1 7.79739000D-02, 2.36643000D-02,-2.00556000D-02,-8.36422000D-02,
     1-1.24538500D-01,-1.50056200D-01,-1.64934200D-01,-1.72345900D-01,
     1-1.73962700D-01,-1.74451700D-01,-1.74474600D-01,-1.74459900D-01,
     1-1.74055800D-01,-1.72853700D-01,-1.68579900D-01,-1.62457000D-01,
     1-1.55067000D-01,-1.46849600D-01,-1.38131200D-01,-1.29156200D-01,
     1-1.20123300D-01,-1.11172500D-01,-1.02412700D-01,-9.39273000D-02,
     1-8.57810000D-02,-7.80164000D-02,-7.06700000D-02,-6.37641000D-02,
     1-5.73118000D-02,-5.13185000D-02,-4.57832000D-02,-4.07003000D-02,
     1-3.60578000D-02,-3.18402000D-02,-2.80272000D-02,-2.45978000D-02/
      DATA (EKW(I), I = 43, 85) /
     1-2.15297000D-02,-1.87967000D-02,-1.63689000D-02,-1.42247000D-02,
     1-1.23371000D-02,-1.06810000D-02,-9.23030000D-03,-7.96820000D-03,
     1-6.87030000D-03,-5.91780000D-03,-5.09230000D-03,-4.37820000D-03,
     1-3.76260000D-03,-3.23090000D-03,-2.77400000D-03,-2.38000000D-03,
     1-2.04230000D-03,-1.75210000D-03,-1.50300000D-03,-1.28990000D-03,
     1-1.10690000D-03,-9.49800000D-04,-8.15000000D-04,-7.00200000D-04,
     1-6.03000000D-04,-5.16200000D-04,-4.46600000D-04,-3.86400000D-04,
     1-3.32800000D-04,-2.90600000D-04,-2.46600000D-04,-2.15400000D-04,
     1-1.88900000D-04,-1.43400000D-04,-1.08600000D-04,-8.68000000D-05,
     1-6.82000000D-05,-5.28000000D-05,-4.04000000D-05,-3.14000000D-05,
     1-1.85000000D-05,-1.21000000D-05,-9.10000000D-06/
      DATA (DEKW(I), I = 1, 42) /     -5.30655810D+00,-4.00127680D+00,
     1-3.07706650D+00,-2.40187940D+00,-1.89616290D+00,-1.50962800D+00,
     1-1.20917840D+00,-9.72328800D-01,-7.83386600D-01,-5.07357200D-01,
     1-3.22463400D-01,-1.95624800D-01,-1.07126900D-01,-4.46910000D-02,
     1-2.06553000D-02,-4.17850000D-03, 7.10000000D-06, 3.25050000D-03,
     1 1.66583000D-02, 3.10060000D-02, 5.31071000D-02, 6.84124000D-02,
     1 7.86685000D-02, 8.51517000D-02, 8.87911000D-02, 9.02855000D-02,
     1 9.01292000D-02, 8.87163000D-02, 8.63479000D-02, 8.32555000D-02,
     1 7.96259000D-02, 7.55929000D-02, 7.12906000D-02, 6.68052000D-02,
     1 6.22279000D-02, 5.76293000D-02, 5.30699000D-02, 4.86005000D-02,
     1 4.42689000D-02, 4.01145000D-02, 3.61640000D-02, 3.24479000D-02/
       DATA (DEKW(I), I = 43, 85) /
     1 2.89726000D-02, 2.57562000D-02, 2.28079000D-02, 2.01113000D-02,
     1 1.76772000D-02, 1.54842000D-02, 1.35072000D-02, 1.17639000D-02,
     1 1.02208000D-02, 8.86280000D-03, 7.66860000D-03, 6.62540000D-03,
     1 5.71640000D-03, 4.92440000D-03, 4.23610000D-03, 3.63840000D-03,
     1 3.12560000D-03, 2.68280000D-03, 2.30080000D-03, 1.97300000D-03,
     1 1.69130000D-03, 1.44420000D-03, 1.23550000D-03, 1.05950000D-03,
     1 9.12000000D-04, 7.64500000D-04, 6.55700000D-04, 5.61100000D-04,
     1 4.77100000D-04, 4.12400000D-04, 3.46200000D-04, 2.95000000D-04,
     1 2.53200000D-04, 1.86000000D-04, 1.36600000D-04, 1.02700000D-04,
     1 7.73000000D-05, 5.56000000D-05, 3.62000000D-05, 2.48000000D-05,
     1 1.75000000D-05, 8.50000000D-06, 5.90000000D-06/
      END
C
C****
C
      SUBROUTINE MXLNEQ(A,NN,IDA,DETT,JRANK,EPS,IN,MMM)
C
C    PROGRAMMED BY R. HOTCHKISS, U. COMP. CTR., U. OF MINN., REVISED OCT.73
C    modified for the VAX by Bruce Garrett, Nov. 1980
C
C   Output from this subprogram is written to the file linked to FORTRAN unit
C   IUNIT.  The variable IUNIT is set in a parameter statement, and has a 
C   default value equal to 6.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MACHEP
      PARAMETER (IPRT = 6)
      DIMENSION A(IDA, 61),IN(IDA)
C
C     DATA ZERO,ONE/0.0D0,1.0D0/,MACHEP/1.D35/
C     REVISED BY N.ABUSALBI TO ALLOW A DET OF THE ORDER OF 10**37
C
      DATA ZERO,ONE/0.0D0,1.0D0/,MACHEP/1.D38/                          2/11/82
      DATA INTMX /400000/
C
C INITIATE SOME LOCAL AND OUTPUT VARIABLES
C
      JRANK = NN
      N = NN
      ID=IDA
      APS=EPS
      OPS = EPS
C
C error checking
C
      IF(N.GT.0 .AND. N.LE.ID) GO TO 1
      WRITE(IPRT,607)
      WRITE(IPRT,600) N,ID
      STOP 'MXLNEQ 1'
    1 IF(ID.GE.N .AND. ID.LT.INTMX) GO TO 2
      WRITE(IPRT,607)
      WRITE(IPRT,601) ID,N,INTMX
      STOP 'MXLNEQ 2'
    2 IF((ONE+OPS) .NE. ONE .AND. OPS.GE.ZERO) GO TO 3
      WRITE(IPRT,607)
      WRITE(IPRT,602) OPS
      STOP 'MXLNEQ 3'
    3 MOD = MMM
      MM = MMM
      NM = IABS(MM)
      IF(NM.LT.INTMX) GO TO 4
      WRITE(IPRT,607)
      WRITE(IPRT,603) MM,INTMX
      STOP 'MXLNEQ 4'
    4 NM = NN + NM
      NMIDA = NM*ID*2
      IF(NMIDA.LT.INTMX) GO TO 5
      WRITE(IPRT,607)
      WRITE(IPRT,604) NMIDA,INTMX
      STOP 'MXLNEQ 5'
    5 CONTINUE
C
C  end of error checking
C CONTINUE TO INITIALIZE
C
      K1=1
      NFLAG=0
      DET=ONE
C
C    MOD IS -1 DET ONLY (this version does not have this option)
C        = or > 0 FOR INV, DET AND 0 OR MORE SETS OF LIN EQNS
C        <-1 FOR LIN EQNS ONLY AND DET
C    MAIN GAUSS JORDAN LOOP BEGINS
C
      DO 75 K=1,N
C
C  SEARCH FOR LARGEST PIVOT CANDIDATE IN 
C  REMAINING LOWER RIGHT SQUARE MATRIX
C
      PIV=ZERO
      L=K
      DO 10 I=K,N
      P=ABS(A(I,K))
      IF(PIV.GE.P) GO TO 10
      PIV=P
      L=I
10    CONTINUE
C
C PIVOT WITH ABS VALUE PIV AND SUBSCRIPTS L AND M HAS BEEN FOUND
C
      PIVOT=A(L,K)
C
C CONTINUE IF PIV LARGER THAN USER EPS
C
      IF(PIV.GT.OPS) GO TO 20
      IF(EPS.EQ.OPS)JRANK=K-1
C
C EPS TEST FAILED, CHECK FOR ZERO PIVOT
C
      IF(PIV.GT.ZERO) GO TO 15
C
C PIVOT IS ZERO, TERMINATE PROGRAM UNLESS MOD=-1,IE, DET ONLY CASE
C      this version does not have the det only mode
C     IF(MOD.NE.(-1)) GO TO 11
C ZERO PIVOT MEANS ZERO DET AND EXIT IF DET ONLY MODE
C
C     DETT=0
C     RETURN
      WRITE(IPRT,607)
      WRITE(IPRT,605) K
      STOP 'MXLNEQ 6'
C
C ISSUE NON-FATAL MESSAGE, PIV .LE. EPS
C
   15 WRITE(IPRT,607)
      WRITE(IPRT,606) K,EPS
C
C SET OPS TO 0 SO SOLUTION MAY CONTINUE AFTER ERROR MESSAGE
C
      OPS=ZERO
      PIV=PIVOT
C
C CALCULATE DETERMINANT AND CHECK FOR OVERFLOW
C
C 20  DET=ONE
   20 DET=PIVOT*DET                                                   DWS11/6/84
      IF(ABS(DET).LT.MACHEP) GO TO 25
      WRITE(IPRT,607)
      WRITE(IPRT,608)K,DET
      STOP 'MXLNEQ 7'
C
C RESET LEADING ROW DO INDEX FOR DET ONLY AND LIN EQN ONLY CASE
C
   25 IF(MOD.LT.0)K1=K
C
C SAVE PIVOT INDEX
C
      IN(K)=L
C
C CHECK FOR ROW INTERCHANGE
C
      IF(L.EQ.K) GO TO 50
      DET=-DET
C
C INTERCHANGE ROW CONTAINING PIVOT AND CURRENT ROW
C ONLY PARTIAL ROWS NEED BE EXCHANGED FOR DET ONLY 
C OR LIN EQN ONLY SINCE LOWER LEFT PARTIALLY FORMED 
C TRIANGLE IS NOT NEEDED
C
      DO 30 J=K1,NM
      Z=A(L,J)
      A(L,J)=A(K,J)
30    A(K,J)=Z
C
C PIVOT ELEMENT IS NOW ON DIAGONAL
C SAVE DIVISION TIME BY USING RECIPROCAL OF PIVOT
C
50    PIVOT=ONE/PIVOT
C
C PRE-DIVIDE NECESSARY PORTION OF PIVOT ROW
C
      DO 55 J=K1,NM
55    A(K,J)=A(K,J)*PIVOT
C
C SET PIVOT ELEMENT TO ZERO SO MAIN REDUCTION 
C STEP DOESNT OPERATE ON PIVOT ROW
C
      A(K,K)=ZERO
C
C SWEEP THROUGH ALL OR PART OF MATRIX USING 
C KTH ROW, PIVOT ROW, TO REDUCE THE MATRIX
C
      DO 70 I=K1,N
      Z=A(I,K)
      IF(Z.EQ.ZERO) GO TO 70
C
C THIS CHECK NOT ONLY PREVENTS OPERATING ON PIVOT ROW BUT CATCHES
C  OTHER ZEROES IN PIVOT COLUMNS. THESE OTHER ZEROES WOULD LEAVE JTH
C  ROW UNCHANGED IN FOLLOWING LOOP SO CONSIDERABLE TIME MAY BE SAVED BY
C  SKIPPING OPERATION
C
      DO 65 J=K1,NM
65    A(I,J)=A(I,J)-Z*A(K,J)
C
C  THE INVERSE IS CREATED IN PLACE BY SUBSTITUTING AN IDENTITY MATRIX
C  COL BY COL, SINCE WE ARE SUBT. THE PIVOT ROW FROM OFF DIAGONAL 0
C  ELEMENTS AT THIS POINT, WE NOW PLACE -A(I,K)/A(K,K) AT THIS POINT IN
C  THE PIVOT COL
C
      A(I,K)=-Z*PIVOT
70    CONTINUE
C
C SIMILARLY DIVIDING PIVOT ROW BY THE PIVOT IS EQUIVALENT 
C TO PLACING ONE/A(K,K) AT THE PIVOT POINT FOR THE INVERSE
C
75    A(K,K)=PIVOT
      IF(N.EQ.1) GO TO 110
      IF(MOD.GE.0) GO TO 77
C
C BACK SUBSTITUTION FOR LIN EQN ONLY CASE
C
      K1=K1+1
      DO 125 K=K1,NM
      I=N
      DO 125 L=2,N
      I1=I
      I=I-1
      Z=ZERO
      DO 120 J=I1,N
120   Z=Z+A(I,J)*A(J,K)
125   A(I,K)=A(I,K)-Z
      GO TO 110
C
C FINAL REORDERING OF MATRIX
C
77    K=N
C
C SKIP LAST STEP SINCE NO INTERCHANGE COULD OCCUR THERE
C
      DO 105 J=2,N
C
C PERFORM INTERCHANGES IN EXACT REVERSE ORDER OF PREVIOUS EXECUTION
C
      K=K-1
C
C ROW INTERCHANGE DURING INVERSION IMPLIES COL INTERCHANGE HERE
C
      M=IN(K)
      IF(M.EQ.K) GO TO 105
C
C COL INTERCHANGE
C
      DO 85 I=1,N
      Z=A(I,K)
      A(I,K)=A(I,M)
85    A(I,M)=Z
105   CONTINUE
110   DETT=DET
      RETURN
  600 FORMAT(/,2X,T5,'arg 2, n = ',I10,T40,'id = ',I10,
     1       /,2X,T5,'n must be .ge. 1 and .le. id')
  601 FORMAT(/,2X,T5,'arg 3, id = ',I10,T40,'n = ',I10,
     1       /,2X,T5,'id must be .ge. N and .le. ',I10)
  602 FORMAT(/,2X,T5,'arg 6, eps = ',1PE13.5,
     1       /,2X,T5,'eps must be .ge. zero and finite')
  603 FORMAT(/,2X,T5,'arg 8, m = ',I10,
     1       /,2X,T5,'abs(m) must be .le. ',I10 )
  604 FORMAT(/,2X,T5,'size = ',I10,
     1       /,2X,T5,'size = id*(n+abs(m))*2 must be .le. ',I10)
  605 FORMAT(/,2X,T5,'k = ',I10
     1       /,2X,T5,'at step k a gauss-jordan pivot value was zero')
  606 FORMAT(/,2X,T5,'k = ',I10,T40,'eps = ',1PE13.5,
     1      /,2X,T5,'at step k a gauss-jordan pivot value was .le. eps')
  607 FORMAT(/,2X,T5,'*** mxlneq ***')
  608 FORMAT(/,2X,T5,'k = ',I10,T40,'det = ',1PE13.5,
     1       /,2X,T5,'at step k det is too large')
      END
