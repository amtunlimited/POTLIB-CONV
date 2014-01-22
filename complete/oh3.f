c**************************************************** 
c POTLIB 2001: A potential energy surface library for chemical systems
c R. J. Duchovic, Y. L. Volobuev, G. C. Lynch, D. G. Truhlar, T. C. Allison,
c   A. F. Wagner, B. C. Garrett, and J. C. Corchado
c Computer Physics Communications 144, 169-187 (2002), 156, 319-322(E) (2004)
c
c System: OH3
c Name: OH + H2              
c Int Coords:
c Special Features: One derivitive

c Original code written by:
c       G. C. Schatz and H. Elgersma
C       Chem. Phys. Lett. 73, 21 (1980).
c Transcribed to the POTLIB Standard by:
c       Aaron Tagliaboschi <aaron.tagliaboschi@gmail.com>
c       Western Kentucky University, Department of Chemestry
c         Dr. Jeremy B. Maddox
c****************************************************         
      SUBROUTINE PREPOT
      IMPLICIT REAL*8 (A-H,O-Z)
      character text1,text2,text3,text4,text5
      CHARACTER*75 REF(5)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF
      COMMON /CONCOM/ XDE(3), XBETA(3), XRE(3), SATO, GAM(3), REOH, 
     *                REHH, CON(7), ALP(4), CLAM(4), ACON(2)
      DIMENSION DE(3), BETA(3), RE(3), Z(3), ZPO(3), OP3Z(3), ZP3(3), 
     *          TZP3(3), TOP3Z(3), DO4Z(3), B(3), X(3), COUL(3), EXCH(3) 
C
      PARAMETER (R2 = 1.41421356D0)           
c****************************************************         
c any data or labeled commons for PTPACM or POT               
c****************************************************         
      REF(1)= "G. C. Schatz and H. Elgersma"
      REF(2)= "Chem. Phys. Lett. 73, 21 (1980)."
      REF(3)= text3
      REF(4)= text4
      REF(5)= text5
      IRCTNT = 3
      INDEXES(1) = 8
      INDEXES(2) = 1
      INDEXES(3) = 1
      INDEXES(4) = 1
      CALL POTINFO
      CALL ANCVRT
c****************************************************         
c any assignments for labeled common share with pot           
c****************************************************       
C
C   Initialize the flag for the potential calculation
C
         NDER = 1
C
C   Set up the values of DE, BETA, and RE for the three-body LEPS potential.
C
         DE(1)   = XDE(1)
         BETA(1) = XBETA(1)
         RE(1)   = XRE(1)
         DE(2)   = XDE(1)
         BETA(2) = XBETA(1)
         RE(2)   = XRE(1)
         DE(3)   = XDE(3)
         BETA(3) = XBETA(3)
         RE(3)   = XRE(3)
C
C   Echo the potential energy surface parameters to unit 6.
C
      WRITE (6,1000) XDE, XBETA, XRE, SATO
      WRITE (6,1200) (GAM(I),I=1,3),REOH,REHH
      WRITE (6,1300) CON
      WRITE (6,1400) ALP,CLAM,ACON
C
      DO 10 I = 1, 3
         Z(I) = SATO
C   Compute useful constants.
         ZPO(I)  = 1.0D0+Z(I)
         OP3Z(I)  = 1.0D0+3.0D0*Z(I)
         TOP3Z(I) = 2.0D0*OP3Z(I)
         ZP3(I)   = Z(I)+3.0D0
         TZP3(I)  = 2.0D0*ZP3(I)
         DO4Z(I)  = DE(I)/4.0D0/ZPO(I)
         B(I)     = BETA(I)*DO4Z(I)*2.0D0
   10 CONTINUE
C
1000  FORMAT (/,2X,T5,'OH + H2 potential energy function',
     *        //, 2X, T5, 'Potential energy surface parameters ',
     *                    'in hartree atomic units:',
     *        /,2X,T5,'Morse and LEPS parameters:',
     *        /,2X,T5,'Dissociation energies:', T31,1P,3E13.6,
     *        /,2X,T5,'Morse betas:', T31,1P,3E13.6,
     *        /,2X,T5,'Equilibrium bond lengths:',T31,1P,3E13.6,
     *        /,2X,T5,'Sato parameter:', T31, 1PE13.6)
1200  FORMAT (/,2X,T5,'GAM:',T20,1P,3E13.6,
     *        /,2X,T5,'REOH, REHH:',T20,1P,2E13.6)
1300  FORMAT (/,2X,T5,'CON:',T20,1P,4E13.6,(/,2X,T20,1P,4E13.6))
1400  FORMAT (/,2X,T5,'ALP:',T20,1P,4E13.6,
     *        /,2X,T5,'CLAM:',T20,1P,4E13.6, 
     *        /,2X,T5,'ACON:',T20,1P,2E13.6)
C
      RETURN
      END

      SUBROUTINE POT
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF            
      COMMON /USROCM/ PENGYGS,PENGYES(ISURF),PENGYIJ(JSURF),  
     X                DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     X                DIJCART(NATOM,3,JSURF)                  
      COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)          
      COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)       
      COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
      COMMON /POTCM/ RS(6), VTOT, DVDR(6)
      
      VMOR(D,B,T,RR) = D*(1.0D0-EXP(-B*(RR-T)))**2
      DVMOR(D,B,T,RR) = 2.0D0*B*D*(1.0D0-EXP(-B*(RR-T)))*EXP(-B*(RR-T))
             
c****************************************************         
c any labeled common statements for PTPACM or POT             
c****************************************************         
      CALL CARTOU
      CALL CARTTOR
c****************************************************         
c code for V,dVdr as functions of int coord R                 
c or for polyatomics as functions of cartesian coord          
c****************************************************

c Start of original code

C
C   Zero the array which will contain the derivatives of the energy 
C   with respect to the internal coordinates.
C
         IF (NDER .EQ. 1) THEN
             DO 10 I = 1, 6
                   DEGSDR(I) = 0.0D0
10           CONTINUE
         ENDIF
C
C   Calculate the Morse part of the potential energy.
      VTOT = VMOR(DE(1),BETA(1),RE(1),RS(1))+VMOR(DE(2),BETA(2),RE(2),
     *  RS(4))+VMOR(DE(2),BETA(2),RE(2),R(5))
C   Calculate the derivatives of the Morse part of the potential energy.
         IF (NDER .EQ. 1) THEN
             DVDR(1) = DVDR(1)+DVMOR(DE(1),BETA(1),RE(1),RS(1))
             DVDR(4) = DVDR(4)+DVMOR(DE(2),BETA(2),RE(2),RS(4))
             DVDR(5) = DVDR(5)+DVMOR(DE(2),BETA(2),RE(2),RS(5))
         ENDIF
C   Initialize the coordinates for the three-body LEPS part of the potential.
      R(1) = RS(2)
      R(2) = RS(3)
      R(3) = RS(6)
C
C   Calculate the three-body LEPS portion of the potential and update the 
C   energy term.
      CALL POTEN
      VTOT = VTOT+ENGYGS
C   Update the array containing the derivatives with the LEPS derivatives.
         IF (NDER .EQ. 1) THEN
             DVDR(2) = DVDR(2)+DEGSDR(1)
             DVDR(3) = DVDR(3)+DEGSDR(2)
             DVDR(6) = DVDR(6)+DEGSDR(3)
         ENDIF
C   Initialize the coordinates for the H2O part of the potential for H1 and H2.
      R(1) = RS(1)
      R(2) = RS(2)
      R(3) = RS(4)
C   Calculate the H2O part of the potential and update the energy term.
      CALL VH2O
      VTOT = VTOT+ENGYGS
C   Update the array containing the derivatives with the H2O derivatives.
         IF (NDER .EQ. 1) THEN
             DVDR(1) = DVDR(1)+DEGSDR(1)
             DVDR(2) = DVDR(2)+DEGSDR(2)
             DVDR(4) = DVDR(4)+DEGSDR(3)
         ENDIF
C   Initialize the coordinates for the H2O part of the potential for H1 and H3.
      R(1) = RS(1)
      R(2) = RS(3)
      R(3) = RS(5)
C   Calculate the H2O part of the potential and update the energy term.
      CALL VH2O
      VTOT = VTOT+ENGYGS
C   Update the array containing the derivatives with the H2O derivatives.
         IF (NDER .EQ. 1) THEN
             DVDR(1) = DVDR(1)+DEGSDR(1)
             DVDR(3) = DVDR(3)+DEGSDR(2)
             DVDR(5) = DVDR(5)+DEGSDR(3)
         ENDIF
C   Initialize the coordinates for the four-body part of the potential.
      R(1) = RS(2)
      R(2) = RS(3)
      R(3) = RS(4)
      R(4) = RS(5)
C   Calculate the four-body part of the potential and update the energy term.
      CALL V4POT
      VTOT = VTOT+ENGYGS
C   Update the array containing the derivatives.
         IF (NDER .EQ. 1) THEN
             DVDR(2) = DVDR(2)+DEGSDR(1)
             DVDR(3) = DVDR(3)+DEGSDR(2)
             DVDR(4) = DVDR(4)+DEGSDR(3)
             DVDR(5) = DVDR(5)+DEGSDR(4)
         ENDIF
C
C   Adjust the potential to the correct zero, which corresponds to
C   OH(R=RE) and H2(R=RE) at infinity.
C
      VTOT = VTOT-2.0D0*DE(2)+DE(3)

c End of original code        
      CALL EUNITZERO
      IF(NDER.NE.0) THEN
         CALL RTOCART
         IF(NFLAG(1)+NFLAG(2).NE.0) CALL DEDCOU
      ENDIF
      RETURN
      END
      
      BLOCK DATA PTPACM
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
      PARAMETER (ia=1,ib=1,ic=1,id=1)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF  
      COMMON /CONCOM/ XDE(3),XBETA(3),XRE(3),SATO,GAM(3),REOH,REHH,
     *                CON(7),ALP(4),CLAM(4),ACON(2)          
c****************************************************         
c all other labeled commons for PREPOT or POT                 
c****************************************************         
      DATA NDER,ANUZERO,NFLAG,NULBL/1,0.0,1,1,15*0,6,0,0,25*0/
      DATA NATOMS,ICARTR,MDER,MSURF/4,2,1,0/
      DATA NASURF/1,35*0/
c****************************************************         
c all other data statements                                   
c**************************************************** 
      DATA XDE / 0.148201D0, 0.0275690D0, 0.151548D0 /  
      DATA XBETA / 1.260580D0, 0.924180D0, 1.068620D0 /  
      DATA XRE / 1.863300D0, 2.907700D0, 1.428600D0 /
      DATA SATO / 0.10D0 /
      DATA GAM / 2.399700D0, 1.058350D0, 2.399700D0 /
      DATA REOH / 1.808090D0 /
      DATA REHH / 2.861590D0 /
      DATA CON / -.0015920D0, 0.026963D0, 0.0014689D0, 0.080011D0,
     *           0.085816D0, -0.063179D0, 0.101380D0 /
      DATA ALP / 4.773D0, 7.14D0, 2.938D0, 5.28D0 /
      DATA CLAM / 0.10D0, 0.10D0, 0.20D0, 0.03D0 /
      DATA ACON / 0.10D0, 0.009D0 /        
      END
      
c Extra subroutines from the original code

C*****
      SUBROUTINE V4POT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      COMMON /CONCOM/ DUM(22),A(4),C(4),COF(2)
      COMMON /PT1CM/  R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)  
      COMMON /POTCCM/ NSURF, NDER, NDUM(8)
C
      T1 = EXP(-C(1)*(R(1)-A(1))**2-C(1)*(R(2)-A(1))**2-C(3)*(R(3)-A(3))
     *   **2-C(3)*(R(4)-A(3))**2)*COF(1)
      T2 = EXP(-C(2)*(R(1)-A(2))**2-C(2)*(R(2)-A(2))**2-C(4)*(R(3)-A(4))
     *   **2-C(4)*(R(4)-A(4))**2)*COF(2)
      ENGYGS = T1+T2
         IF (NDER .EQ. 1) THEN
            DEGSDR(1) = -2.0D0*(T1*C(1)*(R(1)-A(1))+T2*C(2)*(R(1)-A(2)))
            DEGSDR(2) = -2.0D0*(T1*C(1)*(R(2)-A(1))+T2*C(2)*(R(2)-A(2)))
            DEGSDR(3) = -2.0D0*(T1*C(3)*(R(3)-A(3))+T2*C(4)*(R(3)-A(4)))
            DEGSDR(4) = -2.0D0*(T1*C(3)*(R(4)-A(3))+T2*C(4)*(R(4)-A(4)))
         ENDIF
C
      RETURN
      END
C*****
C
      SUBROUTINE VH2O
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      DIMENSION S(3),Q(3),DQ(3),X(3),DP(3)
      COMMON /CONCOM/ DUM1(10),GAM(3),REOH,REHH,C(7),DUM2(10)
      COMMON /PT1CM/  R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)
      COMMON /POTCCM/ NSURF, NDER, NDUM(8)
C
C     XMAX1 FOR WHEN TANH SET=1.0 ON VAX
C     XMAX2 FOR PREVENTING OVERFLOWS ON VAX
C
      DATA XMAX1,XMAX2 / 15.0D0,43.0D0 /
C
C     STATEMENT FUNCTION
C
      TANLG(XX) = 2.0D0/(1.0D0+EXP(2.0D0*XX))
      S(1) = R(1)-REOH
      S(3) = R(2)-REOH
      S(2) = R(3)-REHH
      DO 30 I = 1, 3
         X(I) = 0.5D0*GAM(I)*S(I)
         Q(I) = 1.0D0-TANH(X(I))
         IF (X(I).LT.XMAX1) GO TO 20
         IF (X(I).LT.XMAX2) GO TO 10
         Q(I) = 0.0D0
         DQ(I) = 0.0D0
         GO TO 30
   10    Q(I) = TANLG(X(I))
   20    DQ(I) = -0.5D0*GAM(I)/COSH(X(I))**2
   30 CONTINUE
      P = C(1)+C(2)*(S(1)+S(3))+C(3)*S(2)+0.5D0*C(4)*(S(1)*S(1)+S(3)*S(
     *    3))+0.5D0*C(5)*S(2)*S(2)+C(6)*S(2)*(S(1)+S(3))+C(7)*S(1)*S(3)
      ENGYGS = Q(1)*Q(2)*Q(3)*P
      IF (NDER .EQ. 1) THEN 
          DP(1) = C(2)+C(4)*S(1)+C(6)*S(2)+C(7)*S(3)
          DP(2) = C(3)+C(5)*S(2)+C(6)*(S(1)+S(3))
          DP(3) = C(2)+C(4)*S(3)+C(6)*S(2)+C(7)*S(1)
          DO 40 I = 1, 3
                TRM1 = 0.0D0
                IF (Q(I).EQ.0.0D0) GO TO 40
                TRM1 = DQ(I)/Q(I)
   40           DEGSDR(I) = E*(TRM1+(DP(I)/P))
          TEMP = DEGSDR(2)
          DEGSDR(2) = DEGSDR(3)
          DEGSDR(3) = TEMP
      ENDIF
C
      RETURN
      END
      
      SUBROUTINE POTEN
      DIMENSION DE(3), BETA(3), RE(3), Z(3), ZPO(3), OP3Z(3), ZP3(3), 
     *          TZP3(3), TOP3Z(3), DO4Z(3), B(3), X(3), COUL(3), 
     *          EXCH(3)
C
C   Initialize the variable used for storing the energy.
C
      ENGYGS = 0.D0
C
      DO 20 I = 1, 3
         X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))
         COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)
         EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)
         ENGYGS  = ENGYGS + COUL(I)
   20 CONTINUE
      RAD = SQRT((EXCH(1)-EXCH(2))**2+(EXCH(2)-EXCH(3))**2+(EXCH(3)-EXCH
     *   (1))**2)
      ENGYGS = ENGYGS - RAD/R2 
C
C   Compute the derivatives of the energy with respect to the internal 
C   coordinates.
         IF (NDER .EQ. 1) THEN
             S = EXCH(1) + EXCH(2) + EXCH(3)
             DO 30 I = 1, 3
                   DEDR(I) = B(I)*X(I)*((3.0D0*EXCH(I)-S)/R2*
     *                       (OP3Z(I)*X(I)-ZP3(I))/RAD-
     *                       ZP3(I)*X(I)+OP3Z(I))
30           CONTINUE
         ENDIF
C
      RETURN
      END
