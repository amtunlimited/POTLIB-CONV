c****************************************************         
C System:         Modified Kauppi-Halonen surface for H2O  
C Common Name:    h20mkh
C Reference:      E. Kauppi and L. Halonen, J.Phys.Chem. 94, 5779 (1990)
c Functionalorm:  Expansion in MO for distances and bond angle
c Definitions:    Bond distances and bond angles
c Notes:          Primarily calibrated to near equilibrium
c****************************************************         
      SUBROUTINE PREPOT
      IMPLICIT REAL*8 (A-H,O-Z)
      character text1,text2,text3,text4,text5
      CHARACTER*75 REF(5)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+2)/2)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF            
c****************************************************         
c any data or labeled commons for PTPACM or POT               
c****************************************************         
      REF(1)= 'H2Se by Kauppi and Halonen'
      REF(2)= 'E. Kauppi and L. Halonen, J.Phys.Chem. 94, 5779 (1990)'
      REF(3)= 'Primarily calibrated by vibrational levels of H2Se'
      REF(4)= ' '
      REF(5)= ' '
      IRCTNT = 2
      INDEXES(1) = 1
      INDEXES(2) = 34
      INDEXES(NATOMS) = 1
      CALL POTINFO
      CALL ANCVRT
c****************************************************         
c any assignments for labeled common share with pot           
c****************************************************         
      RETURN
      END
      SUBROUTINE POT
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+2)/2)
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
c****************************************************         
c any labeled common statements for PTPACM or POT             
c****************************************************         
      common /data1/de,a,t1,t2,alpha,beta,gamma,delta,epsilon,rkappa,
     +              rlambda,rmu,rnu,rho,re,phie
      CALL CARTOU
      CALL CARTTOR
c****************************************************         
c code for V,dVdr as functions of int coord R                 
c or for polyatomics as functions of cartesian coord          
c****************************************************         
      r1 = r(1)
      r2 = r(2)
      cosphi = (r1**2+r2**2-r(3)**2)/2.d0/r1/r2
      phi=acos(cosphi)
c
c  The potential is expanded in terms of y1, y2, theta in
c  the paper by Kauppi and Halonen, as below.
c
      y1=1.0d00-dexp(-a*(r1-re))
      y2=1.0d00-dexp(-a*(r2-re))
      theta=phi-phie
c
c  for efficiency, next we use the greek letter parameters
c  as "force constants" and calculate the potential energy.
c
      y1sq=y1*y1
      y2sq=y2*y2
      y1y2=y1*y2
      tsq=theta*theta
c
      ENGYGS = (y1sq+y2sq)*(rmu*tsq+rnu*theta+de)       
     +       + (y1+y2)*(rkappa*theta+rlambda*tsq+beta*y1y2)
     +       + y1y2*(alpha+rho*theta)
     +       + tsq*(gamma+delta*theta+epsilon*tsq)
     +       + t2*(y1sq*y1sq+y2sq*y2sq)
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
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF            
c****************************************************         
c all other labeled commons for PREPOT or POT                 
c****************************************************         
      common /data1/de,a,t1,t2,alpha,beta,gamma,delta,epsilon,rkappa,
     +              rlambda,rmu,rnu,rho,re,phie
      DATA NDER,ANUZERO,NFLAG,NULBL/0,0.0,1,1,15*0,6,0,0,25*0/
      DATA NATOMS,ICARTR,MDER,MSURF/3,3,0,0/
      DATA NASURF/1,35*0/
c****************************************************         
c all other data statements                                   
c****************************************************         
c
c these are the potential parameters: t=theta, r=r, rp=r' (see Kauppi).
c
	data      de,a,t1,t2,
     +            alpha,beta,gamma,delta,epsilon,rkappa,
     +            rlambda,rmu,rnu,rho,
     +            re,phie/
     +    0.144760654398338406     ,       !de
     +    0.881450520988266817     ,       !a
     +    0.000000000000000000e+00 ,       !t1
     +    0.110900893492383024e-01 ,       !t2
     +   -0.179227352954911988e-02 ,       !alpha
     +   -0.896136764774559938e-03 ,       !beta
     +    0.802328447657812416e-01 ,       !gamma
     +   -0.615096649861561571e-02 ,       !delta
     +    0.000000000000000000e+00 ,       !epsilon
     +    0.110148223195307097e-01 ,       !rkappa
     +   -0.188996679643331431e-01 ,       !rlambda
     +   -0.250743588903988915e-01 ,       !rmu
     +    0.270905047360925745e-02 ,       !rnu
     +   -0.274545233066495094e-01 ,       !rho
     +    2.75900001813446805     ,        !re
     +    1.58074469820494201     /        !thetae
      END
