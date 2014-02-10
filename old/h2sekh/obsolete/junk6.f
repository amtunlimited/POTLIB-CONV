c****************************************************         
c System: H2Se
c Common Name: h2sekh,
c Functional Form: Expanded into 2 bond length and 1 bond angle               
c Internal Coords: 2 bond length and 1 bond angle
c Features: Calibrated for the equilibrium state          
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
      COMMON /DATACM/ de,a,t1,t2,                                                          
     +            alpha,beta,gamma,delta,epsilon,rkappa,                        
     +            rlambda,rmu,rnu,rho,                                          
     +            re,phie               
c****************************************************         
c any data or labeled commons for PTPACM or POT               
c****************************************************         
      REF(1)= 'Potential Energy for H2Se'
      REF(2)= ' '
      REF(3)= ' '
      REF(4)= ' '
      REF(5)= ' '
      IRCTNT = 3
      INDEXES(1) = 1
      INDEXES(2) = 34                                   
      INDEXES(3) = 1
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
      COMMON /DATACM/ de,a,t1,t2,                                                          
     +            alpha,beta,gamma,delta,epsilon,rkappa,                        
     +            rlambda,rmu,rnu,rho,                                          
     +            re,phie          
c****************************************************         
c any labeled common statements for PTPACM or POT             
c****************************************************         
      CALL CARTOU
      CALL CARTTOR
      
      cosphi=((R(1)**2 + R(2)**2 - R(3)**2) / (2.d0 * R(1) * R(2)))                                                 
      phi=acos(cosphi)  
c****************************************************         
c code for V,dVdr as functions of int coord R                 
c or for polyatomics as functions of cartesian coord          
c****************************************************
c START OF ORIGINAL CODE
                           
                                                       
                                                                                
c                                                                               
c  The potential is expanded in terms of y1, y2, theta in                       
c  the paper by Kauppi and Halonen, as below.                                   
c                                                                               
	y1=1.0d00-dexp(-a*(R(1)-re))                                                     
	y2=1.0d00-dexp(-a*(R(2)-re))                                                     
	theta=phi-phie                                                                 
c                                                                               
c  for efficiency, next we use the greek letter parameters                      
C  as "force constants" and calculate the potential energy.                     
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
c END OF ORIGINAL CODE
         
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
c     PARAMETER (ia=1,ib=1,ic=1,id=1)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF
      COMMON /DATACM/ de,a,t1,t2,                                                          
     +            alpha,beta,gamma,delta,epsilon,rkappa,                        
     +            rlambda,rmu,rnu,rho,                                          
     +            re,phie          
c****************************************************         
c all other labeled commons for PREPOT or POT                 
c****************************************************         
      DATA NDER,ANUZERO,NFLAG,NULBL/0,0.0,1,1,15*0,6,0,0,25*0/
      DATA NATOMS,ICARTR,MDER,MSURF/3,3,0,0/
      DATA NASURF/1,35*0/
c****************************************************         
c all other data statements                                   
c****************************************************
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
