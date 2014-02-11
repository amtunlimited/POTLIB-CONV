c****************************************************         
c System: H2O
c Common Name: h2okh
c Functional Form: Expansion in 2 Bond Distances, 1 Bond Angle
c Internal Coords: 2 Bond Distances, 1 Bond Angle
c Special Features: Designed for near-equalibrium vibrational levels
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
      COMMON /DATACM/ de,a,t1,t2,alpha,beta,gamma,delta,epsilon,rkappa,              
     &            rlambda,rmu,rnu,rho,re,phie           
c****************************************************         
c any data or labeled commons for PTPACM or POT               
c****************************************************         
      REF(1)= 'H20 Potencial from Kauppi-Halonen'
      REF(2)= ' '
      REF(3)= ' '
      REF(4)= ' '
      REF(5)= ' '
      IRCTNT = 3
      INDEXES(1) = 1
      INDEXES(2) = 8                                           
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
      COMMON /DATACM/ de,a,t1,t2,alpha,beta,gamma,delta,epsilon,rkappa,              
     &            rlambda,rmu,rnu,rho,re,phie       
c****************************************************         
c any labeled common statements for PTPACM or POT             
c****************************************************         
      CALL CARTOU
      CALL CARTTOR
                                  
       cosphi=(R(1)**2 + R(2)**2 -R(3)**2) / (2.d0 * R(1) * R(2))                                                 
       phi=acos(cosphi)                                                         
c****************************************************         
c code for V,dVdr as functions of int coord R                 
c or for polyatomics as functions of cartesian coord          
c****************************************************  

c     Kauppi-Halonen surface for H2O 
                         
        y1=1.d0-dexp(-a*(R(1)-re))                                                
        y2=1.d0-dexp(-a*(R(2)-re))                                                
        theta=phi-phie                                                          
c                                                                               
        y1sq=y1*y1                                                              
        y2sq=y2*y2                                                              
        y1y2=y1*y2                                                              
        tsq=theta*theta                                                         
c                                                                               
      ENGYGS = (y1sq+y2sq)*(rmu*tsq+rnu*theta+de)                                 
     &       + (y1+y2)*(rkappa*theta+rlambda*tsq+beta*y1y2)                     
     &       + y1y2*(alpha+rho*theta)                                           
     &       + tsq*(gamma+delta*theta+epsilon*tsq)                              
     &       + t2*(y1sq*y1sq+y2sq*y2sq)                                         
                                                                                         
       
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
      DATA NDER,ANUZERO,NFLAG,NULBL/0,0.0,1,1,15*0,6,0,0,25*0/
      DATA NATOMS,ICARTR,MDER,MSURF/3,3,0,0/
      DATA NASURF/1,35*0/
c****************************************************         
c all other data statements                                   
c****************************************************
	COMMON /DATACM/ de,a,t1,t2,alpha,beta,gamma,delta,epsilon,rkappa,              
     &            rlambda,rmu,rnu,rho,re,phie
     
       data      de,a,t1,t2,alpha,beta,gamma,delta,epsilon,rkappa,              
     &            rlambda,rmu,rnu,rho,re,phie/                                  
     &    0.195286504346741885d0   , 1.17789560286810513d0   ,                  
     &    0.000000000000000000d+00 , 0.127759668062261463d-01 ,                 
     &   -0.394936461650084985d-02 ,-0.609580676207029798d-02 ,                 
     &    0.822960364451488097d-01 ,-0.249937976602641748d-01 ,                 
     &   -0.104937253831824621d-01 , 0.375501181199095618d-01 ,                 
     &   -0.149932556445900381d-01 ,-0.995022445892074572d-02 ,                 
     &    0.187750590599547809d-01 ,-0.146428799417284827d-01 ,                 
     &    1.80941259980574709d0   ,  1.82404363854352902d0   /                           
      END
