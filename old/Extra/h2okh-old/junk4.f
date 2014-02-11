c****************************************************         
c chemical system: H+OH on the ground state water surface
c common name:  of potential: kauppi and halonen H2O PES
c functional form of PES: nothing special
c int coord: valence (2 bond distances, 1 bond angle
c special features: design to represent bottom of H2O PES, 
c                   not H+O2 asymptot          
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
      REF(1)= 'Kauppi and Halonen H2O PES for near equil region'
      REF(2)= 'J Phys Chem 94, 5779 (1990)'
      REF(3)= ' '
      REF(4)= ' '
      REF(5)= ' '
      IRCTNT = 2
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
      common /data1/de,a,t1,t2,alpha,beta,gamma,delta,epsilon,rkappa,
     &              rlambda,rmu,rnu,rho,re,phie                      
c****************************************************         
c any labeled common statements for PTPACM or POT             
c****************************************************         
      CALL CARTOU
      CALL CARTTOR
c      write(6,*) 'POT R:'
c      write(6,101) (r(i),i=1,3)
c101   format(8f10.5)
c****************************************************         
c code for V,dVdr as functions of int coord R                 
c or for polyatomics as functions of cartesian coord          
c****************************************************         
c
c construct 2 OH bond lengths(au), HOH angle(radians) used by PES
c
      r1 = r(1)
      r2 = r(2)
      cosphi = (r1**2+r2**2-r(3)**2)/2.d0/r1/r2
      phi = acos(cosphi)
c      write(6,*) 'POT R,phi:'
c      write(6,101) r1,r2,phi
c
c calculate the energy in atomic units as expansion in terms of y1,y2,theta
c
      y1=1.d0-dexp(-a*(r1-re))                                                
      y2=1.d0-dexp(-a*(r2-re))                                                
      theta=phi-phie                                                          
      y1sq=y1*y1                                                              
      y2sq=y2*y2                                                              
      y1y2=y1*y2                                                              
      tsq=theta*theta                                                         
      engygs = (y1sq+y2sq)*(rmu*tsq+rnu*theta+de)                                 
     &       + (y1+y2)*(rkappa*theta+rlambda*tsq+beta*y1y2)                     
     &       + y1y2*(alpha+rho*theta)                                           
     &       + tsq*(gamma+delta*theta+epsilon*tsq)                              
     &       + t2*(y1sq*y1sq+y2sq*y2sq)                                         
c      write(6,*) 'POT energy'
c      write(6,101) engygs
c      return
      CALL EUNITZERO
      IF(NDER.NE.0) THEN
         CALL RTOCART
         IF(NFLAG(1)+NFLAG(2).NE.0) CALL DEDCOU
      ENDIF
      RETURN
      END
      subroutine surf(ndof,cartx,v,n)                                           
c     Kauppi-Halonen surface for H2O                                            
      implicit double precision(a-h,o-z)                                        
      dimension cartx(ndof+3,*),v(*)                                            
      common /data1/de,a,t1,t2,alpha,beta,gamma,delta,epsilon,rkappa,
     &              rlambda,rmu,rnu,rho,re,phie
      do 8888 i=1,n                                                             
       c3mc9=cartx(3,i)-cartx(9,i)                                              
       c2mc8=cartx(2,i)-cartx(8,i)                                              
       c1mc7=cartx(1,i)-cartx(7,i)                                              
       c6mc9=cartx(6,i)-cartx(9,i)                                              
       c5mc8=cartx(5,i)-cartx(8,i)                                              
       c4mc7=cartx(4,i)-cartx(7,i)                                              
       r1=dsqrt(c3mc9*c3mc9+c2mc8*c2mc8+c1mc7*c1mc7)                            
       r2=dsqrt(c6mc9*c6mc9+c5mc8*c5mc8+c4mc7*c4mc7)                            
       r1dotr2=(c1mc7*c4mc7+c2mc8*c5mc8+c3mc9*c6mc9)                            
       cosphi=r1dotr2/(r1*r2)                                                   
       phi=acos(cosphi)                                                         
c       write(6,*) 'SURF R,phi:'
c       write(6,101) r1,r2,phi
c101    format(8f10.5)
c      input are the two o-h bond lengths and the h-o-h angle; r1,r2,theta.     
c      all calculations in atomic units; angles in radians.                     
c  the potential is expanded in terms of y1,y2,theta in the paper by kauppi     
c  and halonen, as below.                                                       
        y1=1.d0-dexp(-a*(r1-re))                                                
        y2=1.d0-dexp(-a*(r2-re))                                                
        theta=phi-phie                                                          
c                                                                               
        y1sq=y1*y1                                                              
        y2sq=y2*y2                                                              
        y1y2=y1*y2                                                              
        tsq=theta*theta                                                         
c                                                                               
      v(i) = (y1sq+y2sq)*(rmu*tsq+rnu*theta+de)                                 
     &       + (y1+y2)*(rkappa*theta+rlambda*tsq+beta*y1y2)                     
     &       + y1y2*(alpha+rho*theta)                                           
     &       + tsq*(gamma+delta*theta+epsilon*tsq)                              
     &       + t2*(y1sq*y1sq+y2sq*y2sq)                                         
8888  continue         
c      write(6,*) 'SURF energy'
c      write(6,101) v(1)
      return                                                                    
      end                                                                       
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
     &              rlambda,rmu,rnu,rho,re,phie                      
      DATA NDER,ANUZERO,NFLAG,NULBL/0,0.0,1,1,15*0,6,0,0,25*0/
      DATA NATOMS,ICARTR,MDER,MSURF/3,3,0,0/
      DATA NASURF/1,35*0/
c****************************************************         
c all other data statements                                   
c****************************************************         
c     energies in atomic units; angles in radians.                     
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
