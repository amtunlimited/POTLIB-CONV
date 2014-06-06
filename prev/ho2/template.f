c***********************************************************************
c POTLIB 2001: A potential energy surface library for chemical systems
c R. J. Duchovic, Y. L. Volobuev, G. C. Lynch, D. G. Truhlar, T. C. Allison,
c   A. F. Wagner, B. C. Garrett, and J. C. Corchado
c Computer Physics Communications 144, 169-187 (2002), 156, 319-322(E) (2004)
c
c System:
c Name:               
c Int Coords:
c Special Features:
c
c Original code written by:
c       $REFERENCE$
c Transcribed to the POTLIB Standard by:
c       Aaron Tagliaboschi <aaron.tagliaboschi@gmail.com>
c       Western Kentucky University, Department of Chemestry
c         Dr. Jeremy B. Maddox
c***********************************************************************
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
c     Extra common blocks

      REF(1)= text1
      REF(2)= text2
      REF(3)= text3
      REF(4)= text4
      REF(5)= text5
      IRCTNT = aa
      INDEXES(1) = bb1
c     ...
      INDEXES(NATOMS) = bbnatoms
      CALL POTINFO
      CALL ANCVRT
c     Initializations and onetime calculations
         
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
c     Extra common blocks  
       
      CALL CARTOU
      CALL CARTTOR
c     Start of original code

c     End of original code        
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
c     Extra common blocks
        
      DATA NDER,ANUZERO,NFLAG,NULBL/0,0.0,1,1,15*0,6,0,0,25*0/
      DATA NATOMS,ICARTR,MDER,MSURF/ia,ib,ic,id/
      DATA NASURF/1,35*0/
c     Data statments from the original code
         
      END
c     This space is for other subroutines from the original code

