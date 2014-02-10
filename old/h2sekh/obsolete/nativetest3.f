       PROGRAM test
       
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (NATOM=25,N3ATOM=3*NATOM)
       PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
       COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER
       COMMON /USROCM/ PENGYGS,PENGYES(ISURF),PENGYIJ(JSURF),  
     X                DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     X                DIJCART(NATOM,3,JSURF)
       COMMON /PT31CM/ CARTX(9), ENERGY, DEDR(3)
       
       do i=1,9
         read (5,*) CARTX(i)
       enddo
       
       call PREPOT
       
       call NPREPOT
       call NPOT
       
       write (6,*) 'ENERGY: ', ENERGY
       
       end
