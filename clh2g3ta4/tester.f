      Implicit Real*8 (A-H,O-Z)
C
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
C READ INPUT POSITIONS,DERIVATIVES, EXCITED STATE, NFLAG1-5 INFO
C
      read(5,*) ncart
      n = ncart/3
      do i = 1,n
         do j = 1,3
            read(5,*) cart(i,j)
         enddo
      enddo
      read(5,*) natoms, nder, nasurf, (nflag(i),i=1,5)
C
C CALL POTENTIAL
C
      Call prepot
      Call pot
C
C LIST RESULTS
C
C GROUND ELECTRONIC STATE OUTPUT
      write(6,100) pengygs           
100   format(2x,'GS Energy',f25.14,/)
      if(nder.eq.1) then
         write(6,200) ((dgscart(i,j),j=1,3),i=1,natoms)
200      format(2x, 'Derivative',5x,3f18.11)
      endif
C
C EXCITED ELECTRONIC STATE OUTPUT
C
      do i = 2,isurf+1
         if(nasurf(i,i).ne.0) then
            write(6,101) i-1,pengyes(i-1)           
101         format(2x,i1,'th Excited State Energy',g25.15)
            if(nder.eq.1) 
     x         write(6,200) ((descart(m,j,i-1),j=1,3),m=1,natoms)
         endif
      enddo
C
C COUPLING ELECTRONIC STATE OUTPUT
C
      if(nflag(4).eq.0) then
         l = 0
         do i = 1,isurf
            do j = i+1,isurf+1
               l = l + 1
               if(nasurf(i,j)+nasurf(j,i).ne.0) then
                  write(6,102) i,j,pengyij(l)           
102               format(2x,i1,'-',i1,' Coupling Energy     ',g25.15)
                  if(nder.eq.1) 
     x               write(6,200) ((dijcart(m,n,l),n=1,3),m=1,natoms)
               endif
            enddo
         enddo
      else
         write(6,112) (pengyij(i),i=1,3)
112      format(2x,' nflag(4) controled coupling energies'/3g25.15)
      endif
      stop
      end




