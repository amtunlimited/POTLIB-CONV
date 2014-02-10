c to read in bohr the cartesian coordinates of nh3, call nh3, and obtain
c energy and derivatives in return
c
      implicit real*8 (a-h,o-z)
      dimension xb(12),dxb(12)
      PARAMETER (NATOM=25,ISURF=5,JSURF=ISURF*(ISURF+1)/2)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER
      COMMON /USROCM/ PENGYGS,PENGYES(ISURF),PENGYIJ(JSURF),
     X                DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     X                DIJCART(NATOM,3,JSURF)
c
c read in input
C X(I) = CARTESIAN CORDINATES
C     I = 1 - 3 :   X, Y, Z, OF N
C     I = 4 - 6 :   X, Y, Z, OF H1
C     I = 7 - 9 :   X, Y, Z, OF H2
C     I = 10 - 12 : X, Y, Z, OF H3
c
      n3tm = 12
      do i = 1,n3tm
         read(5,*) xb(i)
      enddo
      write(6,*) 'x,y,z for n,h1,h2,h3 in bohr'
      do i = 1,10,3
         write(6,101) (xb(j),j=i,i+2)
101      format(8f10.5)
      enddo
      do i = 1,10,3
         j = 1+i/3
         l = 0
         do k = i,i+2
            l = l+1
            cart(j,l) = xb(k)
         enddo
         write(6,103) j
103      format('cart(',i1,',l),l=1,2,3')
         write(6,101) (cart(j,l),l=1,3)
      enddo
c
c call nh3 potential
c
      call prepot
      call pot
      write(6,*) 'POT energy in au,kcal'
      write(6,101) pengygs,pengygs*627.5096d0
      write(6,*) 'POT derivatives in au'
      do i = 1,4
         write(6,101) (dgscart(i,j),j=1,3)
      enddo
      call setup(n3tm)
      call surf(v,xb,dxb,n3tm)
      write(6,*) 'energy in au,kcal'
      write(6,101) v,v*627.5096d0
      write(6,*) 'derivatives in au'
      do i = 1,10,3
         write(6,101) (dxb(j),j=i,i+2)
      enddo
      stop
      end
