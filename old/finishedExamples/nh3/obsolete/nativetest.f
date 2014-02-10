c to read in bohr the cartesian coordinates of nh3, call nh3, and obtain
c energy and derivatives in return
c
      implicit real*8 (a-h,o-z)
      dimension xb(12),dxb(12)
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
c
c call nh3 potential
c
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
