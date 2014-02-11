       PROGRAM test
       
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       COMMON /PT31CM/ CARTX(9), ENERGY, DEDR(3)
       
       do i=1,9
         read (5,*) CARTX(i)
       enddo
       
       call NPREPOT
       call NPOT
       
       write (6,*) 'ENERGY: ', ENERGY
       
       end
