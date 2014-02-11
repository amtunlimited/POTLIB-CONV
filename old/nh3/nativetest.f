       PROGRAM TEST
       
       dimension XB(12), DXB(12)
       
c      read in the file
       do i=1,12
         read (5,*) XB(i)
       enddo
       
c      native calls
       call setup(12)
       call surf(V, XB, DXB, 12)
       
c      write out the results
       write (6,*) "Potential: ",V
       write (6,*) "Derivatives:"
       do i=1,12
         write (6,*) "  ",DXB(i)
       enddo
