      Program scrptwrtr
      Implicit Real*8 (A-H,O-Z)
      character*1 name(23),blank
      data blank/' '/
      Read(5,15) ntitle
15    format(i3)
      write(6,*) '# Script to examine first changes potlib routines'
      write(6,*) '#'
      write(6,1)
1     format(' echo " "')
         write(6,*) 'date'
         write(6,1)
      write(6,1)
      do i = 1,ntitle
         read(5,5) name
 5       format(40a1)
         do j = 1,23
            if(name(j).eq.blank) then
               jend = j-1
               go to 10
            endif
         enddo
 10      continue
         write(6,101) name
101      format(' echo "  **** doing ',23a1,' ****  "')
         write(6,1)
         write(6,*) 
     +      ' cd /thry/potlib/potlib-2006/new'
         write(6,102) name
102      format(' xlf ',23a1,
     +  ' /thry/potlib/potlib-2006/examnew/examine.o',
     +  ' -o /thry/potlib/potlib-2006/examnew/examine.x')
         write(6,*) 
     +      ' cd /thry/potlib/potlib-2006/examnew'
         write(6,*) ' examine.x > examine.out'
         write(6,1)
         write(6,*) ' more examine.out'
         write(6,1)
         write(6,103) name
103      format(' echo "  **** thru ',23a1,' ****  "')
         write(6,1)
      enddo
      stop
      end
