      Program scrptwrtr
      Implicit Real*8 (A-H,O-Z)
      character*20 name(25)
      Read(5,15) ntitle
15    format(i2)
      write(6,*) '# Script to clean and copy directories'
      write(6,*) '#'
      write(6,1)
1     format(' echo " "')
      write(6,*) 'date'
      write(6,1)
      do i = 1,ntitle
         read(5,5) name(i)
 5       format(a20)
         write(6,101) name(i)
 101     format(' echo "  **** cleaning ',a20,' in now ****  "')
         write(6,102) name(i)
 102     format(' cd /thry/potlib/potlib-2006/now/',a20)
         write(6,*) ' rm *.*'
         write(6,104) name(i)
 104     format(' echo "  **** thru ',a20,' ****  "')
         write(6,1)
      enddo
         write(6,*) ' echo "  **** copying new to now ****  "'
         write(6,1)
         write(6,*) ' cp -pR /thry/potlib/potlib-2006/new/*',
     +              ' /thry/potlib/potlib-2006/now'
      write(6,1)
      do i = 1,ntitle
         write(6,121) name(i)
 121     format(' echo "  **** cleaning ',a20,' in new ****  "')
         write(6,122) name(i)
 122     format(' cd /thry/potlib/potlib-2006/new/',a20)
         write(6,*) ' rm *.*'
         write(6,124) name(i)
 124     format(' echo "  **** thru ',a20,' ****  "')
         write(6,1)
      enddo
      stop
      end
