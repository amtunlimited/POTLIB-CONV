c first stage modifier of a native potential subroutine to make it
c potlib compliant.  Based on Figs 1 - 3 in CPC 144, 169 (2002)
c
      Implicit Real*8 (A-H,O-Z)
      character title(80)
      character*62 title1(2),title2(2),title3(3),title4(3),title5
      character*62 title6(4),title7(3),title8(3),title9(3),title10(4)
      character*62 title11(3),title12(3),title13
      data title1/
     x '      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     ',
     x '     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  '/
      data title2/
     x '      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  ',
     x '     X                NATOMS,ICARTR,MDER,MSURF,REF            '/
      data title3/
     x '      COMMON /USROCM/ PENGYGS,PENGYES(ISURF),PENGYIJ(JSURF),  ',
     x '     X                DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),',
     x '     X                DIJCART(NATOM,3,JSURF)                  '/
      data title4/
     x '      COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)          ',
     x '      COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)       ',
     x '      COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)       '/
      data title5/
     x '      DATA NDER,ANUZERO,NFLAG,NULBL/0,0.0,1,1,15*0,6,0,0,25*0/'/
      data title6/
     x 'c****************************************************         ',
     x 'c ID chemical system, common name of potential,               ',
     x 'c functional form of PES, int coor, special features          ',
     x 'c****************************************************         '/
      data title7/
     x 'c****************************************************         ',
     x 'c any data or labeled commons for PTPACM or POT               ',
     x 'c****************************************************         '/
      data title8/
     x 'c****************************************************         ',
     x 'c any assignments for labeled common share with pot           ',
     x 'c****************************************************         '/
      data title9/
     x 'c****************************************************         ',
     x 'c any labeled common statements for PTPACM or POT             ',
     x 'c****************************************************         '/
      data title10/
     x 'c****************************************************         ',
     x 'c code for V,dVdr as functions of int coord R                 ',
     x 'c or for polyatomics as functions of cartesian coord          ',
     x 'c****************************************************         '/
      data title11/
     x 'c****************************************************         ',
     x 'c all other labeled commons for PREPOT or POT                 ',
     x 'c****************************************************         '/
      data title12/
     x 'c****************************************************         ',
     x 'c all other data statements                                   ',
     x 'c****************************************************         '/
      data title13/
     x 'c     ...                                                     '/
c
c place prepot at the head of the new routine (Fig.2 of CPC paper)
c
      write(6,111) title6
111   format(a62)
      write(6,*) '     SUBROUTINE PREPOT'
      write(6,*) '     IMPLICIT REAL*8 (A-H,O-Z)'
      write(6,*) '     character text1,text2,text3,text4,text5'
      write(6,*) '     CHARACTER*75 REF(5)'
      write(6,*) '     PARAMETER (NATOM=25,N3ATOM=3*NATOM)'
      write(6,*) '     PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)'
      write(6,111) title1
      write(6,111) title2
      write(6,111) title7
      write(6,*) '     REF(1)= text1'
      write(6,*) '     REF(2)= text2'
      write(6,*) '     REF(3)= text3'
      write(6,*) '     REF(4)= text4'
      write(6,*) '     REF(5)= text5'
      write(6,*) '     IRCTNT = aa'
      write(6,*) '     INDEXES(1) = bb1'
      write(6,*) '     INDEXES(2) = bb2'
      write(6,111) title13
      write(6,*) '     INDEXES(NATOMS) = bbnatoms'
      write(6,*) '     CALL POTINFO'
      write(6,*) '     CALL ANCVRT'
      write(6,111) title8
      write(6,*) '     RETURN'
      write(6,*) '     END'
c
c create a pot entry right after prepot (Fig. 3 of CPC paper)
c
      write(6,*) '     SUBROUTINE POT'
      write(6,*) '     IMPLICIT REAL*8 (A-H,O-Z)'
      write(6,*) '     CHARACTER*75 REF(5)'
      write(6,*) '     PARAMETER (NATOM=25,N3ATOM=3*NATOM)'
      write(6,*) '     PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)'
      write(6,111) title1
      write(6,111) title2
      write(6,111) title3
      write(6,111) title4
      write(6,111) title9
      write(6,*) '     CALL CARTOU'
      write(6,*) '     CALL CARTTOR'
      write(6,111) title10
      write(6,*) '     CALL EUNITZERO'
      write(6,*) '     IF(NDER.NE.0) THEN'
      write(6,*) '        CALL RTOCART'
      write(6,*) '        IF(NFLAG(1)+NFLAG(2).NE.0) CALL DEDCOU'
      write(6,*) '     ENDIF'
      write(6,*) '     RETURN'
      write(6,*) '     END'
c
c read native routine
c
 1    continue
      read(5,2,end=99) title
      write(6,2) title
 2    format(80a1)
      go to 1
 99   continue
c
c write out blockdata (Fig. 1 of CPC paper)
c
      write(6,*) '     BLOCK DATA PTPACM'
      write(6,*) '     IMPLICIT REAL*8 (A-H,O-Z)'
      write(6,*) '     CHARACTER*75 REF(5)'
      write(6,*) '     PARAMETER (NATOM=25,N3ATOM=3*NATOM)'
      write(6,*) '     PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)'
      write(6,*) '     PARAMETER (ia=1,ib=1,ic=1,id=1)'
      write(6,111) title1
      write(6,111) title2
      write(6,111) title11
      write(6,111) title5
      write(6,*) '     DATA NATOMS,ICARTR,MDER,MSURF/ia,ib,ic,id/'
      write(6,*) '     DATA NASURF/1,35*0/'
      write(6,111) title12
      write(6,*) '     END'
      stop
      end
