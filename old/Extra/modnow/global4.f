      Program global
      Implicit Real*8 (A-H,O-Z)
C
      PARAMETER(IRPL=9)
      PARAMETER(ITL=80)
C
      character*50  title1
      character*50  title2
      character*50  title3
      character*50  title4
      character*50  title5
      character*50  title6
      character*50  title7
      character*50  title8
      character*50  title9
      character*50  title10
      character*50  title11
      character*62  title12
      character*50  title13
      character*50  title14
      character*62  title15
      character*62  title16
      character*50  title17
      character*50  title18
      character*62  title19
      character*50  title20
      character*50  title21
      character*62  title22
      character*62  title23
      character*50  title24
C
      character*1 title(80)
C
      DATA title1 /'C                                                 '/
      DATA title2 /'      CHARACTER*75 REF(5)                         '/
      DATA title3 /'C                                                 '/
      DATA title4 /'      PARAMETER (I3ATOMS = 3)                     '/
      DATA title5 /'      PARAMETER (ISURF = 1)                       '/
      DATA title6 /'      PARAMETER (JSURF = ISURF*(ISURF+1)/2)       '/
      DATA title7 /'      PARAMETER (PI = 3.141592653589793D0)        '/
      DATA title8 /'      PARAMETER (ANUZERO=0.0D0)                   '/
      DATA title9 /'      PARAMETER (NATOMS = 3)                      '/
      DATA title10/'      PARAMETER (IREORDER = 0)                    '/
      DATA title11/'C                                                 '/
      DATA 
     +title12
     +/'      COMMON /PV1CM/ R(I3ATOMS),ENGYGS,DEGSDR(I3ATOMS)        '/
      DATA title13/'      COMMON /PV2CM/ NSURF, NDER, NFLAG(20)       '/
      DATA title14/'      COMMON /PV3CM/ EZERO                        '/
      DATA 
     +title15
     +/'      COMMON /PV4CM/ ENGYES(ISURF),DEESDR(I3ATOMS,ISURF)      '/
      DATA
     +title16
     +/'      COMMON /PV5CM/ ENGYIJ(JSURF),DEIJDR(I3ATOMS,JSURF)      '/
      DATA title17/'      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC         '/
      DATA title18/'C                                                 '/ 
      DATA
     +title19
     +/'      COMMON /INFOCM/ CART(NATOMS,3),DGSCART(NATOMS,3),       '/
      DATA title20/'     +                DESCART(NATOMS,3,ISURF),    '/
      DATA title21/'     +                DIJCART(NATOMS,3,JSURF),    '/
      DATA
     +title22
     +/'     +                NULBL(NATOMS),INDEXES(NATOMS),          '/
      DATA
     +title23
     +/'     +                IAMAN(NATOMS,2),IASYM(NATOMS,2),REF     '/
      DATA title24/'C                                                 '/
C
 1    continue
      DO KK=1,80
         TITLE(KK) = ' '
      END DO
      read(5,2,end=99) title
 2    format(80a1)
      do i = 1,69
C*****************************************
C                                        *
C   The first "IF" check results in      *
C   Comment Lines remaining unchanged.   *
C   Otherwise, the large "IF" check      *
C   can cause a Comment line to be       *
C   deleted.                             *
C                                        *
C*****************************************
         if(  (title(1).eq.' ').and.
     +       ((title(i).eq.'C').or.(title(i).eq.'c')).and.
     +       ((title(i+1).eq.'O').or.(title(i+1).eq.'o')).and.
     +       ((title(i+2).eq.'M').or.(title(i+2).eq.'m')).and.
     +       ((title(i+3).eq.'M').or.(title(i+3).eq.'m')).and.
     +       ((title(i+4).eq.'O').or.(title(i+4).eq.'o')).and.
     +       ((title(i+5).eq.'N').or.(title(i+5).eq.'n')).and.
     +       (title(i+6).eq.' ').and.
     +       (title(i+7).eq.'/').and.
     +       (title(i+8).eq.'P').and.
     +       (title(i+9).eq.'T').and.
     +       ((title(i+10).eq.'3').or.(title(i+10).eq.'4')).and.
     +       ((title(i+11).eq.'1').or.(title(i+11).eq.'2').or.
     +        (title(i+11).eq.'4').or.(title(i+11).eq.'5')) 
     +          ) go to 10
      end do
      write(6,2) title
      do i = 1,73
         if( ((title(i).eq.'I').or.(title(i).eq.'i')).and.
     +       ((title(i+1).eq.'M').or.(title(i+1).eq.'m')).and.
     +       ((title(i+2).eq.'P').or.(title(i+2).eq.'p')).and.
     +       ((title(i+3).eq.'L').or.(title(i+3).eq.'l')).and.
     +       ((title(i+4).eq.'I').or.(title(i+4).eq.'i')).and.
     +       ((title(i+5).eq.'C').or.(title(i+5).eq.'c')).and.
     +       ((title(i+6).eq.'I').or.(title(i+6).eq.'i')).and.
     +       ((title(i+7).eq.'T').or.(title(i+7).eq.'t')) ) then
               WRITE(6,3) TITLE1
               WRITE(6,3) TITLE2
               WRITE(6,3) TITLE3
               WRITE(6,3) TITLE4
               WRITE(6,3) TITLE5
               WRITE(6,3) TITLE6
               WRITE(6,3) TITLE7
               WRITE(6,3) TITLE8
               WRITE(6,3) TITLE9
               WRITE(6,3) TITLE10
               WRITE(6,3) TITLE11
               WRITE(6,4) TITLE12
               WRITE(6,3) TITLE13
               WRITE(6,3) TITLE14
               WRITE(6,4) TITLE15
               WRITE(6,4) TITLE16
               WRITE(6,3) TITLE17
               WRITE(6,3) TITLE18
               WRITE(6,4) TITLE19
               WRITE(6,3) TITLE20
               WRITE(6,3) TITLE21
               WRITE(6,4) TITLE22
               WRITE(6,4) TITLE23
               WRITE(6,3) TITLE24
 3             format(a50)
 4             format(a62)
         endif
      enddo
 10   continue
      go to 1
 99   continue
      stop
      end
