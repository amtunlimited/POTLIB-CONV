      Program global
      Implicit Real*8 (A-H,O-Z)
C
      PARAMETER(ITN1=19)
      PARAMETER(ITN2=25)
      PARAMETER(ITN3=25)
      PARAMETER(ITN4=25)
      PARAMETER(ITN5=25)
      PARAMETER(IRPL=23)
      PARAMETER(ITL=80)
C
      character*1 title1(ITN1),title2(ITN2),title3(ITN3),
     +            title4(ITN4),title5(ITN5)
      character*1 title(ITL)
      character*1 replace1(IRPL)
C
      DATA replace1 /'N','A','S','U','R','F','(','I',
     +               'S','U','R','F','+','1',',','I',
     +               'S','U','R','F','+','1',')'/
      DATA title1   /'N','A','S','U','R','F','(','I',
     +               'S','U','R','F',',','I',
     +               'S','U','R','F',')'/

      IFORM=0
 1    continue
      DO KK=1,ITL
         TITLE(KK)=' '
      END DO
      read(5,2,end=99) title
 2    format(80a1)
      IF( ((title(1).eq.'C').or.(title(1).eq.'c')).and.
     +     (IFORM.eq.1)           ) IFORM=0
      IF( ((title(1).ne.'C').or.(title(1).ne.'c')).and.
     +     (title(6).eq.' ').and.
     +     (IFORM.eq.1)           ) IFORM=0
      IF( (title(1).eq.'C').or.(title(1).eq.'c')) then
          write(6,2) title
          go to 10
      ENDIF
C************************************
C                                   *
C   This piece of code will leave   *
C   FORMAT statements unchanged.    *
C                                   *
C************************************
      if( (IFORM.eq.1).and.(title(6).ne.' ') ) then
          write(6,2) title
          go to 10
      end if
      DO JJ=1,65
         if( ((title(JJ).eq.'F').or.(title(JJ).eq.'f')).and.
     +       ((title(JJ+1).eq.'O').or.(title(JJ+1).eq.'o')).and.
     +       ((title(JJ+2).eq.'R').or.(title(JJ+2).eq.'r')).and.
     +       ((title(JJ+3).eq.'M').or.(title(JJ+3).eq.'m')).and.
     +       ((title(JJ+4).eq.'A').or.(title(JJ+4).eq.'a')).and.
     +       ((title(JJ+5).eq.'T').or.(title(JJ+5).eq.'t')).and.
     +       ((title(JJ+6).eq.'(').or.((title(JJ+6).eq.' ').and.
     +        (title(JJ+7).eq.'('))) )THEN
             write(6,2) title
             IFORM=1
             go to 10
         end if
      END DO
      iflag = 0
      IFLAG2 = 0
      IF( (title(1).eq.'C').or.(title(1).eq.'c')) then
          IFIRST=80
      ELSE
          IFIRST=72
      ENDIF
      do i = 1,IFIRST-ITN1+1
         if( (title(i).eq.title1(ITN1-18)).and.
     +       (title(i+1).eq.title1(ITN1-17)).and.
     +       (title(i+2).eq.title1(ITN1-16)).and.
     +       (title(i+3).eq.title1(ITN1-15)).and.
     +       (title(i+4).eq.title1(ITN1-14)).and.
     +       (title(i+5).eq.title1(ITN1-13)).and.
     +       (title(i+6).eq.title1(ITN1-12)).and.
     +       (title(i+7).eq.title1(ITN1-11)).and.
     +       (title(i+8).eq.title1(ITN1-10)).and.
     +       (title(i+9).eq.title1(ITN1-9)).and.
     +       (title(i+10).eq.title1(ITN1-8)).and.
     +       (title(i+11).eq.title1(ITN1-7)).and.
     +       (title(i+12).eq.title1(ITN1-6)).and.
     +       (title(i+13).eq.title1(ITN1-5)).and.
     +       (title(i+14).eq.title1(ITN1-4)).and.
     +       (title(i+15).eq.title1(ITN1-3)).and.
     +       (title(i+16).eq.title1(ITN1-2)).and.
     +       (title(i+17).eq.title1(ITN1-1)).and.
     +       (title(i+18).eq.title1(ITN1)) ) THEN
              IEND=IFIRST-ITN1+IRPL
                 IF(IEND.LE.IFIRST) THEN
                    write(6,2) (title(k),k=1,i-1),
     +                         replace1,
     +                         (title(k),k=i+ITN1,IEND)
                 ELSE
                    write(6,2) (title(k),k=1,i-1),
     +                         replace1,
     +                         (title(k),k=i+ITN1,IFIRST-IRPL+ITN1)
                    ISTART=IFIRST-ITN1+1
                    DO JJ=ISTART,IFIRST
                       IF(TITLE(JJ).NE.' ') THEN
                         IF((TITLE(1).EQ.'C').OR.
     +                      (TITLE(1).EQ.'c')) THEN
                             WRITE(6,3) 
     +                        (title(k),k=IFIRST-IRPL+ITN1+1,IFIRST)
 3                           FORMAT('C     ',80a1)
                         ELSE
                            WRITE(6,3) 
     +                       (title(k),k=IFIRST-IRPL+ITN1+1,IFIRST)
 4                          FORMAT('     +',80a1)
                         ENDIF
                         IFLAG2=1
                       ENDIF
                       if(iflag2.eq.1) go to 10
                    END DO
                 ENDIF
            iflag = 1
         endif
         if(iflag.eq.1) go to 10
      enddo
      if(iflag.eq.0) write(6,2) title
 10   continue
      go to 1
 99   continue
      stop
      end


