c***********************************************************************
c POTLIB 2001: A potential energy surface library for chemical systems
c R. J. Duchovic, Y. L. Volobuev, G. C. Lynch, D. G. Truhlar, T. C. Allison,
c   A. F. Wagner, B. C. Garrett, and J. C. Corchado
c Computer Physics Communications 144, 169-187 (2002), 156, 319-322(E) (2004)
c
c System: ClH2
c Name: ClH2g3ta3
c Int Coords: R12 R23 R13
c Special Features: Has three different energy calculations, chosen by NFLAG(9)
c
c Original code written by:
c       Unknown
c Transcribed to the POTLIB Standard by:
c       Aaron Tagliaboschi <aaron.tagliaboschi@gmail.com>
c       Western Kentucky University, Department of Chemestry
c         Dr. Jeremy B. Maddox
c***********************************************************************
      SUBROUTINE PREPOT
      IMPLICIT REAL*8 (A-H,O-Z)
      character text1,text2,text3,text4,text5
      CHARACTER*75 REF(5)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF            
c     Extra common blocks
      COMMON /POTCM/  de(3),re(3),b(3),s(3),pf(3),spl(3),rj,ap,at,bp,
     X                q2,q4,ckcau,cangau,c1,c2,esep,zl,wvnsep,sqrt3
     
      REF(1)= "Potential code called for ClH2 surface TA3."
      REF(2)= ""
      REF(3)= ""
      REF(4)= ""
      REF(5)= ""
      IRCTNT = 2
      INDEXES(1) = 17
      INDEXES(2) = 1
      INDEXES(3) = 1
      CALL POTINFO
      CALL ANCVRT
c     Initializations and onetime calculations

      data b      / 1.8674d0, 1.9413d0, 1.8674d0 /
      data de     / 106.447d0, 109.458d0, 106.447d0 /
      data re     / 1.2732d0, 0.74127d0, 1.2732d0 /
      data ckcau  / 627.5095d0 /
      data cangau / 0.529177249d0 /
      data s      / 0.1835d0, 0.167d0, 0.1835d0 /
      data rj     / 0.0758016022d0 /
      data ap     / 0.0008627355d0 /
      data at     / 0.2981969860d0 /
      data bp     / 0.1439106543d0 /
      data q2     / 0.6940323070d0 /
      data q4     / 1.6963092005d0 /
      data c1     / 3.1053877618397071175758280546d+0 /
      data c2     / -1.8942608491155350536449828878d+0 /
      data spl    / 2.6477d0, 1.8701d0, 4.5178d0 /
      data wvnsep / 100.d0 /

      do 10 i=1,3
        de(i)=de(i)/ckcau
        re(i)=re(i)/cangau
        b(i)=b(i)*cangau
        pf(i)=0.5d0*de(i)*((1.d0-s(i))/(1.d0+s(i)))
10    continue

      esep=136.748116495694973d0/ckcau
      
      zl=(esep-100.d0*2.859144d-3/ckcau)/esep
c      zl=0.0d0
      
      sqrt3=dsqrt(3.d0)

      RETURN
      END

      SUBROUTINE POT
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
      dimension e(3)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF            
      COMMON /USROCM/ PENGYGS,PENGYES(ISURF),PENGYIJ(JSURF),  
     X                DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     X                DIJCART(NATOM,3,JSURF)                  
      COMMON /PT1CM/  R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)          
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)       
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)       
c     Extra common blocks  
      COMMON /POTCM/  de(3),re(3),b(3),s(3),pf(3),spl(3),rj,ap,at,bp,
     X                q2,q4,ckcau,cangau,c1,c2,esep,zl,wvnsep,sqrt3
       
      CALL CARTOU
      CALL CARTTOR
c     Start of original code
      
      f1=dexp(-ap*(r(1)+r(3))**4)
      f2=dexp(-at*(r(1)-r(3))**2)
      f3=dexp(-bp*(r(1)+r(3)-r(2))**2)
      ctheta=(r(1)**2+r(3)**2-r(2)**2)/(2.d0*r(1)*r(3))
      if (ctheta .lt. -1.d0) then
        ctheta=-1.d0
      else if (ctheta .gt. 1.d0) then
        ctheta=1.d0
      end if
      stheta=dsin(dacos(ctheta))
      f4=1.d0+q2*stheta**2+q4*stheta**4
      v3b=f4*(rj*f1*f2*f3)

      x=(r(1)-r(3))/r(2)
      g=(0.5d0*(1.d0+ctheta))**6
      fmod=1.d0+c1*g*(1.d0-x**2)+c2*g*(1.d0-x**4)

      y1=dexp(-b(1)*(r(1)-re(1)))
      vs1=de(1)*(y1-2.d0)*y1
      vt1=pf(1)*(fmod*y1+2.d0)*y1
      y2=dexp(-b(2)*(r(2)-re(2)))
      vs2=de(2)*(y2-2.d0)*y2
      vt2=pf(2)*(y2+2.d0)*y2
      y3=dexp(-b(3)*(r(3)-re(3)))
      vs3=de(3)*(y3-2.d0)*y3
      vt3=pf(3)*(fmod*y3+2.d0)*y3
c
c Calculate original Hamiltonian matrix and obtain eigenvalues and
c eigenvectors
c
      h11=de(2)+vs2+0.25d0*(vs1+vs3)+0.75d0*(vt1+vt3)+v3b
      h22=de(2)+0.75d0*(vs1+vs3)+vt2+0.25d0*(vt1+vt3)+v3b
      h12=sqrt3*0.25d0*(vs1-vs3+vt3-vt1)
      h12p=-dsqrt(h12**4/((1.227d-2*dexp(-2.d0*b(2)*r(2))**2)+h12**2))
      call dlaev2(h11,h12p,h22,rt1,rt2,cs1,sn1)
c
c Calculate new value of e2 and backtransform to new Hamiltonian matrix
c
      z=(rt1-rt2)/esep
      r10=4.0d0
      r20=0.5d0
      cfd=2.d0*r(1)**2+2.d0*r(3)**2-r(2)**2
      if (cfd .lt. 0.d0) cfd=0.d0
      if (cfd .eq. 0.d0) then
        coschi2=1.d0
      else
        coschi2=((r(3)**2-r(1)**2)/(r(2)*dsqrt(cfd)))**2
      end if
      if (coschi2 .gt. 1.d0) coschi2=1.d0
      if (coschi2 .lt. -1.d0) coschi2=-1.d0
      sinchi4=(1.d0-coschi2)**2
      e2new=rt1-(1.d0-sinchi4)*0.5d0*zl*(rt1-rt2)*((z+1.d0)
     &      *dexp(0.5d0*(1.d0-z)))*dexp(-((r(1)-spl(1)+r(3)
     &      -spl(3))**4/r10**4+(r(2)-spl(2))**4/r20**4))

      tp11=rt2*sn1*sn1+e2new*cs1*cs1
      tp12=(e2new-rt2)*sn1*cs1
      tp22=rt2*cs1*cs1+e2new*sn1*sn1
c
c Calculate scaling factor for coupling and transform matrix
c
      chi=0.5d0*(1.d0-dtanh(b(2)*(r(2)-5.3d0)))
      if (chi .lt. 0.d0) chi=0.d0
      if (chi .gt. 1.d0) chi=1.d0
      ab=tp11-tp22
      ac=ab*ab+4.d0*(1.d0+chi)*tp12*tp12
      ad=ab*dsqrt(ab*ab+4.d0*(1.d0-chi*chi)*tp12*tp12)
      ae1=ac+ad
      if (ae1 .lt. 0.d0) ae1=0.d0
      ae2=ac-ad
      if (ae2 .lt. 0.d0) ae2=0.d0
      af=8.d0*tp12*tp12+2.d0*ab*ab
      ax1=dsqrt(ae1/af)
      if (ax1 .gt. 1.d0) ax1=1.d0
      ax2=dsqrt(ae2/af)
      if (ax2 .gt. 1.d0) ax2=1.d0
      ax=max(ax1,ax2)
      if (tp12 .ge. 0.d0) then
        ay=-dsqrt(1.d0-ax*ax)
      else
        ay=dsqrt(1.d0-ax*ax)
      end if
c
c Return value of appropriate surface
c
      e(1)=tp11*ax*ax-2.d0*tp12*ax*ay+tp22*ay*ay
      e(2)=(tp11-tp22)*ax*ay+tp12*(ax*ax-ay*ay)
      e(3)=tp11*ay*ay+2.d0*tp12*ax*ay+tp22*ax*ax
      
      ENGYGS=e(1)
      
c     End of original code        
      CALL EUNITZERO
      IF(NDER.NE.0) THEN
         CALL RTOCART
         IF(NFLAG(1)+NFLAG(2).NE.0) CALL DEDCOU
      ENDIF
      RETURN
      END
      
      BLOCK DATA PTPACM
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (NATOM=25,N3ATOM=3*NATOM)
      PARAMETER (ISURF=5,JSURF=ISURF*(ISURF+1)/2)
      PARAMETER (ia=1,ib=1,ic=1,id=1)
      COMMON /USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),     
     X                NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER  
      COMMON /INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),IRCTNT,  
     X                NATOMS,ICARTR,MDER,MSURF,REF 
      COMMON /POTCM/  de(3),re(3),b(3),s(3),pf(3),spl(3),rj,ap,at,bp,
     X                q2,q4,ckcau,cangau,c1,c2,esep,zl,wvnsep,sqrt3     
c     Extra common blocks
        
      DATA NDER,ANUZERO,NFLAG,NULBL/0,0.0,1,1,15*0,6,0,0,25*0/
      DATA NATOMS,ICARTR,MDER,MSURF/3,3,0,0/
      DATA NASURF/1,35*0/
c     Data statments from the original code

         
      END
c     This space is for other subroutines from the original code

      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
c
c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
c
c     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
c     ..
c
c  Purpose
c  =======
c
c  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
c     [  A   B  ]
c     [  B   C  ].
c  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
c  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
c  eigenvector for RT1, giving the decomposition
c
c     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
c     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
c
c  Arguments
c  =========
c
c  A       (input) DOUBLE PRECISION
c          The (1,1) element of the 2-by-2 matrix.
c
c  B       (input) DOUBLE PRECISION
c          The (1,2) element and the conjugate of the (2,1) element of
c          the 2-by-2 matrix.
c
c  C       (input) DOUBLE PRECISION
c          The (2,2) element of the 2-by-2 matrix.
c
c  RT1     (output) DOUBLE PRECISION
c          The eigenvalue of larger absolute value.
c
c  RT2     (output) DOUBLE PRECISION
c          The eigenvalue of smaller absolute value.
c
c  CS1     (output) DOUBLE PRECISION
c  SN1     (output) DOUBLE PRECISION
c          The vector (CS1, SN1) is a unit right eigenvector for RT1.
c
c  Further Details
c  ===============
c
c  RT1 is accurate to a few ulps barring over/underflow.
c
c  RT2 may be inaccurate if there is massive cancellation in the
c  determinant A*C-B*B; higher precision or correctly rounded or
c  correctly truncated arithmetic would be needed to compute RT2
c  accurately in all cases.
c
c  CS1 and SN1 are accurate to a few ulps barring over/underflow.
c
c  Overflow is possible only if RT1 is within a factor of 5 of overflow.
c  Underflow is harmless if the input data is 0 or exceeds
c     underflow_threshold / macheps.
c
c =====================================================================
c
c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
c     ..
c     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
c
c     .. Executable Statements ..
c
c     Compute the eigenvalues
c
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
c
c       Includes case AB=ADF=0
c
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
c
c        Order of execution important.
c        To get fully accurate smaller eigenvalue,
c        next line needs to be executed in higher precision.
c
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1

c        Order of execution important.
c        To get fully accurate smaller eigenvalue,
c        next line needs to be executed in higher precision.
c
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
c
c        Includes case RT1 = RT2 = 0
c
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
c
c     Compute the eigenvector
c
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
c
c     End of DLAEV2
c
      END
