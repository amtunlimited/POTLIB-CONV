      subroutine prepot
      implicit double precision (a-h,o-z)
      dimension array(256),x(3)
      common /pt1cm/ r(3),energy,dedr(3)
      common /pt2cm/ nder,nflag(20)
      external f
        call prepes(1,array)
      return
      entry pot
        tmp=r(2)
        r(2)=r(3)
        r(3)=tmp
        if (nder .eq. 1) then
          n=3
          tol=1.d-15
          maxit=50
          do 20 iwrt=1,3
            h=0.01d0
            x(1)=r(1)
            x(2)=r(2)
            x(3)=r(3)
            if ((iwrt .eq. 1) .or. (iwrt .eq. 3)) then
              xa=x(iwrt)
              xb=x(iwrt)+2.d0*h
            else
              xa=x(iwrt)-2.d0*h
              xb=x(iwrt)
            end if
            call dderiv(f,n,x,iwrt,xa,xb,h,tol,maxit,iflag,fval,dfdx,
     &                  errest,nitrs)
            dedr(iwrt)=dfdx
   20     continue
        end if
        call surfac(v,r,1)
        energy=v
        tmp=r(2)
        r(2)=r(3)
        r(3)=tmp
        tmp=dedr(2)
        dedr(2)=dedr(3)
        dedr(3)=tmp
      return
      end
C
      double precision function f(x)
      implicit double precision (a-h,o-z)
      dimension x(3)
        call surfac(v,x,1)
        f=v
      return
      end
C
      subroutine surfac(v,r,iref)
      implicit double precision (a-h,o-z)
      dimension r(3)
      data pi/3.141592653/
      data az/133.45/,b1,b2,b3,b4,b5/-14.70,0.975,90.289,3.6174,0.45/
      data fz/-28.14/,rz/2.8214/
c     cl+hcl scaled fit to polci calculation with quartic bend due
c     to garrett,truhlar,wagner and dunning,jcp 78,4400(1983)
c     modified to add in bondi et al leps potential for angular displacements
c     away from clhcl collinear and to make morse hcl parameters on leps
c     consistent with polci fit.
      rab=r(1)*0.529177
      rbc=r(3)*0.529177
      call potsp(rab,rbc,vcoll)
      phi=-1.
      if(r(1).ne.0..and.r(3).ne.0.)
     xphi=(r(1)**2+r(3)**2-r(2)**2)/(2.*r(1)*r(3))
      if(phi.gt.1.0)  phi=1.0
      if(phi.lt.-1.0) phi=-1.0
      phi=acos(phi)
      rb=max(r(1),r(3))
      delt=max(rb-rz,0.d0)
      arg=b4*sqrt(delt**3)
      if(arg.gt.60.) arg=60.
      aphi=24.*(b3/cosh(arg)/cosh(b5*delt))
      fphi=2.*(fz+b1*delt**2)*exp(-b2*delt)
      vbend=0.5*fphi*(phi-pi)**2+aphi*(phi-pi)**4/24.
      del=phi-pi
      switch=tanh(20.*(phi-2.60))
      call fpot(r,vleps)
      v=(vcoll+vbend)/627.51*0.5*(1.+switch)+0.5*(1.-switch)*vleps
c     write(6,*) vcoll, vbend,vleps
c     v=v*27.211
      return
      end
c this is an ap version of the clhcl collinear potential
c
      subroutine prepes(ncyc, array)
      implicit double precision (a-h,o-z)
      dimension title(5),array(128)
      dimension csp(48,4),scrap(48),atemp(46),ctemp(48)
      dimension ebnd(46),stuff(205)
      common /proinf/r1s,r2s,vsp,bsp(46,3),phi(46),xbc(3),xab(3),
     1   cbc(3),cab(3),abc(3),aab(3),nphi
      common/data/depr,derg,betpr,betrg,repr,rerg
      equivalence (stuff(1),r1s)
      data detrad/.0174532e0/
c
c indirect disk input into potential routine
c
c      open(unit=17,recl=133,status='old',access='sequential',
c     1   form='formatted')
       open(unit=17,file='pot.dat')
      rewind 17
      read(17,843) title
  843 format(5a4)
      read(17,*) stuff,nphi
      write(6,131) title
  131 format( ' spline coef for nonbending pot stored in file:'
     1   ,5a4)
      write(6,101)
  101 format( ' potential info from the disk:'/)
      write(6,105) r1s,r2s,vsp
  105 format( ' r1s,r2s,vsp = ',3f10.5)
      write(6,111) (xbc(i),i= 1,3)
  111 format( ' bc: de,re,be = ',10x,3f10.5)
      write(6,113) (phi(i),(bsp(i,j),j = 1,3),i = 1,nphi)
  113 format( ' phi,de,re,be = ',4f10.5)
      write(6,115) (xab(i),i = 1,3)
  115 format( ' ab: de,re,be = ',10x,3f10.5)
      write(6,117) (cbc(i),abc(i),i = 1,3),(cab(i),aab(i),i = 1,3)
  117 format( ' for asymptote, mo parameter value = (asymp. value) + c*e
     1xp(-a*(ri-ris))'/' (cbc,abc,i = 1,3) = ',6f10.5/
     1   ' (cab,aab,i = 1,3) = ',6f10.5/)
      aabinf = xbc(1)-xab(1)
c
c cock splines
c
      do 20 k = 1,3
      dis  = spline(1,nphi,phi,bsp(1,k),csp(1,k),scrap,angle )
   20 continue
c
c set up recommended value of abc and aab
c the recommendation is based on making the derivative wrt r1 at r1s
c border or r2 at the r2s border equal at the bottom of the well
c along the border.  on either side of the border, v has the form
c   d(1-exp(-b(r-re)))**2 + a
c where d,b,re, and a are either functions of r1 or r2 or the swing
c angle.  the derivative of the first term is 0.0 at the bottom of the
c well.  therefore the dertermination of alfbc or alfab is that
c   d[ainf - (ainf-a(bdry))*exp(-alf(r-rs))]/dr
c      = da(phi)/dr = d(phi)/dr*da/dphi
c simplified this is
c   alf = (da/dphi)/{(ainf-a(bdry))*d(phi)/dr}
c       = (da/dphi)/{(ainf-a(bdry))*le*sign}
c where the sign depends on whether it is the 0 or 90 bdry.
c
      do 22 i = 1,nphi
      angle = phi(i)
      dis=spline(2,nphi,phi,bsp(1,1),csp(1,1),scrap,angle)
      req=spline(2,nphi,phi,bsp(1,2),csp(1,2),scrap,angle)
      bet=spline(2,nphi,phi,bsp(1,3),csp(1,3),scrap,angle)
      atemp(i) = vsp-dis*(1.-exp(-bet*req))**2
   22 continue
      a90 =spline(1,nphi,phi,atemp,ctemp,scrap,90.0d0)
      a90 =spline(3,nphi,phi,atemp,ctemp,scrap,90.0d0)
      dis =spline(4,nphi,phi,atemp,ctemp,scrap,90.0d0)
      req=spline(2,nphi,phi,bsp(1,2),csp(1,2),scrap,90.0d0)
      alfab = dis/(aabinf-a90)/req
      a0 =spline(3,nphi,phi,atemp,ctemp,scrap,0.0d0)
      dis =spline(4,nphi,phi,atemp,ctemp,scrap,0.0d0)
      req=spline(2,nphi,phi,bsp(1,2),csp(1,2),scrap,0.0d0)
      alfbc = dis/a0/req
      write(6,129) alfbc,alfab
  129 format(' abc and aab which makes potential',
     1   ' derivative cont. at mep on swing angle bdry'/
     1   10x,2f10.5)
      write(6,119)
  119 format(//)
c
c setup any constants
c
      dbc = xbc(1)
      derg = xbc(1)/23.061e0
      rerg = (r2s-xbc(2))/.529177e0
      betrg = xbc(3)*.529177e0
      depr = xab(1)/23.061e0
      repr = (r1s-xab(2))/.529177e0
      betpr = xab(3)*.529177e0
      return
c
c entry into potential calculator
c
      entry potsp(r1,r2,v)
      if(r1.gt.r1s) go to 100
      if(r2.gt.r2s) go to 102
c
c interaction region of polar coordinates
c
      dr1=r1s-r1
      dr2=r2s-r2
      angle=atan(dr1/dr2)
      angle=angle/detrad
      r=sqrt(dr1*dr1+dr2*dr2)
      dis=spline(2,nphi,phi,bsp(1,1),csp(1,1),scrap,angle)
      req=spline(2,nphi,phi,bsp(1,2),csp(1,2),scrap,angle)
      bet=spline(2,nphi,phi,bsp(1,3),csp(1,3),scrap,angle)
   30 v = dis*(1.d0-exp(bet*(r-req)))**2 + vsp-dis*(1.-exp(
     1   -bet*req))**2
      return
c
c large r1 asymptotic recion
c
  100 if(r2.gt.r2s) go to 104
      r = r2s-r2
      sep = r1-r1s
      req = xbc(2) + cbc(2)*dexp(-dmin1(abc(2)*sep,25.d0))
      dis = xbc(1) + cbc(1)*dexp(-dmin1(abc(1)*sep,25.d0))
      bet = xbc(3) + cbc(3)*dexp(-dmin1(abc(3)*sep,25.d0))
      v = dis*(1.d0-dexp(bet*(r-req)))**2 + a0*
     1   dexp(-dmin1(abc(1)*sep,25.d0))
      return
c
c large r2 asymptotic region
c
  102 continue
      r = r1s-r1
      sep = r2-r2s
      dis = xab(1) + cab(1)*exp(-min(aab(1)*sep,25.d0))
      req = xab(2) +cab(2)*exp(-min(aab(2)*sep,25.d0))
      bet = xab(3) + cab(3)*exp(-min(aab(3)*sep,25.d0))
      v = dis*(1.d0-exp(bet*(r-req)))**2 + aabinf - (aabinf-a90)*
     1   exp(-min(aab(1)*sep,25.d0))
      return
c
c large r1,r2 3 body break-up region
c
  104 continue
c      write(6,800)
c  800 format('0',1x,'r1 and r2 both beyond swing point stop')
      v=dis
c      stop
      return
      end
      function spline(isw,nn,x,y,c,d,xpt)
      implicit double precision (a-h,o-z)
      parameter (nc=nn+2, nd=nn+2)
      dimension x(nn),y(nn)
      dimension c(nc),d(nd)
c
c
c  this is a subroutine for fitting data with a cubic spline
c  polynomial and evaluating that polynomial at a given point
c  or its derivative at a given point
c
c
c  calling sequence .......
c
c     isw ... control option
c         isw=1  if a cubic spline is to be fitted to the set of knots
c                defined by the arrays x and y.  the spline coefficients
c                are stored in the array c.
c         isw=2  if the spline defined by the coefficient array 'c' is
c                to be evaluated (interpolated) at the point defined by
c                the parameter 'xpt'.
c         isw=3  as in isw=2, only the derivative/3.d0 is also calculated at xpt
c         isw=4  the derivative calculated by the last use of spline with isw=3
c                is returned.
c
c     nn ... the number of knots (data points) to which the spline is to
c            be fitted
c
c     x,y ... the arrays defining the knots.  the x-values must be in
c             increasing order.  the arrays must be dimensioned at least
c             nn.
c
c     c ... the array that contains the cubic spline coefficients.
c           must be dimensioned at least nn+2 .
c
c     d ... a work space.  must be dimensioned at least nn+2 .
c
c     xpt ... the point at which the interpolation is desired (if isw is
c              set to 2).  the value of spline is set to the
c              interpolated value.
c
c
c *****  user notes  *****
c
c     interpolation involves at least two steps .......
c
c       a.  call spline with the knots.  this sets up the
c           coefficient array c.
c           eg.  dumy=spline(1,nn,x,y,c,d,xpt)
c
c       b.  call spline with the array c which was defined by the
c           previous call and will be used to find the value at the
c           point 'xpt' .
c           eg.   value=spline(2,nn,x,y,c,d,xpt)
c
c
c     step 'a' need be executed only once for a given set of knots.
c     step b may be executed as many times as necessary.
c
c
c
2     n=nn
      np1=n+1
      np2=n+2
      z=xpt
24    go to (4,5,6,7),isw
4     c(1)=y(1)
      d(1)=1.0e0
      c(np1)=0.0e0
      d(np1)=0.0e0
      c(np2)=0.0e0
      d(np2)=0.0e0
      do 41 i=2,n
      c(i)=y(i)-y(1)
41    d(i)=x(i)-x(1)
      do 410 i=3,np2
      if(d(i-1).ne.0)go to 43
      write(6,1001)
      call exit
43    pivot=1.0e0/d(i-1)
      if(i.ge.np2)go to 45
      supd=x(i-1)-x(i-2)
      if(supd.ge.0)go to 44
      write(6,1000)
      call exit
44    supd=supd*supd*supd
      go to 46
45    supd=1.0e0
46    dfact=supd*pivot
      cfact=c(i-1)*pivot
      if(i.gt.n)go to 48
      do 47 j=i,n
      v=x(j)-x(i-2)
      c(j)=c(j)-d(j)*cfact
47    d(j)=v*v*v-d(j)*dfact
48    continue
      if(i.ge.np2)go to 49
      c(np1)=c(np1)-d(np1)*cfact
      d(np1)=1.0e0-d(np1)*dfact
49    c(np2)=c(np2)-d(np2)*cfact
410   d(np2)=x(i-2)-d(np2)*dfact
      do 411 i=1,n
      j=np2-i
      if(j.ne.np1)go to 413
      v=1.0e0
      go to 414
413   v=x(j)-x(j-1)
      v=v*v*v
414   if(d(j+1).ne.0)go to 415
      write(6,1001)
      call exit
415   c(j+1)=c(j+1)/d(j+1)
411   c(j)=c(j)-c(j+1)*v
      if(d(2).ne.0)go to 416
      write(6,1001)
      call exit
416   c(2)=c(2)/d(2)
      return
5     spline=c(1)+c(2)*(z-x(1))
      do 51 i=1,n
      v=z-x(i)
      if(v.gt.0)go to 51
      return
51    spline=spline+c(i+2)*v*v*v
      return
    6 continue
      spline=c(1)+c(2)*(z-x(1))
      deriv = c(2)/3.e0
      do 53 i = 1,n
      v=z-x(i)
      if(v.le.0) return
      v2 = v*v
      spline = spline + c(i+2)*v2*v
      deriv = deriv + c(i+2)*v2
   53 continue
      return
    7 continue
      spline = 3.e0*deriv
      return
1000  format(5x,'***** error in spline ... unordered x-values  *****')
1001  format(5x,' ***** error in spline ... divide fault *****')
      end
      subroutine fpot(ri,v)
      implicit double precision (a-h,o-z)
c     leps potential for cl+hcl from bondi,connor,manz,romelt(mol phys,50,467
c     (1983) modified to make morse hcl parameters consistent with polci
      dimension r(3),re(3),del(3),bet(3),q(3)
      dimension a(3),de(3),ri(3)
c      data de/0.1697216,0.0924230,0.1697216/
c      data re/2.409447,3.756848,2.409447/
c      data bet/0.988484,1.059392,0.988484/
c      parameters appropriate to polci asymptotic potential
      data de/0.169797,0.0924230,0.169797/
      data re/2.40846,3.756848,2.40846/
      data bet/0.988185,1.059392,0.988185/
      data del/0.115,0.115,0.115/
      r(1)=ri(3)
      r(2)=ri(2)
      r(3)=ri(1)
      do 11 i=1,3
      ex=exp(-bet(i)*(r(i)-re(i)))
      q(i)=0.25e0*de(i)*((3.e0+del(i))*ex**2-(2.e0+6.e0*del(i))*ex)/
     1   (1.e0+del(i))
      a(i)=0.25e0*de(i)*((1.e0+3.e0*del(i))*ex**2-(6.e0+2.e0*del(i))
     1   *ex)/ (1.e0+del(i))
   11 continue
      u=q(1)+q(2)+q(3)-sqrt(a(1)**2+a(2)**2+a(3)**2-a(1)*a(2)-
     1   a(2)*a(3)-a(1)*a(3))
      v=(u+de(1))
      return
      end
      subroutine dderiv(f,n,x,iwrt,xa,xb,h,tol,maxit,iflag,fval,dfdx,
     &  errest,nitrs)
c
c**********************************************************************
c                                                                     *
c  Numerical Differentiation                                          *
c  Tom Allison, 28 October 1994                                       *
c    last modified 4 March 1996 - TCA                                 *
c                                                                     *
c**********************************************************************
c
c  DERIV is a subroutine which will return the numerical derivative of
c  a given (double precision) function at the selected point.  It will
c  automatically select either centered or one-sided derivatives based
c  on the input parameters.  The derivative subroutine uses three-point
c  derivative formulas and a Richardson's Extrapolation scheme.  An
c  estimate of the error, the number of iterations required to produce
c  the result, and a result flag are returned.
c
c  The formulas used in this code and an excellent explanation of the
c  methods used may be found in Sections 4.1 and 4.2 of "Numerical
c  Analysis" by Richard L. Burden and J. Douglas Faires, Fourth Edition,
c  PWS-Kent Publishing Company, Boston, 1989.
c
c  Some of the unique features of this code are its small memory requirements
c  due to the particular implementation of the Richardson's Extrapolation
c  (see notes below) and its ability to reuse function evaluations to reduce
c  the overall work required.  By reusing function evaluations, one-sided
c  differences require n+2 function evaluations, three-point central
c  differences require 2n+1 function evaluations, and five-point central
c  differences requite 2n+3 function evaluations, where n is the number
c  of iterations required to compute the derivative to the desired tolerance.
c
c  note:  the term one-sided differences used in the comments refers to the
c         use of a forward-difference formula if hs > 0, and a backward-
c         difference formula if hs < 0, where hs is the step size.
c
c  f       function (double precision) - input
c          f must be declared external in the calling program.
c          f must be of the form f(x) where x is a double precision vector
c            containing the point at which the derivative is to be evaluated.
c  n       scalar (integer) - input
c          n is the dimension of the function whose derivative is to be
c            evaluated.
c  x       vector (double precision) - input
c          x must be a vector of length n.
c          x is the point at which the derivative is to be evaluated.
c  iwrt    scalar (integer) - input
c          iwrt is the index of the variable [x(iwrt)] that the derivative
c            is being evaluated with respect to.
c  xa,xb   scalar (double precision) - input
c          bounds for the derivative approximation of the form [xa,xb]
c            subject to xa <= x(iwrt) <= xb, where x(iwrt) is the
c            variable the derivative is being evaluated with respect to.
c  h       scalar (double precision) - input/output
c          h is the initial step size used in the derivative calculation.
c            if h = 0.0 the routine will use a value of h that is equal
c            to half the size of the interval [xa,xb].  the absolute value
c            of h is used by the program.  on output, the value of h will
c            be the signed value of h used for the initial step size
c  tol     scalar (double precision)
c          tol is the tolerance that the routine will attempt to achieve
c            in its calculation of the derivative.  if the tolerance
c            criterion can not be met, the routine will return an error
c            flag (see iflag below).
c  maxit   scalar (integer) - input
c          maximum number of allowed iterations involved in evaluating
c            the derivative.  if the maximum number of iterations is
c            exceeded, the routine will return an error flag (see iflag
c            below).
c  iflag   scalar (integer) - output
c          iflag is used to return the status of the derivative evaluation.
c          iflag = 0  status ok, continue iteration (used internally)
c          iflag = 1  tolerance satisfied in derivative evaluation.
c          iflag = 2  maximum number of allowed iterations exceeded
c                     (the routine returns its best guess).
c          iflag = 3  no likely improvement.  the routines error estimate
c                     has exceeded the error estimate on the last iteration
c                     so the routine exits and returns its best guess.
c          iflag = 4  the number of variables in the function has exceeded
c                     an internal limit
c          iflag = 5  the index variable iwrt is greater than n
c          iflag = 6  the requested maximum number of iterations has exceeded
c                     an internal limit
c          iflag = 7  x(iwrt) is outside the range [xa,xb]
c          iflag = 8  the step size (h) is larger than the interval
c          iflag = 9  xa is greater than xb
c  fval    scalar (double precision) - output
c          value of the function at the point where the derivative is to
c          be evaluated.
c  dfdx    scalar (double precision) - output
c          dfdx is the numerical derivative of the function at the selected
c            point.
c  errest  scalar (double precision) - output
c          errest is an estimate of the error in the evaluation of the
c            derivative.
c  nitrs   scalar (integer) - output
c          number of iterations performed in evaluation of the derivative.
c
      implicit double precision(a-h,o-z)

      external f

      double precision f

c nvars is the maximum number of allowed variables in the function
      parameter(nvars=10)
c minit is the minimum number of iterations which must be performed
c before we give up and return a "best guess"
      parameter(minit=5)
c naits is the maximum number of allowed iterations (internal limit)
      parameter(naits=50)

      dimension x(n),p(nvars),dn(naits)

c assign value of step size
      if (h .eq. 0.d0) h=0.5d0*(xb-xa)
      hs=dabs(h)
      if (x(iwrt) .eq. xb) hs=-dabs(h)
      h=hs

c decide on one-sided (idtyp=1), central (idtyp=2) three-point,
c or central (idtyp=3) five-point differences
      idtyp=0
      xia=dabs(x(iwrt)-xa)
      xib=dabs(x(iwrt)-xb)
      if (x(iwrt) .eq. xa .or. x(iwrt) .eq. xb) then
        idtyp=1
      else
        if (xia .ge. dabs(2.d0*hs) .and. xib .ge. dabs(2.d0*hs)) then
          idtyp=3
        else if (xia .ge. dabs(hs) .and. xib .ge. dabs(hs)) then
          idtyp=2
        endif
      endif
c if the range [xa,xb] is too small on one side to do central differences,
c then attempt to use one-sided differences before returning an error
      if (idtyp .eq. 0) then
        if (xia .gt. xib) then
          xb=x(iwrt)
        else if (xib .ge. xia) then
          xa=x(iwrt)
        end if
        idtyp=1
      endif

c set initial values of input parameters in case there are errors
      iflag=0
      dfdx=0.d0
      errest=0.d0
      nitrs=0

c check input parameters and report errors if necessary
      if (n .gt. nvars) iflag=4
      if (iwrt .gt. n) iflag=5
      if (maxit .gt. naits) iflag=6
      if (x(iwrt) .lt. xa .or. x(iwrt) .gt. xb) iflag=7
      if (dabs(xb-xa) .lt. dabs(hs)) iflag=8
      if (xb .lt. xa) iflag=9

c if there are errors, jump to subroutine exit
      if (iflag .ne. 0) goto 100

c store the point (x) at which the derivative is to be taken
      do 5 i=1,n
5       p(i)=x(i)

c variables:
c md    - counts the iterations
c fterr - error estimate from the last iteration
c fx    - value of f(x)        one-sided derivative
c fxh   - value of f(x+hs)     one-sided derivative
c fx2h  - value of f(x+2*hs)   one-sided derivative
c fxph  - value of f(x+hs)     central three-point and five-point derivative
c fxmh  - value of f(x-hs)     central three-point and five-point derivative
c fxp2h - value of f(x+2*hs)   central five-point derivative
c fxm2h - value of f(x-2*hs)   central five-point derivative
c dn()  - storage for the extrapolation table
c hs    - step size used in the differencing scheme
c dcomp - value of the derivative computed in the last iteration
c dnew  - new value of derivative (to the order of the extrapolation)
c iflag - error flag.  iflag is equal to zero during iteration, and is set
c         to an appropriate exit code (1..3) when the computation is complete
c         also stores codes (4..9) for input errors
c x()   - point at which the function is to be evaluated
c terr  - current error estimate
c bn    - 4**(j-1) term which occurs in the Richardson extrapolation scheme
c dtemp - temporary storage for the newest element in the extrapolation table

c initialize variables and perform first evaluation of derivative
c first iteration
      md=1
c set initial error ridiculously high
      fterr=1.d+100

c evaluate function at the point where the derivative is to be taken
      fval=f(x)

c perform first evaluation of the derivative and save appropriate function
c evaluations for the second iteration
      if (idtyp .eq. 1) then
c one-sided differences, three-point formula
        fx=fval
        x(iwrt)=p(iwrt)+hs
        fxh=f(x)
        x(iwrt)=p(iwrt)+2.d0*hs
        fx2h=f(x)
        dn(1)=(0.5d0/hs)*(-3.d0*fx+4.d0*fxh-fx2h)
        fx2h=fxh
      else if (idtyp .eq. 2) then
c central differences, three-point formula
        x(iwrt)=p(iwrt)+hs
        fxph=f(x)
        x(iwrt)=p(iwrt)-hs
        fxmh=f(x)
        dn(1)=(0.5d0/hs)*(fxph-fxmh)
      else
c central differences, five-point formula
        x(iwrt)=p(iwrt)+hs
        fxph=f(x)
        x(iwrt)=p(iwrt)-hs
        fxmh=f(x)
        x(iwrt)=p(iwrt)+2.d0*hs
        fxp2h=f(x)
        x(iwrt)=p(iwrt)-2.d0*hs
        fxm2h=f(x)
        dn(1)=(1.d0/(12.d0*hs))*(fxm2h-8.d0*fxmh+8.d0*fxph-fxp2h)
        fxp2h=fxph
        fxm2h=fxmh
      end if

c update step size for next iteration
      hs=hs*0.5d0

c loop for further iterations
10    continue

c a little houskeeping
        dcomp=dn(md)
        md=md+1
        if (md .ge. maxit) iflag=2

c evaluate the function at the appropriate points and compute the new value
c of the derivative (dnew).  save appropriate function evaluations for use
c on the next interation
        if (idtyp .eq. 1) then
c one-sided differences, three-point formula
          x(iwrt)=p(iwrt)+hs
          fxh=f(x)
          dnew=(0.5d0/hs)*(-3.d0*fx+4.d0*fxh-fx2h)
          fx2h=fxh
        else if (idtyp .eq. 2) then
c central differences, three-point formula
          x(iwrt)=p(iwrt)+hs
          fxph=f(x)
          x(iwrt)=p(iwrt)-hs
          fxmh=f(x)
          dnew=(0.5d0/hs)*(fxph-fxmh)
        else
c central differences, five-point formula
          x(iwrt)=p(iwrt)+hs
          fxph=f(x)
          x(iwrt)=p(iwrt)-hs
          fxmh=f(x)
          dnew=(1.d0/(12.d0*hs))*(fxm2h-8.d0*fxmh+8.d0*fxph-fxp2h)
          fxp2h=fxph
          fxm2h=fxmh
        endif

c update the step size for the next iteration
        hs=hs*0.5d0

c Richardson extrapolation
c   the next few lines of code update the extrapolation table
c   this is a slight modification of "standard" Richardson extrapolation
c   only the current row of the extrapolation table is saved since the
c   previous rows are no longer needed.  this is done to save memory.
c   on exit the new value of the derivative is dnew and is stored in dn(md)
        bn=4.d0
        do 20 i=2,md
          dtemp=(bn*dnew-dn(i-1))/(bn-1.d0)
          dn(i-1)=dnew
          dnew=dtemp
          bn=bn*4.d0
20      continue
        dn(md)=dnew

c compute a new error estimate based on the current and last computed
c value of the derivative
        terr=dabs(dcomp-dn(md))
c if the error estimate is larger on this iteration than it was on
c the last, set an exit flag, but only if we have performed more than
c the minimum number of iterations
        if (terr .gt. fterr .and. md .ge. minit) iflag=3
c if the current error estimate is less than the tolerance, set an exit flag
        if (terr .le. tol) iflag=1
c save the error estimate for this iteration if we intend to perform
c another iteration
        if (iflag .eq. 0) fterr=terr

c if iflag .eq. 0 then perform another iteration
      if (iflag .eq. 0) goto 10

c if iflag .eq. 3 (error estimate increased on the last iteration) then
c we can only return a "best guess"
      if (iflag .eq. 3) dfdx=dcomp
c if the tolerance criterion has been met (iflag .eq. 1) then return the
c computed value of the derivative.  if the maximum number of iterations
c have been performed, return our best value so far.
      if (iflag .eq. 1 .or. iflag .eq. 2) dfdx=dn(md)
c return the number of iterations performed
      nitrs=md
c return the error estimate
      errest=terr

100   continue

c return the original value of x
      do 110 i=1,n
110     x(i)=p(i)

      return

      end
