      subroutine dfbtr(nvmax,n,nmax,il,iu,a,b,c,f)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Perform the back substitution for a scalar
c     tridiagonal system of equations.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(nvmax,nmax),b(nvmax,nmax),c(nvmax,nmax),f(nvmax,nmax)
c
      il1 = il+1
      is  = il
c
c      f=binv*f
c
cdir$ ivdep
      do 1000 izz=1,n
      f(izz,is) = b(izz,is)*f(izz,is)
 1000 continue
c
c      forward sweep
c
      do 9100 is=il1,iu
      ir = is-1
      it = is+1
c      first row reduction
cdir$ ivdep
      do 1001 izz=1,n
c
c      f=binv*f
c
      f(izz,is) = b(izz,is)*(f(izz,is)-a(izz,is)*f(izz,ir)) 
 1001 continue
 9100 continue
c      back substitution
      do 9180 ii=il1,iu
      is = il+iu-ii
      it = is+1
cdir$ ivdep
      do 1002 izz=1,n
      f(izz,is) = f(izz,is)-c(izz,is)*f(izz,it)
 1002 continue
 9180 continue
      return
      end
