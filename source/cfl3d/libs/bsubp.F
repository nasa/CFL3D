      subroutine bsubp(id1,id2,a,b,c,f,i1,i2,il,iu,g,h)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Performs the back substitution for a block 5x5 tridi-
c     agonal matrix equation (periodic) solution.  The vectorization is
c     over points i1-i2 and the tridiagonal matrix equation spans points
c     il-iu.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(id1,id2,5,5),b(id1,id2,5,5),c(id1,id2,5,5)
      dimension g(id1,id2,5,5),h(id1,id2,5,5)
      dimension f(id1,id2,5)
c
c     solve matrix equation
c
      il1 = il+1
      is  = il
c
c      f=binv*f
c
cdir$ ivdep
      do 100 i=i1,i2
      t1        = b(i,is,1,1)*(f(i,is,1))
      t2        = b(i,is,2,2)*(f(i,is,2)-b(i,is,2,1)*t1)
      t3        = b(i,is,3,3)*(f(i,is,3)-b(i,is,3,1)*t1-b(i,is,3,2)*t2)
      t4        = b(i,is,4,4)*(f(i,is,4)-b(i,is,4,1)*t1-b(i,is,4,2)*t2
     .                                  -b(i,is,4,3)*t3)
      f(i,is,5) = b(i,is,5,5)*(f(i,is,5)-b(i,is,5,1)*t1-b(i,is,5,2)*t2
     .                                -b(i,is,5,3)*t3-b(i,is,5,4)*t4)
      f(i,is,4) = t4-b(i,is,4,5)*f(i,is,5)
      f(i,is,3) = t3-b(i,is,3,5)*f(i,is,5)-b(i,is,3,4)*f(i,is,4)
      f(i,is,2) = t2-b(i,is,2,5)*f(i,is,5)-b(i,is,2,4)*f(i,is,4)
     .              -b(i,is,2,3)*f(i,is,3)
      f(i,is,1) = t1-b(i,is,1,5)*f(i,is,5)-b(i,is,1,4)*f(i,is,4)
     .              -b(i,is,1,3)*f(i,is,3)-b(i,is,1,2)*f(i,is,2)
  100 continue
c
c      forward sweep
c 
      iux = iu-1
      do 200 is=il1,iux
      ir  = is-1
      it  = is+1
c      first row reduction
      do 110 m=1,5
cdir$ ivdep
      do 102 i=i1,i2
      f(i,is,m)= f(i,is,m)-a(i,is,m,1)*f(i,ir,1)
     .                    -a(i,is,m,2)*f(i,ir,2)
     .                    -a(i,is,m,3)*f(i,ir,3)
     .                    -a(i,is,m,4)*f(i,ir,4)
     .                    -a(i,is,m,5)*f(i,ir,5)
  102 continue
  110 continue
c
c      f=binv*f
c
cdir$ ivdep
      do 300 i=i1,i2
      t1        = b(i,is,1,1)*(f(i,is,1))
      t2        = b(i,is,2,2)*(f(i,is,2)-b(i,is,2,1)*t1)
      t3        = b(i,is,3,3)*(f(i,is,3)-b(i,is,3,1)*t1-b(i,is,3,2)*t2)
      t4        = b(i,is,4,4)*(f(i,is,4)-b(i,is,4,1)*t1-b(i,is,4,2)*t2
     .                                  -b(i,is,4,3)*t3)
      f(i,is,5) = b(i,is,5,5)*(f(i,is,5)-b(i,is,5,1)*t1-b(i,is,5,2)*t2
     .                                  -b(i,is,5,3)*t3-b(i,is,5,4)*t4)
      f(i,is,4) = t4-b(i,is,4,5)*f(i,is,5)
      f(i,is,3) = t3-b(i,is,3,5)*f(i,is,5)-b(i,is,3,4)*f(i,is,4)
      f(i,is,2) = t2-b(i,is,2,5)*f(i,is,5)-b(i,is,2,4)*f(i,is,4)
     .              -b(i,is,2,3)*f(i,is,3)
      f(i,is,1) = t1-b(i,is,1,5)*f(i,is,5)-b(i,is,1,4)*f(i,is,4)
     .              -b(i,is,1,3)*f(i,is,3)-b(i,is,1,2)*f(i,is,2)
  300 continue
  200 continue
c
      is = iu
      ir = is-1
      it = is+1
c      first row reduction
      do 500 m=1,5
cdir$ ivdep
      do 490 i=i1,i2
      f(i,is,m) = f(i,is,m)-a(i,is,m,1)*f(i,ir,1)
     .                     -a(i,is,m,2)*f(i,ir,2)
     .                     -a(i,is,m,3)*f(i,ir,3)
     .                     -a(i,is,m,4)*f(i,ir,4)
     .                     -a(i,is,m,5)*f(i,ir,5)
  490 continue
  500 continue
      iu2 = iu-2
      do 600 ix=il,iu2
      do 590 m=1,5
cdir$ ivdep
      do 580 i=i1,i2
      f(i,is,m) = f(i,is,m)-h(i,ix,m,1)*f(i,ix,1)
     .                     -h(i,ix,m,2)*f(i,ix,2)
     .                     -h(i,ix,m,3)*f(i,ix,3)
     .                     -h(i,ix,m,4)*f(i,ix,4)
     .                     -h(i,ix,m,5)*f(i,ix,5)
  580 continue
  590 continue
  600 continue
c
c      f=binv*f
c
cdir$ ivdep
      do 700 i=i1,i2
      t1        = b(i,is,1,1)*(f(i,is,1))
      t2        = b(i,is,2,2)*(f(i,is,2)-b(i,is,2,1)*t1)
      t3        = b(i,is,3,3)*(f(i,is,3)-b(i,is,3,1)*t1-b(i,is,3,2)*t2)
      t4        = b(i,is,4,4)*(f(i,is,4)-b(i,is,4,1)*t1-b(i,is,4,2)*t2
     .                                  -b(i,is,4,3)*t3)
      f(i,is,5) = b(i,is,5,5)*(f(i,is,5)-b(i,is,5,1)*t1-b(i,is,5,2)*t2
     .                                  -b(i,is,5,3)*t3-b(i,is,5,4)*t4)
      f(i,is,4) = t4-b(i,is,4,5)*f(i,is,5)
      f(i,is,3) = t3-b(i,is,3,5)*f(i,is,5)-b(i,is,3,4)*f(i,is,4)
      f(i,is,2) = t2-b(i,is,2,5)*f(i,is,5)-b(i,is,2,4)*f(i,is,4)
     .              -b(i,is,2,3)*f(i,is,3)
      f(i,is,1) = t1-b(i,is,1,5)*f(i,is,5)-b(i,is,1,4)*f(i,is,4)
     .              -b(i,is,1,3)*f(i,is,3)-b(i,is,1,2)*f(i,is,2)
  700 continue
c
c      back substitution
c
      iux = il1
      do 310 iqq=il1,iux
      is  = il+iu-iqq
      it  = is+1
      do 307 m=1,5
cdir$ ivdep
      do 305 i=i1,i2
      f(i,is,m) = f(i,is,m)-c(i,is,m,1)*f(i,it,1)
     .                     -c(i,is,m,2)*f(i,it,2)
     .                     -c(i,is,m,3)*f(i,it,3)
     .                     -c(i,is,m,4)*f(i,it,4)
     .                     -c(i,is,m,5)*f(i,it,5)
  305 continue
  307 continue
  310 continue
c
      il11 = il1+1
      do 900 iqq=il11,iu
      is   = il+iu-iqq
      it   = is+1
      do 840 m=1,5
cdir$ ivdep
      do 820 i=i1,i2
      f(i,is,m) = f(i,is,m)-c(i,is,m,1)*f(i,it,1)
     .                     -c(i,is,m,2)*f(i,it,2)
     .                     -c(i,is,m,3)*f(i,it,3)
     .                     -c(i,is,m,4)*f(i,it,4)
     .                     -c(i,is,m,5)*f(i,it,5)-g(i,is,m,1)*f(i,iu,1)
     .                                           -g(i,is,m,2)*f(i,iu,2)
     .                                           -g(i,is,m,3)*f(i,iu,3)
     .                                           -g(i,is,m,4)*f(i,iu,4)
     .                                           -g(i,is,m,5)*f(i,iu,5)
  820 continue
  840 continue
  900 continue
      return
      end
