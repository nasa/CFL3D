      subroutine dlutrp(nvmax,n,nmax,il,iu,a,b,c,g,h)
c***********************************************************************
c     Purpose:  Perform the LU decomposition for scalar a periodic 
c     tridiagonal system of equations.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(nvmax,nmax),b(nvmax,nmax),c(nvmax,nmax)
      dimension g(nvmax,nmax),h(nvmax,nmax)
c
c      periodic solver  5-3-85
c
c      inversion of block tridiagonal...a,b,c are scalars
c      f is forcing function and solution is output in f
c      solution is by upper triangularization with unity diagonal
c      block inversions use nonpivoted lu decomposition
c      il and iu are starting and finishing indices
c      b,c,and e are overloaded
c
c      set g and h
c
cdir$ ivdep
      do 1000 izz=1,n
      g(izz,il) = a(izz,il)
      h(izz,il) = c(izz,iu)
 1000 continue
c
      il1 = il+1
      i   = il
c
c      l-u decomposition
c
cdir$ ivdep
      do 1001 izz=1,n
      b(izz,i) = 1.e0/b(izz,i)
 1001 continue
c
      if (i.eq.iu) go to 1030
c
c      c=binv*c
c      g=binv*g
c
cdir$ ivdep
      do 1002 izz=1,n
      c(izz,i) = b(izz,i)*(c(izz,i))
      g(izz,i) = b(izz,i)*(g(izz,i))
 1002 continue
 1030 continue
c      forward sweep
      iu2 = iu-2
      if (iu2.eq.il) go to 101
      do 100 i=il1,iu2
      ir = i-1
      it = i+1
c      row reduction
cdir$ ivdep
      do 1003 izz=1,n
      b(izz,i)  =  b(izz,i)-a(izz,i)*c(izz,ir)
      g(izz,i)  = -a(izz,i)*g(izz,ir)
      h(izz,i)  = -h(izz,ir)*c(izz,ir)
      b(izz,iu) =  b(izz,iu)-h(izz,ir)*g(izz,ir)
c
c      l-u decomposition
c
      b(izz,i) = 1.e0/b(izz,i)
c
c      c=binv*c
c
      c(izz,i) = b(izz,i)*(c(izz,i))
      g(izz,i) = b(izz,i)*(g(izz,i))
 1003 continue
 1130 continue
  100 continue
  101 continue
c
      i  = iu-1
      ir = i-1
      it = i+1
c      row reduction
cdir$ ivdep
      do 1004 izz=1,n
      b(izz,i)  = b(izz,i)-a(izz,i)*c(izz,ir)
      c(izz,i)  = c(izz,i)-a(izz,i)*g(izz,ir)
c************************** 9-29-86
      a(izz,iu) = a(izz,iu)-h(izz,ir)*c(izz,ir)
      b(izz,iu) = b(izz,iu)-h(izz,ir)*g(izz,ir)
c
c      l-u decomposition
c
      b(izz,i) = 1.e0/b(izz,i)
c
c      c=binv*c
c
      c(izz,i) = b(izz,i)*(c(izz,i))
 1004 continue
      i  = iu
      ir = i-1
      it = i+1
c      row reduction
cdir$ ivdep
      do 1005 izz=1,n
      b(izz,i) = b(izz,i)-a(izz,i)*c(izz,ir)
c
c      l-u decomposition
c
      b(izz,i) = 1.e0/b(izz,i)
 1005 continue
      return
      end
