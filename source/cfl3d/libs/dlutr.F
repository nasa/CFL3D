      subroutine dlutr(nvmax,n,nmax,il,iu,a,b,c)
c***********************************************************************
c     Purpose:  Perform the scalar tridiagonal (LU) decomposition.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(nvmax,nmax),b(nvmax,nmax),c(nvmax,nmax)
c
c      inversion of scalar tridiagonal...a,b,c are scalars
c      f is forcing function and solution is output in f
c      solution is by upper triangularization with unity diagonal
c      block inversions use nonpivoted lu decomposition
c      il and iu are starting and finishing indices
c      b and c are overloaded
c
      il1 = il+1
      i   = il
c
c      l-u decomposition
c
cdir$ ivdep
      do 1000 izz=1,n
      b(izz,i) = 1.e0/b(izz,i)
 1000 continue
c
      if (i.eq.iu) go to 1030
c
c      c=ainv*c
c
cdir$ ivdep
      do 1001 izz=1,n
      c(izz,i) = b(izz,i)*c(izz,i)
 1001 continue
 1030 continue
c      forward sweep
      do 100 i=il1,iu
      ir = i-1
      it = i+1
c      first row reduction
cdir$ ivdep
      do 1002 izz=1,n
      b(izz,i) = b(izz,i)-a(izz,i)*c(izz,ir)
c
c      l-u decomposition
c
      b(izz,i) = 1.e0/b(izz,i)
 1002 continue
      if (i.eq.iu) go to 100
c
c      c=ainv*c
c
cdir$ ivdep
      do 1003 izz=1,n
      c(izz,i) = b(izz,i)*c(izz,i)
 1003 continue
  100 continue
      return
      end
