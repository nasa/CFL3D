      subroutine triv(jdim,kdim,jl,ju,kl,ku,x,a,b,c,f)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Solve a scalar tridiagonal system of equations
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(jdim,kdim),b(jdim,kdim),c(jdim,kdim)
      dimension x(jdim,kdim),f(jdim,kdim)
c
      do 10 j=jl,ju
      x(j,kl)=c(j,kl)/b(j,kl)
      f(j,kl)=f(j,kl)/b(j,kl)
 10   continue
      klp1 = kl +1
      do 1 i=klp1,ku
      do 20 j=jl,ju
         z=1./(b(j,i)-a(j,i)*x(j,i-1))
         x(j,i)=c(j,i)*z
         f(j,i)=(f(j,i)-a(j,i)*f(j,i-1))*z
 20   continue
1     continue
c
      kupkl=ku+kl
      do 2 i1=klp1,ku
         i=kupkl-i1
         do 30 j=jl,ju
         f(j,i)=f(j,i)-x(j,i)*f(j,i+1)
 30      continue
2     continue
c
      return
      end
