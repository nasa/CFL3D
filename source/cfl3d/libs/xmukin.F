      subroutine xmukin(n,temp,visc,tinf)
c***********************************************************************
c     Purpose:  Computes Sutherland's formula.  Note that this routine
c     is only called by ctime.  (Mu is generally computed in-line
c     throughout the code using this formula.)
c     Also note that an older version of this routine used a linear
c     law for low temperatures, but it was not consistently used
c     throughout the whole code, and so has been taken out.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension temp(n),visc(n)
c
      common /fluid2/ pr,prt,cbar
c
c      molecular viscosity
c
      c2b    = cbar/tinf
      c2bp   = c2b+1.e0
c
      do 1000 j=1,n
         t5      =  temp(j)
         t6      =  sqrt(t5)
         visc(j) =  c2bp*t5*t6/(c2b+t5)
 1000 continue
      return
      end
