      subroutine transp(mdim,ndim,jmax1,kmax1,msub1,l,x1,y1,z1,
     .                  dx,dy,dz,intmx,int)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Translate "from" block to provide complete coverage for
c     interpolation for cases in which the complete physical domain is 
c     not modeled.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x1(mdim,ndim,msub1),y1(mdim,ndim,msub1),
     .          z1(mdim,ndim,msub1),dx(intmx,msub1),dy(intmx,msub1),
     .          dz(intmx,msub1)
c
      do 10 j=1,jmax1
      do 10 k=1,kmax1
      x1(j,k,l) = x1(j,k,l) + dx(int,l)
      y1(j,k,l) = y1(j,k,l) + dy(int,l)
      z1(j,k,l) = z1(j,k,l) + dz(int,l)
 10   continue
      return
      end
