      subroutine project(x0,y0,z0,xc,yc,zc,a1,a2,a3,xcp,ycp,zcp)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Project point xc,yc,zc into plane containing the point
c     x0,y0,z0 and which has the unit normal with components a1,a2,a3
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      a=(xc-x0)*a1+(yc-y0)*a2+(zc-z0)*a3
      xcp = xc-a*a1
      ycp = yc-a*a2
      zcp = zc-a*a3
      return
      end
