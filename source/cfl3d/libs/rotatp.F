      subroutine rotatp(mdim,ndim,jmax1,kmax1,msub1,l,x1,y1,z1,
     .                  dthetx,dthety,dthetz,xorig,yorig,zorig,
     .                  mbl,nn,intmx,int)
c***********************************************************************
c     Purpose:  Rotate "from" block to provide complete coverage for
c     interpolation for cases in which the complete physical domain is
c     not modeled.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x1(mdim,ndim,msub1),y1(mdim,ndim,msub1),
     .          z1(mdim,ndim,msub1),dthetx(intmx,msub1),
     .          dthety(intmx,msub1),dthetz(intmx,msub1),
     .          xorig(nn),yorig(nn),zorig(nn)
c
      xorg = xorig(mbl)
      yorg = yorig(mbl)
      zorg = zorig(mbl)
      dthx = dthetx(int,l)
      dthy = dthety(int,l)
      dthz = dthetz(int,l)
c
      if (abs(real(dthx)) .gt. 0.) then
c 
c        rotation about a line parallel to the x-axis
c
         ca = cos(dthx)
         sa = sin(dthx)
         do 10 j=1,jmax1
         do 10 k=1,kmax1
         ytemp = (y1(j,k,l) - yorg)*ca - (z1(j,k,l) - zorg)*sa + yorg
         ztemp = (y1(j,k,l) - yorg)*sa + (z1(j,k,l) - zorg)*ca + zorg
         y1(j,k,l) = ytemp
         z1(j,k,l) = ztemp
 10      continue
      end if
c
      if (abs(real(dthy)) .gt. 0.) then
c
c        rotation about a line parallel to the y-axis
c
         ca = cos(dthy)
         sa = sin(dthy)
         do 20 j=1,jmax1
         do 20 k=1,kmax1
         xtemp =  (x1(j,k,l) - xorg)*ca + (z1(j,k,l) - zorg)*sa + xorg
         ztemp = -(x1(j,k,l) - xorg)*sa + (z1(j,k,l) - zorg)*ca + zorg
         x1(j,k,l) = xtemp
         z1(j,k,l) = ztemp
 20      continue
      end if
c
      if (abs(real(dthz)) .gt. 0.) then
c
c        rotation about a line parallel to the z-axis
c
         ca = cos(dthz)
         sa = sin(dthx)
         do 30 j=1,jmax1
         do 30 k=1,kmax1
         xtemp = (x1(j,k,l) - xorg)*ca - (y1(j,k,l) - yorg)*sa + xorg
         ytemp = (x1(j,k,l) - xorg)*sa + (y1(j,k,l) - yorg)*ca + yorg
         x1(j,k,l) = xtemp
         y1(j,k,l) = ytemp
 30      continue
      end if
c
      return
      end
