      subroutine grdmove(nbl,jdim,kdim,idim,x,y,z,xorig,yorig,
     .                  zorig,xorg,yorg,zorg,thetax,thetay,
     .                  thetaz)
c***********************************************************************
c     Purpose:  Move the grid from one position to another. On input,
c     (x,y,z) are the grid coordinates in the old position, with origin
c     (xorig,yorig,zorig). Given the new origin (xorg,yorg,zorg) and
c     rotational displacements (thetax,thetay,thetaz) the grid is first
c     rotated and then translated to the new position; on return,
c     (x,y,z) correspond to the new position. No translation occurs if
c     (xorg,yorg,zorg) = (xorig,yorig,zorig); no rotation occurs if
c     (thetax,thetay,thetaz) = (0,0,0)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c
c     rotation about x-axis
      sa = sin(thetax)
      ca = cos(thetax)
      do 20 j=1,jdim
      do 20 k=1,kdim
      do 20 i=1,idim
      yy = (y(j,k,i)-yorig)*ca-(z(j,k,i)-zorig)*sa+yorig
      zz = (y(j,k,i)-yorig)*sa+(z(j,k,i)-zorig)*ca+zorig
      y(j,k,i) = yy
      z(j,k,i) = zz
   20 continue
c     
c     rotation about y-axis
      sa = sin(thetay)
      ca = cos(thetay)
      do 30 j=1,jdim
      do 30 k=1,kdim
      do 30 i=1,idim
      xx =  (x(j,k,i)-xorig)*ca+(z(j,k,i)-zorig)*sa+xorig
      zz = -(x(j,k,i)-xorig)*sa+(z(j,k,i)-zorig)*ca+zorig
      x(j,k,i) = xx
      z(j,k,i) = zz
   30 continue
c     
c     rotation about z-axis
      sa = sin(thetaz)
      ca = cos(thetaz)
      do 40 j=1,jdim
      do 40 k=1,kdim
      do 40 i=1,idim
      xx = (x(j,k,i)-xorig)*ca-(y(j,k,i)-yorig)*sa+xorig
      yy = (x(j,k,i)-xorig)*sa+(y(j,k,i)-yorig)*ca+yorig
      x(j,k,i) = xx
      y(j,k,i) = yy
   40 continue
c
c     translation
      do 10 j=1,jdim
      do 10 k=1,kdim
      do 10 i=1,idim
      x(j,k,i) = x(j,k,i) + xorg - xorig
      y(j,k,i) = y(j,k,i) + yorg - yorig
      z(j,k,i) = z(j,k,i) + zorg - zorig
  10  continue
c
      return
      end
