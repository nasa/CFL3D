      subroutine direct(x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,
     .                   a1,a2,a3,imaxa,nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute (normalized) directed area components, 
c     or equivalently, components of unit normal to cell face
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      xxie = x8-x7
      yxie = y8-y7
      zxie = z8-z7
      xeta = x6-x5
      yeta = y6-y5
      zeta = z6-z5
      a1 = yxie*zeta-zxie*yeta
      a2 = zxie*xeta-xxie*zeta
      a3 = xxie*yeta-yxie*xeta
      d  = sqrt(a1*a1+a2*a2+a3*a3)
      if(real(d) .le. 0.) then
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),*) ' WARNINING: a cell with zero area has ',
     .             'been detected in subroutine direct'
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),*) '  - severe problem'
        d=1.
      end if
      a1 = a1/d
      a2 = a2/d
      a3 = a3/d
c
c     find coordinate direction with maximum area component
c
      amax  = ccabs(a1)
      imaxa = 1
      if (abs(real(a2)).gt.real(amax)) then
         amax  = ccabs(a2)
         imaxa = 2
      end if
      if (abs(real(a3)).gt.real(amax)) then
         amax  = ccabs(a3)
         imaxa = 3
      end if
      return
      end      
