      subroutine rotatmc(irot,rfreqr,omegx,omegy,omegz,xorg,yorg,zorg,
     .                thetax,thetay,thetaz,xmc,ymc,zmc,iupdat,time2)
c***********************************************************************
c     Purpose:  Determines increment to moment center due to rotation
c
c     irot....modulation for rotational motion
c             = 0 no rotation
c             = 1 constant rotation speed
c             = 2 sinusoidal variation of rotational displacement
c             = 3 smooth increase in rotational displacement,               
c                 asypmtotically reaching a fixed rotational displacement
c
c     iupdat..flag to update moment center position
c             = 0 don't update position
c             > 0 update position
c
c     rotations are taken with positive angular displacement following
c     the right hand rule
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
c     ft modulates the rotation
c
      if (irot .eq. 0)  then
         return
      else if (irot .eq. 1)  then
         ft     = time2
      else if (irot .eq. 2)  then
         ft     = sin(rfreqr*time2)
      else if (irot .eq. 3)  then
         expt   = exp(-rfreqr*time2)
         ft     = 1.-expt
      end if
c
      if (abs(real(omegx)) .gt. 0.0) then

c***********************************************************************
c        rotate about an axis parallel to the x-axis
c***********************************************************************
c
         if (time2 .ne. 0.) then
            theold = thetax
         else
            theold = 0.e0
         end if
c
         theta    = omegx*ft
         dtheta   = theta - theold
         ca = cos(dtheta)
         sa = sin(dtheta)
c
         if (iupdat .gt. 0) then
            thetax = theta
            yml = ymc
            zml = zmc
            xmc = xmc
            ymc = (yml-yorg)*ca-(zml-zorg)*sa+yorg
            zmc = (yml-yorg)*sa+(zml-zorg)*ca+zorg
         end if
c
      else if (abs(real(omegy)) .gt. 0.0) then
c
c***********************************************************************
c        rotate about an axis parallel to the y-axis
c***********************************************************************
c
         if (time2 .ne. 0.) then
            theold = thetay
         else
            theold = 0.e0
         end if
c
         theta    = omegy*ft
         dtheta   = theta - theold
         ca = cos(dtheta)
         sa = sin(dtheta)
c
         if (iupdat .gt. 0) then
            thetay = theta
            xml = xmc
            zml = zmc
            xmc =  (xml-xorg)*ca+(zml-zorg)*sa+xorg
            ymc =   ymc
            zmc = -(xml-xorg)*sa+(zml-zorg)*ca+zorg
         end if
c
      else if (abs(real(omegz)) .gt. 0.0) then
c
c***********************************************************************
c        rotate about an axis parallel to the z-axis
c***********************************************************************
c
         if (time2 .ne. 0.) then
            theold = thetaz
         else
            theold = 0.e0
         end if
c
         theta    = omegz*ft
         dtheta   = theta - theold
         ca = cos(dtheta)
         sa = sin(dtheta)
c
         if (iupdat .gt. 0) then
            thetaz = theta
            xml = xmc
            yml = ymc
            xmc = (xml-xorg)*ca-(yml-yorg)*sa+xorg
            ymc = (xml-xorg)*sa+(yml-yorg)*ca+yorg
            zmc = zmc
         end if
c
      end if
c
      return
      end
