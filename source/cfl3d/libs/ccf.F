      subroutine ccf(x,y,ca,sa,cl,xm,uf,wf,cf,pi)
c***********************************************************************
c     Purpose:  Modify uf, wf, and cf (velocities and speed of sound) at
c     the farfield, based on point vortex correction (used when i2d=-1 and
c     farfield bctype=1003).
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
c   far field boundary condition (point vortex)
c
      xe= ca*x+sa*y
      ye=-sa*x+ca*y
      beta=sqrt(1.e0-xm*xm)
      re=xe**2+(ye*beta)**2
      fact=cl*beta*xm/(re*4.e0*pi)
      ue=xm+fact*ye
      ve=  -fact*xe
      uf= ca*ue - sa*ve
      wf= sa*ue + ca*ve
      cf=sqrt(1.e0+0.2e0*(xm*xm-uf*uf-wf*wf))
      return
      end
