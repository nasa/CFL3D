      subroutine initnonin(nbl,jdim,kdim,idim,q,qj0,qk0,qi0,
     .                tj0,tk0,ti0,
     .                vol,volj0,volk0,voli0,
     .                x,y,z)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Increment the initial conditions on a mesh to include
c     the rotational component for noninertial calculations.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension tj0(kdim,idim-1,2,4),tk0(jdim,idim-1,2,4),
     .          ti0(jdim,kdim  ,2,4)
      dimension volj0(kdim,idim-1,4),volk0(jdim,idim-1,4),
     .          voli0(jdim,kdim  ,4),vol(jdim,kdim,idim-1)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
c
      wx = xrotrate
      wy = yrotrate
      wz = zrotrate
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      do 1000 j=1,jdim1
        do 1000 k=1,kdim1
          do 1000 i=1,idim1
 
c    compute cell centers 
 
      cx = 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      cy = 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      cz = 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      cx = cx - xcentrot
      cy = cy - ycentrot
      cz = cz - zcentrot
 
c     increase velocity with rigid body rotaion component:       
 
c     Uinf + r x omega is the same as Uinf - omega x r 
 
      q(j,k,i,2) = q(j,k,i,2) + ( cy * wz - cz * wy )
      q(j,k,i,3) = q(j,k,i,3) + ( cz * wx - cx * wz )
      q(j,k,i,4) = q(j,k,i,4) + ( cx * wy - cy * wx )
 
 1000 continue
c
      return
      end
