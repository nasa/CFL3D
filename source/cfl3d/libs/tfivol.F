      subroutine tfivol(idim,jdim,kdim,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,
     .                  arci,arcj,arck,nou,bou,nbuf,ibufdim,myid)
c***********************************************************************
c     Purpose: compute transfinite interpolation in a 3d volume, using
c              arc-length blending functions
c*********************************************************************** 
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension x1(jdim,kdim,idim),y1(jdim,kdim,idim),z1(jdim,kdim,idim)
      dimension x2(jdim,kdim,idim),y2(jdim,kdim,idim),z2(jdim,kdim,idim)
      dimension x3(jdim,kdim,idim),y3(jdim,kdim,idim),z3(jdim,kdim,idim)
      dimension arci(jdim,kdim,idim),arcj(jdim,kdim,idim),
     .          arck(jdim,kdim,idim)
c 
      dimension phi(2),psi(2),omg(2)
c
c     eta projector  -  1st step 
c
      do j=1,jdim
         do k=1,kdim
            do i=1,idim
               eta = arcj(j,k,i)
               psi(1) = eta
               psi(2) = 1.-eta
               x1(j,k,i) = psi(2)*x(1,k,i) 
     .                   + psi(1)*x(jdim,k,i)
               y1(j,k,i) = psi(2)*y(1,k,i) 
     .                   + psi(1)*y(jdim,k,i)
               z1(j,k,i) = psi(2)*z(1,k,i) 
     .                   + psi(1)*z(jdim,k,i)
            end do
         end do
      end do
c 
c     zeta projector  -  2nd step
c
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               zeta = arck(j,k,i)
               omg(1) = zeta
               omg(2) = 1.-zeta
               x2(j,k,i) = omg(2)*(x(j,1,i)    - x1(j,1,i))
     .                   + omg(1)*(x(j,kdim,i) - x1(j,kdim,i))
               y2(j,k,i) = omg(2)*(y(j,1,i)    - y1(j,1,i))
     .                   + omg(1)*(y(j,kdim,i) - y1(j,kdim,i))
               z2(j,k,i) = omg(2)*(z(j,1,i)    - z1(j,1,i))
     .                   + omg(1)*(z(j,kdim,i) - z1(j,kdim,i))
            end do
         end do
      end do
c 
c     xi projector  -  3rd step
c
      do i=1,idim
         do j=1,jdim
            do k=1,kdim
               xi = arci(j,k,i)
               phi(1) = xi
               phi(2) = 1.-xi
               x3(j,k,i) = phi(2)*(x(j,k,1)    - x1(j,k,1)
     .                   - x2(j,k,1))
     .                   + phi(1)*(x(j,k,idim) - x1(j,k,idim)
     .                   - x2(j,k,idim))
               y3(j,k,i) = phi(2)*(y(j,k,1)    - y1(j,k,1)
     .                   - y2(j,k,1))
     .                   + phi(1)*(y(j,k,idim) - y1(j,k,idim)
     .                   - y2(j,k,idim))
               z3(j,k,i) = phi(2)*(z(j,k,1)    - z1(j,k,1)
     .                   - z2(j,k,1))
     .                   + phi(1)*(z(j,k,idim) - z1(j,k,idim)
     .                   - z2(j,k,idim))
            end do
         end do
      end do
c
c     summation of xi, eta, and zeta projectors  -  4th step
c
      do i=1,idim
         do j=1,jdim
            do k=1,kdim
               x(j,k,i)  = x1(j,k,i) + x2(j,k,i) + x3(j,k,i)
               y(j,k,i)  = y1(j,k,i) + y2(j,k,i) + y3(j,k,i)
               z(j,k,i)  = z1(j,k,i) + z2(j,k,i) + z3(j,k,i)
            end do
         end do
      end do
c
      return
      end 
