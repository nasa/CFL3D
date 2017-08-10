c  ---------------------------------------------------------------------------
c  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
c  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
c  with point-matched, patched, overset, or embedded connectivities. Both
c  multigrid and mesh sequencing are available in time-accurate or
c  steady-state modes.
c
c  Copyright 2001 United States Government as represented by the Administrator
c  of the National Aeronautics and Space Administration. All Rights Reserved.
c 
c  The CFL3D platform is licensed under the Apache License, Version 2.0 
c  (the "License"); you may not use this file except in compliance with the 
c  License. You may obtain a copy of the License at 
c  http://www.apache.org/licenses/LICENSE-2.0. 
c 
c  Unless required by applicable law or agreed to in writing, software 
c  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
c  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
c  License for the specific language governing permissions and limitations 
c  under the License.
c  ---------------------------------------------------------------------------
c
      subroutine tfiface(idim,jdim,kdim,x,y,z,x1,y1,z1,x2,y2,z2,
     .                   arci,arcj,arck,i1,i2,j1,j2,k1,k2,nou,bou,nbuf,
     .                   ibufdim,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose: compute transfinite interpolation on block faces, using
c              arc-length blending functions
c
c     this subroutine expects one and only pair of (i1,i2), (j1,j2), or
c     (k1,k2) to contain identical indicies; the identical indicies
c     determine which face the TFI is carried out on.
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
      dimension arci(jdim,kdim,idim),arcj(jdim,kdim,idim),
     .          arck(jdim,kdim,idim)
c 
      dimension phi(2),psi(2),omg(2)
      common /zero/ iexp
c
c     tolerance for switch to linear blending function
c     (10.**(-iexp) is machine zero)
c
      tol = max(1.e-07,10.**(-iexp+1))
c
c
      if (i1.eq.i2) then
c
         i = i1
c
c        check that j and k ranges span a logically 2d face
c
         if (j1.eq.j2 .or. k1.eq.k2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping in tfiface...error'',
     .      '' in j,k range for face i ='',i4)') i1
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
c        eta projector  -  1st step 
c
         do j=j1,j2
            do k=k1,k2
               denom = (arcj(j2,k,i) - arcj(j1,k,i))
               if(real(denom).lt.real(tol)) then 
                 eta = 0.
               else
                 eta = (arcj(j,k,i)  - arcj(j1,k,i))
     .                  / denom
               end if
               psi(1) = eta
               psi(2) = 1.-eta
               x1(j,k,i) = psi(2)*x(j1,k,i) 
     .                   + psi(1)*x(j2,k,i)
               y1(j,k,i) = psi(2)*y(j1,k,i) 
     .                   + psi(1)*y(j2,k,i)
               z1(j,k,i) = psi(2)*z(j1,k,i) 
     .                   + psi(1)*z(j2,k,i)
            end do
         end do
c 
c        zeta projector  -  2nd step
c
         do k=k1,k2
            do j=j1,j2
               denom = (arck(j,k2,i) - arck(j,k1,i))
               if(real(denom).lt.real(tol)) then 
                 zeta = 0.
               else
                 zeta = (arck(j,k,i)  - arck(j,k1,i))
     .                  / denom
               end if
               omg(1) = zeta
               omg(2) = 1.-zeta
               x2(j,k,i) = omg(2)*(x(j,k1,i) - x1(j,k1,i))
     .                   + omg(1)*(x(j,k2,i) - x1(j,k2,i))
               y2(j,k,i) = omg(2)*(y(j,k1,i) - y1(j,k1,i))
     .                   + omg(1)*(y(j,k2,i) - y1(j,k2,i))
               z2(j,k,i) = omg(2)*(z(j,k1,i) - z1(j,k1,i))
     .                   + omg(1)*(z(j,k2,i) - z1(j,k2,i))
            end do
         end do
c
c        summation of eta and zeta projectors  -  3rd step
c
         do j=j1,j2
            do k=k1,k2
               x(j,k,i)  = x1(j,k,i) + x2(j,k,i)
               y(j,k,i)  = y1(j,k,i) + y2(j,k,i)
               z(j,k,i)  = z1(j,k,i) + z2(j,k,i)
            end do
         end do
c
      end if
c
      if (j1.eq.j2) then
c
         j = j1
c
c        check that i and k ranges span a logically 2d face
c
         if (i1.eq.i2 .or. k1.eq.k2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping in tfiface...error'',
     .      '' in i,k range for face j ='',i4)') j1
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
c        xi projector  -  1st step
c
         do i=i1,i2
            do k=k1,k2
               denom = (arci(j,k,i2) - arci(j,k,i1))
               if(real(denom).lt.real(tol)) then 
                 xi = 0.
               else
                 xi = (arci(j,k,i)  - arci(j,k,i1))
     .                / denom
               end if
               phi(1) = xi
               phi(2) = 1.-xi
               x1(j,k,i) = phi(2)*x(j,k,i1)
     .                   + phi(1)*x(j,k,i2)
               y1(j,k,i) = phi(2)*y(j,k,i1)
     .                   + phi(1)*y(j,k,i2)
               z1(j,k,i) = phi(2)*z(j,k,i1)
     .                   + phi(1)*z(j,k,i2)
            end do
         end do
c
c        zeta projector  -  2nd step
c
         do k=k1,k2
            do i=i1,i2
               denom = (arck(j,k2,i) - arck(j,k1,i))
               if(real(denom).lt.real(tol)) then 
                 zeta = 0.
               else
                 zeta = (arck(j,k,i)  - arck(j,k1,i))
     .                  / denom
               end if
               omg(1) = zeta
               omg(2) = 1.-zeta
               x2(j,k,i) = omg(2)*(x(j,k1,i) - x1(j,k1,i))
     .                   + omg(1)*(x(j,k2,i) - x1(j,k2,i))
               y2(j,k,i) = omg(2)*(y(j,k1,i) - y1(j,k1,i))
     .                   + omg(1)*(y(j,k2,i) - y1(j,k2,i))
               z2(j,k,i) = omg(2)*(z(j,k1,i) - z1(j,k1,i))
     .                   + omg(1)*(z(j,k2,i) - z1(j,k2,i))
            end do
         end do
c
c        summation of xi and zeta projectors  -  3rd step
c
         do i=i1,i2
               do k=k1,k2
               x(j,k,i)  = x1(j,k,i) + x2(j,k,i)
               y(j,k,i)  = y1(j,k,i) + y2(j,k,i)
               z(j,k,i)  = z1(j,k,i) + z2(j,k,i)
            end do
         end do
c
      end if
c
      if (k1.eq.k2) then
c
         k = k1
c
c        check that i and j ranges span a logically 2d face
c
         if (i1.eq.i2 .or. j1.eq.j2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping in tfiface...error'',
     .      '' in i,j range for face k ='',i4)') k1
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
c        xi projector  -  1st step
c
         do i=i1,i2
            do j=j1,j2
               denom = (arci(j,k,i2) - arci(j,k,i1))
               if(real(denom).lt.real(tol)) then 
                 xi = 0.
               else
                 xi = (arci(j,k,i)  - arci(j,k,i1))
     .                / denom
               end if
               phi(1) = xi
               phi(2) = 1.-xi
               x1(j,k,i) = phi(2)*x(j,k,i1)
     .                   + phi(1)*x(j,k,i2)
               y1(j,k,i) = phi(2)*y(j,k,i1)
     .                   + phi(1)*y(j,k,i2)
               z1(j,k,i) = phi(2)*z(j,k,i1)
     .                   + phi(1)*z(j,k,i2)
            end do
         end do
c
c       eta projector  -  2nd step
c
         do j=j1,j2
            do i=i1,i2
               denom = (arcj(j2,k,i) - arcj(j1,k,i))
               if(real(denom).lt.real(tol)) then 
                 eta = 0.
               else
                 eta = (arcj(j,k,i)  - arcj(j1,k,i))
     .               / denom
               end if
               psi(1) = eta
               psi(2) = 1.-eta
               x2(j,k,i) = psi(2)*(x(j1,k,i) - x1(j1,k,i))
     .                   + psi(1)*(x(j2,k,i) - x1(j2,k,i))
               y2(j,k,i) = psi(2)*(y(j1,k,i) - y1(j1,k,i))
     .                   + psi(1)*(y(j2,k,i) - y1(j2,k,i))
               z2(j,k,i) = psi(2)*(z(j1,k,i) - z1(j1,k,i))
     .                   + psi(1)*(z(j2,k,i) - z1(j2,k,i))
            end do
         end do
c
c        summation of xi and eta projectors  -  3rd step
c
         do i=i1,i2
            do j=j1,j2
               x(j,k,i)  = x1(j,k,i) + x2(j,k,i)
               y(j,k,i)  = y1(j,k,i) + y2(j,k,i)
               z(j,k,i)  = z1(j,k,i) + z2(j,k,i)
            end do
         end do
c
      end if
c
      return
      end 
