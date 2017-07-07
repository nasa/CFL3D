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
      subroutine tfivol(idim,jdim,kdim,iskp,jskp,kskp,iskmax,jskmax,
     .                  kskmax,iskip,jskip,kskip,isktyp,x,y,z,
     .                  x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,
     .                  z5,x6,y6,z6,x7,y7,z7,arci,arcj,arck,nou,bou,
     .                  nbuf,ibufdim,myid,maxbl,nbl)
c
c     $Id$
c
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
      dimension x4(jdim,kdim,idim),y4(jdim,kdim,idim),z4(jdim,kdim,idim)
      dimension x5(jdim,kdim,idim),y5(jdim,kdim,idim),z5(jdim,kdim,idim)
      dimension x6(jdim,kdim,idim),y6(jdim,kdim,idim),z6(jdim,kdim,idim)
      dimension x7(jdim,kdim,idim),y7(jdim,kdim,idim),z7(jdim,kdim,idim)
      dimension arci(jdim,kdim,idim),arcj(jdim,kdim,idim),
     .          arck(jdim,kdim,idim)
      dimension iskip(maxbl,500)
      dimension jskip(maxbl,500)
      dimension kskip(maxbl,500)
      dimension iskmax(maxbl),jskmax(maxbl),kskmax(maxbl)
      dimension phi(2),psi(2),omg(2)
c
      common /zero/ iexp
c
c     tolerance for switch to linear blending function
c     (10.**(-iexp) is machine zero)
c
      tol = max(1.e-07,10.**(-iexp+1))
c
      if(abs(isktyp).eq.1) then
       do i1= 1,idim-iskp,iskp
        is = i1
        ie = i1+iskp
        do j1= 1,jdim-jskp,jskp
         js = j1
         je = j1+jskp
         do k1= 1,kdim-kskp,kskp
          ks = k1
          ke = k1+kskp
          do j=js,je
           do k=ks,ke
            do i=is,ie
               denomj= (arcj(je,k,i)- arcj(js,k,i))
               if(real(denomj).lt.real(tol)) then
                 eta = 0.
               else
                 eta = (arcj(j,k,i) - arcj(js,k,i))
     .                  / denomj
               end if
               psi(1) = eta
               psi(2) = 1.-eta
               denomk= (arck(j,ke,i)- arck(j,ks,i))
               if(real(denomk).lt.real(tol)) then
                 zeta = 0.
               else
                 zeta = (arck(j,k,i) - arck(j,ks,i))
     .                  / denomk
               end if
               omg(1) = zeta
               omg(2) = 1.-zeta
               denomi = (arci(j,k,ie)- arci(j,k,is))
               if(real(denomi).lt.real(tol)) then
                 xi = 0.
               else
                 xi = (arci(j,k,i) - arci(j,k,is))
     .                  / denomi
               end if
               phi(1) = xi
               phi(2) = 1.-xi
               x1(j,k,i) = psi(2)*x(js,k,i) + psi(1)*x(je,k,i)
               y1(j,k,i) = psi(2)*y(js,k,i) + psi(1)*y(je,k,i)
               z1(j,k,i) = psi(2)*z(js,k,i) + psi(1)*z(je,k,i)
               x2(j,k,i) = omg(2)*x(j,ks,i) + omg(1)*x(j,ke,i)
               y2(j,k,i) = omg(2)*y(j,ks,i) + omg(1)*y(j,ke,i)
               z2(j,k,i) = omg(2)*z(j,ks,i) + omg(1)*z(j,ke,i)
               x3(j,k,i) = phi(2)*x(j,k,is) + phi(1)*x(j,k,ie)
               y3(j,k,i) = phi(2)*y(j,k,is) + phi(1)*y(j,k,ie)
               z3(j,k,i) = phi(2)*z(j,k,is) + phi(1)*z(j,k,ie)
               x4(j,k,i) = psi(2)*omg(2)     *x(js,ks,i)
     .                    +psi(2)*omg(1)     *x(js,ke,i)
     .                    +psi(1)*omg(2)     *x(je,ks,i)
     .                    +psi(1)*omg(1)     *x(je,ke,i)
               y4(j,k,i) = psi(2)*omg(2)     *y(js,ks,i)
     .                    +psi(2)*omg(1)     *y(js,ke,i)
     .                    +psi(1)*omg(2)     *y(je,ks,i)
     .                    +psi(1)*omg(1)     *y(je,ke,i)
               z4(j,k,i) = psi(2)*omg(2)     *z(js,ks,i)
     .                    +psi(2)*omg(1)     *z(js,ke,i)
     .                    +psi(1)*omg(2)     *z(je,ks,i)
     .                    +psi(1)*omg(1)     *z(je,ke,i)
               x5(j,k,i) = psi(2)*phi(2)     *x(js,k,is)
     .                    +psi(2)*phi(1)     *x(js,k,ie)
     .                    +psi(1)*phi(2)     *x(je,k,is)
     .                    +psi(1)*phi(1)     *x(je,k,ie)
               y5(j,k,i) = psi(2)*phi(2)     *y(js,k,is)
     .                    +psi(2)*phi(1)     *y(js,k,ie)
     .                    +psi(1)*phi(2)     *y(je,k,is)
     .                    +psi(1)*phi(1)     *y(je,k,ie)
               z5(j,k,i) = psi(2)*phi(2)     *z(js,k,is)
     .                    +psi(2)*phi(1)     *z(js,k,ie)
     .                    +psi(1)*phi(2)     *z(je,k,is)
     .                    +psi(1)*phi(1)     *z(je,k,ie)
               x6(j,k,i) = omg(2)*phi(2)     *x(j,ks,is)
     .                    +omg(2)*phi(1)     *x(j,ks,ie)
     .                    +omg(1)*phi(2)     *x(j,ke,is)
     .                    +omg(1)*phi(1)     *x(j,ke,ie)
               y6(j,k,i) = omg(2)*phi(2)     *y(j,ks,is)
     .                    +omg(2)*phi(1)     *y(j,ks,ie)
     .                    +omg(1)*phi(2)     *y(j,ke,is)
     .                    +omg(1)*phi(1)     *y(j,ke,ie)
               z6(j,k,i) = omg(2)*phi(2)     *z(j,ks,is)
     .                    +omg(2)*phi(1)     *z(j,ks,ie)
     .                    +omg(1)*phi(2)     *z(j,ke,is)
     .                    +omg(1)*phi(1)     *z(j,ke,ie)
               x7(j,k,i) = psi(2)*omg(2)*phi(2)*x(js,ks,is)
     .                    +psi(2)*omg(2)*phi(1)*x(js,ks,ie)
     .                    +psi(2)*omg(1)*phi(2)*x(js,ke,is)
     .                    +psi(2)*omg(1)*phi(1)*x(js,ke,ie)
     .                    +psi(1)*omg(2)*phi(2)*x(je,ks,is)
     .                    +psi(1)*omg(2)*phi(1)*x(je,ks,ie)
     .                    +psi(1)*omg(1)*phi(2)*x(je,ke,is)
     .                    +psi(1)*omg(1)*phi(1)*x(je,ke,ie)
               y7(j,k,i) = psi(2)*omg(2)*phi(2)*y(js,ks,is)
     .                    +psi(2)*omg(2)*phi(1)*y(js,ks,ie)
     .                    +psi(2)*omg(1)*phi(2)*y(js,ke,is)
     .                    +psi(2)*omg(1)*phi(1)*y(js,ke,ie)
     .                    +psi(1)*omg(2)*phi(2)*y(je,ks,is)
     .                    +psi(1)*omg(2)*phi(1)*y(je,ks,ie)
     .                    +psi(1)*omg(1)*phi(2)*y(je,ke,is)
     .                    +psi(1)*omg(1)*phi(1)*y(je,ke,ie)
               z7(j,k,i) = psi(2)*omg(2)*phi(2)*z(js,ks,is)
     .                    +psi(2)*omg(2)*phi(1)*z(js,ks,ie)
     .                    +psi(2)*omg(1)*phi(2)*z(js,ke,is)
     .                    +psi(2)*omg(1)*phi(1)*z(js,ke,ie)
     .                    +psi(1)*omg(2)*phi(2)*z(je,ks,is)
     .                    +psi(1)*omg(2)*phi(1)*z(je,ks,ie)
     .                    +psi(1)*omg(1)*phi(2)*z(je,ke,is)
     .                    +psi(1)*omg(1)*phi(1)*z(je,ke,ie)
            end do
           end do
          end do
          do i=is+1,ie-1
           do j=js+1,je-1
            do k=ks+1,ke-1
               x(j,k,i)  = x1(j,k,i) + x2(j,k,i) + x3(j,k,i)
     .                   - x4(j,k,i) - x5(j,k,i) - x6(j,k,i)
     .                   + x7(j,k,i)
               y(j,k,i)  = y1(j,k,i) + y2(j,k,i) + y3(j,k,i)
     .                   - y4(j,k,i) - y5(j,k,i) - y6(j,k,i)
     .                   + y7(j,k,i)
               z(j,k,i)  = z1(j,k,i) + z2(j,k,i) + z3(j,k,i)
     .                   - z4(j,k,i) - z5(j,k,i) - z6(j,k,i)
     .                   + z7(j,k,i)
            end do
           end do
          end do
         enddo
        enddo
       enddo
      else
       do i1=1,iskmax(nbl)-1
        is= iskip(nbl,i1)
        ie= iskip(nbl,i1+1)
        do j1=1,jskmax(nbl)-1
         js= jskip(nbl,j1)
         je= jskip(nbl,j1+1)
         do k1=1,kskmax(nbl)-1
          ks= kskip(nbl,k1)
          ke= kskip(nbl,k1+1)
c
          do j=js,je
           do k=ks,ke
            do i=is,ie
               denomj= (arcj(je,k,i)- arcj(js,k,i))
               if(real(denomj).lt.real(tol)) then
                 eta = 0.
               else
                 eta = (arcj(j,k,i) - arcj(js,k,i))
     .                  / denomj
               end if
               psi(1) = eta
               psi(2) = 1.-eta
               denomk= (arck(j,ke,i)- arck(j,ks,i))
               if(real(denomk).lt.real(tol)) then
                 zeta = 0.
               else
                 zeta = (arck(j,k,i) - arck(j,ks,i))
     .                  / denomk
               end if
               omg(1) = zeta
               omg(2) = 1.-zeta
               denomi = (arci(j,k,ie)- arci(j,k,is))
               if(real(denomi).lt.real(tol)) then
                 xi = 0.
               else
                 xi = (arci(j,k,i) - arci(j,k,is))
     .                  / denomi
               end if
               phi(1) = xi
               phi(2) = 1.-xi
               x1(j,k,i) = psi(2)*x(js,k,i) + psi(1)*x(je,k,i)
               y1(j,k,i) = psi(2)*y(js,k,i) + psi(1)*y(je,k,i)
               z1(j,k,i) = psi(2)*z(js,k,i) + psi(1)*z(je,k,i)
               x2(j,k,i) = omg(2)*x(j,ks,i) + omg(1)*x(j,ke,i)
               y2(j,k,i) = omg(2)*y(j,ks,i) + omg(1)*y(j,ke,i)
               z2(j,k,i) = omg(2)*z(j,ks,i) + omg(1)*z(j,ke,i)
               x3(j,k,i) = phi(2)*x(j,k,is) + phi(1)*x(j,k,ie)
               y3(j,k,i) = phi(2)*y(j,k,is) + phi(1)*y(j,k,ie)
               z3(j,k,i) = phi(2)*z(j,k,is) + phi(1)*z(j,k,ie)
               x4(j,k,i) = psi(2)*omg(2)     *x(js,ks,i)
     .                    +psi(2)*omg(1)     *x(js,ke,i)
     .                    +psi(1)*omg(2)     *x(je,ks,i)
     .                    +psi(1)*omg(1)     *x(je,ke,i)
               y4(j,k,i) = psi(2)*omg(2)     *y(js,ks,i)
     .                    +psi(2)*omg(1)     *y(js,ke,i)
     .                    +psi(1)*omg(2)     *y(je,ks,i)
     .                    +psi(1)*omg(1)     *y(je,ke,i)
               z4(j,k,i) = psi(2)*omg(2)     *z(js,ks,i)
     .                    +psi(2)*omg(1)     *z(js,ke,i)
     .                    +psi(1)*omg(2)     *z(je,ks,i)
     .                    +psi(1)*omg(1)     *z(je,ke,i)
               x5(j,k,i) = psi(2)*phi(2)     *x(js,k,is)
     .                    +psi(2)*phi(1)     *x(js,k,ie)
     .                    +psi(1)*phi(2)     *x(je,k,is)
     .                    +psi(1)*phi(1)     *x(je,k,ie)
               y5(j,k,i) = psi(2)*phi(2)     *y(js,k,is)
     .                    +psi(2)*phi(1)     *y(js,k,ie)
     .                    +psi(1)*phi(2)     *y(je,k,is)
     .                    +psi(1)*phi(1)     *y(je,k,ie)
               z5(j,k,i) = psi(2)*phi(2)     *z(js,k,is)
     .                    +psi(2)*phi(1)     *z(js,k,ie)
     .                    +psi(1)*phi(2)     *z(je,k,is)
     .                    +psi(1)*phi(1)     *z(je,k,ie)
               x6(j,k,i) = omg(2)*phi(2)     *x(j,ks,is)
     .                    +omg(2)*phi(1)     *x(j,ks,ie)
     .                    +omg(1)*phi(2)     *x(j,ke,is)
     .                    +omg(1)*phi(1)     *x(j,ke,ie)
               y6(j,k,i) = omg(2)*phi(2)     *y(j,ks,is)
     .                    +omg(2)*phi(1)     *y(j,ks,ie)
     .                    +omg(1)*phi(2)     *y(j,ke,is)
     .                    +omg(1)*phi(1)     *y(j,ke,ie)
               z6(j,k,i) = omg(2)*phi(2)     *z(j,ks,is)
     .                    +omg(2)*phi(1)     *z(j,ks,ie)
     .                    +omg(1)*phi(2)     *z(j,ke,is)
     .                    +omg(1)*phi(1)     *z(j,ke,ie)
               x7(j,k,i) = psi(2)*omg(2)*phi(2)*x(js,ks,is)
     .                    +psi(2)*omg(2)*phi(1)*x(js,ks,ie)
     .                    +psi(2)*omg(1)*phi(2)*x(js,ke,is)
     .                    +psi(2)*omg(1)*phi(1)*x(js,ke,ie)
     .                    +psi(1)*omg(2)*phi(2)*x(je,ks,is)
     .                    +psi(1)*omg(2)*phi(1)*x(je,ks,ie)
     .                    +psi(1)*omg(1)*phi(2)*x(je,ke,is)
     .                    +psi(1)*omg(1)*phi(1)*x(je,ke,ie)
               y7(j,k,i) = psi(2)*omg(2)*phi(2)*y(js,ks,is)
     .                    +psi(2)*omg(2)*phi(1)*y(js,ks,ie)
     .                    +psi(2)*omg(1)*phi(2)*y(js,ke,is)
     .                    +psi(2)*omg(1)*phi(1)*y(js,ke,ie)
     .                    +psi(1)*omg(2)*phi(2)*y(je,ks,is)
     .                    +psi(1)*omg(2)*phi(1)*y(je,ks,ie)
     .                    +psi(1)*omg(1)*phi(2)*y(je,ke,is)
     .                    +psi(1)*omg(1)*phi(1)*y(je,ke,ie)
               z7(j,k,i) = psi(2)*omg(2)*phi(2)*z(js,ks,is)
     .                    +psi(2)*omg(2)*phi(1)*z(js,ks,ie)
     .                    +psi(2)*omg(1)*phi(2)*z(js,ke,is)
     .                    +psi(2)*omg(1)*phi(1)*z(js,ke,ie)
     .                    +psi(1)*omg(2)*phi(2)*z(je,ks,is)
     .                    +psi(1)*omg(2)*phi(1)*z(je,ks,ie)
     .                    +psi(1)*omg(1)*phi(2)*z(je,ke,is)
     .                    +psi(1)*omg(1)*phi(1)*z(je,ke,ie)
            end do
           end do
          end do
c
          do i=is+1,ie-1
           do j=js+1,je-1
            do k=ks+1,ke-1
               x(j,k,i)  = x1(j,k,i) + x2(j,k,i) + x3(j,k,i)
     .                   - x4(j,k,i) - x5(j,k,i) - x6(j,k,i)
     .                   + x7(j,k,i)
               y(j,k,i)  = y1(j,k,i) + y2(j,k,i) + y3(j,k,i)
     .                   - y4(j,k,i) - y5(j,k,i) - y6(j,k,i)
     .                   + y7(j,k,i)
               z(j,k,i)  = z1(j,k,i) + z2(j,k,i) + z3(j,k,i)
     .                   - z4(j,k,i) - z5(j,k,i) - z6(j,k,i)
     .                   + z7(j,k,i)
            end do
           end do
          end do
         enddo
        enddo
       enddo
 
      end if

      return
      end 
