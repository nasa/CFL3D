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
      subroutine genforce(jdim,kdim,idim,sk,sj,si,cz,cy,cx,bcj,bck,bci,
     .                    blank,nbl,qj0,qk0,qi0,maxbl,maxseg,n2,xmdj,
     .                    xmdk,xmdi,aesrfdat,nmds,maxaes,maxsegdg,
     .                    nsegdfrm,jcsi,jcsf,kcsi,kcsf,icsi,icsf,
     .                    idfrmseg,iaes,iaesurf)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Integrate the generalized forces on the body for 
c               aeroelastic analysis. Only pressure forces are 
c               included. The computations of the generalized 
c               forces are consistant with current methodolgy used
c               calculating "regular" forces in subroutine force. 
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xmdj(kdim,idim,6,nmds),xmdk(jdim,idim,6,nmds),
     .          xmdi(jdim,kdim,6,nmds)
      dimension aesrfdat(5,maxaes)
      dimension sk(jdim,kdim,idim-1,5),sj(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5)
      dimension blank(jdim,kdim,idim)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),
     .          qi0(jdim,kdim,5,4)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg),
     .          iaesurf(maxbl,maxsegdg)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /igrdtyp/ ip3dgrd,ialph
c
c**************************************************************************
c
c     sign conventions for generalized forces:
c
c       cx.....x-component of force
c              positive for force in +x direction
c       cy.....y-component of force
c              positive for force in +y direction
c       cz.....z-component of force
c              positive for force in +z direction
c
c     xmdj,xmdk,xmdi are the modal deflections for the
c     j,k,i faces of the current block. The 4th index of
c     these arrays indicates both the face number and the
c     x,y,z component, e.g.:
c
c       xmdj(k,i,ll,m), ll = 1...x-component of mode m on j=1
c                       ll = 2...y-component of mode m on j=1
c                       ll = 3...z-component of mode m on j=1
c                       ll = 4...x-component of mode m on j=jdim
c                       ll = 5...y-component of mode m on j=jdim
c                       ll = 6...z-component of mode m on j=jdim
c
c**************************************************************************
c
      cpc   = 2.e0/(gamma*xmach*xmach)
c
      cz    = 0.e0
      cy    = 0.e0 
      cx    = 0.e0
c
c     loop over all aeroelastic surface segments (those with idfrmseg=99)
c
      do iseg=1,nsegdfrm(nbl)
c
         if (iaesurf(nbl,iseg).eq.iaes .and.
     .       idfrmseg(nbl,iseg).eq.99) then
c
c******************************************************************
c           generalized forces on k=constant surfaces
c******************************************************************
c
            if (kcsi(nbl,iseg).eq.kcsf(nbl,iseg)) then
c
               k = kcsi(nbl,iseg)
               if (k.eq.1) then
                  kk   = 1
                  kc   = 1
                  l    = 0
               else if (k.eq.kdim) then
                  kk   = 2
                  kc   = kdim - 1
                  l    = 3
               end if
c
               ist = icsi(nbl,iseg)
               ifn = icsf(nbl,iseg) - 1
               jst = jcsi(nbl,iseg)
               jfn = jcsf(nbl,iseg) - 1
c
               sgn = 1.0
               if(kk .gt. 1) sgn = -1.0
c
               do 10 i=ist,ifn 
               cxl   = 0.e0
               cyl   = 0.e0
               czl   = 0.e0
c
               do 15 j=jst,jfn
c
               xmdavg = (xmdk(j,i  ,l+1,n2)+xmdk(j+1,i  ,l+1,n2)
     .                +  xmdk(j,i+1,l+1,n2)+xmdk(j+1,i+1,l+1,n2))*.25
               ymdavg = (xmdk(j,i  ,l+2,n2)+xmdk(j+1,i  ,l+2,n2)
     .                +  xmdk(j,i+1,l+2,n2)+xmdk(j+1,i+1,l+2,n2))*.25
               zmdavg = (xmdk(j,i  ,l+3,n2)+xmdk(j+1,i  ,l+3,n2)
     .                +  xmdk(j,i+1,l+3,n2)+xmdk(j+1,i+1,l+3,n2))*.25
c          
               dcp = -(qk0(j,i,5,kk+kk-1)/p0-1.e0)*cpc*sk(j,k,i,4)
               dcx = dcp*sk(j,k,i,1)*sgn 
               dcy = dcp*sk(j,k,i,2)*sgn 
               dcz = dcp*sk(j,k,i,3)*sgn 
c
c              only use contributions from points with interface 
c              (solid surface) boundary conditions, and only those
c              points not blanked out
c
               fact = bck(j,i,kk)*blank(j,kc,i)
               dcx  = dcx*fact
               dcy  = dcy*fact
               dcz  = dcz*fact
c
               cxl   = cxl+dcx*xmdavg 
               cyl   = cyl+dcy*ymdavg
               czl   = czl+dcz*zmdavg
c
   15          continue
c
c              integrated values
               cz    = cz+czl
               cy    = cy+cyl
               cx    = cx+cxl 
c
   10          continue
c
            end if
c
c******************************************************************
c           generalized forces on j=constant surfaces
c******************************************************************
c
            if (jcsi(nbl,iseg).eq.jcsf(nbl,iseg)) then
c
               j = jcsi(nbl,iseg)
               if (j.eq.1) then
                  jj   = 1
                  jc   = 1
                  l    = 0
               else if (j.eq.jdim) then
                  jj   = 2
                  jc   = jdim - 1
                  l    = 3
               end if
c
               ist = icsi(nbl,iseg)
               ifn = icsf(nbl,iseg) - 1
               kst = kcsi(nbl,iseg)
               kfn = kcsf(nbl,iseg) - 1
c
               sgn = 1.0
               if(jj .gt. 1) sgn = -1.0
c
               do 40 i=ist,ifn
               cxl   = 0.e0
               cyl   = 0.e0
               czl   = 0.e0
c
               do 45 k=kst,kfn
c
               xmdavg = (xmdj(k,i  ,l+1,n2)+xmdj(k+1,i  ,l+1,n2)
     .                +  xmdj(k,i+1,l+1,n2)+xmdj(k+1,i+1,l+1,n2))*.25
               ymdavg = (xmdj(k,i  ,l+2,n2)+xmdj(k+1,i  ,l+2,n2)
     .                +  xmdj(k,i+1,l+2,n2)+xmdj(k+1,i+1,l+2,n2))*.25
               zmdavg = (xmdj(k,i  ,l+3,n2)+xmdj(k+1,i  ,l+3,n2)
     .                +  xmdj(k,i+1,l+3,n2)+xmdj(k+1,i+1,l+3,n2))*.25
c
               dcp = -(qj0(k,i,5,jj+jj-1)/p0-1.e0)*cpc*sj(j,k,i,4)
               dcx = dcp*sj(j,k,i,1)*sgn
               dcy = dcp*sj(j,k,i,2)*sgn
               dcz = dcp*sj(j,k,i,3)*sgn
c
c              only use contributions from points with interface 
c              (solid surface) boundary conditions, and only those
c              points not blanked out
c
               fact = bcj(k,i,jj)*blank(jc,k,i)
               dcx  = dcx*fact
               dcy  = dcy*fact
               dcz  = dcz*fact
c
               cxl   = cxl+dcx*xmdavg
               cyl   = cyl+dcy*ymdavg
               czl   = czl+dcz*zmdavg
c
   45          continue
c
c              integrated values
               cz    = cz+czl
               cy    = cy+cyl
               cx    = cx+cxl
c
   40          continue
c
            end if
c
c******************************************************************
c           generalized forces on i=constant surfaces
c******************************************************************
c
            if (icsi(nbl,iseg).eq.icsf(nbl,iseg)) then
c
               i = icsi(nbl,iseg)
               if (i.eq.1) then
                  ii   = 1
                  ic   = 1
                  l    = 0
               else if (i.eq.idim) then
                  ii   = 2
                  ic   = idim - 1
                  l    = 3
               end if
c
               jst = jcsi(nbl,iseg)
               jfn = jcsf(nbl,iseg) - 1
               kst = kcsi(nbl,iseg)
               kfn = kcsf(nbl,iseg) - 1
c
               sgn = 1.0
               if(ii .gt. 1) sgn = -1.0
c
               do 70 j=jst,jfn
               cxl   = 0.e0
               cyl   = 0.e0
               czl   = 0.e0
c
               do 75 k=kst,kfn
c
               xmdavg = (xmdi(j,k  ,l+1,n2)+xmdi(j+1,k  ,l+1,n2)
     .                +  xmdi(j,k+1,l+1,n2)+xmdi(j+1,k+1,l+1,n2))*.25
               ymdavg = (xmdi(j,k  ,l+2,n2)+xmdi(j+1,k  ,l+2,n2)
     .                +  xmdi(j,k+1,l+2,n2)+xmdi(j+1,k+1,l+2,n2))*.25
               zmdavg = (xmdi(j,k  ,l+3,n2)+xmdi(j+1,k  ,l+3,n2)
     .                +  xmdi(j,k+1,l+3,n2)+xmdi(j+1,k+1,l+3,n2))*.25
c
               dcp = -(qi0(j,k,5,ii+ii-1)/p0-1.e0)*cpc*si(j,k,i,4)
               dcx = dcp*si(j,k,i,1)*sgn
               dcy = dcp*si(j,k,i,2)*sgn
               dcz = dcp*si(j,k,i,3)*sgn
c
c              only use contributions from points with interface 
c              (solid surface) boundary conditions, and only those
c              points not blanked out
c
               fact = bci(j,k,ii)*blank(j,k,ic)
               dcx  = dcx*fact
               dcy  = dcy*fact
               dcz  = dcz*fact
c
               cxl   = cxl+dcx*xmdavg
               cyl   = cyl+dcy*ymdavg
               czl   = czl+dcz*zmdavg
c
   75          continue
c
c              integrated values
               cz    = cz+czl
               cy    = cy+cyl
               cx    = cx+cxl
c
   70          continue
c
            end if
c
         end if
c
      end do
c
      return
      end 
