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
      subroutine update(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,vol,dtj,
     .                  vist3d,blank,x,y,z,res,wk0,vmuk,vmuj,vmui,wk,
     .                  nwork,nbl,iover,vk0,bcj,bck,bci,nou,bou,nbuf,
     .                  ibufdim,myid,mblk2nd,maxbl,volk0,tursav,
     .                  tk0,cmuv,iadvance,nummem,ux)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Update the solution in time.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension wk(nwork),wk0(idim*jdim*22)
      dimension q(jdim,kdim,idim,5),res(jdim,kdim,idim-1,5)
      dimension vmuk(jdim-1,idim-1,2),vmuj(kdim-1,idim-1,2),
     .          vmui(jdim-1,kdim-1,2)
      dimension x(jdim*kdim*idim),y(jdim*kdim*idim),z(jdim*kdim*idim)
      dimension qj0(kdim*(idim-1),5),qk0(jdim*(idim-1),5),
     .          qi0(jdim*kdim,5)
      dimension sk(jdim*kdim*(idim-1),5),si(jdim*kdim*idim,5),
     .          sj(jdim*kdim*(idim-1),5)
      dimension vol(jdim*kdim*(idim-1)),dtj(jdim*kdim*(idim-1)),
     .          vist3d(jdim,kdim,idim),blank(jdim,kdim,idim)
     .         ,vk0(jdim,idim-1,1,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension mblk2nd(maxbl),volk0(jdim,idim-1,4)
      dimension iadvance(maxbl)
      dimension tursav(jdim,kdim,idim,nummem),tk0(jdim,idim-1,nummem,4),
     .          cmuv(idim-1,kdim-1,idim-1)
      dimension ux(jdim-1,kdim-1,idim-1,9)
c
      common /chk/ ichk
      common /sklton/ isklton
      common /twod/ i2d
      common /mms/ iexact_trunc,iexact_disc,iexact_ring
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
      common /iupdate/ iupdatemean
c
      if (iupdatemean .eq. 0) return
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      term  = 1.0
      if (i2d.eq.1 .or. iforcev0.eq.1) term = 0.0
c
      alpq = -0.2
      phiq = 1./0.5
      betq = 1.0 + alpq*phiq
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),17) nbl
      end if
   17 format(1x,46hupdating with 3-D 3-factor AF scheme for block,i6)
c
c     assemble and solve matrix equations
c
      call af3f(nbl,jdim,kdim,idim,q,vol,qj0,qk0,qi0,dtj,sj,sk,si,
     .          res,vist3d,x,y,z,blank,vmuk,resd,wk,nwork,wk0,1,iover,
     .          vk0,bcj,bck,bci,nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl,
     .          volk0,tursav,tk0,cmuv,iadvance,nummem,ux)
c
c      update 3-d solution
c
      nt    = jdim*kdim
      nplq  = min(idim1,999000/nt)
      npl   = nplq
      do 4100 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      n = nt*npl - jdim -1
cdir$ ivdep 
      do 1002 izz=1,n 
c
c      ensure positivity of density
c
      t1             = res(izz,1,i,1)/q(izz,1,i,1)
      t2             = res(izz,1,i,1)/( betq + ccabs(t1)*phiq )
      res(izz,1,i,1) = ccvmgt(t2,res(izz,1,i,1),
     .                 (real(t1).lt.real(alpq)))
c
c      ensure positivity of pressure
c
      t1             = res(izz,1,i,5)/q(izz,1,i,5)
      t2             = res(izz,1,i,5)/( betq + ccabs(t1)*phiq )
      res(izz,1,i,5) = ccvmgt(t2,res(izz,1,i,5),
     .                 (real(t1).lt.real(alpq)))
c
      q(izz,1,i,1) =  q(izz,1,i,1)+res(izz,1,i,1)
      q(izz,1,i,2) =  q(izz,1,i,2)+res(izz,1,i,2)
      q(izz,1,i,3) = (q(izz,1,i,3)+res(izz,1,i,3))*term
      q(izz,1,i,4) =  q(izz,1,i,4)+res(izz,1,i,4)
      q(izz,1,i,5) =  q(izz,1,i,5)+res(izz,1,i,5)
c
 1002 continue
 4100 continue
      if (iexact_ring .eq. 1) then
c       overwrite ring of exact MMS values in outer 2 rows of grid
        call exact_flow_q_ring(jdim,kdim,idim,x,y,z,q,
     +      iexact_trunc,iexact_disc)
      end if
c
      if (ichk.eq.1) then
         epsz = 1.0e-05
         epss = 1.0e+03
         do 4200 i=1,idim1
         do 4200 k=1,kdim1
         do 4200 j=1,jdim1
         if (real(q(j,k,i,5)).lt.real(epsz) .or.
     .       real(q(j,k,i,1)).lt.real(epsz) .or.
     .       real(q(j,k,i,5)).gt.real(epss) .or.
     .       real(q(j,k,i,1)).gt.real(epss)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in update on block ',nbl
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1500)j,k,i,(real(q(j,k,i,m)),m=1,5)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
 4200    continue
      end if
 1500 format(1x,32h *neg. (or large) d/p*(j,k,i,q)=,3i5,5e12.5)
      return
      end
