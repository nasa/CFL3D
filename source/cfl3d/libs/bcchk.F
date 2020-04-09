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
      subroutine bcchk(idim,jdim,kdim,q,qi0,qj0,qk0,blank,ibcflg,nbl,
     .                 nou,bou,nbuf,ibufdim,myid,istop,igridg,maxbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determine if boundary conditions have been set.  Also,
c     fill endpoints for safety.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5),qi0(jdim,kdim,5,4)
      dimension qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension blank(jdim,kdim,idim),igridg(maxbl)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      if (ibcflg.eq.0) then
c
c     don't reset boundary values from initial values for checking if
c     adjacent interior point is blanked out...otherwise this routine (below)
c     will think there is a problem when in fact there is not
c
      do 10 i=1,idim1
      do 20 k=1,kdim
      qj0(k,i,1,1) = (1.0-blank(1,k,i))*qj0(k,i,1,1)
     .              - 1.0*blank(1,k,i)
      qj0(k,i,5,1) = (1.0-blank(1,k,i))*qj0(k,i,5,1)
     .              - 1.0*blank(1,k,i)
      qj0(k,i,1,3) = (1.0-blank(jdim1,k,i))*qj0(k,i,1,3)
     .              - 1.0*blank(jdim1,k,i)
      qj0(k,i,5,3) = (1.0-blank(jdim1,k,i))*qj0(k,i,5,3)
     .              - 1.0*blank(jdim1,k,i)
   20 continue
      do 30 j=1,jdim
      qk0(j,i,1,1) = (1.0-blank(j,1,i))*qk0(j,i,1,1) 
     .              - 1.0*blank(j,1,i)
      qk0(j,i,5,1) = (1.0-blank(j,1,i))*qk0(j,i,5,1) 
     .              - 1.0*blank(j,1,i)
      qk0(j,i,1,3) = (1.0-blank(j,kdim1,i))*qk0(j,i,1,3) 
     .              - 1.0*blank(j,kdim1,i)
      qk0(j,i,5,3) = (1.0-blank(j,kdim1,i))*qk0(j,i,5,3) 
     .              - 1.0*blank(j,kdim1,i)
   30 continue
   10 continue
c
      do 40 k=1,kdim
      do 40 j=1,jdim
      qi0(j,k,1,1) = (1.0-blank(j,k,1))*qi0(j,k,1,1)
     .              - 1.0*blank(j,k,1)
      qi0(j,k,5,1) = (1.0-blank(j,k,1))*qi0(j,k,5,1)
     .              - 1.0*blank(j,k,1)
      qi0(j,k,1,3) = (1.0-blank(j,k,idim1))*qi0(j,k,1,3)
     .              - 1.0*blank(j,k,idim1)
      qi0(j,k,5,3) = (1.0-blank(j,k,idim1))*qi0(j,k,5,3)
     .              - 1.0*blank(j,k,idim1)
   40 continue
c
      else if (ibcflg.eq.1) then
c
c     fill endpoints for safety with multi-plane vectorization
c
      do 300 l=1,5
      do 320 i=1,idim1
      q(jdim,kdim,i,l) = q(jdim1,kdim1,i,l)
c
      qj0(kdim,i,l,1)  = q(1,kdim1,i,l)
      qj0(kdim,i,l,2)  = q(1,kdim1,i,l)
      qj0(kdim,i,l,3)  = q(jdim1,kdim1,i,l)
      qj0(kdim,i,l,4)  = q(jdim1,kdim1,i,l)
c
      qk0(jdim,i,l,1)  = qk0(jdim1,i,l,1)
      qk0(jdim,i,l,2)  = qk0(jdim1,i,l,2)
      qk0(jdim,i,l,3)  = qk0(jdim1,i,l,3)
      qk0(jdim,i,l,4)  = qk0(jdim1,i,l,4)
  320 continue
cdir$ ivdep
      do 330 izz=1,jdim1
      q(izz,kdim,idim,l)  = q(izz,kdim1,idim1,l)
c
      qi0(izz,kdim,l,1)   = q(izz,kdim1,1,l)
      qi0(izz,kdim,l,2)   = q(izz,kdim1,1,l)
      qi0(izz,kdim,l,3)   = q(izz,kdim1,idim1,l)
      qi0(izz,kdim,l,4)   = q(izz,kdim1,idim1,l)
  330 continue
c
      q(jdim,kdim,idim,l) = q(jdim1,kdim1,idim1,l)
c
      qi0(jdim,kdim,l,1)  = q(jdim1,kdim1,1,l)
      qi0(jdim,kdim,l,2)  = q(jdim1,kdim1,1,l)
      qi0(jdim,kdim,l,3)  = q(jdim1,kdim1,idim1,l)
      qi0(jdim,kdim,l,4)  = q(jdim1,kdim1,idim1,l)
c
      do 350 k=1,kdim1
      qi0(jdim,k,l,1)  = qi0(jdim1,k,l,1)
      qi0(jdim,k,l,2)  = qi0(jdim1,k,l,2)
      qi0(jdim,k,l,3)  = qi0(jdim1,k,l,3)
      qi0(jdim,k,l,4)  = qi0(jdim1,k,l,4)
c
      q(jdim,k,idim,l) = q(jdim1,k,idim,l)
  350 continue
  300 continue
c
      else if (ibcflg.eq.2) then
c
c     check for negative (or large) densities and pressures
c
      epsz = 1.0e-10
      epss = 1.0e+05
      do 110 i=1,idim1
      do 120 k=1,kdim
c 
      ib = int(blank(1,k,i))
c
      if ((real(qj0(k,i,1,1)).le.real(epsz) .or.
     .     real(qj0(k,i,5,1)).le.real(epsz) .or.
     .     real(qj0(k,i,1,1)).ge.real(epss) .or. 
     .     real(qj0(k,i,5,1)).ge.real(epss)) .and. 
     .   ib.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),200)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),15) nbl,igridg(nbl),k,i,real(qj0(k,i,1,1)),
     .                           real(qj0(k,i,5,1))
         if (qj0(k,i,1,1).eq.-1.0 .and. qj0(k,i,5,1).eq.-1.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         else
            if (istop.gt.0) then
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
      end if
   15 format(6h block,i4,7h (grid ,i4,1h),25h  on j=1 boundary at k,i=,
     .2i5,15h qj0(1),qj0(5)=,2e12.5)
c
      ib = int(blank(jdim1,k,i))
c
      if ((real(qj0(k,i,1,3)).le.real(epsz) .or. 
     .     real(qj0(k,i,5,3)).le.real(epsz) .or.
     .     real(qj0(k,i,1,3)).ge.real(epss) .or. 
     .     real(qj0(k,i,5,3)).ge.real(epss)) .and.
     .   ib.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),200)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),25) nbl,igridg(nbl),k,i,real(qj0(k,i,1,3)),
     .                           real(qj0(k,i,5,3))
         if (qj0(k,i,1,3).eq.-1.0 .and. qj0(k,i,5,3).eq.-1.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         else
            if (istop.gt.0) then
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
      end if
   25 format(6h block,i4,7h (grid ,i4,1h),23h  on j=jdim boundary at,
     .5h k,i=,2i5,15h qj0(1),qj0(5)=,2e12.5)
  120 continue
c
      do 130 j=1,jdim
c 
      ib = int(blank(j,1,i))
c
      if ((real(qk0(j,i,1,1)).le.real(epsz) .or. 
     .     real(qk0(j,i,5,1)).le.real(epsz) .or.
     .     real(qk0(j,i,1,1)).ge.real(epss) .or. 
     .     real(qk0(j,i,5,1)).ge.real(epss)) .and.
     .   ib.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),200)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),35) nbl,igridg(nbl),j,i,real(qk0(j,i,1,1)),
     .                           real(qk0(j,i,5,1))
         if (qk0(j,i,1,1).eq.-1.0 .and. qk0(j,i,5,1).eq.-1.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         else
            if (istop.gt.0) then
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
      end if
   35 format(6h block,i4,7h (grid ,i4,1h),25h  on k=1 boundary at j,i=,
     .2i5,15h qk0(1),qk0(5)=,2e12.5)
c
      ib = int(blank(j,kdim1,i))
c
      if ((real(qk0(j,i,1,3)).le.real(epsz) .or. 
     .     real(qk0(j,i,5,3)).le.real(epsz) .or.
     .     real(qk0(j,i,1,3)).ge.real(epss) .or. 
     .     real(qk0(j,i,5,3)).ge.real(epss)) .and.
     .   ib.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),200)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),45) nbl,igridg(nbl),j,i,real(qk0(j,i,1,3)),
     .                           real(qk0(j,i,5,3))
         if (real(qk0(j,i,1,3)).eq.-1.0 .and. 
     .       real(qk0(j,i,5,3)).eq.-1.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         else
            if (istop.gt.0) then
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
      end if
   45 format(6h block,i4,7h (grid ,i4,1h),23h  on k=kdim boundary at,
     .5h j,i=,2i5,15h qk0(1),qk0(5)=,2e12.5)
  130 continue
  110 continue
c
      do 140 k=1,kdim
      do 140 j=1,jdim
c
      ib = int(blank(j,k,1))
c
      if ((real(qi0(j,k,1,1)).le.real(epsz) .or.
     .     real(qi0(j,k,5,1)).le.real(epsz) .or.
     .     real(qi0(j,k,1,1)).ge.real(epss) .or. 
     .     real(qi0(j,k,5,1)).ge.real(epss)) .and.
     .   ib.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),200)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),55) nbl,igridg(nbl),j,k,real(qi0(j,k,1,1)),
     .                           real(qi0(j,k,5,1))
         if (real(qi0(j,k,1,1)).eq.-1.0 .and. 
     .       real(qi0(j,k,5,1)).eq.-1.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         else
            if (istop.gt.0) then
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
      end if
   55 format(6h block,i4,7h (grid ,i4,1h),25h  on i=1 boundary at j,k=,
     .2i5,15h qi0(1),qi0(5)=,2e12.5)
c
      ib = int(blank(j,k,idim1))
c
      if ((real(qi0(j,k,1,3)).le.real(epsz) .or. 
     .     real(qi0(j,k,5,3)).le.real(epsz) .or.
     .     real(qi0(j,k,1,3)).ge.real(epss) .or. 
     .     real(qi0(j,k,5,3)).ge.real(epss)) .and.
     .   ib.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),200)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),65) nbl,igridg(nbl),j,k,real(qi0(j,k,1,3)),
     .                           real(qi0(j,k,5,3))
         if (real(qi0(j,k,1,3)).eq.-1.0 .and. 
     .       real(qi0(j,k,5,3)).eq.-1.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         else
            if (istop.gt.0) then
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
      end if
   65 format(6h block,i4,7h (grid ,i4,1h),23h  on i=idim boundary at,
     .5h j,k=,2i5,15h qi0(1),qi0(5)=,2e12.5)
  140 continue
      end if
c
  100 format(40h After analysis, most probable cause is:,
     .23h boundary data not set.)
  200 format(52h boundary conditions resulted in negative (or large),
     .24h density and/or pressure)
c
      return
      end
