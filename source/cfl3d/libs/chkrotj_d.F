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
      subroutine chkrotj_d(nbl,jdim,kdim,idim,x,y,z,nblp,jdimp,kdimp,
     .                  idimp,nface,bcdata,xyzjp,ista,iend,jsta,jend,
     .                  ksta,kend,mdim,ndim,lcnt,xorig,yorig,zorig,
     .                  maxbl,period_miss,lbcprd,nou,bou,nbuf,ibufdim,
     .                  myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Check to make sure that the proper rotation angle for
c     periodic boundary conditions has been set. The check is done by
c     rotating the periodic block face through the specified angle, then
c     checking for point-to-point match with the current block face.
c
c     if j=1 is the periodic face in block nbl, it is assumed periodic 
c     with j=jdim in block nblp. if j=jdim is the periodic face in block
c     nbl, it is assumed periodic with j=1  in block nblp. furthermore,
c     it is assumed that i and k run in the same direction in blocks
c     nbl and nblp. similar assumptions hold when i or k faces are 
c     periodic. nbl and nblp may be identical, but they need not be.
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
      dimension xyzjp(kdim,idim,3)
      dimension bcdata(mdim,ndim,2,12),period_miss(lbcprd)
      dimension xorig(maxbl),yorig(maxbl),zorig(maxbl)
c
      common /conversion/ radtodeg
c
      ip = 1
c
      dthtx = bcdata(1,1,ip,2)/radtodeg
      dthty = bcdata(1,1,ip,3)/radtodeg
      dthtz = bcdata(1,1,ip,4)/radtodeg
c
c     check that only one of dthtx,dthty,dthtz is non-zero
c
      if (dthtx.ne.0) then
	 if (dthty.ne.0. .or. dthtz.ne.0.) then 
            nou(1) = min(nou(1)+1,ibufdim)
	    write(bou(nou(1),1),101)
	    call termn8(myid,-1,ibufdim,nbuf,bou,nou)
	 end if
      end if
      if (dthty.ne.0) then
	 if (dthtx.ne.0. .or. dthtz.ne.0.) then
            nou(1) = min(nou(1)+1,ibufdim)
	    write(bou(nou(1),1),101)
	    call termn8(myid,-1,ibufdim,nbuf,bou,nou)
	 end if
      end if
      if (dthtz.ne.0) then
	 if (dthtx.ne.0. .or. dthty.ne.0.) then
            nou(1) = min(nou(1)+1,ibufdim)
	    write(bou(nou(1),1),101)
	    call termn8(myid,-1,ibufdim,nbuf,bou,nou)
	 end if
      end if
c
      eps = 0
c
c     j=1 boundary
c
      if (nface.eq.3) then
c
         j  = 1
         jp = jdimp
c
c        check that idim = idimp and kdim = kdimp
c
         if (idim.ne.idimp .or. kdim.ne.kdimp) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),102)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),104) idim,kdim,idimp,kdimp
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
c        rotate periodic block face (entire face)
c
         jdum = 1
         kdum = kdimp
         idum = idimp
         call grdmove(nblp,jdum,kdum,idum,xyzjp(1,1,1),xyzjp(1,1,2),
     .                xyzjp(1,1,3),xorig(nblp),yorig(nblp),zorig(nblp),
     .                xorig(nblp),yorig(nblp),zorig(nblp),dthtx,dthty,
     .                dthtz)
c
c        check for point match
c
         do 3010 i=ista,iend
         do 3010 k=ksta,kend
         xt1 = x(j,k,i)
         yt1 = y(j,k,i)
         zt1 = z(j,k,i)
         xr1 = xyzjp(k,i,1)
         yr1 = xyzjp(k,i,2)
         zr1 = xyzjp(k,i,3)
         eps = ccmax(eps,sqrt((xr1-xt1)**2+(yr1-yt1)**2+(zr1-zt1)**2))
 3010    continue
c
         period_miss(lcnt) = eps
c
      end if
c
c     j=jdim boundary
c
      if (nface.eq.4) then
c
         j  = jdim
         jp = 1
c
c        check that idim = idimp and kdim = kdimp
c
         if (idim.ne.idimp .or. kdim.ne.kdimp) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),102)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),104) idim,kdim,idimp,kdimp
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
c        rotate periodic block face (entire face)
c
         jdum = 1
         kdum = kdimp
         idum = idimp
         call grdmove(nblp,jdum,kdum,idum,xyzjp(1,1,1),xyzjp(1,1,2),
     .                xyzjp(1,1,3),xorig(nblp),yorig(nblp),zorig(nblp),
     .                xorig(nblp),yorig(nblp),zorig(nblp),dthtx,dthty,
     .                dthtz)
c
c        check for point match
c
         do 4010 i=ista,iend
         do 4010 k=ksta,kend
         xt1 = x(j,k,i)
         yt1 = y(j,k,i)
         zt1 = z(j,k,i)
         xr1 = xyzjp(k,i,1)
         yr1 = xyzjp(k,i,2)
         zr1 = xyzjp(k,i,3)
         eps = ccmax(eps,sqrt((xr1-xt1)**2+(yr1-yt1)**2+(zr1-zt1)**2))
 4010    continue
c
         period_miss(lcnt) = eps
c
      end if
c
 101  format(1x,41h2 of the 3 dtht values currently must = 0)
 102  format(1x,'periodic block face must be of the same 2 ',
     .          'dimensions (and orientation) as the current ',
     .          'block face')
 104  format(1x,24hidim,kdim,idimp,kdimp = ,4i5)
      return
      end
