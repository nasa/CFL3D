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
      subroutine xupdt (q,qj0,qk0,qi0,jdim,kdim,idim,nbl,ldim,bcj0,
     .                  bck0,bci0,maxbl,iitot,iibg,kkbg,jjbg,ibcg,
     .                  lbg,ibpntsg,qb,nou,bou,nbuf,ibufdim,int_updt)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Update fringe points of overlapped grids with boundary
c     values which have been interpolated from other grids to provide 
c     the mechanism for coupling the various grids. 
c
c     int_updt = 0 prevents interior (fringe) points from being 
c                  updated. this switch was needed with the new qout
c                  routine of version 6 in order to match cell-center
c                  output data with version 5
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,ldim), qi0(jdim,kdim,ldim,4),
     .          qj0(kdim,idim-1,ldim,4),qk0(jdim,idim-1,ldim,4)
      dimension bci0(jdim,kdim,2),bcj0(kdim,idim-1,2),
     .          bck0(jdim,idim-1,2)
      dimension iibg(iitot),kkbg(iitot),jjbg(iitot),ibcg(iitot),
     .          lbg(maxbl),ibpntsg(maxbl,4),qb(iitot,5,3)
c
      common /sklton/ isklton
c
c     temporarily set isklton = 0 to supress output messages...
c     the messages are now output from bc_info
c
      isklt_sav = isklton
      isklton   = 0
c
c   Set index for loading into appropriate qb array 
c        (1 for q, 2 for vk0,vj0,vi0, 3 for tj0,tk0,ti0)
      if (ldim .eq. 5) then
        iset=1
      else if (ldim .eq. 1) then
        iset=2
      else
        iset=3
      end if
c
c     if (isklton.eq.1)
c    . write(15,*) ' chimera grid interpolations, block ',nbl
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
c
      lsta = lbg(nbl)
      lend = lsta-1
c
      if (ibpntsg(nbl,1).gt.0 .and. int_updt.ne.0) then
         lend = lsta+ibpntsg(nbl,1)-1
         do 11 l=lsta,lend
         do 10 ll=1,ldim
         q(jjbg(l),kkbg(l),iibg(l),ll) = qb(ibcg(l),ll,iset) 
   10    continue
   11    continue
         if (isklton.eq.1) then
c           write(15,*)' interior values updated = ',ibpntsg(nbl,1)
            if (iset.eq.1) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1600)
            end if
         end if
      end if
c
      if (ibpntsg(nbl,2).gt.0) then
         mp0   = 0
         mpdim = 0
         lsta = lend+1
         lend = lsta+ibpntsg(nbl,2)-1
         do 21 l=lsta,lend
         mp = max(0,1+jjbg(l)-jdim1) + 1 + max(0,-jjbg(l))
         mpp = max(0,3-mp)
         mp0 = mp0 + min(1,mpp)
         mpp = max(0,mp-2)
         mpdim = mpdim + min(1,mpp)
         mm = mp/3 + 1
         do 20 ll=1,ldim
         qj0(kkbg(l),iibg(l),ll,mp) = qb(ibcg(l),ll,iset)
         bcj0(kkbg(l),iibg(l),mm) = 0.0
   20    continue
   21    continue
         if (isklton.eq.1) then
c           write(15,*)' QJ0 values updated = ',ibpntsg(nbl,2)
            if (mp0.gt.0 .and. iset.eq.1) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1603)
            end if
            if (mpdim.gt.0 .and. iset.eq.1) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1604)
            end if
         end if
      end if
c
      if (ibpntsg(nbl,3).gt.0) then
         mp0   = 0
         mpdim = 0
         lsta = lend+1
         lend = lsta+ibpntsg(nbl,3)-1
         do 31 l=lsta,lend
         mp = max(0,1+kkbg(l)-kdim1) + 1 + max(0,-kkbg(l))
         mpp = max(0,3-mp)
         mp0 = mp0 + min(1,mpp)
         mpp = max(0,mp-2)
         mpdim = mpdim + min(1,mpp)
         mm = mp/3 + 1
         do 30 ll=1,ldim
         qk0(jjbg(l),iibg(l),ll,mp) = qb(ibcg(l),ll,iset)
         bck0(jjbg(l),iibg(l),mm) = 0.0
   30    continue
   31    continue
         if (isklton.eq.1) then
c           write(15,*)' QK0 values updated = ',ibpntsg(nbl,3)
            if (mp0.gt.0 .and. iset.eq.1) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1605)
            end if
            if (mpdim.gt.0 .and. iset.eq.1) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1606)
            end if
         end if
      end if
c
      if (ibpntsg(nbl,4).gt.0) then
         mp0   = 0
         mpdim = 0
         lsta = lend+1
         lend = lsta+ibpntsg(nbl,4)-1
         do 41 l=lsta,lend
         mp = max(0,1+iibg(l)-idim1) + 1 + max(0,-iibg(l))
         mpp = max(0,3-mp)
         mp0 = mp0 + min(1,mpp)
         mpp = max(0,mp-2)
         mpdim = mpdim + min(1,mpp)
         mm = mp/3 + 1
         do 40 ll=1,ldim
         qi0(jjbg(l),kkbg(l),ll,mp) = qb(ibcg(l),ll,iset)
         bci0(jjbg(l),kkbg(l),mm) = 0.0
   40    continue
   41    continue
         if (isklton.eq.1) then
c           write(15,*)' QI0 values updated = ',ibpntsg(nbl,4)
            if (mp0.gt.0 .and. iset.eq.1) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1601)
            end if
            if (mpdim.gt.0 .and. iset.eq.1) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1602)
            end if
         end if
      end if
c
      isklton = isklt_sav
c
 1600 format(' ','  hole    chimera grid interpolation     type    0')
 1601 format(' ','  i=   1  chimera grid interpolation     type    0')
 1602 format(' ','  i=idim  chimera grid interpolation     type    0')
 1603 format(' ','  j=   1  chimera grid interpolation     type    0')
 1604 format(' ','  j=jdim  chimera grid interpolation     type    0')
 1605 format(' ','  k=   1  chimera grid interpolation     type    0')
 1606 format(' ','  k=kdim  chimera grid interpolation     type    0')
c
      return
      end 
