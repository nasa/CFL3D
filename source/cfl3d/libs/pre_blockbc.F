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
      subroutine pre_blockbc(nbl,lw,icount,idimg,jdimg,kdimg,isav_blk,
     .                       nblk,nbli,limblk,isva,nblon,mxbli,nou,bou,
     .                       nbuf,ibufdim,myid,maxbl,ierrflg)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Set up data arrays for 1-1 bc message passing
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension nblk(2,mxbli),limblk(2,6,mxbli),isva(2,2,mxbli),
     .          nblon(mxbli),idimg(maxbl),jdimg(maxbl),kdimg(maxbl),
     .          isav_blk(2*mxbli,17),lw(65,maxbl)
c
      if(abs(nbli).gt.0) then 
c
c     loop over all 1-1 interfaces
c
      do 100 n=1,abs(nbli)
c
c        check to see if block interface n is turned on
c
         if (nblon(n) .ge. 0) then
c
c           if block interface n is a boundary of current block then
c     .     set up blocking parameters
c
            if (nbl.eq.nblk(1,n) .or. nbl.eq.nblk(2,n)) then
c
c              it = 1  if nblk(1,n) is the block being advanced
c              it = 2  if nblk(2,n) is the block being advanced
c              ir = 1  if nblk(1,n) is the neighboring block
c              ir = 2  if nblk(2,n) is the neighboring block
c
               it = 1
               ir = 2
               if (nbl.eq.nblk(2,n)) it = 2
               if (nbl.eq.nblk(2,n)) ir = 1
c
c              allow for 1-1 blocking in same grid
c
               itime = 1
               if (nblk(1,n).eq.nblk(2,n)) itime = 2
               do 101 iti = 1,itime
                  if (iti.gt.1) then
                     it = 1
                     ir = 2
                  end if
c
c                 define current and neighbor blocks
c
                  ic_blk = nblk(it,n)
                  in_blk = nblk(ir,n)
c
c                 set dimensions of current blocks involved
c
                  idimn = idimg(in_blk)
                  jdimn = jdimg(in_blk)
                  kdimn = kdimg(in_blk)
                  idimc = idimg(ic_blk)
                  jdimc = jdimg(ic_blk)
                  kdimc = kdimg(ic_blk)
c
                  if (isva(ir,1,n)+isva(ir,2,n) .eq. 3) then
                     jface = 3
                     jedge = 1
                     if (limblk(ir,3,n).ne.1) jedge = kdimg(in_blk)
                  end if
                  if (isva(ir,1,n)+isva(ir,2,n) .eq. 4) then
                     jface = 2
                     jedge = 1
                     if (limblk(ir,2,n).ne.1) jedge = jdimg(in_blk)
                  end if
                  if (isva(ir,1,n)+isva(ir,2,n) .eq. 5) then
                     jface = 1
                     jedge = 1
                     if (limblk(ir,1,n).ne.1) jedge = idimg(in_blk)
                  end if
c
c                 k = constant interface
c
                  if (isva(it,1,n)+isva(it,2,n) .eq. 3) then
                     lwt   = lw(3,ic_blk)
                     iedge = 1
                     if (limblk(it,3,n).ne.1) then
                        lwt   = lwt + jdimc*(idimc-1)*5*2
                        iedge = 2
                     end if
                     iss = limblk(it,1,n)
                     ise = limblk(it,4,n)
                     jss = limblk(it,2,n)
                     jse = limblk(it,5,n)
                     kss = limblk(it,3,n)
                     kse = limblk(it,6,n)
                     if (iss.gt.ise) then
                        iss=iss+1
                     else
                        ise=ise+1
                     end if
                     if (jss.gt.jse) then
                        jss=jss+1
                     else
                        jse=jse+1
                     end if
c
c                 j = constant interface 
c
                  else if (isva(it,1,n)+isva(it,2,n) .eq. 4) then
                     lwt   = lw(2,ic_blk)
                     iedge = 1
                     if (limblk(it,2,n).ne.1) then
                        lwt   = lwt + kdimc*(idimc-1)*5*2
                        iedge = 2
                     end if
                     iss = limblk(it,1,n)
                     ise = limblk(it,4,n)
                     jss = limblk(it,2,n)
                     jse = limblk(it,5,n)
                     kss = limblk(it,3,n)
                     kse = limblk(it,6,n)
                     if (iss.gt.ise) then
                        iss=iss+1
                     else
                        ise=ise+1
                     end if
                     if (kss.gt.kse) then
                        kss=kss+1
                     else
                        kse=kse+1
                     end if
c
c                 i = constant interface
c
                  else if (isva(it,1,n)+isva(it,2,n) .eq. 5) then
                     lwt   = lw(4,ic_blk)
                     iedge = 1
                     if (limblk(it,1,n).ne.1) then
                        lwt   = lwt + jdimc*kdimc*5*2
                        iedge = 2
                     end if
                     iss = limblk(it,1,n)
                     ise = limblk(it,4,n)
                     kss = limblk(it,3,n)
                     kse = limblk(it,6,n)
                     jss = limblk(it,2,n)
                     jse = limblk(it,5,n)
                     if (kss.gt.kse) then
                        kss=kss+1
                     else
                        kse=kse+1
                     end if
                     if (jss.gt.jse) then
                        jss=jss+1
                     else
                        jse=jse+1
                     end if
                  end if
c
c                 put appropriate values into isav_blk array
c
                  icount = icount + 1
                  if (icount.gt.2*mxbli) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*) 'err: icount in',
     .               ' pre_bcblk = ',icount
                     call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
                  end if
                  isav_blk(icount,1)  = n
                  isav_blk(icount,2)  = it
                  isav_blk(icount,3)  = ir
                  isav_blk(icount,4)  = ic_blk
                  isav_blk(icount,5)  = in_blk
                  isav_blk(icount,6)  = jface 
                  isav_blk(icount,7)  = jedge 
                  isav_blk(icount,8)  = isva(it,1,n)+isva(it,2,n)
                  isav_blk(icount,9)  = lwt
                  isav_blk(icount,10) = iedge
                  isav_blk(icount,11) = iss
                  isav_blk(icount,12) = ise
                  isav_blk(icount,13) = jss
                  isav_blk(icount,14) = jse
                  isav_blk(icount,15) = kss
                  isav_blk(icount,16) = kse
                  isav_blk(icount,17) = iti
c
  101          continue
            end if
         end if
  100 continue
c
      end if
c
      return
      end
