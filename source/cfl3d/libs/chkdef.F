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
      subroutine chkdef(nbl,idim,jdim,kdim,bci,bcj,bck,icsi,icsf,
     .                  jcsi,jcsf,kcsi,kcsf,nsegdfrm,idfrmseg,maxbl,
     .                  maxsegdg,nou,bou,nbuf,ibufdim,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Make sure all deforming surface segments do not have
c              interface type bc flags - note: this is not a 100%
c              check that the deforming surface segments have a solid
c              wall bc (as they must), since at one bc (1012) uses
c              the wall-type bc flag and is not a solid surface bc.
c              However, if of interface type, then there is definitely 
c              a problem
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg)
c
      do iseg=1,nsegdfrm(nbl)
c
c        check i=const surfaces
c
         if (icsi(nbl,iseg) .eq. icsf(nbl,iseg)) then
            if (icsi(nbl,iseg) .eq. 1) then
               ll = 1
            else
               ll = 2
            end if
            do j=jcsi(nbl,iseg),jcsf(nbl,iseg)-1
               do k=kcsi(nbl,iseg),kcsf(nbl,iseg)-1
                 if (bci(j,k,ll) .eq. 0.) then
                    nou(1) = min(nou(1)+1,ibufdim)
                    write(bou(nou(1),1),100) nbl,iseg
                    nou(1) = min(nou(1)+1,ibufdim)
                    write(bou(nou(1),1),101) j,k
                    call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                 end if
               end do
            end do
         end if
c
c        check j=const surfaces
c
         if (jcsi(nbl,iseg) .eq. jcsf(nbl,iseg)) then
            if (jcsi(nbl,iseg) .eq. 1) then
               ll = 1
            else
               ll = 2
            end if
            do i=icsi(nbl,iseg),icsf(nbl,iseg)-1
               do k=kcsi(nbl,iseg),kcsf(nbl,iseg)-1
                 if (bcj(k,i,ll) .eq. 0.) then
                    nou(1) = min(nou(1)+1,ibufdim)
                    write(bou(nou(1),1),100) nbl,iseg
                    nou(1) = min(nou(1)+1,ibufdim)
                    write(bou(nou(1),1),102) k,i
                    call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                 end if
               end do
            end do
         end if
c
c        check k=const surfaces
c
         if (kcsi(nbl,iseg) .eq. kcsf(nbl,iseg)) then
            if (kcsi(nbl,iseg) .eq. 1) then
               ll = 1
            else
               ll = 2
            end if
            do j=jcsi(nbl,iseg),jcsf(nbl,iseg)-1
               do i=icsi(nbl,iseg),icsf(nbl,iseg)-1
                 if (bck(j,i,ll) .eq. 0.) then
                    nou(1) = min(nou(1)+1,ibufdim)
                    write(bou(nou(1),1),100) nbl,iseg
                    nou(1) = min(nou(1)+1,ibufdim)
                    write(bou(nou(1),1),103) j,i
                    call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                 end if
               end do
            end do
         end if
      end do
c
 100  format('stopping...deforming face has interface type bc',
     .       ' block',i4,' segment',i3)
 101  format('  at j,k',2i4)
 102  format('  at k,i',2i4)
 103  format('  at j,i',2i4)
c
      return
      end 
