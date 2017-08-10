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
      subroutine chkrap(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .                 ibcinfo,jbcinfo,kbcinfo,nface,idim,jdim,kdim,
     .                 nn,maxxseg,ista,iend,jsta,jend,ksta,kend,iwrap,
     .                 jwrap,kwrap,nou,bou,nbuf,ibufdim,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Check for wraparound in j, k or i  directions in order 
c     to apply bctype 1012 (singular axis - full plane)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension nbci0(nn),nbcidim(nn),nbcj0(nn),nbcjdim(nn),
     .           nbck0(nn),nbckdim(nn),ibcinfo(nn,maxxseg,7,2),
     .           jbcinfo(nn,maxxseg,7,2),kbcinfo(nn,maxxseg,7,2)
c
      iokrap = 0
      iwrap = 0
      jwrap = 0
      kwrap = 0
c
c***********************************************************************
c     k=1 or k=kdim is a singular axis
c***********************************************************************
c
      if(nface.eq.5.or.nface.eq.6) then
c
c       check for wrap around in j direction, adjacent to k-axis, 
c       within same i-range as the singular axis segment
c
        jwrap = 1
        do 10 nseg=1,nbcj0(nbl)
        ista1 = jbcinfo(nbl,nseg,2,1)
        iend1 = jbcinfo(nbl,nseg,3,1)
        ksta1 = jbcinfo(nbl,nseg,4,1)
        kend1 = jbcinfo(nbl,nseg,5,1)
        if((nface.eq.5 .and. ksta1.eq.1) .or.
     .  (nface.eq.6 .and. kend1.eq.kdim)) then 
          if(ista1.lt.iend .and. iend1.gt.ista) then
c           this segment on the j=1 face must have a wrap around in 
c           order to apply bctype 1012 - turn OFF wraparound flag if not
            jbc = jbcinfo(nbl,nseg,1,1)
            if(jbc.ne.0) jwrap = 0
c
c           if(jbc.ne.0)
c    .        write(61,*) 'block ',nbl,' turning off jwrap on j0',
c    .        ' segment ', nseg
c
          end if
        end if
   10   continue
c
c       check for wrap around in j direction, adjacent to k-axis, 
c       within same i-range as the singular axis segment
c
        do 11 nseg=1,nbcjdim(nbl)
        ista1 = jbcinfo(nbl,nseg,2,2)
        iend1 = jbcinfo(nbl,nseg,3,2)
        ksta1 = jbcinfo(nbl,nseg,4,2)
        kend1 = jbcinfo(nbl,nseg,5,2)
        if((nface.eq.5 .and. ksta1.eq.1) .or.
     .  (nface.eq.6 .and. kend1.eq.kdim)) then
          if(ista1.lt.iend .and. iend1.gt.ista) then
c           this segment on the j=jdim face must have a wrap around in 
c           order to apply bctype 1012 - turn OFF wraparound flag if not
            jbc = jbcinfo(nbl,nseg,1,2)
            if(jbc.ne.0) jwrap = 0 
c  
c           if(jbc.ne.0)
c    .        write(61,*) 'block ',nbl,' turning off jwrap on jdim',
c    .        ' segment ', nseg
c
          end if
        end if
   11   continue
c
c       check for wrap around in i direction, adjacent to k-axis,
c       within same j-range as the singular axis segment
c
        iwrap = 1
        do 12 nseg=1,nbci0(nbl)
        jsta1 = ibcinfo(nbl,nseg,2,1)
        jend1 = ibcinfo(nbl,nseg,3,1)
        ksta1 = ibcinfo(nbl,nseg,4,1)
        kend1 = ibcinfo(nbl,nseg,5,1)
        if((nface.eq.5 .and. ksta1.eq.1) .or.
     .  (nface.eq.6 .and. kend1.eq.kdim)) then
          if(jsta1.lt.jend .and. jend1.gt.jsta) then
c           this segment on the i=1 face must have a wrap around in
c           order to apply bctype 1012 - turn OFF wraparound flag if not
            ibc = ibcinfo(nbl,nseg,1,1)
            if(ibc.ne.0) iwrap = 0
c   
c           if(ibc.ne.0)
c    .        write(61,*) 'block ',nbl,' turning off iwrap on i0',
c    .        ' segment ', nseg
c
          end if
        end if
   12   continue
c
c       check for wrap around in i direction, adjacent to k-axis,
c       within same j-range as the singular axis segment
c
        do 13 nseg=1,nbcidim(nbl)
        jsta1 = ibcinfo(nbl,nseg,2,2)
        jend1 = ibcinfo(nbl,nseg,3,2)
        ksta1 = ibcinfo(nbl,nseg,4,2)
        kend1 = ibcinfo(nbl,nseg,5,2)
        if((nface.eq.5 .and. ksta1.eq.1) .or.
     .  (nface.eq.6 .and. kend1.eq.kdim)) then
          if(jsta1.lt.jend .and. jend1.gt.jsta) then
c           this segment on the i=idim face must have a wrap around in
c           order to apply bctype 1012 - turn OFF wraparound flag if not
            ibc = ibcinfo(nbl,nseg,1,2)
            if(ibc.ne.0) iwrap = 0
c   
c           if(ibc.ne.0)
c    .        write(61,*) 'block ',nbl,' turning off iwrap on idim',
c    .        ' segment ', nseg
c
          end if
        end if
   13   continue
c
c       MUST have proper wrap around set in either i or j
c       directions in order to use bctype 1012 on k=1 or k=kdim
c
        if(jwrap.gt.0 .or. iwrap.gt.0) then
          iokrap = 1
        else
          if(nface.eq.5) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) ' error in applying bctype 1012 on ',
     .      'k=1 boundary:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) '    must have wrap around in i ',
     .      'or j directions'
          else 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) ' error in applying bctype 1012 on ',
     .      'k=kdim boundary:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) '    must have wrap around in i ',
     .      'or j directions'
          end if
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
c
      end if
c
c***********************************************************************
c     j=1 or j=jdim is a singular axis
c***********************************************************************
c
      if(nface.eq.3 .or. nface.eq.4) then
c
c       check for wrap around in k direction, adjacent to j-axis, 
c       within same i-range as the singular axis segment
c
        kwrap = 1
        do 30 nseg=1,nbck0(nbl)
        ista1 = kbcinfo(nbl,nseg,2,1)
        iend1 = kbcinfo(nbl,nseg,3,1)
        jsta1 = kbcinfo(nbl,nseg,4,1)
        jend1 = kbcinfo(nbl,nseg,5,1)
        if((nface.eq.3 .and. jsta1.eq.1) .or.
     .  (nface.eq.4 .and. jend1.eq.jdim)) then
          if(ista1.lt.iend .and. iend1.gt.ista) then
c           this segment on the k=1 face must have a wrap around in 
c           order to apply bctype 1012 - turn OFF wraparound flag if not
            kbc = kbcinfo(nbl,nseg,1,1)
            if(kbc.ne.0) kwrap = 0 
c
c           if(kbc.ne.0)
c    .        write(61,*) 'block ',nbl,' turning off kwrap on k0',
c    .        ' segment ', nseg
c
          end if
        end if
   30   continue
c
c       check for wrap around in k direction, adjacent to j-axis, 
c       within same i-range as the singular axis segment
c
        do 31 nseg=1,nbckdim(nbl)
        ista1 = kbcinfo(nbl,nseg,2,2)
        iend1 = kbcinfo(nbl,nseg,3,2)
        jsta1 = kbcinfo(nbl,nseg,4,2)
        jend1 = kbcinfo(nbl,nseg,5,2)
        if((nface.eq.3 .and. jsta1.eq.1) .or.
     .  (nface.eq.4 .and. jend1.eq.jdim)) then
          if(ista1.lt.iend .and. iend1.gt.ista) then
c           this segment on the k=1 face must have a wrap around in 
c           order to apply bctype 1012 - turn OFF wraparound flag if not
            kbc = kbcinfo(nbl,nseg,1,2)
            if(kbc.ne.0) kwrap = 0
c   
c           if(kbc.ne.0)
c    .         write(61,*) 'block ',nbl,' turning off kwrap on kdim',
c    .         ' segment ', nseg
c
          end if
        end if
   31   continue
c
c       check for wrap around in i direction, adjacent to j-axis,
c       within same k-range as the singular axis segment
c
        iwrap = 1
        do 32 nseg=1,nbci0(nbl)
        jsta1 = ibcinfo(nbl,nseg,2,1)
        jend1 = ibcinfo(nbl,nseg,3,1)
        ksta1 = ibcinfo(nbl,nseg,4,1)
        kend1 = ibcinfo(nbl,nseg,5,1)
        if((nface.eq.3 .and. jsta1.eq.1) .or.
     .  (nface.eq.4 .and. jend1.eq.jdim)) then
          if(ksta1.lt.kend .and. kend1.gt.ksta) then
c           this segment on the i=1 face must have a wrap around in
c           order to apply bctype 1012 - turn OFF wraparound flag if not
            ibc = ibcinfo(nbl,nseg,1,1)
            if(ibc.ne.0) iwrap = 0
c   
c           if(ibc.ne.0)
c    .        write(61,*) 'block ',nbl,' turning off iwrap on i0',
c    .        ' segment ', nseg
c
          end if
        end if
   32   continue
c
c       check for wrap around in i direction, adjacent to j-axis,
c       within same k-range as the singular axis segment
c
        do 33 nseg=1,nbcidim(nbl)
        jsta1 = ibcinfo(nbl,nseg,2,2)
        jend1 = ibcinfo(nbl,nseg,3,2)
        ksta1 = ibcinfo(nbl,nseg,4,2)
        kend1 = ibcinfo(nbl,nseg,5,2)
        if((nface.eq.3 .and. jsta1.eq.1) .or.
     .  (nface.eq.4 .and. jend1.eq.jdim)) then
          if(ksta1.lt.kend .and. kend1.gt.ksta) then
c           this segment on the i=1 face must have a wrap around in
c           order to apply bctype 1012 - turn OFF wraparound flag if not
            ibc = ibcinfo(nbl,nseg,1,2)
            if(ibc.ne.0) iwrap = 0
c   
c           if(ibc.ne.0)
c    .        write(61,*) 'block ',nbl,' turning off iwrap on idim',
c    .        ' segment ', nseg
c
          end if
        end if
   33   continue
c
c       MUST have proper wrap around set in either i or k
c       directions in order to use bctype 1012 on j=1 or j=jdim
c
        if(kwrap.gt.0 .or. iwrap.gt.0) then
          iokrap = 1
        else
          if(nface.eq.3) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) ' error in applying bctype 1012 on ',
     .      'j=1 boundary:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) '    must have wrap around in i ',
     .      'or k directions'
          else   
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) ' error in applying bctype 1012 on ',
     .      'j=jdim boundary:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) '    must have wrap around in i ',
     .      'or k directions'
          end if
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
c
      end if
c
c***********************************************************************
c     i=1 or i=idim is a singular axis
c***********************************************************************
c
      if(nface.eq.1 .or. nface.eq.2) then
c
c       check for wrap around in j direction, adjacent to i-axis, 
c       within same k-range as the singular axis segment
c
        jwrap = 1
        do 50 nseg=1,nbcj0(nbl)
        ista1 = jbcinfo(nbl,nseg,2,1)
        iend1 = jbcinfo(nbl,nseg,3,1)
        ksta1 = jbcinfo(nbl,nseg,4,1)
        kend1 = jbcinfo(nbl,nseg,5,1)
        if((nface.eq.1 .and. ista1.eq.1) .or.
     .  (nface.eq.2 .and. iend1.eq.idim)) then
          if(ksta1.lt.kend .and. kend1.gt.ksta) then
c           this segment on the j=1 face must have a wrap around in 
c           order to apply bctype 1012 - turn OFF wraparound flag if not
            jbc = jbcinfo(nbl,nseg,1,1)
            if(jbc.ne.0) jwrap = 0
c
c           if(jbc.ne.0)
c    .        write(61,*) 'block ',nbl,' turning off jwrap on j0',
c    .        ' segment ', nseg
c
          end if
        end if
   50   continue
c
c       check for wrap around in j direction, adjacent to i-axis, 
c       within same k-range as the singular axis segment
c
        do 51 nseg=1,nbcjdim(nbl)
        ista1 = jbcinfo(nbl,nseg,2,2)
        iend1 = jbcinfo(nbl,nseg,3,2)
        ksta1 = jbcinfo(nbl,nseg,4,2)
        kend1 = jbcinfo(nbl,nseg,5,2)
        if((nface.eq.1 .and. ista1.eq.1) .or.
     .  (nface.eq.2 .and. iend1.eq.idim)) then
         if(ksta1.lt.kend .and. kend1.gt.ksta) then
c          this segment on the j=jdim face must have a wrap around in 
c          order to apply bctype 1012 - turn OFF wraparound flag if not
           jbc = jbcinfo(nbl,nseg,1,2)
           if(jbc.ne.0) jwrap = 0 
c   
c          if(jbc.ne.0)
c    .        write(61,*) 'block ',nbl,' turning off jwrap on jdim',
c    .        ' segment ', nseg
c
         end if
        end if
   51   continue
c
c       check for wrap around in k direction, adjacent to i-axis,
c       within same j-range as the singular axis segment
c
        kwrap = 1
        do 52 nseg=1,nbck0(nbl)
        ista1 = kbcinfo(nbl,nseg,2,1)
        iend1 = kbcinfo(nbl,nseg,3,1)
        jsta1 = kbcinfo(nbl,nseg,4,1)
        jend1 = kbcinfo(nbl,nseg,5,1)
        if((nface.eq.1 .and. ista1.eq.1) .or.
     .  (nface.eq.2 .and. iend1.eq.idim)) then
          if(jsta1.lt.jend .and. jend1.gt.jsta) then
c           this segment on the k=1 face must have a wrap around in
c           order to apply bctype 1012 - turn OFF wraparound flag if not
            kbc = kbcinfo(nbl,nseg,1,1)
            if(kbc.ne.0) kwrap = 0
c   
c           if(kbc.ne.0)
c    .        write(61,*) 'block ',nbl,' turning off kwrap on k0',
c    .        ' segment ', nseg
c
          end if
        end if
   52   continue
c
c       check for wrap around in k direction, adjacent to i-axis,
c       within same j-range as the singular axis segment
c
        do 53 nseg=1,nbckdim(nbl)
        ista1 = kbcinfo(nbl,nseg,2,2)
        iend1 = kbcinfo(nbl,nseg,3,2)
        jsta1 = kbcinfo(nbl,nseg,4,2)
        jend1 = kbcinfo(nbl,nseg,5,2)
        if((nface.eq.1 .and. ista1.eq.1) .or.
     .  (nface.eq.2 .and. iend1.eq.idim)) then
          if(jsta1.lt.jend .and. jend1.gt.jsta) then
c           this segment on the k=kdim face must have a wrap around in
c           order to apply bctype 1012 - turn OFF wraparound flag if not
            kbc = kbcinfo(nbl,nseg,1,2)
            if(kbc.ne.0) kwrap = 0
c   
c           if(kbc.ne.0)
c    .        write(61,*) 'block ',nbl,' turning off kwrap on kdim',
c    .        ' segment ', nseg
c
          end if
        end if
   53   continue
c
c       MUST have proper wrap around set in either k or j
c       directions in order to use bctype 1012 on i=1 or i=idim
c
        if(jwrap.gt.0 .or. kwrap.gt.0) then
          iokrap = 1
        else
          if(nface.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) ' error in applying bctype 1012 on ',
     .      'i=1 boundary:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) '    must have wrap around in k ',
     .      'or j directions'
          else   
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) ' error in applying bctype 1012 on ',
     .      'i=idim boundary:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) '    must have wrap around in k ',
     .      'or j directions'
          end if
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
c
      end if
c
      return
      end
