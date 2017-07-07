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
      subroutine rpatch0(intmx,msub1,iindex,ninter)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Read in generalized-coordinate interpolation data
c     from a file; only the data needed to evaluate sizing 
c     requirements are stored at this time (this is a modified
c     version of subroutine rpatch)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      real dum
c
      common /params/ lmaxgr,lmaxbl,lmxseg,lmaxcs,lnplts,lmxbli,lmaxxe,
     .                lnsub1,lintmx,lmxxe,liitot,isum,lncycm,
     .                isum_n,lminnode,isumi,isumi_n,lmptch,
     .                lmsub1,lintmax,libufdim,lnbuf,llbcprd,
     .                llbcemb,llbcrad,lnmds,lmaxaes,lnslave,lmxsegdg,
     .                lnmaster
c
      dimension iindex(intmx,2*msub1+9)
c
c     read patched-grid connectivity file
c
      read(22) ninter
      lintmax = ninter
      lnsub1 = 0
      do 1500 n=1,abs(ninter)
      read(22) iindex(n,1)
      nfb = iindex(n,1)
      lnsub1 = max(lnsub1,nfb)
      read(22) (iindex(n,1+ll),ll=1,nfb)
      read(22) iindex(n,nfb+2)
      idum1 = iindex(n,nfb+2)
      read(22) (iindex(n,nfb+2+ll),ll=1,nfb)
      read(22) iindex(n,2*nfb+3)
      read(22) iindex(n,2*nfb+4)
      idum2 = iindex(n,2*nfb+4)
      read(22) iindex(n,2*nfb+5)
      idum3 = iindex(n,2*nfb+5)
      read(22) iindex(n,2*nfb+6)
      read(22) iindex(n,2*nfb+7)
      read(22) iindex(n,2*nfb+8)
      read(22) iindex(n,2*nfb+9)
      nbl = idum1
      lst  =idum3
      len = lst + idum2 - 1
      read(22) (ndum,nnn=lst,len)
      read(22) ((dum,nnn=lst,len),ll=1,2)
 1500 continue
c
c     don't need to read angular displacements for
c     sizing purposes
c
      lmaxxe = len
c
c     set ninter = -ninter to allow full read of
c     the patch file in subroutine setup
c
      ninter = -ninter
c
      return
      end
