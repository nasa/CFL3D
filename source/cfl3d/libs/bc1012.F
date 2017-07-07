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
      subroutine bc1012(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,iwrap,jwrap,kwrap,iuns,
     .                  nou,bou,nbuf,ibufdim,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set singular axis - full plane boundary conditions
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
c
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c     works for full-plane only - assume checks for appropriateness
c     of this boundary condition have been made PRIOR to entering this routine
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=1012 *
c            * * * * * * * * * * * * * * * * * * * * * *
c******************************************************************************
c      j=1 boundary            singular axis - full plane           bctype 1012
c******************************************************************************
c
      if (nface.eq.3) then
c
c     wraparound in i-direction
c
      if(iwrap.gt.0) then
        do 10 k=ksta,kend1
        do 10 l=1,5
        sumq = 0.0
        do 11 i=ista,iend1
        sumq = sumq+q(1,k,i,l)
11      continue
        sumq = sumq/float(iend1-ista+1)
        do 10 i=ista,iend1
        qj0(k,i,l,1) = sumq
        qj0(k,i,l,2) = (q(1,k,i,l)-sumq)*2.0
        bcj(k,i,1)   = 1.0
10      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 101 k=ksta,kend1
          sumvis = 0.0
          do 102 i=ista,iend1
          sumvis = sumvis+vist3d(1,k,i)
102       continue
          sumvis = sumvis/float(iend1-ista+1)
          do 101 i=ista,iend1
            vj0(k,i,1,1) = sumvis
            vj0(k,i,1,2) = (vist3d(1,k,i)-sumvis)*2.0
101       continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 104 k=ksta,kend1
          do 104 l=1,nummem
          sumtur = 0.0
          do 105 i=ista,iend1
          sumtur = sumtur + tursav(1,k,i,l)
105       continue
          sumtur = sumtur/float(iend1-ista+1)
          do 104 i=ista,iend1
            tj0(k,i,l,1) = sumtur
            tj0(k,i,l,2) = (tursav(1,k,i,l)-sumtur)*2.0
104       continue
        end if
        end if
c
      end if
c
c     wraparound in k-direction
c
      if(kwrap.gt.0) then
        do 13 i=ista,iend1
        do 13 l=1,5
        sumq = 0.0
        do 14 k=ksta,kend1
        sumq = sumq+q(1,k,i,l)
14      continue
        sumq = sumq/float(kend1-ksta+1)
        do 13 k=ksta,kend1
        qj0(k,i,l,1) = sumq
        qj0(k,i,l,2) = (q(1,k,i,l)-sumq)*2.0
        bcj(k,i,1)   = 1.0
13      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 106 i=ista,iend1
          sumvis = 0.0
          do 107 k=ksta,kend1
            sumvis = sumvis+vist3d(1,k,i)
107       continue
          sumvis = sumvis/float(kend1-ksta+1)
          do 106 k=ksta,kend1
            vj0(k,i,1,1) = sumvis
            vj0(k,i,1,2) = (vist3d(1,k,i)-sumvis)*2.0
106       continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 108 i=ista,iend1
          do 108 l=1,nummem
          sumtur = 0.0
          do 109 k=ksta,kend1
          sumtur = sumtur + tursav(1,k,i,l)
109       continue
          sumtur = sumtur/float(kend1-ksta+1)
          do 108 k=ksta,kend1
            tj0(k,i,l,1) = sumtur
            tj0(k,i,l,2) = (tursav(1,k,i,l)-sumtur)*2.0
108       continue
        end if
        end if
c
      end if
c
      end if
c******************************************************************************
c      j=jdim boundary         singular axis - full plane           bctype 1012
c******************************************************************************
c
      if (nface.eq.4) then
c
c     wraparound in i-direction
c
      if(iwrap.gt.0) then
        do 20 k=ksta,kend1
        do 20 l=1,5
        sumq = 0.0
        do 21 i=ista,iend1
        sumq = sumq+q(jdim1,k,i,l)
21      continue
        sumq = sumq/float(iend1-ista+1)
        do 20 i=ista,iend1
        qj0(k,i,l,3) = sumq
        qj0(k,i,l,4) = (sumq-q(jdim1,k,i,l))*2.0
        bcj(k,i,2)   = 1.0
20      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 201 k=ksta,kend1
          sumvis = 0.0
          do 202 i=ista,iend1
          sumvis = sumvis+vist3d(jdim1,k,i)
202       continue
          sumvis = sumvis/float(iend1-ista+1)
          do 201 i=ista,iend1
            vj0(k,i,1,3) = sumvis
            vj0(k,i,1,4) = (sumvis-vist3d(jdim1,k,i))*2.0
201       continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 204 k=ksta,kend1
          do 204 l=1,nummem
          sumtur = 0.0
          do 205 i=ista,iend1
          sumtur = sumtur + tursav(jdim1,k,i,l)
205       continue
          sumtur = sumtur/float(iend1-ista+1)
          do 204 i=ista,iend1
            tj0(k,i,l,3) = sumtur
            tj0(k,i,l,4) = (sumtur-tursav(jdim1,k,i,l))*2.0
204       continue
        end if
        end if
c
      end if
c
c     wraparound in k-direction
c
      if(kwrap.gt.0) then
        do 23 i=ista,iend1
        do 23 l=1,5
        sumq = 0.0
        do 24 k=ksta,kend1
        sumq = sumq+q(jdim1,k,i,l)
24      continue
        sumq = sumq/float(kend1-ksta+1)
        do 23 k=ksta,kend1
        qj0(k,i,l,3) = sumq
        qj0(k,i,l,4) = (sumq-q(jdim1,k,i,l))*2.0
        bcj(k,i,2)   = 1.0
23      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 206 i=ista,iend1
          sumvis = 0.0
          do 207 k=ksta,kend1
            sumvis = sumvis+vist3d(jdim1,k,i)
207       continue
          sumvis = sumvis/float(kend1-ksta+1)
          do 206 k=ksta,kend1
            vj0(k,i,1,3) = sumvis
            vj0(k,i,1,4) = (sumvis-vist3d(jdim1,k,i))*2.0
206       continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 208 i=ista,iend1
          do 208 l=1,nummem
          sumtur = 0.0
          do 209 k=ksta,kend1
          sumtur = sumtur + tursav(jdim1,k,i,l)
209       continue
          sumtur = sumtur/float(kend1-ksta+1)
          do 208 k=ksta,kend1
            tj0(k,i,l,3) = sumtur
            tj0(k,i,l,4) = (sumtur-tursav(jdim1,k,i,l))*2.0
208       continue
        end if
        end if
c
      end if
c
      end if
c******************************************************************************
c      k=1 boundary            singular axis - full plane           bctype 1012
c******************************************************************************
c
      if (nface.eq.5) then
c
c     wraparound in i-direction
c
      if(iwrap.gt.0) then
        do 30 j=jsta,jend1
        do 30 l=1,5
        sumq = 0.0
        do 31 i=ista,iend1
        sumq = sumq+q(j,1,i,l)
31      continue
        sumq = sumq/float(iend1-ista+1)
        do 30 i=ista,iend1
        qk0(j,i,l,1) = sumq
        qk0(j,i,l,2) = (q(j,1,i,l)-sumq)*2.0
        bck(j,i,1)   = 1.0
30      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 301 j=jsta,jend1
          sumvis = 0.0
          do 302 i=ista,iend1
          sumvis = sumvis+vist3d(j,1,i)
302       continue
          sumvis = sumvis/float(iend1-ista+1)
          do 301 i=ista,iend1
            vk0(j,i,1,1) = sumvis
            vk0(j,i,1,2) = (vist3d(j,1,i)-sumvis)*2.0
301       continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 304 j=jsta,jend1
          do 304 l=1,nummem
          sumtur = 0.0
          do 305 i=ista,iend1
          sumtur = sumtur + tursav(j,1,i,l)
305       continue
          sumtur = sumtur/float(iend1-ista+1)
          do 304 i=ista,iend1
            tk0(j,i,l,1) = sumtur
            tk0(j,i,l,2) = (tursav(j,1,i,l)-sumtur)*2.0
304       continue
        end if
        end if
c
      end if
c
c     wraparound in j-direction
c
      if(jwrap.gt.0) then
        do 33 i=ista,iend1
        do 33 l=1,5
        sumq = 0.0
        do 34 j=jsta,jend1
        sumq = sumq+q(j,1,i,l)
34      continue
        sumq = sumq/float(jend1-jsta+1)
        do 33 j=jsta,jend1
        qk0(j,i,l,1) = sumq
        qk0(j,i,l,2) = (q(j,1,i,l)-sumq)*2.0
        bck(j,i,1)   = 1.0
33      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 306 i=ista,iend1
          sumvis = 0.0
          do 307 j=jsta,jend1
            sumvis = sumvis+vist3d(j,1,i)
307       continue
          sumvis = sumvis/float(jend1-jsta+1)
          do 306 j=jsta,jend1
            vk0(j,i,1,1) = sumvis
            vk0(j,i,1,2) = (vist3d(j,1,i)-sumvis)*2.0
306       continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 308 i=ista,iend1
          do 308 l=1,nummem
          sumtur = 0.0
          do 309 j=jsta,jend1
          sumtur = sumtur + tursav(j,1,i,l)
309       continue
          sumtur = sumtur/float(jend1-jsta+1)
          do 308 j=jsta,jend1
            tk0(j,i,l,1) = sumtur
            tk0(j,i,l,2) = (tursav(j,1,i,l)-sumtur)*2.0
308       continue
        end if
        end if
c
      end if
c
      end if
c******************************************************************************
c      k=kdim boundary            singular axis - full plane        bctype 1012
c******************************************************************************
c
      if (nface.eq.6) then
c
c     wraparound in i-direction
c
      if(iwrap.gt.0) then
        do 40 j=jsta,jend1
        do 40 l=1,5
        sumq = 0.0
        do 41 i=ista,iend1
        sumq = sumq+q(j,kdim1,i,l)
41      continue
        sumq = sumq/float(iend1-ista+1)
        do 40 i=ista,iend1
        qk0(j,i,l,3) = sumq
        qk0(j,i,l,4) = (sumq-q(j,kdim1,i,l))*2.0
        bck(j,i,2)   = 1.0
40      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 401 j=jsta,jend1
          sumvis = 0.0
          do 402 i=ista,iend1
          sumvis = sumvis+vist3d(j,kdim1,i)
402       continue
          sumvis = sumvis/float(iend1-ista+1)
          do 401 i=ista,iend1
            vk0(j,i,1,3) = sumvis
            vk0(j,i,1,4) = (sumvis-vist3d(j,kdim1,i))*2.0
401       continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 404 j=jsta,jend1
          do 404 l=1,nummem
          sumtur = 0.0
          do 405 i=ista,iend1
          sumtur = sumtur + tursav(j,kdim1,i,l)
405       continue
          sumtur = sumtur/float(iend1-ista+1)
          do 404 i=ista,iend1
            tk0(j,i,l,3) = sumtur
            tk0(j,i,l,4) = (sumtur-tursav(j,kdim1,i,l))*2.0
404       continue
        end if
        end if
c
      end if
c
c     wraparound in j-direction
c
      if(jwrap.gt.0) then
        do 43 i=ista,iend1
        do 43 l=1,5
        sumq = 0.0
        do 44 j=jsta,jend1
        sumq = sumq+q(j,kdim1,i,l)
44      continue
        sumq = sumq/float(jend1-jsta+1)
        do 43 j=jsta,jend1
        qk0(j,i,l,3) = sumq
        qk0(j,i,l,4) = (sumq-q(j,kdim1,i,l))*2.0
        bck(j,i,2)   = 1.0
43      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 406 i=ista,iend1
          sumvis = 0.0
          do 407 j=jsta,jend1
            sumvis = sumvis+vist3d(j,kdim1,i)
407       continue
          sumvis = sumvis/float(jend1-jsta+1)
          do 406 j=jsta,jend1
            vk0(j,i,1,3) = sumvis
            vk0(j,i,1,4) = (sumvis-vist3d(j,kdim1,i))*2.0
406       continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 408 i=ista,iend1
          do 408 l=1,nummem
          sumtur = 0.0
          do 409 j=jsta,jend1
          sumtur = sumtur + tursav(j,kdim1,i,l)
409       continue
          sumtur = sumtur/float(jend1-jsta+1)
          do 408 j=jsta,jend1
            tk0(j,i,l,3) = sumtur
            tk0(j,i,l,4) = (sumtur-tursav(j,kdim1,i,l))*2.0
408       continue
        end if
        end if
c
      end if
c
      end if
c******************************************************************************
c      i=1 boundary            singular axis - full plane           bctype 1012
c******************************************************************************
c
      if (nface.eq.1) then
c
c     wraparound in j-direction
c
      if(jwrap.gt.0) then
        do 50 k=ksta,kend1
        do 50 l=1,5
        sumq = 0.0
        do 51 j=jsta,jend1
        sumq = sumq+q(j,k,1,l)
51      continue
        sumq = sumq/float(jend1-jsta+1)
        do 50 j=jsta,jend1
        qi0(j,k,l,1) = sumq
        qi0(j,k,l,2) = (q(j,k,1,l)-sumq)*2.0
        bci(j,k,1)   = 1.0
50      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 501 k=ksta,kend1
          sumvis = 0.0
          do 502 j=jsta,jend1
          sumvis = sumvis+vist3d(j,k,1)
502       continue
          sumvis = sumvis/float(jend1-jsta+1)
          do 501 j=jsta,jend1
            vi0(j,k,1,1) = sumvis
            vi0(j,k,1,2) = (vist3d(j,k,1)-sumvis)*2.0
501       continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 504 k=ksta,kend1
          do 504 l=1,nummem
          sumtur = 0.0
          do 505 j=jsta,jend1
          sumtur = sumtur + tursav(j,k,1,l)
505       continue
          sumtur = sumtur/float(jend1-jsta+1)
          do 504 j=jsta,jend1
            ti0(j,k,l,1) = sumtur
            ti0(j,k,l,2) = (tursav(j,k,1,l)-sumtur)*2.0
504       continue
        end if
        end if
c
      end if
c
c     wraparound in k-direction
c
      if(kwrap.gt.0) then
        do 53 j=jsta,jend1
        do 53 l=1,5
        sumq = 0.0
        do 54 k=ksta,kend1
        sumq = sumq+q(j,k,1,l)
54      continue
        sumq = sumq/float(kend1-ksta+1)
        do 53 k=ksta,kend1
        qi0(j,k,l,1) = sumq
        qi0(j,k,l,2) = (q(j,k,1,l)-sumq)*2.0
        bci(j,k,1)   = 1.0
53      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 506 j=jsta,jend1
          sumvis = 0.0
          do 507 k=ksta,kend1
            sumvis = sumvis+vist3d(j,k,1)
507       continue
          sumvis = sumvis/float(jend1-jsta+1)
          do 506 k=ksta,kend1
            vi0(j,k,1,1) = sumvis
            vi0(j,k,1,2) = (vist3d(j,k,1)-sumvis)*2.0
506       continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 508 j=jsta,jend1
          do 508 l=1,nummem
          sumtur = 0.0
          do 509 k=ksta,kend1
          sumtur = sumtur + tursav(j,k,1,l)
509       continue
          sumtur = sumtur/float(jend1-jsta+1)
          do 508 k=ksta,kend1
            ti0(j,k,l,1) = sumtur
            ti0(j,k,l,2) = (tursav(j,k,1,l)-sumtur)*2.0
508       continue
        end if
        end if
c
      end if
c
      end if
c******************************************************************************
c      i=idim boundary         singular axis - full plane           bctype 1012
c******************************************************************************
c
      if (nface.eq.2) then
c
c     wraparound in j-direction
c
      if(jwrap.gt.0) then
        do 60 k=ksta,kend1
        do 60 l=1,5
        sumq = 0.0
        do 61 j=jsta,jend1
        sumq = sumq+q(j,k,idim1,l)
61      continue
        sumq = sumq/float(jend1-jsta+1)
        do 60 j=jsta,jend1
        qi0(j,k,l,3) = sumq
        qi0(j,k,l,4) = (sumq-q(j,k,idim1,l))*2.0
        bci(j,k,2)   = 1.0
60      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 601 k=ksta,kend1
          sumvis = 0.0
          do 602 j=jsta,jend1
          sumvis = sumvis+vist3d(j,k,idim1)
602       continue
          sumvis = sumvis/float(jend1-jsta+1)
          do 601 j=jsta,jend1
            vi0(j,k,1,3) = sumvis
            vi0(j,k,1,4) = (sumvis-vist3d(j,k,idim1))*2.0
601       continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 604 k=ksta,kend1
          do 604 l=1,nummem
          sumtur = 0.0
          do 605 j=jsta,jend1
          sumtur = sumtur + tursav(j,k,idim1,l)
605       continue
          sumtur = sumtur/float(jend1-jsta+1)
          do 604 j=jsta,jend1
            ti0(j,k,l,3) = sumtur
            ti0(j,k,l,4) = (sumtur-tursav(j,k,idim1,l))*2.0
604       continue
        end if
        end if
c
      end if
c
c     wraparound in k-direction
c
      if(kwrap.gt.0) then
        do 63 j=jsta,jend1
        do 63 l=1,5
        sumq = 0.0
        do 64 k=ksta,kend1
        sumq = sumq+q(j,k,idim1,l)
64      continue
        sumq = sumq/float(kend1-ksta+1)
        do 63 k=ksta,kend1
        qi0(j,k,l,3) = sumq
        qi0(j,k,l,4) = (sumq-q(j,k,idim1,l))*2.0
        bci(j,k,2)   = 1.0
63      continue
c
65      continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 606 j=jsta,jend1
          sumvis = 0.0
          do 607 k=ksta,kend1
            sumvis = sumvis+vist3d(j,k,idim1)
607       continue
          sumvis = sumvis/float(jend1-jsta+1)
          do 606 k=ksta,kend1
            vi0(j,k,1,3) = sumvis
            vi0(j,k,1,4) = (sumvis-vist3d(j,k,idim1))*2.0
606     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 608 j=jsta,jend1
          do 608 l=1,nummem
          sumtur = 0.0
          do 609 k=ksta,kend1
          sumtur = sumtur + tursav(j,k,idim1,l)
609       continue
          sumtur = sumtur/float(jend1-jsta+1)
          do 608 k=ksta,kend1
            ti0(j,k,l,3) = sumtur
            ti0(j,k,l,4) = (sumtur-tursav(j,k,idim1,l))*2.0
608       continue
        end if
        end if
c
      end if
      end if
c
      return
      end
