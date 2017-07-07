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
      subroutine setup(lw,lw2,w,mgwk,itest,jtest,ktest,maxbl,maxgr,
     .                 ngrid,ncgg,iemg,nblg)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Read in the grid calulate some preliminary information,
c     such as 1D array indicies
c***********************************************************************
c
c     maxbl   - maximum number of blocks
c     maxgr   - maximum number of grids
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile
c
      dimension w(mgwk),lw(65,maxbl),lw2(43,maxbl)
      dimension itest(maxgr),jtest(maxgr),ktest(maxgr)
      dimension ncgg(maxgr),iemg(maxgr),nblg(maxgr)
c
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /igrdtyp/ ip3dgrd,ialph
c
      iwr  = 1
c
      write(25,92)
   92 format(/,1x,16hREADING IN GRIDS)
      do 1007 igrid=1,ngrid
      iem  = iemg(igrid)
      nbl  = nblg(igrid)
      call lead(nbl,lw,lw2,maxbl)
c
c      grids
c
      if(ip3dgrd.eq.0) then
c      
c       cfl3d type
c      
        read(1) jdum,kdum,idum
        write(25,93)igrid,idum,jdum,kdum
   93   format(/,13h reading grid,i3,24h of dimensions (I/J/K) :,3i5)
        if (jdum.ne.jdim .or. kdum.ne.kdim .or. idum.ne.idim) then
           write(25,*) ' stopping....inconsistency in grid data file'
           stop
        end if
        read(1) (w(i),i=lx,lvis-1)
c      
      else
c
c       plot3d type
c
        if (igrid.eq.1) then
           read(1) ndum
           if (ndum.ne.ngrid) then
              write(25,*) ' stopping....ngrid = ',ngrid,
     .        ' but grid file contains ',ndum,' grids'
              stop
           end if
           read(1) (itest(ll),jtest(ll),ktest(ll),ll=1,ngrid)
        end if
        write(25,93)igrid,itest(igrid),jtest(igrid),ktest(igrid)
        if (jtest(igrid).ne.jdim .or. ktest(igrid).ne.kdim .or.
     .  itest(igrid).ne.idim) then
          write(25,*) ' stopping....inconsistency in grid data file'
          stop
        end if
c
c       set ialph=1 to be same as in old ronnie
c
        ialph = 1
c
        irr=0
        call rp3d(w(lx),w(ly),w(lz),jdim,kdim,idim,igrid,ialph,0,irr)
        if (irr .ne. 0) then
          write(11,*) ' Stopping... error reading grid...'
          write(11,*) ' (Common error:  grid not written in same',
     .    ' precision that RONNIE was compiled)'
          stop
        end if
c
      end if
c
c      coarser global grids
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
        do 10 m=1,ncg
        nbl   = nbl+1
        lxc   = lw(10,nbl)
        lyc   = lw(11,nbl)
        lzc   = lw(12,nbl)
        call collx(w(lx),w(ly),w(lz),w(lxc),w(lyc),w(lzc),jdim,kdim,
     .  idim,jj2,kk2,ii2)
        call lead(nbl,lw,lw2,maxbl)
   10   continue
      end if
1007  continue
      return
      end
