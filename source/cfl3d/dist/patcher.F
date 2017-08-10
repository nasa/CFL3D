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
      subroutine patcher(nbl,lw,w,mgwk,wk,nwork,ncall,ioutpt,it_thro,
     .                   maxbl,msub1,intmx,mxxe,mptch,jdimg,kdimg,
     .                   idimg,windx,nintr,iindx,llimit,iitmax,
     .                   mmcxie,mmceta,ncheck,iifit,mblkpt,iic0,
     .                   iiorph,iitoss,ifiner,factjlo,factjhi,
     .                   factklo,factkhi,dx,dy,dz,dthetx,dthety,
     .                   dthetz,isav_dpat,isav_dpat_b,
     .                   xte,yte,zte,xmi,ymi,zmi,xmie,ymie,zmie,
     .                   jjmax1,kkmax1,jimage,kimage,xorig,yorig,zorig,
     .                   jte,kte,sxie,seta,sxie2,seta2,xie2s,eta2s,
     .                   temp,x2,y2,z2,nblk1,nblk2,jmm,kmm,x1,y1,z1,
     .                   lout,xif1,xif2,etf1,etf2,ireq_ar,
     .                   myid,myhost,mycomm,mblk2nd,nou,bou,nbuf,
     .                   ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Calculate patched-grid interpolation coefficients 
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
#   ifdef DBLE_PRECSN
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_DOUBLE_COMPLEX
#      else
#        define MY_MPI_REAL MPI_DOUBLE_PRECISION
#      endif
#   else
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_COMPLEX
#      else
#        define MY_MPI_REAL MPI_REAL
#      endif
#   endif
#   ifdef BUILD_MPE
#     include "mpef.h"
#   endif
      dimension istat(MPI_STATUS_SIZE)
#endif
c
c     maxbl - maximum number of blocks
c     maxgr - maximum number of grids
c     intmx - maximum number of block interpolations
c     mptch - maximum dimension of any block face involved in a patch
c     mxxe  - size of 1-d array used to store interpolation coefficients for
c             all patch interfaces, including those on coarser blocks
c     msub1 - maximum number of blocks a single patch face may be
c                interpolated from
c***********************************************************************
c
c      For each interpolation (int):
c   
c      iindx(int,1) - number of sub-faces which make up the "from" side
c                      of the patch surface
c      
c      iindx(int,2)
c           .
c           .       - global block number(s) of the blocks which contain
c           .         the sub-faces on the "from' side of the patch surface
c           .         note: must have iindx(int,1) values
c           .  
c      iindx(int,iindx(int,1)+1)    
c                      
c      iindx(int,iindx(int,1)+2))  - global block number of the block which
c                                    contains the "to" side of the patch 
c                                    surface
c
c      iindx(int,iindx(int,1)+3))
c                  .                
c                  .                - topology of patch in "from" (sub)-block(s)
c                  .                  of the form mn (see below)
c                  .
c      iindx(int,2*iindx(int,1)+2) 
c
c                       
c      iindx(int,2*iindx(int,1)+3)  - topology of the "to" side of the
c                                      patch surface
c
c      iindx(int,2*iindx(int,1)+4)  - number of grid points in cross-section 
c                                      of the "to" grid
c      iindx(int,2*iindx(int,1)+5)  - starting location of generalized 
c                      coordinate data for interpolation to the "to" side
c                      of the patch surface from the "from" side of the 
c                      patch surface
c
c      iindx(int,2*iindx(int,1)+6)  - starting location in xie for the
c                                     "to" side of the patch interface 
c
c      iindx(int,2*iindx(int,1)+7)  - ending location in xie for the
c                                     "to" side of the patch interface
c
c      iindx(int,2*iindx(int,1)+8)  - starting location in eta for the
c                                     "to" side of the patch interface
c
c      iindx(int,2*iindx(int,1)+9)  - ending location in eta for the
c                                     "to" side of the patch interface
c
c
c      iindx(int,2*iindx(int,1)+9+1)
c                  .
c                  .                - starting location for the search
c                  .                  range in xie on the "from" side
c                  .                  of the patch interface
c                  .
c      iindx(int,3*iindx(int,1)+9)
c
c
c      iindx(int,3*iindx(int,1)+9+1)
c                  .
c                  .                - ending location for the search
c                  .                  range in xie on the "from" side
c                  .                  of the patch interface
c                  .
c      iindx(int,4*iindx(int,1)+9)
c
c
c      iindx(int,4*iindx(int,1)+9+1)
c                  .
c                  .                - starting location for the search
c                  .                  range in eta on the "from" side
c                  .                  of the patch interface
c                  .
c      iindx(int,5*iindx(int,1)+9)
c
c
c      iindx(int,5*iindx(int,1)+9+1)
c                  .
c                  .               - ending location for the search
c                  .                 range in eta on the "from" side
c                  .                 of the patch interface
c                  .
c      iindx(int,6*iindx(int,1)+9)
c
c***********************************************************************
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,avgg,avgq
      character*120 bou(ibufdim,nbuf)
c
      integer xi1,xi2,et1,et2,xif1l,xif2l,etf1l,etf2l,
     .        xi1f,xi2f,et1f,et2f
      integer lout(msub1),xif1(msub1),xif2(msub1),etf1(msub1),
     .        etf2(msub1)
c
      dimension nou(nbuf)
      dimension jjmax2(1),kkmax2(1)
      dimension w(mgwk),wk(nwork),lw(65,maxbl),mblk2nd(maxbl)
      dimension xorig(maxbl),yorig(maxbl),zorig(maxbl)
      dimension xte(mptch+2,mptch+2,msub1),
     .          yte(mptch+2,mptch+2,msub1),
     .          zte(mptch+2,mptch+2,msub1)
      dimension xmi(mptch+2,mptch+2,msub1),
     .          ymi(mptch+2,mptch+2,msub1),
     .          zmi(mptch+2,mptch+2,msub1)
      dimension xmie(mptch+2,mptch+2,msub1),
     .          ymie(mptch+2,mptch+2,msub1),
     .          zmie(mptch+2,mptch+2,msub1)
      dimension jjmax1(msub1),kkmax1(msub1)
      dimension jimage(msub1,mptch+2,mptch+2),
     .          kimage(msub1,mptch+2,mptch+2)
      dimension jte(msub1),kte(msub1)
      dimension sxie(mptch+2,mptch+2,msub1),
     .          seta(mptch+2,mptch+2,msub1),
     .          sxie2(mptch+2,mptch+2),
     .          seta2(mptch+2,mptch+2)
      dimension xie2s(mptch+2,mptch+2),eta2s(mptch+2,mptch+2)
      dimension temp((mptch+2)*(mptch+2))
      dimension x2(mptch+2,mptch+2),y2(mptch+2,mptch+2),
     .          z2(mptch+2,mptch+2)
      dimension nblk1(mptch+2),nblk2(mptch+2)
      dimension jmm(mptch+2),kmm(mptch+2)
      dimension jdimg(maxbl),kdimg(maxbl),idimg(maxbl)
      dimension x1(mptch+2,mptch+2),y1(mptch+2,mptch+2),
     .          z1(mptch+2,mptch+2)
      dimension windx(mxxe,2),iindx(intmx,6*msub1+9),
     .          llimit(intmx),iitmax(intmx),mmcxie(intmx),
     .          mmceta(intmx),ncheck(maxbl),iifit(intmx),
     .          mblkpt(mxxe),iic0(intmx),iiorph(intmx),iitoss(intmx),
     .          ifiner(intmx)
      dimension factjlo(intmx,msub1),factjhi(intmx,msub1),
     .          factklo(intmx,msub1),factkhi(intmx,msub1)
      dimension dx(intmx,msub1),dy(intmx,msub1),dz(intmx,msub1),
     .          dthetx(intmx,msub1),dthety(intmx,msub1),
     .          dthetz(intmx,msub1)
      dimension isav_dpat(intmx,17),isav_dpat_b(intmx,msub1,6)
      dimension ireq_ar(intmx*3)
c
      common /tol/ epsc,epsc0,epsreen,epscoll
      common /sklt1/ isklt1
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
      common /zero/ iexp
      common /save/ locv,ireq
      common /is_dpatch/ maxdcnt
c
c******* set iavg flag for coarser level data ********
c
c     iavg = 0, use finer-level averages only if search on
c               coarser level fails
c            1, always use finer-level averages for coarser
c               level interpolation data
c
      iavg = 0
c
c********** set tolerances, etc. **********************
c
c     (10.**(-iexp) is machine zero)
c
c     expansion factors for "from" grids
      factj  = 0.01
      factk  = 0.01
c     factj  = 0.002
c     factk  = 0.002
c
c     convergence tolerance for generalized coordinates
      epsc = max(1.e-07,10.**(-iexp+1))
c
c     threshold for collapsed boundaries
      epscoll = max(1.0e-10,10.**(-iexp+1))
c
c     threshold for reentrant boundaries
      epsreen = max(1.0e-09,10.**(-iexp+1))
c
c     threshold for C-0 continuous interfaces
      epsc0   = max(1.0e-07,10.**(-iexp+1))
c
c******************************************************
c
      isklt1 = 0
      if (ioutpt .gt. 0) isklt1 = 1
      if (it_thro.gt.1) go to 909
c
#if defined DIST_MPI
c
c     set baseline tag values
c
      itag_ptch = 1
c
c     post a bunch of receives first (for non-buffering implementations)
c     set the request index and index for wk
c
      locv = 1
      ireq = 1
      itag = 0
c
      do lcnt = 1,maxdcnt
c        nbll is the current (to) block
c        mbl is the source (from) block
         nbll    = isav_dpat(lcnt,1)
         nd_dest = mblk2nd(nbll) 
         lmax1   = isav_dpat(lcnt,2)
         do ll = 1, lmax1
            itag = itag + 1
            mbl = isav_dpat_b(lcnt,ll,1)
            nd_srce = mblk2nd(mbl)
            if (nd_srce.ne.myid) then
               if (nd_dest.eq.myid) then
                  mdim1  = isav_dpat_b(lcnt,ll,3)
                  ndim1  = isav_dpat_b(lcnt,ll,4)
                  jkdim  = 3*mdim1*ndim1 
                  lcheck = locv + jkdim
                  if (lcheck.gt.nwork) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*) ' stopping in patcher....',
     .                              ' work array insufficient'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  mytag = itag_ptch + itag
                  call MPI_IRecv (wk(locv), jkdim,
     .                            MY_MPI_REAL,
     .                            nd_srce,mytag,mycomm,
     .                            ireq_ar(ireq),ierr)
                  locv = lcheck
                  ireq   = ireq + 1
               endif
            endif
         end do
      end do
c
      if (ireq.gt.intmx*3) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),999) ireq,intmx*3
 999     format(' problem in patcher...ireq = ',i4,
     .   ' but max allowable value = intmx*3 = ',i4)
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      endif
c
c     loop over all blocks looking for blocks that need to send out
c     info to other processors
c
      itag = 0
c
      do lcnt = 1,maxdcnt
c        nbll is the current (to) block
c        mbl is the source (from) block
         lmax1 = isav_dpat(lcnt,2)
         do ll = 1, lmax1
            itag = itag + 1
            mbl  = isav_dpat_b(lcnt,ll,1)
            if (mblk2nd(mbl).eq.myid) then
               nbll    = isav_dpat(lcnt,1)
               nd_dest = mblk2nd(nbll)
               if (nd_dest.ne.myid) then
                  lx1    = lw(10,mbl)
                  ly1    = lw(11,mbl)
                  lz1    = lw(12,mbl)
                  icheck = isav_dpat(lcnt,17)
                  jindex = iindx(icheck,ll+lmax1+2)
                  jdim1  = jdimg(mbl)
                  kdim1  = kdimg(mbl)
                  idim1  = idimg(mbl)
                  mdim1  = isav_dpat_b(lcnt,ll,3)
                  ndim1  = isav_dpat_b(lcnt,ll,4)
                  jkdim  = mdim1*ndim1
                  lcheck = locv + 3*jkdim
                  if (lcheck.gt.nwork) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),*)' stopping in patcher....',
     .                                ' work array insufficient'
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  call loadgr(w,mgwk,lx1,ly1,lz1,jindex,
     .                        wk(locv),wk(locv+jkdim),
     .                        wk(locv+2*jkdim),
     .                        mdim1,ndim1,idim1,jdim1,kdim1)
                  nvals = 3*jkdim
                  mytag = itag_ptch + itag
                  call MPI_Send(wk(locv), nvals,
     .                          MY_MPI_REAL,
     .                          nd_dest, mytag, mycomm, ierr)
                  locv = lcheck
               endif
            endif
         end do
      end do
c
      ireq = 1
      locv = 1
#endif
c
 909     continue
c
      do lcnt = 1,maxdcnt
c        nbll is the current (to) block
c        mbl is the source (from) block
c        nbl is block coming into this routine
         nbll    = isav_dpat(lcnt,1)
         if (nbll.eq.nbl) then
            nd_dest = mblk2nd(nbl)
            lmax1   = isav_dpat(lcnt,2)
            icheck  = isav_dpat(lcnt,17)
            lst    = iindx(icheck,2*lmax1+5)
            npt    = iindx(icheck,2*lmax1+4)
            ifit   = iifit(icheck)
            limit0 = llimit(icheck)
            itmax  = iitmax(icheck)
            mcxie  = mmcxie(icheck)
            mceta  = mmceta(icheck)
            ic0    = iic0(icheck)
            iorph  = iiorph(icheck)
            itoss0 = iitoss(icheck)
            xi1    = iindx(icheck,2*lmax1+6)
            xi2    = iindx(icheck,2*lmax1+7)
            et1    = iindx(icheck,2*lmax1+8)
            et2    = iindx(icheck,2*lmax1+9)
            if (mcxie.gt.100) mcxie = 200
            if (mceta.gt.100) mceta = 200
c
            if (isklt1.gt.0) then
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),*)
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),*)
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),721) icheck
721            format(1x,'generalized coordinate interpolation ',
     .                   'number ',i3)
            end if
c
            do 1605 l=1,lmax1
            mbl = iindx(icheck,l+1)
            xif1l = iindx(icheck,2*lmax1+9+l)
            xif2l = iindx(icheck,3*lmax1+9+l)
            etf1l = iindx(icheck,4*lmax1+9+l)
            etf2l = iindx(icheck,5*lmax1+9+l)
c            
c           patch surface to be interpolated from is an i=constant surface
c
            if (iindx(icheck,l+lmax1+2)/10.eq.1) then
               jjmax1(l) = jdimg(mbl)
               kkmax1(l) = kdimg(mbl)
               if (isklt1.gt.0) then
                  if (iindx(icheck,l+lmax1+2).eq.11) then
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1957) 1,mbl,xif1l,xif2l,
     .                                      etf1l,etf2l
                  else
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1957) idimg(mbl),mbl,xif1l,
     .                                      xif2l,etf1l,etf2l
                  end if
               end if
            end if
c      
c           patch surface to be interpolated from is a j=constant surface
c
            if (iindx(icheck,l+lmax1+2)/10.eq.2) then
               jjmax1(l) = kdimg(mbl)
               kkmax1(l) = idimg(mbl)
               if (isklt1.gt.0) then
                  if (iindx(icheck,l+lmax1+2).eq.21) then
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1958) 1,mbl,xif1l,xif2l,
     .                                      etf1l,etf2l
                  else 
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1958) jdimg(mbl),mbl,xif1l,
     .                                      xif2l,etf1l,etf2l
                 end if
               end if
            end if      
c           
c           patch surface to be interpolated from is a k=constant surface
c
            if (iindx(icheck,l+lmax1+2)/10.eq.3) then
               jjmax1(l) = jdimg(mbl)
               kkmax1(l) = idimg(mbl)
               if (isklt1.gt.0) then
                  if (iindx(icheck,l+lmax1+2).eq.31) then
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1959) 1,mbl,xif1l,xif2l,
     .                                      etf1l,etf2l
                  else 
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1959) kdimg(mbl),mbl,xif1l,
     .                                      xif2l,etf1l,etf2l
                  end if
               end if
            end if
c
            if (jjmax1(l).gt.mptch .or. kkmax1(l).gt.mptch) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'('' program terminated in '',
     .         ''dynamic patching routines - see file '',a60)') grdmov
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*) 'stopping ... mptch =',mptch,
     .         ' too small '
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if 
c
 1605       continue
c
c           patch surface to be interpolated to is an i=constant surface
c
            if (iindx(icheck,2*lmax1+3)/10.eq.1) then
               jjmax2(1) = jdimg(nbl)
               kkmax2(1) = kdimg(nbl)
               if (isklt1.gt.0) then
                  if (iindx(icheck,2*lmax1+3).eq.11) then
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1857) 1,nbl,xi1,xi2,et1,et2
                  else 
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1857) idimg(nbl),nbl,xi1,xi2,
     .                                      et1,et2
                  end if
               end if
            end if
c            
c           patch surface to be interpolated to is a j=constant surface
c
            if (iindx(icheck,2*lmax1+3)/10.eq.2) then
               jjmax2(1) = kdimg(nbl)
               kkmax2(1) = idimg(nbl)
               if (isklt1.gt.0) then
                  if (iindx(icheck,2*lmax1+3).eq.21) then
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1858) 1,nbl,xi1,xi2,et1,et2 
                  else 
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1858) jdimg(nbl),nbl,xi1,xi2,
     .                                      et1,et2
                  end if
               end if     
            end if
c            
c           patch surface to be interpolated to is a k=constant surface
c
            if (iindx(icheck,2*lmax1+3)/10.eq.3) then
               jjmax2(1) = jdimg(nbl)
               kkmax2(1) = idimg(nbl)
               if (isklt1.gt.0) then
                  if (iindx(icheck,2*lmax1+3).eq.31) then
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1859) 1,nbl,xi1,xi2,et1,et2
                  else 
                     nou(4) = min(nou(4)+1,ibufdim)
                     write(bou(nou(4),4),1859) kdimg(nbl),nbl,xi1,xi2,
     .                                      et1,et2
                  end if
               end if
            end if
c
            if (jjmax2(1).gt.mptch .or. kkmax2(1).gt.mptch) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'('' program terminated in '',
     .         ''dynamic patching routines - see file '',a60)') grdmov
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*) 'stopping ... mptch =',mptch,
     .         ' too small '
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if 
c
 1857       format(' ',16x,'to surface i =',i3,' of block',i3,2x,'(j: ',
     .             i3,' to ',i3,', k: ',i3,' to ',i3,')')
 1858       format(' ',16x,'to surface j =',i3,' of block',i3,2x,'(k: ',
     .             i3,' to ',i3,', i: ',i3,' to ',i3,')')
 1859       format(' ',16x,'to surface k =',i3,' of block',i3,2x,'(j: ',
     .             i3,' to ',i3,', i: ',i3,' to ',i3,')')
c
 1957       format(' ','interpolation from surface i =',i3,' of block',
     .             i3,2x,'(j: ',i3,' to ',i3,', k: ',i3,' to ',i3,')')
 1958       format(' ','interpolation from surface j =',i3,' of block',
     .             i3,2x,'(k: ',i3,' to ',i3,', i: ',i3,' to ',i3,')')
 1959       format(' ','interpolation from surface k =',i3,' of block',
     .             i3,2x,'(j: ',i3,' to ',i3,', i: ',i3,' to ',i3,')')
c
            if (isklt1.gt.0) then
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),*)
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),99)
            end if
   99       format(' ','note: j and k referred to below are indicies', 
     .             ' local to the patch surface')
c
c           set-up for "to" side of patch interface
c
            jmax2  = jjmax2(1)
            kmax2  = kkmax2(1) 
            lx2    = lw(10,nbl)
            ly2    = lw(11,nbl)
            lz2    = lw(12,nbl)
            idim2g = idimg(nbl)
            jdim2g = jdimg(nbl)
            kdim2g = kdimg(nbl)
            jindex = iindx(icheck,2*lmax1+3)
c
c           load proper grid from 1-d array into 2-d work array.
c
            call loadgr(w,mgwk,lx2,ly2,lz2,jindex,x2,y2,z2,mptch+2,
     .                  mptch+2,idim2g,jdim2g,kdim2g)
c
c           check for collapsed grid lines
c
            if (isklt1.gt.0) then
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),3958) nbl
            end if
 3958       format(' ','   checking for collapsed grid lines on ',
     .      ' "to"  side in block ',i3)
c
            call collapse(mptch+2,mptch+2,jmax2,kmax2,x2,y2,z2,
     .                    nou,bou,nbuf,ibufdim)
c
c           set-up for "from" side of patch interface
c
            iself  = 0
            do 125 l=1,lmax1
            mbl    = iindx(icheck,l+1)
            lx1    = lw(10,mbl)
            ly1    = lw(11,mbl)
            lz1    = lw(12,mbl)
            jmax1  = jjmax1(l)
            kmax1  = kkmax1(l)
            jindex = iindx(icheck,l+lmax1+2)
            idim1  = idimg(mbl)
            jdim1  = jdimg(mbl)
            kdim1  = kdimg(mbl)
c           iself = 1 if a block face communicates with itself
            itest1 = 100*nbl+iindx(icheck,2*lmax1+3)
            itest2 = 100*mbl+iindx(icheck,l+lmax1+2)
            if (itest1 .eq. itest2) iself = 1
c
c           load proper grid from 1-d array into 2-d work array
c
#if       defined DIST_MPI
            if (myid.ne.mblk2nd(mbl)) then
               mdim1 = jmax1
               ndim1 = kmax1
               jkdim = 3*mdim1*ndim1
               call MPI_Wait (ireq_ar(ireq), istat, ierr)
               call mreal(wk(locv),mdim1,ndim1,mptch+2,mptch+2,
     .                    x1,y1,z1)
               locv = locv + jkdim
               ireq = ireq + 1
            else
               call loadgr(w,mgwk,lx1,ly1,lz1,jindex,x1,y1,z1,mptch+2,
     .                     mptch+2,idim1,jdim1,kdim1)
            end if
#else
            call loadgr(w,mgwk,lx1,ly1,lz1,jindex,x1,y1,z1,mptch+2,
     .                  mptch+2,idim1,jdim1,kdim1)
#endif
c
c           check for collapsed grid lines
c
            if (isklt1.gt.0) then
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),3957) mbl
            end if
 3957       format(' ','   checking for collapsed grid lines on ',
     .                    '"from" side in block ',i3)
c
            call collapse(mptch+2,mptch+2,jmax1,kmax1,x1,y1,z1,
     .                    nou,bou,nbuf,ibufdim)
c
c           check for branch cuts 
c
            if (isklt1.gt.0) then
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),98) mbl
            end if
   98       format(' ','   checking for branch cuts',
     .             ' on "from" side in block ',i3) 
c
            xif1l = iindx(icheck,2*lmax1+9+l)
            xif2l = iindx(icheck,3*lmax1+9+l)
            etf1l = iindx(icheck,4*lmax1+9+l)
            etf2l = iindx(icheck,5*lmax1+9+l)
            call rechk(mptch+2,mptch+2,jimage,kimage,msub1,
     .                 jmax1,kmax1,l,x1,y1,z1,xif1l,xif2l,etf1l,etf2l,
     .                 nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl)
c
c           expand "from" grid(s) at boundaries to insure that the 
c           "to" grid is completely covered
c
            if (isklt1.gt.0) then
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),1755) mbl
            end if
 1755       format(' ','   expanding grid boundaries on "from" side',
     .             ' in block ',i3)
c
c           for now, set lo/hi values for factj/factk to single
c           factj/factk values set above...later may add the
c           ability to read individual values on input
c
            factjlo(icheck,l) = factj
            factjhi(icheck,l) = factj
            factklo(icheck,l) = factk
            factkhi(icheck,l) = factk
c
c           should not really need any significant expansion
c           for c-0 grids
c
            if (ic0.gt.0) then
               factjlo(icheck,l) = 1.e-5
               factjhi(icheck,l) = 1.e-5
               factklo(icheck,l) = 1.e-5
               factkhi(icheck,l) = 1.e-5
            end if
c
            call expand(mptch+2,mptch+2,msub1,jmax1,kmax1,l,
     .                  x1,y1,z1,xte,yte,zte,factjlo(icheck,l),
     .                  factjhi(icheck,l),factklo(icheck,l),
     .                  factkhi(icheck,l),jmax2,kmax2,x2,y2,z2)
            jte(l) = jjmax1(l)+2
            kte(l) = kkmax1(l)+2
c
c           translate/rotate "from" blocks as needed to provide
c           sufficient coverage on "from" side when only part
c           of the physical domain is modeled
c
            if (abs(real(dx(icheck,l))).gt.0. .or.
     .          abs(real(dy(icheck,l))).gt.0. .or.
     .          abs(real(dz(icheck,l))).gt.0.) then
                jjte = jte(l)
                kkte = kte(l)
                call transp(mptch+2,mptch+2,jjte,kkte,msub1,l,
     .                      xte,yte,zte,dx,dy,dz,intmx,icheck)
            end if
            if (abs(real(dthetx(icheck,l))).gt.0. .or. 
     .          abs(real(dthety(icheck,l))).gt.0. .or.
     .          abs(real(dthetz(icheck,l))).gt.0.) then
                jjte = jte(l)
                kkte = kte(l)
                call rotatp(mptch+2,mptch+2,jjte,kkte,msub1,l,
     .                      xte,yte,zte,dthetx,dthety,dthetz,
     .                      xorig,yorig,zorig,mbl,maxbl,intmx,icheck)
            end if
c
c           search range on "from" side
c
            xif1(l) = iindx(icheck,2*lmax1+9+l)
            if (xif1(l).gt.1) then
               xif1(l) = xif1(l) + 1
            end if
            xif2(l) = iindx(icheck,3*lmax1+9+l)
            if (xif2(l) .eq. jjmax1(l)) then
               xif2(l) = xif2(l) + 2
            else
               xif2(l) = xif2(l) + 1
            end if
            etf1(l) = iindx(icheck,4*lmax1+9+l)
            if (etf1(l).gt.1) then
               etf1(l) = etf1(l) + 1
            end if
            etf2(l) = iindx(icheck,5*lmax1+9+l) + 1
            if (etf2(l) .eq. kkmax1(l)) then
               etf2(l) = etf2(l) + 2
            else
               etf2(l) = etf2(l) + 1
            end if
c
  125       continue
c
            if (ifiner(icheck).ne.0) then
               icheckf = ifiner(icheck)
               lmax1f  = iindx(icheckf,1)
               lstf    = iindx(icheckf,2*lmax1f+5)
               nptf    = iindx(icheckf,2*lmax1f+4)
               xi1f    = iindx(icheckf,2*lmax1f+6)
               xi2f    = iindx(icheckf,2*lmax1f+7)
               et1f    = iindx(icheckf,2*lmax1f+8)
               et2f    = iindx(icheckf,2*lmax1f+9)
            else
               lstf    = iindx(icheck,2*lmax1+5)
               nptf    = iindx(icheck,2*lmax1+4)
               xi1f    = iindx(icheck,2*lmax1+6)
               xi2f    = iindx(icheck,2*lmax1+7)
               et1f    = iindx(icheck,2*lmax1+8)
               et2f    = iindx(icheck,2*lmax1+9)
            end if
c
c           obtain interpolation coefficients for coarser levels
c           by averaging finer-level coefficients
c
            if (ifiner(icheck).ne.0 .and. iavg .gt.0) then
               call avgint(windx(lst,1),windx(lst,2),mblkpt(lst),npt,
     .                     windx(lstf,1),windx(lstf,2),mblkpt(lstf),
     .                     nptf,xi1,xi2,et1,et2,xi1f,xi2f,et1f,et2f)
            end if
c
c           determine projection of x2,y2,z2 points onto generalized 
c           coordinate system(s) defined by the grid(s) on "from" side
c
            call invert(mptch+2,mptch+2,msub1,1,jte,kte,
     .                  lmax1,xte,yte,zte,xmi,ymi,zmi,xmie,ymie,zmie,
     .                  limit0,jjmax2,kkmax2,x2,y2,z2,windx(lst,1),
     .                  windx(lst,2),mblkpt(lst),temp,jimage,kimage,
     .                  ifit,itmax,sxie,seta,sxie2,seta2,xie2s,eta2s,
     .                  intmx,icheck,nblk1,nblk2,jmm,kmm,mcxie,mceta,
     .                  lout,xi1,xi2,et1,et2,npt,ic0,iorph,itoss0,ncall,
     .                  ioutpt,xif1,xif2,etf1,etf2,iself,ifiner(icheck),
     .                  windx(lstf,1),windx(lstf,2),mblkpt(lstf),nptf,
     .                  xi1f,xi2f,et1f,et2f,iavg,nou,
     .                  bou,nbuf,ibufdim,myid,mblk2nd,maxbl)
c
         end if
      end do
c
      return
      end
