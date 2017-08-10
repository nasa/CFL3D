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
      subroutine updateg(lw,lw2,w,mgwk,wk,nwork,iupdat,iseqr,maxbl,
     .                   maxgr,maxseg,nbci0,nbcj0,nbck0,nbcidim,
     .                   nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                   nblock,levelg,igridg,utrans,vtrans,wtrans,
     .                   omegax,omegay,omegaz,xorig,yorig,zorig,
     .                   thetax,thetay,thetaz,rfreqt,rfreqr,xorig0,
     .                   yorig0,zorig0,time2,thetaxl,thetayl,thetazl,
     .                   itrans,irotat,idefrm,ncgg,iadvance,
     .                   nou,bou,nbuf,ibufdim,myid,myhost,mycomm,
     .                   mblk2nd,irigb,irbtrim,nt)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Update grid to new position, and obtain corresponding
c     grid-boundary velocities for use in the boundary conditions. Also
c     collocate new grid position to coarser levels and obtain grid-
c     boundary velocities on coarser levels, and update moment center
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
      dimension trnsfr(10)
#endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension w(mgwk),lw(65,maxbl),lw2(43,maxbl),wk(nwork)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),
     .          nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          ibcinfo(maxbl,maxseg,7,2),jbcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2)
      dimension levelg(maxbl),igridg(maxbl)
      dimension utrans(maxbl),vtrans(maxbl),wtrans(maxbl),
     .          omegax(maxbl),omegay(maxbl),omegaz(maxbl),
     .          xorig(maxbl),yorig(maxbl),zorig(maxbl),
     .          thetax(maxbl),thetay(maxbl),thetaz(maxbl),
     .          rfreqt(maxbl),rfreqr(maxbl),xorig0(maxbl),
     .          yorig0(maxbl),zorig0(maxbl),time2(maxbl),
     .          thetaxl(maxbl),thetayl(maxbl),thetazl(maxbl),
     .          itrans(maxbl),irotat(maxbl),idefrm(maxbl)
      dimension ncgg(maxgr),iadvance(maxbl),mblk2nd(maxbl)
c
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /sklton/ isklton
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /motionmc/ xmc0,ymc0,zmc0,utransmc,vtransmc,wtransmc,
     .                  omegaxmc,omegaymc,omegazmc,xorigmc,yorigmc,
     .                  zorigmc,xorig0mc,yorig0mc,zorig0mc,thetaxmc,
     .                  thetaymc,thetazmc,dxmxmc,dymxmc,dzmxmc,
     .                  dthxmxmc,dthymxmc,dthzmxmc,rfreqtmc,
     .                  rfreqrmc,itransmc,irotatmc,time2mc
      common /trim/ dmtrmn,dmtrmnm,dlcln,dlclnm,trtol,cmy,cnw,alf0,
     .              alf1,dzdt,thtd0,thtd1,zrg0,zrg1,dtrmsmx,dtrmsmn,
     .              dalfmx,ddtmx,ddtrm0,ddtrm1,itrmt,itrminc,fp(4,4),
     .              tp(4,4),zlfct,epstr,relax,ittrst 
c
      if (isklton.eq.1) then
#if defined DIST_MPI
         if (myid.eq.1) then
#endif
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),400)
  400    format(1x,1h )
#if defined DIST_MPI
         end if
#endif
      end if
c
      icnt  = 0
      do 100 nbl = 1,nblock
      if (myid.eq.mblk2nd(nbl) .and. iadvance(nbl).ge.0 .and.
     .   (levelg(nbl).ge.lglobal .and. 
     .    levelg(nbl).le.levelt(iseqr))) then
c
         call lead(nbl,lw,lw2,maxbl)
         iuns1 = max(irotat(nbl),itrans(nbl))
         iuns  = max(iuns1,irigb,irbtrim)
c
         if (iuns .gt. 0.) then
c
c           zero out temp array used to store incremental grid point
c           velocities and boundary accelerations due to translation 
c           and rotation components of grid motion. NOTE: between this
c           point and the call to subroutine tmetric, the first lt1wk-1
c           elements of the wk array must only be used to  store
c           increments to the grid point velocities/accelerations.
c           (temporarily treat this section of the temporary array 
c           as permanent!)
c
c           temporary storage locations:
c           lvel  = start of grid point velocity array
c           lacci = start of i-boundary point acceleration array
c           laccj = start of j-boundary point acceleration array
c           lacck = start of k-boundary point acceleration array
c           lt1wk = start of work array for subroutines rotate/metric
c           lt2wk = start of work array for subroutine metric
c           lt3wk = start of work array for subroutine metric
c
            lvel  = 1
            lacci = jdim*kdim*idim*3+lvel
            laccj = jdim*kdim*3*2+lacci
            lacck = kdim*idim*3*2+laccj
            lt1wk = jdim*idim*3*2+lacck
            lt2wk = jdim*kdim*idim*3+lt1wk
            lt3wk = jdim*kdim*6+lt2wk
c
c           zero out velocity/acceleration work arrays
c
            mdim  = lt1wk - 1
            do 55 izz=1,mdim
             wk(lvel+izz-1) = 0.e0
   55       continue
c
c           translation/rotation corresponding to rigid-body modes
c           (currently restricted to z-translation and y-rotation only)
c
            if (irigb .gt. 0 .and. nt .gt. 0) then
c
               if (isklton.eq.1) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),605) nbl
  605             format(1x,20htranslating rb-block,i4,
     .                   16h to new position)
               end if
c 
               if (nwork.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),410)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
c
               xold = xorig(nbl)
               yold = yorig(nbl)
               zold = zorig(nbl)
c
c     This formulation assumes small pitch angles.  sin(thetay)=0,
c     cos(thetay) = 1.
c
               zorig(nbl) = zrg1
               utrans(nbl)= 0.
               vtrans(nbl)= 0.
               wtrans(nbl)= (zorig(nbl) - zold)/dt
               call trans(jdim,kdim,idim,wk(lvel),wk(lacci),wk(laccj),
     .                    wk(lacck),w(lx),w(ly),w(lz),99,
     .                    rfreqt(nbl),utrans(nbl),vtrans(nbl),
     .                    wtrans(nbl),xorig(nbl),yorig(nbl),
     .                    zorig(nbl),xold,yold,zold,xorig0(nbl),
     .                    yorig0(nbl),zorig0(nbl),iupdat,time2(nbl))
c
c              rotation corresponding to rigid-body modes
c              (currently restricted to y-rotation only)
c
               if(isklton.eq.1)then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),615) nbl
               end if
  615          format(1x,19hrotating rb-block   ,i4
     .                                ,16h to new position)
c
               nroom = nwork-lt1wk
               mdim  = jdim*kdim*idim*3
               if (nroom.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),420)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
c
               call rotate(jdim,kdim,idim,wk(lvel),wk(lacci),wk(laccj),
     .                     wk(lacck),wk(lt1wk),w(lx),w(ly),w(lz),nbl,
     .                     99,rfreqr(nbl),omegax(nbl),
     .                     omegay(nbl),omegaz(nbl),xorig(nbl),
     .                     yorig(nbl),zorig(nbl),thetax(nbl),
     .                     thetay(nbl),thetaz(nbl),thetaxl(nbl),
     .                     thetayl(nbl),thetazl(nbl),iupdat,time2(nbl),
     .                     nou,bou,nbuf,ibufdim)
            end if
c
c           trim step for rigid-body modes
c           (trim currently restricted to y-rotation only)
c
            if (irbtrim .gt. 0 .and. nt .gt. 0) then
c
               if(isklton.eq.1)then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),625) nbl
               end if
  625          format(1x,19hrotating rb-block  ,i4,8h to trim)
c
               nroom = nwork-lt1wk
               mdim  = jdim*kdim*idim*3
               if (nroom.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),420)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
c
               omegay(nbl)  = epstr*(thetay(nbl)-thetayl(nbl))
               call rotate(jdim,kdim,idim,wk(lvel),wk(lacci),wk(laccj),
     .                     wk(lacck),wk(lt1wk),w(lx),w(ly),w(lz),nbl,
     .                     99,rfreqr(nbl),omegax(nbl),
     .                     omegay(nbl),omegaz(nbl),xorig(nbl),
     .                     yorig(nbl),zorig(nbl),thetax(nbl),
     .                     thetay(nbl),thetaz(nbl),thetaxl(nbl),
     .                     thetayl(nbl),thetazl(nbl),iupdat,time2(nbl),
     .                     nou,bou,nbuf,ibufdim)
            end if
c
c           translation 
c
            if (itrans(nbl) .gt. 0 .and. irbtrim .eq. 0 .and.
     .          irigb .eq. 0 .and. nt .gt. 0) then
c
               if (isklton.eq.1) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),405) nbl
  405             format(1x,17htranslating block,i4,
     .                   16h to new position)
               end if
c
               if (nwork.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),410)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  410          format(43h not enough work space for subroutine trans)
c
               xold = xorig(nbl)
               yold = yorig(nbl)
               zold = zorig(nbl)
               call trans(jdim,kdim,idim,wk(lvel),wk(lacci),wk(laccj),
     .                    wk(lacck),w(lx),w(ly),w(lz),itrans(nbl),
     .                    rfreqt(nbl),utrans(nbl),vtrans(nbl),
     .                    wtrans(nbl),xorig(nbl),yorig(nbl),
     .                    zorig(nbl),xold,yold,zold,xorig0(nbl),
     .                    yorig0(nbl),zorig0(nbl),iupdat,time2(nbl))
            end if
c
c           rotation 
c
            if (irotat(nbl) .gt. 0 .and. irbtrim .eq. 0 .and.
     .          irigb .eq. 0 .and. nt .gt. 0) then
c
               if(isklton.eq.1)then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),415) nbl
               end if
  415          format(1x,17hrotating block   ,i4,16h to new position)
c
               nroom = nwork-lt1wk
               mdim  = jdim*kdim*idim*3
               if (nroom.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),420)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  420          format(44h not enough work space for subroutine rotate)
c
               thetaxl(nbl) = thetax(nbl)
               thetayl(nbl) = thetay(nbl)
               thetazl(nbl) = thetaz(nbl)
               call rotate(jdim,kdim,idim,wk(lvel),wk(lacci),wk(laccj),
     .                     wk(lacck),wk(lt1wk),w(lx),w(ly),w(lz),nbl,
     .                     irotat(nbl),rfreqr(nbl),omegax(nbl),
     .                     omegay(nbl),omegaz(nbl),xorig(nbl),
     .                     yorig(nbl),zorig(nbl),thetax(nbl),
     .                     thetay(nbl),thetaz(nbl),thetaxl(nbl),
     .                     thetayl(nbl),thetazl(nbl),iupdat,time2(nbl),
     .                     nou,bou,nbuf,ibufdim)
            end if
c
c           if the current block will also undergo deformation, defer
c           updating metrics (temporal and spatial) until the changes
c           due to deformation are added.

            if (idefrm(nbl) .eq. 0) then
c
c              calculate face-average values of velocity and acceleration
c              on block boundaries and place in permanent storage for use
c              in boundary condition routines
c
               call xtbatb(jdim,kdim,idim,w(lxtbj),w(lxtbk),w(lxtbi),
     .                     w(latbj),w(latbk),w(latbi),wk(lvel),
     .                     wk(lacci),wk(laccj),wk(lacck))
c
               nroom = nwork-lt3wk
               mdim  = jdim*kdim*idim*5
               if (nroom.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),425)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  425          format(45h not enough work space for metric subroutines)
c
c              calculate spatial metrics for updated grid
c
               iflag = -1
               call metric(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),
     .                     w(lsk),w(lsi),wk(lt2wk),wk(lt3wk),nbl,
     .                     iflag,icnt,nbci0,nbcj0,nbck0,nbcidim,
     .                     nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                     maxbl,maxseg,nou,bou,nbuf,ibufdim,myid,
     .                     mblk2nd)

c
c              calculate temporal metrics for updated grid
c
               call tmetric(jdim,kdim,idim,w(lsj),w(lsk),w(lsi),
     .                      w(lx),w(ly),w(lz),wk(lvel),wk(lt1wk),
     .                      wk(lt2wk),wk(lt3wk),nbl)
c
c              coarser levels
c
               ncg = ncgg(igridg(nbl)) - (mseq-iseqr)
               if (ncg.gt.0 .and. mgflag.gt.0) then
                  nbll = nbl
                  do 1820 m=1,ncg
                  nbll = nbll+1
                  nbllm1 = nbll - 1
                  time2(nbll)  = time2(nbl)
                  xorig(nbll)  = xorig(nbl)
                  yorig(nbll)  = yorig(nbl)
                  zorig(nbll)  = zorig(nbl)
                  thetax(nbll) = thetax(nbl)
                  thetay(nbll) = thetay(nbl)
                  thetaz(nbll) = thetaz(nbl)
                  lvolc  = lw( 8,nbll)
                  lxc    = lw(10,nbll)
                  lyc    = lw(11,nbll)
                  lzc    = lw(12,nbll)
                  lxtbjc = lw(36,nbll)
                  lxtbkc = lw(37,nbll)
                  lxtbic = lw(38,nbll)
                  latbjc = lw(39,nbll)
                  latbkc = lw(40,nbll)
                  latbic = lw(41,nbll)
                  lvelc  = lt1wk
c
                  if(isklton.eq.1)then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),850) nbll,ii2,jj2,kk2
                  end if
  850             format(1x,24h  creating coarser block,i4,
     .            24h of dimensions (I/J/K) :,3i4)
c
c                 collocate xyz
                  call collx(w(lx),w(ly),w(lz),w(lxc),w(lyc),w(lzc),
     .                       jdim,kdim,idim,jj2,kk2,ii2)
c
c                 collocate grid point velocity
                  call collxt(wk(lvel),wk(lvelc),jdim,kdim,idim,
     .                        jj2,kk2,ii2,nbllm1,nou,bou,nbuf,ibufdim)
                  nv = jj2*kk2*ii2*3
                  do 1825 izz = 1,nv
                  wk(lvel+izz-1) = wk(lvelc+izz-1)
 1825             continue
c
c                 collocate i0/idim boundary velocity/acceleration
                  call collxtb(w(lxtbi),w(lxtbic),jdim,kdim,
     .                        jj2,kk2,nbllm1)
                  call collxtb(w(latbi),w(latbic),jdim,kdim,
     .                        jj2,kk2,nbllm1)
c
c                 collocate j0/jdim boundary velocity/acceleration 
                  call collxtb(w(lxtbj),w(lxtbjc),kdim,idim-1,
     .                        kk2,ii2-1,nbllm1)
                  call collxtb(w(latbj),w(latbjc),kdim,idim-1,
     .                        kk2,ii2-1,nbllm1)
c
c                 collocate k0/kdim boundary velocity/acceleration
                  call collxtb(w(lxtbk),w(lxtbkc),jdim,idim-1,
     .                        jj2,ii2-1,nbllm1)
                  call collxtb(w(latbk),w(latbkc),jdim,idim-1,
     .                        jj2,ii2-1,nbllm1)
c
c                 calculate spatial metrics for updated coarser grid
c
                  call lead(nbll,lw,lw2,maxbl)
c
                  lvel  = 1
                  lt1wk = jdim*kdim*idim*3+lvel
                  lt2wk = jdim*kdim*idim*3+lt1wk
                  lt3wk = jdim*kdim*6+lt2wk
c
                  iflag = -1
                  call metric(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),
     .                        w(lsk),w(lsi),wk(lt2wk),wk(lt3wk),nbll,
     .                        iflag,icnt,nbci0,nbcj0,nbck0,nbcidim,
     .                        nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                        maxbl,maxseg,nou,bou,nbuf,ibufdim,myid,
     .                        mblk2nd)

c
c                 calculate temporal metrics for updated coarser grid
c
                  call tmetric(jdim,kdim,idim,w(lsj),w(lsk),w(lsi),
     .                         w(lx),w(ly),w(lz),wk(lvel),wk(lt1wk),
     .                         wk(lt2wk),wk(lt3wk),nbll)
c
 1820             continue
c
                  call lead(nbl,lw,lw2,maxbl)
c
               end if
            end if
         end if
c
      end if
c
      if (isklton.eq.1) then
#   ifdef FASTIO
         call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,32)
#   else
         call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
      end if
c
  100 continue
c
c     update moment center location
c
      if (myid.eq.myhost) then
         if (itransmc .gt. 0) then
c
           if (isklton.eq.1) then
               write(11,505)
  505          format(1x,41htranslating moment center to new position)
            end if
c
            call transmc(itransmc,rfreqtmc,utransmc,vtransmc,wtransmc,
     .                   xorigmc,yorigmc,zorigmc,xorig0mc,yorig0mc,
     .                   zorig0mc,xmc,ymc,zmc,iupdat,time2mc)
c
         end if
c
         if (irotatmc .gt. 0) then
c
            if (isklton.eq.1) then
               write(11,506)
  506          format(1x,39h rotating moment center to new position)
            end if
c
            call rotatmc(irotatmc,rfreqrmc,omegaxmc,omegaymc,omegazmc,
     .                   xorigmc,yorigmc,zorigmc,thetaxmc,thetaymc,
     .                   thetazmc,xmc,ymc,zmc,iupdat,time2mc)
c
         end if
      end if
c
#if defined DIST_MPI
c     broadcast new moment center data to all processors
c
      if (myid.eq.myhost) then
         trnsfr(1)  = xmc
         trnsfr(2)  = ymc
         trnsfr(3)  = zmc
         trnsfr(4)  = xorigmc
         trnsfr(5)  = yorigmc
         trnsfr(6)  = zorigmc
         trnsfr(7)  = thetaxmc
         trnsfr(8)  = thetaymc 
         trnsfr(9)  = thetazmc
         trnsfr(10) = time2mc
      end if
c 
      call MPI_Bcast (trnsfr, 10, MY_MPI_REAL, myhost,
     .                mycomm, ierr)
c
      if (myid.ne.myhost) then
         xmc      = trnsfr(1)
         ymc      = trnsfr(2)
         zmc      = trnsfr(3)
         xorigmc  = trnsfr(4)
         yorigmc  = trnsfr(5)
         zorigmc  = trnsfr(6)
         thetaxmc = trnsfr(7)
         thetaymc = trnsfr(8)
         thetazmc = trnsfr(9)
         time2mc  = trnsfr(10)
      end if
#endif
c
      return
      end
