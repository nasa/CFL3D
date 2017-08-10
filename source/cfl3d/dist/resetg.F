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
      subroutine resetg(lw,lw2,w,mgwk,wk,nwork,iupdat,iseqr,maxbl,
     .                  maxgr,maxseg,nbci0,nbcj0,nbck0,nbcidim,
     .                  nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                  nblock,levelg,igridg,utrans,vtrans,wtrans,
     .                  omegax,omegay,omegaz,xorig,yorig,zorig,
     .                  thetax,thetay,thetaz,rfreqt,rfreqr,xorig0,
     .                  yorig0,zorig0,time2,thetaxl,thetayl,thetazl,
     .                  itrans,irotat,idefrm,ncgg,iadvance,
     .                  dxmx,dymx,dzmx,dthymx,dthzmx,dthxmx,
     .                  iitot,iovrlp,lig,lbg,iipntsg,ibpntsg,
     .                  qb,iibg,kkbg,jjbg,ibcg,nou,bou,nbuf,ibufdim,
     .                  myid,myhost,mycomm,mblk2nd,iresetb,nt)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Check to see if any blocks in the grid have been 
c     translated/rotated out of bounds set by input deck. If so,
c     reset these blocks so they are at or near the initial positions.
c     If the rotational displacement of a block is reset, must also
c     rotate the solution to correspond to the reset position.
c
c     Resetting is allowed only for constant translational speed
c     (itrans = 1) or constant rotational speed (irotat = 1)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
      dimension istat(MPI_STATUS_SIZE)
#endif
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension w(mgwk),wk(nwork),lw(65,maxbl),lw2(43,maxbl)
      dimension iresetb(maxbl) 
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
      dimension dxmx(maxbl),dymx(maxbl),dzmx(maxbl),dthxmx(maxbl),
     .          dthymx(maxbl),dthzmx(maxbl)
      dimension iovrlp(maxbl),lig(maxbl),lbg(maxbl),iipntsg(maxbl),
     .          ibpntsg(maxbl,4),qb(iitot,5,3),iibg(iitot),
     .          kkbg(iitot),jjbg(iitot),ibcg(iitot)
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
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /conversion/ radtodeg
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /motionmc/ xmc0,ymc0,zmc0,utransmc,vtransmc,wtransmc,
     .                  omegaxmc,omegaymc,omegazmc,xorigmc,yorigmc,
     .                  zorigmc,xorig0mc,yorig0mc,zorig0mc,thetaxmc,
     .                  thetaymc,thetazmc,dxmxmc,dymxmc,dzmxmc,
     .                  dthxmxmc,dthymxmc,dthzmxmc,rfreqtmc,
     .                  rfreqrmc,itransmc,irotatmc,time2mc
c
      do 5 nbl = 1,nblock
      iresetb(nbl) = 0
      ireset = 0
      if (iadvance(nbl) .lt. 0) go to 5
      iuns = max(irotat(nbl),itrans(nbl))
      if (levelg(nbl).ge.lglobal .and. iuns.gt.0)  then
#if defined DIST_MPI
         if (myid.eq.mblk2nd(nbl)) then
#endif
c
c        check for need to reset
c
         call lead(nbl,lw,lw2,maxbl)
         adxmx   = ccabs(dxmx(nbl))
         adymx   = ccabs(dymx(nbl))
         adzmx   = ccabs(dzmx(nbl))
         adthxmx = ccabs(dthxmx(nbl))
         adthymx = ccabs(dthymx(nbl))
         adthzmx = ccabs(dthzmx(nbl))
         axorig  = ccabs(xorig(nbl) - xorig0(nbl))
         ayorig  = ccabs(yorig(nbl) - yorig0(nbl))
         azorig  = ccabs(zorig(nbl) - zorig0(nbl))
         athetx  = ccabs(thetax(nbl))
         athety  = ccabs(thetay(nbl))
         athetz  = ccabs(thetaz(nbl))
         if (itrans(nbl) .eq. 1) then 
            if (real(adxmx).gt.0. .and. real(axorig) .gt. real(adxmx))
     .         ireset=1
            if (real(adymx).gt.0. .and. real(ayorig).gt.real(adymx))
     .         ireset=1
            if (real(adzmx).gt.0. .and. real(azorig).gt.real(adzmx))
     .         ireset=1
         else
            if (real(adxmx).gt.0. .and. real(axorig).gt.real(adxmx)) 
     .         ireset= -1
            if (real(adymx).gt.0. .and. real(ayorig).gt.real(adymx))
     .         ireset= -1
            if (real(adzmx).gt.0. .and. real(azorig).gt.real(adzmx))
     .         ireset= -1
         end if
         if (irotat(nbl) .eq. 1) then 
            if (real(adthxmx).gt.0. .and. real(athetx).gt.real(adthxmx))
     .         ireset=1
            if (real(adthymx).gt.0. .and. real(athety).gt.real(adthymx))
     .         ireset=1
            if (real(adthzmx).gt.0. .and. real(athetz).gt.real(adthzmx))
     .         ireset=1
         else
            if (real(adthxmx).gt.0. .and. real(athetx).gt.real(adthxmx))
     .         ireset= -1
            if (real(adthymx).gt.0. .and. real(athety).gt.real(adthymx))
     .         ireset= -1
            if (real(adthzmx).gt.0. .and. real(athetz).gt.real(adthzmx))
     .         ireset= -1
         end if
c
         if (ireset .ge. 0) iresetb(nbl) = ireset
         if (ireset .lt. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),99) nbl
         end if
c
c        reset block nbl
c
         if (ireset .gt. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100) nbl,igridg(nbl)
c
c           reset translation
c
            if (itrans(nbl).gt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),101) real(xorig(nbl)),
     .               real(yorig(nbl)),real(zorig(nbl))
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),102) real(dxmx(nbl)),
     .               real(dymx(nbl)),real(dzmx(nbl))
               xorg = xorig(nbl)
               yorg = yorig(nbl)
               zorg = zorig(nbl)
               if (abs(real(utrans(nbl))) .gt. 0.) then
                  t2max      = ccabs(dxmx(nbl))/ccabs(utrans(nbl))
                  time2(nbl) = time2(nbl)  - t2max - dt
                  xorig(nbl) = xorig0(nbl) + time2(nbl)*utrans(nbl)
               end if
               if (abs(real(vtrans(nbl))) .gt. 0.) then
                  t2max      = ccabs(dymx(nbl))/ccabs(vtrans(nbl))
                  time2(nbl) = time2(nbl)  - t2max - dt
                  yorig(nbl) = yorig0(nbl) + time2(nbl)*vtrans(nbl)
               end if
               if (abs(real(wtrans(nbl))) .gt. 0.) then
                  t2max      = ccabs(dzmx(nbl))/ccabs(wtrans(nbl))
                  time2(nbl) = time2(nbl)  - t2max - dt
                  zorig(nbl) = zorig0(nbl) + time2(nbl)*wtrans(nbl)
               end if
               call grdmove(nbl,jdim,kdim,idim,w(lx),w(ly),w(lz),
     .                      xorg,yorg,zorg,
     .                      xorig(nbl),yorig(nbl),zorig(nbl),
     .                      thetax(nbl),thetay(nbl),thetaz(nbl))
            end if
c
c           reset rotation
c
            if (irotat(nbl).gt.0) then
               thx = thetax(nbl)*radtodeg
               thy = thetay(nbl)*radtodeg
               thz = thetaz(nbl)*radtodeg
               thxmx = dthxmx(nbl)*radtodeg
               thymx = dthymx(nbl)*radtodeg
               thzmx = dthzmx(nbl)*radtodeg
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),103) real(thx),real(thy),real(thz)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),104) real(thxmx),real(thymx),
     .               real(thzmx)
               xorg = xorig(nbl)
               yorg = yorig(nbl)
               zorg = zorig(nbl)
               thtx = thetax(nbl)
               thty = thetay(nbl)
               thtz = thetaz(nbl)
               if (abs(real(omegax(nbl))) .gt. 0.) then
                  t2max       = ccabs(dthxmx(nbl))/ccabs(omegax(nbl))
                  time2(nbl)  = time2(nbl) - t2max - dt
                  thetax(nbl) = time2(nbl)*omegax(nbl)
               end if
               if (abs(real(omegay(nbl))) .gt. 0.) then
                  t2max       = ccabs(dthymx(nbl))/ccabs(omegay(nbl))
                  time2(nbl)  = time2(nbl) - t2max - dt
                  thetay(nbl) = time2(nbl)*omegay(nbl)
               end if
               if (abs(real(omegaz(nbl))) .gt. 0.) then
                  t2max       = ccabs(dthzmx(nbl))/ccabs(omegaz(nbl))
                  time2(nbl)  = time2(nbl) - t2max - dt
                  thetaz(nbl) = time2(nbl)*omegaz(nbl)
               end if
               dthtx  = thetax(nbl) - thtx 
               dthty  = thetay(nbl) - thty
               dthtz  = thetaz(nbl) - thtz
               call grdmove(nbl,jdim,kdim,idim,w(lx),w(ly),w(lz),
     .                      xorg,yorg,zorg,
     .                      xorig(nbl),yorig(nbl),zorig(nbl),
     .                      dthtx,dthty,dthtz)
c
c              We rotated the x,y,zs thru an extra delta t, since we 
c              will be updating them again by delta t at the end of this
c              routine (in the call to updateg).  For the Qs, however,
c              we don't want to rotate them by extra amount:
c
               dthtx  = dthtx + dt*omegax(nbl)
               dthty  = dthty + dt*omegay(nbl)
               dthtz  = dthtz + dt*omegaz(nbl)
c
c              interior values
               lqwk = 1
               call rotateq(jdim,kdim,idim,w(lq),wk(lqwk),1,idim,
     .                      1,jdim,1,kdim,dthtx,dthty,dthtz)
               do 40 l=1,jdim*kdim*idim*5
               w(l+lq-1) = wk(l+lqwk-1)
   40          continue
c              qj0 values
               call rotateq0(kdim,idim-1,w(lqj0),wk(lqwk),1,kdim,
     .                       1,idim-1,dthtx,dthty,dthtz)
               do 45 l=1,kdim*(idim-1)*5*4
               w(l+lqj0-1) = wk(l+lqwk-1)
   45          continue
c              qk0 values
               call rotateq0(jdim,idim-1,w(lqk0),wk(lqwk),1,jdim,
     .                       1,idim-1,dthtx,dthty,dthtz)
               do 50 l=1,jdim*(idim-1)*5*4
               w(l+lqk0-1) = wk(l+lqwk-1)
   50          continue
c              qi0 values
               call rotateq0(jdim,kdim,w(lqi0),wk(lqwk),1,jdim,
     .                       1,kdim,dthtx,dthty,dthtz)
               do 55 l=1,jdim*kdim*5*4
               w(l+lqi0-1) = wk(l+lqwk-1)
   55          continue
c              qc0 values for second-order time advancement
               if (abs(ita).eq.2) then
                  call rotateq(jdim,kdim,idim-1,w(lqc0),wk(lqwk),
     .                         1,idim-1,1,jdim,1,kdim,dthtx,dthty,dthtz)
                  do 60 l=1,jdim*kdim*(idim-1)*5
                  w(l+lqc0-1) = wk(l+lqwk-1)
   60             continue
               end if
c              qb values for chimera scheme
               if (iovrlp(nbl).eq.1) then
                  call rotateqb(nbl,dthtx,dthty,dthtz,maxbl,iitot,ibcg,
     .                          lig,lbg,ibpntsg,iipntsg,qb)
               end if
            end if
c
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),105)
c
#if defined DIST_MPI
         end if
c
#endif
         end if
#   ifdef FASTIO
      call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl,15)
#   else
      call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl)
#   endif
c
      end if
c
    5 continue
c
c     check to see if moment center needs to be reset
c
#if defined DIST_MPI
      if (myid.eq.myhost) then
c
#endif
      iresetmc = 0
c
      if (itransmc.eq.1 .or. irotatmc.eq.1) then 
         adxmx   = ccabs(dxmxmc)
         adymx   = ccabs(dymxmc)
         adzmx   = ccabs(dzmxmc)
         adthxmx = ccabs(dthxmxmc)
         adthymx = ccabs(dthymxmc)
         adthzmx = ccabs(dthzmxmc)
         axorig  = ccabs(xorigmc)
         ayorig  = ccabs(yorigmc)
         azorig  = ccabs(zorigmc)
         athetx  = ccabs(thetaxmc)
         athety  = ccabs(thetaymc)
         athetz  = ccabs(thetazmc)
         if (real(adxmx).gt.0. .and. real(axorig).gt.real(adxmx))
     .       iresetmc = 1
         if (real(adymx).gt.0. .and. real(ayorig).gt.real(adymx))
     .       iresetmc = 1
         if (real(adzmx).gt.0. .and. real(azorig).gt.real(adzmx))
     .       iresetmc = 1
         if (real(adthxmx).gt.0. .and. real(athetx).gt.real(adthxmx))
     .       iresetmc = 1
         if (real(adthymx).gt.0. .and. real(athety).gt.real(adthymx))
     .       iresetmc = 1
         if (real(adthzmx).gt.0. .and. real(athetz).gt.real(adthzmx))
     .       iresetmc = 1
         if (iresetmc .gt. 0) then
            write(11,*)
            write(11,106)
            if (itransmc.gt.0) then
               write(11,101) real(xorigmc),real(yorigmc),real(zorigmc)
               write(11,102) real(dxmxmc),real(dymxmc),real(dzmxmc)
               if (abs(real(utransmc)) .gt. 0.) then
                  t2maxmc = ccabs(dxmxmc)/ccabs(utransmc)
                  time2mc = time2mc  - t2maxmc - dt
                  xorigmc = xorig0mc + time2mc*utransmc
                  xmc     = xmc0     + time2mc*utransmc
               end if
               if (abs(real(vtransmc)) .gt. 0.) then
                  t2maxmc = ccabs(dymxmc)/ccabs(vtransmc)
                  time2mc = time2mc  - t2maxmc - dt
                  yorigmc = yorig0mc + time2mc*vtransmc
                  ymc     = ymc0     + time2mc*vtransmc
               end if
               if (abs(real(wtransmc)) .gt. 0.) then
                  t2maxmc = ccabs(dzmxmc)/ccabs(wtransmc)
                  time2mc = time2mc  - t2maxmc - dt
                  zorigmc = zorig0mc + time2mc*wtransmc
                  zmc     = zmc0     + time2mc*wtransmc
               end if
            end if
            if (irotatmc.gt.0) then
               thx   = thetaxmc*radtodeg
               thy   = thetaymc*radtodeg
               thz   = thetazmc*radtodeg
               thxmx  = dthxmxmc*radtodeg
               thymx  = dthymxmc*radtodeg
               thzmx  = dthzmxmc*radtodeg
               thtxmc = thetaxmc
               thtymc = thetaymc
               thtzmc = thetazmc
               write(11,103) real(thx),real(thy),real(thz)
               write(11,104) real(thxmx),real(thymx),real(thzmx)
               if (abs(real(omegaxmc)) .gt. 0.) then
                  t2maxmc  = ccabs(dthxmxmc)/ccabs(omegaxmc)
                  time2mc  = time2mc  - t2maxmc - dt
                  thetaxmc = time2mc*omegaxmc
                  dthtxmc = thetaxmc - thtxmc
                  yml = ymc
                  zml = zmc
                  ca  = cos(dthtxmc)
                  sa  = sin(dthtxmc)
                  ymc = (yml-yorigmc)*ca-(zml-zorigmc)*sa+yorg
                  zmc = (yml-yorigmc)*sa+(zml-zorigmc)*ca+zorg
               end if
               if (abs(real(omegaymc)) .gt. 0.) then
                  t2maxmc  = ccabs(dthymxmc)/ccabs(omegaymc)
                  time2mc  = time2mc  - t2maxmc - dt
                  thetaymc = time2mc*omegaymc
                  dthtymc = thetaymc - thtymc 
                  xml = xmc
                  zml = zmc
                  ca  = cos(dthtymc)
                  sa  = sin(dthtymc)
                  xmc =  (xml-xorigmc)*ca+(zml-zorigmc)*sa+xorg
                  zmc = -(xml-xorigmc)*sa+(zml-zorigmc)*ca+zorg
               end if
               if (abs(real(omegazmc)) .gt. 0.) then
                  t2maxmc  = ccabs(dthzmxmc)/ccabs(omegazmc)
                  time2mc  = time2mc  - t2maxmc - dt
                  thetazmc = time2mc*omegazmc
                  dthtzmc = thetazmc - thtzmc
                  xml = xmc
                  yml = ymc
                  ca  = cos(dthtzmc)
                  sa  = sin(dthtzmc)
                  xmc = (xml-xorigmc)*ca-(yml-yorigmc)*sa+xorg
                  ymc = (xml-xorigmc)*sa+(yml-yorigmc)*ca+yorg
               end if
            end if
c
            write(11,105)
c
         end if
c
      end if
#if defined DIST_MPI
c
      end if
#endif
c
c     update time, grid and metrics
c
      ireset = 0
      do 80 nbl=1,nblock
#if defined DIST_MPI
      myidchk = mblk2nd(nbl)
      mytag = nbl
      if (myid.eq.myidchk) then
         call MPI_Send (iresetb(nbl), 1, MPI_INTEGER, myhost,
     &                  mytag, mycomm, ierr)
      else if (myid.eq.myhost) then
         call MPI_Recv (iresetb(nbl), 1, MPI_INTEGER, myidchk,
     &                  mytag, mycomm, istat, ierr)
      end if
#endif
      if (iresetb(nbl) .gt. 0) then
         time2(nbl) = time2(nbl) + dt
         ireset = 1
      end if
   80 continue
#if defined DIST_MPI
      if (myid.eq.myhost) then
#endif
      if (iresetmc .gt. 0) then
         time2mc = time2mc + dt
         ireset = 1
      end if
#if defined DIST_MPI
      end if
c
c     broadcast the hosts' reset flag so all processors will call
c     updateg if ireset = 1 on the host
c
      call MPI_Bcast (ireset, 1, MPI_INTEGER, myhost,
     &                mycomm, ierr)
#endif
c
      if (ireset .gt. 0) then
         irigb0   = 0
         irbtrim0 = 0
         call updateg(lw,lw2,w,mgwk,wk,nwork,iupdat,iseqr,maxbl,
     .                maxgr,maxseg,nbci0,nbcj0,nbck0,nbcidim,
     .                nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                nblock,levelg,igridg,utrans,vtrans,wtrans,
     .                omegax,omegay,omegaz,xorig,yorig,zorig,
     .                thetax,thetay,thetaz,rfreqt,rfreqr,xorig0,
     .                yorig0,zorig0,time2,thetaxl,thetayl,thetazl,
     .                itrans,irotat,idefrm,ncgg,iadvance,
     .                nou,bou,nbuf,ibufdim,myid,myhost,mycomm,
     .                mblk2nd,irigb0,irbtrim0,nt)
      end if
c
   99 format(15h WARNING: block,i4,20h will not be reset -,
     .       44h resetting allowed only if itrans/irotat = 1)
  100 format(28h resetting position of block,i4,
     .       6h (grid,i4,1h))
  101 format(39h   current x,y,z displacements:        ,
     .       3(f8.3,2x))
  102 format(39h   max. allowable x,y,z displacements: ,
     .       3(f8.3,2x))
  103 format(27h   current x,y,z rotational,
     .       23h displacements:        ,3(f8.3,2x))
  104 format(34h   max. allowable x,y,z rotational,
     .       16h displacements: ,3(f8.3,2x))
  105 format(' ')
  106 format(24h resetting moment center)
c
      return
      end
