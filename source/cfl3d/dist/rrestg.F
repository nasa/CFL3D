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
      subroutine rrestg(nbl,igrid,jdim,kdim,idim,x,y,z,xnm2,
     .                  ynm2,znm2,deltj,deltk,delti,qc0,
     .                  nflagg,iuns,utrans,vtrans,wtrans,
     .                  omegax,omegay,omegaz,xorig,yorig,zorig,
     .                  dxmx,dymx,dzmx,dthxmx,dthymx,
     .                  dthzmx,thetax,thetay,thetaz,rfreqt,
     .                  rfreqr,xorig0,yorig0,zorig0,time2,
     .                  thetaxl,thetayl,thetazl,itrans,irotat,idefrm,
     .                  utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                  xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                  rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,kcsi,kcsf,
     .                  freq,gmass,damp,x0,gf0,nmds,maxaes,aesrfdat,
     .                  perturb,myhost,myid,mycomm,mblk2nd,maxbl,
     .                  ibufdim,nbuf,bou,nou,nsegdfrm,idfrmseg,
     .                  iaesurf,maxsegdg,wk,nwork,idima,jdima,kdima,
     .                  tursav2,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Read the restart file to get the required info for a 
c     dynamic mesh (nflagg=1); read the restart file to get qc0 data 
c     for 2nd order accurate restart (in time) (nflagg=0).
c
c     For rigid (non-deforming) grids undergoing prescribed motion,
c     only the dynamic mesh parameters are needed, from which the
c     last grid orientation (at t(n)) is reconstructed. 
c
c     For deforming grids, or grids without prescribed motion, a
c     more complete output of the state of the entire grid system 
c     is needed: the last grid definition (at t(n)) for first order
c     in time, or, for second order in time, the grid definition at
c     t(n-1) and the boundary displacements at t(n), from which the
c     the last grid definition (at t(n)) is reconstructed.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
      dimension istat(MPI_STATUS_SIZE)
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
#endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf),wk(nwork)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension xnm2(jdim,kdim,idim),ynm2(jdim,kdim,idim),
     .          znm2(jdim,kdim,idim)
      dimension qc0(jdim,kdim,idim-1,5),tursav2(jdim,kdim,idim,2*nummem)
      dimension deltj(kdim,idim,3,2),deltk(jdim,idim,3,2),
     .          delti(jdim,kdim,3,2)
      dimension utrans(maxbl),vtrans(maxbl),wtrans(maxbl),
     .          omegax(maxbl),omegay(maxbl),omegaz(maxbl),
     .          xorig(maxbl),yorig(maxbl),zorig(maxbl),
     .          thetax(maxbl),thetay(maxbl),thetaz(maxbl),
     .          rfreqt(maxbl),rfreqr(maxbl),xorig0(maxbl),
     .          yorig0(maxbl),zorig0(maxbl),time2(maxbl),
     .          thetaxl(maxbl),thetayl(maxbl),thetazl(maxbl),
     .          itrans(maxbl),irotat(maxbl),idefrm(maxbl)
      dimension dxmx(maxbl),dymx(maxbl),dzmx(maxbl),dthxmx(maxbl),
     .          dthymx(maxbl),dthzmx(maxbl)
      dimension mblk2nd(maxbl)
      dimension utrnsae(maxbl,maxsegdg),vtrnsae(maxbl,maxsegdg),
     .          wtrnsae(maxbl,maxsegdg),omgxae(maxbl,maxsegdg),
     .          omgyae(maxbl,maxsegdg),omgzae(maxbl,maxsegdg),
     .          xorgae(maxbl,maxsegdg),yorgae(maxbl,maxsegdg),
     .          zorgae(maxbl,maxsegdg),thtxae(maxbl,maxsegdg),
     .          thtyae(maxbl,maxsegdg),thtzae(maxbl,maxsegdg),
     .          rfrqtae(maxbl,maxsegdg),rfrqrae(maxbl,maxsegdg)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg),
     .          iaesurf(maxbl,maxsegdg)
      dimension freq(nmds,maxaes),gmass(nmds,maxaes),x0(2*nmds,maxaes),
     .          gf0(2*nmds,maxaes),damp(nmds,maxaes),
     .          perturb(nmds,maxaes,4)
      dimension aesrfdat(5,maxaes)
c
      common /cfl/ dt0,dtold
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /maxiv/ ivmx
      common /motionmc/ xmc0,ymc0,zmc0,utransmc,vtransmc,wtransmc,
     .                  omegaxmc,omegaymc,omegazmc,xorigmc,yorigmc,
     .                  zorigmc,xorig0mc,yorig0mc,zorig0mc,thetaxmc,
     .                  thetaymc,thetazmc,dxmxmc,dymxmc,dzmxmc,
     .                  dthxmxmc,dthymxmc,dthzmxmc,rfreqtmc,
     .                  rfreqrmc,itransmc,irotatmc,time2mc
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /igrdtyp/ ip3dgrd,ialph
      common /cgns/ icgns,iccg,ibase,nzones,nsoluse,irind,jrind,krind
      common /twod/ i2d
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt
c
#if defined DIST_MPI
      nd_dest = mblk2nd(nbl)
c
c     set baseline tag values
c
      ioffset    = maxbl
      itag_qc0   = 1
      itag_tursav= itag_qc0   + ioffset
      itag_dmdat = itag_tursav+ ioffset
      itag_x     = itag_dmdat + ioffset
      itag_y     = itag_x     + ioffset
      itag_z     = itag_y     + ioffset
      itag_xnm2  = itag_z     + ioffset
      itag_ynm2  = itag_xnm2  + ioffset
      itag_znm2  = itag_ynm2  + ioffset
      itag_deltj = itag_znm2  + ioffset
      itag_deltk = itag_deltj + ioffset
      itag_delti = itag_deltk + ioffset
      itag_dtold = itag_delti + ioffset
#endif
c
      jdim1=jdim-1
      kdim1=kdim-1
      idim1=idim-1
c
      if (nflagg .eq. 0) then
c
c *******************************
c       read 2nd order time data
c *******************************
c
      if (myid.eq.myhost) then
      if (icgns .ne. 1) then
      write(11,50) nbl
   50 format(33h reading 2nd order time data from,
     .       23h restart file for block,i4)
c
         read(2) jdum,kdum,idum
         if (jdum.ne.jdim.and.kdum.ne.kdim.and.idum.ne.idim) then
            write(11,55) nbl
   55       format(43h error in 2nd order time section of restart,
     .             15h file for block,i3)
            write(11,65) jdum,kdum,idum,jdim,kdim,idim
   65       format(29h restart file has j x k x i :,3i5,/,
     .             29h should be        j x k x i :,3i5)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         read(2) ((((qc0(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),
     .               l=1,5)
c     itime2read=1=default=YES, read 2nd order time-accurate turb info
c     itime2read=0=NO, do not read (for use with older CFL3D 2nd-order
c                  time-accurate restart files, which was missing this info)
         if (itime2read .ne. 0) then
         read(2) dtold
         if (ivmx .ge. 4) then
         read(2) ((((tursav2(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),
     .               l=1,2)
         else
         read(2) ((((dum,j=1,jdim1),k=1,kdim1),i=1,idim1),l=1,2)
         end if
c   NOTE: for ivmx .ge. 30, currently you CANNOT restart to/from a
c   different ivmx; if this is ever allowed, then the following may 
c   read/write incorrectly!!!
         if (ivmx .ge. 30) then
          read(2) nummem_read
          if (nummem_read .ne. nummem) then
            write(11,'('' ERROR: nummem_read .ne. nummem in rrestg'')')
            write(11,'(''   cannot change ivisc with time-accurate'')')
            write(11,'(''   run for ivisc.ge.30'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
          end if
          read(2) ((((tursav2(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),
     .               l=3,nummem_read)
         end if
         end if
      else
#if defined CGNS
         write(901,'('' reading 2nd order time data from cgns'',
     +    '' file for block '',i4)') nbl
         call rsecord(iccg,ibase,igrid,idima,jdima,kdima,idim,
     +     jdim,kdim,wk,nsoluse,irind,jrind,krind,i2d,dt,
     +     ialph,qc0,tursav2,dtold,ivmx,nummem)
#endif
      end if
      if (dtold .ne. dt) then
        write(11,'('' time step has changed: old='',f12.6,
     .   '', new='',f12.6)') dtold,dt
        write(11,'('' ...1st step will be forced to be 1st order'',
     .   '' in time (to maintain order property)'')')
      end if
c
c        fill in edge values of qc0 array for safety before passing
c        the data to the appropriate node
c
         do l=1,5
            do i=1,idim1
               k = kdim
               do j=1,jdim1
                  qc0(j,k,i,l) = qiv(l)
               end do
               j=jdim
               do k=1,kdim
                  qc0(j,k,i,l) = qiv(l)
               end do
            end do
         end do
c        fill in edge values of tursav2 array for safety before passing
c        the data to the appropriate node
c
         if (ivmx .ge. 4) then
         do i=1,idim
            k = kdim
            do j=1,jdim1
              do l=1,nummem
               tursav2(j,k,i,l) = tur10(l)
               tursav2(j,k,i,nummem+l) = 0.
              enddo
            end do
            j=jdim
            do k=1,kdim
              do l=1,nummem
               tursav2(j,k,i,l) = tur10(l)
               tursav2(j,k,i,nummem+l) = 0.
              end do
            end do
         end do
         end if
c
      end if
#if defined DIST_MPI
c
c     send/receive data to/on the appropriate node
c
      jkim5 = jdim*kdim*idim1*5
      mytag = itag_qc0 + nbl
      if (myid .eq. myhost) then
         call MPI_Send(qc0,jkim5,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(qc0,jkim5,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
      end if
      if (ivmx .ge. 4) then
      jki4 = jdim*kdim*idim*2*nummem
      mytag = itag_tursav + nbl
      if (myid .eq. myhost) then
         call MPI_Send(tursav2,jki4,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(tursav2,jki4,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
      end if
      end if
      mytag = itag_dtold + nbl
      nvals=1
      if (myid .eq. myhost) then
         call MPI_Send(dtold,nvals,MY_MPI_REAL,
     .                 nd_dest,mytag,mycomm,ierr)
      else if (myid .eq. mblk2nd(nbl)) then
         call MPI_Recv(dtold,nvals,MY_MPI_REAL,
     .                 myhost,mytag,mycomm,istat,ierr)
      end if
c
#endif
c
      return
c
      else if (nflagg .eq. 1) then
c
c *******************************
c        read dynamic mesh data   
c *******************************
c
         if (myid.eq.myhost) then
         if (icgns .ne. 1) then
            read(2) iuns
         else
#if defined CGNS
         call readiuns(iccg,ibase,igrid,iuns)
#endif
         end if
         end if
#if defined DIST_MPI
c
c        send iuns to the appropriate node
c
         mytag = itag_dmdat
         if (myid.eq.myhost) then
            call MPI_Send (iuns, 1, MPI_INTEGER,
     .                     nd_dest, mytag, mycomm, ierr)
         else if (myid.eq.nd_dest) then
            call MPI_Recv (iuns, 1, MPI_INTEGER,
     .               myhost, mytag, mycomm, istat, ierr)
         end if
#endif
c
         if (iuns .ne. 0) then
c
            if (myid.eq.myhost) then
            if (icgns .ne. 1) then
               write(11,105) nbl
  105          format(31h reading dynamic mesh data from,
     .                23h restart file for block,i4)
c
               read(2) jdum,kdum,idum
               if (jdum.ne.jdim.and.kdum.ne.kdim.and.idum.ne.idim) then
                  write(11,100) nbl
  100             format(41h error in dynamic mesh section of restart,
     .                   15h file for block,i3)
                  write(11,101) jdum,kdum,idum,jdim,kdim,idim
  101             format(29h restart file has j x k x i :,3i5,/,
     .                   29h should be        j x k x i :,3i5)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
            end if
            end if
#if defined DIST_MPI
c
c           move a copy of the motion common block from the owner
c           of nbl to the host, to allow data read from input file
c           to be checked against data read from the restart file
c
            mytag = itag_dmdat + nbl
            nval  = 58 + 20*nsegdfrm(nbl)
            if (myid.eq.mblk2nd(nbl)) then
               idir = 0
               call mvdat(nbl,idir,maxbl,utrans,vtrans,wtrans,omegax,
     .                    omegay,omegaz,xorig,yorig,zorig,dxmx,dymx,
     .                    dzmx,dthxmx,dthymx,dthzmx,thetax,thetay,
     .                    thetaz,rfreqt,rfreqr,xorig0,yorig0,zorig0,
     .                    time2,thetaxl,thetayl,thetazl,itrans,irotat,
     .                    idefrm,utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,
     .                    omgzae,xorgae,yorgae,zorgae,thtxae,thtyae,
     .                    thtzae,rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                    kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                    aesrfdat,perturb,nsegdfrm,idfrmseg,iaesurf,
     .                    maxsegdg,wk,nwork)
               call MPI_Send (wk, nval, MY_MPI_REAL,
     .                        myhost, mytag, mycomm,ierr)
            else if (myid.eq.myhost) then
               call MPI_Recv (wk, nval, MY_MPI_REAL,
     .                        nd_dest, mytag, mycomm, istat, ierr)
               idir = 1
               call mvdat(nbl,idir,maxbl,utrans,vtrans,wtrans,omegax,
     .                    omegay,omegaz,xorig,yorig,zorig,dxmx,dymx,
     .                    dzmx,dthxmx,dthymx,dthzmx,thetax,thetay,
     .                    thetaz,rfreqt,rfreqr,xorig0,yorig0,zorig0,
     .                    time2,thetaxl,thetayl,thetazl,itrans,irotat,
     .                    idefrm,utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,
     .                    omgzae,xorgae,yorgae,zorgae,thtxae,thtyae,
     .                    thtzae,rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                    kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                    aesrfdat,perturb,nsegdfrm,idfrmseg,iaesurf,
     .                    maxsegdg,wk,nwork)
            end if
#endif
c
            if (myid.eq.myhost) then
            if (icgns .ne. 1) then
               read(2)
     .            itrans1,rfreqt1,xorig1,yorig1,zorig1,xorig01,yorig01,
     .            zorig01,utrans1,vtrans1,wtrans1,dxmx1,dymx1,dzmx1,
     .            itransmc1,rfreqtmc1,xorigmc1,yorigmc1,zorigmc1,
     .            xorig0mc1,yorig0mc1,zorig0mc1,utransmc1,vtransmc1,
     .            wtransmc1,xmc1,ymc1,zmc1,dxmxmc1,dymxmc1,dzmxmc1,
     .            irotat1,rfreqr1,thetax1,thetay1,thetaz1,
     .            omegax1,omegay1,omegaz1,dthxmx1,dthymx1,dthzmx1,
     .            irotatmc1,rfreqrmc1,thetaxmc1,thetaymc1,thetazmc1,
     .            omegaxmc1,omegaymc1,omegazmc1,dthxmxmc1,dthymxmc1,
     .            dthzmxmc1,time21,time2mc1,dt1
            else
#if defined CGNS
            write(901,'('' reading dynamic mesh data from cgns file'',
     +        '' for block'',i4)') nbl
            call rgrdmov(iccg,ibase,igrid,ialph,
     .            itrans1,rfreqt1,xorig1,yorig1,zorig1,xorig01,yorig01,
     .            zorig01,utrans1,vtrans1,wtrans1,dxmx1,dymx1,dzmx1,
     .            itransmc1,rfreqtmc1,xorigmc1,yorigmc1,zorigmc1,
     .            xorig0mc1,yorig0mc1,zorig0mc1,utransmc1,vtransmc1,
     .            wtransmc1,xmc1,ymc1,zmc1,dxmxmc1,dymxmc1,dzmxmc1,
     .            irotat1,rfreqr1,thetax1,thetay1,thetaz1,
     .            omegax1,omegay1,omegaz1,dthxmx1,dthymx1,dthzmx1,
     .            irotatmc1,rfreqrmc1,thetaxmc1,thetaymc1,thetazmc1,
     .            omegaxmc1,omegaymc1,omegazmc1,dthxmxmc1,dthymxmc1,
     .            dthzmxmc1,time21,time2mc1,dt1)
#endif
            end if
c
c              fix-up for change in moving grid parameters: 
c              change the block time counter (time2) to insure
c              continuous displacement of grid at restart; this also
c              insures that the grid resetting option (if used) will 
c              remain in sync. note: the value of time21 from the
c              restart file is first decremented by dt1 since the last 
c              thing that is done in subroutine mgblk before writing 
c              the restart file is to increment time2, i.e the restart 
c              file contains a displacement and a time21 that are 
c              actually out of sync by dt1 (as is always the case
c              *just* before the grid is updated). After time2 is 
c              is scaled to maintain continuous displacement, dt1
c              is added to the new time2.
c
               if (irotat(nbl).eq.1) then
                  if (abs(real(omegax1)).gt.0. .and. 
     .               real(omegax(nbl)).ne.real(omegax1)) then
                     if (real(omegay1).eq.0. .and.
     .                  real(omegaz1).eq.0.) then
                        time21 = time21 - dt1
                        time211 = omegax1/omegax(nbl)*time21
                        write(11,201) nbl,igrid
                        write(11,202) real(omegax(nbl)),real(omegax1)
                        write(11,203) real(time211),real(time21)
                        time21 = time211 + dt1
 201                    format(/,3x,16h WARNING: block ,i3,
     .                         7h (grid ,i3,1h),
     .                         21h rotational speed was,
     .                         22h changed upon restart:)
 202                    format(3x,14h new omegax = ,f10.4,
     .                         3x,14h old omegax = ,f10.4)
 203                    format(3x,31h to compensate, adjusting block,
     .                         21h time counter, time2:,/,
     .                         3x,14h new time2  = ,f10.4,
     .                         3x,14h old time2  = ,f10.4,/)
                     else
                        write(11,204)
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
 204                    format(/,3x,26h stopping...rotation speed,
     .                         34h change in more than one direction,
     .                         /,3x,29h not allowed between restarts) 
                     end if
                  end if
                  if (abs(real(omegay1)).gt.0. .and. 
     .               real(omegay(nbl)).ne.real(omegay1)) then
                     if (real(omegax1).eq.0. .and.
     .                  real(omegaz1).eq.0.) then
                        time21 = time21 - dt1
                        time211 = omegay1/omegay(nbl)*time21
                        write(11,201) nbl,igrid
                        write(11,205) real(omegay(nbl)),real(omegay1)
                        write(11,203) real(time211),real(time21)
                        time21 = time211 + dt1
 205                    format(3x,14h new omegay = ,f10.4,
     .                         3x,14h old omegay = ,f10.4)
                     else
                        write(11,204)
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                  end if
                  if (abs(real(omegaz1)).gt.0. .and. 
     .               real(omegaz(nbl)).ne.real(omegaz1)) then
                     if (real(omegax1).eq.0. .and.
     .                  real(omegay1).eq.0.) then
                        time21 = time21 - dt1
                        time211 = omegaz1/omegaz(nbl)*time21
                        write(11,201) nbl,igrid
                        write(11,206) real(omegaz(nbl)),real(omegaz1)
                        write(11,203) real(time211),real(time21)
                        time21 = time211 + dt1
 206                    format(3x,14h new omegaz = ,f10.4,
     .                         3x,14h old omegaz = ,f10.4)
                     else
                        write(11,204)
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                  end if
               end if
               if (irotatmc.eq.1) then
                  if (abs(real(omegaxmc1)).gt.0. .and.
     .               real(omegaxmc).ne.real(omegaxmc1)) then
                     if (real(omegaymc1).eq.0. .and.
     .                  real(omegazmc1).eq.0.) then
                        time2mc1 = time2mc1 - dt1
                        time2mc11 = omegaxmc1/omegaxmc*time2mc1
                        write(11,301) 
                        write(11,302) real(omegaxmc),real(omegaxmc1)
                        write(11,303) real(time2mc11),real(time2mc1)
                        time2mc1 = time2mc11 + dt
 301                    format(/,3x,23h WARNING: moment center,
     .                         34h rotational speed was changed upon,
     .                          8h restart:)
 302                    format(3x,16h new omegaxmc = ,f10.4,
     .                         3x,16h old omegaxmc = ,f10.4)
 303                    format(3x,32h to compensate, adjusting moment,
     .                         30h center time counter, time2mc:,/,
     .                         3x,16h new time2mc  = ,f10.4,
     .                         3x,16h old time2mc  = ,f10.4,/)
                     else
                        write(11,304)
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
 304                    format(/,3x,25h stopping...moment center,
     .                           25h rotation speed change in,
     .                           28h more than one direction not,
     .                         /,3x,25h allowed between restarts)
                     end if
                  end if
                  if (abs(real(omegaymc1)).gt.0. .and.
     .               real(omegaymc).ne.real(omegaymc1)) then
                     if (real(omegaxmc1).eq.0. .and.
     .                  real(omegazmc1).eq.0.) then
                        time2mc1 = time2mc1 - dt1
                        time2mc11 = omegaymc1/omegaymc*time2mc1
                        write(11,301)
                        write(11,305) real(omegaymc),real(omegaymc1)
                        write(11,303) real(time2mc11),real(time2mc1)
                        time2mc1 = time2mc11 + dt1
 305                    format(3x,16h new omegaymc = ,f10.4,
     .                         3x,16h old omegaymc = ,f10.4)
                     else
                        write(11,304)
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                  end if
                  if (abs(real(omegazmc1)).gt.0. .and. 
     .               real(omegazmc).ne.real(omegazmc1)) then
                     if (real(omegaxmc1).eq.0. .and.
     .                  real(omegaymc1).eq.0.) then
                        time2mc1 = time2mc1 - dt1
                        time2mc11 = omegazmc1/omegazmc*time2mc1
                        write(11,301)
                        write(11,306) real(omegazmc),real(omegazmc1)
                        write(11,303) real(time2mc11),real(time2mc1)
                        time2mc1 = time2mc11 + dt1
 306                    format(3x,16h new omegazmc = ,f10.4,
     .                         3x,16h old omegazmc = ,f10.4)
                     else
                        write(11,304)
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                  end if
               end if
c
c              set old position/time values to current values
c
               time2(nbl)    = time21
               xorig(nbl)    = xorig1
               yorig(nbl)    = yorig1
               zorig(nbl)    = zorig1
               thetax(nbl)   = thetax1
               thetay(nbl)   = thetay1
               thetaz(nbl)   = thetaz1
               time2mc  = time2mc1
               xorigmc  = xorigmc1
               yorigmc  = yorigmc1
               zorigmc  = zorigmc1
               xmc      = xmc1
               ymc      = ymc1
               zmc      = zmc1
               thetaxmc = thetaxmc1
               thetaymc = thetaymc1
               thetazmc = thetazmc1
c
            end if
c
            if (idefrm(nbl).gt.0) then
               if (myid.eq.myhost) then
                  if (icgns .ne. 1) then
                  read(2) idefrm(nbl),nsegdfrm(nbl)
                  do is=1,nsegdfrm(nbl)
                     read(2) idfrmsdum,utrnsae(nbl,is),
     .                       vtrnsae(nbl,is),wtrnsae(nbl,is),
     .                       omgxae(nbl,is),omgyae(nbl,is),
     .                       omgzae(nbl,is),xorgae(nbl,is),
     .                       yorgae(nbl,is),zorgae(nbl,is),
     .                       thtxae(nbl,is),thtyae(nbl,is),
     .                       thtzae(nbl,is),rfrqtae(nbl,is),
     .                       rfrqrae(nbl,is),icsi(nbl,is),
     .                       icsf(nbl,is),jcsi(nbl,is),jcsf(nbl,is),
     .                       kcsi(nbl,is),kcsf(nbl,is)
                  end do
                  if (abs(ita) .gt. 1) then
                     read(2) 
     .                    (((xnm2(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
     .                    (((ynm2(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
     .                    (((znm2(j,k,i),j=1,jdim),k=1,kdim),i=1,idim)
                  else
                     read(2)
     .                    (((xdum,j=1,jdim),k=1,kdim),i=1,idim),
     .                    (((ydum,j=1,jdim),k=1,kdim),i=1,idim),
     .                    (((zdum,j=1,jdim),k=1,kdim),i=1,idim)
                  end if
                  read(2) (((x(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
     .                    (((y(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
     .                    (((z(j,k,i),j=1,jdim),k=1,kdim),i=1,idim)
                  if(ndgrd.ne.0.and.iunst.gt.1) then
                   read(98,*)(((x(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .                       (((y(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .                       (((z(j,k,i),i=1,idim),j=1,jdim),k=1,kdim)
                  end if
                  else
#if defined CGNS
                  write(901,'('' reading deforming grid info from'',
     .             '' cgns file for block '',i4)') nbl
                  call getnsegdfrm(iccg,ibase,igrid,nsegdfrm(nbl))
                  if (abs(ita) .gt. 1) then
                  call rdeform(iccg,ibase,igrid,nsegdfrm(nbl),
     .            wk,jdima,kdima,idima,i2d,
     .            jdim,kdim,idim,ialph,idefrm(nbl),utrnsae(nbl,1),
     .            vtrnsae(nbl,1),wtrnsae(nbl,1),omgxae(nbl,1),
     .            omgyae(nbl,1),omgzae(nbl,1),xorgae(nbl,1),
     .            yorgae(nbl,1),zorgae(nbl,1),thtxae(nbl,1),
     .            thtyae(nbl,1),thtzae(nbl,1),rfrqtae(nbl,1),
     .            rfrqrae(nbl,1),icsi(nbl,1),icsf(nbl,1),jcsi(nbl,1),
     .            jcsf(nbl,1),kcsi(nbl,1),kcsf(nbl,1),
     .            x,y,z,xnm2,ynm2,znm2)
                  else
                  call rdeform(iccg,ibase,igrid,nsegdfrm(nbl),
     .            wk,jdima,kdima,idima,i2d,
     .            jdim,kdim,idim,ialph,idefrm(nbl),utrnsae(nbl,1),
     .            vtrnsae(nbl,1),wtrnsae(nbl,1),omgxae(nbl,1),
     .            omgyae(nbl,1),omgzae(nbl,1),xorgae(nbl,1),
     .            yorgae(nbl,1),zorgae(nbl,1),thtxae(nbl,1),
     .            thtyae(nbl,1),thtzae(nbl,1),rfrqtae(nbl,1),
     .            rfrqrae(nbl,1),icsi(nbl,1),icsf(nbl,1),jcsi(nbl,1),
     .            jcsf(nbl,1),kcsi(nbl,1),kcsf(nbl,1),
     .            x,y,z,wk,wk,wk)
                  end if
#endif
                  end if
               end if
#if defined DIST_MPI
               if (abs(ita) .gt. 1) then
                  jki = jdim*kdim*idim
                  if (myid.eq.myhost) then
                     mytag = itag_xnm2 + nbl
                     call MPI_Send (xnm2, jki, MY_MPI_REAL,
     .                              nd_dest, mytag, mycomm, ierr)
                     mytag = itag_ynm2 + nbl
                     call MPI_Send (ynm2, jki, MY_MPI_REAL,
     .                              nd_dest, mytag, mycomm, ierr)
                     mytag = itag_znm2 + nbl
                     call MPI_Send (znm2, jki, MY_MPI_REAL,
     .                              nd_dest, mytag, mycomm, ierr)
                  else if (myid.eq.mblk2nd(nbl)) then
                     mytag = itag_xnm2 + nbl
                     call MPI_Recv (xnm2, jki, MY_MPI_REAL,
     .                              myhost, mytag, mycomm, istat, ierr)
                     mytag = itag_ynm2 + nbl
                     call MPI_Recv (ynm2, jki, MY_MPI_REAL,
     .                              myhost, mytag, mycomm, istat, ierr)
                     mytag = itag_znm2 + nbl
                     call MPI_Recv (znm2, jki, MY_MPI_REAL,
     .                              myhost, mytag, mycomm, istat, ierr)
                  end if
               end if
               jki = jdim*kdim*idim
               if (myid.eq.myhost) then
                  mytag = itag_x + nbl
                  call MPI_Send (x, jki, MY_MPI_REAL,
     .                           nd_dest, mytag, mycomm, ierr)
                  mytag = itag_y + nbl
                  call MPI_Send (y, jki, MY_MPI_REAL,
     .                           nd_dest, mytag, mycomm, ierr)
                  mytag = itag_z + nbl
                  call MPI_Send (z, jki, MY_MPI_REAL,
     .                           nd_dest, mytag, mycomm, ierr)
               else if (myid.eq.mblk2nd(nbl)) then
                  mytag = itag_x + nbl
                  call MPI_Recv (x, jki, MY_MPI_REAL,
     .                           myhost, mytag, mycomm, istat, ierr)
                  mytag = itag_y + nbl
                  call MPI_Recv (y, jki, MY_MPI_REAL,
     .                           myhost, mytag, mycomm, istat, ierr)
                  mytag = itag_z + nbl
                  call MPI_Recv (z, jki, MY_MPI_REAL,
     .                           myhost, mytag, mycomm, istat, ierr)
               end if
#endif
            end if
#if defined DIST_MPI
c
c           for distributed code, send updated motion block back
c           to the owner
c
            mytag = itag_dmdat + nbl
            nval  = 58 + 20*nsegdfrm(nbl)
            if (myid.eq.myhost) then
               idir = 0
               call mvdat(nbl,idir,maxbl,utrans,vtrans,wtrans,omegax,
     .                    omegay,omegaz,xorig,yorig,zorig,dxmx,dymx,
     .                    dzmx,dthxmx,dthymx,dthzmx,thetax,thetay,
     .                    thetaz,rfreqt,rfreqr,xorig0,yorig0,zorig0,
     .                    time2,thetaxl,thetayl,thetazl,itrans,irotat,
     .                    idefrm,utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,
     .                    omgzae,xorgae,yorgae,zorgae,thtxae,thtyae,
     .                    thtzae,rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                    kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                    aesrfdat,perturb,nsegdfrm,idfrmseg,iaesurf,
     .                    maxsegdg,wk,nwork)
               call MPI_Send (wk, nval, MY_MPI_REAL,
     .                        nd_dest, mytag, mycomm, ierr)
            else if (myid.eq.mblk2nd(nbl)) then
               call MPI_Recv (wk, nval, MY_MPI_REAL,
     .                        myhost, mytag, mycomm, istat, ierr)
               idir = 1
               call mvdat(nbl,idir,maxbl,utrans,vtrans,wtrans,omegax,
     .                    omegay,omegaz,xorig,yorig,zorig,dxmx,dymx,
     .                    dzmx,dthxmx,dthymx,dthzmx,thetax,thetay,
     .                    thetaz,rfreqt,rfreqr,xorig0,yorig0,zorig0,
     .                    time2,thetaxl,thetayl,thetazl,itrans,irotat,
     .                    idefrm,utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,
     .                    omgzae,xorgae,yorgae,zorgae,thtxae,thtyae,
     .                    thtzae,rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                    kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                    aesrfdat,perturb,nsegdfrm,idfrmseg,iaesurf,
     .                    maxsegdg,wk,nwork)
            end if
#endif
c
         end if
c
      end if
c
      return
      end
