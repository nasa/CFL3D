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
      subroutine setup(lw,lw2,w,mgwk,wk,nwork,iwk,iwork,itest,
     .                 jtest,ktest,maxbl,mxbli,maxgr,maxseg,nsub1,
     .                 maxxe,intmax,iitot,ncycmax,lwdat,lig,lbg,iovrlp,
     .                 qb,nblock,iviscg,jdimg,kdimg,idimg,utrans,
     .                 vtrans,wtrans,omegax,omegay,omegaz,xorig,
     .                 yorig,zorig,dxmx,dymx,dzmx,dthxmx,dthymx,
     .                 dthzmx,thetax,thetay,thetaz,rfreqt,rfreqr,
     .                 xorig0,yorig0,zorig0,time2,thetaxl,thetayl,
     .                 thetazl,itrans,irotat,idefrm,bcvali,bcvalj,
     .                 bcvalk,nbci0,nbcidim,nbcj0,nbcjdim,
     .                 nbck0,nbckdim,ibcinfo,jbcinfo,kbcinfo,bcfilei,
     .                 bcfilej,bcfilek,ngrid,ncgg,nblg,iemg,inewgg,
     .                 rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,cmxw,cmyw,
     .                 cmzw,n_clcd,clcd,nblocks_clcd,blocks_clcd,
     .                 chdw,swetw,fmdotw,cfttotw,cftmomw,cftpw,
     .                 cftvw,rmstr,nneg,ntr,windex,
     .                 ninter,iindex,nblkpt,dthetxx,dthetyy,dthetzz,
     .                 iibg,kkbg,jjbg,ibcg,dxintg,dyintg,dzintg,iiig,
     .                 jjig,kkig,ibpntsg,iipntsg,mblk2nd,nou,bou,nbuf,
     .                 ibufdim,ireq_qb,igridg,bcfiles,mxbcfil,
     .                 utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                 xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                 rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                 kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                 aesrfdat,perturb,islavept,nslave,iskip,jskip,
     .                 kskip,bmat,stm,stmi,gforcn,gforcnm,xxn,
     .                 nsegdfrm,idfrmseg,iaesurf,maxsegdg,nmaster,
     .                 aehist,timekeep,inpl3d,nplots,nplot3d,levelg,
     .                 iadvance,xs,gforcs,xorgae0,yorgae0,zorgae0,
     .                 icouple,lfgm,nblk,limblk,isva,nblelst,
     .                 iskmax,jskmax,kskmax,ue,irdrea,nbli,nummem)
#   ifdef CMPLX
#   else
      use module_kwstm, only:kws_init
#   endif
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Read in the grid and restart data and calculate some
c     preliminary information, such as the metrics, for subsequent use.
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
      dimension istat(MPI_STATUS_SIZE)
#endif
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,avgg,avgq
      character*241 avgq2,avgq2pert,clcds,clcdp,output_dir
      character*80 filname
      character*80 bcfiles(mxbcfil)
      character*120 bou(ibufdim,nbuf)
      character*32 basedesired
c
      integer bcfilei,bcfilej,bcfilek
      integer isndrea(2)
      integer irdrea(maxgr)
c
      dimension nou(nbuf)
      dimension w(mgwk),wk(nwork),iwk(iwork),lw(65,maxbl),
     .          lwdat(maxbl,maxseg,6),bcdat(12),intpts(4),lw2(43,maxbl)
      dimension itest(maxgr),jtest(maxgr),ktest(maxgr)
      dimension lig(maxbl),lbg(maxbl),iovrlp(maxbl),qb(iitot,5,3),
     .          iibg(iitot),kkbg(iitot),jjbg(iitot),ibcg(iitot),
     .          dxintg(iitot),dyintg(iitot),dzintg(iitot),iiig(iitot),
     .          jjig(iitot),kkig(iitot),ibpntsg(maxbl,4),iipntsg(maxbl)
      dimension igridg(maxbl),iviscg(maxbl,3),jdimg(maxbl),kdimg(maxbl),
     .          idimg(maxbl)
      dimension utrans(maxbl),vtrans(maxbl),wtrans(maxbl),omegax(maxbl),
     .          omegay(maxbl),omegaz(maxbl),xorig(maxbl),yorig(maxbl),
     .          zorig(maxbl),dxmx(maxbl),dymx(maxbl),dzmx(maxbl),
     .          dthxmx(maxbl),dthymx(maxbl),dthzmx(maxbl),thetax(maxbl),
     .          thetay(maxbl),thetaz(maxbl),rfreqt(maxbl),rfreqr(maxbl),
     .          xorig0(maxbl),yorig0(maxbl),zorig0(maxbl),time2(maxbl),
     .          thetaxl(maxbl),thetayl(maxbl),thetazl(maxbl),
     .          itrans(maxbl),irotat(maxbl),idefrm(maxbl)
      dimension bcvali(maxbl,maxseg,12,2),bcvalj(maxbl,maxseg,12,2),
     .          bcvalk(maxbl,maxseg,12,2),nbci0(maxbl),nbcidim(maxbl),
     .          nbcj0(maxbl),nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          ibcinfo(maxbl,maxseg,7,2),jbcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2)
      dimension bcfilei(maxbl,maxseg,2),bcfilej(maxbl,maxseg,2),
     .          bcfilek(maxbl,maxseg,2)
      dimension ncgg(maxgr),nblg(maxgr),iemg(maxgr),inewgg(maxgr)
      integer blocks_clcd
      dimension rms(ncycmax),clw(ncycmax),cdw(ncycmax),cdpw(ncycmax),
     .          cdvw(ncycmax),cxw(ncycmax),cyw(ncycmax),czw(ncycmax),
     .          cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),chdw(ncycmax),
     .          clcd(2,n_clcd,ncycmax), blocks_clcd(2,nblocks_clcd),
     .          swetw(ncycmax),fmdotw(ncycmax),cfttotw(ncycmax),
     .          cftmomw(ncycmax),cftpw(ncycmax),cftvw(ncycmax),
     .          rmstr(ncycmax,nummem),nneg(ncycmax,nummem),
     .          timekeep(ncycmax)
      dimension windex(maxxe,2),iindex(intmax,6*nsub1+9),nblkpt(maxxe),
     .          dthetxx(intmax,nsub1),dthetyy(intmax,nsub1),
     .          dthetzz(intmax,nsub1)
      dimension mblk2nd(maxbl),ireq_qb(maxbl*8),mytag_qb(8)
      dimension utrnsae(maxbl,maxsegdg),vtrnsae(maxbl,maxsegdg),
     .          wtrnsae(maxbl,maxsegdg),omgxae(maxbl,maxsegdg),
     .          omgyae(maxbl,maxsegdg),omgzae(maxbl,maxsegdg),
     .          xorgae(maxbl,maxsegdg),yorgae(maxbl,maxsegdg),
     .          zorgae(maxbl,maxsegdg),thtxae(maxbl,maxsegdg),
     .          thtyae(maxbl,maxsegdg),thtzae(maxbl,maxsegdg),
     .          rfrqtae(maxbl,maxsegdg),rfrqrae(maxbl,maxsegdg)
      dimension xorgae0(maxbl,maxsegdg),yorgae0(maxbl,maxsegdg),
     .          zorgae0(maxbl,maxsegdg),icouple(maxbl,maxsegdg)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg),
     .          iaesurf(maxbl,maxsegdg)
      dimension freq(nmds,maxaes),gmass(nmds,maxaes),x0(2*nmds,maxaes),
     .          gf0(2*nmds,maxaes),damp(nmds,maxaes),
     .          perturb(nmds,maxaes,4),aehist(ncycmax,3,nmds,maxaes)
      dimension aesrfdat(5,maxaes),islavept(nslave,nmaster,5)
      dimension iskip(maxbl,500),jskip(maxbl,500),kskip(maxbl,500)
      dimension ue(3*nslave)
      dimension bmat(2*nmds,2*nmds,maxaes),stm(2*nmds,2*nmds,maxaes),
     .          stmi(2*nmds,2*nmds,maxaes),xxn(2*nmds,maxaes),
     .          gforcn(2*nmds,maxaes),gforcnm(2*nmds,maxaes),
     .          gforcs(2*nmds,maxaes),xs(2*nmds,maxaes)
      dimension inpl3d(nplots,11),levelg(maxbl),iadvance(maxbl)
      dimension nblk(2,mxbli),limblk(2,6,mxbli),
     .          isva(2,2,mxbli)
      dimension nblelst(maxbl,2)
      dimension iskmax(maxbl)
      dimension jskmax(maxbl)
      dimension kskmax(maxbl)
c
      common /cgns/ icgns,iccg,ibase,nzones,nsoluse,irind,jrind,krind
      common /complx/ xmach_img,alpha_img,beta_img,reue_img,tinf_img,
     .                geom_img,surf_img,xrotrate_img,yrotrate_img,
     .                zrotrate_img
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /ginfo2/ lq2avg,iskip_blocks,inc_2d(3),inc_coarse(3)
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /ncyct/ ncyctot
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /sminn/ isminc,ismincforce
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /wrestq/ irest,irest2
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
      common /filenam2/ avgq2,avgq2pert,clcds,clcdp,output_dir
      common /igrdtyp/ ip3dgrd,ialph
      common /mydist2/ nnodes,myhost,myid,mycomm
      common /elastic/ ndefrm,naesrf
      common /elastic_ss/ idef_ss 
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt 
      common /ghost/ irghost,iwghost
      common /rigidbody/ irigb,irbtrim
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
      common /gridtrans/ roll_angle
c
      if (myid.eq.myhost) then
      if (icgns .ne. 1) then
         nzones = ngrid
         open(unit=2,file=restrt,form='unformatted',status='unknown')
         rewind( 2)
         if( iclcd .eq. 1 .or. iclcd .eq. 2 ) then
            open(unit=102,file=clcds,form='unformatted',
     $           status='unknown')
            rewind(102)
         end if
      else
#if defined CGNS
c   Open database
      basedesired='Base'
      idimdesired=3
      call ropencgns(grid,basedesired,idimdesired,iccg,
     .     ibase,nzones)
#endif
      end if
      end if
c
      ntt     = 0
      ntq     = 0
      ntr     = 0
      time    = 0.
      time2mc = 0.
      xnumavg  = 0.
      xnumavg2  = 0.
      do 800 nbl=1,nblock
      time2(nbl) = 0.
  800 continue
      do nc=1,ncycmax
         timekeep(nc) = 0.
         if (naesrf.gt.0) then
            do iaes=1,naesrf
               nmodes = int(aesrfdat(5,iaes))
               do nm=1,nmodes
                  aehist(nc,1,nm,iaes) = 0.
                  aehist(nc,2,nm,iaes) = 0.
                  aehist(nc,3,nm,iaes) = 0.
               end do
            end do
         end if
      end do
      if(iunst.gt.1) then
        do ns = 1,3*nslave
          ue(ns) = 0.
        enddo
      end if
c
      lig(1) = 1
      lbg(1) = 1
c
      if (myid.eq.mblk2nd(1)) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),10)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
      end if
   10 format(36h***** BEGINNING INITIALIZATION *****)
c
c*******************************************************
c     read in grids
c*******************************************************
#if defined DIST_MPI
c
c     set baseline tag values
c
      ioffset   = maxbl
      itag_x    = 1
      itag_y    = itag_x    + ioffset
      itag_z    = itag_y    + ioffset
      itag_flag = itag_z    + ioffset
#endif
c
c     flag to read complex-valued grid
c
      igeom_img = 0
c
#   ifdef CMPLX
      if (real(geom_img).gt.0.) then
         igeom_img = 1
      end if
c
#   endif
      if (icgns .eq. 1 .and. myid.eq.myhost) then
c   Check number of zones
      if (nzones .ne. ngrid) then
        write(901,'('' inconsistency between cgns # of zones & input.'',
     + '' nzones,ngrid='',2i6)') nzones,ngrid
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      end if
c
      do 1000 igrid=1,ngrid
c
      iem  = iemg(igrid)
      nbl  = nblg(igrid)
      call lead(nbl,lw,lw2,maxbl)
c
      if (icgns .eq. 1) then
#if defined CGNS
        if (myid.eq.myhost) then
        if(nwork .lt. idim*jdim*kdim) then
          write(901,'('' not enough memory for cgns grid read.'')')
          write(901,'('' nwork in wk='',i6,''.  Needed = '',i6)')
     +      nwork,idim*jdim*kdim
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
        call getgrid(iccg,ibase,igrid,idim,jdim,kdim,wk,
     +    ialph,w(lx),w(ly),w(lz))
        end if
#   ifdef FASTIO
         mytag   = itag_flag + nbl
         nd_scrc = mblk2nd(nbl)
         if (mblk2nd(nbl).eq.myid) then
            isndrea(1)=myid
            isndrea(2)=igrid
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100) igrid,idim,jdim,kdim
            call MPI_Send (isndrea, 2, MPI_INTEGER,
     *                     myhost, mytag, mycomm, ierr)
         else if (myid.eq.myhost) then
            call MPI_Recv (isndrea, 2, MPI_INTEGER,
     *                     nd_scrc,mytag,mycomm,istat,ierr)
            irdrea(isndrea(2))=isndrea(1)
         end if
#   else
         if (myid.eq.mblk2nd(nbl)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100)igrid,idim,jdim,kdim
         end if
#   endif
#endif
      else
      if (ip3dgrd.eq.0) then
c	  
c        cfl3d type
c
         if (myid .eq. myhost) then
            read(1) jdum,kdum,idum  
         end if
c
  100       format(13h reading grid,i6,24h of dimensions (I/J/K) :,3i6)
#   ifdef FASTIO
         mytag   = itag_flag + nbl
         nd_scrc = mblk2nd(nbl)
         if (myid.eq.mblk2nd(nbl)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100) igrid,idim,jdim,kdim
            isndrea(1)=myid
            isndrea(2)=igrid
            call MPI_Send (isndrea, 2, MPI_INTEGER,
     *                     myhost, mytag, mycomm, ierr)
         else if (myid.eq.myhost) then
            call MPI_Recv (isndrea, 2, MPI_INTEGER,
     *                     nd_scrc,mytag,mycomm,istat,ierr)
            irdrea(isndrea(2))=isndrea(1)
         end if
#   else
         if (myid.eq.mblk2nd(nbl)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100) igrid,idim,jdim,kdim
         end if
#   endif
c
         if (myid .eq. myhost) then
            if (jdum.ne.jdim .or. kdum.ne.kdim .or. idum.ne.idim) then
               write(11,*) ' stopping....inconsistency in'
               write(11,*) ' grid data file'
               write(11,*) ' for grid',igrid
               write(11,*) ' input file: idim,jdim,kdim',
     .         idim,jdim,kdim
               write(11,*) ' grid  file: idim,jdim,kdim',
     .         idum,jdum,kdum
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
c
c           if it is desired to allow ialph to be used with cfl3d type 
c           grids as well as plot3d type grids (i.e. ialph = 0 means 
c           alpha measured in x-z plane; ialph = 1 means alpha measured
c           in x-y  plane), then comment out the following line in
c           subroutine global:      if (ngrid .gt. 0) ialph = 0
c
            if (ialph.ne.0 .and. iipv.eq.1) then
               write(11,'('' To use farfield point vortex, '',
     .          ''grid must be in x-z plane'')')
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            irr=0
            call rcfl(w(lx),w(ly),w(lz),jdim,kdim,idim,igrid,ialph,
     .                igeom_img,irr)
            if (irr .ne. 0) then
              write(11,*) ' Stopping... error reading grid...'
              write(11,*) ' (Common error:  grid not written in same',
     .        ' precision that CFL3D was compiled)'
              call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
      else
c
c        plot3d type     	  
c
         if (myid .eq.myhost) then
            if (igrid.eq.1) then
               read(1) ndum
               if (ndum.ne.ngrid) then
                  write(11,*) ' stopping....ngrid = ',ngrid,
     .            ' but grid file contains ',ndum,' grids'
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
               read(1) (itest(ll),jtest(ll),ktest(ll),ll=1,ngrid)
            end if
         end if
c
#   ifdef FASTIO
         mytag   = itag_flag + nbl
         nd_scrc = mblk2nd(nbl)
         if (mblk2nd(nbl).eq.myid) then
            isndrea(1)=myid
            isndrea(2)=igrid
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100) igrid,idim,jdim,kdim
            call MPI_Send (isndrea, 2, MPI_INTEGER,
     *                     myhost, mytag, mycomm, ierr)
         else if (myid.eq.myhost) then
            call MPI_Recv (isndrea, 2, MPI_INTEGER,
     *                     nd_scrc,mytag,mycomm,istat,ierr)
            irdrea(isndrea(2))=isndrea(1)
         end if
#   else
         if (myid.eq.mblk2nd(nbl)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),100)igrid,idim,jdim,kdim
         end if
#   endif
c
         if (myid .eq.myhost) then
            if (jtest(igrid).ne.jdim .or. ktest(igrid).ne.kdim .or. 
     .         itest(igrid).ne.idim) then
               write(11,*) ' stopping....inconsistency in',
     .                     ' grid data file'
               write(11,*) ' for grid',igrid
               write(11,*) ' input file: idim,jdim,kdim',
     .         idim,jdim,kdim
               write(11,*) ' grid  file: idim,jdim,kdim',
     .         itest(igrid),jtest(igrid),ktest(igrid)
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            irr=0
            call rp3d(w(lx),w(ly),w(lz),jdim,kdim,idim,igrid,ialph,
     .                igeom_img,irr)
            if (irr .ne. 0) then
              write(11,*) ' Stopping... error reading grid...'
              write(11,*) ' (Common error:  grid not written in same',
     .        ' precision that CFL3D was compiled)'
              call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
      end if
c
      end if 
c
c
      if (myid.eq.myhost) then
         if (igrid .eq. 1) then
            write(11,*)
            write(11,'(''input roll angle = '',f9.4,'' degrees'')') 
     .            real(roll_angle)
            roll_angle = roll_angle*atan(1.0)/45.0   ! radians
         end if
         if (real(ccabs(roll_angle)) .gt. 0.) then
            nbl = nblg(igrid)
            call grdmove(nbl,jdim,kdim,idim,w(lx),w(ly),w(lz),
     .                   xcentrot,ycentrot,zcentrot,
     .                   xcentrot,ycentrot,zcentrot,
     .                   roll_angle,0.0,0.0)
         end if
      end if
c
#if defined DIST_MPI
      nd_dest = mblk2nd(nbl)
      nvals   = jdim*kdim*idim
c
      mytag = itag_x + nbl
      if (myid .eq. myhost) then
         call MPI_Send (w(lx), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (w(lx), nvals, MY_MPI_REAL,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_y + nbl
      if (myid .eq. myhost) then
         call MPI_Send (w(ly), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (w(ly), nvals, MY_MPI_REAL,
     .                 myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_z + nbl
      if (myid .eq. myhost) then
         call MPI_Send (w(lz), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (w(lz), nvals, MY_MPI_REAL,
     .                 myhost, mytag, mycomm, istat, ierr)
      end if
#endif
c
c     store off original grid; grid displacements are added to
c     the original grid if mesh deforms
c
      if (idefrm(nbl).gt.0) then
         if (myid.eq.mblk2nd(nbl)) then
               do lll=lx,lvis-1
                  w(lxnm1+lll-lx) = w(lll)
               end do
            end if
      end if
c
c     collocate grid points on coarser levels; store off original
c     grid if mesh deforms
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 1010 m=1,ncg
         nbl = nbl+1
         if (myid.eq.mblk2nd(nbl)) then
            lxc = lw(10,nbl)
            lyc = lw(11,nbl)
            lzc = lw(12,nbl)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),110) nbl,ii2,jj2,kk2
  110       format(25h   creating coarser block,i6,
     .      24h of dimensions (I/J/K) :,3i6)
            call collx(w(lx),w(ly),w(lz),w(lxc),w(lyc),w(lzc),jdim,
     .                 kdim,idim,jj2,kk2,ii2)
            call lead(nbl,lw,lw2,maxbl)
            if (idefrm(nbl).gt.0) then
               do lll=lx,lvis-1
                  w(lxnm1+lll-lx) = w(lll)
               end do
            end if
         end if
 1010    continue
      end if
c        
      nblout = nblg(igrid)
#   ifdef FASTIO
      call writ_buffast(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl,16)
#   else
      call writ_buf(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl)
#   endif
c
 1000 continue
c
c     initialize slavept data for mesh deformation
c
      if (iunst.gt.1 .or. idef_ss .gt. 0) then
         call setslave(lw,lw2,w,mgwk,wk,nwork,maxbl,mxbli,maxgr,maxseg,
     .                 nbci0,nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,
     .                 ibcinfo,jbcinfo,kbcinfo,nblock,idefrm,
     .                 iadvance,nou,bou,nbuf,ibufdim,myid,myhost,
     .                 mycomm,mblk2nd,icsi,icsf,jcsi,jcsf,kcsi,
     .                 kcsf,islavept,nslave,nsegdfrm,idfrmseg,
     .                 maxsegdg,iwk,iwork,nmaster,ngrid,jskip,kskip,
     .                 iskip,nblg,levelg,lfgm,nblk,limblk,isva,nblelst,
     .                 nnodes,iskmax,jskmax,kskmax,nbli)
      end if
c
      if (idef_ss .gt. 0) then
c
c        read in new surface grid and deform mesh to reflect the
c        new surface
c
#if defined DIST_MPI
         if (myid .eq. myhost) then
#endif
         iunitr = 39
         open (unit=iunitr,file='newsurf.p3d',form='formatted',
     .         status='old')
         read(iunitr,*) ndum
         read(iunitr,*) (idum,jdum,kdum,ll=1,ndum)
#if defined DIST_MPI
         end if
#endif
c
c        read in the new surface data and temporarily store in 
c        the deltj/deltk/delti arrays
c
         do igrid = 1,ngrid
            nbl = nblg(igrid)
            if (idefrm(nbl) .eq. 1) then
               call lead(nbl,lw,lw2,maxbl)
#if defined DIST_MPI
               if (myid .eq. myhost) then
#endif
               do iseg=1,nsegdfrm(nbl)
                  call rsurf(maxbl,maxsegdg,idim,jdim,kdim,w(ldelti),
     .                       w(ldeltj),w(ldeltk),nbl,icsi,icsf,jcsi,
     .                       jcsf,kcsi,kcsf,iseg,nou,bou,nbuf,
     .                       ibufdim,iunitr)
               end do
#if defined DIST_MPI
               end if
c
               nd_dest = mblk2nd(nbl)
c
               nvals   = kdim*idim*3*2
               mytag = itag_x + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(w(ldeltj),nvals,MY_MPI_REAL,
     .                          nd_dest,mytag,mycomm,ierr)
               else if (mblk2nd(nbl).eq.myid) then
                  call MPI_Recv(w(ldeltj),nvals,MY_MPI_REAL,
     .                          myhost,mytag,mycomm,istat,ierr)
               end if
c
               nvals   = jdim*idim*3*2
               mytag = itag_y + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(w(ldeltk),nvals,MY_MPI_REAL,
     .                          nd_dest,mytag,mycomm,ierr)
               else if (mblk2nd(nbl).eq.myid) then
                  call MPI_Recv(w(ldeltk),nvals,MY_MPI_REAL,
     .                          myhost,mytag,mycomm,istat,ierr)
               end if
c
               nvals   = jdim*kdim*3*2
               mytag = itag_z + nbl
               if (myid .eq. myhost) then
                  call MPI_Send(w(ldelti),nvals,MY_MPI_REAL,
     .                          nd_dest,mytag,mycomm,ierr)
               else if (mblk2nd(nbl).eq.myid) then
                  call MPI_Recv(w(ldelti),nvals,MY_MPI_REAL,
     .                          myhost,mytag,mycomm,istat,ierr)
               end if
c
#endif
c
c              calculate delta displacement between new surface and
c              current surface
c 
#if defined DIST_MPI
               if (mblk2nd(nbl).eq.myid) then
#endif
               lwk1 = 1
               lwk2 = lwk1 + kdim*idim*2
               lwk3 = lwk2 + jdim*idim*2
               lwk4 = lwk3 + jdim*kdim*2
               if (nwork.lt.lwk4) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),444)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  444          format(37h not enough work space for subroutine,
     .                 8h getdelt)
c
c              wk(lwk1-lwk4) contain a flag that prevents points common
c              to multiple segments from being updated more than once
c
               do lll=lwk1,lwk4
                  wk(lll) = 1.
               end do
c
               do iseg=1,nsegdfrm(nbl)
                  call getdelt(maxbl,maxsegdg,idim,jdim,kdim,w(ldelti),
     .                        w(ldeltj),w(ldeltk),w(lx),w(ly),w(lz),
     .                        nbl,icsi,icsf,jcsi,jcsf,kcsi,kcsf,
     .                        iseg,nou,bou,nbuf,ibufdim,wk(lwk1),
     .                        wk(lwk2),wk(lwk3))
               end do
#if defined DIST_MPI
               end if
#endif
c
            end if
         end do
c
c        deform volume grid to fit new surface 
c
         iwk1 = 2*maxbl
         iwork1 = iwork - iwk1
         if (iwork1.lt.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...not enough'',
     .           ''integer work space for subroutine updatedg'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         do iii = 1,iwk1
            iwk(iii) = 0
         end do
         iseqr  = mseq
         iupdat = 1
         call updatedg(lw,lw2,w,mgwk,wk,nwork,iupdat,iseqr,maxbl,
     .                 maxgr,maxseg,nbci0,nbcj0,nbck0,nbcidim,
     .                 nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                 nblock,levelg,igridg,idefrm,ncgg,iadvance,nou,
     .                 bou,nbuf,ibufdim,myid,myhost,mycomm,mblk2nd,
     .                 utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                 xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                 rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                 kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                 aesrfdat,perturb,itrans,irotat,islavept,nslave,
     .                 iskip,jskip,kskip,xs,xxn,nsegdfrm,idfrmseg,
     .                 iaesurf,maxsegdg,iwk,nmaster,1,xorig,yorig,
     .                 zorig,xorgae0,yorgae0,zorgae0,icouple,
     .                 iwk(maxbl),nnodes,nblelst,iskmax,jskmax,kskmax,
     .                 ue)
c
c         collocate deformed grid points on coarser levels
c
          do igrid=1,ngrid
             nbl = nblg(igrid)
             call lead(nbl,lw,lw2,maxbl)
             ncg = ncgg(igrid)
             if (ncg.gt.0) then
                do m=1,ncg
                nbl = nbl+1
                if (myid.eq.mblk2nd(nbl)) then
                   lxc = lw(10,nbl)
                   lyc = lw(11,nbl)
                   lzc = lw(12,nbl)
                   nou(1) = min(nou(1)+1,ibufdim)
                   write(bou(nou(1),1),110) nbl,ii2,jj2,kk2
                   call collx(w(lx),w(ly),w(lz),w(lxc),w(lyc),w(lzc),
     .                        jdim,kdim,idim,jj2,kk2,ii2)
                   call lead(nbl,lw,lw2,maxbl)
                end if
                end do
             end if
          end do
c
      end if
c
c*******************************************************
c     compute grid metrics
c*******************************************************
c
      if (myid.eq.mblk2nd(1)) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
      end if
      iflag = 0
      icnt  = 0
      do 1200 igrid=1,ngrid
      iem  = iemg(igrid)
      nbl  = nblg(igrid)
      if (myid.eq.mblk2nd(nbl)) then
         call lead(nbl,lw,lw2,maxbl)
         if (nwork.lt.jdim*kdim*6+jdim*kdim*idim*5) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),300)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
  300    format(50h not enough work space for subroutine metric (must,
     .          16h increase mwork))
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),310) nbl,igrid
  310    format(28h computing metrics for block,i6,7h  (grid,i6,1h))
         call metric(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),w(lsk),
     .               w(lsi),wk(1),wk(jdim*kdim*6+1),nbl,iflag,icnt,
     .               nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .               nbckdim,ibcinfo,jbcinfo,kbcinfo,maxbl,maxseg,
     .               nou,bou,nbuf,ibufdim,myid,mblk2nd)

      end if
c
c     compute grid metrics on coarser levels
c
      ncg = ncgg(igrid)
      iem = iemg(igrid)
      if (ncg.gt.0) then
         isk1    = isklton
         isklton = 2
         do 1210 m=1,ncg
         nbl    = nbl+1
         if (mblk2nd(nbl).eq.myid) then
            call lead(nbl,lw,lw2,maxbl)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),320) nbl
  320       format(38h   computing metrics for coarser block,i6)
            call metric(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),w(lsk),
     .                  w(lsi),wk(1),wk(jdim*kdim*6+1),nbl,iflag,icnt,
     .                  nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                  nbckdim,ibcinfo,jbcinfo,kbcinfo,maxbl,maxseg,
     .                  nou,bou,nbuf,ibufdim,myid,mblk2nd)
         end if
 1210    continue
         isklton = isk1
      end if
c
      nblout = nblg(igrid)
#   ifdef FASTIO
      call writ_buffast(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl,17)
#   else
      call writ_buf(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl)
#   endif
c
#if defined DIST_MPI
         mytag   = itag_flag + nblout
         nd_scrc = mblk2nd(nblout)
c
         if (mblk2nd(nblout).eq.myid) then
            call MPI_Send (iflag, 1, MPI_INTEGER,
     .                     myhost, mytag, mycomm, ierr)
         else if (myid.eq.myhost) then
            call MPI_Recv (iflag1, 1, MPI_INTEGER,
     .                     nd_scrc,mytag,mycomm,istat,ierr)
            iflag = max(iflag, iflag1)
         end if
c
#endif
 1200 continue
c
      if (iflag .gt. 0 .and. myid.eq.myhost) then
        write(11,'('' Fatal error(s) uncovered in grid metrics'')')
        if (negvol .eq. 0) then
           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
      end if
c
c*******************************************************
c     compute cell volumes
c*******************************************************
c
      if (myid.eq.mblk2nd(1)) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
      end if
c
      do 1300 igrid=1,ngrid
      iem  = iemg(igrid)
      nbl  = nblg(igrid)
      if (myid.eq.mblk2nd(nbl)) then
         call lead(nbl,lw,lw2,maxbl)
         if (nwork.lt.jdim*kdim*15) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),400)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
  400    format(51h not enough work space for subroutine cellvol (must,
     .          16h increase mwork))
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),410)nbl,igrid
  410    format(33h computing cell volumes for block,i6,
     .           7h  (grid,i6,1h))
         call cellvol(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),
     .                w(lsk),w(lsi),w(lvol),wk,nou,bou,nbuf,ibufdim,
     .                myid,mblk2nd,maxbl,nbl,0,iflagv,imin,imax,jmin,
     .                jmax,kmin,kmax)
      end if
c
c     collocate cell volumes on coarser levels
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 1310 m=1,ncg
         nbl    = nbl+1
         if (myid.eq.mblk2nd(nbl)) then
            lvolc  = lw( 8,nbl)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),420) nbl
  420       format(43h   computing cell volumes for coarser block,i6)
            call collv(w(lvol),w(lvolc),jdim,kdim,idim,jj2,kk2,ii2)
            call lead(nbl,lw,lw2,maxbl)
         end if
 1310    continue
      end if
c
      nblout = nblg(igrid)
#   ifdef FASTIO
      call writ_buffast(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl,18)
#   else
      call writ_buf(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl)
#   endif
c
 1300 continue
c
c*******************************************************	  
c     read/set auxilliary data arrays for 2000 series
c     boundary conditions
c*******************************************************
#if defined DIST_MPI
c
c     set baseline tag values
c
      ioffset   = maxbl*maxseg
      itag_bci0 = 1
      itag_bcid = itag_bci0 + ioffset
      itag_bcj0 = itag_bcid + ioffset
      itag_bcjd = itag_bcj0 + ioffset
      itag_bck0 = itag_bcjd + ioffset
      itag_bckd = itag_bck0 + ioffset
      itag2k    = 0
#endif
c 
      if (myid.eq.mblk2nd(1)) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
      end if
c
      do 1400 igrid=1,ngrid
      iem  = iemg(igrid)
      nbl  = nblg(igrid)
#if defined DIST_MPI
      nd_dest = mblk2nd(nbl)
#endif
      call lead(nbl,lw,lw2,maxbl)
c
      do 1420 nseg = 1,nbci0(nbl)
      mm = 1
      if (abs(ibcinfo(nbl,nseg,1,mm)).ge.2000 .and.
     .    abs(ibcinfo(nbl,nseg,1,mm)).lt.3000) then
         ldata = lwdat(nbl,nseg,1)
         mdim = ibcinfo(nbl,nseg,3,mm)-ibcinfo(nbl,nseg,2,mm)
         ndim = ibcinfo(nbl,nseg,5,mm)-ibcinfo(nbl,nseg,4,mm)  
#if defined DIST_MPI
         nvals = mdim*ndim*12*2
         itag2k = itag2k + 1
         if (myid .eq. myhost) then
#endif
         filname = bcfiles(bcfilei(nbl,nseg,mm))
         do 1425 l=1,12
 1425    bcdat(l) = bcvali(nbl,nseg,l,mm)
         call readdat(w(ldata),mdim,ndim,filname,bcdat,nou,bou,nbuf,
     .                ibufdim,myid,mblk2nd,maxbl)
#if defined DIST_MPI
         mytag = itag_bci0 + itag2k
         call MPI_Send (w(ldata), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
         else if (mblk2nd(nbl).eq.myid) then
         mytag = itag_bci0 + itag2k
         call MPI_Recv (w(ldata), nvals, MY_MPI_REAL,
     .                  myhost,mytag,mycomm,istat,ierr)
         end if
#endif
      end if
 1420 continue
c
      do 1430 nseg = 1,nbcidim(nbl)
      mm = 2
      if (abs(ibcinfo(nbl,nseg,1,mm)).ge.2000 .and. 
     .    abs(ibcinfo(nbl,nseg,1,mm)).lt.3000) then
         ldata = lwdat(nbl,nseg,2)
         mdim = ibcinfo(nbl,nseg,3,mm)-ibcinfo(nbl,nseg,2,mm)
         ndim = ibcinfo(nbl,nseg,5,mm)-ibcinfo(nbl,nseg,4,mm)
#if defined DIST_MPI
         nvals = mdim*ndim*12*2
         itag2k = itag2k + 1
         if (myid .eq. myhost) then
#endif
         filname = bcfiles(bcfilei(nbl,nseg,mm))
         do 1435 l=1,12
 1435    bcdat(l) = bcvali(nbl,nseg,l,mm)
         call readdat(w(ldata),mdim,ndim,filname,bcdat,nou,bou,nbuf,
     .                ibufdim,myid,mblk2nd,maxbl)
#if defined DIST_MPI
         mytag = itag_bcid + itag2k
         call MPI_Send (w(ldata), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
         else if (mblk2nd(nbl).eq.myid) then
         mytag = itag_bcid + itag2k
         call MPI_Recv (w(ldata), nvals, MY_MPI_REAL,
     .                  myhost,mytag, mycomm, istat, ierr)
         end if
#endif
      end if
 1430 continue
c 
      do 1440 nseg = 1,nbcj0(nbl)
      mm = 1
      if (abs(jbcinfo(nbl,nseg,1,mm)).ge.2000 .and. 
     .    abs(jbcinfo(nbl,nseg,1,mm)).lt.3000) then
         ldata = lwdat(nbl,nseg,3)
         mdim = jbcinfo(nbl,nseg,5,mm)-jbcinfo(nbl,nseg,4,mm)
         ndim = jbcinfo(nbl,nseg,3,mm)-jbcinfo(nbl,nseg,2,mm)
#if defined DIST_MPI
         nvals = mdim*ndim*12*2
         itag2k = itag2k + 1
         if (myid .eq. myhost) then
#endif
         filname = bcfiles(bcfilej(nbl,nseg,mm))
         do 1445 l=1,12
 1445    bcdat(l) = bcvalj(nbl,nseg,l,mm)
         call readdat(w(ldata),mdim,ndim,filname,bcdat,nou,bou,nbuf,
     .                ibufdim,myid,mblk2nd,maxbl)
#if defined DIST_MPI
         mytag = itag_bcj0 + itag2k
         call MPI_Send (w(ldata), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
         else if (mblk2nd(nbl).eq.myid) then
         mytag = itag_bcj0 + itag2k
         call MPI_Recv (w(ldata), nvals, MY_MPI_REAL,
     .                  myhost,mytag, mycomm, istat, ierr)
         end if
#endif
      end if
 1440 continue
c
      do 1450 nseg = 1,nbcjdim(nbl)
      mm = 2
      if (abs(jbcinfo(nbl,nseg,1,mm)).ge.2000 .and. 
     .    abs(jbcinfo(nbl,nseg,1,mm)).lt.3000) then
         ldata = lwdat(nbl,nseg,4)
         mdim = jbcinfo(nbl,nseg,5,mm)-jbcinfo(nbl,nseg,4,mm)
         ndim = jbcinfo(nbl,nseg,3,mm)-jbcinfo(nbl,nseg,2,mm)
#if defined DIST_MPI
         nvals = mdim*ndim*12*2
         itag2k = itag2k + 1
         if (myid .eq. myhost) then
#endif
         filname = bcfiles(bcfilej(nbl,nseg,mm))
         do 1455 l=1,12
 1455    bcdat(l) = bcvalj(nbl,nseg,l,mm)
         call readdat(w(ldata),mdim,ndim,filname,bcdat,nou,bou,nbuf,
     .                ibufdim,myid,mblk2nd,maxbl)
#if defined DIST_MPI
         mytag = itag_bcjd + itag2k
         call MPI_Send (w(ldata), nvals, MY_MPI_REAL,
     .                 nd_dest, mytag, mycomm, ierr)
         else if (mblk2nd(nbl).eq.myid) then
         mytag = itag_bcjd + itag2k
         call MPI_Recv (w(ldata), nvals, MY_MPI_REAL,
     .                  myhost,mytag, mycomm, istat, ierr)
         end if
#endif
      end if
 1450 continue
c 
      do 1460 nseg = 1,nbck0(nbl)
      mm = 1
      if (abs(kbcinfo(nbl,nseg,1,mm)).ge.2000 .and. 
     .    abs(kbcinfo(nbl,nseg,1,mm)).lt.3000) then
         ldata = lwdat(nbl,nseg,5)
         mdim = kbcinfo(nbl,nseg,5,mm)-kbcinfo(nbl,nseg,4,mm)
         ndim = kbcinfo(nbl,nseg,3,mm)-kbcinfo(nbl,nseg,2,mm)
#if defined DIST_MPI
         nvals = mdim*ndim*12*2
         itag2k = itag2k + 1
         if (myid .eq. myhost) then
#endif
         filname = bcfiles(bcfilek(nbl,nseg,mm))
         do 1465 l=1,12
 1465    bcdat(l) = bcvalk(nbl,nseg,l,mm)
         call readdat(w(ldata),mdim,ndim,filname,bcdat,nou,bou,nbuf,
     .                ibufdim,myid,mblk2nd,maxbl)
#if defined DIST_MPI
         mytag = itag_bck0 + itag2k
         call MPI_Send (w(ldata), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
         else if (mblk2nd(nbl).eq.myid) then
         mytag = itag_bck0 + itag2k
         call MPI_Recv (w(ldata), nvals, MY_MPI_REAL,
     .                  myhost,mytag, mycomm, istat, ierr)
         end if
#endif
      end if
 1460 continue
c
      do 1470 nseg = 1,nbckdim(nbl)
      mm = 2
      if (abs(kbcinfo(nbl,nseg,1,mm)).ge.2000 .and. 
     .    abs(kbcinfo(nbl,nseg,1,mm)).lt.3000) then
         ldata = lwdat(nbl,nseg,6)
         mdim = kbcinfo(nbl,nseg,5,mm)-kbcinfo(nbl,nseg,4,mm)
         ndim = kbcinfo(nbl,nseg,3,mm)-kbcinfo(nbl,nseg,2,mm)
#if defined DIST_MPI
         nvals = mdim*ndim*12*2
         itag2k = itag2k + 1
         if (myid .eq. myhost) then
#endif
         filname = bcfiles(bcfilek(nbl,nseg,mm))
         do 1475 l=1,12
 1475    bcdat(l) = bcvalk(nbl,nseg,l,mm)
         call readdat(w(ldata),mdim,ndim,filname,bcdat,nou,bou,nbuf,
     .                ibufdim,myid,mblk2nd,maxbl)
#if defined DIST_MPI
         mytag = itag_bckd + itag2k
         call MPI_Send (w(ldata), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
         else if (mblk2nd(nbl).eq.myid) then
         mytag = itag_bckd + itag2k
         call MPI_Recv (w(ldata), nvals, MY_MPI_REAL,
     .                 myhost,mytag, mycomm, istat, ierr)
         end if
#endif
      end if
 1470 continue
c
c     collocate auxilliary data on coarser grids
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 1480 m=1,ncg
         nbl    = nbl+1
         if (mblk2nd(nbl).eq.myid) then
            lbcdjc = lw(42,nbl)
            lbcdkc = lw(43,nbl)
            lbcdic = lw(44,nbl)
c
            do 1490 nseg = 1,nbci0(nbl-1)
            mm = 1
            if (abs(ibcinfo(nbl,nseg,1,mm)).ge.2000 .and.
     .         abs(ibcinfo(nbl,nseg,1,mm)).lt.3000) then
               ldata = lwdat(nbl-1,nseg,1)
               mdim = ibcinfo(nbl-1,nseg,3,mm)-ibcinfo(nbl-1,nseg,2,mm)
               ndim = ibcinfo(nbl-1,nseg,5,mm)-ibcinfo(nbl-1,nseg,4,mm)
               ldatac = lwdat(nbl,nseg,1)
               mdimc = ibcinfo(nbl,nseg,3,mm)-ibcinfo(nbl,nseg,2,mm)
               ndimc = ibcinfo(nbl,nseg,5,mm)-ibcinfo(nbl,nseg,4,mm)
               call colldat(w(ldata),mdim,ndim,w(ldatac),mdimc,ndimc)
            end if
 1490       continue
c
            do 1500 nseg = 1,nbcidim(nbl)
            mm = 2
            if (abs(ibcinfo(nbl,nseg,1,mm)).ge.2000 .and.
     .         abs(ibcinfo(nbl,nseg,1,mm)).lt.3000) then
               ldata = lwdat(nbl-1,nseg,2)
               mdim = ibcinfo(nbl-1,nseg,3,mm)-ibcinfo(nbl-1,nseg,2,mm)
               ndim = ibcinfo(nbl-1,nseg,5,mm)-ibcinfo(nbl-1,nseg,4,mm)
               ldatac = lwdat(nbl,nseg,2)
               mdimc = ibcinfo(nbl,nseg,3,mm)-ibcinfo(nbl,nseg,2,mm)
               ndimc = ibcinfo(nbl,nseg,5,mm)-ibcinfo(nbl,nseg,4,mm)
               call colldat(w(ldata),mdim,ndim,w(ldatac),mdimc,ndimc)
            end if
 1500       continue
c
            do 1510 nseg = 1,nbcj0(nbl-1)
            mm = 1
            if (abs(jbcinfo(nbl,nseg,1,mm)).ge.2000 .and.
     .         abs(jbcinfo(nbl,nseg,1,mm)).lt.3000) then
               ldata = lwdat(nbl-1,nseg,3)
               mdim = jbcinfo(nbl-1,nseg,5,mm)-jbcinfo(nbl-1,nseg,4,mm)
               ndim = jbcinfo(nbl-1,nseg,3,mm)-jbcinfo(nbl-1,nseg,2,mm)
               ldatac = lwdat(nbl,nseg,3)
               mdimc = jbcinfo(nbl,nseg,5,mm)-jbcinfo(nbl,nseg,4,mm)
               ndimc = jbcinfo(nbl,nseg,3,mm)-jbcinfo(nbl,nseg,2,mm)
               call colldat(w(ldata),mdim,ndim,w(ldatac),mdimc,ndimc)
            end if
 1510       continue
c
            do 1520 nseg = 1,nbcjdim(nbl)
            mm = 2
            if (abs(jbcinfo(nbl,nseg,1,mm)).ge.2000 .and.
     .         abs(jbcinfo(nbl,nseg,1,mm)).lt.3000) then
               ldata = lwdat(nbl-1,nseg,4)
               mdim = jbcinfo(nbl-1,nseg,5,mm)-jbcinfo(nbl-1,nseg,4,mm)
               ndim = jbcinfo(nbl-1,nseg,3,mm)-jbcinfo(nbl-1,nseg,2,mm)
               ldatac = lwdat(nbl,nseg,4)
               mdimc = jbcinfo(nbl,nseg,5,mm)-jbcinfo(nbl,nseg,4,mm)
               ndimc = jbcinfo(nbl,nseg,3,mm)-jbcinfo(nbl,nseg,2,mm)
               call colldat(w(ldata),mdim,ndim,w(ldatac),mdimc,ndimc)
            end if
 1520       continue
c
            do 1530 nseg = 1,nbck0(nbl-1)
            mm = 1
            if (abs(kbcinfo(nbl,nseg,1,mm)).ge.2000 .and.
     .         abs(kbcinfo(nbl,nseg,1,mm)).lt.3000) then
               ldata = lwdat(nbl-1,nseg,5)
               mdim = kbcinfo(nbl-1,nseg,5,mm)-kbcinfo(nbl-1,nseg,4,mm)
               ndim = kbcinfo(nbl-1,nseg,3,mm)-kbcinfo(nbl-1,nseg,2,mm)
               ldatac = lwdat(nbl,nseg,5)
               mdimc = kbcinfo(nbl,nseg,5,mm)-kbcinfo(nbl,nseg,4,mm)
               ndimc = kbcinfo(nbl,nseg,3,mm)-kbcinfo(nbl,nseg,2,mm)
               call colldat(w(ldata),mdim,ndim,w(ldatac),mdimc,ndimc)
            end if
 1530       continue
c
            do 1540 nseg = 1,nbckdim(nbl)
            mm = 2
            if (abs(kbcinfo(nbl,nseg,1,mm)).ge.2000 .and.
     .         abs(kbcinfo(nbl,nseg,1,mm)).lt.3000) then
               ldata = lwdat(nbl-1,nseg,6)
               mdim = kbcinfo(nbl-1,nseg,5,mm)-kbcinfo(nbl-1,nseg,4,mm)
               ndim = kbcinfo(nbl-1,nseg,3,mm)-kbcinfo(nbl-1,nseg,2,mm)
               ldatac = lwdat(nbl,nseg,6)
               mdimc = kbcinfo(nbl,nseg,5,mm)-kbcinfo(nbl,nseg,4,mm)
               ndimc = kbcinfo(nbl,nseg,3,mm)-kbcinfo(nbl,nseg,2,mm)
               call colldat(w(ldata),mdim,ndim,w(ldatac),mdimc,ndimc)
            end if
 1540       continue
c
         end if
 1480    continue
      end if
c
 1400 continue
c
c*******************************************************
c     read overlapped grid interpolation data and set
c     iblank arrays
c*******************************************************
#if defined DIST_MPI
c
c     set baseline tag values
c
      ioffset     = maxbl
      itag_blk    = 1
      itag_lig    = itag_blk    + ioffset
      itag_lbg    = itag_lig    + ioffset
      itag_ibptsg = itag_lbg    + ioffset
      itag_iiptsg = itag_ibptsg + ioffset
      itag_jjbg   = itag_iiptsg + ioffset
      itag_kkbg   = itag_jjbg   + ioffset
      itag_iibg   = itag_kkbg   + ioffset
      itag_ibcg   = itag_iibg   + ioffset
      itag_jjig   = itag_ibcg   + ioffset
      itag_kkig   = itag_jjig   + ioffset
      itag_iiig   = itag_kkig   + ioffset
      itag_dxintg = itag_iiig   + ioffset
      itag_dyintg = itag_dxintg + ioffset
      itag_dzintg = itag_dyintg + ioffset
      itag_data   = itag_dzintg + ioffset
c
#endif
c
      do 1600 igrid=1,ngrid
      iem  = iemg(igrid)
      nbl  = nblg(igrid)
      call lead(nbl,lw,lw2,maxbl)
      if (iovrlp(nbl).eq.1) then
         if (myid .eq. mblk2nd(nbl)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),500)nbl,igrid
  500       format(38h reading overlap information for block,i6,
     .             7h  (grid,i6,1h))
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),101)nbl
  101       format(1x,37hreading overlap information for block,1x,i6)
         end if
         if (myid .eq. myhost) then
            if (nbl.gt.1) then
               if (iovrlp(nbl-1).eq.0) then
                  lig(nbl) = lig(nbl-1)
                  lbg(nbl) = lbg(nbl-1)
               end if
            end if
            call getibk(w(lblk),jdim,kdim,idim,nbl,intpts,nblpts,ibpnts,
     .                  iipnts,iitot,maxbl,iibg,kkbg,jjbg,ibcg,lig,
     .                  lbg,dxintg,dyintg,dzintg,iiig,jjig,kkig,
     .                  ibpntsg,iipntsg,myid,ibufdim,nbuf,bou,nou)
#if defined DIST_MPI
            nd_dest = mblk2nd(nbl)
            mytag = itag_data 
            call MPI_Send (nblpts, 1, MPI_INTEGER,
     .                     nd_dest, mytag, mycomm, ierr)
            mytag = itag_data + 1
            call MPI_Send (ibpnts, 1, MPI_INTEGER,
     .                     nd_dest, mytag, mycomm, ierr)
            mytag = itag_data + 2
            call MPI_Send (iipnts, 1, MPI_INTEGER,
     .                     nd_dest, mytag, mycomm, ierr)
            mytag = itag_data + 3
            call MPI_Send (intpts, 4, MPI_INTEGER,
     .                     nd_dest, mytag, mycomm, ierr)
#endif
         end if
         if (myid .eq. mblk2nd(nbl)) then
c
#if defined DIST_MPI
            mytag = itag_data
            call MPI_Recv (nblpts, 1, MPI_INTEGER,
     .                    myhost,mytag, mycomm, istat, ierr)
            mytag = itag_data + 1
            call MPI_Recv (ibpnts, 1, MPI_INTEGER,
     .                    myhost,mytag, mycomm, istat, ierr)
            mytag = itag_data + 2
            call MPI_Recv (iipnts, 1, MPI_INTEGER,
     .                    myhost,mytag, mycomm, istat, ierr)
            mytag = itag_data + 3
            call MPI_Recv (intpts, 4, MPI_INTEGER,
     .                    myhost,mytag, mycomm, istat, ierr)
c
#endif
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''  ibpnts,intpts,iipnts = '',6i8)')
     .                           ibpnts,intpts,iipnts
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*) '  counting iblank values....'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''  number of blank values .ne. 1 '',
     .                     ''are '',i8)') nblpts
         end if
c
      else
c
         if (myid .eq. mblk2nd(nbl)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),510)nbl,igrid
  510       format(38h setting blank array to 1.0 for block ,i6,
     .             7h  (grid,i6,1h))
         end if
         if (myid.eq.myhost) then
            call setblk(w(lblk),jdim,kdim,idim,nbl)
            if (nbl.gt.1) then
               lig(nbl) = lig(nbl-1)
               lbg(nbl) = lbg(nbl-1)
            end if
         end if
      end if
c
#if defined DIST_MPI
      nd_dest = mblk2nd(nbl)
      mytag = itag_blk + nbl
      nvals = jdim*kdim*idim
      if (myid .eq. myhost) then
         call MPI_Send (w(lblk), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (w(lblk), nvals, MY_MPI_REAL,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_lig + nbl
      if (myid .eq. myhost) then
         call MPI_Send (lig(nbl), 1, MPI_INTEGER, nd_dest,
     .                 mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (lig(nbl), 1, MPI_INTEGER, myhost,
     .                  mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_lbg + nbl
      if (myid .eq. myhost) then
         call MPI_Send (lbg(nbl), 1, MPI_INTEGER, nd_dest,
     .                 mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (lbg(nbl), 1, MPI_INTEGER, myhost,
     .                  mytag, mycomm, istat, ierr)
      end if
c
      if (iovrlp(nbl).eq.1) then
c
      mytag = itag_ibptsg + nbl
      if (myid .eq. myhost) then
         call MPI_Send (ibpntsg(nbl,1), 1, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (ibpntsg(nbl,1), 1, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
      if (myid .eq. myhost) then
         call MPI_Send (ibpntsg(nbl,2), 1, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (ibpntsg(nbl,2), 1, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
      if (myid .eq. myhost) then
         call MPI_Send (ibpntsg(nbl,3), 1, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (ibpntsg(nbl,3), 1, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
      if (myid .eq. myhost) then
         call MPI_Send (ibpntsg(nbl,4), 1, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (ibpntsg(nbl,4), 1, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_iiptsg + nbl
      if (myid .eq. myhost) then
         call MPI_Send (iipntsg(nbl), 1, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (iipntsg(nbl), 1, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_jjbg + nbl
      nvals = ibpntsg(nbl,1) + ibpntsg(nbl,2)
     .      + ibpntsg(nbl,3) + ibpntsg(nbl,4)
      if (myid .eq. myhost) then
         call MPI_Send (jjbg(lbg(nbl)), nvals, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (jjbg(lbg(nbl)), nvals, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_kkbg + nbl
      nvals = ibpntsg(nbl,1) + ibpntsg(nbl,2)
     .      + ibpntsg(nbl,3) + ibpntsg(nbl,4)
      if (myid .eq. myhost) then
         call MPI_Send (kkbg(lbg(nbl)), nvals, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (kkbg(lbg(nbl)), nvals, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_iibg + nbl
      nvals = ibpntsg(nbl,1) + ibpntsg(nbl,2)
     .      + ibpntsg(nbl,3) + ibpntsg(nbl,4)
      if (myid .eq. myhost) then
         call MPI_Send (iibg(lbg(nbl)), nvals, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (iibg(lbg(nbl)), nvals, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_ibcg + nbl
      nvals = ibpntsg(nbl,1) + ibpntsg(nbl,2)
     .      + ibpntsg(nbl,3) + ibpntsg(nbl,4)
      if (myid .eq. myhost) then
         call MPI_Send (ibcg(lbg(nbl)), nvals, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (ibcg(lbg(nbl)), nvals, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_jjig + nbl
      nvals = iipntsg(nbl)
      if (myid .eq. myhost) then
         call MPI_Send (jjig(lig(nbl)), nvals, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (jjig(lig(nbl)), nvals, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_kkig + nbl
      nvals = iipntsg(nbl)
      if (myid .eq. myhost) then
         call MPI_Send (kkig(lig(nbl)), nvals, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (kkig(lig(nbl)), nvals, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_iiig + nbl
      nvals = iipntsg(nbl)
      if (myid .eq. myhost) then
         call MPI_Send (iiig(lig(nbl)), nvals, MPI_INTEGER,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (iiig(lig(nbl)), nvals, MPI_INTEGER,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_dxintg + nbl
      nvals = iipntsg(nbl)
      if (myid .eq. myhost) then
         call MPI_Send (dxintg(lig(nbl)), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (dxintg(lig(nbl)), nvals, MY_MPI_REAL,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_dyintg + nbl
      nvals = iipntsg(nbl)
      if (myid .eq. myhost) then
         call MPI_Send (dyintg(lig(nbl)), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (dyintg(lig(nbl)), nvals, MY_MPI_REAL,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      mytag = itag_dzintg + nbl
      nvals = iipntsg(nbl)
      if (myid .eq. myhost) then
         call MPI_Send (dzintg(lig(nbl)), nvals, MY_MPI_REAL,
     .                  nd_dest, mytag, mycomm, ierr)
      else if (mblk2nd(nbl).eq.myid) then
         call MPI_Recv (dzintg(lig(nbl)), nvals, MY_MPI_REAL,
     .                  myhost, mytag, mycomm, istat, ierr)
      end if
c
      end if
#endif
c
c     set iblank for coarser grids. note: cfl3d currently
c     does not use interpolation data on coarser overlapped
c     grids
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 1610 m=1,ncg
         nbl    = nbl+1
         lblkc  = lw(18,nbl)
         if (iovrlp(nbl-1).eq.0) then
            lig(nbl) = lig(nbl-1)
            lbg(nbl) = lbg(nbl-1)
         end if
         if (mblk2nd(nbl).eq.myid) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),520) nbl
  520       format(41h   setting blank array to 1.0 for coarser,
     .             6h block,i6)
            call setblk(w(lblkc),jj2,kk2,ii2,nbl)
         end if
         call lead(nbl,lw,lw2,maxbl)
 1610    continue
      end if
c
      nblout = nblg(igrid)
#   ifdef FASTIO
      call writ_buffast(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl,19)
#   else
      call writ_buf(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl)
#   endif
c
 1600 continue
c
#if defined DIST_MPI
      call MPI_Bcast (lig, maxbl, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (iipntsg, maxbl, MPI_INTEGER, myhost,
     .                mycomm, ierr)
      mytag_qb(1) = 0
      do ll=2,8
         mytag_qb(ll) = mytag_qb(ll-1) + maxbl
      end do
      nreq = 1
#endif
c********************************************************
c     read patched-grid interpolation data from file.
c     file must contain data for coarser levels if needed
c********************************************************
c
      if (ninter.lt.0) then
         if (myid .eq. myhost) then
            call rpatch(maxbl,maxxe,intmax,nsub1,windex,ninter,iindex,
     .                  nblkpt,dthetxx,dthetyy,dthetzz,jdimg,kdimg,
     .                  idimg,nou,bou,nbuf,ibufdim,myid)
         end if
#if defined DIST_MPI
         call MPI_Bcast (ninter, 1, MPI_INTEGER, myhost,
     .                   mycomm, ierr)
         nvals = intmax*(2*nsub1+9)
         call MPI_Bcast (iindex, nvals, MPI_INTEGER, myhost,
     .                   mycomm, ierr)
         call MPI_Bcast (windex, maxxe*2, MY_MPI_REAL,
     .                            myhost,
     .                   mycomm, ierr)
         call MPI_Bcast (nblkpt, maxxe, MPI_INTEGER, myhost,
     .                   mycomm, ierr)
         nvals = intmax*nsub1
         call MPI_Bcast (dthetxx, nvals, MY_MPI_REAL,
     .                   myhost, mycomm, ierr)
         call MPI_Bcast (dthetyy, nvals, MY_MPI_REAL, 
     .                   myhost, mycomm, ierr)
         call MPI_Bcast (dthetzz, nvals, MY_MPI_REAL, 
     .                   myhost, mycomm, ierr)
#endif
      end if
c
#if defined DIST_MPI
c     initialize qiv array on the host for use in subroutine rrest
c
      if (myid.eq.myhost) then
         call init_mast
      end if
#endif
c
c*******************************************************
c     set initial conditions as freestream on all meshes
c     (subsequently overwritten if restarting). also set
c     boundary condition flags for use in finding smin.
c*******************************************************
c
      if (myid.eq.mblk2nd(1)) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
      end if
c
#   ifdef CMPLX
c     RSM not working in complex mode
#   else
c     Get some settings for RSM, and print some info:
      if (ivmx .eq. 72) then
        call kws_init(nou,bou,nbuf,ibufdim)
      end if
#   endif
      do 1100 igrid=1,ngrid
      nbl = nblg(igrid)-1
      ncg = ncgg(igrid)
      do 1110 n=1,ncg+1
      iem  = iemg(igrid)
      nbl = nbl+1
      if (myid.eq.mblk2nd(nbl)) then
         call lead(nbl,lw,lw2,maxbl)
         isksav=isklton
         isklton=0
         idum=0
         nttuse=max(ntt-1,1)
c        this first call to init is for safety in subsequent call
c        to bc; the call to bc will alter the initial bc data, so
c        then must call init again to make sure everything ends up
c        with freestream data
         call init(nbl,jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .             w(ltj0),w(ltk0),w(lti0),w(lvol),w(lvolj0),
     .             w(lvolk0),w(lvoli0),nummem,w(lx),w(ly),w(lz),
     .             nou,bou,nbuf,ibufdim,0)
c        call the noninertial initization to speed convergence
         if (noninflag.gt.0) then
           call initnonin(nbl,jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),
     .                    w(lqi0),w(lvol),
     .                    w(lvolj0),w(lvolk0),w(lvoli0),
     .                    w(lx),w(ly),w(lz))
         end if
         call bc(idum,nbl,lw,lw2,w,mgwk,wk,nwork,clw(nttuse),
     .           nou,bou,nbuf,ibufdim,maxbl,maxgr,maxseg,itrans,
     .           irotat,idefrm,igridg,nblg,nbci0,nbcj0,nbck0,
     .           nbcidim,nbcjdim,nbckdim,ibcinfo,jbcinfo,
     .           kbcinfo,bcfilei,bcfilej,bcfilek,lwdat,myid,
     .           idimg,jdimg,kdimg,bcfiles,mxbcfil,nummem)
         isklton=isksav
         if (n.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),200) nbl,igrid
  200       format(45h setting initial conditions as freestream for,
     .             6h block,i6,7h  (grid,i6,1h))
            if (noninflag.gt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),201)
  201          format(35h  correcting initial conditions for,
     .                28h NONINERTIAL reference frame)
            end if
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),210) nbl
  210       format(43h   setting initial conditions as freestream,
     .             18h for coarser block,i6)
            if (noninflag.gt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),211)
  211          format(37h    correcting initial conditions for,
     .                28h NONINERTIAL reference frame)
            end if
         end if
         call init(nbl,jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .             w(ltj0),w(ltk0),w(lti0),w(lvol),w(lvolj0),
     .             w(lvolk0),w(lvoli0),nummem,w(lx),w(ly),w(lz),
     .             nou,bou,nbuf,ibufdim,1)
c        call the noninertial initization to speed convergence
         if (noninflag.gt.0) then
           call initnonin(nbl,jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),
     .                    w(lqi0),w(lvol),
     .                    w(lvolj0),w(lvolk0),w(lvoli0),
     .                    w(lx),w(ly),w(lz))
         endif
      end if
 1110 continue
c
      nblout = nblg(igrid)
#   ifdef FASTIO
      call writ_buffast(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl,20)
#   else
      call writ_buf(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl)
#   endif
c
 1100 continue
c
c*******************************************************
c     initialize rigid body dynamics arrays/data
c*******************************************************
c
      if (irigb.gt.0.or.irbtrim.gt.0) then
            call init_rb(maxaes,maxbl,zorig,aesrfdat) 
      end if
      call init_trim
c
c*******************************************************
c     initialize aeroelastic arrays/data and read modal
c     surface definitions for aeroelastic case
c*******************************************************
#if defined DIST_MPI
c
c     set baseline tag values
c
      ioffset     = maxbl
      itag_aesj   = 1
      itag_aesk   = itag_aesj + ioffset
      itag_aesi   = itag_aesk + ioffset
#endif
c
      if (naesrf.gt.0) then
c
         if (myid.eq.myhost) then
c
c           initialize aeroelastic arrays and data
c
            call init_ae(nmds,maxaes,x0,xxn,wk,bmat,stm,stmi,gforcn,
     .                   gforcnm,freq,damp,gmass,aesrfdat,xs,gforcs)
         end if
#if defined DIST_MPI
c
         nnn = 2*nmds*maxaes
         call MPI_Bcast (xxn, nnn, MY_MPI_REAL, myhost,
     .                   mycomm, ierr)
         call MPI_Bcast (gforcn, nnn, MY_MPI_REAL, myhost,
     .                   mycomm, ierr)
         call MPI_Bcast (gforcnm, nnn, MY_MPI_REAL, myhost,
     .                   mycomm, ierr)
c
         nnn = 2*nmds*2*nmds*maxaes
         call MPI_Bcast (bmat, nnn, MY_MPI_REAL, myhost,
     .                   mycomm, ierr)
         call MPI_Bcast (stm, nnn, MY_MPI_REAL, myhost,
     .                   mycomm, ierr)
         call MPI_Bcast (stmi, nnn, MY_MPI_REAL, myhost,
     .                   mycomm, ierr)
#endif
c
c        read in modal surface definitions
c
         do iaes=1,naesrf
            nmodes = int(aesrfdat(5,iaes))
c
c           setup for plot3d output of modal surface definitions
c
            if (myid.eq.myhost) then
               idum    = 1
               jdum    = 1
               kdum    = 1
               mxaedum = 1
               nmddum  = 1
               nbldum  = 1
               lwk1    = 1
               lwk2    = lwk1 + idum*jdum*kdum
               lwk3    = lwk2 + idum*jdum*kdum
               lwk4    = lwk3 + kdum*idum*6*nmddum*mxaedum
               lwk5    = lwk4 + kdum*idum*6*nmddum*mxaedum
               lwk6    = lwk5 + kdum*idum*6*nmddum*mxaedum
               do nm=1,nmodes
                  iunitw  = 200+nm
                  call pltmode(nm,iaes,nblg,ngrid,maxgr,maxbl,
     .                         nsegdfrm,iaesurf,jbcinfo,kbcinfo,
     .                         ibcinfo,nbcj0,nbcjdim,nbck0,nbckdim,
     .                         nbci0,nbcidim,maxseg,maxsegdg,lw,lw2,
     .                         wk(lwk1),wk(lwk2),wk(lwk3),wk(lwk4),
     .                         wk(lwk5),wk(lwk6),jdum,kdum,idum,
     .                         mxaedum,nmddum,nbldum,iunitw,0)
               end do
            end if
c
            do igrid = 1,ngrid
              nbl = nblg(igrid)
              ncg = ncgg(igrid)
              call lead(nbl,lw,lw2,maxbl)
c             initialize modal data
              do ll=lxmdj,lvelg-1
                 w(ll) = 0.
              end do
            enddo
            do igrid = 1,ngrid
               nbl = nblg(igrid)
               ncg = ncgg(igrid)
               call lead(nbl,lw,lw2,maxbl)
               iaesrf = 0
               do is=1,nsegdfrm(nbl)
                  iaesrf = iaesrf + iaesurf(nbl,is)
               end do
#if defined DIST_MPI
c
c              need to send grid back to host for modal surface
c              plot3d data
c
               if (iaesrf.ne.0) then
c
                  nnn = jdim*kdim*idim
                  mytag = itag_x
                  nd_srce = mblk2nd(nbl)
                  if (myid.eq.mblk2nd(nbl)) then
                     call MPI_Send (w(lx), nnn, MY_MPI_REAL,
     .                              myhost, mytag, mycomm, ierr)
                  else if (myid.eq.myhost) then
                     call MPI_Recv (w(lx), nnn, MY_MPI_REAL,
     .                              nd_srce, mytag, mycomm, istat, ierr)
                  end if
c
                  nnn = jdim*kdim*idim
                  mytag = itag_y
                  nd_srce = mblk2nd(nbl)
                  if (myid.eq.mblk2nd(nbl)) then
                     call MPI_Send (w(ly), nnn, MY_MPI_REAL,
     .                              myhost, mytag, mycomm, ierr)
                  else if (myid.eq.myhost) then
                     call MPI_Recv (w(ly), nnn, MY_MPI_REAL,
     .                              nd_srce, mytag, mycomm, istat, ierr)
                  end if
c
                  nnn = jdim*kdim*idim
                  mytag = itag_z
                  nd_srce = mblk2nd(nbl)
                  if (myid.eq.mblk2nd(nbl)) then
                     call MPI_Send (w(lz), nnn, MY_MPI_REAL,
     .                              myhost, mytag, mycomm, ierr)
                  else if (myid.eq.myhost) then
                     call MPI_Recv (w(lz), nnn, MY_MPI_REAL,
     .                              nd_srce, mytag, mycomm, istat, ierr)
                  end if
c
               end if
#endif
               if (iaesrf.ne.0 .and. myid.eq.myhost) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),*)
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),605) nbl,idim,jdim,kdim
  605             format(32h reading modal surface for block,i6,
     .            24h of dimensions (I/J/K) :,3i6)
                  iunit = 33
                  jst = lxmdj
                  kst = lxmdk
                  ist = lxmdi
                  do nm = 1,nmodes
                   call modread(idim,jdim,kdim,nm,nbl,iunit,jbcinfo,
     .                          kbcinfo,ibcinfo,nbcj0,nbcjdim,nbck0,
     .                          nbckdim,nbci0,nbcidim,maxbl,maxseg,
     .                          nmds,w(jst),w(kst),w(ist),iaes,maxaes)
                   iunitw = 200+nm
                   call pltmode(nm,iaes,nblg,ngrid,maxgr,maxbl,
     .                          nsegdfrm,iaesurf,jbcinfo,kbcinfo,
     .                          ibcinfo,nbcj0,nbcjdim,nbck0,nbckdim,
     .                          nbci0,nbcidim,maxseg,maxsegdg,lw,lw2,
     .                          w(lx),w(ly),w(lz),w(jst),
     .                          w(kst),w(ist),jdim,kdim,idim,
     .                          maxaes,nmds,nbl,iunitw,1)
                  end do
c                 close(iunitw)
               end if

#if defined DIST_MPI
c
               nd_dest = mblk2nd(nbl)
c
               nnn = kdim*idim*6*nmodes
               ldj = lxmdj + nnn*(iaes-1)
               mytag = itag_aesj + nbl
               if (myid .eq. myhost) then
                  call MPI_Send (w(lxmdj), nnn, MY_MPI_REAL,
     .                           nd_dest, mytag, mycomm, ierr)
               else if (mblk2nd(nbl).eq.myid) then
                  call MPI_Recv (w(ldj), nnn, MY_MPI_REAL,
     .                           myhost, mytag, mycomm, istat, ierr)
               end if
c
               nnn = jdim*idim*6*nmodes
               ldk = lxmdk + nnn*(iaes-1)
               mytag = itag_aesk + nbl
               if (myid .eq. myhost) then
                  call MPI_Send (w(lxmdk), nnn, MY_MPI_REAL,
     .                           nd_dest, mytag, mycomm, ierr)
               else if (mblk2nd(nbl).eq.myid) then
                  call MPI_Recv (w(ldk), nnn, MY_MPI_REAL,
     .                           myhost, mytag, mycomm, istat, ierr)
               end if
c
               nnn = jdim*kdim*6*nmodes
               ldi = lxmdi + nnn*(iaes-1)
               mytag = itag_aesi + nbl
               if (myid .eq. myhost) then
                  call MPI_Send (w(lxmdi), nnn, MY_MPI_REAL,
     .                           nd_dest, mytag, mycomm, ierr)
               else if (mblk2nd(nbl).eq.myid) then
                  call MPI_Recv (w(ldi), nnn, MY_MPI_REAL,
     .                           myhost, mytag, mycomm, istat, ierr)
               end if
#endif
c
c              create modal surface definitions for coarser blocks
c
               if (ncg.gt.0) then
                  do m=1,ncg
                     nbl = nbl+1
                     if (iaesrf.ne.0 .and. myid.eq.mblk2nd(nbl)) then
                        lxmdjc = lw(52,nbl)
                        lxmdkc = lw(53,nbl)
                        lxmdic = lw(54,nbl)
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),610) nbl
  610                   format(25h   creating modal surface,
     .                  18h for coarser block,i6)
                        jst = lxmdj
                        kst = lxmdk
                        ist = lxmdi
                        call collmod(w(lxmdj),w(lxmdk),w(lxmdi),
     .                               w(lxmdjc),w(lxmdkc),w(lxmdic),
     .                               jdim,kdim,idim,jj2,kk2,ii2,
     .                               nm,nmds,iaes,maxaes)
                        call lead(nbl,lw,lw2,maxbl)
                     end if
                  end do
               end if
c
               nblout = nblg(igrid)
#   ifdef FASTIO
               call writ_buffast(nblout,11,nou,bou,nbuf,ibufdim,myhost,
     .                       myid,mycomm,mblk2nd,maxbl,21)
#   else
               call writ_buf(nblout,11,nou,bou,nbuf,ibufdim,myhost,
     .                       myid,mycomm,mblk2nd,maxbl)
#   endif
c
            end do
            do nm=1,nmodes
               close(200+nm)
            end do
         end do
      end if
c
c********************************************************
c     read restart data on mseq-1 level from finest level
c     i.e. the highest level computed on in the last run
c********************************************************
c
      if (irest.gt.0) then
         if (myid.eq.mblk2nd(1)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
         end if
c
         if (iteravg .eq. 2) then
         if (myid .eq. myhost) then
           read(97,end=1011,err=1011) ngx
           if (ngx .ne. ngrid) goto 1011
           read(97,end=1011,err=1011) (itest(igrid),jtest(igrid),
     .       ktest(igrid),igrid=1,ngx)

           if (ipertavg .ne. 0) then

              do igrid=1,ngrid
                 if (itest(igrid) .ne. idimg(nblg(igrid)+mseq-1) .or. 
     .                jtest(igrid) .ne. jdimg(nblg(igrid)+mseq-1) .or.
     .                ktest(igrid) .ne. kdimg(nblg(igrid)+mseq-1))
     .                goto 1011
              enddo
              read(98,end=1011,err=1011) ngx 
              if (ngx .ne. ngrid) goto 1011 
              read(98,end=1011,err=1011) (itest(igrid),jtest(igrid), 
     .             ktest(igrid),igrid=1,ngx) 
              do igrid=1,ngrid  
                 if (itest(igrid) .ne. idimg(nblg(igrid)+mseq-1) .or. 
     .                jtest(igrid) .ne. jdimg(nblg(igrid)+mseq-1) .or. 
     .                ktest(igrid) .ne. kdimg(nblg(igrid)+mseq-1)) 
     .                goto 1011 
c     
              enddo
           else
              do igrid=1,ngrid
                 if (itest(igrid) .ne. idimg(nblg(igrid)+mseq-1)-1 .or. 
     .                jtest(igrid) .ne. jdimg(nblg(igrid)+mseq-1)-1 .or.
     .                ktest(igrid) .ne. kdimg(nblg(igrid)+mseq-1)-1) 
     .                goto 1011
              enddo
           end if
           goto 1012
 1011      continue
           write(11,'(/,'' stopping... flag iteravg=2, but cannot'',
     .      '' read running-average Q file'')')
           write(11,'('' ... it either does not exist or it is an'',
     .      '' incorrect file'')')
           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
 1012      continue
         end if
         end if
c
         do 1700 igrid=1,ngrid
         iskipz = 0
         if (igrid.eq.1) iskipz = 1
         iem  = iemg(igrid)
         nbl  = nblg(igrid)
         inewg = inewgg(igrid)
         if (inewg.gt.0) go to 1701
         if (iem.ge.1 .and. mseq.gt.1) go to 1701
         if (icgns .eq. 1) then
         call lead(nbl,lw,lw2,maxbl)
         idima=idim
         jdima=jdim
         kdima=kdim
         end if
         nbl = nbl+(mseq-1)
         call lead(nbl,lw,lw2,maxbl)
         if (myid.eq.mblk2nd(nbl)) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),700) nbl,igrid
  700       format(31h reading restart file for block,i6,
     .      7h  (grid,i6,1h))
c
            if (iteravg .eq. 1) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),811) nbl
  811         format(35h running-average to be started from,
     .          18h scratch for block,i6)
            end if
            if (iteravg .eq. 2) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),812) nbl
  812         format(39h reading running-average file for block,i6)
            end if
c
               if (icgns .eq. 1) then
                 if(nwork .lt. (idima+1)*(jdima+1)*(kdima+1)) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),'('' not enough memory for'',
     .             '' cgns Q read.'')')
                 write(bou(nou(1),1),'('' nwork in wk='',i6,
     .             ''.  Needed = '',i6)')
     +              nwork,(idima+1)*(jdima+1)*(kdima+1)
                 call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                 end if
               else
                 idima=idim
                 jdima=jdim
                 kdima=kdim
               end if
         end if
         if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
            call rrest(nbl,jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .                 ncycmax,ntr,rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,
     .                 cmxw,cmyw,cmzw,
     .                 n_clcd, clcd, nblocks_clcd, blocks_clcd,
     $                 fmdotw,cftmomw,cftpw,cftvw,
     .                 cfttotw,rmstr,nneg,iskipz,w(lvis),
     .                 w(lxib),w(lsnk0),w(lsni0),w(lxkb),w(lnbl),
     .                 w(lcmuv),maxbl,mblk2nd,myid,myhost,mycomm,
     .                 nou,bou,nbuf,ibufdim,igrid,
     .                 wk,idima,jdima,kdima,w(lvj0),w(lvk0),w(lvi0),
     .                 w(ltj0),w(ltk0),w(lti0),w(lqavg),w(lq2avg),
     .                 nummem) 
            ntq     = ntt
            ncycchk = ncyctot+ntt
c
c           check to make sure ncycmax is large enough to store all
c           history info
c
            if (myid .eq. myhost) then
               if (ncycchk.gt.ncycmax) then
                  write(11,*)
                  write(11,*)' stopping...ncycmax too small to store',
     .                    ' prior convergence history AND current',
     .                    ' history '
                  write(11,*)' increase value of ncycmax to at LEAST ',
     .            ncycchk, ncyctot, ncycchk, ntt
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
            end if
c
         end if
 1701    continue
c
#   ifdef FASTIO
         call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,22)
#   else
         call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
c
 1700    continue
c
#if defined DIST_MPI
c
c        re-broadcast isminc to reflect change that may occur
c        in subroutine rrest
c
         call MPI_Bcast (isminc, 1, MPI_INTEGER, myhost, mycomm, ierr)
c
         if (myid.eq.myhost) then
           iflagg = 0
           if (icgns .ne. 1) then
             read(2,end=8988) iflagg
           else
#if defined CGNS
             if (irest.gt.0) then
               call getiflagg(iccg,ibase,nsoluse,iflagg)
             end if
#endif
           end if
 8988      continue
c
c          set iflagg = 0 if steady state (dt < 0)
c
           if (real(dt) .lt. 0.) iflagg = 0
c
         end if
c
         call MPI_Bcast (iflagg, 1, MPI_INTEGER, myhost, mycomm, ierr)
c
         if (iflagg .eq. 0) go to 8989
c
#else
         iflagg=0
         if (icgns .ne. 1) then
           read(2,end=8989) iflagg
         else
#if defined CGNS
           if (irest.gt.0) then
             call getiflagg(iccg,ibase,nsoluse,iflagg)
           end if
#endif
         end if
c
c        set iflagg = 0 if steady state (dt < 0)
c
         if (real(dt) .lt. 0.) iflagg = 0
         if (iflagg .eq. 0) go to 8989
c
#endif
c
c***********************************************************
c     read data for 2nd order time advancement
c***********************************************************
c
         if (iflagg .eq. 1 .or. iflagg .eq. 3) then
            do 1810 igrid=1,ngrid
            iem   = iemg(igrid)
            nbl   = nblg(igrid)
            inewg = inewgg(igrid)
            if (irest.eq.0 .or. inewg.gt.0) go to 1811
            if (iem.ge.1 .and. mseq.gt.1) go to 1811
           if (icgns .eq. 1) then
           call lead(nbl,lw,lw2,maxbl)
           idima=idim
           jdima=jdim
           kdima=kdim
           end if
            nbl = nbl+(mseq-1)
            if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
               call lead(nbl,lw,lw2,maxbl)
           if (icgns .eq. 1) then
             if(nwork .lt. (idima+1)*(jdima+1)*(kdima+1)) then
             write(11,'('' not enough memory for cgns '',
     +          ''unsteady Q read.'')')
             write(11,'('' nwork in wk='',i6,''.  Needed = '',i6)')
     +          nwork,(idima+1)*(jdima+1)*(kdima+1)
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
             end if
           else
             idima=idim
             jdima=jdim
             kdima=kdim
           end if
               call rrestg(nbl,igrid,jdim,kdim,idim,w(lx),
     .              w(ly),w(lz),w(lxnm2),w(lynm2),w(lznm2),
     .              w(ldeltj),w(ldeltk),w(ldelti),
     .              w(lqc0),0,iuns,utrans,vtrans,wtrans,
     .              omegax,omegay,omegaz,xorig,yorig,zorig,
     .              dxmx,dymx,dzmx,dthxmx,dthymx,
     .              dthzmx,thetax,thetay,thetaz,rfreqt,
     .              rfreqr,xorig0,yorig0,zorig0,time2,
     .              thetaxl,thetayl,thetazl,itrans,irotat,idefrm,
     .              utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .              xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .              rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .              kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .              aesrfdat,perturb,myhost,myid,mycomm,mblk2nd,
     .              maxbl,ibufdim,nbuf,bou,nou,nsegdfrm,idfrmseg,
     .              iaesurf,maxsegdg,wk,nwork,idima,jdima,kdima,
     .              w(lxib2),nummem)
            end if
 1811       continue
 1810       continue
         else
            if (myid.eq.mblk2nd(1)) then
               if (time .ne. 0.) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),'('' WARNING:  You are '',
     .            ''restarting a 2nd order in time solution with no'')')
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),'(''           qc0 data saved '',
     .            ''in restart. You may see a glitch in the'')')
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),
     .            '(''           restart as a result.'')')
               end if
            end if
         end if
c
c***********************************************************
c     read dynamic mesh data on restart level and reconstruct
c     the last grid position before restart
c***********************************************************
#if defined DIST_MPI
c
c        set baseline tag values
c
         ioffset   = maxbl
         itag_ns   = 1
         itag_wk   = itag_ns + ioffset
#endif
c
         if (real(dt).gt.0. .and. iunst .gt. 0) then
            if (iflagg .eq. 2 .or. iflagg .eq. 3) then
               if (ndgrd .ne. 0 .and. iunst .gt. 1) then
                 if(myid.eq.myhost) then
                   open(98,file='dgplot3d.bin',form='formatted'
     .                 ,status='old')
                   read(98,*) ngridt
                   read(98,*) (idimt,jdimt,kdimt,igrid=1,ngrid)
                 end if
               end if
               do 1800 igrid=1,ngrid
               iem    = iemg(igrid)
               nbl    = nblg(igrid)
               inewg  = inewgg(igrid)
               if (irest.eq.0 .or. inewg.gt.0) go to 1801
               if (iem.ge.1 .and. mseq.gt.1) go to 1801
            if (icgns .eq. 1) then
            call lead(nbl,lw,lw2,maxbl)
            idima=idim
            jdima=jdim
            kdima=kdim
            end if
               nbl = nbl+(mseq-1)
               call lead(nbl,lw,lw2,maxbl)
           if (icgns .eq. 1) then
             if(nwork .lt. (idima+1)*(jdima+1)*(kdima+1)) then
             write(901,'('' not enough memory for cgns '',
     +          ''unsteady Q read.'')')
             write(901,'('' nwork in wk='',i6,''.  Needed = '',i6)')
     +          nwork,(idima+1)*(jdima+1)*(kdima+1)
             call termn8(myid,-1,ibufdim,nbuf,bou,nou)
             end if
           else
             idima=idim
             jdima=jdim
             kdima=kdim
           end if
               if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
                  iunsn = max(itrans(nbl),irotat(nbl),idefrm(nbl))
                  call rrestg(nbl,igrid,jdim,kdim,idim,w(lx),
     .                 w(ly),w(lz),w(lxnm2),w(lynm2),w(lznm2),
     .                 w(ldeltj),w(ldeltk),w(ldelti),
     .                 w(lqc0),1,iuns,utrans,vtrans,wtrans,
     .                 omegax,omegay,omegaz,xorig,yorig,zorig,
     .                 dxmx,dymx,dzmx,dthxmx,dthymx,
     .                 dthzmx,thetax,thetay,thetaz,rfreqt,
     .                 rfreqr,xorig0,yorig0,zorig0,time2,
     .                 thetaxl,thetayl,thetazl,itrans,irotat,idefrm,
     .                 utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                 xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                 rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                 kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                 aesrfdat,perturb,myhost,myid,mycomm,mblk2nd,
     .                 maxbl,ibufdim,nbuf,bou,nou,nsegdfrm,idfrmseg,
     .                 iaesurf,maxsegdg,wk,nwork,idima,jdima,kdima,
     .                 w(lxib2),nummem)
                  if (myid.eq.myhost) then
                     if (iuns .ne. iunsn) then
                        write(11,'('' Stopping: cannot alter type'',
     .                  '' of grid motion between restarts'')')
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
                  end if
               end if
 1801          continue
c
#   ifdef FASTIO
               call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,
     .                       myid,mycomm,mblk2nd,maxbl,23)
#   else
               call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                       mycomm,mblk2nd,maxbl)
#   endif
c
 1800          continue
c
c              reconstruct last grid position for purely rigid grid case
c
               if (iunst.gt.0) then
                  do igrid=1,ngrid
                     nbl = nblg(igrid)
                     nbl = nbl+(mseq-1)
                     call lead(nbl,lw,lw2,maxbl)
                     iuns = max(itrans(nbl),irotat(nbl))
                     if (mblk2nd(nbl).eq.myid) then
                        if (idefrm(nbl).eq.0 .and. iuns.gt.0) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),810) nbl
  810                      format(25h reconstructing last grid,
     .                            29h position from data for block,
     .                            i6)
                           call grdmove(nbl,jdim,kdim,idim,w(lx),w(ly),
     .                          w(lz),xorig0(nbl),yorig0(nbl),
     .                          zorig0(nbl),xorig(nbl),yorig(nbl),
     .                          zorig(nbl),thetax(nbl),thetay(nbl),
     .                          thetaz(nbl))
                        end if
                     end if
#   ifdef FASTIO
                     call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,
     .                             myhost,myid,mycomm,mblk2nd,maxbl,24)
#   else
                     call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,
     .                             myid,mycomm,mblk2nd,maxbl)
#   endif
                  end do
               end if
c 
            end if
         end if
c
         goto 8990
c
 8989    continue
c
         if (myid.eq.mblk2nd(1)) then
            if (real(dt).gt.0. and. abs(ita) .eq. 2) then
              if (time .ne. 0.) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),'('' WARNING:  You are '',
     .           ''restarting a 2nd order in time solution with no'')')
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),'(''           qc0 data saved '',
     .           ''in restart. You may see a glitch in the'')')
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),
     .           '(''           restart as a result.'')')
              end if
            end if
         end if
c
 8990    continue
c
c        read modal data if aeroelastic
c
         if (myid.eq.myhost) then
            if (naesrf.gt.0) then
            if (icgns .ne. 1) then
               read(2,end=8995) (timekeep(nn),nn=1,ntt)
               do iaes=1,naesrf
                  nmodes = aesrfdat(5,iaes)
                  read(2,end=8995) (xxn(n,iaes),n=1,2*nmodes)
                  read(2,end=8995) (gforcn(n,iaes),n=1,2*nmodes)
                  read(2,end=8995) (gforcnm(n,iaes),n=1,2*nmodes)
                  read(2,end=8995) (((aehist(nn,ll,n,iaes),nn=1,ntt),
     .                                ll=1,3),n=1,nmodes)
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),'('' reading generalized force/'',
     .            ''displacement data for aeroelastic surface '',i6)')
     .            iaes
c
c                 if x0 data is non-zero in input file, overwrite xxn with x0
c                 this allows a static ae solution to be perturbed for dynamic
c                 response
c
                  do n=1,2*nmodes
                     if (abs(real(x0(n,iaes))) .ne. 0.) then
                        xxn(n,iaes) = x0(n,iaes)
                     end if
                  end do
               end do
            else
#if defined CGNS
               maxnum=2*nmds*naesrf
               maxnum2=ntt*3*nmds*naesrf
               if ((maxnum+maxnum2) .gt. nwork) then
                 write(901,'('' not enough room in wk for'',
     .             '' writing data in raeromode'')')
                 call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'('' reading generalized force/'',
     .         ''displacement data from cgns file for'')')
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''    all aeroelastic surfaces'')')
               call raeromode(iccg,ibase,ncycmax,nmds,maxaes,
     .         ntt,naesrf,wk,maxnum,wk2,maxnum2,
     .         timekeep,xxn,gforcn,gforcnm,aehist)
#endif
            end if
            end if
 8995       continue
         end if
#if defined DIST_MPI
c
         nnn = 2*nmds*maxaes
         call MPI_Bcast (xxn, nnn, MY_MPI_REAL, myhost,
     .                   mycomm, ierr)
         call MPI_Bcast (gforcn, nnn, MY_MPI_REAL, myhost,
     .                   mycomm, ierr)
         call MPI_Bcast (gforcnm, nnn, MY_MPI_REAL, myhost,
     .                   mycomm, ierr)
#endif
c
      end if
      if (myid.eq.myhost) then
      if (icgns .ne. 1) then
            close (2)
            close (102)
            if(ndgrd .ne. 0 .and. iunst .gt. 1) then
              close (98)
            end if
      else
#if defined CGNS
c   Close data base
      write(901,'('' Closing CGNS database from setup'')')
      call cg_close_f(iccg, ier)
      if (ier .ne. 0) then
         call cg_error_print_f
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
#endif
      end if
      end if
c
c********************************************************
c     compute directed distance function for Baldwin-
c     Lomax turbulence model
c********************************************************
c
      if (myid.eq.mblk2nd(1)) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
      end if
c
      do 1900 igrid=1,ngrid
      nbl   = nblg(igrid)
      call lead(nbl,lw,lw2,maxbl)
      if (myid.eq.mblk2nd(nbl)) then
         if (ivisc(1).eq.2 .or. ivisc(2).eq.2 .or. ivisc(3).eq.2 .or.
     .       ivisc(1).eq.3 .or. ivisc(2).eq.3 .or. ivisc(3).eq.3) then
            inewg = inewgg(igrid)
            numchk = (jdim-1)*(kdim-1)*(idim-1)*4
            if (nwork.lt.numchk) then
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),900)
              call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
  900       format(48h not enough work space for subroutine dird (must,
     .             16h increase mwork))
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),901) nbl,igrid
  901       format(39h computing directed distances for block,i6,
     .             7h  (grid,i6,1h))
            call dird(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),w(lsk),
     .                w(lsi),w(lsni0),w(lsnk0),w(lsni0),w(lnbl),
     .                w(lxkb),w(lnbl),wk,ivisc,nou,bou,nbuf,ibufdim)
c
c           compute directed distances on coarser levels
c
            ncg = ncgg(igrid)
            iem = iemg(igrid)
            if (ncg.gt.0) then
               do 1910 m=1,ncg
               nbl    = nbl+1
               call lead(nbl,lw,lw2,maxbl)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),910) nbl
  910          format(43h   computing directed distances for coarser,
     .                6h block,i6)
               call dird(jdim,kdim,idim,w(lx),w(ly),w(lz),w(lsj),w(lsk),
     .                   w(lsi),w(lsni0),w(lsnk0),w(lsni0),w(lnbl),
     .                   w(lxkb),w(lnbl),wk,ivisc,nou,bou,nbuf,ibufdim)
1910           continue
            end if
         end if
      end if
c
      nblout = nblg(igrid)
#   ifdef FASTIO
      call writ_buffast(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl,25)
#   else
      call writ_buf(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl)
#   endif
c
1900  continue
c
c********************************************************
c     minimum distance function for field equation
c     turbulence models
c********************************************************
c
      if (isminc .ge. 1) then
      if (isminc .eq. 1 .or. isminc .eq. 2) then
c
c        avoid computing min dist function on grid levels 
c        where it's not needed
c
         if(mseq.eq.2 .and. ncyc1(2).eq.0) then
           j1=1
         else if(mseq.eq.3 .and. ncyc1(2).eq.0) then
           j1=2
         else if(mseq.eq.3 .and. ncyc1(3).eq.0) then
           j1=1
         else if(mseq.eq.4 .and. ncyc1(2).eq.0) then
           j1=3
         else if(mseq.eq.4 .and. ncyc1(3).eq.0) then
           j1=2
         else if(mseq.eq.4 .and. ncyc1(4).eq.0) then
           j1=1
         else if(mseq.eq.5 .and. ncyc1(2).eq.0) then
           j1=4
         else if(mseq.eq.5 .and. ncyc1(3).eq.0) then
           j1=3
         else if(mseq.eq.5 .and. ncyc1(4).eq.0) then
           j1=2
         else if(mseq.eq.5 .and. ncyc1(5).eq.0) then
           j1=1
         else
           j1=0
         end if
c
c        for recursive box algorithm, need some bc info
c        on the finest level regardless of the top level
c        actually computed on via mseq, so temporarily set j1=0
c
         j1sav = j1
         j1 = 0
c
c        Initialize smin to large value over all blocks
c
         itry=0
         do 2000 n=1,ngrid
         iset1 = nblg(n) + j1
         iset2 = nblg(n) + ncgg(n)
         do 2000 nbl=iset1,iset2
         call lead(nbl,lw,lw2,maxbl)
         if (ivisc(1).ge.4 .or. ivisc(2).ge.4 .or. 
     .      ivisc(3).ge.4)then
            itry=itry+1
            if (myid.eq.mblk2nd(nbl)) then
               ntot=(jdim-1)*(kdim-1)
               izz=0
               do 2010 i=1,idim-1
               do 2010 m=1,ntot
               izz=izz+1
               w(lsnk0+izz-1)=1.e20
 2010          continue
            end if
         end if
 2000    continue
c  If no walls at all, don't do min distance calculation
         do 5858 n=1,ngrid
            nbl=nblg(n)
            do iset=1,nbci0(nbl)
              if (abs(ibcinfo(nbl,iset,1,1)) .eq. 2004 .or.
     .            abs(ibcinfo(nbl,iset,1,1)) .eq. 2014 .or.
     .            abs(ibcinfo(nbl,iset,1,1)) .eq. 2024 .or.
     .            abs(ibcinfo(nbl,iset,1,1)) .eq. 2034 .or.
     .            abs(ibcinfo(nbl,iset,1,1)) .eq. 2016) goto 5860
            enddo
            do iset=1,nbcidim(nbl)
              if (abs(ibcinfo(nbl,iset,1,2)) .eq. 2004 .or.
     .            abs(ibcinfo(nbl,iset,1,2)) .eq. 2014 .or.
     .            abs(ibcinfo(nbl,iset,1,2)) .eq. 2024 .or.
     .            abs(ibcinfo(nbl,iset,1,2)) .eq. 2034 .or.
     .            abs(ibcinfo(nbl,iset,1,2)) .eq. 2016) goto 5860
            enddo
            do iset=1,nbcj0(nbl)
              if (abs(jbcinfo(nbl,iset,1,1)) .eq. 2004 .or.
     .            abs(jbcinfo(nbl,iset,1,1)) .eq. 2014 .or.
     .            abs(jbcinfo(nbl,iset,1,1)) .eq. 2024 .or.
     .            abs(jbcinfo(nbl,iset,1,1)) .eq. 2034 .or.
     .            abs(jbcinfo(nbl,iset,1,1)) .eq. 2016) goto 5860
            enddo
            do iset=1,nbcjdim(nbl)
              if (abs(jbcinfo(nbl,iset,1,2)) .eq. 2004 .or.
     .            abs(jbcinfo(nbl,iset,1,2)) .eq. 2014 .or.
     .            abs(jbcinfo(nbl,iset,1,2)) .eq. 2024 .or.
     .            abs(jbcinfo(nbl,iset,1,2)) .eq. 2034 .or.
     .            abs(jbcinfo(nbl,iset,1,2)) .eq. 2016) goto 5860
            enddo
            do iset=1,nbck0(nbl)
              if (abs(kbcinfo(nbl,iset,1,1)) .eq. 2004 .or.
     .            abs(kbcinfo(nbl,iset,1,1)) .eq. 2014 .or.
     .            abs(kbcinfo(nbl,iset,1,1)) .eq. 2024 .or.
     .            abs(kbcinfo(nbl,iset,1,1)) .eq. 2034 .or.
     .            abs(kbcinfo(nbl,iset,1,1)) .eq. 2016) goto 5860
            enddo
            do iset=1,nbckdim(nbl)
              if (abs(kbcinfo(nbl,iset,1,2)) .eq. 2004 .or.
     .            abs(kbcinfo(nbl,iset,1,2)) .eq. 2014 .or.
     .            abs(kbcinfo(nbl,iset,1,2)) .eq. 2024 .or.
     .            abs(kbcinfo(nbl,iset,1,2)) .eq. 2034 .or.
     .            abs(kbcinfo(nbl,iset,1,2)) .eq. 2016) goto 5860
            enddo
 5858    continue
         itry=0
 5860    continue
c
         if (itry.ne.0) then
c
c              calculate minimum distance via recursive-box algorithm 
c
               j1 = j1sav
c
c              check storage availability
               nworkf = 0
               nworki = 0
               if (myid.eq.myhost) then
                  call cntsurf(ns2004,maxbl,maxgr,maxseg,ngrid,nblg,
     .                         nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                         nbckdim,ibcinfo,jbcinfo,kbcinfo,2004)
                  call cntsurf(ns2014,maxbl,maxgr,maxseg,ngrid,nblg,
     .                         nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                         nbckdim,ibcinfo,jbcinfo,kbcinfo,2014)
                  call cntsurf(ns2024,maxbl,maxgr,maxseg,ngrid,nblg,
     .                         nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                         nbckdim,ibcinfo,jbcinfo,kbcinfo,2024)
                  call cntsurf(ns2034,maxbl,maxgr,maxseg,ngrid,nblg,
     .                         nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                         nbckdim,ibcinfo,jbcinfo,kbcinfo,2034)
                  call cntsurf(ns2016,maxbl,maxgr,maxseg,ngrid,nblg,
     .                         nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                         nbckdim,ibcinfo,jbcinfo,kbcinfo,2016)
                  nsurf=ns2004+ns2014+ns2024+ns2034+ns2016
               end if
#if defined DIST_MPI
               call MPI_Bcast (nsurf, 1, MPI_INTEGER, myhost,
     .                         mycomm, ierr)
#endif
               minbox=sqrt(float(nsurf))
               minbox=max(minbox,50)
               nbb = 3*nsurf/minbox
               do 2055 igrid = 1,ngrid
               nbl = nblg(igrid)
               ntempf =  9*nsurf + 7*nbb
     .                +  4*jdimg(nbl)*kdimg(nbl)*idimg(nbl)
               if (ivmx.eq.4 .or. ivmx.eq.25) then
                  ntempi = 15*nsurf + 2*nbb
     .                   +    jdimg(nbl)*kdimg(nbl)*idimg(nbl)
               else
                  ntempi = 11*nsurf + 2*nbb
     .                   +    jdimg(nbl)*kdimg(nbl)*idimg(nbl)
               end if
               nworkf = max(nworkf,ntempf)
               nworki = max(nworki,ntempi)
 2055          continue
c              additional memory allocation for temporary storage of 
c              smin at all grid points on all levels
               nadd = 0
               do 2065 nbl= 1,nblock
               if (mblk2nd(nbl).eq.myid) then
                  nadd = nadd + jdimg(nbl)*kdimg(nbl)*idimg(nbl)
               end if
 2065          continue
               nworkf = nworkf  + nadd
c
c              nworkxs, nworkixs needed to store surface data from
c              subroutines getpts and getptsbb; nworkxs and nworkixs
c              must be compatable with dimension statements for arrays
c              xs and ixs in subroutines getpts, getptsbb, collect_surf,
c              and collect_surfbb
c
               nworkxs  = 4*nsurf
               if (ivmx.eq.4 .or. ivmx.eq.25) then
                  nworkixs = 4*nsurf
               else
                  nworkixs = 0
               end if
c
c              misc. additional integer work space for findmin_new
c
               nwklsmin = maxbl
               nwkireq  = maxbl*maxseg*6
c
               nroomf = nwork    - nworkf   - nworkxs
               nroomi = iwork    - nworki   - nworkixs
     .                - nwklsmin - 2*nwkireq
#if defined DIST_MPI
               if (myid.ne.myhost) then
#endif
               if (nroomf.lt.0) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),753)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  753          format(26h not enough work space for,
     .                23h subroutine findmin_new)
#if defined DIST_MPI
               end if
               if (myid.ne.myhost) then
#endif
               if (nroomi.lt.0) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),754)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
  754          format(34h not enough integer work space for,
     .                23h subroutine findmin_new)
c
#if defined DIST_MPI
               end if
               if (myid.eq.mblk2nd(1)) then
#endif
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),992)
  992          format(32h setting up for minimum distance,
     .         42h calculation using recursive-box algorithm)
#if defined DIST_MPI
               end if
               if (myid.ne.myhost) then
#endif
               lwk1 = 1
               lwk2 = lwk1 + nworkf
               iwk1 = 1
               iwk2 = iwk1 + nworki
               iwk3 = iwk2 + nworkixs
               iwk4 = iwk3 + nwklsmin
               iwk5 = iwk4 + nwkireq
               iwk6 = iwk5 + nwkireq
               do iii = 1,iwk6
                  iwk(iii) = 0
               end do
               call findmin_new(lw,lw2,w,mgwk,wk(lwk1),nworkf,
     .                          iwk(iwk1),nworki,nsurf,j1,
     .                          wk(lwk2),iwk(iwk2),iwk(iwk3),
     .                          iwk(iwk4),iwk(iwk5),ngrid,ncgg,
     .                          nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                          nbckdim,jbcinfo,kbcinfo,ibcinfo,
     .                          nblg,nou,bou,nbuf,ibufdim,maxbl,
     .                          maxgr,maxseg,mblk2nd)
#if defined DIST_MPI
               end if
#endif
c
               do nbl=1,nblock
#   ifdef FASTIO
                  call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,
     .                          myid,mycomm,mblk2nd,maxbl,26)
#   else
                  call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                          mycomm,mblk2nd,maxbl)
#   endif
               end do
c
         end if
      end if
c
c        need to initialize turb quantities if there is turbulence
c        on any block. (NOTE: For the case of restarting in a
c        mesh-sequencing fashion from a coarser grid solution,
c        i.e. mseq > 1, we do NOT want to reinitialize these
c        quantities, since they will be interpolated from the 
c        coarser mesh solution saved in the restart file
c
         if (isminc .eq. 1 .or. isminc .eq. 3) then
            do 2070 nbl=1,nblock
            if (mblk2nd(nbl).eq.myid) then
               call lead(nbl,lw,lw2,maxbl)
               if (nbl.eq.1) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),260)
               end if
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),265) nbl
               call initvist(nbl,jdim,kdim,idim,w(lvis),w(lxib),
     .                       w(lsnk0),w(lcmuv),nummem,w(lx),w(ly),w(lz))
            endif
c
#   ifdef FASTIO
            call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                    mycomm,mblk2nd,maxbl,27)
#   else
            call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                    mycomm,mblk2nd,maxbl)
#   endif
c
 2070       continue
c
  260       format(1h )
  265       format(44h setting turbulent initial conditions, block,i6)
         end if
c
      end if
c   make smin negative in regions where want laminar flow
c   (field equation models only)
      if (ivmx .ge. 4) then
      do n=1,ngrid
        llev=0
        do nbl=nblg(n),nblg(n)+ncgg(n)
          llev=llev+1
          if (mblk2nd(nbl).eq.myid) then
            call lead(nbl,lw,lw2,maxbl)
c     J-dir:
            do nseg=1,nbcj0(nbl)
              if(abs(jbcinfo(nbl,nseg,1,1)).eq.2014) then
                ibeg=jbcinfo(nbl,nseg,2,1)
                iend=jbcinfo(nbl,nseg,3,1)
                jbeg=1
                jend=1
                kbeg=jbcinfo(nbl,nseg,4,1)
                kend=jbcinfo(nbl,nseg,5,1)
                mdim=kend-kbeg
                ndim=iend-ibeg
                ldata=lwdat(nbl,nseg,3)
                call lamfix(jdim,kdim,idim,w(lsnk0),mdim,ndim,w(ldata),
     .                 ibeg,iend,jbeg,jend,kbeg,kend,3,llev)
              end if
            enddo
            do nseg=1,nbcjdim(nbl)
              if(abs(jbcinfo(nbl,nseg,1,2)).eq.2014) then
                ibeg=jbcinfo(nbl,nseg,2,2)
                iend=jbcinfo(nbl,nseg,3,2)
                jbeg=jdim
                jend=jdim
                kbeg=jbcinfo(nbl,nseg,4,2)
                kend=jbcinfo(nbl,nseg,5,2)
                mdim=kend-kbeg
                ndim=iend-ibeg
                ldata=lwdat(nbl,nseg,4)
                call lamfix(jdim,kdim,idim,w(lsnk0),mdim,ndim,w(ldata),
     .                 ibeg,iend,jbeg,jend,kbeg,kend,4,llev)
              end if
            enddo
c     K-dir:
            do nseg=1,nbck0(nbl)
              if(abs(kbcinfo(nbl,nseg,1,1)).eq.2014) then
                ibeg=kbcinfo(nbl,nseg,2,1)
                iend=kbcinfo(nbl,nseg,3,1)
                jbeg=kbcinfo(nbl,nseg,4,1)
                jend=kbcinfo(nbl,nseg,5,1)
                kbeg=1
                kend=1
                mdim=jend-jbeg
                ndim=iend-ibeg
                ldata=lwdat(nbl,nseg,5)
                call lamfix(jdim,kdim,idim,w(lsnk0),mdim,ndim,w(ldata),
     .                 ibeg,iend,jbeg,jend,kbeg,kend,5,llev)
              end if
            enddo
            do nseg=1,nbckdim(nbl)
              if(abs(kbcinfo(nbl,nseg,1,2)).eq.2014) then
                ibeg=kbcinfo(nbl,nseg,2,2)
                iend=kbcinfo(nbl,nseg,3,2)
                jbeg=kbcinfo(nbl,nseg,4,2)
                jend=kbcinfo(nbl,nseg,5,2)
                kbeg=kdim
                kend=kdim
                mdim=jend-jbeg
                ndim=iend-ibeg
                ldata=lwdat(nbl,nseg,6)
                call lamfix(jdim,kdim,idim,w(lsnk0),mdim,ndim,w(ldata),
     .                 ibeg,iend,jbeg,jend,kbeg,kend,6,llev)
              end if
            enddo
c     I-dir:
            do nseg=1,nbci0(nbl)
              if(abs(ibcinfo(nbl,nseg,1,1)).eq.2014) then
                ibeg=1
                iend=1
                jbeg=ibcinfo(nbl,nseg,2,1)
                jend=ibcinfo(nbl,nseg,3,1)
                kbeg=ibcinfo(nbl,nseg,4,1)
                kend=ibcinfo(nbl,nseg,5,1)
                mdim=jend-jbeg
                ndim=kend-kbeg
                ldata=lwdat(nbl,nseg,1)
                call lamfix(jdim,kdim,idim,w(lsnk0),mdim,ndim,w(ldata),
     .                 ibeg,iend,jbeg,jend,kbeg,kend,1,llev)
              end if
            enddo
            do nseg=1,nbcidim(nbl)
              if(abs(ibcinfo(nbl,nseg,1,2)).eq.2014) then
                ibeg=idim
                iend=idim
                jbeg=ibcinfo(nbl,nseg,2,2)
                jend=ibcinfo(nbl,nseg,3,2)
                kbeg=ibcinfo(nbl,nseg,4,2)
                kend=ibcinfo(nbl,nseg,5,2)
                mdim=jend-jbeg
                ndim=kend-kbeg
                ldata=lwdat(nbl,nseg,2)
                call lamfix(jdim,kdim,idim,w(lsnk0),mdim,ndim,w(ldata),
     .                 ibeg,iend,jbeg,jend,kbeg,kend,2,llev)
              end if
            enddo
          endif
        enddo
      enddo
      end if
c
c********************************************************
c     use overlapped grid interpolation data to set 
c     conditions in overlapped cells
c********************************************************
c
      iflg = 0
      do 2100 igrid=1,ngrid
      ncg = ncgg(igrid)
      nbl = nblg(igrid)-1
      do 2110 icg=1,ncg+1
      nbl = nbl+1
      if (iovrlp(nbl).eq.1) then
         if (myid.eq.mblk2nd(1)) then
            if (iflg .eq.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)
            end if
            iflg = 1
         end if
         if (mblk2nd(nbl).eq.myid) then
            call lead(nbl,lw,lw2,maxbl)
            if (isklton.eq.1) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),350) nbl
            end if
            ldim = 5
            need  = (idim+1)*(jdim+1)*(kdim+1)*5
            if (need .gt. nwork) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''stopping...not enough work'',
     .         '' space for subroutine intrbc (q)'')')
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''need, have '',2i11)')
     .         need,nwork
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            call intrbc(w(lq),jdim,kdim,idim,nbl,ldim,maxbl,
     .                  iitot,lig,iipntsg,dxintg,dyintg,dzintg,
     .                  iiig,jjig,kkig,qb,w(lqj0),w(lqk0),
     .                  w(lqi0),wk,w(lbcj),w(lbck),w(lbci),nou,
     .                  bou,nbuf,ibufdim,irghost,1)
c           turbulence quantities
            if(iviscg(nbl,1).ge.2 .or. iviscg(nbl,2).ge.2 .or.
     .         iviscg(nbl,3).ge.2) then
               ldim = 1
               call intrbc(w(lvis),jdim,kdim,idim,nbl,ldim,maxbl,
     .                     iitot,lig,iipntsg,dxintg,dyintg,dzintg,
     .                     iiig,jjig,kkig,qb,w(lvj0),w(lvk0),
     .                     w(lvi0),wk,w(lbcj),w(lbck),w(lbci),nou,
     .                     bou,nbuf,ibufdim,irghost,2)
            end if
            if(iviscg(nbl,1).ge.4 .or. iviscg(nbl,2).ge.4 .or.
     .         iviscg(nbl,3).ge.4) then
               ldim = nummem
               call intrbc(w(lxib),jdim,kdim,idim,nbl,ldim,maxbl,
     .                     iitot,lig,iipntsg,dxintg,dyintg,dzintg,
     .                     iiig,jjig,kkig,qb,w(ltj0),w(ltk0),
     .                     w(lti0),wk,w(lbcj),w(lbck),w(lbci),nou,
     .                     bou,nbuf,ibufdim,irghost,3)
            end if
         end if
c
#if defined DIST_MPI
c
c        share qb information with the processors
c
         nvalqb  = iipntsg(nbl)
c
         if (myid.ne.mblk2nd(nbl) .and. myid.ne.myhost) then

c           recieve the data on non-local nodes
c
            nd_srce = mblk2nd(nbl)
c           receive q data 
            do ll=1,5
               mytag  = mytag_qb(ll) + nbl
               call MPI_Irecv (qb(lig(nbl),ll,1), nvalqb,
     .                         MY_MPI_REAL,
     .                         nd_srce, mytag, mycomm,
     .                         ireq_qb(nreq), ierr)
               nreq = nreq + 1
            end do
            if (ivmx.ge.2) then
c              receive vist3d data
               mytag = mytag_qb(6) + nbl
               call MPI_Irecv (qb(lig(nbl),1,2), nvalqb,
     .                         MY_MPI_REAL,
     .                         nd_srce, mytag, mycomm,
     .                         ireq_qb(nreq), ierr)
               nreq = nreq + 1
            end if
            if (ivmx.ge.4) then
c              receive turb. data            
                 do ll=1,nummem
                    mytag = mytag_qb(6+ll) + nbl
                    call MPI_Irecv (qb(lig(nbl),ll,3), nvalqb,
     .                              MY_MPI_REAL,
     .                              nd_srce, mytag, mycomm,
     .                              ireq_qb(nreq), ierr)
                    nreq = nreq + 1
                 end do
            end if
c
         else if (myid.eq.mblk2nd(nbl)) then
c
c           send the data to non-local nodes
c
            do inode = 1,nnodes
               if (inode .ne. myid) then
                  nd_dest = inode
c                 send q data
                  do ll=1,5
                     mytag  = mytag_qb(ll) + nbl
                     call MPI_Send (qb(lig(nbl),ll,1), nvalqb,
     .                              MY_MPI_REAL,
     .                              nd_dest, mytag,
     .                              mycomm, ierr)
                  end do
                  if (ivmx.ge.2) then
c                    send vist3d data
                     mytag = mytag_qb(6) + nbl
                     call MPI_Send (qb(lig(nbl),1,2), nvalqb,
     .                              MY_MPI_REAL,
     .                              nd_dest, mytag,
     .                              mycomm, ierr)
                  end if
                  if (ivmx.ge.4) then
c                    send turb. data
                     do ll=1,nummem
                        mytag = mytag_qb(6+ll) + nbl
                        call MPI_Send (qb(lig(nbl),ll,3), nvalqb,
     .                                 MY_MPI_REAL,
     .                                 nd_dest, mytag,
     .                                 mycomm, ierr)
                     end do
                  end if
               end if
            end do
c
         end if
c
#endif
      end if
 2110 continue
c
  350 format(1x,43hinterpolating for cells overlapped to block,i6)
c
      nblout = nblg(igrid)
#   ifdef FASTIO
      call writ_buffast(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl,28)
#   else
      call writ_buf(nblout,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl)
#   endif
c
 2100 continue
c
#if defined DIST_MPI
      if (nreq-1 .gt. 0) then
         do nnn=1,nreq-1
            call MPI_Wait(ireq_qb(nnn), istat, ierr)
         end do
      end if
c
#endif
c********************************************************
c     print iteration history
c********************************************************
c
      if (myid.eq.myhost) then
         if (irest.gt.0) then
            write(11,701)
  701       format(/,24h SUMMARY OF RESTART DATA)
            write(11,702)
  702       format(1X,9hiteration,5x,8hresidual,9x,2hCL,12x,2hCD,
     .             12x,2hCY,11x,3hCMY)
            do 1710 iter=1,ntt
            if (iter.lt.10 .or.
     .         (iter/50)*50.eq.iter .or. iter.gt.ntt-10) then
                write(11,703) iter,real(rms(iter)),real(clw(iter)),
     .                        real(cdw(iter)),real(cyw(iter)),
     .                        real(cmyw(iter))
            end if
  703       format(1x,i6,3x,5e14.5)
 1710       continue
         end if
      end if
c
      ntt = ntq
c
      if (myid.eq.myhost) then
         write(11,*)
         write(11,20)
      end if
   20 format(33h***** ENDING INITIALIZATION *****)
c
      return
      end
