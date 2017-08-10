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
      subroutine ronnie(maxbl,maxgr,nwork,intmax,nsub1,maxxe,mpatch,
     .                  nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose: establishes zone-to-zone communication file for cfl3d
c     multi-block applications with the generalized patch interface
c     option 
c***********************************************************************
c
      character*120 bou(ibufdim,nbuf)
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist
      character*50 string
c
      integer stats
c
      dimension ie_pat(5)
      dimension is_pat(5)
      dimension nou(nbuf)

      allocatable :: dthetx(:,:)
      allocatable :: dthety(:,:)
      allocatable :: dthetz(:,:)
      allocatable :: dx(:,:)
      allocatable :: dy(:,:)
      allocatable :: dz(:,:)
      allocatable :: eta2s(:,:)
      integer, allocatable :: etf1(:)
      integer, allocatable :: etf2(:)
      allocatable :: factjhi(:,:)
      allocatable :: factjlo(:,:)
      allocatable :: factkhi(:,:)
      allocatable :: factklo(:,:)
      allocatable :: idimg(:)
      allocatable :: iemg(:)
      allocatable :: ifiner(:)
      allocatable :: ifrom(:)
      allocatable :: igridg(:)
      allocatable :: iic0(:)
      allocatable :: iifit(:)
      allocatable :: iiint1(:)
      allocatable :: iiint2(:)
      allocatable :: iindex(:,:)
      allocatable :: iiorph(:)
      allocatable :: iitmax(:)
      allocatable :: iitoss(:)
      allocatable :: ireq_ar(:)
      allocatable :: isav_pat(:,:)
      allocatable :: isav_pat_b(:,:,:)
      allocatable :: itest(:)
      allocatable :: jdimg(:)
      allocatable :: jimage(:,:,:)
      allocatable :: jjmax1(:)
      allocatable :: jmm(:)
      allocatable :: jte(:)
      allocatable :: jtest(:)
      allocatable :: kdimg(:)
      allocatable :: kimage(:,:,:)
      allocatable :: kkmax1(:)
      allocatable :: kmm(:)
      allocatable :: kte(:)
      allocatable :: ktest(:)
      allocatable :: levelg(:)
      allocatable :: llimit(:)
      allocatable :: lout(:)
      allocatable :: lw(:,:)
      allocatable :: lw2(:,:)
      allocatable :: mblk2nd(:)
      allocatable :: mglevg(:)
      allocatable :: mmceta(:)
      allocatable :: mmcxie(:)
      allocatable :: nblcg(:)
      allocatable :: nblg(:)
      allocatable :: nblk1(:)
      allocatable :: nblk2(:)
      allocatable :: nblkpt(:)
      allocatable :: ncgg(:)
      allocatable :: ncheck(:)
      allocatable :: nemgl(:)
      allocatable :: seta(:,:,:)
      allocatable :: seta2(:,:)
      allocatable :: sxie(:,:,:)
      allocatable :: sxie2(:,:)
      allocatable :: temp(:)
      allocatable :: windex(:,:)
      allocatable :: work(:)
      allocatable :: x1(:,:)
      allocatable :: x2(:,:)
      allocatable :: xie2s(:,:)
      integer, allocatable :: xif1(:)
      integer, allocatable :: xif2(:)
      allocatable :: xmi(:,:,:)
      allocatable :: xmie(:,:,:)
      allocatable :: xorig(:)
      allocatable :: xte(:,:,:)
      allocatable :: y1(:,:)
      allocatable :: y2(:,:)
      allocatable :: ymi(:,:,:)
      allocatable :: ymie(:,:,:)
      allocatable :: yorig(:)
      allocatable :: yte(:,:,:)
      allocatable :: z1(:,:)
      allocatable :: z2(:,:)
      allocatable :: zmi(:,:,:)
      allocatable :: zmie(:,:,:)
      allocatable :: zorig(:)
      allocatable :: zte(:,:,:)

c
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist
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
      common /tracer/ itrace
      common /zero/ iexp
      common /unit5/ iunit5
      common /conversion/ radtodeg
      common /igrdtyp/ ip3dgrd,ialph
      common /is_dpatch/ maxdcnt
c
      ierrflg    = -1
c
c***********************************************************************
c     memory allocation
c***********************************************************************
c
      memuse = 0
c
      allocate( dthetx(intmax,nsub1), stat=stats)
      call umalloc(intmax*nsub1,0,'dthetx',memuse,stats)
      allocate( dthety(intmax,nsub1), stat=stats)
      call umalloc(intmax*nsub1,0,'dthety',memuse,stats)
      allocate( dthetz(intmax,nsub1), stat=stats)
      call umalloc(intmax*nsub1,0,'dthetz',memuse,stats)
      allocate( dx(intmax,nsub1), stat=stats)
      call umalloc(intmax*nsub1,0,'dx',memuse,stats)
      allocate( dy(intmax,nsub1), stat=stats)
      call umalloc(intmax*nsub1,0,'dy',memuse,stats)
      allocate( dz(intmax,nsub1), stat=stats)
      call umalloc(intmax*nsub1,0,'dx',memuse,stats)
      allocate( eta2s(mpatch+2,mpatch+2), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2),0,'eta2s',memuse,stats)
      allocate( etf1(nsub1), stat=stats)
      call umalloc(nsub1,1,'etf1',memuse,stats)
      allocate( etf2(nsub1), stat=stats)
      call umalloc(nsub1,1,'etf2',memuse,stats)
      allocate( factjhi(intmax,nsub1), stat=stats)
      call umalloc(intmax*nsub1,0,'factjhi',memuse,stats)
      allocate( factjlo(intmax,nsub1), stat=stats)
      call umalloc(intmax*nsub1,0,'factjlo',memuse,stats)
      allocate( factkhi(intmax,nsub1), stat=stats)
      call umalloc(intmax*nsub1,0,'factkhi',memuse,stats)
      allocate( factklo(intmax,nsub1), stat=stats)
      call umalloc(intmax*nsub1,0,'factklo',memuse,stats)
      allocate( idimg(maxbl), stat=stats)
      call umalloc(maxbl,1,'idimg',memuse,stats)
      allocate( iemg(maxgr), stat=stats)
      call umalloc(maxgr,1,'iemg',memuse,stats)
      allocate( ifiner(intmax), stat=stats)
      call umalloc(intmax,1,'ifiner',memuse,stats)
      allocate( ifrom(nsub1), stat=stats)
      call umalloc(nsub1,0,'ifrom',memuse,stats)
      allocate( igridg(maxbl), stat=stats)
      call umalloc(maxbl,1,'igridg',memuse,stats)
      allocate( iic0(intmax), stat=stats)
      call umalloc(intmax,1,'iic0',memuse,stats)
      allocate( iifit(intmax), stat=stats)
      call umalloc(intmax,1,'iifit',memuse,stats)
      allocate( iiint1(nsub1), stat=stats)
      call umalloc(nsub1,1,'iiint1',memuse,stats)
      allocate( iiint2(nsub1), stat=stats)
      call umalloc(nsub1,1,'iiint2',memuse,stats)
      allocate( iindex(intmax,6*nsub1+9), stat=stats)
      call umalloc(intmax*(6*nsub1+9),1,'iindex',memuse,stats)
      allocate( iiorph(intmax), stat=stats)
      call umalloc(intmax,1,'iiorph',memuse,stats)
      allocate( iitmax(intmax), stat=stats)
      call umalloc(intmax,1,'iitmax',memuse,stats)
      allocate( iitoss(intmax), stat=stats)
      call umalloc(intmax,1,'iitoss',memuse,stats)
      allocate( ireq_ar(intmax*3), stat=stats)
      call umalloc(intmax*3,1,'ireq_ar',memuse,stats)
      allocate( isav_pat(intmax,17), stat=stats)
      call umalloc(intmax*17,1,'isav_pat',memuse,stats)
      allocate( isav_pat_b(intmax,nsub1,6), stat=stats)
      call umalloc(intmax*nsub1*6,1,'isav_pat_b',memuse,stats)
      allocate( itest(maxgr), stat=stats)
      call umalloc(maxgr,1,'itest',memuse,stats)
      allocate( jdimg(maxbl), stat=stats)
      call umalloc(maxbl,1,'jdimg',memuse,stats)
      allocate( jimage(nsub1,mpatch+2,mpatch+2), stat=stats)
      call umalloc(nsub1*(mpatch+2)*(mpatch+2),1,'jimage',memuse,stats)
      allocate( jjmax1(nsub1), stat=stats)
      call umalloc(nsub1,1,'jjmax1',memuse,stats)
      allocate( jmm(mpatch+2), stat=stats)
      call umalloc(mpatch+2,1,'jmm',memuse,stats)
      allocate( jte(nsub1), stat=stats)
      call umalloc(nsub1,1,'jte',memuse,stats)
      allocate( jtest(maxgr), stat=stats)
      call umalloc(maxgr,1,'jtest',memuse,stats)
      allocate( kdimg(maxbl), stat=stats)
      call umalloc(maxbl,1,'kdimg',memuse,stats)
      allocate( kimage(nsub1,mpatch+2,mpatch+2), stat=stats)
      call umalloc(nsub1*(mpatch+2)*(mpatch+2),1,'kimage',memuse,stats)
      allocate( kkmax1(nsub1), stat=stats)
      call umalloc(nsub1,1,'kkmax1',memuse,stats)
      allocate( kmm(mpatch+2), stat=stats)
      call umalloc(mpatch+2,1,'kmm',memuse,stats)
      allocate( kte(nsub1), stat=stats)
      call umalloc(nsub1,1,'kte',memuse,stats)
      allocate( ktest(maxgr), stat=stats)
      call umalloc(maxgr,1,'ktest',memuse,stats)
      allocate( levelg(maxbl), stat=stats)
      call umalloc(maxbl,1,'levelg',memuse,stats)
      allocate( llimit(intmax), stat=stats)
      call umalloc(intmax,1,'llimit',memuse,stats)
      allocate( lout(nsub1), stat=stats)
      call umalloc(nsub1,1,'lout',memuse,stats)
      allocate( lw(65,maxbl), stat=stats)
      call umalloc(65*maxbl,1,'lw',memuse,stats)
      allocate( lw2(43,maxbl), stat=stats)
      call umalloc(43*maxbl,1,'lw2',memuse,stats)
      allocate( mblk2nd(maxbl), stat=stats)
      call umalloc(maxbl,1,'mblk2nd',memuse,stats)
      allocate( mglevg(maxbl), stat=stats)
      call umalloc(maxbl,1,'mglevg',memuse,stats)
      allocate( mmceta(intmax), stat=stats)
      call umalloc(intmax,1,'mmceta',memuse,stats)
      allocate( mmcxie(intmax), stat=stats)
      call umalloc(intmax,1,'mmcxie',memuse,stats)
      allocate( nblcg(maxbl), stat=stats)
      call umalloc(maxbl,1,'nblcg',memuse,stats)
      allocate( nblg(maxgr), stat=stats)
      call umalloc(maxgr,1,'nblg',memuse,stats)
      allocate( nblk1(mpatch+2), stat=stats)
      call umalloc(mpatch+2,1,'nblk1',memuse,stats)
      allocate( nblk2(mpatch+2), stat=stats)
      call umalloc(mpatch+2,1,'nblk2',memuse,stats)
      allocate( nblkpt(maxxe), stat=stats)
      call umalloc(maxxe,1,'nblkpt',memuse,stats)
      allocate( ncgg(maxgr), stat=stats)
      call umalloc(maxgr,1,'ncgg',memuse,stats)
      allocate( ncheck(maxbl), stat=stats)
      call umalloc(maxbl,1,'ncheck',memuse,stats)
      allocate( nemgl(maxbl), stat=stats)
      call umalloc(maxbl,1,'nemgl',memuse,stats)
      allocate( seta(mpatch+2,mpatch+2,nsub1), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2)*nsub1,0,'seta',memuse,stats)
      allocate( seta2(mpatch+2,mpatch+2), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2),0,'seta2',memuse,stats)
      allocate( sxie(mpatch+2,mpatch+2,nsub1), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2)*nsub1,0,'sxie',memuse,stats)
      allocate( sxie2(mpatch+2,mpatch+2), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2),0,'sxie2',memuse,stats)
      allocate( temp((mpatch+2)*(mpatch+2)), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2),0,'temp',memuse,stats)
      allocate( windex(maxxe,2), stat=stats)
      call umalloc(maxxe*2,0,'windex',memuse,stats)
      allocate( work(nwork), stat=stats)
      call umalloc(nwork,0,'work',memuse,stats)
      allocate( x1(mpatch+2,mpatch+2), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2),0,'x1',memuse,stats)
      allocate( x2(mpatch+2,mpatch+2), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2),0,'x2',memuse,stats)
      allocate( xie2s(mpatch+2,mpatch+2), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2),0,'xie2s',memuse,stats)
      allocate( xif1(nsub1), stat=stats)
      call umalloc(nsub1,1,'xif1',memuse,stats)
      allocate( xif2(nsub1), stat=stats)
      call umalloc(nsub1,1,'xif2',memuse,stats)
      allocate( xmi(mpatch+2,mpatch+2,nsub1), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2)*nsub1,0,'xmi',memuse,stats)
      allocate( xmie(mpatch+2,mpatch+2,nsub1), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2)*nsub1,0,'xmie',memuse,stats)
      allocate( xorig(maxbl), stat=stats)
      call umalloc(maxbl,0,' xorig',memuse,stats)
      allocate( xte(mpatch+2,mpatch+2,nsub1), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2)*nsub1,0,'xte',memuse,stats)
      allocate( y1(mpatch+2,mpatch+2), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2),0,'y1',memuse,stats)
      allocate( y2(mpatch+2,mpatch+2), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2),0,'y2',memuse,stats)
      allocate( ymi(mpatch+2,mpatch+2,nsub1), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2)*nsub1,0,'ymi',memuse,stats)
      allocate( ymie(mpatch+2,mpatch+2,nsub1), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2)*nsub1,0,'ymie',memuse,stats)
      allocate( yorig(maxbl), stat=stats)
      call umalloc(maxbl,0,'yorig',memuse,stats)
      allocate( yte(mpatch+2,mpatch+2,nsub1), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2)*nsub1,0,'yte',memuse,stats)
      allocate( z1(mpatch+2,mpatch+2), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2),0,'z1',memuse,stats)
      allocate( z2(mpatch+2,mpatch+2), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2),0,'z2',memuse,stats)
      allocate( zmi(mpatch+2,mpatch+2,nsub1), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2)*nsub1,0,'zmi',memuse,stats)
      allocate( zmie(mpatch+2,mpatch+2,nsub1), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2)*nsub1,0,'zmie',memuse,stats)
      allocate( zorig(maxbl), stat=stats)
      call umalloc(maxbl,0,'zorig',memuse,stats)
      allocate( zte(mpatch+2,mpatch+2,nsub1), stat=stats)
      call umalloc((mpatch+2)*(mpatch+2)*nsub1,0,'zte',memuse,stats)
c
      string = ' '
      call cputim(0,1,string,myhost,myid,mycomm,25)
c
c     determine machine zero for use in setting tolerances
c     (10.**(-iexp) is machine zero)
c
      icount = 0
      compare = 1.0
      do 4000 i = 1,20
        icount = icount + 1
        add = 1.
        do n=1,i
          add = add*.1
        enddo
        x11 = compare + add
        if(x11.eq.compare)then
         iexp = i-1
         x    = .1**iexp
         goto 4010
        end if
 4000 continue
 4010 continue
c
      pi       = 4.*atan(1.0)
      radtodeg = 180.e0/pi
c
c     set dummy values for parallel-related variables
c
      nnodes = 1
      myhost = 0
      myid   = 0
      mycomm = 0
      do nn=1,maxbl
         mblk2nd(nn) = myhost
      end do
c
c     initialize output buffers
c
      do ll=1,nbuf
         nou(ll) = 0
         do mm=1,ibufdim
            bou(mm,ll) = ' '
         end do
      end do
c
c     output banner
c
      write(25,83)
      write(25,83)
      write(25,87)
      write(25,9900)
 9900 format(2(2h *),38h           RONNIE - CFL3D PATCHED-GRID,
     .13h PREPROCESSOR,11x,2(2h *))
      write(25,87)
      write(25,9990)
 9990 format(2(2h *),43h   VERSION 6.7 :  Computational Fluids Lab,,
     .15h Mail Stop 128,,4x,2(2h *),
     ./2(2h *),18x,41hNASA Langley Research Center, Hampton, VA,
     .3x,2(2h *),/2(2h *),18x,33hRelease Date:  February  1, 2017.,
     .11x,2(2h *))
      write(25,87)
      write(25,83)
      write(25,83)
   83 format(35(2h *))
   87 format(2(2h *),62x,2(2h *))
c
#ifdef CRAY_TIME
c     cray_time implies cray (always double precision)
      write(25,12) float(memuse)/1.e6
#else
#   ifdef DBLE_PRECSN
      write(25,12) float(memuse)/1.e6
#   else
      write(25,13) float(memuse)/1.e6
#   endif
#endif
   12 format(/,' memory allocation: ',f12.6,' Mbytes, double precision')
   13 format(/,' memory allocation: ',f12.6,' Mbytes, single precision')
c
      write(25,88)
   88 format(/19hinput/output files:)
c
      read(iunit5,*)
      read(iunit5,'(a60)')grid
      read(iunit5,'(a60)')output
      read(iunit5,'(a60)')patch
c
      write(25,'(''  '',a60)')grid
      write(25,'(''  '',a60)')output
      write(25,'(''  '',a60)')patch
c
      read(iunit5,*)
      read(iunit5,*) ioflag,itrace
      write(25,398) ioflag,itrace
  398 format(/,6hioflag,4x,6hitrace/,i6,4x,i6)
c
      read(iunit5,10)(title(i),i=1,20)
   10 format(20a4)
      write(25,111)
 111  format(/5htitle)
      write(25,11)(title(i),i=1,20)
   11 format(2h  ,20a4)
c
      read(iunit5,*)
      read(iunit5,*) ngrid
      write(25,1638)
 1638 format(/,15hgrid/level data)
      write(25,1639)
 1639 format(1x,5hngrid)
      write(25,36) ngrid
   36 format(13i6)
c
      ip3dgrd = 0
      if (ngrid.lt.0) then
         ip3dgrd = 1
         ngrid =  iabs(ngrid)
      end if
      nchk = maxgr-ngrid
c
c     check maximum number of grids
c
      if (nchk.lt.0) then
         write(25,1492)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 1492 format(55h stopping - insufficient maximum number of grids(maxgr))
c
      read(iunit5,*)
      nbl = 0
      iemtot = 0
      write(25,1631)
 1631 format(3x,3hncg,3x,3hiem,2x,4hidim,2x,4hjdim,2x,4hkdim)
      do 7001 igrid=1,ngrid
      nbl = nbl+1
      read(iunit5,*) ncg,iem,idim,jdim,kdim
      write(25,36) ncg,iem,idim,jdim,kdim
      iemtot      = iemtot+iem
      ncgg(igrid) = ncg
      if(igrid.eq.1) then
        ncgmax = ncg
        iemmax = iem
      else
        ncgmax = max(ncgmax,ncg)
        iemmax = max(iemmax,iem)
      end if
      iemg(igrid)   = iem
      nblg(igrid) = nbl
      idimg(nbl)  = idim
      jdimg(nbl)  = jdim
      kdimg(nbl)  = kdim
c
      if (ncg.gt.0) then
         if (iem.gt.0) then
            write(25,*)' embedded grids must have ncg = 0'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         do 6885 n=1,ncg
         nbl        = nbl+1
         idimg(nbl) = idimg(nbl-1)/2+1
         jdimg(nbl) = jdimg(nbl-1)/2+1
         kdimg(nbl) = kdimg(nbl-1)/2+1
         if (idimg(nbl-1).le.2) then
c        2-d meshes
            idimg(nbl) = idimg(nbl-1)
         end if
         istop=0
         if (float(idimg(nbl-1)/2) .eq. float(idimg(nbl-1))/2. .and.
     .    idim .gt. 2) then
           write(25,'('' Cannot create coarser level for idim past'',
     .      i6)') idimg(nbl-1)
           istop=1
         end if
         if (float(jdimg(nbl-1)/2) .eq. float(jdimg(nbl-1))/2.) then
           write(25,'('' Cannot create coarser level for jdim past'',
     .      i6)') jdimg(nbl-1)
           istop=1
         end if
         if (float(kdimg(nbl-1)/2) .eq. float(kdimg(nbl-1))/2.) then
           write(25,'('' Cannot create coarser level for kdim past'',
     .      i6)') kdimg(nbl-1)
           istop=1
         end if
         if (istop .eq. 1) then
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
 6885    continue
      end if
 7001 continue
c
      mseq = 1
      if (ncgmax.gt.0) mgflag = 1
      if (iemmax.gt.0) mgflag = 2
c
      if (mseq.gt.1) then
         nbl = 0
         do igrid=1,ngrid
            nbl = nbl+1
            ncg = ncgg(igrid)
            if (ncg.gt.0) nbl = nbl+ncg
         end do
      end if
c
      do m=1,mseq
         mglevg(m) = ncgmax + 1
         nemgl(m)  = iemmax
      end do
c
c     check maximum number of blocks
c
      nblock = nbl
      nchk   = maxbl-nblock
      if (nchk.lt.0) then
         write(25,1649) maxbl,nchk
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 1649 format(1x,12h maxbl,nchk=,2i5)
c
c     determine levelt and levelb
c
c     levelt - starting level for multigrid/time cycling
c     levelb - ending level for multigrid/time cycling
c
      do 17 m=1,mseq
      levelt(m) = ncgmax-(mseq-m)+nemgl(m)+1
      levelb(m) = levelt(m)-(mglevg(m)-1)-nemgl(m)
      if (levelb(m).lt.1) then
         write(25,157)m,levelt(m),levelb(m),ncgmax
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
  157 format(1x,42herror in input, m, levelt, levelb, ncgmax=,4i5)
   17 continue
c
      icall = 1
      iunit = 25
      imode = 0
      call global2(maxbl,maxgr,nsub1,ninter,intmax,ngrid,idimg,
     .             jdimg,kdimg,levelg,ncgg,nblg,iindex,llimit,
     .             iitmax,mmcxie,mmceta,ncheck,iifit,iic0,
     .             iiorph,iitoss,ifiner,dx,dy,dz,dthetx,
     .             dthety,dthetz,myid,mpatch,maxxe,icall,iunit,
     .             nou,bou,ibufdim,nbuf,ifrom,xif1,etf1,xif2,
     .             etf2,igridg,iemg,nblock,ioflag,imode)
c
c     data for grids stored in array w
c
c     starting locations for data in w are stored in
c     lw(i,n) where i is type of data and n is block number
c     the lw and lw2 arrays are sized as in cfl3d to permit
c     reuse of certain subroutines, although only locations 10-12
c     are actually required here
c
c     i correspondence to data is as follows
c     i=10  x
c       11  y
c       12  z
c
c
      nstart = 1
c
      do 1000 nbl=1,nblock
      j = jdimg(nbl)
      k = kdimg(nbl)
      i = idimg(nbl)
c
      igrid = igridg(nbl)
      iem   = iemg(igrid)
c
      j1 = j-1
      k1 = k-1
      i1 = i-1
c
c     set some dummy values just so they are defined in pre_patch
c
      do ll=1,9
         lw(ll,nbl) = nstart
      end do
c
c     x(jdim,kdim,idim)
      lw(10,nbl)  = nstart
c
c     y(jdim,kdim,idim)
      ns         = j*k*i
      lw(11,nbl) = lw(10,nbl)+ns
c
c     z(jdim,kdim,idim)          
      ns         = j*k*i
      lw(12,nbl) = lw(11,nbl)+ns
c
      ns         = j*k*i
      nstart = lw(12,nbl)+ns
c
c     set some more dummy values just so they are defined in pre_patch
c
      do ll=13,65
         lw(ll,nbl) = nstart
      end do
c
c     set up lw2 array (istalling dummy values for those entires
c     not used/defined in ronnie
c
      jdim        = jdimg(nbl)
      kdim        = kdimg(nbl)
      idim        = idimg(nbl)
      lw2(1, nbl) = jdimg(nbl)
      lw2(2, nbl) = kdimg(nbl)
      lw2(3, nbl) = idimg(nbl)
c     note: lw2(4, nbl) = nbl not correct for embedded, but of
c     no consequnce here
      lw2(4, nbl) = nbl
      lw2(5, nbl) = jdim/2+1
      lw2(6, nbl) = kdim/2+1
      lw2(7, nbl) = idim/2+1
      lw2(8, nbl) = 0
      lw2(9, nbl) = 0
      lw2(10,nbl) = 0
      lw2(11,nbl) = 0
      lw2(12,nbl) = 0
      lw2(13,nbl) = 0
      lw2(14,nbl) = 0
      lw2(15,nbl) = 0
      lw2(16,nbl) = 0
      lw2(17,nbl) = 0
      lw2(18,nbl) = 0
      lw2(19,nbl) = 0
      lw2(20,nbl) = 0
      lw2(21,nbl) = 0
      lw2(22,nbl) = 0
      lw2(23,nbl) = 0
      lw2(24,nbl) = 0
      lw2(25,nbl) = 0
      lw2(26,nbl) = 0
      lw2(27,nbl) = 0
      lw2(28,nbl) = 0
      lw2(29,nbl) = 0
      lw2(30,nbl) = 0
      lw2(31,nbl) = 0
      lw2(32,nbl) = 0
      lw2(33,nbl) = 0
      lw2(34,nbl) = 0
      lw2(35,nbl) = 0
      lw2(36,nbl) = 0
      lw2(37,nbl) = 0
      lw2(38,nbl) = 0
      lw2(39,nbl) = 0
      lw2(40,nbl) = 0
c
1000  continue
c
      if (nstart.gt.nwork) then
         write(25,*)' must increase parameter nwork to ',nstart
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      call setup(lw,lw2,work,nstart,itest,jtest,ktest,maxbl,maxgr,
     .           ngrid,ncgg,iemg,nblg)
c
c     set up auxiliary arrays for patching subroutines
c
      ierrflg    = -1
      icount_pat = 0
      levt       = levelt(mseq)
      do levl = 1,levt
         is_pat(levl) = icount_pat + 1
         do 6909 nbl=1,nblock
            if (levl.ne.levelg(nbl)) go to 6909
            icount_pat1 = icount_pat
            call pre_patch(nbl,lw,icount_pat,ninter,
     .                     iindex,intmax,nsub1,isav_pat,
     .                     isav_pat_b,jjmax1,kkmax1,
     .                     iiint1,iiint2,maxbl,jdimg,kdimg,idimg,
     .                     ierrflg)
            if (icount_pat .gt. icount_pat1) then
               ie_pat(levl) = icount_pat
            end if
 6909    continue
      end do
c
      maxdcnt = icount_pat
c
c     cycle through blocks, check for and establishing connection
c     information for those blocks with patching
c
      write(25,91)
   91 format(/,/,1x,46hBEGINNING GENERALIZED-COORDINATE INTERPOLATION)
c
      it_thro = 0
      ncall   = 1
      mgwk    = nstart
      ioutpt  = 1
      mwork = nwork - nstart
c
      do 7000 nbl=1,nblock
      call lead(nbl,lw,lw2,maxbl)
      it_thro = it_thro + 1
      call patcher(nbl,lw,work,mgwk,work(mgwk+1),mwork,ncall,ioutpt,
     .             it_thro,maxbl,nsub1,intmax,maxxe,mpatch,jdimg,kdimg,
     .             idimg,windex,ninter,iindex,llimit,iitmax,
     .             mmcxie,mmceta,ncheck,iifit,nblkpt,iic0,
     .             iiorph,iitoss,ifiner,factjlo,factjhi,
     .             factklo,factkhi,dx,dy,dz,dthetx,dthety,
     .             dthetz,isav_pat,isav_pat_b,
     .             xte,yte,zte,xmi,ymi,zmi,xmie,ymie,zmie,
     .             jjmax1,kkmax1,jimage,kimage,xorig,yorig,zorig,
     .             jte,kte,sxie,seta,sxie2,seta2,xie2s,eta2s,
     .             temp,x2,y2,z2,nblk1,nblk2,jmm,kmm,x1,y1,z1,
     .             lout,xif1,xif2,etf1,etf2,ireq_ar,
     .             myid,myhost,mycomm,mblk2nd,nou,bou,nbuf,
     .             ibufdim)
      call writ_buf(nbl,25,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl)
      call writ_buf(nbl,9,nou,bou,nbuf,ibufdim,myhost,myid,
     .              mycomm,mblk2nd,maxbl)
7000  continue
c
c     write out connectivity file
c
      write(2) ninter
      do 1500 n=1,ninter
      write(2) iindex(n,1)
      nfb = iindex(n,1)
      write(2) (iindex(n,1+ll),ll=1,nfb)
      write(2) iindex(n,nfb+2)
      write(2) (iindex(n,nfb+2+ll),ll=1,nfb)
      write(2) iindex(n,2*nfb+3)
      write(2) iindex(n,2*nfb+4)
      write(2) iindex(n,2*nfb+5)
      write(2) iindex(n,2*nfb+6)
      write(2) iindex(n,2*nfb+7)
      write(2) iindex(n,2*nfb+8)
      write(2) iindex(n,2*nfb+9)
      nbl = iindex(n,nfb+2)
      lst    = iindex(n,2*nfb+5)
      len = lst + iindex(n,2*nfb+4) - 1
      write(2) (nblkpt(nnn),nnn=lst,len)
      write(2) ((windex(nnn,ll),nnn=lst,len),ll=1,2)
1500  continue
c
c     write out angular displacements for ioflag=2 input option
c
      if (ioflag .gt. 1) then
         do n=1,ninter
            nfb = iindex(n,1)
            do ll=1,nfb
               write(2) dthetx(n,ll),dthety(n,ll),dthetz(n,ll)
            end do
         end do
      end if
c 
      write(25,*)
      write(25,9998) patch
9998  format(2x,26hwriting connectivity file:,/,2x,a60,/)
      nfb = iindex(ninter,1)
      lst = iindex(ninter,2*nfb+5)
      len = lst + iindex(ninter,2*nfb+4) - 1
      write(25,9997) len
9997  format('  minimum dimension for parameter MAXXE in CFL3D:',i6)
c
      string = '    timing for complete run - time in seconds     '
      call cputim(-1,nnodes,string,myhost,myid,mycomm,25)
c
c     free memory
c
      ifree = 1
      if (ifree.gt.0) then
         deallocate(work)
         deallocate(lw)
         deallocate(lw2)
         deallocate(mglevg)
         deallocate(nemgl)
         deallocate(levelg)
         deallocate(jtest)
         deallocate(ktest)
         deallocate(itest)
         deallocate(nblg)
         deallocate(iemg)
         deallocate(igridg)
         deallocate(jdimg)
         deallocate(kdimg)
         deallocate(idimg)
         deallocate(nblcg)
         deallocate(ncgg)
         deallocate(xorig)
         deallocate(yorig)
         deallocate(zorig)
         deallocate(windex)
         deallocate(iindex)
         deallocate(llimit)
         deallocate(iitmax)
         deallocate(mmcxie)
         deallocate(mmceta)
         deallocate(ncheck)
         deallocate(iifit)
         deallocate(nblkpt)
         deallocate(iic0)
         deallocate(iiorph)
         deallocate(iitoss)
         deallocate(ifiner)
         deallocate(dx)
         deallocate(dy)
         deallocate(dz)
         deallocate(dthetx)
         deallocate(dthety)
         deallocate(dthetz)
         deallocate(lout)
         deallocate(xif1)
         deallocate(xif2)
         deallocate(etf1)
         deallocate(etf2)
         deallocate(jjmax1)
         deallocate(kkmax1)
         deallocate(iiint1)
         deallocate(iiint2)
         deallocate(jimage)
         deallocate(kimage)
         deallocate(jte)
         deallocate(kte)
         deallocate(jmm)
         deallocate(kmm)
         deallocate(nblk1)
         deallocate(nblk2)
         deallocate(xte)
         deallocate(yte)
         deallocate(zte)
         deallocate(xmi)
         deallocate(ymi)
         deallocate(zmi)
         deallocate(xmie)
         deallocate(ymie)
         deallocate(zmie)
         deallocate(sxie)
         deallocate(seta)
         deallocate(sxie2)
         deallocate(seta2)
         deallocate(xie2s)
         deallocate(eta2s)
         deallocate(temp)
         deallocate(x2)
         deallocate(y2)
         deallocate(z2)
         deallocate(x1)
         deallocate(y1)
         deallocate(z1)
         deallocate(factjlo)
         deallocate(factjhi)
         deallocate(factklo)
         deallocate(factkhi)
         deallocate(ifrom)
         deallocate(isav_pat)
         deallocate(isav_pat_b)
         deallocate(mblk2nd)
         deallocate(ireq_ar)
      end if
c
      return
      end
