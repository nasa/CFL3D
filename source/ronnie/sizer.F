      subroutine sizer(maxgr,maxbl,maxxe,nsub1,intmax,mpatch,nwork,
     .                 nbuf0,ibufdim0,maxgr0,maxbl0,maxxe0,nsub10,
     .                 intmax0,mpatch0,imode)
c
c     $Id$
c
c***********************************************************************
c     Purpose: establishes necessary parameters for the "ronnie" 
c     preprocessor for cfl3d; ronnie determines the interpolation
c     coefficients used in the generalized patch interface option
c     for cfl3d
c
c     imode governs whether the routine is being used a a stand-alone
c     preprocessor, or as part of the dynamic memory allocation in
c     ronnie:
c
c     imode = 0  stand-alone preprocessor
c             1  part of the dynamic memory allocation in ronnie
c***********************************************************************
c
c     maxbl  - maximum number of blocks
c     maxgr  - maximum number of grids
c     nwork  - size of 1-d array used to store xyz
c     intmax - maximum number of block interpolations
c     mpatch - maximum dimension of any block face involved in a patch
c     maxxe  - size of 1-d array used to store interpolation coefficients for
c              all patch interfaces, including those on coarser blocks
c     nsub1  - maximum number of blocks a single patch face may be
c              interpolated from
c
      character*120 bou(ibufdim0,nbuf0)
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist
      character*50 string
c
      integer xif1(nsub10),xif2(nsub10),etf1(nsub10),etf2(nsub10)
c
      dimension levelg(maxbl0)
      dimension mglevg(maxbl0),nemgl(maxbl0)
      dimension igridg(maxbl0),
     .          jdimg(maxbl0),kdimg(maxbl0),idimg(maxbl0),nblcg(maxbl0)
      dimension ncgg(maxgr0),nblg(maxgr0),iemg(maxgr0)
      dimension nou(nbuf0)
      dimension jjmax1(nsub10),kkmax1(nsub10)
      dimension iiint1(nsub10),iiint2(nsub10)
      dimension iindex(intmax0,6*nsub10+9),
     .          llimit(intmax0),iitmax(intmax0),mmcxie(intmax0),
     .          mmceta(intmax0),ncheck(maxbl0),iifit(intmax0),
     .          iic0(intmax0),iiorph(intmax0),
     .          iitoss(intmax0),ifiner(intmax0),ifrom(nsub10)
      dimension dx(intmax0,nsub10),dy(intmax0,nsub10),
     .          dz(intmax0,nsub10),dthetx(intmax0,nsub10),
     .          dthety(intmax0,nsub10),dthetz(intmax0,nsub10)
      dimension isav_pat(intmax0,17),isav_pat_b(intmax0,nsub10,6)
      dimension mblk2nd(maxbl0)
      dimension is_pat(5),ie_pat(5)
c
      common /params/ lmaxgr,lmaxbl,lmxseg,lmaxcs,lnplts,lmxbli,lmaxxe,
     .                lnsub1,lintmx,lmxxe,liitot,isum,lncycm,
     .                isum_n,lminnode,isumi,isumi_n,lmptch,
     .                lmsub1,lintmax,libufdim,lnbuf,llbcprd,
     .                llbcemb,llbcrad,lnmds,lmaxaes,lnslave,lmxsegdg,
     .                lnmaster
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /unit5/ iunit5
      common /conversion/ radtodeg
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist
c
      pointer (ip_levelg,levelg),(ip_mglevg,mglevg),(ip_nemgl,nemgl)
      pointer (ip_ncgg,ncgg),(ip_nblg,nblg),(ip_iemg,iemg)
      pointer (ip_igridg,igridg),(ip_nblcg,nblcg)
      pointer (ip_iindex,iindex),(ip_iitoss,iitoss),(ip_ifiner,ifiner),
     .        (ip_llimit,llimit),(ip_iitmax,iitmax),(ip_mmcxie,mmcxie),
     .        (ip_mmceta,mmceta),(ip_ncheck,ncheck),(ip_iifit,iifit),
     .        (ip_iic0,iic0),(ip_iiorph,iiorph)
      pointer (ip_dx,dx),(ip_dy,dy),(ip_dz,dz),
     .        (ip_dthetx,dthetx),(ip_dthety,dthety),
     .        (ip_dthetz,dthetz)
      pointer (ip_xif1,xif1),(ip_xif2,xif2),
     .        (ip_etf1,etf1),(ip_etf2,etf2),(ip_iiint1,iiint1),
     .        (ip_iiint2,iiint2),(ip_jjmax1,jjmax1),
     .        (ip_kkmax1,kkmax1),(ip_ifrom,ifrom)
      pointer (ip_isav_pat,isav_pat),(ip_isav_pat_b,isav_pat_b)
      pointer (ip_idimg,idimg),(ip_jdimg,jdimg),(ip_kdimg,kdimg)
      pointer (ip_mblk2nd,mblk2nd)
c
      ierrflg    = - 99
c
      icall  = 0
      memuse = 0
c
      call umalloc(ip_mglevg,maxbl0,1,icall,memuse)
      call umalloc(ip_nemgl,maxbl0,1,icall,memuse)
      call umalloc(ip_levelg,maxbl0,1,icall,memuse)
      call umalloc(ip_nblg,maxgr0,1,icall,memuse)
      call umalloc(ip_iemg,maxgr0,1,icall,memuse)
      call umalloc(ip_igridg,maxbl0,1,icall,memuse)
      call umalloc(ip_jdimg,maxbl0,1,icall,memuse)
      call umalloc(ip_kdimg,maxbl0,1,icall,memuse)
      call umalloc(ip_idimg,maxbl0,1,icall,memuse)
      call umalloc(ip_nblcg,maxbl0,1,icall,memuse)
      call umalloc(ip_ncgg,maxgr0,1,icall,memuse)
      call umalloc(ip_iindex,intmax0*(6*nsub10+9),1,icall,memuse)
      call umalloc(ip_llimit,intmax0,1,icall,memuse)
      call umalloc(ip_iitmax,intmax0,1,icall,memuse)
      call umalloc(ip_mmcxie,intmax0,1,icall,memuse)
      call umalloc(ip_mmceta,intmax0,1,icall,memuse)
      call umalloc(ip_ncheck,maxbl0,1,icall,memuse)
      call umalloc(ip_iifit,intmax0,1,icall,memuse)
      call umalloc(ip_iic0,intmax0,1,icall,memuse)
      call umalloc(ip_iiorph,intmax0,1,icall,memuse)
      call umalloc(ip_iitoss,intmax0,1,icall,memuse)
      call umalloc(ip_ifiner,intmax0,1,icall,memuse)
      call umalloc(ip_dx,intmax0*nsub10,0,icall,memuse)
      call umalloc(ip_dy,intmax0*nsub10,0,icall,memuse)
      call umalloc(ip_dz,intmax0*nsub10,0,icall,memuse)
      call umalloc(ip_dthetx,intmax0*nsub10,0,icall,memuse)
      call umalloc(ip_dthety,intmax0*nsub10,0,icall,memuse)
      call umalloc(ip_dthetz,intmax0*nsub10,0,icall,memuse)
      call umalloc(ip_xif1,nsub10,1,icall,memuse)
      call umalloc(ip_xif2,nsub10,1,icall,memuse)
      call umalloc(ip_etf1,nsub10,1,icall,memuse)
      call umalloc(ip_etf2,nsub10,1,icall,memuse)
      call umalloc(ip_jjmax1,nsub10,1,icall,memuse)
      call umalloc(ip_kkmax1,nsub10,1,icall,memuse)
      call umalloc(ip_iiint1,nsub10,1,icall,memuse)
      call umalloc(ip_iiint2,nsub10,1,icall,memuse)
      call umalloc(ip_ifrom,nsub10,0,icall,memuse)
      call umalloc(ip_isav_pat,17*intmax0,1,icall,memuse)
      call umalloc(ip_isav_pat_b,6*intmax0*nsub10,1,icall,memuse)
      call umalloc(ip_mblk2nd,maxbl0,1,icall,memuse)
c
      open (unit=66,file='preronnie.out',form='formatted',
     .      status='unknown')
c
      read(iunit5,*)
c
      read(iunit5,'(a60)')grid
c
      read(iunit5,'(a60)')output
c
      read(iunit5,'(a60)')patch
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
      do nn=1,maxbl0
         mblk2nd(nn) = myhost
      end do
c
c     initialize output buffers
c
      do ll=1,nbuf0
         nou(ll) = 0
         do mm=1,ibufdim0
            bou(mm,ll) = ' '
         end do
      end do
c
c     default to dimensions of 1
c
      lmaxgr   = 1
      lmaxbl   = 1
      lmxseg   = 1
      lmaxcs   = 1
      lnplts   = 1
      lmxbli   = 1
      lmaxxe   = 1
      lnsub1   = 1
      lintmx   = 1
      lmxxe    = 1
      liitot   = 1
      isum     = 1
      lncycm   = 1
      isum_h   = 1
      isum_n   = 1
      lminnode = 1
      isumi    = 1
      isumi_h  = 1
      isumi_n  = 1
      lmptch   = 1
      lmsub1   = 1
      lintmax  = 1
      libufdim = 1
      lnbuf    = 1
      llbcprd  = 1
      llbcemb  = 1
      llbcrad  = 1
c
c     output banner
c
      write(66,83)
      write(66,83)
      write(66,87)
      write(66,9900)
 9900 format(2(2h *),47h                PRERONNIE - RONNIE PREPROCESSOR,
     .11h           ,4x,2(2h *))
      write(66,87)
      write(66,9990)
 9990 format(2(2h *),43h   VERSION 6.0 :  Computational Fluids Lab,,
     .15h Mail Stop 128,,4x,2(2h *),
     ./2(2h *),18x,41hNASA Langley Research Center, Hampton, VA,
     .3x,2(2h *),/2(2h *),18x,33hRelease Date:    August  8, 2001.,
     .11x,2(2h *))
      write(66,87)
      write(66,83)
      write(66,83)
   83 format(35(2h *))
   87 format(2(2h *),62x,2(2h *))
c
#ifdef CRAY_TIME
c     cray_time implies cray (always double precision)
      write(66,12) float(memuse)/1.e6
#else
#   ifdef DBLE_PRECSN
      write(66,12) float(memuse)/1.e6
#   else
      write(66,13) float(memuse)/1.e6
#   endif
#endif
   12 format(/,' memory allocation: ',f12.6,' Mbytes, double precision')
   13 format(/,' memory allocation: ',f12.6,' Mbytes, single precision')
c
      write(66,88)
   88 format(/19hinput/output files:)
c
      write(66,'(''  '',a60)')grid
      write(66,'(''  '',a60)')output
      write(66,'(''  '',a60)')patch
c
      read(iunit5,*)
      read(iunit5,*) ioflag,itrace
      write(66,398) ioflag,itrace
  398 format(/,4x,6hioflag,4x,6hitrace/,2i10)
c
      read(iunit5,10)(title(i),i=1,20)
   10 format(20a4)
      write(66,111)
  111 format(/5htitle)
      write(66,11)(title(i),i=1,20)
   11 format(2h  ,20a4)
c
      read(iunit5,*)
      read(iunit5,*) ngrid
      write(66,1638)
 1638 format(/,15hgrid/level data)
      write(66,1639)
 1639 format(1x,5hngrid)
      write(66,36) ngrid
   36 format(13i6)
c
      ip3dgrd = 0
      if (ngrid.lt.0) then
         ip3dgrd = 1
         ngrid =  iabs(ngrid)
      end if
      nchk = maxgr0-ngrid
c
c     check maximum number of grids
c
      if (nchk.lt.0) then
         write(66,1492)
         call termn8(myid,ierrflg,ibufdim0,nbuf0,bou,nou) 
      end if
 1492 format(55h stopping - insufficient maximum number of grids(maxgr))
c
      read(iunit5,*)
      nbl = 0
      iemtot = 0
      write(66,1631)
 1631 format(3x,3hncg,3x,3hiem,2x,4hidim,2x,4hjdim,2x,4hkdim)
      do 7001 igrid=1,ngrid
      nbl = nbl+1
      read(iunit5,*) ncg,iem,idim,jdim,kdim
      write(66,36) ncg,iem,idim,jdim,kdim
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
            write(66,*)' embedded grids must have ncg = 0'
            call termn8(myid,ierrflg,ibufdim0,nbuf0,bou,nou)
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
           write(66,'('' Cannot create coarser level for idim past'',
     .      i6)') idimg(nbl-1)
           istop=1
         end if
         if (float(jdimg(nbl-1)/2) .eq. float(jdimg(nbl-1))/2.) then
           write(66,'('' Cannot create coarser level for jdim past'',
     .      i6)') jdimg(nbl-1)
           istop=1
         end if
         if (float(kdimg(nbl-1)/2) .eq. float(kdimg(nbl-1))/2.) then
           write(66,'('' Cannot create coarser level for kdim past'',
     .      i6)') kdimg(nbl-1)
           istop=1
         end if
         if (istop .eq. 1) then
            call termn8(myid,ierrflg,ibufdim0,nbuf0,bou,nou)
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
      nchk   = maxbl0-nblock
      if (nchk.lt.0) then
         write(66,1649) maxbl0,nchk
         call termn8(myid,ierrflg,ibufdim0,nbuf0,bou,nou)
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
         write(66,157)m,levelt(m),levelb(m),ncgmax
         call termn8(myid,ierrflg,ibufdim0,nbuf0,bou,nou)
      end if
  157 format(1x,42herror in input, m, levelt, levelb, ncgmax=,4i5)
   17 continue
c
      lmaxgr  = ngrid
      lmaxbl  = nblock
c
      icall1  = 0
      iunit1 = 66
      imode1 = 0
      call global2(maxbl0,maxgr0,nsub10,ninter0,intmax0,ngrid,idimg,
     .             jdimg,kdimg,levelg,ncgg,nblg,iindex,llimit,
     .             iitmax,mmcxie,mmceta,ncheck,iifit,iic0,
     .             iiorph,iitoss,ifiner,dx,dy,dz,dthetx,
     .             dthety,dthetz,myid,mpatch0,maxxe0,icall1,iunit1,
     .             nou,bou,ibufdim0,nbuf0,ifrom,xif1,etf1,xif2,
     .             etf2,igridg,iemg,nblock,ioflag,imode1)
c
c     work array sizing (note: a parallel implementation would
c     also need some work space for message passing
c
      nstart = 1
      do nbl=1,nblock
         j = jdimg(nbl)
         k = kdimg(nbl)
         i = idimg(nbl)
         ns = j*k*i
c        x,y,z storage locations
         lw = nstart
         lw = lw + ns
         lw = lw + ns
         nstart = lw + ns
      end do
c
c     set and print out the necessary parameter sizes
c
      nwork  = nstart
c     need to add 1 for dummy storage of message-passing requirements
      nwork   = nwork + 1
      maxgr   = lmaxgr
      maxbl   = lmaxbl
      maxxe   = lmaxxe
      nsub1   = lnsub1
      mpatch  = lmptch
      intmax  = lintmax
      nbuf    = nbuf0
      ibufdim = ibufdim0
c
      write(66,*)
      write(66,*)'***********************************************'
      write(66,*)
      write(66,*)'    PARAMETER SIZES REQUIRED FOR THIS CASE:'
      write(66,*)
      write(66,*)'***********************************************'
      write(66,*)
      write(66,'('' nwork   = '',i10)') nwork
      write(66,'('' intmax  = '',i10)') intmax
      write(66,'('' maxxe   = '',i10)') maxxe
      write(66,'('' nsub1   = '',i10)') nsub1
      write(66,'('' mpatch  = '',i10)') mpatch
      write(66,'('' maxbl   = '',i10)') maxbl
      write(66,'('' maxgr   = '',i10)') maxgr
      write(66,'('' nbuf    = '',i10)') nbuf
      write(66,'('' ibufdim = '',i10)') ibufdim
c
      nstart = 1
      do nbl=1,nblock
         j = jdimg(nbl)
         k = kdimg(nbl)
         i = idimg(nbl)
         ns = j*k*i
c        x,y,z storage locations
         lw = nstart
         lw = lw + ns
         lw = lw + ns
         nstart = lw + ns
      end do

      write(66,*)'  '
      write(66,*)
     .'*******************************************************'
      write(66,*)'  '
      write(66,*)'            SUMMARY OF STORAGE REQUIREMENTS '
      write(66,*)'  '
      write(66,790) int(nstart)
  790 format('      permanent array w  requires ',i9,' (words)')
c
 800  format('          >>> Estimate for nwork = ',i9,' <<<')
c
      write(66,*)'  '
      write(66,*)'  '
      write(66,800) nstart
      write(66,*)'  '
      write(66,*)
     .'*******************************************************'
      write(66,*)'  '
c
      rewind (iunit5)
c
c     free up memory used by sizing routine
c
      ifree = 1
      if (ifree.gt.0) then
         call ufree(ip_mglevg)
         call ufree(ip_nemgl)
         call ufree(ip_levelg)
         call ufree(ip_nblg)
         call ufree(ip_iemg)
         call ufree(ip_igridg)
         call ufree(ip_jdimg)
         call ufree(ip_kdimg)
         call ufree(ip_idimg)
         call ufree(ip_nblcg)
         call ufree(ip_ncgg)
         call ufree(ip_iindex)
         call ufree(ip_llimit)
         call ufree(ip_iitmax)
         call ufree(ip_mmcxie)
         call ufree(ip_mmceta)
         call ufree(ip_ncheck)
         call ufree(ip_iifit)
         call ufree(ip_iic0)
         call ufree(ip_iiorph)
         call ufree(ip_iitoss)
         call ufree(ip_ifiner)
         call ufree(ip_dx)
         call ufree(ip_dy)
         call ufree(ip_dz)
         call ufree(ip_dthetx)
         call ufree(ip_dthety)
         call ufree(ip_dthetz)
         call ufree(ip_xif1)
         call ufree(ip_xif2)
         call ufree(ip_etf1)
         call ufree(ip_etf2)
         call ufree(ip_jjmax1)
         call ufree(ip_kkmax1)
         call ufree(ip_iiint1)
         call ufree(ip_iiint2)
         call ufree(ip_ifrom)
         call ufree(ip_isav_pat)
         call ufree(ip_isav_pat_b)
         call ufree(ip_mblk2nd)
      end if
c
      if (imode .eq. 1) then
         write(66,'(/,'' memory for preronnie has been deallocated'')')
      else
         write(6,'(/,''preronnie has completed successfully'')')
         write(6,'(''preronnie information has been put in'',
     .  '' file preronnie.out'',/)')
      end if
      close(66)
c
      return
      end
