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
      integer lout(nsub1),xif1(nsub1),xif2(nsub1),etf1(nsub1),
     .        etf2(nsub1)
c
      dimension work(nwork),levelg(maxbl)
      dimension mglevg(maxbl),nemgl(maxbl)
      dimension lw(65,maxbl),lw2(43,maxbl),igridg(maxbl),jdimg(maxbl),
     .          kdimg(maxbl),idimg(maxbl),nblcg(maxbl)
      dimension jtest(maxgr),ktest(maxgr),itest(maxgr)
      dimension ncgg(maxgr),nblg(maxgr),iemg(maxgr)
      dimension nou(nbuf)
      dimension xorig(maxbl),yorig(maxbl),zorig(maxbl)
      dimension xte(mpatch+2,mpatch+2,nsub1),
     .          yte(mpatch+2,mpatch+2,nsub1),
     .          zte(mpatch+2,mpatch+2,nsub1)
      dimension xmi(mpatch+2,mpatch+2,nsub1),
     .          ymi(mpatch+2,mpatch+2,nsub1),
     .          zmi(mpatch+2,mpatch+2,nsub1)
      dimension xmie(mpatch+2,mpatch+2,nsub1),
     .          ymie(mpatch+2,mpatch+2,nsub1),
     .          zmie(mpatch+2,mpatch+2,nsub1)
      dimension jjmax1(nsub1),kkmax1(nsub1)
      dimension iiint1(nsub1),iiint2(nsub1)
      dimension jimage(nsub1,mpatch+2,mpatch+2),
     .          kimage(nsub1,mpatch+2,mpatch+2)
      dimension jte(nsub1),kte(nsub1)
      dimension sxie(mpatch+2,mpatch+2,nsub1),
     .          seta(mpatch+2,mpatch+2,nsub1),
     .          sxie2(mpatch+2,mpatch+2),
     .          seta2(mpatch+2,mpatch+2)
      dimension xie2s(mpatch+2,mpatch+2),eta2s(mpatch+2,mpatch+2)
      dimension temp((mpatch+2)*(mpatch+2))
      dimension x2(mpatch+2,mpatch+2),y2(mpatch+2,mpatch+2),
     .          z2(mpatch+2,mpatch+2)
      dimension nblk1(mpatch+2),nblk2(mpatch+2)
      dimension jmm(mpatch+2),kmm(mpatch+2)
      dimension x1(mpatch+2,mpatch+2),y1(mpatch+2,mpatch+2),
     .          z1(mpatch+2,mpatch+2)
      dimension windex(maxxe,2),iindex(intmax,6*nsub1+9),
     .          llimit(intmax),iitmax(intmax),mmcxie(intmax),
     .          mmceta(intmax),ncheck(maxbl),iifit(intmax),
     .          nblkpt(maxxe),iic0(intmax),iiorph(intmax),
     .          iitoss(intmax),ifiner(intmax)
      dimension factjlo(intmax,nsub1),factjhi(intmax,nsub1),
     .          factklo(intmax,nsub1),factkhi(intmax,nsub1),
     .          ifrom(nsub1)
      dimension dx(intmax,nsub1),dy(intmax,nsub1),dz(intmax,nsub1),
     .          dthetx(intmax,nsub1),dthety(intmax,nsub1),
     .          dthetz(intmax,nsub1)
      dimension isav_pat(intmax,17),isav_pat_b(intmax,nsub1,6)
      dimension ireq_ar(intmax*3),mblk2nd(maxbl)
      dimension is_pat(5),ie_pat(5)
c
      pointer (ip_levelg,levelg),(ip_mglevg,mglevg),(ip_nemgl,nemgl)
      pointer (ip_work,work),(ip_lw,lw),(ip_lw2,lw2)
      pointer (ip_xorig,xorig),(ip_yorig,yorig),(ip_zorig,zorig)
      pointer (ip_jtest,jtest),(ip_ktest,ktest),(ip_itest,itest)
      pointer (ip_ncgg,ncgg),(ip_nblg,nblg),(ip_iemg,iemg)
      pointer (ip_igridg,igridg),(ip_nblcg,nblcg),(ip_ireq_ar,ireq_ar)
      pointer (ip_windex,windex),(ip_iindex,iindex),
     .        (ip_llimit,llimit),(ip_iitmax,iitmax),(ip_mmcxie,mmcxie),
     .        (ip_mmceta,mmceta),(ip_ncheck,ncheck),(ip_iifit,iifit),
     .        (ip_nblkpt,nblkpt),(ip_iic0,iic0),(ip_iiorph,iiorph),
     .        (ip_iitoss,iitoss),(ip_ifiner,ifiner)
      pointer (ip_dx,dx),(ip_dy,dy),(ip_dz,dz),
     .        (ip_dthetx,dthetx),(ip_dthety,dthety),
     .        (ip_dthetz,dthetz)
      pointer (ip_lout,lout),(ip_xif1,xif1),(ip_xif2,xif2),
     .        (ip_etf1,etf1),(ip_etf2,etf2),(ip_iiint1,iiint1),
     .        (ip_iiint2,iiint2),(ip_jjmax1,jjmax1),
     .        (ip_kkmax1,kkmax1),(ip_jimage,jimage),
     .        (ip_kimage,kimage),(ip_jte,jte),(ip_kte,kte),
     .        (ip_jmm,jmm),(ip_kmm,kmm),(ip_nblk1,nblk1),
     .        (ip_nblk2,nblk2),(ip_xte,xte),(ip_yte,yte),
     .        (ip_zte,zte),(ip_xmi,xmi),
     .        (ip_ymi,ymi),(ip_zmi,zmi),
     .        (ip_xmie,xmie),(ip_ymie,ymie),
     .        (ip_zmie,zmie),(ip_sxie,sxie),
     .        (ip_seta,seta),(ip_sxie2,sxie2),
     .        (ip_seta2,seta2),(ip_xie2s,xie2s),
     .        (ip_eta2s,eta2s),(ip_temp,temp),
     .        (ip_x2,x2),(ip_y2,y2),
     .        (ip_z2,z2),(ip_x1,x1),
     .        (ip_y1,y1),(ip_z1,z1),
     .        (ip_factjlo,factjlo),(ip_factjhi,factjhi),
     .        (ip_factklo,factklo),(ip_factkhi,factkhi),
     .        (ip_ifrom,ifrom)
      pointer (ip_isav_pat,isav_pat),(ip_isav_pat_b,isav_pat_b)
      pointer (ip_idimg,idimg),(ip_jdimg,jdimg),(ip_kdimg,kdimg)
      pointer (ip_mblk2nd,mblk2nd)
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
     .        lxold,lyold,lzold,lxorg,lyorg,lzorg
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
      icall  = 0
      memuse = 0
c
      call umalloc(ip_work,nwork,0,icall,memuse)
      call umalloc(ip_lw,65*maxbl,1,icall,memuse)
      call umalloc(ip_lw2,43*maxbl,1,icall,memuse)
      call umalloc(ip_mglevg,maxbl,1,icall,memuse)
      call umalloc(ip_nemgl,maxbl,1,icall,memuse)
      call umalloc(ip_levelg,maxbl,1,icall,memuse)
      call umalloc(ip_jtest,maxgr,1,icall,memuse)
      call umalloc(ip_ktest,maxgr,1,icall,memuse)
      call umalloc(ip_itest,maxgr,1,icall,memuse)
      call umalloc(ip_nblg,maxgr,1,icall,memuse)
      call umalloc(ip_iemg,maxgr,1,icall,memuse)
      call umalloc(ip_igridg,maxbl,1,icall,memuse)
      call umalloc(ip_jdimg,maxbl,1,icall,memuse)
      call umalloc(ip_kdimg,maxbl,1,icall,memuse)
      call umalloc(ip_idimg,maxbl,1,icall,memuse)
      call umalloc(ip_nblcg,maxbl,1,icall,memuse)
      call umalloc(ip_ncgg,maxgr,1,icall,memuse)
      call umalloc(ip_xorig,maxbl,0,icall,memuse)
      call umalloc(ip_yorig,maxbl,0,icall,memuse)
      call umalloc(ip_zorig,maxbl,0,icall,memuse)
      call umalloc(ip_windex,2*maxxe,0,icall,memuse)
      call umalloc(ip_iindex,intmax*(6*nsub1+9),1,icall,memuse)
      call umalloc(ip_llimit,intmax,1,icall,memuse)
      call umalloc(ip_iitmax,intmax,1,icall,memuse)
      call umalloc(ip_mmcxie,intmax,1,icall,memuse)
      call umalloc(ip_mmceta,intmax,1,icall,memuse)
      call umalloc(ip_ncheck,maxbl,1,icall,memuse)
      call umalloc(ip_iifit,intmax,1,icall,memuse)
      call umalloc(ip_nblkpt,maxxe,1,icall,memuse)
      call umalloc(ip_iic0,intmax,1,icall,memuse)
      call umalloc(ip_iiorph,intmax,1,icall,memuse)
      call umalloc(ip_iitoss,intmax,1,icall,memuse)
      call umalloc(ip_ifiner,intmax,1,icall,memuse)
      call umalloc(ip_dx,intmax*nsub1,0,icall,memuse)
      call umalloc(ip_dy,intmax*nsub1,0,icall,memuse)
      call umalloc(ip_dz,intmax*nsub1,0,icall,memuse)
      call umalloc(ip_dthetx,intmax*nsub1,0,icall,memuse)
      call umalloc(ip_dthety,intmax*nsub1,0,icall,memuse)
      call umalloc(ip_dthetz,intmax*nsub1,0,icall,memuse)
      call umalloc(ip_lout,nsub1,1,icall,memuse)
      call umalloc(ip_xif1,nsub1,1,icall,memuse)
      call umalloc(ip_xif2,nsub1,1,icall,memuse)
      call umalloc(ip_etf1,nsub1,1,icall,memuse)
      call umalloc(ip_etf2,nsub1,1,icall,memuse)
      call umalloc(ip_jjmax1,nsub1,1,icall,memuse)
      call umalloc(ip_kkmax1,nsub1,1,icall,memuse)
      call umalloc(ip_iiint1,nsub1,1,icall,memuse)
      call umalloc(ip_iiint2,nsub1,1,icall,memuse)
      call umalloc(ip_jimage,nsub1*(mpatch+2)*(mpatch+2),1,icall,memuse)
      call umalloc(ip_kimage,nsub1*(mpatch+2)*(mpatch+2),1,icall,memuse)
      call umalloc(ip_jte,nsub1,1,icall,memuse)
      call umalloc(ip_kte,nsub1,1,icall,memuse)
      call umalloc(ip_jmm,mpatch+2,1,icall,memuse)
      call umalloc(ip_kmm,mpatch+2,1,icall,memuse)
      call umalloc(ip_nblk1,mpatch+2,1,icall,memuse)
      call umalloc(ip_nblk2,mpatch+2,1,icall,memuse)
      call umalloc(ip_xte,nsub1*(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_yte,nsub1*(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_zte,nsub1*(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_xmi,nsub1*(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_ymi,nsub1*(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_zmi,nsub1*(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_xmie,nsub1*(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_ymie,nsub1*(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_zmie,nsub1*(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_sxie,nsub1*(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_seta,nsub1*(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_sxie2,(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_seta2,(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_xie2s,(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_eta2s,(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_temp,(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_x2,(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_y2,(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_z2,(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_x1,(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_y1,(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_z1,(mpatch+2)*(mpatch+2),0,icall,memuse)
      call umalloc(ip_factjlo,intmax*nsub1,0,icall,memuse)
      call umalloc(ip_factjhi,intmax*nsub1,0,icall,memuse)
      call umalloc(ip_factklo,intmax*nsub1,0,icall,memuse)
      call umalloc(ip_factkhi,intmax*nsub1,0,icall,memuse)
      call umalloc(ip_ifrom,nsub1,0,icall,memuse)
      call umalloc(ip_isav_pat,17*intmax,1,icall,memuse)
      call umalloc(ip_isav_pat_b,6*intmax*nsub1,1,icall,memuse)
      call umalloc(ip_mblk2nd,maxbl,1,icall,memuse)
      call umalloc(ip_ireq_ar,intmax*3,1,icall,memuse)
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
 9990 format(2(2h *),43h   VERSION 6.0 :  Computational Fluids Lab,,
     .15h Mail Stop 128,,4x,2(2h *),
     ./2(2h *),18x,41hNASA Langley Research Center, Hampton, VA,
     .3x,2(2h *),/2(2h *),18x,33hRelease Date:    August 14, 2001.,
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
         call ufree(ip_work)
         call ufree(ip_lw)
         call ufree(ip_lw2)
         call ufree(ip_mglevg)
         call ufree(ip_nemgl)
         call ufree(ip_levelg)
         call ufree(ip_jtest)
         call ufree(ip_ktest)
         call ufree(ip_itest)
         call ufree(ip_nblg)
         call ufree(ip_iemg)
         call ufree(ip_igridg)
         call ufree(ip_jdimg)
         call ufree(ip_kdimg)
         call ufree(ip_idimg)
         call ufree(ip_nblcg)
         call ufree(ip_ncgg)
         call ufree(ip_xorig)
         call ufree(ip_yorig)
         call ufree(ip_zorig)
         call ufree(ip_windex)
         call ufree(ip_iindex)
         call ufree(ip_llimit)
         call ufree(ip_iitmax)
         call ufree(ip_mmcxie)
         call ufree(ip_mmceta)
         call ufree(ip_ncheck)
         call ufree(ip_iifit)
         call ufree(ip_nblkpt)
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
         call ufree(ip_lout)
         call ufree(ip_xif1)
         call ufree(ip_xif2)
         call ufree(ip_etf1)
         call ufree(ip_etf2)
         call ufree(ip_jjmax1)
         call ufree(ip_kkmax1)
         call ufree(ip_iiint1)
         call ufree(ip_iiint2)
         call ufree(ip_jimage)
         call ufree(ip_kimage)
         call ufree(ip_jte)
         call ufree(ip_kte)
         call ufree(ip_jmm)
         call ufree(ip_kmm)
         call ufree(ip_nblk1)
         call ufree(ip_nblk2)
         call ufree(ip_xte)
         call ufree(ip_yte)
         call ufree(ip_zte)
         call ufree(ip_xmi)
         call ufree(ip_ymi)
         call ufree(ip_zmi)
         call ufree(ip_xmie)
         call ufree(ip_ymie)
         call ufree(ip_zmie)
         call ufree(ip_sxie)
         call ufree(ip_seta)
         call ufree(ip_sxie2)
         call ufree(ip_seta2)
         call ufree(ip_xie2s)
         call ufree(ip_eta2s)
         call ufree(ip_temp)
         call ufree(ip_x2)
         call ufree(ip_y2)
         call ufree(ip_z2)
         call ufree(ip_x1)
         call ufree(ip_y1)
         call ufree(ip_z1)
         call ufree(ip_factjlo)
         call ufree(ip_factjhi)
         call ufree(ip_factklo)
         call ufree(ip_factkhi)
         call ufree(ip_ifrom)
         call ufree(ip_isav_pat)
         call ufree(ip_isav_pat_b)
         call ufree(ip_mblk2nd)
         call ufree(ip_ireq_ar)
      end if
c
      return
      end
