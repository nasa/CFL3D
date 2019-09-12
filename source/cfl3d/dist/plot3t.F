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
      subroutine plot3t(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,q,x,y,
     .                  z,xw,blank2,blank,xg,vist3d,iover,nblk,
     .                  nmap,sj,sk,si,smin,ux,turre,vol,qj0,qk0,qi0,
     .                  bcj,bck,bci,wt,cmuv,jdw,kdw,idw,
     .                  nplots,jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                  jeg,keg,ieg,ninter,iindex,intmax,nsub1,
     .                  maxxe,nblkk,nbli,limblk,isva,nblon,mxbli,
     .                  thetay,maxbl,maxgr,myid,myhost,mycomm,
     .                  mblk2nd,inpl3d,nblock,nblkpt,volj0,volk0,voli0,
     .                  vormax,ivmax,jvmax,kvmax,nummem,ifuncdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Write turbulent information for the cell-centered points
c     in PLOT3D format.
c
c     outputs grid/turb info in single precision for use with FAST/PLOT3D
c***********************************************************************
c
#if defined ADP_OFF
#   ifdef CMPLX
#     ifdef DBLE_PRECSN
      implicit complex*8(a-h,o-z)
#     else
      implicit complex(a-h,o-z)
#     endif
#   else
#     ifdef DBLE_PRECSN
      implicit real*8 (a-h,o-z)
#     endif
#   endif
#else
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
#endif
c
#if defined DIST_MPI
#     include "mpif.h"
      dimension istat(MPI_STATUS_SIZE)
#   ifdef P3D_SINGLE
#     define MY_MPI_REAL MPI_REAL
#   else
#     define MY_MPI_REAL MPI_DOUBLE_PRECISION
#   endif
#endif
c
#if defined P3D_SINGLE
      real*4    xw(jdw,kdw,idw,ifuncdim),xg(jdw,kdw,idw,4)
      real*4    xmachw,alphww,reuew,timew
#else
      real xw(jdw,kdw,idw,ifuncdim),xg(jdw,kdw,idw,4)
#endif
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension q(jdim,kdim,idim,5),vist3d(jdim,kdim,idim)
      dimension blank(jdim,kdim,idim),blank2(jdim,kdim,idim)
      dimension sk(jdim,kdim,idim-1,5),sj(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5)
      dimension smin(jdim-1,kdim-1,idim-1),ux(jdim-1,kdim-1,idim-1,9)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension qj0(kdim,idim-1,5,4),
     .          qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4)
      dimension wt(jdim,kdim,9),turre(jdim,kdim,idim,nummem)
      dimension cmuv(jdim-1,kdim-1,idim-1)
      dimension nmap(maxbl),jdimg(maxbl),kdimg(maxbl),idimg(maxbl),
     .          nblcg(maxbl),jsg(maxbl),ksg(maxbl),isg(maxbl),
     .          jeg(maxbl),keg(maxbl),ieg(maxbl),
     .          iindex(intmax,6*nsub1+9),nblkpt(maxxe),nblkk(2,mxbli),
     .          limblk(2,6,mxbli),isva(2,2,mxbli),nblon(mxbli),
     .          inpl3d(nplots,11),thetay(maxbl),mblk2nd(maxbl)
      dimension vol(jdim,kdim,idim-1),volj0(kdim,idim-1,4),
     .          volk0(jdim,idim-1,4),voli0(jdim,kdim,4)
      dimension vormax(maxbl),ivmax(maxbl),jvmax(maxbl),kvmax(maxbl)
c
      common /bin/ ibin,iblnk,iblnkfr,ip3dgrad
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /twod/ i2d
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /lam/ ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi,
     .        i_lam_forcezero
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /reyue/ reue,tinf,ivisc(3)
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /igrdtyp/ ip3dgrd,ialph
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
      common /conversion/ radtodeg
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /konew/ ikoprod,isstdenom,pklimterm,ibeta8kzeta,i_bsl,
     .        keepambient,re_thetat0,i_wilcox06,i_wilcox06_chiw,
     .        i_turbprod_kterm,i_catris_kw,prod2d3dtrace,
     .        i_compress_correct,isstsf,i_wilcox98,i_wilcox98_chiw,
     .        isst2003
      common /plot3dtyp/ ifunct
      common /constit/ i_nonlin,c_nonlin,snonlin_lim,i_tauijs,i_qcr2000,
     .                 i_qcr2013,i_qcr2013v
      common /easmv/ c10,c11,c2,c3,c4,c5,sigk1,cmuc1,ieasm_type
      common /writestuff/ ifort50write, j_ifort50write, i_ifort50write
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset = maxbl
      itag_grid = 1
      itag_q    = itag_grid + ioffset
#endif
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      re=reue/xmach
      ccr2=2.5
c
c
c     initialize xw and xg arrays
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
      do j=1,jw
         do k=1,kw
            do i=1,iw
               do l= 1,ifuncdim
               xw(j,k,i,l) = 0.
               end do
               do l= 1,4
               xg(j,k,i,l) = 0.
               end do
            end do
         end do
      end do
c
c     assign single precision scalars
c
      alphaw = radtodeg*(alpha+thetay(nblk))
      xmachw = xmach
      alphww = alphaw
      reuew  = reue
      timew  = time
c
c***********************************************************************
c     get ux values
c***********************************************************************
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
      call delv(jdim,kdim,idim,q,sj,sk,si,vol,ux,wt,
     .          blank,iover,qj0,qk0,qi0,bcj,bck,bci,nblk,
     .          volj0,volk0,voli0,maxbl,vormax,ivmax,jvmax,kvmax)
c
#if defined DIST_MPI
      end if
c
#endif
c***********************************************************************
c     average grid points to get cell centers and load iinto 1st 3 
c     locations of the single precision xg array
c***********************************************************************
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
      iw = 0
      do 2000 i=i1,i2,i3
      iw = iw + 1
      kw = 0
      do 2000 k=k1,k2,k3
      kw = kw + 1
      jw = 0
      do 2000 j=j1,j2,j3
      jw = jw + 1
      xg(jw,kw,iw,1) = .125*(   x(j,k,i)    +x(j,k,i+1)
     .                        + x(j+1,k,i)  +x(j+1,k,i+1)
     .                        + x(j,k+1,i)  +x(j,k+1,i+1)
     .                        + x(j+1,k+1,i)+x(j+1,k+1,i+1))
c
      xg(jw,kw,iw,2) = .125*(   y(j,k,i)    +y(j,k,i+1)
     .                        + y(j+1,k,i)  +y(j+1,k,i+1)
     .                        + y(j,k+1,i)  +y(j,k+1,i+1)
     .                        + y(j+1,k+1,i)+y(j+1,k+1,i+1))
c
      xg(jw,kw,iw,3) = .125*(   z(j,k,i)    +z(j,k,i+1)
     .                        + z(j+1,k,i)  +z(j+1,k,i+1)
     .                        + z(j,k+1,i)  +z(j,k+1,i+1)
     .                        + z(j+1,k+1,i)+z(j+1,k+1,i+1))
 2000 continue
c
#if defined DIST_MPI
      end if
c
#endif
c***********************************************************************
c     set iblank (blank2) array
c***********************************************************************
c
c     don't compute blanking info if iblnk = 0
c
      if (iblnk .gt. 0) then
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
c     assign default iblank (blank2) array
c
      do 2049 i=1,idim
      do 2049 j=1,jdim
      do 2049 k=1,kdim
      blank2(j,k,i) = 1.
 2049 continue
c
c     iblank (blank2) array for generalized patch interface boundaries 
c
      if (abs(ninter).gt.0) then
         do 1600 inter = 1,abs(ninter)
         lmax1  = iindex(inter,1)
         nbl    = iindex(inter,lmax1+2)
         if (nbl.ne.nblk) go to 1600
         lst    = iindex(inter,2*lmax1+5)
         lcoord = iindex(inter,2*lmax1+3)/10
         lend   = iindex(inter,2*lmax1+3)-lcoord*10
         j21    = iindex(inter,2*lmax1+6)
         j22    = iindex(inter,2*lmax1+7)
         k21    = iindex(inter,2*lmax1+8)
         k22    = iindex(inter,2*lmax1+9)
         if (lcoord.eq.1) then
           if (lend.eq.1) i = 1
           if (lend.eq.2) i = idim1
           do 1610 j = j21,j22-1
           do 1610 k = k21,k22-1
           ll = lst + (j22-j21)*(k-k21) + (j-j21)
           mblk         = iindex(inter,nblkpt(ll)+1)
           blank2(j,k,i) = -float(nmap(mblk))
 1610      continue
         end if
         if (lcoord.eq.2) then
           if (lend.eq.1) j = 1
           if (lend.eq.2) j = jdim1
           do 1620 i = k21,k22-1
           do 1620 k = j21,j22-1
           ll = lst + (j22-j21)*(i-k21) + (k-j21)
           mblk         = iindex(inter,nblkpt(ll)+1)
           blank2(j,k,i) = -float(nmap(mblk))
 1620      continue
         end if
         if (lcoord.eq.3) then
            if (lend.eq.1) k = 1
            if (lend.eq.2) k = kdim1
            do 1630 i = k21,k22-1
            do 1630 j = j21,j22-1
            ll = lst + (j22-j21)*(i-k21) + (j-j21)
            mblk         = iindex(inter,nblkpt(ll)+1)
            blank2(j,k,i) = -float(nmap(mblk))
 1630       continue
         end if
 1600    continue
      end if
c
c     iblank (blank2) array for 1:1 interface boundaries
c
      if(nbli.gt.0) then
        do 100 n=1,abs(nbli)
        if(nblon(n).ge.0) then
          if(nblk.eq.nblkk(1,n) .or. nblk.eq.nblkk(2,n)) then
            it = 1
            ir = 2
            if(nblk.eq.nblkk(2,n)) then
              it = 2
              ir = 1
            end if
c
c           allow for 1-1 blocking in same grid
c
            itime = 1
            if (nblkk(1,n).eq.nblkk(2,n)) itime = 2
            do 101 iti = 1,itime
            if (iti.gt.1) then
               it = 1
               ir = 2
            end if
c
            is = limblk(it,1,n)
            ie = limblk(it,4,n)
            js = limblk(it,2,n)
            je = limblk(it,5,n)
            ks = limblk(it,3,n)
            ke = limblk(it,6,n)
c
            is1 = min(is,ie)
            ie1 = max(is,ie)
            js1 = min(js,je)
            je1 = max(js,je)
            ks1 = min(ks,ke)
            ke1 = max(ks,ke)
c
            ie1 = min(ie1,idim1)
            je1 = min(je1,jdim1)
            ke1 = min(ke1,kdim1)
            is1 = min(is1,idim1)
            js1 = min(js1,jdim1)
            ks1 = min(ks1,kdim1)
c
            do 110 i = is1,ie1
            do 110 j = js1,je1
            do 110 k = ks1,ke1
            blank2(j,k,i) = -float(nmap(nblkk(ir,n)))
  110       continue
  101       continue
          end if
        end if
  100   continue
      end if
c
c     iblank (blank2) array for embedded grids - the underlying
c     coarse grid areas are blanked out if the parameter ibembed > 0
c
      ibembed = 1
c
      if (ibembed.gt.0) then
         do 7500 nblc=1,nblock
         if (nblk.eq.nblc) go to 7500
         nblcc    = nblcg(nblc)
         if (nblcc.eq.nblk) then
            js       = jsg(nblc)
            ks       = ksg(nblc)
            is       = isg(nblc)
            je       = jeg(nblc)
            ke       = keg(nblc)
            ie       = ieg(nblc)
            if (je.gt.js) je = je - 1
            if (ke.gt.ks) ke = ke - 1
            if (ie.gt.is) ie = ie - 1
            do 7501 i=is,ie
            do 7501 j=js,je
            do 7501 k=ks,ke
            blank2(j,k,i) = 0.
 7501       continue
         end if
 7500    continue
      end if
c
c     iblank (blank2) array for overlapped grids
c
      if (iover.eq.1) then
         do 413 i=1,idim1
         do 413 j=1,jdim1
         do 413 k=1,kdim1
         if (real(blank(j,k,i)).eq.0.0) blank2(j,k,i) = 0.
  413    continue
      end if
c
#if defined DIST_MPI
      end if
c
#endif
      end if
c
c***********************************************************************
c      plot3d data
c***********************************************************************
c
      if (lhdr.gt.0  .and. myid.eq.myhost) then
         if (i2d .eq. 1) then
            if (iblnk .eq. 0) then
               write(11,'(''writing turb plot3d file for JDIM X'',
     .         '' KDIM ='',i5,'' x '',i5,'' grid'')') jdim,kdim
               write(11,'(''   plot3d files to be read with 2d/mgrid'',
     .         '' qualifiers'')')
            else
               write(11,'(''writing turb plot3d file for JDIM X'',
     .         '' KDIM ='',i5,'' x '',i5,'' grid'')') jdim,kdim
               write(11,'(''   plot3d files to be read with 2d/mgrid'',
     .         ''/blank qualifiers'')')
            end if
         else
            write(11,93)idim,jdim,kdim
   93       format(48h writing turb plot3d file for IDIM X JDIM X KDIM,
     .      2h =,i5,3h x ,i5,3h x ,i5,5h grid)
            if (iblnk .eq. 0) then
               write(11,2042)
 2042          format(3x,35hplot3d files to be read with /mgrid,
     .         11h qualifiers)
            else
               write(11,2142)
 2142          format(3x,41hplot3d files to be read with /mgrid/blank,
     .         11h qualifiers)
            end if
         end if
         if (ifunct .eq. 0) then
           write(11,'(''   plot3d file is Q file of turbulence'',
     .      '' quantities'')')
           if (i2d .eq. 1) then
              write(11,
     .         '(''  Defaults are production-term, uw, uu, ww'',
     .         '' for 2D'')') 
           else
              write(11,
     .         '(''   Defaults are production-term, uw, uu, ww,'',
     .         '' Sk/eps for 3D'')')
           end if
         else
           write(11,'(''   plot3d file is function file of'',i5,
     .      '' turbulence quantities'')') ifunct
         end if
#        ifdef CMPLX
         if (ip3dgrad .ne. 0) then
            write(11,'(''   NOTE: contains d(turb data)/d(), rather '',
     .                 ''than turb data'')')
            call gradinfo(delh,11)
         end if
#        endif
      end if
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
c     load blank2 into 4th location of the single precision xg array
c
      iw = 0
      do 9000 i=i1,i2,i3
      iw = iw + 1
      jw = 0
      do 9000 j=j1,j2,j3
      jw = jw + 1
      kw = 0
      do 9000 k=k1,k2,k3
      kw = kw + 1
      xg(jw,kw,iw,4) = blank2(j,k,i)
 9000 continue
#if defined DIST_MPI
c
      end if
#endif
#if defined DIST_MPI
c
c     send/receive plot3d grid data
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
c
      nng = jw*kw*iw*4
      nd_srce = mblk2nd(nblk)
      mytag = itag_grid + nblk
c
      if (mblk2nd(nblk).eq.myid) then
         call MPI_Send (xg, nng, MY_MPI_REAL, myhost, mytag,
     &                  mycomm, ierr)
      else if (myid.eq.myhost) then
         call MPI_Recv (xg, nng, MY_MPI_REAL, nd_srce,
     &                  mytag, mycomm, istat, ierr)
      end if
c
#endif
c
      if (myid.eq.myhost .and. iblnk.gt.0) then
c
c     ialph > 0 for a grid that was read in plot3d format with alpha measured
c               in the xy plane (TLNS3D convention)
c
      if (ialph.ne.0 .and. i2d.ne.1) then
         do i=1,iw
            do j=1,jw
               do k=1,kw
                  temp        = xg(j,k,i,2)
                  xg(j,k,i,2) = xg(j,k,i,3)
                  xg(j,k,i,3) = -temp
               end do
            end do
         end do
      end if
c
c     output grid
c
      if (ibin.eq.0) then
         if (i2d.eq.0) then
            if (iblnk .eq. 0) then
               write(3,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(3,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
               write(3,'(10i8)')
     .               (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         else
            if (iblnk .eq. 0) then
               write(3,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(3,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
               write(3,'(10i8)')
     .               (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         end if
      else
         if (i2d.eq.0) then
            if (iblnk .eq. 0) then
               write(3)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(3)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .             (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         else
            if (iblnk .eq. 0) then
               write(3)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(3)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .            (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         end if
      end if
c
      end if
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
c     turbulent quantities
c
c     There are a variety of turbulent quantities that can be computed and printed
c     out.  By default, only 5 (or 4, for 2-D) can be printed out at a time,
c     into a PLOT3D Q-type file.  However, if the keyword ifunct is employed, 
c     a function file is written instead, with "ifunct" number of variables.  
c     Modify the section that assigns quantities to the "xw" arrays in order to 
c     get what you want.  (See "Modify Here" section below.)
c 
      uprime=0.
      vprime=0.
      wprime=0.
      uuprime=0.
      vvprime=0.
      wwprime=0.
      uvprime=0.
      uwprime=0.
      vwprime=0.
      ske=0.
      prod=0.
      fd=0.
      povere=0.
      turbre=0.
      def=0.
      tint=0.
      b11=0.
      b22=0.
      b33=0.
      b12=0.
      b13=0.
      b23=0.
      xii=0.
      xiii=0.
      ff=0.
      zk=0.
      qinvar=0.
      uuu=0.
      vvv=0.
      www=0.
      sigma=0.
      r11slow=0.
      r22slow=0.
      r33slow=0.
      r12slow=0.
      r13slow=0.
      r23slow=0.
      r11rapid=0.
      r22rapid=0.
      r33rapid=0.
      r12rapid=0.
      r13rapid=0.
      r23rapid=0.
      r11=0.
      r22=0.
      r33=0.
      r12=0.
      r13=0.
      r23=0.
c   optional printout to fort.50 (specifically for k-line; k=1 must be body)
c   writes "plus" turbulent values; to be used only for single block runs
c   (if multiple blocks, you will only get the last one (others get overwritten))
      if (ifort50write .eq. 1) then
      write(50,'(''variables="y","uuplus","wwplus","vvplus","uwplus",'',
     + ''"uvplus","vwplus","u_scaled","kplus","yplus","uplus",'',
     + ''"mu_turbulent/mu_lam_ref","eplus","cmu"'')')
      write(50,'(''#'',i5,''  CFL3D scaled turb results along k-line'',
     + '' at j,i='',2i4, '' x='',f16.6)') k2-k1+1,j_ifort50write,
     +  i_ifort50write,0.5*(x(j_ifort50write,1,i_ifort50write)+
     +  x(j_ifort50write+1,1,i_ifort50write))
      write(50,'(''# NOTE: eplus applies to k-eps models only'')')
      end if
c
      iw = 0
      do 4000 i=i1,i2,i3
      iw = iw + 1
      kw = 0
      do 4000 k=k1,k2,k3
      kw = kw + 1
      jw = 0
      do 4000 j=j1,j2,j3
      jw = jw + 1
c
           if (i.ge.ilamlo .and. i.lt.ilamhi .and.
     .         j.ge.jlamlo .and. j.lt.jlamhi .and.
     .         k.ge.klamlo .and. k.lt.klamhi) then
             cutoff=0.
           else
             cutoff=1.
           end if
c          Determine Sij and Wij values:
            s11 = ux(j,k,i,1)
            s22 = ux(j,k,i,5)
            s33 = ux(j,k,i,9)
            s11star = s11 - ((s11+s22+s33)/3.)
            s22star = s22 - ((s11+s22+s33)/3.)
            s33star = s33 - ((s11+s22+s33)/3.)
            s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
            s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
            s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
            w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
            w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
            w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
c           xis=SijSij, wis=WijWij
            xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
            wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
            c1=2.*c10
c           c1=6.8
c           c2=0.36
c           c3=1.25
c           c4=0.40
c           c5=1.88
            gg=1./(c1/2.+c5-1.)
            if (ivmx .lt. 6) then
c   The following gets an estimate for turb kinetic energy for
c   models that do not compute it directly, in case it is needed
c   (from AIAA 92-0439, where structure parameter taken as a1=0.5*(0.31))
              zk=vist3d(j,k,i)/q(j,k,i,1)*sqrt(2.*xis)/(re*.31)
            else
c   Otherwise, use exact computed turb kinetic energy
              zk=turre(j,k,i,2)
            end if
c   If using QCR2013, then k is already added, so do not add here:
            if (i_qcr2013 .eq. 1 .or. i_qcr2013v .eq. 1) then
              zk=0.
            end if
c
            if ((ivmx .le. 7 .or. ivmx .eq. 10 .or. ivmx .eq. 15
     +          .or. ivmx .eq. 16)
     +          .and. i_nonlin .eq. 0) then
              t11 = -2.*vist3d(j,k,i)*s11star
              t22 = -2.*vist3d(j,k,i)*s22star
              t33 = -2.*vist3d(j,k,i)*s33star
c   Add 2/3*k to the diagonal stress terms for eddy viscosity
c   models, for better approximation to normal stresses
              t11=t11+2.*q(j,k,i,1)*zk*re/3.
              t22=t22+2.*q(j,k,i,1)*zk*re/3.
              t33=t33+2.*q(j,k,i,1)*zk*re/3.
c
              t12 = -2.*vist3d(j,k,i)*s12
              t13 = -2.*vist3d(j,k,i)*s13
              t23 = -2.*vist3d(j,k,i)*s23
              if (i_qcr2000 .eq. 1 .or. i_qcr2013 .eq. 1 .or.
     +            i_qcr2013v .eq. 1) then
                denom=sqrt(ux(j,k,i,1)**2 + ux(j,k,i,2)**2 + 
     +                     ux(j,k,i,3)**2 + ux(j,k,i,4)**2 +
     +                     ux(j,k,i,5)**2 + ux(j,k,i,6)**2 +
     +                     ux(j,k,i,7)**2 + ux(j,k,i,8)**2 +
     +                     ux(j,k,i,9)**2)+1.e-20
                o11=0.
                o22=0.
                o33=0.
                o12=2.*w12/denom
                o13=2.*w13/denom
                o23=2.*w23/denom
                t11p=t11-0.3*( o11*t11+o12*t12+o13*t13+
     +                         o11*t11+o12*t12+o13*t13)
                t22p=t22-0.3*(-o12*t12+o22*t22+o23*t23+
     +                        -o12*t12+o22*t22+o23*t23)
                t33p=t33-0.3*(-o13*t13-o23*t23+o33*t33+
     +                        -o13*t13-o23*t23+o33*t33)
                t12p=t12-0.3*( o11*t12+o12*t22+o13*t23+
     +                        -o12*t11+o22*t12+o23*t13)
                t13p=t13-0.3*( o11*t13+o12*t23+o13*t33+
     +                        -o13*t11-o23*t12+o33*t13)
                t23p=t23-0.3*(-o12*t13+o22*t23+o23*t33+
     +                        -o13*t12-o23*t22+o33*t23)
                t11=t11p
                t22=t22p
                t33=t33p
                t12=t12p
                t13=t13p
                t23=t23p
              end if
              if (i_qcr2013 .eq. 1) then
                xis2013 = s11star*s11star + s22star*s22star + 
     +                    s33star*s33star +
     +                    2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
                t11=t11+ccr2*vist3d(j,k,i)*sqrt(2.*xis2013)
                t22=t22+ccr2*vist3d(j,k,i)*sqrt(2.*xis2013)
                t33=t33+ccr2*vist3d(j,k,i)*sqrt(2.*xis2013)
              end if
              if (i_qcr2013v .eq. 1) then
                t11=t11+ccr2*vist3d(j,k,i)*sqrt(2.*wis)
                t22=t22+ccr2*vist3d(j,k,i)*sqrt(2.*wis)
                t33=t33+ccr2*vist3d(j,k,i)*sqrt(2.*wis)
              end if
            else if(ivmx .eq. 12) then
c            Find nonlinear tauij values:
                alpa1=(2.-c4)/2.*gg
                alpa2=(2.-c3)*gg
              omegatemp=ccmax(turre(j,k,i,1),tur10)
c     Get "factre" term:  changes numerator of cmu for nonlinear terms only
              eta=(2.-c3)**2*(gg*gg/4.)*xis/(omegatemp*re)**2
              squig=(2.-c4)**2*(gg*gg/4.)*wis/(omegatemp*re)**2
              eta=ccmincr(eta,10.)
              squig=ccmincr(squig,10.)
              factre=(3.*(1.+eta)+.2e-8*(eta*eta*eta +
     +          squig*squig*squig))/
     +               (3.*(1.+eta)+   .2*(eta*eta*eta +
     +          squig*squig*squig))
              t11 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre/omegatemp*
     +                 (-s12*w12 - s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              t22 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre/omegatemp*
     +                 (s12*w12 - s23*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              t33 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre/omegatemp*
     +                 (s23*w23 + s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              t12 =-2.*vist3d(j,k,i)*s12
     +           -2.*alpa1*vist3d(j,k,i)*factre/omegatemp*
     +                 (s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*vist3d(j,k,i)*s13
     +           -2.*alpa1*vist3d(j,k,i)*factre/omegatemp*
     +                 (s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*vist3d(j,k,i)*s23
     +           -2.*alpa1*vist3d(j,k,i)*factre/omegatemp*
     +                 (s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
            else if(ivmx .eq. 11) then
               alpa1=(2.-c4)/2.*gg
               alpa2=(2.-c3)*gg
              eta=(2.-c3)**2*(gg*gg/4.)*xis*
     +            turre(j,k,i,2)**2/(turre(j,k,i,1)*re)**2
              squig=(2.-c4)**2*(gg*gg/4.)*wis*
     +            turre(j,k,i,2)**2/(turre(j,k,i,1)*re)**2
              eta=ccmincr(eta,10.)
              squig=ccmincr(squig,10.)
              factre=(3.*(1.+eta)+.2e-8*(eta*eta*eta +
     +          squig*squig*squig))/
     +               (3.*(1.+eta)+   .2*(eta*eta*eta +
     +          squig*squig*squig))
c            Find nonlinear tauij values:
              t11 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              t22 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              t33 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              t12 =-2.*vist3d(j,k,i)*s12
     +           -2.*alpa1*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*vist3d(j,k,i)*s13
     +           -2.*alpa1*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*vist3d(j,k,i)*s23
     +           -2.*alpa1*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
            else if(ivmx .eq. 9 .or. ivmx .eq. 13) then
               al10 = 0.25*c1 -1.
               al1 = 3.8
c          2-line improvement for wake
               if (ieasm_type .eq. 0) then
                 al10=al10+1.8864
                 al1=al1-2.
               end if
               al3 = 0.5*c3 -1.
               al4 = 0.5*c4 -1.
              eta1=xis*turre(j,k,i,2)**2/(turre(j,k,i,1)*re)**2
              alpa1 = -al4/(al10-eta1*al1*cmuv(j,k,i))
              alpa2 = -2.*al3/(al10-eta1*cmuv(j,k,i)*al1)
               if (ivmx .eq. 9) then
                 alpa1=0.
                 alpa2=0.
               end if
c            Find nonlinear tauij values:
              t11 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              t22 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              t33 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              t12 =-2.*vist3d(j,k,i)*s12
     +           -2.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*vist3d(j,k,i)*s13
     +           -2.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*vist3d(j,k,i)*s23
     +           -2.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
            else if(ivmx .eq. 8 .or. ivmx .eq. 14) then
               al10 = 0.25*c1 -1.
               al1 = 3.8
c          2-line improvement for wake
               if (ieasm_type .eq. 0 .or. ieasm_type .eq. 3 .or.
     +             ieasm_type .eq. 4) then
                 al10=al10+1.8864
                 al1=al1-2.
               end if
               al3 = 0.5*c3 -1.
               al4 = 0.5*c4 -1.
c              Girimaji JFM 2000 fix to c4
               if (ieasm_type .eq. 4) then
                 eta1_girimaji=xis/(xis+wis)
                 if (real(eta1_girimaji) .lt. 0.5) then
                   c4new=2.0-((2.0-c4)*
     +                   (eta1_girimaji/(1.-eta1_girimaji))**0.75)
                 else
                   c4new=c4
                 end if
                 al4 = 0.5*c4new -1.
               end if
              eta1=xis/(turre(j,k,i,1)*re)**2
              alpa1 = -al4/(al10-eta1*al1*cmuv(j,k,i))
              alpa2 = -2.*al3/(al10-eta1*cmuv(j,k,i)*al1)
               if (ivmx .eq. 8) then
                 alpa1=0.
                 alpa2=0.
               end if
c            Find nonlinear tauij values:
              t11 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)/
     +                 turre(j,k,i,1)*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              t22 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)/
     +                 turre(j,k,i,1)*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              t33 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)/
     +                 turre(j,k,i,1)*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              t12 =-2.*vist3d(j,k,i)*s12
     +           -2.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*vist3d(j,k,i)*s13
     +           -2.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*vist3d(j,k,i)*s23
     +           -2.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
            else if(i_nonlin .ne. 0) then
             denom=ux(j,k,i,1)**2+ux(j,k,i,2)**2+ux(j,k,i,3)**2+
     +             ux(j,k,i,4)**2+ux(j,k,i,5)**2+ux(j,k,i,6)**2+
     +             ux(j,k,i,7)**2+ux(j,k,i,8)**2+ux(j,k,i,9)**2
             denom=sqrt(denom)
             denom=ccmax(denom,snonlin_lim)
             s11t=s11-(s11+s22+s33)/3.
             s22t=s22-(s11+s22+s33)/3.
             s33t=s33-(s11+s22+s33)/3.
c            Find nonlinear tauij values (include estimate for k!):
              t11 = 2.*q(j,k,i,1)*zk*re/3.
     +             -2.*vist3d(j,k,i)*s11t
     +             -8.*c_nonlin*vist3d(j,k,i)/denom*(-s12*w12 - s13*w13)
              t22 = 2.*q(j,k,i,1)*zk*re/3.
     +             -2.*vist3d(j,k,i)*s22t
     +             -8.*c_nonlin*vist3d(j,k,i)/denom*( s12*w12 - s23*w23)
              t33 = 2.*q(j,k,i,1)*zk*re/3.
     +             -2.*vist3d(j,k,i)*s33t
     +             -8.*c_nonlin*vist3d(j,k,i)/denom*( s23*w23 + s13*w13)
              t12 =-2.*vist3d(j,k,i)*s12
     +           -4.*c_nonlin*vist3d(j,k,i)/denom*
     +                 ( s11t*w12 - s22t*w12 - s13*w23 - s23*w13)
              t13 =-2.*vist3d(j,k,i)*s13
     +           -4.*c_nonlin*vist3d(j,k,i)/denom*
     +                 ( s11t*w13 + s12*w23 - s23*w12 - s33t*w13)
              t23 =-2.*vist3d(j,k,i)*s23
     +           -4.*c_nonlin*vist3d(j,k,i)/denom*
     +                 ( s12*w13 + s22t*w23 + s13*w12 - s33t*w23)
            end if
c            Calculate 2-eqn-specific production term (approx form in 
c            certain circumstances). Default to vorticity-type:
            if (ivmx.eq.6 .or. ivmx.eq.7 .or. ivmx.eq.10 .or.
     +          ivmx.eq.16) then
              if (ikoprod .eq. 1) then
                pk=2.*vist3d(j,k,i)*xis/(q(j,k,i,1)*re)
              else
                pk=2.*vist3d(j,k,i)*wis/(q(j,k,i,1)*re)
              end if
            else if (ivmx.eq.8 .or. ivmx.eq.9 .or. ivmx.eq.13 .or.
     +          ivmx.eq.14) then
              if (iturbprod .eq. 1) then
                pk = -(t11*ux(j,k,i,1) + t22*ux(j,k,i,5)
     +            +t33*ux(j,k,i,9) + t12*(ux(j,k,i,2)+ux(j,k,i,4))
     +                             + t13*(ux(j,k,i,3)+ux(j,k,i,7))
     +                             + t23*(ux(j,k,i,6)+ux(j,k,i,8)))/
     +            (q(j,k,i,1)*re)
              else
                pk=2.*vist3d(j,k,i)*xis/(q(j,k,i,1)*re)
              end if
            else if (ivmx.eq.11 .or. ivmx.eq.12) then
              pk = -(t11*ux(j,k,i,1) + t22*ux(j,k,i,5)
     +            +t33*ux(j,k,i,9) + t12*(ux(j,k,i,2)+ux(j,k,i,4))
     +                             + t13*(ux(j,k,i,3)+ux(j,k,i,7))
     +                             + t23*(ux(j,k,i,6)+ux(j,k,i,8)))/
     +            (q(j,k,i,1)*re)
            else if (ivmx .eq. 15) then
              if (ikoprod .eq. 2) then
                pk = -(t11*ux(j,k,i,1) + t22*ux(j,k,i,5)
     +            +t33*ux(j,k,i,9) + t12*(ux(j,k,i,2)+ux(j,k,i,4))
     +                             + t13*(ux(j,k,i,3)+ux(j,k,i,7))
     +                             + t23*(ux(j,k,i,6)+ux(j,k,i,8)))/
     +            (q(j,k,i,1)*re)
              else if (ikoprod .eq. 1) then
                pk=2.*vist3d(j,k,i)*xis/(q(j,k,i,1)*re)
              else
                pk=2.*vist3d(j,k,i)*wis/(q(j,k,i,1)*re)
              end if
            else
              pk=2.*vist3d(j,k,i)*wis/(q(j,k,i,1)*re)
            end if
            pk=pk*cutoff
c
c   optional printout to fort.50 (specifically for k-line; k=1 must be body)
c   writes "plus" turbulent values; to be used only for single block runs
c   (if multiple blocks, you will only get the last one (others get overwritten))
      if (ifort50write .eq. 1) then
      if (j .eq. j_ifort50write .and. i .eq. i_ifort50write) then
        yget=ccabs(smin(j,k,i))
        qset=sqrt(q(j,1,i,2)**2+q(j,1,i,3)**2+q(j,1,i,4)**2)
        tt=gamma*q(j,1,i,5)/q(j,1,i,1)
        c2b=198.6/tinf
        c2bp=c2b+1.0
        fnuw=c2bp*tt*sqrt(tt)/(c2b+tt)
        utau=sqrt(fnuw*qset/(ccabs(smin(j,1,i))*
     +     q(j,1,i,1)*re))
        ypl=re*q(j,1,i,1)*utau*ccabs(smin(j,k,i))/fnuw
        upl=sqrt(q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)**2)
     +        /utau
        zkplus=zk/(utau**2)
        uuplus=abs(t11)/(re*utau*utau)
        vvplus=abs(t22)/(re*utau*utau)
        wwplus=abs(t33)/(re*utau*utau)
        uwplus=t13/(re*utau*utau)
        vwplus=t23/(re*utau*utau)
        uvplus=t12/(re*utau*utau)
        uloc=sqrt(q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)**2)/xmach
        ttx=gamma*q(j,k,i,5)/q(j,k,i,1)
        fnu=c2bp*ttx*sqrt(ttx)/(c2b+ttx)
        if (ivmx .eq. 9 .or. ivmx .eq. 10 .or. ivmx .eq. 11 .or.
     +      ivmx .eq. 13) then
        eplus=turre(j,k,i,1)*fnu/(q(j,k,i,1)*utau*utau*utau*utau)
        else
        eplus=0.
        end if
        if(ivmx .eq. 8 .or. ivmx .eq. 9 .or. ivmx .ge. 13) then
          cmuxx=-cmuv(j,k,i)
        else
          cmuxx=.09
        end if
        write(50,'(14e15.5)') yget,uuplus,wwplus,vvplus,uwplus,
     +   uvplus,vwplus,uloc,zkplus,ypl,upl,vist3d(j,k,i),eplus,
     +   cmuxx
      end if
      end if
c
c      u'/u_inf:
      uprime=sqrt(abs(t11)/(q(j,k,i,1)*re))/xmach
c      v'/u_inf:
      if (ialph .eq. 0) then
        vprime=sqrt(abs(t22)/(q(j,k,i,1)*re))/xmach
      else
        vprime=sqrt(abs(t33)/(q(j,k,i,1)*re))/xmach
      end if
c      w'/u_inf:
      if (ialph .eq. 0) then
        wprime=sqrt(abs(t33)/(q(j,k,i,1)*re))/xmach
      else
        wprime=sqrt(abs(t22)/(q(j,k,i,1)*re))/xmach
      end if
c      u'u'/(u_inf)**2:
      uuprime=t11/(q(j,k,i,1)*re)/xmach**2
c      v'v'/(u_inf)**2:
      if (ialph .eq. 0) then
        vvprime=t22/(q(j,k,i,1)*re)/xmach**2
      else
        vvprime=t33/(q(j,k,i,1)*re)/xmach**2
      end if
c      w'w'/(u_inf)**2:
      if (ialph .eq. 0) then
        wwprime=t33/(q(j,k,i,1)*re)/xmach**2
      else
        wwprime=t22/(q(j,k,i,1)*re)/xmach**2
      end if
c      u'v'/(u_inf)**2:
      if (ialph .eq. 0) then
        uvprime=t12/(q(j,k,i,1)*re)/xmach**2
      else
        uvprime=t13/(q(j,k,i,1)*re)/xmach**2
      end if
c      u'w'/(u_inf)**2:
      if (ialph .eq. 0) then
        uwprime=t13/(q(j,k,i,1)*re)/xmach**2
      else
        uwprime=-t12/(q(j,k,i,1)*re)/xmach**2
      end if
c      v'w'/(u_inf)**2:
      if (ialph .eq. 0) then
        vwprime=t23/(q(j,k,i,1)*re)/xmach**2
      else
        vwprime=-t23/(q(j,k,i,1)*re)/xmach**2
      end if
c      Laminar Viscosity:
      c2b=cbar/tinf
      c2bp=c2b+1.0
      tt=gamma*q(j,k,i,5)/q(j,k,i,1)
      fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
c      S*k/epsilon:
      if (ivmx.eq.6 .or. ivmx.eq.7) then
        ske=sqrt(2.*xis)/(turre(j,k,i,1)*re*0.09)
      else if (ivmx.eq.8 .or. ivmx.eq.12 .or. ivmx .eq. 14) then
        ske=sqrt(2.*xis)/(turre(j,k,i,1)*re)
      else if (ivmx.eq.9 .or. ivmx.eq.10 .or. ivmx.eq.11 .or. 
     .         ivmx.eq.13) then
        ske=sqrt(2.*xis)*turre(j,k,i,2)/(turre(j,k,i,1)*re)
      else if (ivmx.eq.15) then
        ske=sqrt(2.*xis)*turre(j,k,i,2)*q(j,k,i,1)/
     +    (fnu*turre(j,k,i,1)*re)
      else if (ivmx.eq.16) then
        epsill=0.09**0.75*turre(j,k,i,2)**2.5/turre(j,k,i,1)
        ske=sqrt(2.*xis)*turre(j,k,i,2)/(epsill*re)
      end if
c      P(k)*Lr/(u_inf)**3  (where Lr is reference length corresponding to
c      unit 1 of grid)
      prod=pk/(xmach*xmach*xmach)
c      Fd blending function for DDES
      dist = ccabs(smin(j,k,i))
      velterm=ux(j,k,i,1)**2 + ux(j,k,i,2)**2 + ux(j,k,i,3)**2 +
     +        ux(j,k,i,4)**2 + ux(j,k,i,5)**2 + ux(j,k,i,6)**2 +
     +        ux(j,k,i,7)**2 + ux(j,k,i,8)**2 + ux(j,k,i,9)**2
      rd = (vist3d(j,k,i)+fnu)/(q(j,k,i,1)*sqrt(velterm)*0.41*0.41*
     +     dist*dist*re)
      fd = 1.0 - cctanh((8.0*rd)*(8.0*rd)*(8.0*rd)) 
c      P(k)/epsilon:
      if (ivmx.eq.6 .or. ivmx.eq.7) then
        povere=pk/(turre(j,k,i,1)*turre(j,k,i,2)*re*0.09)
      else if (ivmx.eq.8 .or. ivmx.eq.12 .or. ivmx .eq. 14) then
        povere=pk/(turre(j,k,i,1)*turre(j,k,i,2)*re)
      else if (ivmx.eq.9 .or. ivmx.eq.10 .or. ivmx.eq.11 .or.
     .         ivmx.eq.13) then
        povere=pk/(turre(j,k,i,1)*re)
      else if (ivmx.eq.15) then
        povere=pk*q(j,k,i,1)/(fnu*turre(j,k,i,1)*re)
      else if (ivmx.eq.16) then
        epsill=0.09**0.75*turre(j,k,i,2)**2.5/turre(j,k,i,1)
        povere=pk/(epsill*re)
      end if
c      Turb Reynolds number (4/9*k^2/(nu*epsilon), from NASA-CR 198221, 1995:
      if (ivmx.eq.6 .or. ivmx.eq.7) then
        turbre=4.*q(j,k,i,1)*turre(j,k,i,2)/(9.*fnu*turre(j,k,i,1)*0.09)
      else if (ivmx.eq.8 .or. ivmx.eq.12 .or. ivmx .eq. 14) then
        turbre=4.*q(j,k,i,1)*turre(j,k,i,2)/(9.*fnu*turre(j,k,i,1))
      else if (ivmx.eq.9 .or. ivmx.eq.10 .or. ivmx.eq.11 .or.
     .         ivmx.eq.13) then
        turbre=4.*q(j,k,i,1)*turre(j,k,i,2)**2/(9.*fnu*turre(j,k,i,1))
      else if (ivmx.eq.15) then
        turbre=4.*q(j,k,i,1)*q(j,k,i,1)*turre(j,k,i,2)**2/
     +    (9.*fnu*fnu*turre(j,k,i,1))
      else if (ivmx.eq.16) then
        epsill=0.09**0.75*turre(j,k,i,2)**2.5/turre(j,k,i,1)
        turbre=4.*q(j,k,i,1)*turre(j,k,i,2)**2/(9.*fnu*epsill)
      end if
c      Mean deformation (-1: pure solid body rotation, 0: pure shear,
c      +1: pure strain):
      def=(xis-wis)/(xis+wis+1.e-20)
c      Turbulence intensity:
      if (ivmx .ge. 6) then
        tint=sqrt(2.*turre(j,k,i,2)/3.)/xmach
      end if
c      anisotropy tensor component b11:
      b11=t11/(2.*q(j,k,i,1)*zk*re)-.33333333
c      anisotropy tensor component b22:
      if (ialph .eq. 0) then
        b22=t22/(2.*q(j,k,i,1)*zk*re)-.33333333
      else
        b22=t33/(2.*q(j,k,i,1)*zk*re)-.33333333
      end if
c      anisotropy tensor component b33:
      if (ialph .eq. 0) then
        b33=t33/(2.*q(j,k,i,1)*zk*re)-.33333333
      else
        b33=t22/(2.*q(j,k,i,1)*zk*re)-.33333333
      end if
c      anisotropy tensor component b12:
      if (ialph .eq. 0) then
        b12=t12/(2.*q(j,k,i,1)*zk*re)
      else
        b12=t13/(2.*q(j,k,i,1)*zk*re)
      end if
c      anisotropy tensor component b13:
      if (ialph .eq. 0) then
        b13=t13/(2.*q(j,k,i,1)*zk*re)
      else
        b13=-t12/(2.*q(j,k,i,1)*zk*re)
      end if
c      anisotropy tensor component b23:
      if (ialph .eq. 0) then
        b23=t23/(2.*q(j,k,i,1)*zk*re)
      else
        b23=-t23/(2.*q(j,k,i,1)*zk*re)
      end if
c      scalar invariant II:
      if (i2d .eq. 1) then
        if (ialph .eq. 0) then
        xii=-0.5*(b11**2+b22**2+b33**2+2.*b13**2)
        else
        xii=-0.5*(b11**2+b22**2+b33**2+2.*b12**2)
        end if
      else
        xii=-0.5*(b11**2+b22**2+b33**2+2.*b12**2+2.*b13**2+2.*b23**2)
      end if
c      scalar invariant III:
      if (i2d .eq. 1) then
        if (ialph .eq. 0) then
        xiii=0.333333*(b11**3+b22**3+b33**3+3.*b11*b13**2+
     +                                      3.*b33*b13**2) 
        else
        xiii=0.333333*(b11**3+b22**3+b33**3+3.*b11*b12**2+
     +                                      3.*b22*b12**2) 
        end if
      else
        xiii=0.333333*(b11**3+b22**3+b33**3+
     +       3.*b11*b12**2+3.*b11*b13**2+
     +       3.*b22*b12**2+3.*b22*b23**2+
     +       3.*b33*b13**2+3.*b33*b23**2+6.*b12*b13*b23) 
      end if
c      scalar invariant F ("flatness parameter"):
      ff=1.+9.*xii+27.*xiii
c      Q-criterion 2nd invariant of velocity gradient tensor
c      (not sure if this variable makes sense in 2-d)
      if (i2d .eq. 1) then
        if (ialph .eq. 0) then
          qinvar=-0.5*(ux(j,k,i,1)**2+ux(j,k,i,9)**2+
     +                 2.*ux(j,k,i,3)*ux(j,k,i,7))
        else
          qinvar=-0.5*(ux(j,k,i,1)**2+ux(j,k,i,5)**2+
     +                 2.*ux(j,k,i,2)*ux(j,k,i,4))
        end if
      else
        qinvar=-0.5*(ux(j,k,i,1)**2+ux(j,k,i,5)**2+ux(j,k,i,9)**2+
     +               2.*ux(j,k,i,2)*ux(j,k,i,4)+ 
     +               2.*ux(j,k,i,3)*ux(j,k,i,7)+
     +               2.*ux(j,k,i,6)*ux(j,k,i,8))
      end if
c      velocities by Uref
      uuu=q(j,k,i,2)/xmach
      if (ialph .eq. 0) then
        vvv=q(j,k,i,3)/xmach
        www=q(j,k,i,4)/xmach
      else
        vvv=q(j,k,i,4)/xmach
        www=-q(j,k,i,3)/xmach
      end if
c      dissipation*Lr/(u_inf)**3  (where Lr is reference length corresponding to
c      unit 1 of grid)
      if (ivmx.eq.6 .or. ivmx.eq.7) then
        epsilon=turre(j,k,i,1)*turre(j,k,i,2)*re*0.09/
     +          (xmach*xmach*xmach)
      else if (ivmx.eq.8 .or. ivmx.eq.12 .or. ivmx .eq. 14) then
        epsilon=turre(j,k,i,1)*turre(j,k,i,2)*re/(xmach*xmach*xmach)
      else if (ivmx.eq.9 .or. ivmx.eq.10 .or. ivmx.eq.11 .or.
     .         ivmx.eq.13) then
        epsilon=turre(j,k,i,1)*re/(xmach*xmach*xmach)
      else if (ivmx.eq.15) then
        epsilon=fnu*turre(j,k,i,1)*re/q(j,k,i,1)/(xmach*xmach*xmach)
      else if (ivmx.eq.16) then
        epsill=0.09**0.75*turre(j,k,i,2)**2.5/turre(j,k,i,1)
        epsilon=epsill*re/(xmach*xmach*xmach)
      end if
c      sigma is Sk/e, using Wallin-Johansson defn of S:
      if (ivmx.eq.6 .or. ivmx.eq.7) then
        ske=sqrt(xis/2.)/(turre(j,k,i,1)*re*0.09)
      else if (ivmx.eq.8 .or. ivmx.eq.12 .or. ivmx .eq. 14) then
        ske=sqrt(xis/2.)/(turre(j,k,i,1)*re)
      else if (ivmx.eq.9 .or. ivmx.eq.10 .or. ivmx.eq.11 .or.
     .         ivmx.eq.13) then
        ske=sqrt(xis/2.)*turre(j,k,i,2)/(turre(j,k,i,1)*re)
      else if (ivmx.eq.15) then
        ske=sqrt(xis/2.)*turre(j,k,i,2)*q(j,k,i,1)/
     +    (fnu*turre(j,k,i,1)*re)
      else if (ivmx.eq.16) then
        epsill=0.09**0.75*turre(j,k,i,2)**2.5/turre(j,k,i,1)
        ske=sqrt(xis/2.)*turre(j,k,i,2)/(epsill*re)
      end if
c      pressure-strain (only applies for EASMs):
      if (ivmx.eq.8 .or. ivmx.eq.9 .or. ivmx.eq.13 .or.
     +    ivmx.eq.14) then
      s1kbk1=s11*b11+s12*b12+s13*b13
      s2kbk2=s12*b12+s22*b22+s23*b23
      s3kbk3=s13*b13+s23*b23+s33*b33
      s1kbk2=s11*b12+s12*b22+s13*b23
      s1kbk3=s11*b13+s12*b23+s13*b33
      s2kbk3=s12*b13+s22*b23+s23*b33
      b1ksk1=b11*s11+b12*s12+b13*s13
      b2ksk2=b12*s12+b22*s22+b23*s23
      b3ksk3=b13*s13+b23*s23+b33*s33
      b1ksk2=b11*s12+b12*s22+b13*s23
      b1ksk3=b11*s13+b12*s23+b13*s33
      b2ksk3=b12*s13+b22*s23+b23*s33
      sklblk=s11*b11+s22*b22+s33*b33+
     +    2.*s12*b12+2.*s13*b13+2.*s23*b23
      w11=0.
      w22=0.
      w33=0.
      w1kbk1=w11*b11+w12*b12+w13*b13
      w2kbk2=-w12*b12+w22*b22+w23*b23
      w3kbk3=-w13*b13-w23*b23+w33*b33
      w1kbk2=w11*b12+w12*b22+w13*b23
      w1kbk3=w11*b13+w12*b23+w13*b33
      w2kbk3=-w12*b13+w22*b23+w23*b33
      b1kwk1=b11*w11-b12*w12-b13*w13
      b2kwk2=b12*w12+b22*w22-b23*w23
      b3kwk3=b13*w13+b23*w23+b33*w33
      b1kwk2=b11*w12+b12*w22-b13*w23
      b1kwk3=b11*w13+b12*w23+b13*w33
      b2kwk3=b12*w13+b22*w23+b23*w33
      arg=4.*sigma-10.
      if (ivmx .eq. 14 .and. (ieasm_type .eq. 3 .or.
     +    ieasm_type .eq.4)) then
        f=0.5*(1.+cctanh(arg))
      else
        f=1.0
      end if
      c2star=1.2+(c2-1.2)*f+0.84*povere*(f-1.)
c      pstrain*Lr/(u_inf)**3  (where Lr is reference length corresponding to
c      unit 1 of grid)
      r11slow=-(3.4*epsilon+1.8*prod)*b11
      r22slow=-(3.4*epsilon+1.8*prod)*b22
      r33slow=-(3.4*epsilon+1.8*prod)*b33
      r12slow=-(3.4*epsilon+1.8*prod)*b12
      r13slow=-(3.4*epsilon+1.8*prod)*b13
      r23slow=-(3.4*epsilon+1.8*prod)*b23
      r11rapid=(c2star*zk*s11+zk*1.25*(s1kbk1+b1ksk1-0.66667*sklblk)+
     +  zk*0.4*(w1kBk1-b1kWk1))/(xmach*xmach*xmach)
      r22rapid=(c2star*zk*s22+zk*1.25*(s2kbk2+b2ksk2-0.66667*sklblk)+
     +  zk*0.4*(w2kBk2-b2kWk2))/(xmach*xmach*xmach)
      r33rapid=(c2star*zk*s33+zk*1.25*(s3kbk3+b3ksk3-0.66667*sklblk)+
     +  zk*0.4*(w3kBk3-b3kWk3))/(xmach*xmach*xmach)
      r12rapid=(c2star*zk*s12+zk*1.25*(s1kbk2+b1ksk2)+
     +  zk*0.4*(w1kBk2-b1kWk2))/(xmach*xmach*xmach)
      r13rapid=(c2star*zk*s13+zk*1.25*(s1kbk3+b1ksk3)+
     +  zk*0.4*(w1kBk3-b1kWk3))/(xmach*xmach*xmach)
      r23rapid=(c2star*zk*s23+zk*1.25*(s2kbk3+b2ksk3)+
     +  zk*0.4*(w2kBk3-b2kWk3))/(xmach*xmach*xmach)
      r11=r11slow+r11rapid
      r22=r22slow+r22rapid
      r33=r33slow+r33rapid
      r12=r12slow+r12rapid
      r13=r13slow+r13rapid
      r23=r23slow+r23rapid
      end if
c
c     load appropriate data into single precision array
c
#     ifdef CMPLX
c
c     get derivative setp size on the local processor
c
      call gradinfo(delh,-1)
#     endif
c
c     Modify Here:
c
c     By default:
c     - if 3-D, prints out all 5 quantities
c     - if 2-D, only prints out the first 4
c     If keyword ifunct is employed:
c     - writes out "ifunct" number of variables to function file
c     Note that the u, v, and w prime quantities are aligned with the x, y, and
c     z axes OF THE GRID, respectively.  You must transform these quantities 
c     appropriately if they are desired in a different coordinate system.
c
#    ifdef CMPLX
      if (ip3dgrad .eq. 0) then
         xw(jw,kw,iw,1) = real(prod)
         xw(jw,kw,iw,2) = real(uwprime)
         xw(jw,kw,iw,3) = real(uuprime)
         xw(jw,kw,iw,4) = real(wwprime)
         xw(jw,kw,iw,5) = real(ske)
c        can add more here, but only if ifunct keyword is big enough!
         if (ifunct .ge. 6) xw(jw,kw,iw,6) = real(turre(j,k,i,1))
         if (ifunct .ge. 7) xw(jw,kw,iw,7) = real(turre(j,k,i,2))
         if (ifunct .ge. 8) xw(jw,kw,iw,8) = real(vist3d(j,k,i))
         if (ifunct .ge. 9) xw(jw,kw,iw,9) = real(xii)
         if (ifunct .ge.10) xw(jw,kw,iw,10) = real(xiii)
         if (ifunct .ge.11) xw(jw,kw,iw,11) = real(povere)
         if (ifunct .ge.12) xw(jw,kw,iw,12) = real(cmuv(j,k,i))
         if (ifunct .ge.13) xw(jw,kw,iw,13) = real(b11)
         if (ifunct .ge.14) xw(jw,kw,iw,14) = real(b22)
         if (ifunct .ge.15) xw(jw,kw,iw,15) = real(b33)
         if (ifunct .ge.16) xw(jw,kw,iw,16) = real(b13)
         if (ifunct .ge.17) xw(jw,kw,iw,17) = real(uuu)
         if (ifunct .ge.18) xw(jw,kw,iw,18) = real(www)
         if (ifunct .ge.19) xw(jw,kw,iw,19) = real(r11)
         if (ifunct .ge.20) xw(jw,kw,iw,20) = real(r33)
         if (ifunct .ge.21) xw(jw,kw,iw,21) = real(r13)
         if (ifunct .ge.22) xw(jw,kw,iw,22) = real(vvprime)
         if (ifunct .ge.23) xw(jw,kw,iw,23) = real(uvprime)
         if (ifunct .ge.24) xw(jw,kw,iw,24) = real(vwprime)
         if (ifunct .ge.25) xw(jw,kw,iw,25) = real(turre(j,k,i,3))
         if (ifunct .ge.26) xw(jw,kw,iw,26) = real(turre(j,k,i,4))
         if (ifunct .ge.27) xw(jw,kw,iw,27) = real(ux(j,k,i,1))
         if (ifunct .ge.28) xw(jw,kw,iw,28) = real(ux(j,k,i,3))
         if (ifunct .ge.29) xw(jw,kw,iw,29) = real(ux(j,k,i,7))
         if (ifunct .ge.30) xw(jw,kw,iw,30) = real(ux(j,k,i,9))
         if (ifunct .ge.31) xw(jw,kw,iw,31) = real(smin(j,k,i))
      else
         xw(jw,kw,iw,1) = imag(prod)/real(delh)
         xw(jw,kw,iw,2) = imag(uwprime)/real(delh)
         xw(jw,kw,iw,3) = imag(uuprime)/real(delh)
         xw(jw,kw,iw,4) = imag(wwprime)/real(delh)
         xw(jw,kw,iw,5) = imag(ske)/real(delh)
c        can add more here, but only if ifunct keyword is big enough!
         if (ifunct .ge. 6) xw(jw,kw,iw,6) = imag(turre(j,k,i,1))/
     +     real(delh)
         if (ifunct .ge. 7) xw(jw,kw,iw,7) = imag(turre(j,k,i,2))/
     +     real(delh)
         if (ifunct .ge. 8) xw(jw,kw,iw,8) = imag(vist3d(j,k,i))/
     +     real(delh)
         if (ifunct .ge. 9) xw(jw,kw,iw,9) = imag(xii)/real(delh)
         if (ifunct .ge.10) xw(jw,kw,iw,10) = imag(xiii)/real(delh)
         if (ifunct .ge.11) xw(jw,kw,iw,11) = imag(povere)/real(delh)
         if (ifunct .ge.12) xw(jw,kw,iw,12) = imag(cmuv(j,k,i))/
     +     real(delh)
         if (ifunct .ge.13) xw(jw,kw,iw,13) = imag(b11)/real(delh)
         if (ifunct .ge.14) xw(jw,kw,iw,14) = imag(b22)/real(delh)
         if (ifunct .ge.15) xw(jw,kw,iw,15) = imag(b33)/real(delh)
         if (ifunct .ge.16) xw(jw,kw,iw,16) = imag(b13)/real(delh)
         if (ifunct .ge.17) xw(jw,kw,iw,17) = imag(uuu)/real(delh)
         if (ifunct .ge.18) xw(jw,kw,iw,18) = imag(www)/real(delh)
         if (ifunct .ge.19) xw(jw,kw,iw,19) = imag(r11)/real(delh)
         if (ifunct .ge.20) xw(jw,kw,iw,20) = imag(r33)/real(delh)
         if (ifunct .ge.21) xw(jw,kw,iw,21) = imag(r13)/real(delh)
         if (ifunct .ge.22) xw(jw,kw,iw,22) = imag(vvprime)/real(delh)
         if (ifunct .ge.23) xw(jw,kw,iw,23) = imag(uvprime)/real(delh)
         if (ifunct .ge.24) xw(jw,kw,iw,24) = imag(vwprime)/real(delh)
         if (ifunct .ge.25) xw(jw,kw,iw,25) = imag(turre(j,k,i,3))/
     +     real(delh)
         if (ifunct .ge.26) xw(jw,kw,iw,26) = imag(turre(j,k,i,4))/
     +     real(delh)
         if (ifunct .ge.27) xw(jw,kw,iw,27) = imag(ux(j,k,i,1))/
     +     real(delh)
         if (ifunct .ge.28) xw(jw,kw,iw,28) = imag(ux(j,k,i,3))/
     +     real(delh)
         if (ifunct .ge.29) xw(jw,kw,iw,29) = imag(ux(j,k,i,7))/
     +     real(delh)
         if (ifunct .ge.30) xw(jw,kw,iw,30) = imag(ux(j,k,i,9))/
     +     real(delh)
         if (ifunct .ge.31) xw(jw,kw,iw,31) = imag(smin(j,k,i))/
     +     real(delh)
      end if
#     else
      xw(jw,kw,iw,1) = real(prod)
      xw(jw,kw,iw,2) = real(uwprime)
      xw(jw,kw,iw,3) = real(uuprime)
      xw(jw,kw,iw,4) = real(wwprime)
      xw(jw,kw,iw,5) = real(ske)
c     can add more here, but only if ifunct keyword is big enough!
      if (ifunct .ge. 6) xw(jw,kw,iw,6) = real(turre(j,k,i,1))
      if (ifunct .ge. 7) xw(jw,kw,iw,7) = real(turre(j,k,i,2))
      if (ifunct .ge. 8) xw(jw,kw,iw,8) = real(vist3d(j,k,i))
      if (ifunct .ge. 9) xw(jw,kw,iw,9) = real(xii)
      if (ifunct .ge.10) xw(jw,kw,iw,10) = real(xiii)
      if (ifunct .ge.11) xw(jw,kw,iw,11) = real(povere)
      if (ifunct .ge.12) xw(jw,kw,iw,12) = real(cmuv(j,k,i))
      if (ifunct .ge.13) xw(jw,kw,iw,13) = real(b11)
      if (ifunct .ge.14) xw(jw,kw,iw,14) = real(b22)
      if (ifunct .ge.15) xw(jw,kw,iw,15) = real(b33)
      if (ifunct .ge.16) xw(jw,kw,iw,16) = real(b13)
      if (ifunct .ge.17) xw(jw,kw,iw,17) = real(uuu)
      if (ifunct .ge.18) xw(jw,kw,iw,18) = real(www)
      if (ifunct .ge.19) xw(jw,kw,iw,19) = real(r11)
      if (ifunct .ge.20) xw(jw,kw,iw,20) = real(r33)
      if (ifunct .ge.21) xw(jw,kw,iw,21) = real(r13)
      if (ifunct .ge.22) xw(jw,kw,iw,22) = real(vvprime)
      if (ifunct .ge.23) xw(jw,kw,iw,23) = real(uvprime)
      if (ifunct .ge.24) xw(jw,kw,iw,24) = real(vwprime)
      if (ifunct .ge.25) xw(jw,kw,iw,25) = real(turre(j,k,i,3))
      if (ifunct .ge.26) xw(jw,kw,iw,26) = real(turre(j,k,i,4))
      if (ifunct .ge.27) xw(jw,kw,iw,27) = real(ux(j,k,i,1))
      if (ifunct .ge.28) xw(jw,kw,iw,28) = real(ux(j,k,i,3))
      if (ifunct .ge.29) xw(jw,kw,iw,29) = real(ux(j,k,i,7))
      if (ifunct .ge.30) xw(jw,kw,iw,30) = real(ux(j,k,i,9))
      if (ifunct .ge.31) xw(jw,kw,iw,31) = real(smin(j,k,i))
#     endif
 4000 continue
c
#if defined DIST_MPI
      end if
c
c     send/receive plot3d q data
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
c
      nnq = jw*kw*iw*ifuncdim
      nd_srce = mblk2nd(nblk)
      mytag = itag_q + nblk
c
      if (mblk2nd(nblk).eq.myid) then
         call MPI_Send (xw, nnq, MY_MPI_REAL, myhost, mytag,
     &                  mycomm, ierr)
      else if (myid .eq. myhost) then
         call MPI_Recv (xw, nnq, MY_MPI_REAL, nd_srce,
     &                  mytag, mycomm, istat, ierr)
      end if
#endif
c
      if (myid .eq. myhost) then
c
c        output turb data in plot3d q file
c
         if (ibin.eq.0) then
           if (ifunct .eq. 0) then
            if (i2d .eq. 0) then
               write(4,'(5e14.6)') xmachw,alphww,reuew,timew
               write(4,'(5e14.6)')
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
            else
               write(4,'(5e14.6)') xmachw,alphww,reuew,timew
               write(4,'(5e14.6)') 
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,4)
            end if
           else
               write(4,'(5e14.6)')
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,ifunct)
           end if
         else
           if (ifunct .eq. 0) then
            if (i2d .eq. 0) then
               write(4) xmachw,alphww,reuew,timew
               write(4)
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
            else
               write(4) xmachw,alphww,reuew,timew
               write(4)
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,4)
            end if
           else
               write(4)
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,ifunct)
           end if
         end if
      end if
c
      return
      end
