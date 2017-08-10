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
      subroutine plot3c(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,q,
     .                  qi0,qj0,qk0,x,y,z,xw,blank2,blank,xg,iflag,
     .                  vist3d,vi0,vj0,vk0,iover,nblk,nmap,
     .                  smin,ifunc,iplot,jdw,kdw,idw,
     .                  nplots,jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                  jeg,keg,ieg,ninter,iindex,intmax,nsub1,
     .                  maxxe,nblkk,nbli,limblk,isva,nblon,mxbli,
     .                  thetay,maxbl,maxgr,myid,myhost,mycomm,
     .                  mblk2nd,inpl3d,nblock,nblkpt,ipsurf)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Write the output file for the cell-centered points in
c     PLOT3D format.
c
c     outputs grid/solution in single precision for use with FAST/PLOT3D
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
      real*4    xw(jdw,kdw,idw,5),xg(jdw,kdw,idw,4)
      real*4    xmachw,alphww,reuew,timew
#else
      real xw(jdw,kdw,idw,5),xg(jdw,kdw,idw,4)
#endif
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension q(jdim,kdim,idim,5),vist3d(jdim,kdim,idim)
      dimension qi0(jdim,kdim,5,4),qj0(kdim,idim-1,5,4),
     .          qk0(jdim,idim-1,5,4)
      dimension vi0(jdim,kdim,1,4),vj0(kdim,idim-1,1,4),
     .          vk0(jdim,idim-1,1,4)
      dimension blank(jdim,kdim,idim),blank2(jdim,kdim,idim)
      dimension smin(jdim-1,kdim-1,idim-1)
      dimension nmap(maxbl),jdimg(maxbl),kdimg(maxbl),idimg(maxbl),
     .          nblcg(maxbl),jsg(maxbl),ksg(maxbl),isg(maxbl),
     .          jeg(maxbl),keg(maxbl),ieg(maxbl),
     .          iindex(intmax,6*nsub1+9),nblkpt(maxxe),nblkk(2,mxbli),
     .          limblk(2,6,mxbli),isva(2,2,mxbli),nblon(mxbli),
     .          inpl3d(nplots,11),thetay(maxbl),mblk2nd(maxbl)
c
      common /bin/ ibin,iblnk,iblnkfr,ip3dgrad
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /reyue/ reue,tinf,ivisc(3)
      common /twod/ i2d
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /igrdtyp/ ip3dgrd,ialph
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
      common /conversion/ radtodeg
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /fluid2/ pr,prt,cbar
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset = maxbl
      itag_grid = 1
      itag_q    = itag_grid + ioffset
#endif
c
c     output cell centers instead of cell-face centers
c
c     ipsurf = 0
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
c     initialize xw and xg arrays
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
      do j=1,jw
         do k=1,kw
            do i=1,iw
               do l= 1,5
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
c     average grid points to get cell centers and load into first 3
c     locations of the single precision xg array
c***********************************************************************
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
      if (ipsurf.gt. 0) then
c        surface face centers
         if (i1.eq.i2) then
            iw = 1
            if (i1.eq.1) then
               i = 1
            else
               i = idim
            end if
            kw = 0
            do 1500 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 1500 j=j1,j2,j3
            jw = jw + 1
            xg(jw,kw,1,1) = .25*(   x(j,k,i)   + x(j+1,k,i) 
     .                            + x(j,k+1,i) + x(j+1,k+1,i))
c
            xg(jw,kw,1,2) = .25*(   y(j,k,i)   + y(j+1,k,i) 
     .                            + y(j,k+1,i) + y(j+1,k+1,i))
c
            xg(jw,kw,1,3) = .25*(   z(j,k,i)   + z(j+1,k,i) 
     .                            + z(j,k+1,i) + z(j+1,k+1,i))
 1500       continue
         else if (j1.eq.j2) then
            jw = 1
            if (j1.eq.1) then
               j = 1
            else
               j = jdim
            end if
            iw = 0
            do 1510 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 1510 k=k1,k2,k3
            kw = kw + 1
            xg(1,kw,iw,1) = .25*(   x(j,k,i)   + x(j,k,i+1)
     .                            + x(j,k+1,i) + x(j,k+1,i+1))
c
            xg(1,kw,iw,2) = .25*(   y(j,k,i)   + y(j,k,i+1)
     .                            + y(j,k+1,i) + y(j,k+1,i+1))
c
            xg(1,kw,iw,3) = .25*(   z(j,k,i)   + z(j,k,i+1)
     .                            + z(j,k+1,i) + z(j,k+1,i+1))
 1510       continue
         else if (k1.eq.k2) then
            kw = 1
            if (k1.eq.1) then
               k = 1
            else
               k = kdim
            end if
            iw = 0
            do 1520 i=i1,i2,i3
            iw = iw + 1
            jw = 0
            do 1520 j=j1,j2,j3
            jw = jw + 1
            xg(jw,1,iw,1) = .25*(   x(j,k,i)   + x(j,k,i+1)
     .                            + x(j+1,k,i) + x(j+1,k,i+1))
c
            xg(jw,1,iw,2) = .25*(   y(j,k,i)   + y(j,k,i+1)
     .                            + y(j+1,k,i) + y(j+1,k,i+1))
c
            xg(jw,1,iw,3) = .25*(   z(j,k,i)   + z(j,k,i+1)
     .                            + z(j+1,k,i) + z(j+1,k,i+1))
 1520       continue
         end if
      else
c        cell centers
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
     .                           + x(j+1,k,i)  +x(j+1,k,i+1)
     .                           + x(j,k+1,i)  +x(j,k+1,i+1)
     .                           + x(j+1,k+1,i)+x(j+1,k+1,i+1))
c
         xg(jw,kw,iw,2) = .125*(   y(j,k,i)    +y(j,k,i+1)
     .                           + y(j+1,k,i)  +y(j+1,k,i+1)
     .                           + y(j,k+1,i)  +y(j,k+1,i+1)
     .                           + y(j+1,k+1,i)+y(j+1,k+1,i+1))
c
         xg(jw,kw,iw,3) = .125*(   z(j,k,i)    +z(j,k,i+1)
     .                           + z(j+1,k,i)  +z(j+1,k,i+1)
     .                           + z(j,k+1,i)  +z(j,k+1,i+1)
     .                           + z(j+1,k+1,i)+z(j+1,k+1,i+1))
 2000 continue
      end if
c
#if defined DIST_MPI
      end if
c
#endif
c***********************************************************************
c     set iblank (blank2) array
c***********************************************************************
c
c     currently don't need blanking data for printout option (iflag=2),
c     so only compute blank2 array for plot3d output option
c     also, don't compute blanking info if iblnk = 0
c
      if (iflag.eq.1 .and. iblnk.gt.0) then
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
      if (iflag.eq.1) then
c
c***********************************************************************
c      plot3d data
c***********************************************************************
c
      if (lhdr.gt.0 .and. myid.eq.myhost) then
         if (i2d .eq. 1) then
            write(11,'(''writing plot3d file for JDIM X KDIM ='',i5,
     .      '' x '',i5,'' grid'')') jdim,kdim
         else
            write(11,93)idim,jdim,kdim
   93       format(45h writing plot3d file for IDIM X JDIM X KDIM =,
     .      i5,3h x ,i5,3h x ,i5,5h grid)
         end if
         if (inpl3d(iplot,2).eq.1) then
            if (ipsurf.eq.0) then
               write(11,'(''   plot3dg file is an xyz'',
     .         '' file at cell centers'')')
               write(11,'(''   plot3dq file is a    q'',
     .         '' file at cell centers'')')
            else
               write(11,'(''   plot3dg file is an xyz'',
     .         '' file at cell face centers'')')
               write(11,'(''   plot3dq file is a    q'',
     .         '' file at cell face centers'')')
            end if
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then
               write(11,'(''   NOTE: contains dq/d(), rather than q'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (ifunc.eq.3) then
            if (ipsurf.eq.0) then
               write(11,'(''   plot3dg file is an xyz'',
     .         '' file at cell centers'')')
               write(11,'(''   plot3dq file is a function'',
     .         '' file of minimum distance at cell centers'')')
            else
               write(11,'(''   plot3dg file is an xyz'',
     .         '' file at cell face centers'')')
               write(11,'(''   plot3dq file is a function'',
     .         '' file of minimum distance at cell face centers'')')
            end if
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then
               write(11,'(''   NOTE: contains dsmin/d(), rather than'',
     .                    '' smin'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (ifunc.eq.4) then 
            if (ipsurf.eq.0) then
               write(11,'(''    plot3dg file is an xyz'',
     .         '' file at cell centers'')')
               write(11,'(''    plot3dq file is a function'',
     .         '' file of mu_t/mu_lam  at cell centers'')')
            else
               write(11,'(''   plot3dg file is an xyz'',
     .         '' file at cell face centers'')')
               write(11,'(''   plot3dq file is a function'',
     .         '' file of mu_t/mu_lam at cell face centers'')')
            end if
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then
               write(11,'(''   NOTE: contains d(dmu_t/mu_lam)/d(),'',
     .                    '' rather than mu_t/mu_lam'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (ifunc.eq.5) then
            if (ipsurf.eq.0) then
               write(11,'(''   plot3dg file is an xyz file'',
     .         '' at cell centers'')')
               write(11,'(''   plot3dq file is a function'',
     .         '' file of Cp at cell centers'')')
            else
               write(11,'(''   plot3dg file is an xyz file'',
     .         '' at cell face centers'')')
               write(11,'(''   plot3dq file is a function'',
     .         '' file of Cp at cell face centers'')')
            end if
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then
               write(11,'(''   NOTE: contains dCp/d(), rather than'',
     .                    '' Cp'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (ifunc.eq.6) then
            if (ipsurf.eq.0) then
               write(11,'(''   plot3dg file is an xyz file'',
     .         '' at cell centers'')')
               write(11,'(''   plot3dq file is a function'',
     .         '' file of P/Pinf at cell centers'')')
            else
               write(11,'(''   plot3dg file is an xyz file'',
     .         '' at cell face centers'')')
               write(11,'(''   plot3dq file is a function'',
     .         '' file of P/Pinf at cell face centers'')')
            end if
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then
               write(11,'(''   NOTE: contains d(P/Pinf)/d(), rather'',
     .                    '' than P/Pinf'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (ifunc.eq.7) then
            if (ipsurf.eq.0) then
               write(11,'(''   plot3dg file is an xyz file'',
     .         '' at cell centers'')')
               write(11,'(''   plot3dq file is a function'',
     .         '' file of turbulence index at the wall'')')
            else
               write(11,'(''   plot3dg file is an xyz file'',
     .         '' at cell face centers'')')
               write(11,'(''   plot3dq file is a function'',
     .         '' file of turbulence index at the wall'')')
            end if
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then
               write(11,'(''   NOTE: ip3dgrad .ne. 0 has no'',
     .                    '' effect on this'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (i2d .eq. 1) then
            if (iblnk .eq. 0) then
               write(11,2041)
 2041          format(3x,35hplot3d files to be read with /mgrid,
     .         14h/2d qualifiers)
            else
               write(11,2141)
 2141          format(3x,41hplot3d files to be read with /mgrid/blank,
     .         14h/2d qualifiers)
            end if
         else
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
      if (myid.eq.myhost .and. icall1.eq.0) then
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
            if (iblnk.eq.0) then
               write(3,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(3,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .               (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         else
            if (iblnk.eq.0) then
               write(3,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(3,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .               (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         end if
      else
        if (i2d.eq.0) then
            if (iblnk.eq.0) then
               write(3)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(3)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .            (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
        else
            if (iblnk.eq.0) then
               write(3)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(3)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .             (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
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
c
c     load appropriate data into single precision array
c   
c     ifunc...flag to determine whether a plot3d q file or plot3d function
c             file is written. function files contain only one variable to
c             be compatable with the current version of FAST
c           = 0 plot3d q file is written at cell centers
c           > 2 function file is written
c             3 minimum distance to a solid surface (smin) at cell centers
c             4 turbulent viscosity (vist3d) at cell centers
c             5 pressure coefficient (cp) at cell centers
c             6 p/pinf at cell centers
c             7 turbulence index at cell centers, when output on body only
c
c            Note: if nplot3d < 0 in the input file, then ipsurf is
c            set to 1 and solid surface, cell face center data is output,
c            rather than cell center data for nplot3d > 0 (ipsurf=0)
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
c
#     ifdef CMPLX
c
c     get derivative setp size on the local processor
c
      call gradinfo(delh,-1)
#     endif
c
      if (ifunc.eq.0) then
c
c     load solution (q) in xw array
c
      if (ipsurf .eq. 0) then
         iw = 0
         do 4000 i=i1,i2,i3
         iw = iw + 1
         kw = 0
         do 4000 k=k1,k2,k3
         kw = kw + 1
         jw = 0
         do 4000 j=j1,j2,j3
         jw = jw + 1
#        ifdef CMPLX
         if (ip3dgrad .eq. 0) then
            xw(jw,kw,iw,1) = real(q(j,k,i,1))
            xw(jw,kw,iw,5) = real(q(j,k,i,5)/gm1+0.5*(q(j,k,i,2)**2
     .                     + q(j,k,i,3)**2+q(j,k,i,4)**2)*q(j,k,i,1))
            xw(jw,kw,iw,2) = real(q(j,k,i,1)*q(j,k,i,2))
            xw(jw,kw,iw,3) = real(q(j,k,i,1)*q(j,k,i,3))
            xw(jw,kw,iw,4) = real(q(j,k,i,1)*q(j,k,i,4))
         else
            xw(jw,kw,iw,1) = imag(q(j,k,i,1))/real(delh)
            xw(jw,kw,iw,5) = imag(q(j,k,i,5)/gm1+0.5*(q(j,k,i,2)**2
     .                     + q(j,k,i,3)**2+q(j,k,i,4)**2)*q(j,k,i,1))
     .                     / real(delh)
            xw(jw,kw,iw,2) = imag(q(j,k,i,1)*q(j,k,i,2))/real(delh)
            xw(jw,kw,iw,3) = imag(q(j,k,i,1)*q(j,k,i,3))/real(delh)
            xw(jw,kw,iw,4) = imag(q(j,k,i,1)*q(j,k,i,4))/real(delh)
         end if
#        else
         xw(jw,kw,iw,1) = q(j,k,i,1)
         xw(jw,kw,iw,5) = q(j,k,i,5)/gm1+0.5*(q(j,k,i,2)**2
     .                  + q(j,k,i,3)**2+q(j,k,i,4)**2)*q(j,k,i,1)
         xw(jw,kw,iw,2) = q(j,k,i,1)*q(j,k,i,2)
         xw(jw,kw,iw,3) = q(j,k,i,1)*q(j,k,i,3)
         xw(jw,kw,iw,4) = q(j,k,i,1)*q(j,k,i,4)
#        endif
 4000    continue
      else
         if (i1.eq.i2) then
            if (i1.eq.1) then
               mm = 2
            else
               mm = 4
            end if
            iw = 0
            do 4001 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 4001 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 4001 j=j1,j2,j3
            jw = jw + 1
#           ifdef CMPLX
            if (ip3dgrad .eq. 0) then
               xw(jw,kw,iw,1) = real(qi0(j,k,1,mm))
               xw(jw,kw,iw,5) = real(qi0(j,k,5,mm)/gm1
     .                        + 0.5*(qi0(j,k,2,mm)**2
     .                        + qi0(j,k,3,mm)**2
     .                        + qi0(j,k,4,mm)**2)*qi0(j,k,1,mm))
               xw(jw,kw,iw,2) = real(qi0(j,k,1,mm)*qi0(j,k,2,mm))
               xw(jw,kw,iw,3) = real(qi0(j,k,1,mm)*qi0(j,k,3,mm))
               xw(jw,kw,iw,4) = real(qi0(j,k,1,mm)*qi0(j,k,4,mm))
            else
               xw(jw,kw,iw,1) = imag(qi0(j,k,1,mm))/real(delh)
               xw(jw,kw,iw,5) = imag(qi0(j,k,5,mm)/gm1
     .                        + 0.5*(qi0(j,k,2,mm)**2
     .                        + qi0(j,k,3,mm)**2
     .                        + qi0(j,k,4,mm)**2)*qi0(j,k,1,mm))
     .                        / real(delh)
               xw(jw,kw,iw,2) = imag(qi0(j,k,1,mm)*qi0(j,k,2,mm))
     .                        / real(delh)
               xw(jw,kw,iw,3) = imag(qi0(j,k,1,mm)*qi0(j,k,3,mm))
     .                        / real(delh)
               xw(jw,kw,iw,4) = imag(qi0(j,k,1,mm)*qi0(j,k,4,mm))
     .                        / real(delh)
            end if
#           else
            xw(jw,kw,iw,1) = qi0(j,k,1,mm)
            xw(jw,kw,iw,5) = qi0(j,k,5,mm)/gm1+0.5*(qi0(j,k,2,mm)**2
     .                     + qi0(j,k,3,mm)**2
     .                     + qi0(j,k,4,mm)**2)*qi0(j,k,1,mm)
            xw(jw,kw,iw,2) = qi0(j,k,1,mm)*qi0(j,k,2,mm)
            xw(jw,kw,iw,3) = qi0(j,k,1,mm)*qi0(j,k,3,mm)
            xw(jw,kw,iw,4) = qi0(j,k,1,mm)*qi0(j,k,4,mm)
#           endif
 4001       continue
         else if (j1.eq.j2) then
            if (j1.eq.1) then
               mm = 2
            else
               mm = 4
            end if
            iw = 0
            do 4002 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 4002 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 4002 j=j1,j2,j3
            jw = jw + 1
#           ifdef CMPLX
            if (ip3dgrad .eq. 0) then
               xw(jw,kw,iw,1) = real(qj0(k,i,1,mm))
               xw(jw,kw,iw,5) = real(qj0(k,i,5,mm)/gm1
     .                        + 0.5*(qj0(k,i,2,mm)**2
     .                        + qj0(k,i,3,mm)**2
     .                        + qj0(k,i,4,mm)**2)*qj0(k,i,1,mm))
               xw(jw,kw,iw,2) = real(qj0(k,i,1,mm)*qj0(k,i,2,mm))
               xw(jw,kw,iw,3) = real(qj0(k,i,1,mm)*qj0(k,i,3,mm))
               xw(jw,kw,iw,4) = real(qj0(k,i,1,mm)*qj0(k,i,4,mm))
            else
               xw(jw,kw,iw,1) = imag(qj0(k,i,1,mm))/real(delh)
               xw(jw,kw,iw,5) = imag(qj0(k,i,5,mm)/gm1
     .                        + 0.5*(qj0(k,i,2,mm)**2
     .                        + qj0(k,i,3,mm)**2
     .                        + qj0(k,i,4,mm)**2)*qj0(k,i,1,mm))
     .                        / real(delh)
               xw(jw,kw,iw,2) = imag(qj0(k,i,1,mm)*qj0(k,i,2,mm))
     .                        / real(delh)
               xw(jw,kw,iw,3) = imag(qj0(k,i,1,mm)*qj0(k,i,3,mm))
     .                        / real(delh)
               xw(jw,kw,iw,4) = imag(qj0(k,i,1,mm)*qj0(k,i,4,mm))
     .                        / real(delh)
            end if
#           else
            xw(jw,kw,iw,1) = qj0(k,i,1,mm)
            xw(jw,kw,iw,5) = qj0(k,i,5,mm)/gm1+0.5*(qj0(k,i,2,mm)**2
     .                     + qj0(k,i,3,mm)**2
     .                     + qj0(k,i,4,mm)**2)*qj0(k,i,1,mm)
            xw(jw,kw,iw,2) = qj0(k,i,1,mm)*qj0(k,i,2,mm)
            xw(jw,kw,iw,3) = qj0(k,i,1,mm)*qj0(k,i,3,mm)
            xw(jw,kw,iw,4) = qj0(k,i,1,mm)*qj0(k,i,4,mm)
#           endif
 4002       continue
         else if (k1.eq.k2) then
            if (k1.eq.1) then
                  mm = 2
            else
                  mm = 4
               end if
            iw = 0
            do 4003 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 4003 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 4003 j=j1,j2,j3
            jw = jw + 1
#           ifdef CMPLX
            if (ip3dgrad .eq. 0) then
               xw(jw,kw,iw,1) = real(qk0(j,i,1,mm))
               xw(jw,kw,iw,5) = real(qk0(j,i,5,mm)/gm1
     .                        + 0.5*(qk0(j,i,2,mm)**2
     .                        + qk0(j,i,3,mm)**2
     .                        + qk0(j,i,4,mm)**2)*qk0(j,i,1,mm))
               xw(jw,kw,iw,2) = real(qk0(j,i,1,mm)*qk0(j,i,2,mm))
               xw(jw,kw,iw,3) = real(qk0(j,i,1,mm)*qk0(j,i,3,mm))
               xw(jw,kw,iw,4) = real(qk0(j,i,1,mm)*qk0(j,i,4,mm))
            else
               xw(jw,kw,iw,1) = imag(qk0(j,i,1,mm))/real(delh)
               xw(jw,kw,iw,5) = imag(qk0(j,i,5,mm)/gm1
     .                        + 0.5*(qk0(j,i,2,mm)**2
     .                        + qk0(j,i,3,mm)**2
     .                        + qk0(j,i,4,mm)**2)*qk0(j,i,1,mm))
     .                        / real(delh)
               xw(jw,kw,iw,2) = imag(qk0(j,i,1,mm)*qk0(j,i,2,mm))
     .                        / real(delh)
               xw(jw,kw,iw,3) = imag(qk0(j,i,1,mm)*qk0(j,i,3,mm))
     .                        / real(delh)
               xw(jw,kw,iw,4) = imag(qk0(j,i,1,mm)*qk0(j,i,4,mm))
     .                        / real(delh)
            end if
#           else
            xw(jw,kw,iw,1) = qk0(j,i,1,mm)
            xw(jw,kw,iw,5) = qk0(j,i,5,mm)/gm1+0.5*(qk0(j,i,2,mm)**2
     .                     + qk0(j,i,3,mm)**2
     .                     + qk0(j,i,4,mm)**2)*qk0(j,i,1,mm)
            xw(jw,kw,iw,2) = qk0(j,i,1,mm)*qk0(j,i,2,mm)
            xw(jw,kw,iw,3) = qk0(j,i,1,mm)*qk0(j,i,3,mm)
            xw(jw,kw,iw,4) = qk0(j,i,1,mm)*qk0(j,i,4,mm)
#           endif
 4003       continue
         end if
      end if
c
      else
c
c     load desired function in xw array
c
      if (ifunc.eq.3) then
         if (ipsurf .eq. 0) then
            iw = 0
            do 4010 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 4010 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 4010 j=j1,j2,j3
            jw = jw + 1
#           ifdef CMPLX
            if (ip3dgrad .eq. 0) then
               xw(jw,kw,iw,1) = real(ccabs(smin(j,k,i)))
            else
               xw(jw,kw,iw,1) = imag(ccabs(smin(j,k,i)))/real(delh)
            end if
#           else
            xw(jw,kw,iw,1) = real(ccabs(smin(j,k,i)))
#           endif
 4010       continue
         else
c           at cell face centers in a wall, smin is by definition zero
            iw = 0
            do 4011 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 4011 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 4011 j=j1,j2,j3
            jw = jw + 1
            xw(jw,kw,iw,1) = 0.0
 4011       continue
         end if
      end if
      if (ifunc.eq.4) then
         if (ipsurf .eq. 0) then
            iw = 0
            do 4020 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 4020 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 4020 j=j1,j2,j3
            jw = jw + 1
#           ifdef CMPLX
            if (ip3dgrad .eq. 0) then
               xw(jw,kw,iw,1) = real(vist3d(j,k,i))
            else
               xw(jw,kw,iw,1) = imag(vist3d(j,k,i))/real(delh)
            end if
#           else
            xw(jw,kw,iw,1) = vist3d(j,k,i)
#           endif
 4020       continue
         else
            if (i1.eq.i2) then
               if (i1.eq.1) then
                  mm = 2
               else
                  mm = 4
               end if
               iw = 0
               do 4021 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 4021 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 4021 j=j1,j2,j3
               jw = jw + 1
#              ifdef CMPLX
               if (ip3dgrad .eq. 0) then
                  xw(jw,kw,iw,1) = real(vi0(j,k,1,mm))
               else
                  xw(jw,kw,iw,1) = imag(vi0(j,k,1,mm))/real(delh)
               end if
#              else
               xw(jw,kw,iw,1) = vi0(j,k,1,mm)
#              endif
 4021          continue
            else if (j1.eq.j2) then
               if (j1.eq.1) then
                  mm = 2
               else
                  mm = 4
               end if
               iw = 0
               do 4022 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 4022 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 4022 j=j1,j2,j3
               jw = jw + 1
#              ifdef CMPLX
               if (ip3dgrad .eq. 0) then
                  xw(jw,kw,iw,1) = real(vj0(k,i,1,mm))
               else
                  xw(jw,kw,iw,1) = imag(vj0(k,i,1,mm))/real(delh)
               end if
#              else
               xw(jw,kw,iw,1) = vj0(k,i,1,mm)
#              endif
 4022          continue
            else if (k1.eq.k2) then
               if (k1.eq.1) then
                     mm = 2
               else
                     mm = 4
                  end if
               iw = 0
               do 4023 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 4023 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 4023 j=j1,j2,j3
               jw = jw + 1
#              ifdef CMPLX
               if (ip3dgrad .eq. 0) then
                  xw(jw,kw,iw,1) = real(vk0(j,i,1,mm))
               else
                  xw(jw,kw,iw,1) = imag(vk0(j,i,1,mm))/real(delh)
               end if
#              else
               xw(jw,kw,iw,1) = vk0(j,i,1,mm)
#              endif
 4023          continue
            end if
         end if
      end if
      if (ifunc.eq.5) then
         cpc = 2.e0/(gamma*xmach*xmach)
         if (ipsurf.eq.0) then
            iw = 0
            do 4030 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 4030 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 4030 j=j1,j2,j3
            jw = jw + 1
#           ifdef CMPLX
            if (ip3dgrad .eq. 0) then
               xw(jw,kw,iw,1) = real((q(j,k,i,5)/p0-1)*cpc)
            else
               xw(jw,kw,iw,1) = imag((q(j,k,i,5)/p0-1)*cpc)/real(delh)
            end if
#           else
            xw(jw,kw,iw,1) = (q(j,k,i,5)/p0-1)*cpc
#           endif
 4030       continue
         else
            if (i1.eq.i2) then
               if (i1.eq.1) then
                  mm = 2
               else
                  mm = 4
               end if
               iw = 0
               do 4031 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 4031 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 4031 j=j1,j2,j3
               jw = jw + 1
#              ifdef CMPLX
               if (ip3dgrad .eq. 0) then
                  xw(jw,kw,iw,1) = real((qi0(j,k,5,mm)/p0-1)*cpc)
               else
                  xw(jw,kw,iw,1) = imag((qi0(j,k,5,mm)/p0-1)*cpc)
     .                           / real(delh)
               end if
#              else
               xw(jw,kw,iw,1) = (qi0(j,k,5,mm)/p0-1)*cpc
#              endif
 4031          continue
            else if (j1.eq.j2) then
               if (j1.eq.1) then
                  mm = 2
               else
                  mm = 4
               end if
               iw = 0
               do 4032 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 4032 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 4032 j=j1,j2,j3
               jw = jw + 1
#              ifdef CMPLX
               if (ip3dgrad .eq. 0) then
                  xw(jw,kw,iw,1) = real((qj0(k,i,5,mm)/p0-1)*cpc)
               else
                  xw(jw,kw,iw,1) = imag((qj0(k,i,5,mm)/p0-1)*cpc)
     .                           / real(delh)
               end if
#              else
               xw(jw,kw,iw,1) = (qj0(k,i,5,mm)/p0-1)*cpc
#              endif
 4032          continue
            else if (k1.eq.k2) then
               if (k1.eq.1) then
                     mm = 2
               else
                     mm = 4
                  end if
               iw = 0
               do 4033 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 4033 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 4033 j=j1,j2,j3
               jw = jw + 1
#              ifdef CMPLX
               if (ip3dgrad .eq. 0) then
                  xw(jw,kw,iw,1) = real((qk0(j,i,5,mm)/p0-1)*cpc)
               else
                  xw(jw,kw,iw,1) = imag((qk0(j,i,5,mm)/p0-1)*cpc)
     .                           / real(delh)
               end if
#              else
               xw(jw,kw,iw,1) = (qk0(j,i,5,mm)/p0-1)*cpc
#              endif
 4033          continue
            end if
         end if
      end if
      if (ifunc.eq.6) then
         if (ipsurf.eq.0) then
            iw = 0
            do 4040 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 4040 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 4040 j=j1,j2,j3
            jw = jw + 1
#           ifdef CMPLX
            if (ip3dgrad .eq. 0) then
               xw(jw,kw,iw,1) = real(q(j,k,i,5)/p0)
            else
               xw(jw,kw,iw,1) = imag(q(j,k,i,5)/p0)/real(delh)
            end if
#           else
            xw(jw,kw,iw,1) = q(j,k,i,5)/p0
#           endif
 4040       continue
         else
            if (i1.eq.i2) then
               if (i1.eq.1) then
                  mm = 2
               else
                  mm = 4
               end if
               iw = 0
               do 4041 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 4041 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 4041 j=j1,j2,j3
               jw = jw + 1
#              ifdef CMPLX
               if (ip3dgrad .eq. 0) then
                  xw(jw,kw,iw,1) = real(qi0(j,k,5,mm)/p0)
               else
                  xw(jw,kw,iw,1) = imag(qi0(j,k,5,mm)/p0)/real(delh)
               end if
#              else
               xw(jw,kw,iw,1) = qi0(j,k,5,mm)/p0
#              endif
 4041          continue
            else if (j1.eq.j2) then
               if (j1.eq.1) then
                  mm = 2
               else
                  mm = 4
               end if
               iw = 0
               do 4042 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 4042 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 4042 j=j1,j2,j3
               jw = jw + 1
#              ifdef CMPLX
               if (ip3dgrad .eq. 0) then
                  xw(jw,kw,iw,1) = real(qj0(k,i,5,mm)/p0)
               else
                  xw(jw,kw,iw,1) = imag(qj0(k,i,5,mm)/p0)/real(delh)
               end if
#              else
               xw(jw,kw,iw,1) = qj0(k,i,5,mm)/p0
#              endif
 4042          continue
            else if (k1.eq.k2) then
               if (k1.eq.1) then
                  mm = 2
               else
                  mm = 4
               end if
               iw = 0
               do 4043 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 4043 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 4043 j=j1,j2,j3
               jw = jw + 1
#              ifdef CMPLX
               if (ip3dgrad .eq. 0) then
                  xw(jw,kw,iw,1) = real(qk0(j,i,5,mm)/p0)
               else
                  xw(jw,kw,iw,1) = imag(qk0(j,i,5,mm)/p0)/real(delh)
               end if
#              else
               xw(jw,kw,iw,1) = qk0(j,i,5,mm)/p0
#              endif
 4043          continue
            end if
         end if
      end if
      if (ifunc.eq.7) then
         iw = 0
         do 5768 i=i1,i2,i3
         iw = iw + 1
         kw = 0
         do 5768 k=k1,k2,k3
         kw = kw + 1
         jw = 0
         do 5768 j=j1,j2,j3
         jw = jw + 1
c        vorticity at wall is approximate
         vor=sqrt(q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)**2)/
     +       (ccabs(smin(j,k,i)))
         tt=gamma*q(j,k,i,5)/q(j,k,i,1)
         c2b=cbar/tinf
         c2bp=c2b+1.0
         xmu=c2bp*tt*sqrt(tt)/(c2b+tt)
c        iterate to get nuwiggle
         xmut=vist3d(j,k,i)
         xnuwiggle=500.
         do n=1,50
           delta=(xnuwiggle**4 - xmut*xnuwiggle**3-xmut*(7.1**3))/
     +           (4.*xnuwiggle**3-3.*xmut*xnuwiggle**2)
           xnuwiggle=xnuwiggle-delta
c          iteration may not be working
           xnuwiggle=ccmaxcr(xnuwiggle,1.e-20)
           xnuwiggle=ccmincr(xnuwiggle,5000.)
         enddo
#        ifdef CMPLX
         xw(jw,kw,iw,1) = real(xnuwiggle)/real(ccabs(smin(j,k,i)))*
     +    sqrt(real(q(j,k,i,1)))*sqrt(xmach/reue)/(0.41*
     +    sqrt(real(xmu)*real(vor)))
#        else
         xw(jw,kw,iw,1) = xnuwiggle/(ccabs(smin(j,k,i)))*
     +    sqrt(q(j,k,i,1)*xmach/reue)/(0.41*sqrt(xmu*vor))
#        endif
c        Note Spalart / Strelets formula for SST (not used):
c        (need to bring in turre and nummem if want to use it)
c        xw(jw,kw,iw,1) = 7.4*turre(j,k,i,2)**0.31/(ccabs(smin(j,k,i)))*
c    +    xmu/q(j,k,i,1)*(q(j,k,i,1)/(xmu*real(vor)))**0.81*
c    +    (xmach/reue)**0.19
c
c        limit maximum to 1 (represents turbulent)
         if (xw(jw,kw,iw,1) .gt. 1.0) then
           xw(jw,kw,iw,1)=1.0
         end if
 5768    continue
      end if
c
      end if
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
      if (ifunc.eq.0) then
         nnq = jw*kw*iw*5
      else if (ifunc.gt.0) then
         nnq = jw*kw*iw*1
      end if
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
c     output solution
c
      if (ifunc .eq. 0) then
         if (ibin.eq.0) then
            write(4,'(5e14.6)') xmachw,alphww,reuew,timew
            if (i2d.eq.0) then
               if (ialph.eq.0) then
                  write(4,'(5e14.6)')
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
               else
                  do i=1,iw
                     do j=1,jw
                        do k=1,kw
                           xw(j,k,i,3) = -xw(j,k,i,3)
                        end do
                     end do
                  end do
                  write(4,'(5e14.6)')
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .             (((xw(j,k,i,4),i=1,iw),j=1,jw),k=1,kw),
     .             (((xw(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .             (((xw(j,k,i,5),i=1,iw),j=1,jw),k=1,kw)
               end if
            else
               write(4,'(5e14.6)')
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,5)
            end if
         else
            write(4) xmachw,alphww,reuew,timew
            if (i2d.eq.0) then
               if (ialph.eq.0) then
                  write(4)
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
               else
                  do i=1,iw
                     do j=1,jw
                        do k=1,kw
                           xw(j,k,i,3) = -xw(j,k,i,3)
                        end do
                     end do
                  end do
                  write(4)
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .             (((xw(j,k,i,4),i=1,iw),j=1,jw),k=1,kw),
     .             (((xw(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .             (((xw(j,k,i,5),i=1,iw),j=1,jw),k=1,kw)
                end if
            else
                write(4)
     .          ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .          ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,5)
            end if
         end if
      else
         if (ibin.eq.0) then
            write(4,'(5e14.6)')
     .      (((xw(j,k,i,1),i=1,iw),j=1,jw),k=1,kw)
         else
            write(4)
     .      (((xw(j,k,i,1),i=1,iw),j=1,jw),k=1,kw)
         end if
      end if
c
      end if
c
      end if
c
      if (iflag.eq.2) then
c
c***********************************************************************
c     print out data
c***********************************************************************
c
c     NOTE: the turb. viscosity data are stored in the 4th position of the
c           xq array as follows
c
c           xq(j,k,i,4) = vist3d
c
c           the solution q is stored in the xw array as follows
c
c           xw(j,k,i,1) = q(1) (rho)
c           xw(j,k,i,2) = q(2) (u)
c           xw(j,k,i,3) = q(3) (v)
c           xw(j,k,i,4) = q(4) (w)
c           xw(j,k,i,5) = q(5) (p)
c
      if (lhdr.gt.0 .and. myid.eq.myhost) then
         write(11,94)idim,jdim,kdim
   94    format(47h writing printout file for IDIM X JDIM X KDIM =,
     .   i5,3h x ,i5,3h x ,i5,5h grid)
      end if
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
#endif
c     load q and turb. viscostiy into single precision arrays
c
      jdw = (j2-j1)/j3 + 1
      kdw = (k2-k1)/k3 + 1
      idw = (i2-i1)/i3 + 1
c
      if (ipsurf.eq.0) then
         iw = 0
         do 9100 i=i1,i2,i3
         iw = iw + 1
         jw = 0
         do 9100 j=j1,j2,j3
         jw = jw + 1
         kw = 0
         do 9100 k=k1,k2,k3
         kw = kw + 1
         do 9100 l=1,5
         xw(jw,kw,iw,l) = q(j,k,i,l)
 9100 continue
      else
         if (i1.eq.i2) then
            if (i1.eq.1) then
               mm = 2
            else
               mm = 4
            end if
            iw = 0
            do 9101 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 9101 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 9101 j=j1,j2,j3
            jw = jw + 1
            do 9101 l=1,5
            xw(jw,kw,iw,l) = qi0(j,k,l,mm)
 9101       continue
         else if (j1.eq.j2) then
            if (j1.eq.1) then
               mm = 2
            else
               mm = 4
            end if
            iw = 0
            do 9102 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 9102 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 9102 j=j1,j2,j3
            jw = jw + 1
            do 9102 l=1,5
            xw(jw,kw,iw,l) = qj0(k,i,l,mm)
 9102       continue
         else if (k1.eq.k2) then
            if (k1.eq.1) then
               mm = 2
            else
               mm = 4
            end if
            iw = 0
            do 9103 i=i1,i2,i3
            iw = iw + 1
            kw = 0
            do 9103 k=k1,k2,k3
            kw = kw + 1
            jw = 0
            do 9103 j=j1,j2,j3
            jw = jw + 1
            do 9103 l=1,5
            xw(jw,kw,iw,l) = qk0(j,i,l,mm)
 9103       continue
         end if
      end if
c
      if (ivisc(3).gt.1 .or. ivisc(2).gt.1 .or. ivisc(1).gt.1) then
         if (ipsurf.eq.0) then
            iw = 0
            do 9200 i=i1,i2,i3
            iw = iw + 1
            jw = 0
            do 9200 j=j1,j2,j3
            jw = jw + 1
            kw = 0
            do 9200 k=k1,k2,k3
            kw = kw + 1
            xg(jw,kw,iw,4) = vist3d(j,k,i)
 9200       continue
         else
            if (i1.eq.i2) then
               if (i1.eq.1) then
                  mm = 2
               else
                  mm = 4
               end if
               iw = 0
               do 9201 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 9201 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 9201 j=j1,j2,j3
               jw = jw + 1
               xg(jw,kw,iw,4) = vi0(j,k,1,mm)
 9201          continue
            else if (j1.eq.j2) then
               if (j1.eq.1) then
                  mm = 2
               else
                  mm = 4
               end if
               iw = 0
               do 9202 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 9202 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 9202 j=j1,j2,j3
               jw = jw + 1
               xg(jw,kw,iw,4) = vj0(k,i,1,mm)
 9202          continue
            else if (k1.eq.k2) then
               if (k1.eq.1) then
                  mm = 2
               else
                  mm = 4
               end if
               iw = 0
               do 9203 i=i1,i2,i3
               iw = iw + 1
               kw = 0
               do 9203 k=k1,k2,k3
               kw = kw + 1
               jw = 0
               do 9203 j=j1,j2,j3
               jw = jw + 1
               xg(jw,kw,iw,4) = vk0(j,i,1,mm)
 9203          continue
            end if
         end if
      end if
c
#if defined DIST_MPI
      end if
c
c     send/receive x,y,z,q and vist3d data
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
c
      nnq = jw*kw*iw*5
      nd_srce = mblk2nd(nblk)
      mytag = itag_q + nblk
c
      if (mblk2nd(nblk).eq.myid) then
         call MPI_Send (xw, nnq, MY_MPI_REAL, myhost, mytag,
     &                  mycomm, ierr)
      else if (myid.eq.myhost) then
         call MPI_Recv (xw, nnq, MY_MPI_REAL, nd_srce,
     &                  mytag, mycomm, istat, ierr)
      end if
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
      if (myid .eq. myhost) then
c
c        icp=1 prints out cp's; icp=0 prints out pitot pressures
c
         icp=1
c
         if(icp .eq. 1) then
           write(17,77)
   77      format(/4x,1hI,4x,1hJ,4x,1hK,9x,1hX,16x,1hY,16x,1hZ,
     .     14x,6hU/Uinf,11x,6hV/Vinf,11x,6hW/Winf,11x,6hP/Pinf,
     .     11x,6hT/Tinf,12x,4hMACH,10x,6h    cp,12x,9htur. vis.)
         else
           write(17,7)
    7      format(/4x,1hI,4x,1hJ,4x,1hK,9x,1hX,16x,1hY,16x,1hZ,
     .     14x,6hU/Uinf,11x,6hV/Vinf,11x,6hW/Winf,11x,6hP/Pinf,
     .     11x,6hT/Tinf,12x,4hMACH,10x,6hpitotp,12x,9htur. vis.)
         end if
c
         term3 = 1./( (1.+0.5*gm1*xmach*xmach)**(gamma/gm1) )
         iw = 0
         do 5000 i=i1,i2,i3
         iw = iw + 1
         kw = 0
         do 5000 k=k1,k2,k3
         kw = kw + 1
         jw = 0
         do 5000 j=j1,j2,j3
         jw = jw + 1
         q1    = xw(jw,kw,iw,1)
         q2    = xw(jw,kw,iw,2)/xmach
         q3    = xw(jw,kw,iw,3)/xmach
         q4    = xw(jw,kw,iw,4)/xmach
         q5    = gamma*xw(jw,kw,iw,5)
         t1    = q5/q1
         xm1   = sqrt(xmach**2*(q2**2+q3**2+q4**2)/t1)
c
c        turbulent viscosity
c
         edvis    = xg(jw,kw,iw,4)
c
c        cp or pitot pressure
c
         if(icp .eq. 1) then
            pitot = 2.*(q5-1.)/(gamma*xmach*xmach)
         else
            if (real(xm1).gt.1.0) then
               term1 = (0.5*gp1*xm1*xm1)**(gamma/gm1)
               term2 = (2.*gamma*xm1*xm1/gp1 - gm1/gp1)**(1./gm1)
               pitot = q5*term1*term3/term2
            else
               term1 = (1.0+0.5*gm1*xm1*xm1)**(gamma/gm1)
               pitot = q5*term1*term3
            end if
         end if
         if (ialph.eq.0) then
            write(17,29)i,j,k,xg(jw,kw,iw,1),
     .                  xg(jw,kw,iw,2),xg(jw,kw,iw,3),
     .                  real(q2),real(q3),real(q4),real(q5),real(t1),
     .                  real(xm1),real(pitot),real(edvis)
         else
            xg(jw,kw,iw,2) = -xg(jw,kw,iw,2)
            q3 = -q3
            write(17,29)i,j,k,xg(jw,kw,iw,1),xg(jw,kw,iw,3),
     .                  xg(jw,kw,iw,2),real(q2),real(q4),
     .                  real(q3),real(q5),real(t1),real(xm1),
     .                  real(pitot),real(edvis)
         end if
 5000    continue
c
   29    format(3i5,11(1x,e16.9))
c
      end if
c
      end if
c
      return
      end
