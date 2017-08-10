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
      subroutine plot3d_2d(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,
     .                  k1,k2,k3,q,
     .                  qi0,qj0,qk0,x,y,z,xw,blank2,blank,xg,iflag,
     .                  vist3d,iover,nblk,nmap,bcj,bck,bci,
     .                  vj0,vk0,vi0,ifunc,iplot,jdw,kdw,idw,
     .                  nplots,jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                  jeg,keg,ieg,ninter,iindex,intmax,nsub1,
     .                  maxxe,nblkk,nbli,limblk,isva,nblon,mxbli,
     .                  thetay,maxbl,maxgr,myid,myhost,mycomm,
     .                  mblk2nd,inpl3d,nblock,nblkpt,xv,
     .                  sj,sk,si,vol,nset)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Write the output file at the grid points in PLOT3D
c     format and print solution data.
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
#if defined P3D_SINGLE
      real*4    xw(jdw,kdw,idw,nset),xg(jdw,kdw,idw,4)
      real*4    xmachw,alphww,reuew,timew
#else
      real xw(jdw,kdw,idw,nset),xg(jdw,kdw,idw,4)
#endif
c
      dimension xv(jdw,kdw,idw,5)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension vist3d(jdim,kdim,idim)
      dimension vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4)
      dimension blank(jdim,kdim,idim),blank2(jdim,kdim,idim,2)
      dimension nmap(maxbl),jdimg(maxbl),kdimg(maxbl),idimg(maxbl),
     .          nblcg(maxbl),jsg(maxbl),ksg(maxbl),isg(maxbl),
     .          jeg(maxbl),keg(maxbl),ieg(maxbl),
     .          iindex(intmax,6*nsub1+9),nblkpt(maxxe),nblkk(2,mxbli),
     .          limblk(2,6,mxbli),isva(2,2,mxbli),nblon(mxbli), 
     .          inpl3d(nplots,11),thetay(maxbl),mblk2nd(maxbl)
      dimension sk(jdim,kdim,idim-1,5),sj(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
c
      common /bin/ ibin,iblnk,iblnkfr,ip3dgrad
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /reyue/ reue,tinf,ivisc(3)
      common /twod/ i2d
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /wallfun/ iwf(3)
      common /igrdtyp/ ip3dgrd,ialph
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
      common /moovcrs2d/icallcrs,icall2d
      common /conversion/ radtodeg
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
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
c     initialize xw, xv and xg arrays
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
      do j=1,jw
         do k=1,kw
            do i=1,iw
               do l= 1,5
               xw(j,k,i,l) = 0.
               xv(j,k,i,l) = 0.
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
c     load grid into first 3 locations of the single precision xg array
c***********************************************************************
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
      iw = 0
      do 2000 i=i1,i2,i3
      iw = iw + 1
      jw = 0
      do 2000 j=j1,j2,j3
      jw = jw + 1
      kw = 0
      do 2000 k=k1,k2,k3
      kw = kw + 1
      xg(jw,kw,iw,1) = x(j,k,i)
      xg(jw,kw,iw,2) = y(j,k,i)
      xg(jw,kw,iw,3) = z(j,k,i)
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
      blank2(j,k,i,1) = 1.
      blank2(j,k,i,2) = 1.
 2049 continue
c
c     zero out edges and corners if desired (izero > 0)
c
      izero = 0
c     
      if (izero.gt.0) then
         do 2050 i=1,idim
         do 2050 j=1,jdim,jdim1
         do 2050 k=1,kdim,kdim1
         blank2(j,k,i,1) = 0.
 2050    continue
c
         if (idim.gt.2) then
            do 2051 i=1,idim,idim1
            do 2151 j=1,jdim
            blank2(j,1,i,1)    = 0.
            blank2(j,kdim,i,1) = 0.
 2151       continue
            do 2251 k=1,kdim
            blank2(1,k,i,1)    = 0.
            blank2(jdim,k,i,1) = 0.
 2251       continue
 2051       continue
         else
            do 2052 j=1,jdim
            do 2052 k=1,kdim
            blank2(j,k,idim,1) = 0.
 2052       continue
         end if
c
         if (jdim.eq.2) then
            do 2054 i=1,idim
            do 2054 k=1,kdim
            blank2(jdim,k,i,1) = 0.
 2054       continue
         end if
c
         if (kdim.eq.2) then
            do 2056 i=1,idim
            do 2056 k=1,kdim
            blank2(j,kdim,i,1) = 0.
 2056       continue
         end if
      end if
c
c     solid surface iblank (blank2) values (iblank=2 for solid surface)
c
      j11 = 1
      j22 = jdim1
      if(jdim.eq.2) then
        j11 = 1
        j22 = 1
      end if
      k11 = 1
      k22 = kdim1
      if(kdim.eq.2) then
        k11 = 1
        k22 = 1
      end if
      i11 = 1
      i22 = idim1
      if(idim.eq.2) then
        i11 = 1
        i22 = 1
      end if
      i = 1
      do 1190 nnn=1,2
      do 1191 j=j11,j22
      do 1191 k=k11,k22
      kk = min(k+1,k22)
      jj = min(j+1,j22)
      blank2(j,k,i,1) = 1. + ccmax4( bci(j,k,nnn),  bci(jj,k,nnn),
     .                               bci(j,kk,nnn), bci(jj,kk,nnn) )
 1191 continue
      i = idim
 1190 continue
      j = 1
      do 2190 nnn=1,2
      do 2191 i=i11,i22
      do 2191 k=k11,k22
      kk = min(k+1,k22)
      ii = min(i+1,i22)
      blank2(j,k,i,1) = 1. + ccmax4( bcj(k,i,nnn),  bcj(kk,i,nnn),
     .                               bcj(k,ii,nnn), bcj(kk,ii,nnn) )
 2191 continue
      j = jdim
 2190 continue
      k = 1
      do 3190 nnn=1,2
      do 3191 j=j11,j22
      do 3191 i=i11,i22
      jj = min(j+1,j22)
      ii = min(i+1,i22)
      blank2(j,k,i,1) = 1. + ccmax4( bck(j,i,nnn),  bck(jj,i,nnn),
     .                               bck(j,ii,nnn), bck(jj,ii,nnn) )
 3191 continue
      k = kdim
 3190 continue
c
c     iblank (blank2) array for patch interface boundaries
c  
c     Set blank2 array on block edges from nblkpt array.  Nblkpt array
c     corresponds to cell center grid.  For ambiguous points which arise
c     when going from cell-center data to grid point data (e.g. along block
c     edges, corners, or for interior points on a patch interface for which 
c     surrounding cell-centers lie in different blocks), blank2 value 
c     remains equal to previously set values.
c
      if (abs(ninter).gt.0) then
         do 1600 inter=1,abs(ninter)
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
c
         if (lcoord.eq.1) then
            if (lend.eq.1) i = 1
            if (lend.eq.2) i = idim
            if (jdim.gt.2 .and. kdim.gt.2) then
               do 1610 j=j21+1,j22-1
               do 1610 k=k21+1,k22-1
c              connecting block numbers of the four cell-centered points 
c              surrounding the current grid point
               ll = lst + (j22-j21)*(k-1-k21) + (j-1-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21-1) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21)
               mblk3 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21-1)
               mblk4 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2 .and. mblk1.eq.mblk3 .and.
     .             mblk1.eq.mblk4) blank2(j,k,i,1) = -float(nmap(mblk1))
 1610          continue
            else if (jdim.gt.2) then
c              connecting block numbers of the two cell-centered points
c              surrounding the current grid point
               k = k21
               do 1611 j=j21+1,j22-1
               ll = lst + (j22-j21)*(k-k21) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21-1)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1611          continue
            else if (kdim.gt.2) then
               j = j21
               do 1612 k=k21+1,k22-1
               ll = lst + (j22-j21)*(k-k21-1) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1612          continue
            end if
         end if
c
         if (lcoord.eq.2) then
            if (lend.eq.1) j = 1
            if (lend.eq.2) j = jdim
            if (idim.gt.2 .and .kdim.gt.2) then
               do 1620 i=k21+1,k22-1
               do 1620 k=j21+1,j22-1
c              connecting block numbers of the four cell-centered points 
c              surrounding the current grid point
               ll = lst + (j22-j21)*(i-1-k21) + (k-1-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21-1) + (k-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21)
               mblk3 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21-1)
               mblk4 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2 .and. mblk1.eq.mblk3 .and.
     .             mblk1.eq.mblk4) blank2(j,k,i,1) = -float(nmap(mblk1))
 1620          continue
            else if (idim.gt.2) then
               k = j21
               do 1621 i=k21+1,k22-1
               ll = lst + (j22-j21)*(i-k21-1) + (k-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1621          continue
            else if (kdim.gt.2) then
               i = k21
               do 1622 k=j21+1,j22-1
               ll = lst + (j22-j21)*(i-k21) + (k-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21-1)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1622          continue
            end if
         end if
c
         if (lcoord.eq.3) then
            if (lend.eq.1) k = 1
            if (lend.eq.2) k = kdim
            if (idim.gt.2 .and. jdim.gt.2) then
               do 1630 i=k21+1,k22-1
               do 1630 j=j21+1,j22-1
c              connecting block numbers of the four cell-centered points 
c              surrounding the current grid point
               ll = lst + (j22-j21)*(i-1-k21) + (j-1-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21-1) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21)
               mblk3 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21-1)
               mblk4 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2 .and. mblk1.eq.mblk3 .and.
     .             mblk1.eq.mblk4) blank2(j,k,i,1) = -float(nmap(mblk1))
 1630          continue
            else if (idim.gt.2) then
               j = j21
               do 1631 i=k21+1,k22-1
               ll = lst + (j22-j21)*(i-k21-1) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1631          continue
            else if (jdim.gt.2) then
               i = k21
               do 1632 j=j21+1,j22-1
               ll = lst + (j22-j21)*(i-k21) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21-1)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1632          continue
            end if
         end if
c
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
c           cell center indicies ---> grid point indicies
c
            if(isva(it,1,n)+isva(it,2,n).eq.5) then
c             i = constant interface
              if(js.gt.je) js = js+1
              if(js.lt.je) je = je+1
              if(ks.gt.ke) ks = ks+1
              if(ks.lt.ke) ke = ke+1
c             2d cases
              if(jdim.eq.2) then
                js = 1
                je = 2
              end if
              if(kdim.eq.2) then
                ks = 1
                ke = 2
              end if
            end if
            if(isva(it,1,n)+isva(it,2,n).eq.4) then
c             j = constant interface
              if(is.gt.ie) is = is+1
              if(is.lt.ie) ie = ie+1
              if(ks.gt.ke) ks = ks+1
              if(ks.lt.ke) ke = ke+1
c             2d cases
              if(idim.eq.2) then
                is = 1
                ie = 2
              end if
              if(kdim.eq.2) then
                ks = 1
                ke = 2
              end if
            end if
            if(isva(it,1,n)+isva(it,2,n).eq.3) then
c             k = constant interface
              if(js.gt.je) js = js+1
              if(js.lt.je) je = je+1
              if(is.gt.ie) is = is+1
              if(is.lt.ie) ie = ie+1
c             2d cases
              if(jdim.eq.2) then
                js = 1
                je = 2
              end if
              if(idim.eq.2) then
                is = 1
                ie = 2
              end if
            end if
c
            is1 = min(is,ie)
            ie1 = max(is,ie)
            js1 = min(js,je)
            je1 = max(js,je)
            ks1 = min(ks,ke)
            ke1 = max(ks,ke)
c
            do 110 i = is1,ie1
            do 110 j = js1,je1
            do 110 k = ks1,ke1
            blank2(j,k,i,1) = -float(nmap(nblkk(ir,n)))
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
            js = jsg(nblc)
            if (js.lt.jdimg(nblcc) .and. js.gt.1) js = js + 1
            ks = ksg(nblc)
            if (ks.lt.kdimg(nblcc) .and. ks.gt.1) ks = ks + 1
            is = isg(nblc)
            if (is.lt.idimg(nblcc) .and. is.gt.1) is = is + 1
            je = jeg(nblc)
            if (je.gt.2 .and. je.lt.jdimg(nblcc)) je = je - 1
            ke = keg(nblc)
            if (ke.gt.2 .and. ke.lt.kdimg(nblcc)) ke = ke - 1
            ie = ieg(nblc)
            if (ie.gt.2 .and. ie.lt.idimg(nblcc)) ie = ie - 1
            do 7501 i=is,ie
            do 7501 j=js,je
            do 7501 k=ks,ke
            blank2(j,k,i,1) = 0.
 7501       continue
         end if
 7500    continue
      end if
c
c     iblank (blank2) array for overlapped grids
c
      if (iover.eq.1) then
c
c        interior of faces (i=1 and i=idim)
c
         do 413 j=2,jdim1
         do 413 k=2,kdim1
         blank2(j,k,1,2)    = ccmin4(blank(j,k,1),
     .                               blank(j-1,k,1),
     .                               blank(j,k-1,1),
     .                               blank(j-1,k-1,1))
         blank2(j,k,idim,2) = ccmin4(blank(j,k,idim1),  
     .                               blank(j-1,k,idim1),
     .                               blank(j,k-1,idim1),
     .                               blank(j-1,k-1,idim1))
  413    continue
c
c        edges and corners (i=1 and i=idim)
c
         do 414 m=1,2
         k  = 1
         kk = 1
         if (m.eq.2) then
            k  = kdim
            kk = kdim1
         end if
         do 414 j=2,jdim1
         blank2(j,k,1,2)    = ccmin(blank(j,kk,1),
     .                              blank(j-1,kk,1))
         blank2(j,k,idim,2) = ccmin(blank(j,kk,idim1),
     .                              blank(j-1,kk,idim1))
  414    continue
         do 415 m=1,2
         j  = 1
         jj = 1
         if (m.eq.2) then
            j  = jdim
            jj = jdim1
         end if
         do 415 k=2,kdim1
         blank2(j,k,1,2)    = ccmin(blank(jj,k,1),
     .                              blank(jj,k-1,1))
         blank2(j,k,idim,2) = ccmin(blank(jj,k,idim1),
     .                              blank(jj,k-1,idim1))
  415    continue
         blank2(jdim,kdim,1,2)    = blank2(jdim-1,kdim-1,1,2)
         blank2(1,1,1,2)          = blank2(2,2,1,2)
         blank2(jdim,kdim,idim,2) = blank2(jdim-1,kdim-1,idim,2)
         blank2(1,1,idim,2)       = blank2(2,2,idim,2)
c
c        interior of faces (j=1 and j=jdim)
c
         do 513 i=2,idim1
         do 513 k=2,kdim1
         blank2(1,k,i,2)    = ccmin4(blank(1,k,i),
     .                               blank(1,k,i-1),
     .                               blank(1,k-1,i),
     .                               blank(1,k-1,i-1))
         blank2(jdim,k,i,2) = ccmin4(blank(jdim1,k,i),
     .                               blank(jdim1,k,i-1),
     .                               blank(jdim1,k-1,i),
     .                               blank(jdim1,k-1,i-1))
  513    continue
c
c        edges and corners (j=1 and j=jdim)
c  
         do 514 m=1,2
         k  = 1
         kk = 1
         if (m.eq.2) then
            k  = kdim
            kk = kdim1
         end if
         do 514 i=2,idim1
         blank2(1,k,i,2)    = ccmin(blank(1,kk,i),
     .                              blank(1,kk,i-1))
         blank2(jdim,k,i,2) = ccmin(blank(jdim1,kk,i),
     .                              blank(jdim1,kk,i-1))
  514    continue
         do 515 m=1,2
         i  = 1
         ii = 1
         if (m.eq.2) then
            i  = idim
            ii = idim1
         end if
         do 515 k=2,kdim1
         blank2(1,k,i,2)    = ccmin(blank(1,k,ii),
     .                              blank(1,k-1,ii))
         blank2(jdim,k,i,2) = ccmin(blank(jdim1,k,ii),
     .                              blank(jdim1,k-1,ii))
  515    continue
         blank2(1,kdim,idim,2)    = blank2(1,kdim-1,idim-1,2)
         blank2(1,1,1,2)          = blank2(1,2,2,2)
         blank2(jdim,kdim,idim,2) = blank2(jdim,kdim-1,idim-1,2)
         blank2(jdim,1,1,2)       = blank2(jdim,2,2,2)
c
c        interior of faces (k=1 and k=kdim)
c
         do 613 i=2,idim1
         do 613 j=2,jdim1
         blank2(j,1,i,2)    = ccmin4(blank(j,1,i),
     .                               blank(j,1,i-1),
     .                               blank(j-1,1,i),
     .                               blank(j-1,1,i-1))
         blank2(j,kdim,i,2) = ccmin4(blank(j,kdim1,i),
     .                               blank(j,kdim1,i-1),
     .                               blank(j-1,kdim1,i),
     .                               blank(j-1,kdim1,i-1))
  613    continue
c
c        edges and corners (k=1 and k=kdim)
c
         do 614 m=1,2
         j  = 1
         jj = 1
         if (m.eq.2) then
            j  = jdim
            jj = jdim1
         end if
         do 614 i=2,idim1
         blank2(j,1,i,2)    = ccmin(blank(jj,1,i),
     .                              blank(jj,1,i-1))
         blank2(j,kdim,i,2) = ccmin(blank(jj,kdim1,i),
     .                              blank(jj,kdim1,i-1))
  614    continue
         do 615 m=1,2
         i  = 1
         ii = 1
         if (m.eq.2) then
            i  = idim
            ii = idim1
         end if
         do 615 j=2,jdim1
         blank2(j,1,i,2)    = ccmin(blank(j,1,ii),
     .                              blank(j-1,1,ii))
         blank2(j,kdim,i,2) = ccmin(blank(j,kdim1,ii),
     .                              blank(j-1,kdim1,ii))
  615   continue
         blank2(jdim,1,idim,2)    = blank2(jdim-1,1,idim-1,2)
         blank2(1,1,1,2)          = blank2(2,1,2,2)
         blank2(jdim,kdim,idim,2) = blank2(jdim-1,kdim,idim-1,2)
         blank2(1,kdim,1,2)       = blank2(2,kdim,2,2)
c
c       interior cells
c
         do 713 i=2,idim1
         do 713 j=2,jdim1
         do 713 k=2,kdim1
         blank2(j,k,i,2) = ccmin8(blank(j,k,i),
     .                            blank(j-1,k,i),
     .                            blank(j,k-1,i),
     .                            blank(j-1,k-1,i),
     .                            blank(j,k,i-1),
     .                            blank(j-1,k,i-1),
     .                            blank(j,k-1,i-1),
     .                            blank(j-1,k-1,i-1))
  713    continue
c
c        combining topology and overset grid blankings
c
         do 813 i=1,idim
         do 813 j=1,jdim
         do 813 k=1,kdim
         blank2(j,k,i,1) = ccmin(blank2(j,k,i,1),blank2(j,k,i,2))
  813    continue
c
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
         if (inpl3d(iplot,2).eq.0) then
            write(11,'(''   plot3dg file is an xyz file'',
     .      '' at grid points'')')
            write(11,'(''   plot3dq file is a    q file'',
     .      '' at grid points'')')
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then 
               write(11,'(''   NOTE: contains dq/d(), rather than q'')') 
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (ifunc.eq.-5) then
            write(11,'(''   plot3dg file is an xyz file'',
     .      '' at grid points'')')
            write(11,'(''   plot3dq file is a function'',
     .      '' file of Cp at grid points'')')
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then 
               write(11,'(''   NOTE: contains dCp/d(), rather than'',
     .                    '' Cp'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (ifunc.eq.-6) then
            write(11,'(''   plot3dg file is an xyz file'',
     .      '' at grid points'')')
            write(11,'(''   plot3dq file is a function'',
     .      '' file of P/Pinf at grid points'')')
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then 
               write(11,'(''   NOTE: contains d(P/Pinf)/d(), rather'',
     .                    '' than P/Pinf'')')
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
      xg(jw,kw,iw,4) = blank2(j,k,i,1)
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
      if (myid.eq.myhost .and. icall2d.eq.0) then
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
               write(93,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
                 write(93,'(5e14.6)')
     .                     (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                     (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                     (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .                 (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         else
            if (iblnk.eq.0) then
               write(93,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(93,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .               (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         end if
      else
         if (i2d.eq.0) then
            if (iblnk.eq.0) then
               write(93)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(93)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .             (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         else
            if (iblnk.eq.0) then
               write(93)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw) 
            else
               write(93)
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
c     determine q values at grid points and load into xv array
c
      jdw = (j2-j1)/j3 + 1
      kdw = (k2-k1)/k3 + 1
      idw = (i2-i1)/i3 + 1
      ldw = 5
c
      call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,q,
     .            qi0,qj0,qk0,jdw,kdw,idw,xv,ldw) 
c
c     ifunc...flag to determine whether a plot3d q file or plot3d function
c             file is written. function files contain only one variable to
c             be compatable with the current version of FAST
c           = 0 plot3d q file is written at grid points
c           < 0 function file is written
c            -5 pressure coefficient (cp) at grid points
c            -6 p/pinf at grid points
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
c        load solution (q) in xv array
c
         do 4000 i=1,iw
         do 4000 k=1,kw
         do 4000 j=1,jw
         xv(j,k,i,5) = xv(j,k,i,5)/gm1+0.5*(xv(j,k,i,2)**2
     .               + xv(j,k,i,3)**2+xv(j,k,i,4)**2)*xv(j,k,i,1)
         xv(j,k,i,2) = xv(j,k,i,1)*xv(j,k,i,2)
         xv(j,k,i,3) = xv(j,k,i,1)*xv(j,k,i,3)
         xv(j,k,i,4) = xv(j,k,i,1)*xv(j,k,i,4)
 4000    continue
c
      else
c
c        load desired function in xv array
c
         if (ifunc.eq.-5) then
            cpc = 2.e0/(gamma*xmach*xmach)
            do 4030 i=1,iw
            do 4030 k=1,kw
            do 4030 j=1,jw
            xv(j,k,i,1) = (xv(j,k,i,5)/p0-1)*cpc
 4030       continue
         end if
c
         if (ifunc.eq.-6) then
            do 4032 i=1,iw
            do 4032 k=1,kw
            do 4032 j=1,jw
            xv(j,k,i,1) = xv(j,k,i,5)/p0
 4032       continue
         end if
c
      end if
c
c     load xv into single precision xw array
c     (for complex code, load either real or complex part
c     of xv into single precision xw array, depending
c     on whether solution or derivative is to be output)
c
      leq = 5
      if (ifunc .ne. 0) leq = 1
c
#     ifdef CMPLX
c
      if (ip3dgrad .eq. 1) then
c        output derivative of solution
         do l=1,leq
            do i=1,iw
               do k=1,kw
                  do j=1,jw
                     xw(j,k,i,l) = imag(xv(j,k,i,l))/real(delh)
                  end do
               end do
            end do
         end do
      else
c        output solution
         do l=1,leq
            do i=1,iw
               do k=1,kw
                  do j=1,jw
                     xw(j,k,i,l) = real(xv(j,k,i,l))
                  end do
               end do
            end do
         end do
      end if
#     else
      do l=1,leq
         do i=1,iw
            do k=1,kw
               do j=1,jw
                  xw(j,k,i,l) = xv(j,k,i,l)
               end do
            end do
         end do
      end do
#     endif
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
      else if (ifunc.lt.0) then
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
            write(94,'(5e14.6)') xmachw,alphww,reuew,timew
            if (i2d.eq.0) then
               if (ialph.eq.0) then
                  write(94,'(5e14.6)')
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
               else
                  do i=1,iw
                     do j=1,jw
                        do k=1,kw
                           xw(j,k,i,3) = -xw(j,k,i,3)
                        end do
                     end do
                  end do
                  write(94,'(5e14.6)')
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .             (((xw(j,k,i,4),i=1,iw),j=1,jw),k=1,kw),
     .             (((xw(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .             (((xw(j,k,i,5),i=1,iw),j=1,jw),k=1,kw) 
               end if
            else
               write(94,'(5e14.6)')
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,5)
            end if
         else
            write(94)  xmachw,alphww,reuew,timew
            if (i2d.eq.0) then
               if (ialph.eq.0) then
                  write(94)
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
               else
                  do i=1,iw
                     do j=1,jw
                        do k=1,kw
                           xw(j,k,i,3) = -xw(j,k,i,3)
                        end do
                     end do
                  end do
                  write(94)
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .             (((xw(j,k,i,4),i=1,iw),j=1,jw),k=1,kw),
     .             (((xw(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .              (((xw(j,k,i,5),i=1,iw),j=1,jw),k=1,kw) 
               end if
            else
               write(94)
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,5)
            end if
         end if
      else
         if (ibin.eq.0) then
            write(94,'(5e14.6)')
     .      (((xw(j,k,i,1),i=1,iw),j=1,jw),k=1,kw)
         else
            write(94)
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
c     NOTE: for the first section of printed output (up to the end of
c           the 5000 loop), turb. viscosity data are stored in the 4th
c           position of the xg array as follows
c
c           xg(j,k,i,4) = vist3d
c
c           the solution q is stored in the xw array as follows
c
c           xw(j,k,i,1) = q(1) (rho)
c           xw(j,k,i,2) = q(2) (u)
c           xw(j,k,i,3) = q(3) (v)
c           xw(j,k,i,4) = q(4) (w)
c           xw(j,k,i,5) = q(5) (p)
c
c           for the subsequent sections of printed output (loops 6002,
c           7002 and 8002), the following storage is used:
c
c           xg(j,k,i,4) = dn (normal spacing away from a surface)
c
c           xw(j,k,i,1) = pres1 (pressure)
c           xw(j,k,i,2) = temp1 (temperature)
c           xw(j,k,i,3) = Cf    (skin friction coeff.)
c           xw(j,k,i,4) = Ch    (heat transfer coeff.)
c           xw(j,k,i,5) = yplus (turbulent length scale y+)
c
c           in all sections xg(1-3) contain x,y,z, as set at the top
c           of the subroutine
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
c     determine q values at grid points and load into single
c     precision array
c
      jdw = (j2-j1)/j3 + 1
      kdw = (k2-k1)/k3 + 1
      idw = (i2-i1)/i3 + 1
      ldw = 5
      call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,q,
     .            qi0,qj0,qk0,jdw,kdw,idw,xv,ldw)          
      do l=1,ldw
         do i=1,idw
            do k=1,kdw
               do j=1,jdw
                  xw(j,k,i,l) = xv(j,k,i,l)
               end do
            end do
         end do
      end do
c
c     determine eddy viscosity values at grid points and load
c     into single precision array
c
      if (ivisc(3).gt.1 .or. ivisc(2).gt.1 .or. ivisc(1).gt.1) then
         ldw = 1
         call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,vist3d,
     .               vi0,vj0,vk0,jdw,kdw,idw,xv(1,1,1,4),ldw)
      do i=1,idw
         do k=1,kdw
            do j=1,jdw
               xg(j,k,i,4) = xv(j,k,i,4)
            end do
         end do
      end do 
      end if
c
c     old way for eddy viscosity ...incorrect for surface data,
c     but keep for checking with version 5.0
c
      iold = 0
      if (iold .gt. 0) then
      if (ivisc(3).gt.1 .or. ivisc(2).gt.1 .or. ivisc(1).gt.1) then
         iw = 0
         do 9100 i=i1,i2,i3
         iw = iw + 1
         jw = 0
         do 9100 j=j1,j2,j3
         jw = jw + 1
         kw = 0
         do 9100 k=k1,k2,k3
         kw = kw + 1
         jv    = j
         kv    = k
         iv    =    i
         jv    = min(jdim-1,jv)
         kv    = min(kdim-1,kv)
         iv    = min(idim-1,iv)
         jx    = j-1
         kx    = k-1
         ix    = i-1
         jx    = max(1,jx)
         kx    = max(1,kx)
         ix    = max(1,ix)
         xg(jw,kw,iw,4) = 0.125*( vist3d(jv,kv,iv) + vist3d(jv,kv,ix)
     .                  + vist3d(jx,kv,iv)   + vist3d(jx,kv,ix)
     .                  + vist3d(jv,kx,iv)   + vist3d(jv,kx,ix)
     .                  + vist3d(jx,kx,iv)   + vist3d(jx,kx,ix) )
 9100    continue
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
   77      format(/3x,1hI,3x,1hJ,3x,1hK,9x,1hX,16x,1hY,16x,1hZ,
     .     14x,6hU/Uinf,11x,6hV/Vinf,11x,6hW/Winf,11x,6hP/Pinf,
     .     11x,6hT/Tinf,12x,4hMACH,10x,6h    cp,12x,9htur. vis.)
         else
           write(17,7)
    7      format(/3x,1hI,3x,1hJ,3x,1hK,9x,1hX,16x,1hY,16x,1hZ,
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
            write(17,29)i,j,k,xg(jw,kw,iw,1),xg(jw,kw,iw,2),
     .                  xg(jw,kw,iw,3),real(q2),real(q3),real(q4),
     .                  real(q5),real(t1),real(xm1),real(pitot),
     .                  real(edvis)
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
   29    format(3i4,11(1x,e16.9))
c
      end if
c
      if (ivisc(3).ne.0) then
c
c     output viscous-flow data (skin friction, heat transfer, etc) on 
c     k=1 and/or k=kdim surfaces
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
      iw = 0
      do 6002 i=i1,i2,i3
      iw = iw + 1
      id   = i
      id1  = id-1
      if (id1.lt.1)   id1 = 1
      if (id.gt.idim1) id = idim1
      kw = 0
      do 6001 k=k1,k2,k3
      kw = kw + 1
      jw = 0
      if (k.ne.1 .and. k.ne.kdim) go to 6001
      do 6000 j=j1,j2,j3
      jw = jw + 1
      if (j.eq.1 .or. j.eq.jdim)  go to 6000
      if (k.eq.1) then
         kd  = 1
         kd1 = kd+1
         kdx = 1
         m   = 2
         sgn = 1.
      else
         kd  = kdim1
         kd1 = kd-1
         kdx = kdim
         m   = 4
         sgn = -1.
      end if
c
c     wall value - average in J and I (after QFACE, m=2 & 4 are interface
c     (wall) values
c
      q1k1 =       .25*(qk0(j,id,1,m)  +qk0(j-1,id,1,m)
     .                + qk0(j,id1,1,m) +qk0(j-1,id1,1,m))
      q2k1 =       .25*(qk0(j,id,2,m)  +qk0(j-1,id,2,m)
     .                + qk0(j,id1,2,m) +qk0(j-1,id1,2,m))
      q3k1 =       .25*(qk0(j,id,3,m)  +qk0(j-1,id,3,m)
     .                + qk0(j,id1,3,m) +qk0(j-1,id1,3,m))
      q4k1 =       .25*(qk0(j,id,4,m)  +qk0(j-1,id,4,m)
     .                + qk0(j,id1,4,m) +qk0(j-1,id1,4,m))
      q5k1 = gamma*.25*(qk0(j,id,5,m)  +qk0(j-1,id,5,m)
     .                + qk0(j,id1,5,m) +qk0(j-1,id1,5,m))
      t1k1 = q5k1/q1k1
c
c     first cell center location - average in J and I
c
      q1k2 =       .25*(q(j,kd,id,1)   +q(j-1,kd,id,1)
     .                + q(j,kd,id1,1)  +q(j-1,kd,id1,1))
      q2k2 =       .25*(q(j,kd,id,2)   +q(j-1,kd,id,2)
     .                + q(j,kd,id1,2)  +q(j-1,kd,id1,2))
      q3k2 =       .25*(q(j,kd,id,3)   +q(j-1,kd,id,3)
     .                + q(j,kd,id1,3)  +q(j-1,kd,id1,3))
      q4k2 =       .25*(q(j,kd,id,4)   +q(j-1,kd,id,4)
     .                + q(j,kd,id1,4)  +q(j-1,kd,id1,4))
      q5k2 = gamma*.25*(q(j,kd,id,5)   +q(j-1,kd,id,5)
     .                + q(j,kd,id1,5)  +q(j-1,kd,id1,5))
      t1k2 = q5k2/q1k2
c
      if (sk(j  ,kdx,id ,4) .eq. 0. .or. 
     +    sk(j-1,kdx,id ,4) .eq. 0. .or.
     +    sk(j  ,kdx,id1,4) .eq. 0. .or.
     +    sk(j-1,kdx,id1,4) .eq. 0.) then
      dx = x(j,kd+1,i)-x(j,kd,i)
      dy = y(j,kd+1,i)-y(j,kd,i)
      dz = z(j,kd+1,i)-z(j,kd,i)
      dn = sqrt(dx**2+dy**2+dz**2)
      else
      dn = (vol(j  ,kd,id )/sk(j  ,kdx,id ,4)+
     +      vol(j-1,kd,id )/sk(j-1,kdx,id ,4)+
     +      vol(j  ,kd,id1)/sk(j  ,kdx,id1,4)+
     +      vol(j-1,kd,id1)/sk(j-1,kdx,id1,4))/4.
      end if
c
c     Get turb viscosity at wall (0 unless wall fn used)
      if (iwf(3) .eq. 1) then
      avgmut =       .25*(vk0(j,id,1,m)  +vk0(j-1,id,1,m)
     .                  + vk0(j,id1,1,m) +vk0(j-1,id1,1,m))
      else
      avgmut = 0.
      end if
      emuka = (t1k1**1.5)*((1.0+cbar/tinf)/(t1k1+cbar/tinf))
c
c     Use component of velocity parallel to wall (with 2-point formula)
      urel = q2k2-q2k1
      vrel = q3k2-q3k1
      wrel = q4k2-q4k1
      sk1  = (sk(j  ,kdx,id ,1)+sk(j-1,kdx,id ,1)+
     +        sk(j  ,kdx,id1,1)+sk(j-1,kdx,id1,1))/4.
      sk2  = (sk(j  ,kdx,id ,2)+sk(j-1,kdx,id ,2)+
     +        sk(j  ,kdx,id1,2)+sk(j-1,kdx,id1,2))/4.
      sk3  = (sk(j  ,kdx,id ,3)+sk(j-1,kdx,id ,3)+
     +        sk(j  ,kdx,id1,3)+sk(j-1,kdx,id1,3))/4.
      vnorm = (urel*sk1 + vrel*sk2 + wrel*sk3)*sgn
      upar = urel-vnorm*sk1*sgn
      vpar = vrel-vnorm*sk2*sgn
      wpar = wrel-vnorm*sk3*sgn
      Cf   = 2.*sqrt(upar**2 + vpar**2 + wpar**2)
      Cf = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn
      if (real(q2k2).lt.0.) cf = -cf
c     Cf vector
      cfx  = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn
      cfy  = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn
      cfz  = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn
c
c     Note:  Ch definition not standard
c
      eps  = 1.0e-6
      tty1 = 1.0 +gm1*0.5*xmach*xmach
      t1m1 = t1k1-tty1
      if (abs(real(t1m1)).le.real(eps)) then
         Ch = 999.99999
      else
c        2-point formula
            Ch = 2.*(t1k2-t1k1)
            Ch = (emuka+avgmut)/(reue*pr*(t1k1-tty1))*Ch/dn
      end if
c
      pres1 = q5k1
      temp1 = t1k1
      yplus = 0.
c     yplus is value at first GRIDPOINT away from wall
      if (abs(real(Cf)).gt.0.) then
         yplus = dn*reue*q1k1*sqrt(ccabs(Cf*0.5/q1k1))/emuka
      end if
c
c     store data in single precision arrays
c
      xg(jw,kw,iw,1) = x(j,k,i)
      xg(jw,kw,iw,2) = y(j,k,i)
      xg(jw,kw,iw,3) = z(j,k,i)
      xg(jw,kw,iw,4) = dn
c
      xw(jw,kw,iw,1) = pres1
      xw(jw,kw,iw,2) = temp1
      xw(jw,kw,iw,3) = Cf
      xw(jw,kw,iw,4) = Ch
      xw(jw,kw,iw,5) = yplus
      xw(jw,kw,iw,6) = cfx
      xw(jw,kw,iw,7) = cfy
      xw(jw,kw,iw,8) = cfz
c
 6000 continue
 6001 continue
 6002 continue
c
#if defined DIST_MPI
      end if
c
c     send/receive printout data
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
c
      nnq = jw*kw*iw*8
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
#endif
c
c     write printout data to file
c
      if (myid.eq.myhost) then
         isay = 0
         iw = 0
         do 6102 i=i1,i2,i3
         iw = iw + 1
         kw = 0
         do 6101 k=k1,k2,k3
         kw = kw + 1
         jw = 0
         if (k.ne.1 .and. k.ne.kdim) go to 6101
         do 6100 j=j1,j2,j3
         jw = jw + 1
         if (j.eq.1 .or. j.eq.jdim)  go to 6100
         isay = isay+1
         if (isay.eq.1) write(17,53)
   53    format(/3x,1hI,3x,1hJ,3x,1hK,9x,1hX,16x,1hY,16x,1hZ,16x,
     .   2hdn,12x,6hP/Pinf,11x,6hT/Tinf,13x,2hCf,15x,2hCh,13x,5hyplus,
     .   13x,3hCfx,14x,3hCfy,14x,3hCfz)
         if (ialph.eq.0) then
            write(17,22)i,j,k,xg(jw,kw,iw,1),xg(jw,kw,iw,2),
     .                  xg(jw,kw,iw,3),xg(jw,kw,iw,4),
     .                  xw(jw,kw,iw,1),xw(jw,kw,iw,2),
     .                  xw(jw,kw,iw,3),xw(jw,kw,iw,4),
     .                  xw(jw,kw,iw,5),xw(jw,kw,iw,6),
     .                  xw(jw,kw,iw,7),xw(jw,kw,iw,8)
         else
           xg(jw,kw,iw,2) = -xg(jw,kw,iw,2)
           xw(jw,kw,iw,7) = -xw(jw,kw,iw,7)
            write(17,22)i,j,k,xg(jw,kw,iw,1),xg(jw,kw,iw,3),
     .                  xg(jw,kw,iw,2),xg(jw,kw,iw,4),
     .                  xw(jw,kw,iw,1),xw(jw,kw,iw,2),
     .                  xw(jw,kw,iw,3),xw(jw,kw,iw,4),
     .                  xw(jw,kw,iw,5),xw(jw,kw,iw,6),
     .                  xw(jw,kw,iw,8),xw(jw,kw,iw,7)
         end if
 6100    continue
 6101    continue
 6102    continue
      end if
c
      end if
c
      if (ivisc(2).ne.0) then
c
c     output viscous-flow data (skin friction, heat transfer, etc) on
c     j=1 and/or j=jdim surfaces
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
      iw = 0
      do 7002 i=i1,i2,i3
      iw = iw + 1
      id   = i
      id1  = id-1
      if (id1.lt.1)   id1 = 1
      if (id.gt.idim1) id = idim1
      jw = 0
      do 7001 j=j1,j2,j3
      jw = jw + 1
      kw = 0
      if (j.ne.1 .and. j.ne.jdim) go to 7001
      do 7000 k=k1,k2,k3
      kw = kw + 1
      if (k.eq.1 .or. k.eq.kdim)  go to 7000
      if (j.eq.1) then
         jd  = 1
         jd1 = jd+1
         jdx = 1
         m   = 2
         sgn = 1.
      else
         jd  = jdim1
         jd1 = jd-1
         jdx = jdim
         m   = 4
         sgn = -1.
      end if
c
c     wall value - average in K and I (after QFACE, m=2 & 4 are interface
c     (wall) values
c
      q1j1 =       .25*(qj0(k,id,1,m)  +qj0(k-1,id,1,m)
     .                + qj0(k,id1,1,m) +qj0(k-1,id1,1,m))
      q2j1 =       .25*(qj0(k,id,2,m)  +qj0(k-1,id,2,m)
     .                + qj0(k,id1,2,m) +qj0(k-1,id1,2,m))
      q3j1 =       .25*(qj0(k,id,3,m)  +qj0(k-1,id,3,m)
     .                + qj0(k,id1,3,m) +qj0(k-1,id1,3,m))
      q4j1 =       .25*(qj0(k,id,4,m)  +qj0(k-1,id,4,m)
     .                + qj0(k,id1,4,m) +qj0(k-1,id1,4,m))
      q5j1 = gamma*.25*(qj0(k,id,5,m)  +qj0(k-1,id,5,m)
     .                + qj0(k,id1,5,m) +qj0(k-1,id1,5,m))
      t1j1 = q5j1/q1j1
c
c     first cell center location - average in K and I
c
      q1j2 =       .25*(q(jd,k,id,1)   +q(jd,k-1,id,1)
     .                + q(jd,k,id1,1)  +q(jd,k-1,id1,1))
      q2j2 =       .25*(q(jd,k,id,2)   +q(jd,k-1,id,2)
     .                + q(jd,k,id1,2)  +q(jd,k-1,id1,2))
      q3j2 =       .25*(q(jd,k,id,3)   +q(jd,k-1,id,3)
     .                + q(jd,k,id1,3)  +q(jd,k-1,id1,3))
      q4j2 =       .25*(q(jd,k,id,4)   +q(jd,k-1,id,4)
     .                + q(jd,k,id1,4)  +q(jd,k-1,id1,4))
      q5j2 = gamma*.25*(q(jd,k,id,5)   +q(jd,k-1,id,5)
     .                + q(jd,k,id1,5)  +q(jd,k-1,id1,5))
      t1j2 = q5j2/q1j2
c
      if (sj(jdx,k  ,id ,4) .eq. 0. .or. 
     +    sj(jdx,k-1,id ,4) .eq. 0. .or.
     +    sj(jdx,k  ,id1,4) .eq. 0. .or.
     +    sj(jdx,k-1,id1,4) .eq. 0.) then
      dx = x(jd+1,k,i)-x(jd,k,i)
      dy = y(jd+1,k,i)-y(jd,k,i)
      dz = z(jd+1,k,i)-z(jd,k,i)
      dn = sqrt(dx**2+dy**2+dz**2)
      else
      dn = (vol(jd,k  ,id )/sj(jdx,k  ,id ,4)+
     +      vol(jd,k-1,id )/sj(jdx,k-1,id ,4)+
     +      vol(jd,k  ,id1)/sj(jdx,k  ,id1,4)+
     +      vol(jd,k-1,id1)/sj(jdx,k-1,id1,4))/4.
      end if
c
c     Get turb viscosity at wall (0 unless wall fn used)
      if (iwf(2) .eq. 1) then
      avgmut =       .25*(vj0(k,id,1,m)  +vj0(k-1,id,1,m)
     .                  + vj0(k,id1,1,m) +vj0(k-1,id1,1,m))
      else
      avgmut = 0.
      end if
      emuka = (t1j1**1.5)*((1.0+cbar/tinf)/(t1j1+cbar/tinf))
c
c     Use component of velocity parallel to wall
      urel = q2j2-q2j1
      vrel = q3j2-q3j1
      wrel = q4j2-q4j1
      sj1  = (sj(jdx,k  ,id ,1)+sj(jdx,k-1,id ,1)+
     +        sj(jdx,k  ,id1,1)+sj(jdx,k-1,id1,1))/4.
      sj2  = (sj(jdx,k  ,id ,2)+sj(jdx,k-1,id ,2)+
     +        sj(jdx,k  ,id1,2)+sj(jdx,k-1,id1,2))/4.
      sj3  = (sj(jdx,k  ,id ,3)+sj(jdx,k-1,id ,3)+
     +        sj(jdx,k  ,id1,3)+sj(jdx,k-1,id1,3))/4.
      vnorm = (urel*sj1 + vrel*sj2 + wrel*sj3)*sgn
      upar = urel-vnorm*sj1*sgn
      vpar = vrel-vnorm*sj2*sgn
      wpar = wrel-vnorm*sj3*sgn
      Cf   = 2.*sqrt(upar**2 + vpar**2 + wpar**2)
      Cf = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn
      if (real(q2j2).lt.0.) cf = -cf
c     Cf vector
      cfx  = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn
      cfy  = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn
      cfz  = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn
c
c     Note:  Ch definition not standard
c
      eps  = 1.0e-6
      tty1 = 1.0 +gm1*0.5*xmach*xmach
      t1m1 = t1j1-tty1
      if (abs(real(t1m1)).le.real(eps)) then
         Ch = 999.99999
      else
c        2-point formula
            Ch = 2.*(t1j2-t1j1)
            Ch = (emuka+avgmut)/(reue*pr*(t1j1-tty1))*Ch/dn
      end if
c
      pres1 = q5j1
      temp1 = t1j1
      yplus = 0.
c     yplus is valus at first GRIDPOINT away from wall
      if (abs(real(Cf)).gt.0.) then
         yplus = dn*reue*q1j1*sqrt(ccabs(Cf*0.5/q1j1))/emuka
      end if
c
c     store data in single precision arrays
c
      xg(jw,kw,iw,1) = x(j,k,i)
      xg(jw,kw,iw,2) = y(j,k,i)
      xg(jw,kw,iw,3) = z(j,k,i)
      xg(jw,kw,iw,4) = dn
c
      xw(jw,kw,iw,1) = pres1
      xw(jw,kw,iw,2) = temp1
      xw(jw,kw,iw,3) = Cf
      xw(jw,kw,iw,4) = Ch
      xw(jw,kw,iw,5) = yplus
      xw(jw,kw,iw,6) = cfx
      xw(jw,kw,iw,7) = cfy
      xw(jw,kw,iw,8) = cfz
c
 7000 continue
 7001 continue
 7002 continue
c
#if defined DIST_MPI
      end if
c
c     send/receive printout data
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
c
      nnq = jw*kw*iw*8
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
#endif
c     write printout data to file
c
      if (myid.eq.myhost) then
         isay = 0
         iw = 0
         do 7102 i=i1,i2,i3
         iw = iw + 1
         jw = 0
         do 7101 j=j1,j2,j3
         jw = jw + 1
         kw = 0
         if (j.ne.1 .and. j.ne.jdim) go to 7101
         do 7100 k=k1,k2,k3
         kw = kw + 1
         if (k.eq.1 .or. k.eq.kdim)  go to 7100
         isay = isay+1
         if (isay.eq.1) write(17,53)
         if (ialph.eq.0) then
            write(17,22)i,j,k,xg(jw,kw,iw,1),xg(jw,kw,iw,2),
     .                  xg(jw,kw,iw,3),xg(jw,kw,iw,4),
     .                  xw(jw,kw,iw,1),xw(jw,kw,iw,2),
     .                  xw(jw,kw,iw,3),xw(jw,kw,iw,4),
     .                  xw(jw,kw,iw,5),xw(jw,kw,iw,6),
     .                  xw(jw,kw,iw,7),xw(jw,kw,iw,8)
         else
            xg(jw,kw,iw,2) = -xg(jw,kw,iw,2)
            xw(jw,kw,iw,7) = -xw(jw,kw,iw,7)
            write(17,22)i,j,k,xg(jw,kw,iw,1),xg(jw,kw,iw,3),
     .                  xg(jw,kw,iw,2),xg(jw,kw,iw,4),
     .                  xw(jw,kw,iw,1),xw(jw,kw,iw,2),
     .                  xw(jw,kw,iw,3),xw(jw,kw,iw,4),
     .                  xw(jw,kw,iw,5),xw(jw,kw,iw,6),
     .                  xw(jw,kw,iw,8),xw(jw,kw,iw,7)
         end if
 7100    continue
 7101    continue
 7102    continue
      end if
c
      end if
c
      if (idim.eq.2) return
      if (ivisc(1).ne.0) then
c
c     output viscous-flow data (skin friction, heat transfer, etc) on
c     i=1 and/or i=idim surfaces
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
      jw = 0
      do 8002 j=j1,j2,j3
      jw = jw + 1
      iw = 0
      if (j.eq.1 .or. j.eq.jdim)  go to 8002
      do 8001 i=i1,i2,i3
      iw = iw + 1
      kw = 0
      if (i.ne.1 .and. i.ne.idim) go to 8001
      do 8000 k=k1,k2,k3
      kw = kw + 1
      if (k.eq.1 .or. k.eq.kdim)  go to 8000
      if (i.eq.1) then
         id  = 1
         id1 = id+1
         idx = 1
         m   = 2
         sgn = 1.
      else
         id  = idim1
         id1 = id-1
         idx = idim
         m   = 4
         sgn = -1.
      end if
c
c     wall value - average in J and k (after QFACE, m=2 & 4 are interface
c     (wall) values
c
      q1i1 =       .25*(qi0(j,k,1,m)   +qi0(j-1,k,1,m)
     .                + qi0(j,k-1,1,m) +qi0(j-1,k-1,1,m))
      q2i1 =       .25*(qi0(j,k,2,m)   +qi0(j-1,k,2,m)
     .                + qi0(j,k-1,2,m) +qi0(j-1,k-1,2,m))
      q3i1 =       .25*(qi0(j,k,3,m)   +qi0(j-1,k,3,m)
     .                + qi0(j,k-1,3,m) +qi0(j-1,k-1,3,m))
      q4i1 =       .25*(qi0(j,k,4,m)   +qi0(j-1,k,4,m)
     .                + qi0(j,k-1,4,m) +qi0(j-1,k-1,4,m))
      q5i1 = gamma*.25*(qi0(j,k,5,m)   +qi0(j-1,k,5,m)
     .                + qi0(j,k-1,5,m) +qi0(j-1,k-1,5,m))
      t1i1 = q5i1/q1i1
c
c     first cell center location - average in J and K
c
      q1i2 =       .25*(q(j,k,id,1)    +q(j,k-1,id,1)
     .                + q(j-1,k,id,1)  +q(j-1,k-1,id,1))
      q2i2 =       .25*(q(j,k,id,2)    +q(j,k-1,id,2)
     .                + q(j-1,k,id,2)  +q(j-1,k-1,id,2))
      q3i2 =       .25*(q(j,k,id,3)    +q(j,k-1,id,3)
     .                + q(j-1,k,id,3)  +q(j-1,k-1,id,3))
      q4i2 =       .25*(q(j,k,id,4)    +q(j,k-1,id,4)
     .                + q(j-1,k,id,4)  +q(j-1,k-1,id,4))
      q5i2 = gamma*.25*(q(j,k,id,5)    +q(j,k-1,id,5)
     .                + q(j-1,k,id,5)  +q(j-1,k-1,id,5))
      t1i2 = q5i2/q1i2
c
      if (si(j  ,k  ,idx,4) .eq. 0. .or.                 
     +    si(j  ,k-1,idx,4) .eq. 0. .or.
     +    si(j-1,k  ,idx,4) .eq. 0. .or.
     +    si(j-1,k-1,idx,4) .eq. 0.) then
      dx = x(j,k,id+1)-x(j,k,id)
      dy = y(j,k,id+1)-y(j,k,id)
      dz = z(j,k,id+1)-z(j,k,id)
      dn = sqrt(dx**2+dy**2+dz**2)
      else
      dn = (vol(j  ,k  ,id)/si(j  ,k  ,idx,4)+
     +      vol(j  ,k-1,id)/si(j  ,k-1,idx,4)+
     +      vol(j-1,k  ,id)/si(j-1,k  ,idx,4)+
     +      vol(j-1,k-1,id)/si(j-1,k-1,idx,4))/4.
      end if
c
c     Get turb viscosity at wall (0 unless wall fn used)
      if (iwf(1) .eq. 1) then
      avgmut =       .25*(vi0(j,k,1,m)   +vi0(j-1,k,1,m)
     .                  + vi0(j,k-1,1,m) +vi0(j-1,k-1,1,m))
      else
      avgmut = 0.
      end if
      emuka = (t1i1**1.5)*((1.0+cbar/tinf)/(t1i1+cbar/tinf))
c
c     Use component of velocity parallel to wall
      urel = q2i2-q2i1
      vrel = q3i2-q3i1
      wrel = q4i2-q4i1
      si1  = (si(j  ,k  ,idx,1)+si(j  ,k-1,idx,1)+
     +        si(j-1,k  ,idx,1)+si(j-1,k-1,idx,1))/4.
      si2  = (si(j  ,k  ,idx,2)+si(j  ,k-1,idx,2)+
     +        si(j-1,k  ,idx,2)+si(j-1,k-1,idx,2))/4.
      si3  = (si(j  ,k  ,idx,3)+si(j  ,k-1,idx,3)+
     +        si(j-1,k  ,idx,3)+si(j-1,k-1,idx,3))/4.
      vnorm = (urel*si1 + vrel*si2 + wrel*si3)*sgn
      upar = urel-vnorm*si1*sgn
      vpar = vrel-vnorm*si2*sgn
      wpar = wrel-vnorm*si3*sgn
      Cf   = 2.*sqrt(upar**2 + vpar**2 + wpar**2)
      Cf = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn
      if (real(q2i2).lt.0.) cf = -cf
c     Cf vector
      cfx  = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn
      cfy  = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn
      cfz  = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn
c
c     Note:  Ch definition not standard
c
      eps  = 1.0e-6
      tty1 = 1.0 +gm1*0.5*xmach*xmach
      t1m1 = t1i1-tty1
      if (abs(real(t1m1)).le.real(eps)) then
         Ch = 999.99999
      else
c        2-point formula
            Ch = 2.*(t1i2-t1i1)
            Ch = (emuka+avgmut)/(reue*pr*(t1i1-tty1))*Ch/dn
      end if
c
      pres1 = q5i1
      temp1 = t1i1
      yplus = 0.
c     yplus is valus at first GRIDPOINT away from wall
      if (abs(real(Cf)).gt.0.) then
         yplus = dn*reue*q1i1*sqrt(ccabs(Cf*0.5/q1i1))/emuka
      end if
c
c     store data in single precision arrays
c
      xg(jw,kw,iw,1) = x(j,k,i)
      xg(jw,kw,iw,2) = y(j,k,i)
      xg(jw,kw,iw,3) = z(j,k,i)
      xg(jw,kw,iw,4) = dn
c
      xw(jw,kw,iw,1) = pres1
      xw(jw,kw,iw,2) = temp1
      xw(jw,kw,iw,3) = Cf
      xw(jw,kw,iw,4) = Ch
      xw(jw,kw,iw,5) = yplus
      xw(jw,kw,iw,6) = cfx
      xw(jw,kw,iw,7) = cfy
      xw(jw,kw,iw,8) = cfz
c
 8000 continue
 8001 continue
 8002 continue
#if defined DIST_MPI
      end if
c
c     send/receive printout data
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
c
      nnq = jw*kw*iw*8
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
#endif
c
c     write printout data to file
c
      if (myid.eq.myhost) then
         isay = 0
         jw = 0
         do 8102 j=j1,j2,j3
         jw = jw + 1
         iw = 0
         if (j.eq.1 .or. j.eq.jdim)  go to 8102
         do 8101 i=i1,i2,i3
         iw = iw + 1
         kw = 0
         if (i.ne.1 .and. i.ne.idim) go to 8101
         do 8100 k=k1,k2,k3
         kw = kw + 1
         if (k.eq.1 .or. k.eq.kdim)  go to 8100
         isay = isay+1
         if (isay.eq.1) write(17,53)
         if (ialph.eq.0) then
            write(17,22)i,j,k,xg(jw,kw,iw,1),xg(jw,kw,iw,2),
     .                  xg(jw,kw,iw,3),xg(jw,kw,iw,4),
     .                  xw(jw,kw,iw,1),xw(jw,kw,iw,2),
     .                  xw(jw,kw,iw,3),xw(jw,kw,iw,4),
     .                  xw(jw,kw,iw,5),xw(jw,kw,iw,6),
     .                  xw(jw,kw,iw,7),xw(jw,kw,iw,8)
         else
            xg(jw,kw,iw,2) = -xg(jw,kw,iw,2)
            xw(jw,kw,iw,7) = -xw(jw,kw,iw,7)
            write(17,22)i,j,k,xg(jw,kw,iw,1),xg(jw,kw,iw,3),
     .                  xg(jw,kw,iw,2),xg(jw,kw,iw,4),
     .                  xw(jw,kw,iw,1),xw(jw,kw,iw,2),
     .                  xw(jw,kw,iw,3),xw(jw,kw,iw,4),
     .                  xw(jw,kw,iw,5),xw(jw,kw,iw,6),
     .                  xw(jw,kw,iw,8),xw(jw,kw,iw,7)
         end if
 8100    continue
 8101    continue
 8102    continue
      end if
c
      end if
c
   22 format(3i4,12(1x,e16.9))
c
      end if
c
      return
      end
