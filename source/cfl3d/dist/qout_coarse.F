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
      subroutine qout_coarse(iseq,lw,lw2,w,mgwk,wk,nwork,nplots,
     .                iovrlp,iibg,                              
     .                kkbg,jjbg,ibcg,lbg,ibpntsg,qb,lwdat,nbci0,
     .                nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,jbcinfo,
     .                kbcinfo,ibcinfo,bcfilei,bcfilej,bcfilek,
     .                itrans,irotat,idefrm,nblock,levelg,igridg,iviscg,
     .                jdimg,kdimg,idimg,nblg,clw,ncycmax,nplot3d,
     .                inpl3d,ip3dsurf,nprint,inpr,iadvance,mycomm,
     .                myid,myhost,mblk2nd,nou,bou,nbuf,ibufdim,maxbl,
     .                maxgr,maxseg,iitot,jsg,ksg,isg,jeg,keg,ieg,
     .                ninter,windex,iindex,nblkpt,intmax,nsub1,maxxe,
     .                nblk,nbli,limblk,isva,nblon,mxbli,thetay,
     .                ip3ddim,nmap,iwk,iwork,xorig,yorig,zorig,
     .                period_miss,geom_miss,epsc0,epsrot,isav_blk,
     .                isav_pat,isav_pat_b,isav_emb,isav_prd,
     .                lbcprd,lbcemb,
     .                dthetxx,dthetyy,dthetzz,nblcg,lfgm,istat2_bl,
     .                istat2_pa,istat2_pe,istat2_em,istat_size,
     .                vormax,ivmax,jvmax,kvmax,bcfiles,mxbcfil,iprnsurf,
     .                nt,mov_coarse,iinc_coarse,jinc_coarse,kinc_coarse,
     .                nummem)

c   Using existing plot3d machinery involves arrays ip3ddim(..,nplots)
c             and, therefore, nplots needs to be .ge. number of blocks
c             Also, iinc, jinc, etc. are global constants, ie, same for all blks
c             Except, we halve it for blocks where that incr is not permissible
c             => increments should be multiples of 2
c
c     NOTE: pass through inpl3d uotouched to plot3d; this will result in a 
c           duplicate message of the type "plot3dg file written as grid points"
c           But, that's probably acceptable
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Output data for plotting or printing.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
#if defined DIST_MPI
#     include "mpif.h"
#   ifdef DBLE_PRECSN
#     define MY_MPI_REAL MPI_DOUBLE_PRECISION
#   else
#     define MY_MPI_REAL MPI_REAL
#   endif
#endif
c
c     maxbl   - maximum number of blocks
c     maxgr   - maximum number of grids
c     nplots  - maximum number of data sets to output via PLOT3D or print
c     ncycmax - maximum number of time-steps/cycles
c
      character*80  bcfiles(mxbcfil)
      character*120 bou(ibufdim,nbuf)
c
      integer bcfilei,bcfilej,bcfilek
c
      dimension nou(nbuf)
      dimension istat2_bl(istat_size,mxbli*5),
     .          istat2_pa(istat_size,intmax*nsub1*3),
     .          istat2_em(istat_size,lbcemb*3),
     .          istat2_pe(istat_size,lbcprd*5)
      dimension vormax(maxbl),ivmax(maxbl),jvmax(maxbl),kvmax(maxbl)
      dimension w(mgwk),wk(nwork),lw(65,maxbl),lwdat(maxbl,maxseg,6)
      dimension nmap(maxbl),mblk2nd(maxbl),iwk(iwork),lw2(43,maxbl)
      dimension iovrlp(maxbl),iibg(iitot),kkbg(iitot),jjbg(iitot),
     .          ibcg(iitot),lbg(maxbl),ibpntsg(maxbl,4),qb(iitot,5,3)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2),
     .          bcfilei(maxbl,maxseg,2),bcfilej(maxbl,maxseg,2),
     .          bcfilek(maxbl,maxseg,2)
      dimension itrans(maxbl),irotat(maxbl),levelg(maxbl),igridg(maxbl),
     .          iviscg(maxbl,3),jdimg(maxbl),kdimg(maxbl),idimg(maxbl),
     .          nblcg(maxbl),nblg(maxgr),clw(ncycmax),inpl3d(nplots,11),
     .          inpr(nplots,11),iadvance(maxbl),thetay(maxbl),
     .          idefrm(maxbl)
      dimension jsg(maxbl),ksg(maxbl),isg(maxbl),jeg(maxbl),keg(maxbl),
     .          ieg(maxbl)
      dimension windex(maxxe,2),iindex(intmax,6*nsub1+9),nblkpt(maxxe)
      dimension dthetxx(intmax,nsub1),dthetyy(intmax,nsub1),
     .          dthetzz(intmax,nsub1)
      dimension nblk(2,mxbli),limblk(2,6,mxbli),isva(2,2,mxbli),
     .          nblon(mxbli)
      dimension geom_miss(2*mxbli),period_miss(lbcprd)
      dimension isav_prd(lbcprd,12)
      dimension isav_blk(2*mxbli,17)
      dimension isav_emb(lbcemb,12)
      dimension isav_pat(intmax,17),isav_pat_b(intmax,nsub1,6)
      dimension xorig(maxbl),yorig(maxbl),zorig(maxbl)
      dimension ip3ddim(3,nplots)
c
      common /bin/ ibin,iblnk,iblnkfr,ip3dgrad
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
      common /maxiv/ ivmx
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /twod/ i2d
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
      common /moovcrs2d/icallcrs,icall2d
      common /conversion/ radtodeg
c
c***********************************************************************
c     write output file in PLOT3D format
c     in this routine, iptype = inpl3d(n,2), where n is the block number 
c     iptype = 0...output q at grid points 
c              1...output q cell centers 
c              2...turbulence data at cell centers, output in place of
c                  the standard plot3d q vector 
c              3...smin at cell centers 
c              4...vist3d at cell centers 
c              5...cp at cell centers 
c             -5...cp at grid points 
c              6...p/pinf at cell centers 
c             -6...p/pinf at grid points 
c***********************************************************************
c
c
      if (1.gt.0) then                !  do this regardless of nplot3d
c
c        update all bc's for level lglobal (and above, if embeded),
c        even if not all blocks at that level get output to plot3d or
c        printout files - this is different from older versions of the
c        code, in which bc's are only called if needed to get the latest
c        values for output.
c
         level_sav = level
         lglobal = lfgm-(mseq-iseq)
c
         do level=lglobal,levelt(iseq)
c

            if (iipv.gt.0) then
               nttuse=max(ntt-1,1)
               clwuse = clw(nttuse)
            end if

#if defined(DIST_MPI)
            call MPI_Bcast(clwuse,1,MY_MPI_REAL,myhost,
     .                     mycomm,ierr)
            if (myid.ne.myhost) then
c
#endif
            do nbl=1,nblock
c              need to call even for blocks not advanced
               if (level.eq.levelg(nbl)) then
                  if (mblk2nd(nbl).eq.myid) then
                     call lead(nbl,lw,lw2,maxbl)
                     call bc(1,nbl,lw,lw2,w,mgwk,wk,nwork,
     .                       clwuse,nou,bou,nbuf,ibufdim,maxbl,
     .                       maxgr,maxseg,itrans,irotat,idefrm,
     .                       igridg,nblg,nbci0,nbcj0,nbck0,nbcidim,
     .                       nbcjdim,nbckdim,ibcinfo,jbcinfo,
     .                       kbcinfo,bcfilei,bcfilej,bcfilek,
     .                       lwdat,myid,idimg,jdimg,kdimg,bcfiles,
     .                       mxbcfil,nummem)
                  end if
               end if
            end do
c
c           update periodic boundary conditions
c
            lres   = 1
            nsafe  = nwork-lres+1
            mneed  = lbcprd*5
            iwk1   = 1
            iwk2   = iwk1 + mneed
            iwk3   = iwk2 + mneed
            iwk4   = iwk3 + mneed
            iwk5   = iwk4 + mneed*2
            iwk6   = iwk5 + mneed
            iwork1 = iwork - iwk6
            if (iwork1.lt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''stopping...not enough integer '',
     .                        ''work space for subroutine bc_period'')')
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''have, need = '',2i12)') 
     .         iwork,iwk6
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            do iii = 1,iwk6
               iwk(iii) = 0
            end do
            call bc_period(1,nbl,lw,lw2,w,mgwk,wk(lres),nsafe,maxbl,
     .                     maxgr,maxseg,iadvance,bcfilei,bcfilej,
     .                     bcfilek,lwdat,xorig,yorig,zorig,jdimg,kdimg,
     .                     idimg,lbcprd,isav_prd,
     .                     period_miss,epsrot,iwk(iwk1),iwk(iwk2),
     .                     iwk(iwk3),iwk(iwk4),iwk(iwk5),myid,myhost,
     .                     mycomm,mblk2nd,nou,bou,nbuf,ibufdim,
     .                     istat2_pe,istat_size,bcfiles,mxbcfil,nummem)
c
c           update embeded-grid boundary conditions
c
            lres   = 1
            nsafe  = nwork-lres+1
            mneed  = lbcemb*3
            iwk1   = 1
            iwk2   = iwk1 + mneed
            iwk3   = iwk2 + mneed
            iwk4   = iwk3 + mneed
            iwk5   = iwk4 + mneed*2
            iwk6   = iwk5 + mneed
            iwork1 = iwork - iwk6
            if (iwork1.lt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''stopping...not enough integer '',
     .                        ''work space for subroutine bc_embed'')')
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''have, need = '',2i12)')
     .         iwork,iwk6
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            do iii = 1,iwk6
               iwk(iii) = 0
            end do
            call bc_embed(1,nbl,lw,lw2,w,mgwk,wk(lres),nsafe,maxbl,
     .                    maxgr,lbcemb,iadvance,idimg,jdimg,
     .                    kdimg,isav_emb,iwk(iwk1),
     .                    iwk(iwk2),iwk(iwk3),iwk(iwk4),iwk(iwk5),
     .                    myid,myhost,mycomm,mblk2nd,nou,bou,nbuf,
     .                    ibufdim,iviscg,istat2_em,istat_size,nummem)
c
c           update 1-1 block boundary conditions
c
            lres   = 1
            nsafe  = nwork-lres+1
            mneed  = mxbli*5
            iwk1   = 1
            iwk2   = iwk1 + mneed
            iwk3   = iwk2 + mneed
            iwk4   = iwk3 + mneed
            iwk5   = iwk4 + mneed*2
            iwk6   = iwk5 + mneed
            iwork1 = iwork - iwk6
            if (iwork1.lt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''stopping...not enough integer '',
     .                        ''work space for subroutine bc_blkint'')')
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''have, need = '',2i12)')
     .         iwork,iwk6
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            do iii = 1,iwk6
               iwk(iii) = 0
            end do
            call bc_blkint(1,nbl,lw,lw2,w,mgwk,wk(lres),nsafe,maxbl,
     .                     maxgr,mxbli,iadvance,geom_miss,epsc0,nblk,
     .                     nbli,limblk,isva,nblon,jdimg,kdimg,idimg,
     .                     mblk2nd,isav_blk,iwk(iwk1),
     .                     iwk(iwk2),iwk(iwk3),iwk(iwk4),iwk(iwk5),
     .                     nou,bou,nbuf,ibufdim,myid,myhost,mycomm,
     .                     istat2_bl,istat_size,nummem)
c
c           update patch-grid boundary conditions
c
            lres   = 1
            nsafe  = nwork-lres+1
            mneed  = intmax*nsub1*3
            iwk1   = 1
            iwk2   = iwk1 + mneed
            iwk3   = iwk2 + mneed
            iwk4   = iwk3 + mneed
            iwk5   = iwk4 + mneed*2
            iwk6   = iwk5 + mneed*2
            iwork1 = iwork - iwk6
            if (iwork1.lt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''stopping...not enough integer'',
     .                       '' work space for subroutine bc_patch'')')
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''have, need = '',2i12)')
     .         iwork,iwk6
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            do iii = 1,iwk6
               iwk(iii) = 0
            end do
            call bc_patch(1,nbl,lw,lw2,w,mgwk,wk(lres),nsafe,maxbl,
     .                    maxgr,intmax,nsub1,maxxe,iadvance,jdimg,kdimg,
     .                    idimg,ninter,windex,iindex,nblkpt,dthetxx,
     .                    dthetyy,dthetzz,isav_pat,isav_pat_b,
     .                    iwk(iwk1),iwk(iwk2),iwk(iwk3),
     .                    iwk(iwk4),iwk(iwk5),myid,myhost,mycomm,
     .                    mblk2nd,nou,bou,nbuf,ibufdim,
     .                    istat2_pa,istat_size,nummem)
c
c           update chimera boundary conditions
c
c           don't update interior (fringe) points if cell center data 
c           is to be output (this is needed to get cell center output
c           identical to version 5)
c
            int_updt = 1
            !do n=1,nplot3d                   ! only do grid point output
            !   if (inpl3d(n,2).gt.0) then
            !      int_updt = 0
            !   end if
            !end do
c
            do nbl=1,nblock
               if (iadvance(nbl).ge.0) then
                  if (level.eq.levelg(nbl)) then
                     if (mblk2nd(nbl).eq.myid) then
                        call lead(nbl,lw,lw2,maxbl)
                        lres  = 1
                        nsafe = nwork-lres+1
                        call bc_xmera(1,nbl,lw,lw2,w,mgwk,wk(lres),
     .                                nsafe,maxbl,iitot,iviscg,iovrlp,
     .                                lbg,ibpntsg,qb,iibg,kkbg,jjbg,
     .                                ibcg,nou,bou,nbuf,ibufdim,
     .                                int_updt,nummem)
                     end if
                  end if
               end if
            end do
c
c           call qface to install face-center values in 
c           the qi0/qj0/qk0 arrays
c
            do nbl=1,nblock
               if (level.eq.levelg(nbl)) then
                  if (mblk2nd(nbl).eq.myid) then
                     call lead(nbl,lw,lw2,maxbl)
                     ldim = 5
                     call qface(jdim,kdim,idim,w(lq),w(lqj0),
     .                          w(lqk0),w(lqi0),w(lbcj),w(lbck),
     .                          w(lbci),w(lblk),ldim)
                     if (ivmx.ge.2) then
                        ldim = 1
                        call qface(jdim,kdim,idim,w(lvis),w(lvj0),
     .                             w(lvk0),w(lvi0),w(lbcj),w(lbck),
     .                             w(lbci),w(lblk),ldim)
                     end if
                  end if
               end if
            end do
c
#if defined(DIST_MPI)
            end if
c
#endif
         end do
c
         level = level_sav
c
      end if
c
      !if (nplot3d.le.0) go to 231                                  
c
      ncount    = 0
      np3d      = nblock/lglobal                                    
      !np3d      = nplot3d                                          
      ifunc     = 0
c
c     if zone has function file output, all must, and all
c     must have the same function output
c
      !if (abs(inpl3d(1,2)).gt.2) ifunc = inpl3d(1,2)
c
      do 60 n=1,np3d                                                
      !do 60 n=1,nplot3d                                            
      if (n.eq.1) then
         if (myid .eq. myhost) then
            if (ibin.eq.0) then
               if (icallcrs.eq.0) write(101,'(3i5)') np3d           
               write(92,'(3i5)') np3d
            else
               if (icallcrs.eq.0) write(101) np3d                   
               write(92) np3d
            end if
         end if
      end if
c
      !nbl = inpl3d(n,1)                                            
      nbl = lglobal*(n-1)+1                                         
c
      if (nbl.gt.nblock) then
         if (myid .eq. myhost) then
            write(11,777)nbl
  777       format(6h Block,i3,18h does not exist.  ,
     .             25hNo plot3d output printed.)
         end if
         go to 60
      end if
c
      ncount = ncount+1
      call lead(nbl,lw,lw2,maxbl)
      !if incr is too large, allow 1 level of halving it any given block 
      i1 = 1                                       !inpl3d(n,3)     
      i2 = idim                                    !inpl3d(n,4)     
      i3 = iinc_coarse                             !inpl3d(n,5)     
      if((iinc_coarse.gt.1).and.(mod(idim-1,iinc_coarse).ne.0))
     .  i3=iinc_coarse/2                                            
      j1 = 1                                       !inpl3d(n,6)     
      j2 = jdim                                    !inpl3d(n,7)     
      j3 = jinc_coarse                             !inpl3d(n,8)     
      if((jinc_coarse.gt.1).and.(mod(jdim-1,jinc_coarse).ne.0))
     .  j3=jinc_coarse/2                                            
      k1 = 1                                       !inpl3d(n,9)     
      k2 = kdim                                    !inpl3d(n,10)    
      k3 = kinc_coarse                             !inpl3d(n,11)    
      if((kinc_coarse.gt.1).and.(mod(kdim-1,kinc_coarse).ne.0))
     .  k3=kinc_coarse/2                                            
      !if (inpl3d(n,2).gt.0) then     ! only use grid pts for coarse movie
c     !   cell center dimensions      
      !   call lead(nbl,lw,lw2,maxbl) 
      !   i2 = min(idim-1,i2)         
      !   j2 = min(jdim-1,j2)         
      !   k2 = min(kdim-1,k2)         
      !   i1 = min(idim-1,i1)         
      !   j1 = min(jdim-1,j1)         
      !   k1 = min(kdim-1,k1)         
      !end if                         
      ip3ddim(1,ncount) = (i2-i1)/i3+1
      ip3ddim(2,ncount) = (j2-j1)/j3+1
      ip3ddim(3,ncount) = (k2-k1)/k3+1
   60 continue
c
      if (myid .eq. myhost) then
         if (ibin.eq.0) then
            if (i2d.eq.0) then
               if (icallcrs.eq.0) 
     .         write(101,'(3i5)') ((ip3ddim(i,n),i=1,3),n=1,ncount)    
               if (ifunc.eq.0) then
                  write(92,'(3i5)')((ip3ddim(i,n),i=1,3),n=1,ncount)
               else
                  write(92,'(3i5)')((ip3ddim(i,n),i=1,3),1,n=1,ncount)
               end if
            else
               if (icallcrs.eq.0) 
     .         write(101,'(3i5)') ((ip3ddim(i,n),i=2,3),n=1,ncount)    
               if (ifunc.eq.0) then
                  write(92,'(3i5)') ((ip3ddim(i,n),i=2,3),n=1,ncount)
               else
                  write(92,'(3i5)') ((ip3ddim(i,n),i=2,3),1,n=1,ncount)
               end if
            end if
         else
            if (i2d.eq.0) then
               if (icallcrs.eq.0) 
     .         write(101) ((ip3ddim(i,n),i=1,3),n=1,ncount)            
               if (ifunc.eq.0) then
                  write(92) ((ip3ddim(i,n),i=1,3),n=1,ncount)
               else
                  write(92) ((ip3ddim(i,n),i=1,3),1,n=1,ncount)
               end if
            else
               if (icallcrs.eq.0) 
     .         write(101) ((ip3ddim(i,n),i=2,3),n=1,ncount)            
               if (ifunc.eq.0) then
                  write(92) ((ip3ddim(i,n),i=2,3),n=1,ncount)
               else
                  write(92) ((ip3ddim(i,n),i=2,3),1,n=1,ncount)
               end if
            end if
         end if
      end if
c
c     correspondence between global block number and plot3d block
c     number stored in nmap(n) for n=1,nblock
c
      do 1097 n = 1,nblock
      nmap(n) = 1
      do 1098 nnn = 1,np3d                                         
      !do 1098 nnn = 1,nplot3d                                     
      m = lglobal*(nnn-1)+1                                        
      !m = inpl3d(nnn,1)                                           
      if (n.eq.m) nmap(n) = nnn
 1098 continue
 1097 continue             
c
      if (myid.eq.myhost) then
         if (lhdr .gt. 0) write(11,1096)
      end if
 1096 format(1h )
c
      do 70 n=1,np3d                                               
      !do 70 n=1,nplot3d                                           
c
      nbl = lglobal*(n-1)+1                                        
      !nbl = inpl3d(n,1)                                           
c
      if (nbl.gt.nblock .or. nbl.le. 0) go to 70
c
      if (iblnkfr .eq. 0) then
c
c        temporarily set blank values at fringe points (not holes)
c        to 1 for plotting purposes (helps reduce gaps in plots)
c
         call lead(nbl,lw,lw2,maxbl)
         if (mblk2nd(nbl).eq.myid .and. iovrlp(nbl).ne.0) then
            call blnkfr(nbl,iibg,kkbg,jjbg,ibpntsg,lbg,iitot,w(lblk),
     .                  jdim,kdim,idim,maxbl,1.)
         end if
      end if
c
      ncount = ncount+1
      call lead(nbl,lw,lw2,maxbl)
      i1 = 1                                       !inpl3d(n,3)     
      i2 = idim                                    !inpl3d(n,4)     
      i3 = iinc_coarse                             !inpl3d(n,5)     
      if((iinc_coarse.gt.1).and.(mod(idim-1,iinc_coarse).ne.0))
     .  i3=iinc_coarse/2                                            
      j1 = 1                                       !inpl3d(n,6)     
      j2 = jdim                                    !inpl3d(n,7)     
      j3 = jinc_coarse                             !inpl3d(n,8)     
      if((jinc_coarse.gt.1).and.(mod(jdim-1,jinc_coarse).ne.0))
     .  j3=jinc_coarse/2                                            
      k1 = 1                                       !inpl3d(n,9)     
      k2 = kdim                                    !inpl3d(n,10)    
      k3 = kinc_coarse                             !inpl3d(n,11)    
      if((kinc_coarse.gt.1).and.(mod(kdim-1,kinc_coarse).ne.0))
     .  k3=kinc_coarse/2                                            
c
      iflag = 1
c
      if (0.le.0) then              ! use grid pt option for coarse movie
      !if (inpl3d(n,2).le.0) then   
c
c        grid point data
c
         if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
c
            call lead(nbl,lw,lw2,maxbl)
c
            jdw = (j2-j1)/j3 + 1
            kdw = (k2-k1)/k3 + 1
            idw = (i2-i1)/i3 + 1
c
c           check storage availability
c
            ixwk = 1
            nset = 5
            ibwk = ixwk + jdw*kdw*idw*nset
            ixgk = ibwk + jdim*kdim*idim*2
            ixvk = ixgk + jdw*kdw*idw*4
            nroom = nwork - (ixvk+jdw*kdw*idw*5-1)
            if (nroom .lt. 0.) then
               if (myid.eq.myhost) then
                  write(11,'('' not enough memory for plot3d_coarse'')')
                  write(11,'('' have, need = '',2i12)') nwork,
     .            nwork-nroom
                  write(11,'('' Aborting -not calling plot3d_coarse'')')
#if defined DIST_MPI
                  call MPI_ABORT(MPI_COMM_WORLD, myid, mpierror)
#else
                  stop
#endif
               end if
            else
               call plot3d_coarse(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,
     .                     k1,k2,k3,
     .                     w(lq),w(lqi0),w(lqj0),w(lqk0),w(lx),
     .                     w(ly),w(lz),wk(ixwk),wk(ibwk),
     .                     w(lblk),wk(ixgk),iflag,w(lvis),iovrlp(nbl),
     .                     nbl,nmap,w(lbcj),w(lbck),w(lbci),
     .                     w(lvj0),w(lvk0),w(lvi0),ifunc,n,jdw,kdw,idw,
     .                     nplots,jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                     jeg,keg,ieg,ninter,iindex,intmax,nsub1,
     .                     maxxe,nblk,nbli,limblk,isva,nblon,mxbli,
     .                     thetay,maxbl,maxgr,myid,myhost,mycomm,
     .                     mblk2nd,inpl3d,nblock,nblkpt,wk(ixvk),
     .                     w(lsj),w(lsk),w(lsi),w(lvol),nset)
            end if
c
         end if
c
      else if (0.eq.1 .or. 0.gt.2) then    ! never use cell center here
      !else if (inpl3d(n,2).eq.1 .or. inpl3d(n,2).gt.2) then          
c
c        cell center or face center data
c
         if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
c
            call lead(nbl,lw,lw2,maxbl)
c
            i2 = min(idim-1,i2)
            j2 = min(jdim-1,j2)
            k2 = min(kdim-1,k2)
            i1 = min(idim-1,i1)
            j1 = min(jdim-1,j1)
            k1 = min(kdim-1,k1)
            jdw = (j2-j1)/j3 + 1
            kdw = (k2-k1)/k3 + 1
            idw = (i2-i1)/i3 + 1
c
c           check storage availability
c
            ixwk = 1
            ibwk = ixwk + jdw*kdw*idw*5
            ixgk = ibwk + jdim*kdim*idim
            nroom = nwork - (ixgk+jdw*kdw*idw*4-1)
            if (nroom .lt. 0.) then
               if (myid.eq.myhost) then
                  write(11,'('' not enough memory for plot3c'')')
                  write(11,'('' have, need = '',2i12)') nwork,
     .            nwork-nroom
                  write(11,'('' not writing out plot3d files'')')
               end if
            else
               call plot3c(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,
     .                     w(lq),w(lqi0),w(lqj0),w(lqk0),w(lx),w(ly),
     .                     w(lz),wk(ixwk),wk(ibwk),w(lblk),wk(ixgk),
     .                     iflag,w(lvis),w(lvi0),w(lvj0),w(lvk0),
     .                     iovrlp(nbl),nbl,nmap,w(lsnk0),ifunc,n,
     .                     jdw,kdw,idw,
     .                     nplots,jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                     jeg,keg,ieg,ninter,iindex,intmax,nsub1,
     .                     maxxe,nblk,nbli,limblk,isva,nblon,mxbli,
     .                     thetay,maxbl,maxgr,myid,myhost,mycomm,
     .                     mblk2nd,inpl3d,nblock,nblkpt,ip3dsurf)
            end if
         end if
c
      elseif(0.eq.2) then    ! never write turb vars on cc
      !else                  
c
c        cell center turbulence data (plot3d q file format)
c
         if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
c
            if (ivmx .gt. 1) then
c
               call lead(nbl,lw,lw2,maxbl)
c
               i2 = min(idim-1,i2)
               j2 = min(jdim-1,j2)
               k2 = min(kdim-1,k2)
               i1 = min(idim-1,i1)
               j1 = min(jdim-1,j1)
               k1 = min(kdim-1,k1)
               jdw = (j2-j1)/j3 + 1
               kdw = (k2-k1)/k3 + 1
               idw = (i2-i1)/i3 + 1
c
c              check storage availability
c
               if (ivmx.eq.8 .or. ivmx.eq.9 .or. ivmx.ge.11) then
c                 have permanent storage for ux, starting at lux
                  ixwk  = 1
                  ibwk  = ixwk + jdw*kdw*idw*5
                  ixgk  = ibwk + jdim*kdim*idim
                  ibwk3 = ixgk + jdw*kdw*idw*4
                  nroom=nwork-(ibwk3+jdim*kdim*9)
                  if (nroom .lt. 0.) then
                     if (myid.eq.myhost) then
                        write(11,'('' not enough memory for plot3t'')')
                        write(11,'('' have, need = '',2i12)') nwork,
     .                  nwork-nroom
                        write(11,'('' not writing out plot3d files'')')
                     end if
                  else
                     call plot3t(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,
     .                           k1,k2,k3,w(lq),w(lx),w(ly),w(lz),
     .                           wk(ixwk),wk(ibwk),w(lblk),wk(ixgk),
     .                           w(lvis),iovrlp(nbl),nbl,nmap,w(lsj),
     .                           w(lsk),w(lsi),w(lsnk0),w(lux),
     .                           w(lxib),w(lvol),w(lqj0),w(lqk0),
     .                           w(lqi0),w(lbcj),w(lbck),w(lbci),
     .                           wk(ibwk3),w(lcmuv),jdw,kdw,idw,nplots,
     .                           jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                           ieg,jeg,keg,ninter,iindex,intmax,nsub1,
     .                           maxxe,nblk,nbli,limblk,isva,nblon,
     .                           mxbli,thetay,maxbl,maxgr,myid,myhost,
     .                           mycomm,mblk2nd,inpl3d,nblock,nblkpt,
     .                           w(lvolj0),w(lvolk0),w(lvoli0),vormax,
     .                           ivmax,jvmax,kvmax,nummem,5)
                  end if
               else
c                 need to grab ux storage from temporary work array,
c                 starting at location ibwk2
                  ixwk  = 1
                  ibwk  = ixwk + jdw*kdw*idw*5
                  ixgk  = ibwk + jdim*kdim*idim
                  ibwk2 = ixgk + jdw*kdw*idw*4
                  ibwk3 = ibwk2 +(jdim-1)*(kdim-1)*(idim-1)*9
                  nroom=nwork-(ibwk3+jdim*kdim*9)
                  if (nroom .lt. 0.) then
                     if (myid.eq.myhost) then
                        write(11,'('' not enough memory for plot3t'')')
                        write(11,'('' have, need = '',2i12)') nwork,
     .                  nwork-nroom
                        write(11,'('' not writing out plot3d files'')')
                     end if
                  else
                     call plot3t(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,
     .                           k1,k2,k3,w(lq),w(lx),w(ly),w(lz),
     .                           wk(ixwk),wk(ibwk),w(lblk),wk(ixgk),
     .                           w(lvis),iovrlp(nbl),nbl,nmap,w(lsj),
     .                           w(lsk),w(lsi),w(lsnk0),wk(ibwk2),
     .                           w(lxib),w(lvol),w(lqj0),w(lqk0),
     .                           w(lqi0),w(lbcj),w(lbck),w(lbci),
     .                           wk(ibwk3),w(lcmuv),jdw,kdw,idw,nplots,
     .                           jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                           jeg,keg,ieg,ninter,iindex,intmax,nsub1,
     .                           maxxe,nblk,nbli,limblk,isva,nblon,
     .                           mxbli,thetay,maxbl,maxgr,myid,myhost,
     .                           mycomm,mblk2nd,inpl3d,nblock,nblkpt,
     .                           w(lvolj0),w(lvolk0),w(lvoli0),vormax,
     .                           ivmax,jvmax,kvmax,nummem,5)
                  end if
               end if
c
            end if
c
         end if
c
      end if
c
      if (iblnkfr .eq. 0) then
c
c        reset blank values at fringe points to 0.
c
         if (mblk2nd(nbl).eq.myid .and. iovrlp(nbl).ne.0) then
            call blnkfr(nbl,iibg,kkbg,jjbg,ibpntsg,lbg,iitot,w(lblk),
     .                  jdim,kdim,idim,maxbl,0.)
         end if
      end if
c
   70 continue
c
      if (myid.eq.myhost) then
         if (lhdr .gt. 0) write(11,1096)
      endif
c
c     for stationary grid cases, set icallcrs flag to prevent output
c     of grid to plot3d file on subsequent calls (which occur only if 
c     abs(movie) > 0).  For dynamic grid cases, grid is output every
c     time plot3d routine is is called.
c
      if (iunst.eq.0) then
         icallcrs = 1
      else
         icallcrs = 0
      end if
c
  231 continue
c
c***********************************************************************
c     Print solution data.
c***********************************************************************
c
      if (0.gt.0) then       ! never write to print file here
      !if (nprint.gt.0) then 
c
      if (myid.eq.myhost) then
         if (lhdr .gt. 0 .and. np3d .le. 0) write(11,1096)            
         !if (lhdr .gt. 0 .and. nplot3d .le. 0) write(11,1096)        
      end if
c
      if (myid .eq. myhost) then
         write(17,111)(real(title(i)),i=1,20)
  111    format(/,20a4)
         write(17,21)
   21    format(6x,4hMach,5x,5halpha,6x,4hbeta,6x,4hReUe,3x,7hTinf,dR,
     .   6x,4htime)
         alphaw=radtodeg*alpha
         betaw =radtodeg*beta
         write(17,20)real(xmach),real(alphaw),real(betaw),real(reue),
     .               real(tinf),real(time)
   20    format(3f10.5,e10.3,2f10.5)
      end if
c
      do 80 n=1,nprint
c
      nbl = inpr(n,1)
c
      if (nbl.gt.nblock) then
         if (myid .eq. myhost) then
            write(11,77)nbl
         end if
   77    format(6h Block,i3,36h does not exist.  No output printed.)
         go to 80
      end if
c
      i1 = inpr(n,3)
      i2 = inpr(n,4)
      i3 = inpr(n,5)
      j1 = inpr(n,6)
      j2 = inpr(n,7)
      j3 = inpr(n,8)
      k1 = inpr(n,9)
      k2 = inpr(n,10)
      k3 = inpr(n,11)
c  
      iflag = 2
c
      if (inpr(n,2).eq.0) then
c
c        grid point data
c
         if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
c
            call lead(nbl,lw,lw2,maxbl)
c
            if (myid.eq.myhost) then
               write(17,348) nbl,igridg(nbl),idim,jdim,kdim
  348          format(//1x,5hBLOCK,i4,2x,5h(GRID,i4,1h),5x,
     .               15hIDIM,JDIM,KDIM=,3I5)
               write(17,349)
  349          format(1x,32hNOTE: endpts may not be reliable)
            end if
c
            jdw = (j2-j1)/j3 + 1
            kdw = (k2-k1)/k3 + 1
            idw = (i2-i1)/i3 + 1
c
c           check storage availability
c
            ixwk = 1
            nset = 8
            ibwk = ixwk + jdw*kdw*idw*nset
            ixgk = ibwk + jdim*kdim*idim*2
            ixvk = ixgk + jdw*kdw*idw*4
            nroom = nwork - (ixvk+jdw*kdw*idw*5-1)
            if (nroom .lt. 0.) then
               if (myid.eq.myhost) then
                  write(11,'('' not enough memory for plot3d'')')
                  write(11,'('' have, need = '',2i12)') nwork,
     .            nwork-nroom
                  write(11,'('' not writing out prout file'')')
               end if
            else
               call plot3d_coarse(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,
     .                     k1,k2,k3,
     .                     w(lq),w(lqi0),w(lqj0),w(lqk0),w(lx),
     .                     w(ly),w(lz),wk(ixwk),wk(ibwk),
     .                     w(lblk),wk(ixgk),iflag,w(lvis),iovrlp(nbl),
     .                     nbl,nmap,w(lbcj),w(lbck),w(lbci),
     .                     w(lvj0),w(lvk0),w(lvi0),ifunc,n,jdw,kdw,idw,
     .                     nplots,jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                     jeg,keg,ieg,ninter,iindex,intmax,nsub1,
     .                     maxxe,nblk,nbli,limblk,isva,nblon,mxbli,
     .                     thetay,maxbl,maxgr,myid,myhost,mycomm,
     .                     mblk2nd,inpl3d,nblock,nblkpt,wk(ixvk),
     .                     w(lsj),w(lsk),w(lsi),w(lvol),nset)
            end if
c
         end if
c
      else
c
c        cell center data
c
         if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
c
            call lead(nbl,lw,lw2,maxbl)
c
            if (myid.eq.myhost) then
               write(17,348)nbl,igridg(nbl),idim,jdim,kdim
            end if
c
            i2 = min(idim-1,i2)
            j2 = min(jdim-1,j2)
            k2 = min(kdim-1,k2)
            i1 = min(idim-1,i1)
            j1 = min(jdim-1,j1)
            k1 = min(kdim-1,k1)
            jdw = (j2-j1)/j3 + 1
            kdw = (k2-k1)/k3 + 1
            idw = (i2-i1)/i3 + 1
c
c           check storage availability
c
            ixwk = 1
            ibwk = ixwk + jdw*kdw*idw*5
            ixgk = ibwk + jdim*kdim*idim
            nroom = nwork - (ixgk+jdw*kdw*idw*4-1)
            if (nroom .lt. 0.) then
               if (myid.eq.myhost) then
                  write(11,'('' not enough memory for plot3c'')')
                  write(11,'('' have, need = '',2i12)') nwork,
     .            nwork-nroom
                  write(11,'('' not writing out plot3d files'')')
               end if
            else
               call plot3c(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,
     .                     w(lq),w(lqi0),w(lqj0),w(lqk0),w(lx),w(ly),
     .                     w(lz),wk(ixwk),wk(ibwk),w(lblk),wk(ixgk),
     .                     iflag,w(lvis),w(lvi0),w(lvj0),w(lvk0),
     .                     iovrlp(nbl),nbl,nmap,w(lsnk0),ifunc,n,
     .                     jdw,kdw,idw,
     .                     nplots,jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                     jeg,keg,ieg,ninter,iindex,intmax,nsub1,
     .                     maxxe,nblk,nbli,limblk,isva,nblon,mxbli,
     .                     thetay,maxbl,maxgr,myid,myhost,mycomm,
     .                     mblk2nd,inpl3d,nblock,nblkpt,iprnsurf)
            end if
         end if
c
      end if
c
   80 continue
c
      if (myid.eq.myhost) then
         if (lhdr .gt. 0) write(11,1096)
      end if
c
      end if
      return
      end
