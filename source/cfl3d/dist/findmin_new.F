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
      subroutine findmin_new(lw,lw2,w,mgwk,wk,nwork,iwk,iwork,nsurf,j1,
     .                       xs,ixs,lsminn,ireq_xs,ireq_bb,ngrid,
     .                       ncgg,nbci0,nbcj0,nbck0,
     .                       nbcidim,nbcjdim,nbckdim,jbcinfo,
     .                       kbcinfo,ibcinfo,nblg,nou,bou,nbuf,ibufdim,
     .                       maxbl,maxgr,maxseg,mblk2nd)
c
c     $Id$
c
c***********************************************************************
c     Purpose: To serve as a driver routine for computing distance to
c              closest solid surface (note that search region involves
c              information across all the blocks for all solid segments)
c
c              This routine is much (MUCH!) faster than the old findmin
c              routine for large grids ( > 10e6 points). For small or 
c              medium grids, it is slightly slower, although the result
c              in all cases is believed to be more accurate, especially 
c              where the grid is skewed.
c
c              Completely new version of the tlns3d routine disdr 
c              introduced by Larry Wigton (Boeing) May 1995. Uses
c              routines written by Mac Ice (Boeing) based on a 
c              recursive-box algorithm. This speeds up nearest 
c              point calculation.
c
c              This routine assumes the minumum distance at every point
c              has previously been initialized to a large value, and 
c              that the number of viscous surface points (nsurf) has
c              been counted
c
c     coding history: L. Wigton.....original C version for tlns3d
c                     E. Parlette...converted tlns3d C version to fortran
c                     R. Biedron....converted fortran version for cfl3d
c                                   added mpi calls
c                                   added mods for Baldwin-Barth model
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
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension w(mgwk),lw(65,maxbl),lw2(43,maxbl)
      dimension wk(nwork),xs(4*nsurf)
      dimension iwk(iwork),ixs(4*nsurf)
      dimension lsminn(maxbl)
      dimension ireq_xs(maxbl*maxseg*6),ireq_bb(maxbl*maxseg*6)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension ncgg(maxgr),nblg(maxgr),mblk2nd(maxbl)
c
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /maxiv/ ivmx
      common /mydist2/ nnodes,myhost,myid,mycomm
c
      integer vlist,bbdef,surf
      integer surfi
c
c
c***************************************************************************
c     memory allocation based on nsurf
c***************************************************************************
c
      call initi(iwork)
      call initf(nwork)
      surf  = ifalloc(nou,bou,nbuf,ibufdim,myid,4*nsurf)
      ntri  = iialloc(nou,bou,nbuf,ibufdim,myid,nsurf)
      iptri = iialloc(nou,bou,nbuf,ibufdim,myid,8*nsurf)
      isurf = 0
      surfi = 0 
      if (ivmx.eq.4 .or. ivmx.eq.25) 
     .  surfi = iialloc(nou,bou,nbuf,ibufdim,myid,4*nsurf)
c
c***************************************************************************
c     collect all the viscous surface points into an array xsurf(nsurf,3)
c***************************************************************************
c
c     each interior surface point has 8 neighbors, which help define
c     surface triangles:
c
c            p1 --- p2 --- p3
c            |    / |  \   |
c            |   /  |   \  |
c            |  /   |    \ |
c            p4 --- sp --- p5
c            |  \   |    / |
c            |   \  |   /  |
c            |    \ |  /   |
c            p6 --- p7 --- p8
c
c     Points on a boundary will have fewer neighbors (as few as three
c     for a corner point.)
c     ntri(1,i) == the number of neighbor triangles for point i.
c     ntri(2,i) == the index in array iptri at which these neighbor
c                  points are defined.
c     iptri(ntri(2,i))   == p4
c     iptri(ntri(2,i)+1) == p7
c     iptri(ntri(2,i)+2) == p5
c     iptri(ntri(2,i)+3) == p7
c     iptri(ntri(2,i)+4) == p5
c     iptri(ntri(2,i)+5) == p2
c     iptri(ntri(2,i)+6) == p4
c     iptri(ntri(2,i)+7) == p2
c
#if defined DIST_MPI
      nreq    = 1
      nreq_bb = 1
c
c     set baseline tag values
c
      ioffset  = maxgr*maxseg*6
      itag     = 0
      itag_xs  = 1
      itag_ixs = itag_xs + ioffset
#endif
c
c     assemble surface points from all blocks into xs array, and
c     for Baldwin-Barth mode, assemble surface indicies from all
c     blocks into ixs array. the xs and ixs arrays are sent/recieved
c     to from all blocks.
c
      do igrid = 1,ngrid
         nbl  = nblg(igrid)
         call lead(nbl,lw,lw2,maxbl)
c
c        i=constant surfaces
c
         do m = 1,2
         if (m.eq.1) then
            ns    = nbci0(nbl)
            nface = 1
         else
            ns    = nbcidim(nbl)
            nface = 2
         end if
            do iseg = 1,ns
               nbctype  =  ibcinfo(nbl,iseg,1,m)
               n1beg    =  ibcinfo(nbl,iseg,2,m)
               n1end    =  ibcinfo(nbl,iseg,3,m)
               n2beg    =  ibcinfo(nbl,iseg,4,m)
               n2end    =  ibcinfo(nbl,iseg,5,m)
               npts     =  (n1end-n1beg+1)*(n2end-n2beg+1)
               nvalxs   =  4*npts
               nvalbb   =  4*npts
c
               if(abs(nbctype).eq.2004 .or. abs(nbctype).eq.2014 .or.
     .            abs(nbctype).eq.2024 .or. abs(nbctype).eq.2034 .or.
     .            abs(nbctype).eq.2016) then
c
                  if (myid .eq. mblk2nd(nbl)) then
                     if (ivmx .eq.4 .or. ivmx.eq.25) then
                        call getptsbb(1,1,1,idim,jdim,kdim,idim,
     .                  jdim,kdim,w(lx),w(ly),w(lz),nface,n1beg,n1end,
     .                  n2beg,n2end,xs,ixs,nsurf,nbl,
     .                  w(lblk),w(lbci),w(lbcj),w(lbck))
                     else
                        call getpts(1,1,1,idim,jdim,kdim,idim,
     .                  jdim,kdim,w(lx),w(ly),w(lz),nface,n1beg,n1end,
     .                  n2beg,n2end,xs,nsurf,nbl,
     .                  w(lblk),w(lbci),w(lbcj),w(lbck))
                     end if
                  end if
c
#if defined DIST_MPI
c
                  itag = itag + 1
                  if (myid .ne. mblk2nd(nbl)) then

c                    recieve the data on non-local nodes
c
                     nd_srce = mblk2nd(nbl)
                     mytag   = itag_xs + itag
                     call MPI_Irecv (xs, nvalxs,
     .                               MY_MPI_REAL,
     .                               nd_srce,mytag, mycomm, 
     .                               ireq_xs(nreq), ierr)
                     nreq = nreq + 1
c
                     if (ivmx .eq.4 .or. ivmx.eq.25) then
                        nd_srce = mblk2nd(nbl)
                        mytag   = itag_ixs + itag
                        call MPI_Irecv(ixs, nvalbb, 
     .                                 MPI_INTEGER,
     .                                 nd_srce,mytag,mycomm,
     .                                 ireq_bb(nreq_bb),ierr)
                        nreq_bb = nreq_bb + 1
                     end if
c
                  else 
c
c                    send the data to non-local nodes
c
                     do inode = 1,nnodes
                        if (inode .ne. myid) then
                           nd_dest = inode
                           mytag   = itag_xs + itag
                           call MPI_Send (xs, nvalxs,
     .                                    MY_MPI_REAL,
     .                                    nd_dest, mytag, 
     .                                    mycomm, ierr)
                           if (ivmx .eq.4 .or. ivmx.eq.25) then
                              nd_dest = inode
                              mytag   = itag_ixs + itag
                              call MPI_Send (ixs, nvalbb,
     .                                      MPI_INTEGER,
     .                                      nd_dest, mytag,
     .                                      mycomm, ierr)
                           end if
                        end if
                     end do
c
                  end if
c
                  if (nreq-1 .gt. 0) then
                     call MPI_Wait(ireq_xs(nreq-1), istat, ierr)
                  end if
                  if (ivmx.eq.4 .or. ivmx.eq.25) then
                     if (nreq_bb-1 .gt. 0) then
                        call MPI_Wait(ireq_bb(nreq_bb-1), istat, ierr)
                     end if
                  end if
c
#endif
c
                  if (ivmx .eq.4 .or. ivmx.eq.25) then
                     call collect_surfbb(1,1,1,idim,jdim,kdim,idim,
     .               jdim,kdim,xs,ixs,nface,n1beg,
     .               n1end,n2beg,n2end,nsurf,wk(surf),isurf,iwk(ntri),
     .               iwk(iptri),iwk(surfi),nbl)
                  else
                     call collect_surf(1,1,1,idim,jdim,kdim,idim,
     .               jdim,kdim,xs,nface,n1beg,n1end,
     .               n2beg,n2end,nsurf,wk(surf),isurf,iwk(ntri),
     .               iwk(iptri))
                  end if
c
                end if
            end do
         end do
c
c        j=constant surfaces
c
         do m = 1,2
         if (m.eq.1) then
            ns    = nbcj0(nbl)
            nface = 3
         else
            ns    = nbcjdim(nbl)
            nface = 4
         end if
            do iseg = 1,ns
               nbctype  =  jbcinfo(nbl,iseg,1,m)
               n1beg    =  jbcinfo(nbl,iseg,4,m)
               n1end    =  jbcinfo(nbl,iseg,5,m)
               n2beg    =  jbcinfo(nbl,iseg,2,m)
               n2end    =  jbcinfo(nbl,iseg,3,m)
               npts     =  (n1end-n1beg+1)*(n2end-n2beg+1)
               nvalxs   =  4*npts
               nvalbb   =  4*npts
c
               if(abs(nbctype).eq.2004 .or. abs(nbctype).eq.2014 .or.
     .            abs(nbctype).eq.2024 .or. abs(nbctype).eq.2034 .or.
     .            abs(nbctype).eq.2016) then
c
                  if (myid .eq. mblk2nd(nbl)) then
                     if (ivmx .eq.4 .or. ivmx.eq.25) then
                        call getptsbb(1,1,1,idim,jdim,kdim,idim,
     .                  jdim,kdim,w(lx),w(ly),w(lz),nface,n1beg,n1end,
     .                  n2beg,n2end,xs,ixs,nsurf,nbl,
     .                  w(lblk),w(lbci),w(lbcj),w(lbck))
                     else
                        call getpts(1,1,1,idim,jdim,kdim,idim,
     .                  jdim,kdim,w(lx),w(ly),w(lz),nface,n1beg,n1end,
     .                  n2beg,n2end,xs,nsurf,nbl,
     .                  w(lblk),w(lbci),w(lbcj),w(lbck))
                     end if
                  end if
c
#if defined DIST_MPI
c
                  itag = itag + 1
                  if (myid .ne. mblk2nd(nbl)) then

c                    recieve the data on non-local nodes
c
                     nd_srce = mblk2nd(nbl)
                     mytag   = itag_xs + itag
                     call MPI_Irecv (xs, nvalxs,
     .                               MY_MPI_REAL,
     .                               nd_srce,mytag, mycomm, 
     .                               ireq_xs(nreq), ierr)
                     nreq = nreq + 1
c
                     if (ivmx .eq.4 .or. ivmx.eq.25) then
                        nd_srce = mblk2nd(nbl)
                        mytag   = itag_ixs + itag
                        call MPI_Irecv(ixs, nvalbb, 
     .                                 MPI_INTEGER,
     .                                 nd_srce,mytag,mycomm,
     .                                 ireq_bb(nreq_bb),ierr)
                        nreq_bb = nreq_bb + 1
                     end if
c
                  else 
c
c                    send the data to non-local nodes
c
                     do inode = 1,nnodes
                        if (inode .ne. myid) then
                           nd_dest = inode
                           mytag   = itag_xs + itag
                           call MPI_Send (xs, nvalxs,
     .                                    MY_MPI_REAL,
     .                                    nd_dest, mytag, 
     .                                    mycomm, ierr)
                           if (ivmx .eq.4 .or. ivmx.eq.25) then
                              nd_dest = inode
                              mytag   = itag_ixs + itag
                              call MPI_Send (ixs, nvalbb,
     .                                      MPI_INTEGER,
     .                                      nd_dest, mytag,
     .                                      mycomm, ierr)
                           end if
                        end if
                     end do
c
                  end if
c
                  if (nreq-1 .gt. 0) then
                     call MPI_Wait(ireq_xs(nreq-1), istat, ierr)
                  end if
                  if (ivmx.eq.4 .or. ivmx.eq.25) then
                     if (nreq_bb-1 .gt. 0) then
                        call MPI_Wait(ireq_bb(nreq_bb-1), istat, ierr)
                     end if
                  end if
c
#endif
c
                  if (ivmx .eq.4 .or. ivmx.eq.25) then
                     call collect_surfbb(1,1,1,idim,jdim,kdim,idim,
     .               jdim,kdim,xs,ixs,nface,n1beg,
     .               n1end,n2beg,n2end,nsurf,wk(surf),isurf,iwk(ntri),
     .               iwk(iptri),iwk(surfi),nbl)
                  else
                     call collect_surf(1,1,1,idim,jdim,kdim,idim,
     .               jdim,kdim,xs,nface,n1beg,n1end,
     .               n2beg,n2end,nsurf,wk(surf),isurf,iwk(ntri),
     .               iwk(iptri))
                  end if
c
                end if
            end do
         end do
c
c        k=constant surfaces
c
         do m = 1,2
         if (m.eq.1) then
            ns    = nbck0(nbl)
            nface = 5
         else
            ns    = nbckdim(nbl)
            nface = 6
         end if
            do iseg = 1,ns
               nbctype  =  kbcinfo(nbl,iseg,1,m)
               n1beg    =  kbcinfo(nbl,iseg,2,m)
               n1end    =  kbcinfo(nbl,iseg,3,m)
               n2beg    =  kbcinfo(nbl,iseg,4,m)
               n2end    =  kbcinfo(nbl,iseg,5,m)
               npts     =  (n1end-n1beg+1)*(n2end-n2beg+1)
               nvalxs   =  4*npts
               nvalbb   =  4*npts
c
               if(abs(nbctype).eq.2004 .or. abs(nbctype).eq.2014 .or.
     .            abs(nbctype).eq.2024 .or. abs(nbctype).eq.2034 .or.
     .            abs(nbctype).eq.2016) then
c
                  if (myid .eq. mblk2nd(nbl)) then
                     if (ivmx .eq.4 .or. ivmx.eq.25) then
                        call getptsbb(1,1,1,idim,jdim,kdim,idim,
     .                  jdim,kdim,w(lx),w(ly),w(lz),nface,n1beg,n1end,
     .                  n2beg,n2end,xs,ixs,nsurf,nbl,
     .                  w(lblk),w(lbci),w(lbcj),w(lbck))
                     else
                        call getpts(1,1,1,idim,jdim,kdim,idim,
     .                  jdim,kdim,w(lx),w(ly),w(lz),nface,n1beg,n1end,
     .                  n2beg,n2end,xs,nsurf,nbl,
     .                  w(lblk),w(lbci),w(lbcj),w(lbck))
                     end if
                  end if
c
#if defined DIST_MPI
c
                  itag = itag + 1
                  if (myid .ne. mblk2nd(nbl)) then

c                    recieve the data on non-local nodes
c
                     nd_srce = mblk2nd(nbl)
                     mytag   = itag_xs + itag
                     call MPI_Irecv (xs, nvalxs,
     .                               MY_MPI_REAL,
     .                               nd_srce,mytag, mycomm, 
     .                               ireq_xs(nreq), ierr)
                     nreq = nreq + 1
c
                     if (ivmx .eq.4 .or. ivmx.eq.25) then
                        nd_srce = mblk2nd(nbl)
                        mytag   = itag_ixs + itag
                        call MPI_Irecv(ixs, nvalbb, 
     .                                 MPI_INTEGER,
     .                                 nd_srce,mytag,mycomm,
     .                                 ireq_bb(nreq_bb),ierr)
                        nreq_bb = nreq_bb + 1
                     end if
c
                  else
c
c                    send the data to non-local nodes
c
                     do inode = 1,nnodes
                        if (inode .ne. myid) then
                           nd_dest = inode
                           mytag   = itag_xs + itag
                           call MPI_Send (xs, nvalxs,
     .                                    MY_MPI_REAL,
     .                                    nd_dest, mytag, 
     .                                    mycomm, ierr)
                           if (ivmx .eq.4 .or. ivmx.eq.25) then
                              nd_dest = inode
                              mytag   = itag_ixs + itag
                              call MPI_Send (ixs, nvalbb,
     .                                      MPI_INTEGER,
     .                                      nd_dest, mytag,
     .                                      mycomm, ierr)
                           end if
                        end if
                     end do
c
                  end if
c
                  if (nreq-1 .gt. 0) then
                     call MPI_Wait(ireq_xs(nreq-1), istat, ierr)
                  end if
                  if (ivmx.eq.4 .or. ivmx.eq.25) then
                     if (nreq_bb-1 .gt. 0) then
                        call MPI_Wait(ireq_bb(nreq_bb-1), istat, ierr)
                     end if
                  end if
c
#endif
c
                  if (ivmx .eq.4 .or. ivmx.eq.25) then
                     call collect_surfbb(1,1,1,idim,jdim,kdim,idim,
     .               jdim,kdim,xs,ixs,nface,n1beg,
     .               n1end,n2beg,n2end,nsurf,wk(surf),isurf,iwk(ntri),
     .               iwk(iptri),iwk(surfi),nbl)
                  else
                     call collect_surf(1,1,1,idim,jdim,kdim,idim,
     .               jdim,kdim,xs,nface,n1beg,n1end,
     .               n2beg,n2end,nsurf,wk(surf),isurf,iwk(ntri),
     .               iwk(iptri))
                  end if
c
               end if
            end do
         end do
c
      end do
c
#if defined DIST_MPI
c
c     check tags
c
      if (itag_xs + itag .gt. itag_ixs) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)'not enough tags between itag_xs ',
     .                 'and itag_ixs'
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)'itag,itag_xs,itag_ixs = ',
     .                  itag,itag_xs,itag_ixs
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
#endif
c
c***********************************************************************
c     determine minimum distance
c***********************************************************************
c
c     to help locate identical points on boundaries, we sort with
c     respect to x coordinate.  Better to use a stable sort on
c     x,y,z and ntri consecutively, but this is good enough
c     for time being.
c
      if (ivmx .eq. 4 .or. ivmx.eq.25) then
         call sort_xbb(nsurf,wk(surf),iwk(ntri),iwk(iptri),iwk(surfi),
     .                 wk,iwk,nou,bou,nbuf,ibufdim,myid)
      else
         call sort_x(nsurf,wk(surf),iwk(ntri),iwk(iptri),wk,iwk,
     .               nou,bou,nbuf,ibufdim,myid)
      end if
c
c     run Mac Ice routines for making boxes that surround
c     groups of viscous surface points
c
      minbox=sqrt(float(nsurf))
      minbox=max(minbox,50)
      nbb = 3*nsurf/minbox
c
c     bbdef(1..6,i) == xmin,xmax,ymin,ymax,zmin,zmax of i-th
c                      bounding box.
c     ipv(1,i) == number of points in i-th bounding box.
c     ipv(2,i) == location in vlist of points in i-th bounding box
c     vlist(ipv(2,i),...,vlist(ipv(2,i)+ipv(1,i)-1) ==
c           location in array surf of points in i-th bounding box
c
      bbdef = ifalloc(nou,bou,nbuf,ibufdim,myid,6*nbb)
      ipv = iialloc(nou,bou,nbuf,ibufdim,myid,2*nbb)
      vlist = iialloc(nou,bou,nbuf,ibufdim,myid,nsurf)
      call makebb(minbox,10,nsurf,wk(surf),wk(surf+nsurf),
     . wk(surf+2*nsurf),nbb,ntotv,wk(bbdef),iwk(ipv),iwk(vlist),iwk,
     . nou,bou,nbuf,ibufdim,myid)
c
c    calculate distance from each point in the field to points on
c    the viscous surfaces and then look at neighboring triangles
c
c     The nearest point on any surface to a field point is found
c     as follows:
c
c     1. The distance from the field point to each bounding
c        box is calculated.
c     2. The nearest bounding box is found, and the the nearest
c        point in that box is found (by exhaustive search.)
c     3. Go to the next nearest bounding box. If it is farther
c        away than the nearest point found so far, then we have
c        the nearest point.
c     4. Find the nearest point in this box, and compare with the
c        previous nearest point. Return to step 3.
c
c     because the permanent storage space for for the minimum 
c     distance array in cfl3d is sized for cell centers, and the
c     distance calculated in calc_dist is to grid points, the
c     returned minimum distance function from calc_dist is
c     stored in the work array wk, starting at lsmin.
c
      do igrid=1,ngrid
c
c        fine level
c
         nbl  = nblg(igrid)+j1
         call lead(nbl,lw,lw2,maxbl)
         if (myid .eq.  mblk2nd(nbl)) then
            lsmin = ifalloc(nou,bou,nbuf,ibufdim,myid,jdim*kdim*idim)
            lsminn(nbl) = lsmin
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),993) nbl
            if (ivmx .eq. 4 .or. ivmx.eq.25) then
               call calc_distbb(1,1,1,idim,jdim,kdim,idim,jdim,kdim,
     .         idim,jdim,kdim,w(lx),w(ly),w(lz),wk(lsmin),
     .         nsurf,wk(surf),nbb,wk(bbdef),iwk(ipv),
     .         iwk(vlist),iwk(ntri),iwk(iptri),iwk(surfi),w(lsni0),
     .         w(lxkb),w(lxib),w(lnbl),wk,iwk,j1,
     .         nou,bou,nbuf,ibufdim,myid)
            else
               call calc_dist(1,1,1,idim,jdim,kdim,idim,jdim,kdim,
     .         idim,jdim,kdim,w(lx),w(ly),w(lz),wk(lsmin),
     .         nsurf,wk(surf),nbb,wk(bbdef),iwk(ipv),
     .         iwk(vlist),iwk(ntri),iwk(iptri),wk,iwk,
     .         nou,bou,nbuf,ibufdim,myid)
            end if
         end if
c
c        coarser levels 
c
         ncg = ncgg(igrid)-j1
         if (ncg.gt.0) then
            do 105 m=1,ncg
            nbl    = nbl+1
            if (myid .eq.  mblk2nd(nbl)) then
               lsminc = ifalloc(nou,bou,nbuf,ibufdim,myid,jj2*kk2*ii2)
               lsminn(nbl) = lsminc
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),993) nbl
               call distcg(jdim,kdim,idim,wk(lsmin),
     .                     jj2,kk2,ii2,wk(lsminc))
               call lead(nbl,lw,lw2,maxbl)
               lsmin = lsminc
            end if
 105        continue
         end if
c
         if (ncg.gt.0 .and. (ivmx.eq.4 .or. ivmx.eq.25)) then
            nblc  = nblg(igrid)+j1
            do 1049 m=1,ncg
            call lead(nblc,lw,lw2,maxbl)
            nbl    = nblc
            nblc   = nblc+1
            if (myid .eq.  mblk2nd(nbl)) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)'getting bbarth data on ',
     .         'coarser block ',nblc
               jdimc = jdim/2+1
               kdimc = kdim/2+1
               idimc = idim/2+1
               call bbarthcg(jdim,kdim,idim,jdimc,kdimc,idimc,
     .         w(lw(15,nbl)),w(lw(26,nbl)),w(lw(19,nbl)),w(lw(27,nbl)),
     .         w(lw(15,nblc)),w(lw(26,nblc)),w(lw(19,nblc)),
     .         w(lw(27,nblc)))
            end if
 1049       continue
         end if
c
      end do
c
c     generate distance to cell-centers and store in permanent array w
c
      do 1050 igrid = 1,ngrid
      nbl  = nblg(igrid)+j1
      if (myid .eq.  mblk2nd(nbl)) then
         call lead(nbl,lw,lw2,maxbl)
         lsmin = lsminn(nbl)
         call distcc(jdim,kdim,idim,wk(lsmin),w(lsnk0))
         ncg = ncgg(igrid)-j1
         if (ncg.gt.0) then
            do 106 m=1,ncg
            nbl    = nbl+1
            call lead(nbl,lw,lw2,maxbl)
            lsmin = lsminn(nbl)
            call distcc(jdim,kdim,idim,wk(lsmin),w(lsnk0))
  106       continue
         end if
      end if
 1050 continue
c
  993 format(45h computing smin (min dist func) for field eqn,
     .       21h turb model for block,i6)
c
      return
      end
      subroutine getpts (imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                   x,y,z,nface,n1beg,n1end,n2beg,n2end,
     .                   xs,nsurf,nbl,blank,bci,bcj,bck)
c***********************************************************************
c     Purpose:  this subroutine collects the solid surface points from
c     a given segment on a block into the array xs. the points are put 
c     in the array as triplets (x,y,z) so that the addressing is
c     xs(3,npts).
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(jmp1,kmp1,imp1),y(jmp1,kmp1,imp1),z(jmp1,kmp1,imp1),
     .          xs(4,nsurf)
      dimension bcj(kmp1,imp1-1,2),bck(jmp1,imp1-1,2),
     .          bci(jmp1,kmp1,2),blank(jmp1,kmp1,imp1)
c
      ipts    = n1end - n1beg + 1
c
c slarge...large number to add at points that: 1) have a bctype=2004/14/24/34/16
c          but have an interface flag (bci/bcj/bck) of 0 rather than
c          the standard viscous wall value of 1. This occurs in certain
c          patching applications where part of a solid wall bc is
c          subsequently overwritten in the patching routine - those patched
c          cells must not be treated as lying on a solid surface; 2) surface
c          points that lie in a hole (chimera applications). adding the
c          large number to the distance to such points effectively excludes
c          them from becoming tagged as the minimum distance point to any
c          field point.
c
      slarge = 1.e30
      hole   = 1.
c
c---- branch based on block face number
c
c
      mm = 1
      if (nface .eq. 2 .or. nface .eq. 4 .or. nface .eq. 6) mm = 2
c
      go to (100,200,300,400,500,600) nface
c
  100 iw = imn
      iw1 = imn
      go to 201
  200 iw = imx
      iw1 = imx-1
  201 do kw=n2beg,n2end
        kk = min(kw,n2end-1)
        k1 = max(kw-1,n2beg)
        do jw=n1beg,n1end
          jj = min(jw,n1end-1)
          j1 = max(jw-1,n1beg)
          i = 1 +jw-n1beg +ipts*(kw-n2beg)
          bciuse1=ccmax(bci(jj,kk,mm),bci(j1,kk,mm))
          bciuse2=ccmax(bci(jj,k1,mm),bci(j1,k1,mm))
          bciuse =ccmax(bciuse1,bciuse2)
          blankuse1=ccmax(blank(jj,kk,iw1),blank(j1,kk,iw1))
          blankuse2=ccmax(blank(jj,k1,iw1),blank(j1,k1,iw1))
          blankuse =ccmax(blankuse1,blankuse2)
          factr = slarge*((1. - blankuse)*hole +
     .                    (1. -   bciuse))
          xs(1,i) = x(jw,kw,iw)
          xs(2,i) = y(jw,kw,iw)
          xs(3,i) = z(jw,kw,iw)
          xs(4,i) = factr
        enddo
      enddo
      return
c
  300 jw = jmn
      jw1 = jmn
      go to 401
  400 jw = jmx
      jw1 = jmx-1
  401 do iw=n2beg,n2end
        ii = min(iw,n2end-1)
        i1 = max(iw-1,n2beg)
        do kw=n1beg,n1end
          kk = min(kw,n1end-1)
          k1 = max(kw-1,n1beg)
          i = 1 +kw-n1beg +ipts*(iw-n2beg)
          bcjuse1=ccmax(bcj(kk,ii,mm),bcj(k1,ii,mm))
          bcjuse2=ccmax(bcj(kk,i1,mm),bcj(k1,i1,mm))
          bcjuse =ccmax(bcjuse1,bcjuse2)
          blankuse1=ccmax(blank(jw1,kk,ii),blank(jw1,k1,ii))
          blankuse2=ccmax(blank(jw1,kk,i1),blank(jw1,k1,i1))
          blankuse =ccmax(blankuse1,blankuse2)
          factr = slarge*((1. - blankuse)*hole +
     .                    (1. -   bcjuse))
          xs(1,i) = x(jw,kw,iw)
          xs(2,i) = y(jw,kw,iw)
          xs(3,i) = z(jw,kw,iw)
          xs(4,i) = factr
        enddo
      enddo
      return
c
  500 kw = kmn
      kw1 = kmn
      go to 601
  600 kw = kmx
      kw1 = kmx-1
  601 do jw=n2beg,n2end
        jj = min(jw,n2end-1) 
        j1 = max(jw-1,n2beg)
        do iw=n1beg,n1end
          ii = min(iw,n1end-1)
          i1 = max(iw-1,n1beg)
          i = 1 +iw-n1beg +ipts*(jw-n2beg)
          bckuse1=ccmax(bck(jj,ii,mm),bck(j1,ii,mm))
          bckuse2=ccmax(bck(jj,i1,mm),bck(j1,i1,mm))
          bckuse =ccmax(bckuse1,bckuse2)
          blankuse1=ccmax(blank(jj,kw1,ii),blank(j1,kw1,ii))
          blankuse2=ccmax(blank(jj,kw1,i1),blank(j1,kw1,i1))
          blankuse =ccmax(blankuse1,blankuse2)
          factr = slarge*((1. - blankuse)*hole +
     .                    (1. -   bckuse))
          xs(1,i) = x(jw,kw,iw)
          xs(2,i) = y(jw,kw,iw)
          xs(3,i) = z(jw,kw,iw)
          xs(4,i) = factr
        enddo
      enddo
      return
      end
      subroutine getptsbb(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                    x,y,z,nface,n1beg,n1end,n2beg,n2end,
     .                    xs,ixs,nsurf,nbl,blank,bci,bcj,bck)
c***********************************************************************
c     Purpose:  this subroutine collects the solid surface points from 
c     a given segment on a block into the array xs. also puts
c     index and block info into array ixs for use with Baldwin-Barth
c     turbulence model.
c
c     the points are put in the array as triplets (x,y,z) so
c     that the addressing is xs(3,npts); similarly, the index
c     and block infor are put in the array as i,j,k,nbl so
c     that the addressing is ixs(4,npts)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(jmp1,kmp1,imp1),y(jmp1,kmp1,imp1),z(jmp1,kmp1,imp1),
     .          xs(4,nsurf),ixs(4,nsurf)
      dimension bcj(kmp1,imp1-1,2),bck(jmp1,imp1-1,2),
     .          bci(jmp1,kmp1,2),blank(jmp1,kmp1,imp1)
c
      ipts    = n1end - n1beg + 1
      jpts    = n2end - n2beg + 1
c
c slarge...large number to add at points that: 1) have a bctype=2004/14/24/34/16
c          but have an interface flag (bci/bcj/bck) of 0 rather than
c          the standard viscous wall value of 1. This occurs in certain
c          patching applications where part of a solid wall bc is
c          subsequently overwritten in the patching routine - those patched
c          cells must not be treated as lying on a solid surface; 2) surface
c          points that lie in a hole (chimera applications). adding the
c          large number to the distance to such points effectively excludes
c          them from becoming tagged as the minimum distance point to any
c          field point.
c
      slarge = 1.e30
      hole   = 1.
c
c---- branch based on block face number
c
      mm = 1
      if (nface .eq. 2 .or. nface .eq. 4 .or. nface .eq. 6) mm = 2
c
      go to (100,200,300,400,500,600) nface
c
  100 iw = imn
      iw1 = imn
      go to 201
  200 iw = imx
      iw1 = imx - 1
  201 do kw=n2beg,n2end
        kk = min(kw,n2end-1)
        k1 = max(kw-1,n2beg)
        do jw=n1beg,n1end
          jj = min(jw,n1end-1)
          j1 = max(jw-1,n1beg)
          i = 1 +jw-n1beg +ipts*(kw-n2beg)
          bciuse1=ccmax(bci(jj,kk,mm),bci(j1,kk,mm))
          bciuse2=ccmax(bci(jj,k1,mm),bci(j1,k1,mm))
          bciuse =ccmax(bciuse1,bciuse2)
          blankuse1=ccmax(blank(jj,kk,iw1),blank(j1,kk,iw1))
          blankuse2=ccmax(blank(jj,k1,iw1),blank(j1,k1,iw1))
          blankuse =ccmax(blankuse1,blankuse2)
          factr = slarge*((1. - blankuse)*hole +
     .                    (1. -   bciuse))
          xs(1,i) = x(jw,kw,iw)
          xs(2,i) = y(jw,kw,iw)
          xs(3,i) = z(jw,kw,iw)
          xs(4,i) = factr
          ixs(1,i) = iw
          ixs(2,i) = jw
          ixs(3,i) = kw
          ixs(4,i) = nbl
        enddo
      enddo
      return
c
  300 jw = jmn
      jw1 = jmn
      go to 401
  400 jw = jmx
      jw1 = jmx - 1
  401 do iw=n2beg,n2end
        ii = min(iw,n2end-1)
        i1 = max(iw-1,n2beg)
        do kw=n1beg,n1end
          kk = min(kw,n1end-1)
          k1 = max(kw-1,n1beg)
          i = 1 +kw-n1beg +ipts*(iw-n2beg)
          bcjuse1=ccmax(bcj(kk,ii,mm),bcj(k1,ii,mm))
          bcjuse2=ccmax(bcj(kk,i1,mm),bcj(k1,i1,mm))
          bcjuse =ccmax(bcjuse1,bcjuse2)
          blankuse1=ccmax(blank(jw1,kk,ii),blank(jw1,k1,ii))
          blankuse2=ccmax(blank(jw1,kk,i1),blank(jw1,k1,i1))
          blankuse =ccmax(blankuse1,blankuse2)
          factr = slarge*((1. - blankuse)*hole +
     .                    (1. -   bcjuse))
          xs(1,i) = x(jw,kw,iw)
          xs(2,i) = y(jw,kw,iw)
          xs(3,i) = z(jw,kw,iw)
          xs(4,i) = factr
          ixs(1,i) = iw
          ixs(2,i) = jw
          ixs(3,i) = kw
          ixs(4,i) = nbl
        enddo
      enddo
      return
c
  500 kw = kmn
      kw1 = kmn
      go to 601
  600 kw = kmx
      kw1 = kmx - 1
  601 do jw=n2beg,n2end
        jj = min(jw,n2end-1) 
        j1 = max(jw-1,n2beg)
        do iw=n1beg,n1end
          ii = min(iw,n1end-1)
          i1 = max(iw-1,n1beg)
          i = 1 +iw-n1beg +ipts*(jw-n2beg)
          bckuse1=ccmax(bck(jj,ii,mm),bck(j1,ii,mm))
          bckuse2=ccmax(bck(jj,i1,mm),bck(j1,i1,mm))
          bckuse =ccmax(bckuse1,bckuse2)
          blankuse1=ccmax(blank(jj,kw1,ii),blank(j1,kw1,ii))
          blankuse2=ccmax(blank(jj,kw1,i1),blank(j1,kw1,i1))
          blankuse =ccmax(blankuse1,blankuse2)
          factr = slarge*((1. - blankuse)*hole +
     .                    (1. -   bckuse))
          xs(1,i) = x(jw,kw,iw)
          xs(2,i) = y(jw,kw,iw)
          xs(3,i) = z(jw,kw,iw)
          xs(4,i) = factr
          ixs(1,i) = iw
          ixs(2,i) = jw
          ixs(3,i) = kw
          ixs(4,i) = nbl
        enddo
      enddo
      return
      end
      subroutine bbdist(ng,grid,nsurf,surf,nbb,bbdef,ipv,vlist,dist,
     .                  idist,ncalc,wk3d5,nou,bou,nbuf,ibufdim,myid)
c***********************************************************************
c     Purpose:  Driver routine for determining the nearest bounding
c     box to each field point
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension grid(ng,3),surf(nsurf,4),dist(ng),idist(ng)
      dimension bbdef(6,nbb),ipv(2,nbb)
      dimension wk3d5(*)
      integer   vlist(nsurf)
      integer wrk,tsurf,test
      wrk = ifalloc(nou,bou,nbuf,ibufdim,myid,nbb)
      tsurf = ifalloc(nou,bou,nbuf,ibufdim,myid,4*nsurf)
      test = ifalloc(nou,bou,nbuf,ibufdim,myid,nsurf)
      call bbdst1(ng,grid,nsurf,surf,
     .           nbb,bbdef,ipv,vlist,dist,idist,ncalc,
     .           wk3d5(wrk),wk3d5(tsurf),wk3d5(test))
      call ffree(nou,bou,nbuf,ibufdim,myid,nsurf)
      call ffree(nou,bou,nbuf,ibufdim,myid,4*nsurf)
      call ffree(nou,bou,nbuf,ibufdim,myid,nbb)
c
      return
      end
      subroutine bbdst1(ng,grid,nsurf,surf,nbb,bbdef,ipv,vlist,dist,
     .                  idist,ncalc,wrk,tsurf,test)
c***********************************************************************
c     Purpose:  To identify nearest bounding box to each field point
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension grid(ng,3),surf(nsurf,4),dist(ng),idist(ng)
      dimension bbdef(6,nbb),ipv(2,nbb)
      dimension wrk(nbb), tsurf(nsurf,4), test(nsurf)
      integer   vlist(nsurf)
      do i=1,nsurf
      kk=vlist(i)
      tsurf(i,1)=surf(kk,1)
      tsurf(i,2)=surf(kk,2)
      tsurf(i,3)=surf(kk,3)
      tsurf(i,4)=surf(kk,4)
      end do
      ncalc=0
      do 100 i=1,ng
      x=grid(i,1)
      y=grid(i,2)
      z=grid(i,3)
      smin=1.0e34
      jp=0
c   build table of distances from grid to bounding boxes
      do 120 j=1,nbb
      px=x
      if (real(px) .le. real(bbdef(1,j))) px=bbdef(1,j)
      if (real(px) .ge. real(bbdef(2,j))) px=bbdef(2,j)
      py=y
      if (real(py) .le. real(bbdef(3,j))) py=bbdef(3,j)
      if (real(py) .ge. real(bbdef(4,j))) py=bbdef(4,j)
      pz=z
      if (real(pz) .le. real(bbdef(5,j))) pz=bbdef(5,j)
      if (real(pz) .ge. real(bbdef(6,j))) pz=bbdef(6,j)
      wrk(j)=(x-px)**2+(y-py)**2+(z-pz)**2
 120  continue
      do 200 j=1,nbb
c   find nearest bounding box that has not been searched
      bbmin=wrk(1)
      do ii=2,nbb
      bbmin=ccmin(bbmin,wrk(ii))
      end do
      jj=isrcheq(nbb,wrk,1,bbmin)
      wrk(jj)=2.0e34
c   stop searching when nearest bounding box is too far away
      if (real(bbmin) .gt. real(smin)) goto 201
      n=ipv(1,jj)
      l=ipv(2,jj)
      testmin=1.0e34
      do 160 k=1,n
      xs=tsurf(l,1)
      ys=tsurf(l,2)     
      zs=tsurf(l,3)
c     tsurf(l,4) contains 1) a large value for surface pts. that
c     should not be minimum distance pts (see note in subroutine
c     (collect_surf) or 2) zero for all regular surface pts.
      test(k)=(x-xs)**2+(y-ys)**2+(z-zs)**2+tsurf(l,4)
      testmin=ccmin(testmin,test(k))
      l=l+1
 160  continue
      ncalc=ncalc+n
      if (real(testmin) .lt. real(smin)) then
        kmin=isrcheq(n,test,1,testmin)
        smin=testmin
        jp=vlist(ipv(2,jj)+kmin-1)
      end if
 200  continue
 201  continue
      dist(i)=sqrt(smin)
      idist(i)=jp
 100  continue
      return
      end
      subroutine calc_dist(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                     imp2,jmp2,kmp2,x,y,z,smin,nsurf,surf,nbb,
     .                     bbdef,ipv,vlist,ntri,iptri,wk3d5,iwrk,
     .                     nou,bou,nbuf,ibufdim,myid)
c***********************************************************************
c     Purpose:  To find minimum distance to field points using the
c     recursive-box algorithm
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension x(jmp1,kmp1,imp1),smin(jmp2,kmp2,imp2)
      dimension y(jmp1,kmp1,imp1),z(jmp1,kmp1,imp1)
      dimension wk3d5(*), iwrk(*)
      integer vlist,grid,dist
      dimension surf(nsurf,4),bbdef(6,nbb),ipv(2,nbb),vlist(nsurf)
      dimension ntri(nsurf),iptri(nsurf,8)
      ng=(imx-imn+1)*(jmx-jmn+1)*(kmx-kmn+1)
      grid = ifalloc(nou,bou,nbuf,ibufdim,myid,3*ng)
      dist = ifalloc(nou,bou,nbuf,ibufdim,myid,ng)
      idist = iialloc(nou,bou,nbuf,ibufdim,myid,ng)
      ig=0
      do i=imn,imx
      do k=kmn,kmx
CDIR$ IVDEP
      do j=jmn,jmx
        ig=ig+1
        wk3d5(grid+ig-1)=x(j,k,i)
        wk3d5(grid+ig+ng-1)=y(j,k,i)
        wk3d5(grid+ig+2*ng-1)=z(j,k,i)
      end do
      end do
      end do
      call bbdist(ng,wk3d5(grid),nsurf,surf,
     .            nbb,bbdef,ipv,vlist,wk3d5(dist),iwrk(idist),
     .            ncalc,wk3d5,nou,bou,nbuf,ibufdim,myid)
c
c  put in calculation to triangles
c
      do i=1,ng
      isurf=iwrk(idist+i-1)
      xp=wk3d5(grid+i-1)
      yp=wk3d5(grid+i+ng-1)
      zp=wk3d5(grid+i+2*ng-1)
      numtri=ntri(isurf)
      if (numtri .eq. 4) then
        do itri=1,numtri
        it1=iptri(isurf,2*itri-1)
        it2=iptri(isurf,2*itri)
        call triang(wk3d5(dist+i-1),xp,yp,zp,
     .     surf(isurf,1),surf(isurf,2),surf(isurf,3),surf(isurf,4),
     .     surf(it1,1),surf(it1,2),surf(it1,3),surf(it1,4),
     .     surf(it2,1),surf(it2,2),surf(it2,3),surf(it2,4))
        end do
      else
        xps=surf(isurf,1)
        yps=surf(isurf,2)
        zps=surf(isurf,3)
        do ixmin=isurf,1,-1
        if (real(surf(ixmin,1)) .lt. real(xps)) go to 10
        end do
        ixmin=0
 10     continue
        ixmin=ixmin+1
        do ixmax=isurf,nsurf
        if (real(surf(ixmax,1)) .gt. real(xps)) go to 20
        end do
        ixmax=nsurf+1
 20     continue
        ixmax=ixmax-1
        do isurf=ixmin,ixmax
          if (real(surf(isurf,2)) .eq. real(yps) .and. 
     .        real(surf(isurf,3)) .eq. real(zps) .and.
     .        ntri(isurf) .ne. 4) then
            numtri=ntri(isurf)
            do itri=1,numtri
            it1=iptri(isurf,2*itri-1)
            it2=iptri(isurf,2*itri)
            call triang(wk3d5(dist+i-1),xp,yp,zp,
     .         surf(isurf,1),surf(isurf,2),surf(isurf,3),surf(isurf,4),
     .         surf(it1,1),surf(it1,2),surf(it1,3),surf(it1,4),
     .         surf(it2,1),surf(it2,2),surf(it2,3),surf(it2,4))
            end do
          end if
        end do
      end if
      end do
      ig=0
      do i=imn,imx
      do k=kmn,kmx
      do j=jmn,jmx
      ig=ig+1
      smin(j,k,i)=wk3d5(dist+ig-1)
      end do
      end do
      end do
      call ifree(nou,bou,nbuf,ibufdim,myid,ng)
      call ffree(nou,bou,nbuf,ibufdim,myid,ng)
      call ffree(nou,bou,nbuf,ibufdim,myid,3*ng)
      return
      end
      subroutine collect_surf(imn,jmn,kmn,imxs,jmxs,kmxs,imp1s,jmp1s,
     .                        kmp1s,xs,nface,n1beg,n1end,n2beg,n2end,
     .                        nsurf,surf,isurf,ntri,iptri)
c***********************************************************************
c     Purpose:  To store coordinates of surface points in the surf
c     array and to identify the neigbors of each surface pt.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension     xs(4,nsurf),surf(nsurf,4)
      dimension     ntri(nsurf),iptri(nsurf,8)
c
c Collect surface points.
c Now add collection of local triangles.
c Important side effect is that isurf is incremented properly
c
      isurf_beg=isurf
      ipts = n1end - n1beg + 1
c
      do k=n2beg,n2end
         do j=n1beg,n1end
            isurf=isurf+1
            i = 1 + (j-n1beg) + ipts*(k-n2beg)
            surf(isurf,1)    = xs(1,i)
            surf(isurf,2)    = xs(2,i)
            surf(isurf,3)    = xs(3,i)
            surf(isurf,4)    = xs(4,i)
         end do
      end do
c
c     The surf array (see above) contains the surface points.
c     Each interior surface point has 8 neighbors, which help define
c     surface triangles:
c
c            p1 --- p2 --- p3
c            |  \   |    / |
c            |   \  |   /  |
c            |    \ |  /   |
c            p4 --- sp --- p5
c            |    / |  \   |
c            |   /  |   \  |
c            |  /   |    \ |
c            p6 --- p7 --- p8
c
c     Points on a boundary will have fewer neighbors (as few as three
c     for a corner point.)
c     ntri(1,i) == the number of neighbor points for point i.
c     ntri(2,i) == the index in array iptri at which these neighbor
c                  points are defined.
c     iptri(ntri(2,i)),...iptri(ntri(2,i)+ntri(1,i)-1) == the indices
c           in the surf array of the neighbors of point i.
c
      isurf=isurf_beg
      n1inc=1
      n2inc=n1end-n1beg+1
      do j=n2beg,n2end
      do i=n1beg,n1end
        isurf=isurf+1
        itri=0
        do jj=max(j-1,n2beg),min(j+1,n2end)
        do ii=max(i-1,n1beg),min(i+1,n1end)
          if (ii .ne. i .and. jj .ne. j) then
          itri=itri+1
          iptri(isurf,2*itri-1) = isurf+(ii-i)*n1inc
          iptri(isurf,2*itri)   = isurf+(jj-j)*n2inc
          end if
        end do
        end do
        ntri(isurf)=itri
      end do
      end do
      return
      end
      subroutine distcc(jdim,kdim,idim,smingp,smincc) 
c***********************************************************************
c     Purpose: To average the distance function based at grid points 
c              to one based at cell centers
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension smingp(jdim,kdim,idim),smincc(jdim-1,kdim-1,idim-1)
c
      iw = 0
      do 2000 i=1,idim-1
      iw = iw + 1
      kw = 0
      do 2000 k=1,kdim-1
      kw = kw + 1
      jw = 0
      do 2000 j=1,jdim-1
      jw = jw + 1
      smincc(jw,kw,iw) = .125*(smingp(j,k,i)    +smingp(j,k,i+1)
     .                       + smingp(j+1,k,i)  +smingp(j+1,k,i+1)
     .                       + smingp(j,k+1,i)  +smingp(j,k+1,i+1)
     .                       + smingp(j+1,k+1,i)+smingp(j+1,k+1,i+1))
 2000 continue
c
      return
      end
      subroutine distcg(jdimf,kdimf,idimf,sminf,jdim,kdim,idim,sminc) 
c***********************************************************************
c     Purpose: To collocate the fine-grid minimum distance function
c              on a fine grid to a coarser grid
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension sminf(jdimf,kdimf,idimf),sminc(jdim,kdim,idim)
c
      ii   = 0
      iinc = 2
      if (idimf.eq.2) iinc = 1
      do 10 i=1,idimf,iinc
      ii   = ii+1
      kk   = 0
      do 10 k=1,kdimf,2
      kk   = kk+1
      jj   = 0
      do 10 j=1,jdimf,2
      jj   = jj+1
      sminc(jj,kk,ii) = sminf(j,k,i)
   10 continue
      return
      end
      subroutine initf(isize)
c***********************************************************************
c     Purpose:  To initialize pointers to floating-point variables in
c     recursive-box algorithm
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /alloc/ iptr,imax,ifptr,ifmax
c
      ifptr = 1
      ifmax = isize
c
      return
      end
      subroutine initi(isize)
c***********************************************************************
c     Purpose:  To initialize pointers to interger variables in
c     recursive-box algorithm
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /alloc/ iptr,imax,ifptr,ifmax
c
      iptr = 1
      imax = isize
c
      return
      end
      subroutine ifree(nou,bou,nbuf,ibufdim,myid,isize)
c***********************************************************************
c     Purpose:  To free up pointers to interger variables in
c     recursive-box algorithm
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /alloc/ iptr,imax,ifptr,ifmax
c
      if (isize.ge.iptr) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*) 'stopping...ifree failed: ',isize,iptr
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      endif
      iptr = iptr - isize
c
      return
      end
      subroutine ffree(nou,bou,nbuf,ibufdim,myid,isize)
c***********************************************************************
c     Purpose:  To free up pointers to floating-point variables in
c     recursive-box algorithm
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /alloc/ iptr,imax,ifptr,ifmax
c
      if (isize.ge.ifptr) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*) 'stopping...ffree failed: ',isize,ifptr
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      endif
      ifptr = ifptr - isize
c
      return
      end
       subroutine makebb(minpb0,maxlv0,nwv0,wx,wy,wz,
     .                 nbb0,bblen,bbdef,iv,vlist,iwrk,
     .                 nou,bou,nbuf,ibufdim,myid)
c***********************************************************************
c     Purpose:  Driver routine for generating the bounding boxes for the
c     recursive-box algorithm
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
        character*120 bou(ibufdim,nbuf)
c
        dimension nou(nbuf)
	dimension wx(*), wy(*), wz(*)
	dimension iwrk(*)
        dimension bbdef(6,*)
	integer bblen,iv(2,*), vlist(*)

	integer totbbv
	common /bbcom/ nbbv,minpbb,maxlev,nwv,nbb,totbbv,
     .                 bbmin(3), bbmax(3)
	integer bbi

	minpbb = minpb0
	maxlev = maxlv0
	nwv = nwv0
	nbb = 0
	totbbv = 0
c
c       The array bbi is a list of the surface points in wx,wy,wz.
c
        bbi = iialloc(nou,bou,nbuf,ibufdim,myid,nwv)
	do 10 i=1,nwv
	   iwrk(bbi+i-1) = i
 10	continue
c
c       Create a bounding box containing all the surface points.
c
	call calcbb(bbmin,bbmax,nwv,iwrk(bbi),wx,wy,wz)

c
c       Now subdivide the bounding box. Subdivision occurs maxlev
c       times, or until there are less than minpbb points in the
c       bounding box, whichever comes first.
c
	level = 0
	call spltbb(level,bbdef,iv,vlist,iwrk(bbi),wx,wy,wz)
	call ifree(nou,bou,nbuf,ibufdim,myid,nwv)
        nbb0 = nbb
	bblen = totbbv
c
	return
	end
        subroutine calcbb(bbmin,bbmax,nvi,vi,wx,wy,wz)
c***********************************************************************
c       Purpose:  Calculate bounding boxes
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
	parameter(bval=1.0e20, sval=-1.0e20)
c
	dimension bbmin(3), bbmax(3)
	integer vi(*)
	dimension wx(*), wy(*), wz(*)
c
c       vi is an array of nvi indices into the wx,wy,wz arrays.
c       These indices represent nvi points. This subroutine
c       creates a bounding box that contains these nvi points.
c
	do 10 i=1,3
	   bbmin(i) = bval
	   bbmax(i) = sval
 10	continue

	do 20 j=1,nvi
	   i = vi(j)
	   if(real(wx(i)).lt.real(bbmin(1))) bbmin(1) = wx(i)
	   if(real(wy(i)).lt.real(bbmin(2))) bbmin(2) = wy(i)
	   if(real(wz(i)).lt.real(bbmin(3))) bbmin(3) = wz(i)
	   if(real(wx(i)).gt.real(bbmax(1))) bbmax(1) = wx(i)
	   if(real(wy(i)).gt.real(bbmax(2))) bbmax(2) = wy(i)
	   if(real(wz(i)).gt.real(bbmax(3))) bbmax(3) = wz(i)
 20	continue

	return
	end
        subroutine getvrt(bbmin,bbmax,nwv,nbbv,bbi,wx,wy,wz)
c***********************************************************************
c       Purpose:  To search the wx,wy,wz arrays to find all points that
c       fall within a bounding box
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
	dimension bbmin(3), bbmax(3)
	integer bbi(*)
	dimension wx(*), wy(*), wz(*)
c
	nbbv = 0
	do 10 i=1,nwv
	   if((real(wx(i)).ge.real(bbmin(1))) .and.
     .        (real(wx(i)).le.real(bbmax(1))) .and.
     .        (real(wy(i)).ge.real(bbmin(2))) .and.
     .        (real(wy(i)).le.real(bbmax(2))) .and.
     .        (real(wz(i)).ge.real(bbmin(3))) .and.
     .        (real(wz(i)).le.real(bbmax(3)))) then
              nbbv = nbbv + 1
              bbi(nbbv) = i
	   endif
 10	continue

	return
	end
        subroutine shells(n,is,s)
c***********************************************************************
c       Purpose:  Shell sort (see Knuth)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
	dimension is(*), s(*)
c
	integer gap,temp
c
	gap = n / 2
 100	if(real(gap).le.0) goto 999
	   do 10 i=gap+1,n
		j = i - gap
 101		if (j.le.0) goto 10
		if (real(s(is(j))).le.real(s(is(j+gap)))) goto 10
		   temp = is(j)
		   is(j) = is(j+gap)
		   is(j+gap) = temp
		   j = j - gap
		   goto 101
 10	   continue
	   gap = gap / 2
	   goto 100
 999	continue

	return
	end
        subroutine spltbb(level0,bbdef,iv,vlist,bbi,wx,wy,wz)
c***********************************************************************
c       Purpose:  To subdivide bounding boxes
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
	parameter(MAXBB=200)
c
	dimension bbdef(6,*)
	integer bbi(*), iv(2,*), vlist(*)
	dimension wx(*), wy(*), wz(*)
c
	integer totbbv
	common /bbcom/ nbbv,minpbb,maxlev,nwv,nbb,totbbv,
     .                 bbmin(3), bbmax(3)
	dimension stack(7,MAXBB)
	dimension parmin(3), parmax(3)
	dimension plmin(3), plmax(3)
	dimension prmin(3), prmax(3)
c
	isptr = 1
	call push(stack,isptr,bbmin,bbmax,level0)
	k = 1
 100	if (isptr.le.1) goto 999
	   call pop(stack,isptr,parmin,parmax,level)
	   call getvrt(parmin,parmax,nwv,nbbv,bbi,wx,wy,wz)
	   if (nbbv.ge.minpbb.and.level.lt.maxlev) then
		bx = parmax(1) - parmin(1)
		by = parmax(2) - parmin(2)
		bz = parmax(3) - parmin(3)
		if (real(bx).gt.real(by).and.real(bx).gt.real(bz)) then
		   call shells(nbbv,bbi,wx)
		   iright = nbbv/2
		   i1 = bbi(iright)
		   i2 = bbi(iright+1)
 101               if(real(wx(i1)).ne.real(wx(i2)) .or. 
     .                iright.ge.nbbv)goto 102
			iright = iright + 1
                        if (iright .ge. nbbv) goto 102
		        i1 = bbi(iright)
		        i2 = bbi(iright+1)
			goto 101
 102		   continue
		elseif (real(by).gt.real(bx) .and. 
     .                  real(by).gt.real(bz)) then
		   call shells(nbbv,bbi,wy)
		   iright = nbbv/2
		   i1 = bbi(iright)
		   i2 = bbi(iright+1)
 103               if(real(wy(i1)).ne.real(wy(i2)) .or.
     .                iright.ge.nbbv)goto 104
			iright = iright + 1
                        if (iright .ge. nbbv) goto 104
		        i1 = bbi(iright)
		        i2 = bbi(iright+1)
			goto 103
 104		   continue
		else
		   call shells(nbbv,bbi,wz)
		   iright = nbbv/2
		   i1 = bbi(iright)
		   i2 = bbi(iright+1)
 105               if(real(wz(i1)).ne.real(wz(i2)) .or.
     .                iright.ge.nbbv)goto 106
			iright = iright + 1
                        if (iright .ge. nbbv) goto 106
		        i1 = bbi(iright)
		        i2 = bbi(iright+1)
			goto 105
 106		   continue
		endif
		call calcbb(plmin,plmax,iright,bbi,wx,wy,wz)
		call calcbb(prmin,prmax,nbbv-iright,bbi(iright+1),
     .			    wx,wy,wz)
		call push(stack,isptr,prmin,prmax,level+1)
		call push(stack,isptr,plmin,plmax,level+1)
	   else
		if (nbbv.gt.0) then
		   nbb = nbb + 1
		   bbdef(1,nbb) = parmin(1)
		   bbdef(2,nbb) = parmax(1)
		   bbdef(3,nbb) = parmin(2)
		   bbdef(4,nbb) = parmax(2)
		   bbdef(5,nbb) = parmin(3)
		   bbdef(6,nbb) = parmax(3)
		   iv(1,nbb) = nbbv
		   iv(2,nbb) = k
		   do 10 i=1,nbbv
			vlist(k) = bbi(i)
			k = k + 1
			totbbv = totbbv + 1
 10		   continue
		endif
	   endif
	   goto 100
 999	continue

	return
	end
        subroutine push(st,is,bbmin,bbmax,lev)
c***********************************************************************
c       Purpose:  Place an item into the stack, and adust the pointer
c       accordingly
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
	dimension st(7,*), bbmin(3),bbmax(3)
c
	do 10 i=1,3
	   ii = i+3
	   st(i,is) = bbmin(i)
	   st(ii,is) = bbmax(i)
 10	continue
	st(7,is) = lev + 0.1
	is = is + 1
c
	return
	end
        subroutine pop(st,is,bbmin,bbmax,lev)
c***********************************************************************
c       Purpose:  Remove an item from the stack, and adust the pointer
c       accordingly
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
	dimension st(7,*), bbmin(3),bbmax(3)
c
	is = is - 1
	do 10 i=1,3
	   ii = i+3
	   bbmin(i) = st(i,is)
	   bbmax(i) = st(ii,is)
 10	continue
	lev = st(7,is)
c
	return
	end
      subroutine sort_x(nsurf,surf,ntri,iptr,wk3d5,iwrk,
     .                  nou,bou,nbuf,ibufdim,myid)
c***********************************************************************
c     Purpose: To  sort surface points with respect to x coordinate
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension surf(nsurf,4),ntri(nsurf),iptr(nsurf,8)
      dimension wk3d5(*), iwrk(*)
      integer temp
      iperm = iialloc(nou,bou,nbuf,ibufdim,myid,nsurf)
      temp = ifalloc(nou,bou,nbuf,ibufdim,myid,nsurf)
      call heap_sort(nsurf,surf,iwrk(iperm))
      do j=1,4
         call move_real(nsurf,surf(1,j),iwrk(iperm),wk3d5(temp))
      end do
      call ffree(nou,bou,nbuf,ibufdim,myid,nsurf)
      itemp = iialloc(nou,bou,nbuf,ibufdim,myid,nsurf)
      call move_integer(nsurf,ntri,iwrk(iperm),iwrk(itemp))
      do j=1,8
         call move_integer(nsurf,iptr(1,j),iwrk(iperm),iwrk(itemp))
      end do
      call ifree(nou,bou,nbuf,ibufdim,myid,nsurf)
      inv_iperm = iialloc(nou,bou,nbuf,ibufdim,myid,nsurf)
      do i=1,nsurf
	 ii = iwrk(iperm+i-1)
         iwrk(inv_iperm+ii-1)=i
      end do
      do i=1,nsurf
         numtri=ntri(i)
         do j=1,2*numtri
	    ii = iptr(i,j)
            iptr(i,j)=iwrk(inv_iperm+ii-1)
         end do
      end do
      call ifree(nou,bou,nbuf,ibufdim,myid,nsurf)
      call ifree(nou,bou,nbuf,ibufdim,myid,nsurf)
      return
      end
      subroutine move_real(n,x,iperm,temp)
c***********************************************************************
c     Purpose:  Rearrange items in the real array x based on the
c     pointer iperm
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(n),iperm(n),temp(n)
      do i=1,n
      temp(i)=x(iperm(i))
      end do
      do i=1,n
      x(i)=temp(i)
      end do
      return
      end
      subroutine move_integer(n,ix,iperm,itemp)
c***********************************************************************
c     Purpose:  Rearrange items in the integer array x based on the
c     pointer iperm
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension ix(n),iperm(n),itemp(n)
      do i=1,n
      itemp(i)=ix(iperm(i))
      end do
      do i=1,n
      ix(i)=itemp(i)
      end do
      return
      end
      subroutine triang(dist,pp1,pp2,pp3,
     .aa1,aa2,aa3,aa4,bb1,bb2,bb3,bb4,cc1,cc2,cc3,cc4)
c***********************************************************************
c     Purpose:  To find the closest distance from a field point to the
c     actual surface (i.e. not simply the closest discrete surface
c     point), using local trangulation of the surface.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      if(real(aa4).ne.0..or.real(bb4).ne.0..or.real(cc4).ne.0.) return
c
      a1=bb1-aa1
      a2=bb2-aa2
      a3=bb3-aa3
      b1=cc1-aa1
      b2=cc2-aa2
      b3=cc3-aa3
      p1=pp1-aa1
      p2=pp2-aa2
      p3=pp3-aa3
      daa=a1*a1+a2*a2+a3*a3
      dab=a1*b1+a2*b2+a3*b3
      dbb=b1*b1+b2*b2+b3*b3
      den=dab**2-daa*dbb
      if (real(den) .eq. 0) go to 100
      dap=a1*p1+a2*p2+a3*p3
      dbp=b1*p1+b2*p2+b3*p3
      s=(dab*dbp-dbb*dap)/den
      t=(dab*dap-daa*dbp)/den
      if (real(s) .lt. 0 .or. real(t) .lt. 0 .or. 
     .    real(t+s) .gt. 1) go to 100
      r1=p1-s*a1-t*b1
      r2=p2-s*a2-t*b2
      r3=p3-s*a3-t*b3
      dist=ccmin(dist,sqrt(r1**2+r2**2+r3**2))
      return
 100  continue
      dsq=dist*dist
      if (real(daa) .eq. 0) go to 200
      dap=a1*p1+a2*p2+a3*p3
      t=dap/daa
      if (real(t) .lt. 0 .or. real(t) .gt. 1) go to 200
      r1=p1-t*a1
      r2=p2-t*a2
      r3=p3-t*a3
      dsq=ccmin(dsq,r1**2+r2**2+r3**2)
 200  continue
      if (real(dbb) .eq. 0) go to 300
      dbp=b1*p1+b2*p2+b3*p3
      t=dbp/dbb
      if (real(t) .lt. 0. .or. real(t) .gt. 1) go to 300
      r1=p1-t*b1
      r2=p2-t*b2
      r3=p3-t*b3
      dsq=ccmin(dsq,r1**2+r2**2+r3**2)
 300  continue
      p1=pp1-bb1
      p2=pp2-bb2
      p3=pp3-bb3
      a1=cc1-bb1
      a2=cc2-bb2
      a3=cc3-bb3
      daa=a1*a1+a2*a2+a3*a3
      if (real(daa) .eq. 0) go to 400
      dap=a1*p1+a2*p2+a3*p3
      t=dap/daa
      if (real(t) .lt. 0 .or. real(t) .gt. 1) go to 400
      r1=p1-t*a1
      r2=p2-t*a2
      r3=p3-t*a3
      dsq=ccmin(dsq,r1**2+r2**2+r3**2)
400   continue
      if (real(dsq) .lt. real(dist)*real(dist)) dist=sqrt(dsq)
      return
      end
      subroutine heap_sort(n,arrin,indx)
c***********************************************************************
c     Purpose:  To sort a list of points 
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension arrin(n)
      integer indx(n)
c      
      do 11 j=1,n
         indx(j) = j
 11   continue
c
      if(n .eq. 1)return
c      
      l = n/2 + 1
      ir= n
 10   continue
      if(l.gt.1)then
         l=l-1
         indxt = indx(l)
         q     = arrin(indxt)
      else
         indxt = indx(ir)
         q     = arrin(indxt)
         indx(ir) = indx(1)
         ir = ir - 1
         if(ir.eq.1)then
            indx(1) = indxt
            return
         endif
      endif
      i=l
      j=l+l
 20   if(j.le.ir)then
         if(j.lt.ir)then
            if(real(arrin(indx(j))).lt.real(arrin(indx(j+1))))j=j+1
         endif
         if(real(q).lt.real(arrin(indx(j))))then
            indx(i) = indx(j)
            i = j
            j = j+j
         else
            j = ir+1
         endif
         go to 20
      endif
      indx(i) = indxt
      go to 10
      end
      subroutine bbarthcg(jdim,kdim,idim,jdimc,kdimc,idimc,xjb,xkb,xib,
     .                    blnum,xjbc,xkbc,xibc,blnumc)
c***********************************************************************
c     Purpose:  To set the i,j,k and block number indicies needed for 
c     the Baldwin-Barth model on coarser levels
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xjb(jdim-1,kdim-1,idim-1),xkb(jdim-1,kdim-1,idim-1),
     .          xib(jdim,kdim,idim,2),blnum(jdim-1,kdim-1,idim-1),
     .          xjbc(jdimc-1,kdimc-1,idimc-1),
     .          xkbc(jdimc-1,kdimc-1,idimc-1),
     .          xibc(jdimc,kdimc,idimc,2),
     .          blnumc(jdimc-1,kdimc-1,idimc-1) 
c
      jdim1 = jdim - 1
      kdim1 = kdim - 1
      idim1 = idim - 1
c
      iinc = 2
      if (idim.eq.2) iinc = 1
      ii = 0
      do 10 i=1,idim1,iinc
      ii = ii + 1
      kk = 0
      do 10 k=1,kdim1,2
      kk = kk + 1
      jj = 0
      do 10 j=1,jdim1,2
      jj = jj + 1
         blnumc(jj,kk,ii) = blnum(j,k,i) + 1
         xjbc(jj,kk,ii)   = int(xjb(j,k,i)/2 + 1)
         xkbc(jj,kk,ii)   = int(xkb(j,k,i)/2 + 1)
         xibc(jj,kk,ii,2) = int(xib(j,k,i,2)/2 + 1)
  10  continue
c
      return
      end
      subroutine calc_distbb(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                       imp2,jmp2,kmp2,x,y,z,smin,nsurf,surf,nbb,
     .                       bbdef,ipv,vlist,ntri,iptri,ibbarth,
     .                       xjb,xkb,xib,blnum,wk3d5,iwrk,j1,
     .                       nou,bou,nbuf,ibufdim,myid)
c***********************************************************************
c     Purpose:  To find minimum distance to field points using the
c     recursive-box algorithm; this routine also sets the corresponding
c     indicies i,j,k,blnum for use with the Baldwin-Barth turb. model
c***********************************************************************
c
c     this routine is a modified version of calc_dist for
c     the baldwin-barth model - mods by r. t. biedron
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension ibbarth(nsurf,4),xjb(jmp1-1,kmp1-1,imp1-1),
     .          xkb(jmp1-1,kmp1-1,imp1-1),xib(jmp1,kmp1,imp1,2),
     .          blnum(jmp1-1,kmp1-1,imp1-1)
      dimension x(jmp1,kmp1,imp1),smin(jmp2,kmp2,imp2)
      dimension y(jmp1,kmp1,imp1),z(jmp1,kmp1,imp1)
      dimension wk3d5(*), iwrk(*)
      integer vlist,grid,dist
      dimension surf(nsurf,4),bbdef(6,nbb),ipv(2,nbb),vlist(nsurf)
      dimension ntri(nsurf),iptri(nsurf,8)
      ng=(imx-imn+1)*(jmx-jmn+1)*(kmx-kmn+1)
      grid = ifalloc(nou,bou,nbuf,ibufdim,myid,3*ng)
      dist = ifalloc(nou,bou,nbuf,ibufdim,myid,ng)
      idist = iialloc(nou,bou,nbuf,ibufdim,myid,ng)
      ig=0
      do i=imn,imx
      do k=kmn,kmx
CDIR$ IVDEP
      do j=jmn,jmx
        ig=ig+1
        wk3d5(grid+ig-1)=x(j,k,i)
        wk3d5(grid+ig+ng-1)=y(j,k,i)
        wk3d5(grid+ig+2*ng-1)=z(j,k,i)
      end do
      end do
      end do
      call bbdist(ng,wk3d5(grid),nsurf,surf,
     .            nbb,bbdef,ipv,vlist,wk3d5(dist),iwrk(idist),
     .            ncalc,wk3d5,nou,bou,nbuf,ibufdim,myid)
c
c     note: min. distance and ibbarth pointers are calculated/set
c     for grid points. what cfl3d ultimately needs are cell-center
c     values. the min. distance value is averaged elsewhere (by 
c     calling subroutine distcc in findmin_new). for the xib, xjb
c     xkb and blnum pointers, we can't really average, so just take
c     the easy way out and say that the cell center value of these
c     pointers are the same as the grid point values with the same
c     index. note also that iwrk(i) contains imx1*jmx1*kmx1 grid-
c     point values, whereas xib....blnum can only be set for 
c     (imx1-1)*(jmx1-1)*(kmx1-1) values, so to keep the conversion
c     straight, must loop over all j,k but skip when at the maximum
c
      ifactr = 2**j1
      imx1 = imx - 1
      jmx1 = jmx - 1
      kmx1 = kmx - 1
      ig=0
      do i=imn,imx1
      do k=kmn,kmx
CDIR$ IVDEP
      do j=jmn,jmx
        ig=ig+1
        if (j.eq.jmx .or. k.eq.kmx) cycle
        isurf=iwrk(idist+ig-1)
        xib(j,k,i,2) = int((ibbarth(isurf,1)-1)/ifactr + 1)
        xjb(j,k,i)   = int((ibbarth(isurf,2)-1)/ifactr + 1)
        xkb(j,k,i)   = int((ibbarth(isurf,3)-1)/ifactr + 1)
        blnum(j,k,i) = ibbarth(isurf,4) + j1
      end do
      end do
      end do
c
c  put in calculation to triangles
c
      do i=1,ng
      isurf=iwrk(idist+i-1)
      xp=wk3d5(grid+i-1)
      yp=wk3d5(grid+i+ng-1)
      zp=wk3d5(grid+i+2*ng-1)
      numtri=ntri(isurf)
      if (numtri .eq. 4) then
        do itri=1,numtri
        it1=iptri(isurf,2*itri-1)
        it2=iptri(isurf,2*itri)
        call triang(wk3d5(dist+i-1),xp,yp,zp,
     .      surf(isurf,1),surf(isurf,2),surf(isurf,3),surf(isurf,4),
     .      surf(it1,1),surf(it1,2),surf(it1,3),surf(it1,4),
     .      surf(it2,1),surf(it2,2),surf(it2,3),surf(it2,4))
        end do
      else
        xps=surf(isurf,1)
        yps=surf(isurf,2)
        zps=surf(isurf,3)
        do ixmin=isurf,1,-1
        if (real(surf(ixmin,1)) .lt. real(xps)) go to 10
        end do
        ixmin=0
 10     continue
        ixmin=ixmin+1
        do ixmax=isurf,nsurf
        if (real(surf(ixmax,1)) .gt. real(xps)) go to 20
        end do
        ixmax=nsurf+1
 20     continue
        ixmax=ixmax-1
        do isurf=ixmin,ixmax
          if (real(surf(isurf,2)) .eq. real(yps) .and. 
     .        real(surf(isurf,3)) .eq. real(zps) .and.
     .        ntri(isurf) .ne. 4) then
              numtri=ntri(isurf)
              do itri=1,numtri
              it1=iptri(isurf,2*itri-1)
              it2=iptri(isurf,2*itri)
              call triang(wk3d5(dist+i-1),xp,yp,zp,
     .          surf(isurf,1),surf(isurf,2),surf(isurf,3),surf(isurf,4),
     .          surf(it1,1),surf(it1,2),surf(it1,3),surf(it1,4),
     .          surf(it2,1),surf(it2,2),surf(it2,3),surf(it2,4))
            end do
          end if
        end do
      end if
      end do
      ig=0
      do i=imn,imx
      do k=kmn,kmx
      do j=jmn,jmx
      ig=ig+1
      smin(j,k,i)=wk3d5(dist+ig-1)
      end do
      end do
      end do
      call ifree(nou,bou,nbuf,ibufdim,myid,ng)
      call ffree(nou,bou,nbuf,ibufdim,myid,ng)
      call ffree(nou,bou,nbuf,ibufdim,myid,3*ng)
      return
      end
      subroutine collect_surfbb(imn,jmn,kmn,imxs,jmxs,kmxs,imp1s,jmp1s,
     .                          kmp1s,xs,ixs,nface,n1beg,n1end,n2beg,
     .                          n2end,nsurf,surf,isurf,ntri,iptri,
     .                          ibbarth,nbl)
c***********************************************************************
c     Purpose:  To store coordinates of surface points in the surf
c     array and to identify the neigbors of each surface pt. this 
c     routine also stores the additional data needed for the Baldwin-
c     Barth turbulence model.
c***********************************************************************
c
c     this routine is a modified version of collect_surf for
c     the baldwin-barth model - mods by r. t. biedron
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension     xs(4,nsurf),ixs(4,nsurf),surf(nsurf,4)
      dimension     ntri(nsurf),iptri(nsurf,8)
      dimension     ibbarth(nsurf,4)
c
c Collect surface points.
c Now add collection of local triangles.
c Important side effect is that isurf is incremented properly
c
      isurf_beg=isurf
      ipts = n1end - n1beg + 1
c
      do k=n2beg,n2end
         do j=n1beg,n1end
            isurf=isurf+1
            i = 1 + (j-n1beg) + ipts*(k-n2beg)
            surf(isurf,1)    = xs(1,i)
            surf(isurf,2)    = xs(2,i)
            surf(isurf,3)    = xs(3,i)
            surf(isurf,4)    = xs(4,i)
            ibbarth(isurf,1) = ixs(1,i)
            ibbarth(isurf,2) = ixs(2,i)
            ibbarth(isurf,3) = ixs(3,i)
            ibbarth(isurf,4) = ixs(4,i)
         end do
      end do
c
c
c     The surf array (see above) contains the surface points.
c     Each interior surface point has 8 neighbors, which help define
c     surface triangles:
c
c            p1 --- p2 --- p3
c            |  \   |    / |
c            |   \  |   /  |
c            |    \ |  /   |
c            p4 --- sp --- p5
c            |    / |  \   |
c            |   /  |   \  |
c            |  /   |    \ |
c            p6 --- p7 --- p8
c
c     Points on a boundary will have fewer neighbors (as few as three
c     for a corner point.)
c     ntri(1,i) == the number of neighbor points for point i.
c     ntri(2,i) == the index in array iptri at which these neighbor
c                  points are defined.
c     iptri(ntri(2,i)),...iptri(ntri(2,i)+ntri(1,i)-1) == the indices
c           in the surf array of the neighbors of point i.
c
      isurf=isurf_beg
      n1inc=1
      n2inc=n1end-n1beg+1
      do j=n2beg,n2end
      do i=n1beg,n1end
        isurf=isurf+1
        itri=0
        do jj=max(j-1,n2beg),min(j+1,n2end)
        do ii=max(i-1,n1beg),min(i+1,n1end)
          if (ii .ne. i .and. jj .ne. j) then
          itri=itri+1
          iptri(isurf,2*itri-1) = isurf+(ii-i)*n1inc
          iptri(isurf,2*itri)   = isurf+(jj-j)*n2inc
          end if
        end do
        end do
        ntri(isurf)=itri
      end do
      end do
      return
      end
      subroutine sort_xbb(nsurf,surf,ntri,iptr,ibbarth,wk3d5,iwrk,
     .                    nou,bou,nbuf,ibufdim,myid)
c***********************************************************************
c     Purpose: To  sort surface points with respect to x coordinate
c***********************************************************************
c
c     this routine is a modified version of sort_x for
c     the baldwin-barth model - mods by r. t. biedron
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension surf(nsurf,4),ntri(nsurf),iptr(nsurf,8)
      dimension ibbarth(nsurf,4)
      dimension wk3d5(*), iwrk(*)
      integer temp
      iperm = iialloc(nou,bou,nbuf,ibufdim,myid,nsurf)
      temp = ifalloc(nou,bou,nbuf,ibufdim,myid,nsurf)
      call heap_sort(nsurf,surf,iwrk(iperm))
      do j=1,4
         call move_real(nsurf,surf(1,j),iwrk(iperm),wk3d5(temp))
      end do
      call ffree(nou,bou,nbuf,ibufdim,myid,nsurf)
      itemp = iialloc(nou,bou,nbuf,ibufdim,myid,nsurf)
      call move_integer(nsurf,ntri,iwrk(iperm),iwrk(itemp))
      do j=1,8
         call move_integer(nsurf,iptr(1,j),iwrk(iperm),iwrk(itemp))
      end do
      do j=1,4
      call move_integer(nsurf,ibbarth(1,j),iwrk(iperm),iwrk(itemp))
      end do 
      call ifree(nou,bou,nbuf,ibufdim,myid,nsurf)
      inv_iperm = iialloc(nou,bou,nbuf,ibufdim,myid,nsurf)
      do i=1,nsurf
	 ii = iwrk(iperm+i-1)
         iwrk(inv_iperm+ii-1)=i
      end do
      do i=1,nsurf
         numtri=ntri(i)
         do j=1,2*numtri
	    ii = iptr(i,j)
            iptr(i,j)=iwrk(inv_iperm+ii-1)
         end do
      end do
      call ifree(nou,bou,nbuf,ibufdim,myid,nsurf)
      call ifree(nou,bou,nbuf,ibufdim,myid,nsurf)
      return
      end
      integer function iialloc(nou,bou,nbuf,ibufdim,myid,isize)
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /alloc/ iptr,imax,ifptr,ifmax
c
      iialloc = iptr
      iptr = iptr + isize
      if (iptr.gt.(imax+1)) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*) 'stopping...iialloc failed: ',
     .   isize,(iptr-isize)
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      endif
c
      return
      end
      integer function ifalloc(nou,bou,nbuf,ibufdim,myid,isize)
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /alloc/ iptr,imax,ifptr,ifmax
c
      ifalloc = ifptr
      ifptr = ifptr + isize
      if (ifptr.gt.(ifmax+1)) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*) 'stopping...ifalloc failed: ',
     .   isize,(ifptr-isize)
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      endif
c
      return
      end
