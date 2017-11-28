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
      subroutine pltmode(nm,iaes,nblg,ngrid,maxgr,maxbl,
     .                   nsegdfrm,iaesurf,jbcinfo,kbcinfo,
     .                   ibcinfo,nbcj0,nbcjdim,nbck0,nbckdim,
     .                   nbci0,nbcidim,maxseg,maxsegdg,lw,lw2,
     .                   x,y,z,xmdj,xmdk,xmdi,jdim1,kdim1,idim1,
     .                   maxaes,nmds,nbl1,iunitw,iflag)
c
c     $Id$
c
c***********************************************************************
c     Purpose: To write the user-input modal surface out to a plot3d 
c     file that can subsequently be viewed in order to help sort out
c     any problems with the modal surface data.
c***********************************************************************
      character*21 titlemode
      dimension nblg(maxgr),nbci0(maxbl),nbcidim(maxbl),
     .          nbcj0(maxbl),nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          ibcinfo(maxbl,maxseg,7,2),jbcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2),nsegdfrm(maxbl),
     .          iaesurf(maxbl,maxsegdg),lw(65,maxbl),lw2(43,maxbl)
      dimension x(jdim1,kdim1,idim1),y(jdim1,kdim1,idim1),
     .          z(jdim1,kdim1,idim1)
      dimension xmdj(kdim1,idim1,6,nmds,maxaes),
     .          xmdk(jdim1,idim1,6,nmds,maxaes),
     .          xmdi(jdim1,kdim1,6,nmds,maxaes)
      dimension jdm(500),kdm(500),idm(500)
c
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /elastic/ ndefrm,naesrf
      common /igrdtyp/ ip3dgrd,ialph
c 
      if (iaes.gt.9) then
         if (nm.gt.99) then
            len1 = 20
            write (titlemode,'("aesurf",i2,"_mode",i3,".plt")') iaes,nm
         else if(nm.gt.9) then
            len1 = 19
            write (titlemode,'("aesurf",i2,"_mode",i2,".plt")') iaes,nm
         else
            len1 = 18
            write (titlemode,'("aesurf",i2,"_mode",i1,".plt")') iaes,nm
         end if
      else
         if (nm.gt.99) then
            len1 = 19
            write (titlemode,'("aesurf",i1,"_mode",i3,".plt")') iaes,nm
         else if(nm.gt.9) then
            len1 = 18
            write (titlemode,'("aesurf",i1,"_mode",i2,".plt")') iaes,nm
         else
            len1 = 17
            write (titlemode,'("aesurf",i1,"_mode",i1,".plt")') iaes,nm
         end if
      end if
      do i = len1+1, 21
         titlemode(i:i) = ' '
      end do
      open(iunitw,file=titlemode(1:len1),form='formatted',
     .     status='unknown')
      write(iunitw,23110)
23110             format(' title="tecplot modeshape file"',/,
     .              ' variables="x","y","z","phix","phiy","phiz"')
23111             format(' ZONE I = ',i5,' J = ',i5)
c
      if (iflag.eq.0) then
c
c        get number of zones and zone dimensions for plot3d header
c
         rewind(iunitw)
         ngd = 0
         do igrid = 1,ngrid
            nbl = nblg(igrid)
            iaesrf = 0
            do is=1,nsegdfrm(nbl)
               iaesrf = iaesrf + iaesurf(nbl,is)
            end do
            if (iaesrf.ne.0) then
               call lead(nbl,lw,lw2,maxbl)
               do m = 1,2
                  if (m.eq.1) then
                     j    = 1
                     l    = 0
                     nseg = nbcj0(nbl)
                  else
                     j    = jdim
                     l    = 3
                     nseg = nbcjdim(nbl)
                  end if
                  do ns = 1,nseg
                     lbc = abs(jbcinfo(nbl,ns,1,m))
                     iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
     .                                *(lbc-2024)*(lbc-2034)*(lbc-2016)
                     if (iflg.eq.0) then
                        ngd = ngd + 1
                        ist = jbcinfo(nbl,ns,2,m)
                        ifn = jbcinfo(nbl,ns,3,m)
                        kst = jbcinfo(nbl,ns,4,m)
                        kfn = jbcinfo(nbl,ns,5,m)
                        idm(ngd) = ifn-ist+1
                        kdm(ngd) = kfn-kst+1
                        jdm(ngd) = 1
                     end if
                  end do
               end do
               do m = 1,2
                  if (m.eq.1) then
                     k    = 1
                     l    = 0
                     nseg = nbck0(nbl)
                  else
                     k    = kdim
                     l    = 3
                     nseg = nbckdim(nbl)
                  end if
                  do ns = 1,nseg
                    lbc = abs(kbcinfo(nbl,ns,1,m))
                    iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
     .                               *(lbc-2024)*(lbc-2034)*(lbc-2016)
                    if (iflg.eq.0) then
                        ngd = ngd + 1
                        ist = kbcinfo(nbl,ns,2,m)
                        ifn = kbcinfo(nbl,ns,3,m)
                        jst = kbcinfo(nbl,ns,4,m)
                        jfn = kbcinfo(nbl,ns,5,m)
                        idm(ngd) = ifn-ist+1
                        kdm(ngd) = 1
                        jdm(ngd) = jfn-jst+1
                    end if
                 end do
               end do
               do m = 1,2
                  if (m.eq.1) then
                     i    = 1
                     l    = 0
                     nseg = nbci0(nbl)
                  else
                     i    = idim
                     l    = 3
                     nseg = nbcidim(nbl)
                  end if
                  do ns = 1,nseg
                     lbc = abs(ibcinfo(nbl,ns,1,m))
                     iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
     .                                *(lbc-2024)*(lbc-2034)*(lbc-2016)
                     if (iflg.eq.0) then
                        ngd = ngd + 1
                        jst = ibcinfo(nbl,ns,2,m)
                        jfn = ibcinfo(nbl,ns,3,m)
                        kst = ibcinfo(nbl,ns,4,m)
                        kfn = ibcinfo(nbl,ns,5,m)
                        idm(ngd) = 1
                        kdm(ngd) = kfn-kst+1
                        jdm(ngd) = jfn-jst+1
                    end if
                 end do
               end do
            end if
         end do
c
      else if (iflag.gt.0) then
c
c        write out x,y,z of modal surface
c
         ngd = 0
         do m = 1,2
            if (m.eq.1) then
               j    = 1
               l    = 0
               nseg = nbcj0(nbl1)
            else
               j    = jdim1
               l    = 3
               nseg = nbcjdim(nbl1)
            end if
            do ns = 1,nseg
               lbc = abs(jbcinfo(nbl1,ns,1,m))
               iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
     .                          *(lbc-2024)*(lbc-2034)*(lbc-2016)
               if (iflg.eq.0) then
                  ngd = ngd + 1
                  ist = jbcinfo(nbl1,ns,2,m)
                  ifn = jbcinfo(nbl1,ns,3,m)
                  kst = jbcinfo(nbl1,ns,4,m)
                  kfn = jbcinfo(nbl1,ns,5,m)
                  write(iunitw,23111) ifn-ist+1,kfn-kst+1
                  if(ialph.eq.0) then
                    do k = kst,kfn 
                      do i = ist,ifn
                        write(iunitw,'(6(1x,e16.8))')
     .                      x(j,k,i),y(j,k,i),z(j,k,i),
     .                      real(xmdj(k,i,l+1,nm,iaes)),
     .                      real(xmdj(k,i,l+2,nm,iaes)),
     .                      real(xmdj(k,i,l+3,nm,iaes)) 
                      enddo
                    enddo
                  else
                    do k = kst,kfn 
                      do i = ist,ifn
                        write(iunitw,'(6(1x,e16.8))')
     .                      x(j,k,i),z(j,k,i),-y(j,k,i),
     .                      real(xmdj(k,i,l+1,nm,iaes)),
     .                      real(xmdj(k,i,l+3,nm,iaes)),
     .                     -real(xmdj(k,i,l+2,nm,iaes)) 
                      enddo
                    enddo
                  end if
               end if
            end do
         end do
         do m = 1,2
            if (m.eq.1) then
               k    = 1
               l    = 0
               nseg = nbck0(nbl1)
            else
               k    = kdim1
               l    = 3
               nseg = nbckdim(nbl1)
            end if
            do ns = 1,nseg
               lbc = abs(kbcinfo(nbl1,ns,1,m))
               iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
     .                          *(lbc-2024)*(lbc-2034)*(lbc-2016)
               if (iflg.eq.0) then
                  ngd = ngd + 1
                  ist = kbcinfo(nbl1,ns,2,m)
                  ifn = kbcinfo(nbl1,ns,3,m)
                  jst = kbcinfo(nbl1,ns,4,m)
                  jfn = kbcinfo(nbl1,ns,5,m)
                  write(iunitw,23111) ifn-ist+1,jfn-jst+1
                  if(ialph.eq.0) then
                    do j = jst,jfn
                      do i = ist,ifn
                        write(iunitw,'(6(1x,e16.8))')
     .                      x(j,k,i),y(j,k,i),z(j,k,i),
     .                      real(xmdk(j,i,l+1,nm,iaes)),
     .                      real(xmdk(j,i,l+2,nm,iaes)),
     .                      real(xmdk(j,i,l+3,nm,iaes))
                      enddo
                    enddo
                  else
                    do j = jst,jfn
                      do i = ist,ifn
                        write(iunitw,'(6(1x,e16.8))')
     .                      x(j,k,i),z(j,k,i),-y(j,k,i),
     .                      real(xmdk(j,i,l+1,nm,iaes)),
     .                      real(xmdk(j,i,l+3,nm,iaes)),
     .                     -real(xmdk(j,i,l+2,nm,iaes))
                      enddo
                    enddo
                  end if
              end if
           end do
         end do
         do m = 1,2
            if (m.eq.1) then
               i    = 1
               l    = 0
               nseg = nbci0(nbl1)
            else
               i    = idim1
               l    = 3
               nseg = nbcidim(nbl1)
            end if
            do ns = 1,nseg
               lbc = abs(ibcinfo(nbl1,ns,1,m))
               iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
     .                          *(lbc-2024)*(lbc-2034)*(lbc-2016)
               if (iflg.eq.0) then
                  ngd = ngd + 1
                  jst = ibcinfo(nbl1,ns,2,m)
                  jfn = ibcinfo(nbl1,ns,3,m)
                  kst = ibcinfo(nbl1,ns,4,m)
                  kfn = ibcinfo(nbl1,ns,5,m)
                  write(iunitw,23111) jfn-jst+1,kfn-kst+1
                  if(ialph.eq.0) then
                    do k = kst,kfn
                      do j = jst,jfn
                        write(iunitw,'(6(1x,e16.8))')
     .                      x(j,k,i),y(j,k,i),z(j,k,i),
     .                      real(xmdi(j,k,l+1,nm,iaes)),
     .                      real(xmdi(j,k,l+2,nm,iaes)),
     .                      real(xmdi(j,k,l+3,nm,iaes))
                      enddo
                    enddo
                  else
                    do k = kst,kfn
                      do j = jst,jfn
                        write(iunitw,'(6(1x,e16.8))')
     .                      x(j,k,i),z(j,k,i),-y(j,k,i),
     .                      real(xmdi(j,k,l+1,nm,iaes)),
     .                      real(xmdi(j,k,l+2,nm,iaes)),
     .                     -real(xmdi(j,k,l+3,nm,iaes))
                      enddo
                    enddo
                  end if
              end if
           end do
         end do
c
      end if
c
      return
      end
