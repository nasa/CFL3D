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
      dimension bcvali(maxbl,maxseg,7,2),bcvalj(maxbl,maxseg,7,2),
     .          nblg(maxgr),nbci0(maxbl),nbcidim(maxbl),
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
c 
      if (iaes.gt.9) then
         if (nm.gt.99) then
            len1 = 20
            write (titlemode,'("aesurf",i2,"_mode",i3,".p3d")') iaes,nm
         else if(nm.gt.9) then
            len1 = 19
            write (titlemode,'("aesurf",i2,"_mode",i2,".p3d")') iaes,nm
         else
            len1 = 18
            write (titlemode,'("aesurf",i2,"_mode",i1,".p3d")') iaes,nm
         end if
      else
         if (nm.gt.99) then
            len1 = 19
            write (titlemode,'("aesurf",i1,"_mode",i3,".p3d")') iaes,nm
         else if(nm.gt.9) then
            len1 = 18
            write (titlemode,'("aesurf",i1,"_mode",i2,".p3d")') iaes,nm
         else
            len1 = 17
            write (titlemode,'("aesurf",i1,"_mode",i1,".p3d")') iaes,nm
         end if
      end if
      do i = len1+1, 21
         titlemode(i:i) = ' '
      end do
      open(iunitw,file=titlemode(1:len1),form='formatted',
     .     status='unknown')
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
                     lbc = iabs(jbcinfo(nbl,ns,1,m))
                     if (naesrf.gt.0) lbc =
     .                  abs(jbcinfo(nbl,ns,1,m))*jbcinfo(nbl,ns,6,m)
                     iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
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
                     if (naesrf.gt.0) lbc =
     .                  abs(kbcinfo(nbl,ns,1,m))*kbcinfo(nbl,ns,6,m)
                    iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
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
                     if (naesrf.gt.0) lbc =
     .                  abs(ibcinfo(nbl,ns,1,m))*ibcinfo(nbl,ns,6,m)
                     iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
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
         write(iunitw,*) ngd
         write(iunitw,*) (idm(nn),jdm(nn),kdm(nn),nn=1,ngd)
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
               lbc = iabs(jbcinfo(nbl1,ns,1,m))
               if (naesrf.gt.0) lbc =
     .            abs(jbcinfo(nbl1,ns,1,m))*jbcinfo(nbl1,ns,6,m)
               iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
               if (iflg.eq.0) then
                  ngd = ngd + 1
                  ist = jbcinfo(nbl1,ns,2,m)
                  ifn = jbcinfo(nbl1,ns,3,m)
                  kst = jbcinfo(nbl1,ns,4,m)
                  kfn = jbcinfo(nbl1,ns,5,m)
                  write(iunitw,*)
     .                 ((x(j,k,i),i=ist,ifn),k=kst,kfn),
     .                 ((y(j,k,i),i=ist,ifn),k=kst,kfn),
     .                 ((real(xmdj(k,i,l+3,nm,iaes)),i=ist,ifn),
     .                                               k=kst,kfn)
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
               if (naesrf.gt.0) lbc =
     .            abs(kbcinfo(nbl1,ns,1,m))*kbcinfo(nbl1,ns,6,m)
               iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
               if (iflg.eq.0) then
                  ngd = ngd + 1
                  ist = kbcinfo(nbl1,ns,2,m)
                  ifn = kbcinfo(nbl1,ns,3,m)
                  jst = kbcinfo(nbl1,ns,4,m)
                  jfn = kbcinfo(nbl1,ns,5,m)
                  write(iunitw,*)
     .            ((x(j,k,i),i=ist,ifn),j=jst,jfn),
     .            ((y(j,k,i),i=ist,ifn),j=jst,jfn),
     .            ((real(xmdk(j,i,l+3,nm,iaes)),i=ist,ifn),
     .                                          j=jst,jfn)
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
               if (naesrf.gt.0) lbc =
     .            abs(ibcinfo(nbl1,ns,1,m))*ibcinfo(nbl1,ns,6,m)
               iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
               if (iflg.eq.0) then
                  ngd = ngd + 1
                  jst = ibcinfo(nbl1,ns,2,m)
                  jfn = ibcinfo(nbl1,ns,3,m)
                  kst = ibcinfo(nbl1,ns,4,m)
                  kfn = ibcinfo(nbl1,ns,5,m)
                  write(iunitw,*)
     .                 ((x(j,k,i),j=jst,jfn),k=kst,kfn),
     .                 ((y(j,k,i),j=jst,jfn),k=kst,kfn),
     .                 ((real(xmdi(j,k,l+3,nm,iaes)),j=jst,jfn),
     .                                               k=kst,kfn)
              end if
           end do
         end do
c
      end if
c
      return
      end
