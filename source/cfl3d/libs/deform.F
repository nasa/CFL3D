      subroutine deform(nbl,idim,jdim,kdim,x,y,z,xnm2,ynm2,znm2,
     .                  xnm1,ynm1,znm1,deltj,deltk,delti,
     .                  vel,icsi,icsf,jcsi,jcsf,kcsi,kcsf,
     .                  maxbl,time,dt,ita,nou,bou,nbuf,ibufdim,myid,
     .                  idefrm,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,
     .                  nbckdim,ibcinfo,jbcinfo,kbcinfo,maxseg,wk,nsurf,
     .                  irst,iflag,ivert,slavept,nslave,iskip,jskip,
     .                  kskip,nsegdfrm,idfrmseg,iaesurf,maxsegdg,
     .                  nmaster,iseq)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Compute mesh deformation via Transfinite Interpolation
c              using arc-length blending functions, compute grid speeds
c              at grid points, store of the current grid as the one at
c              the previous time step, and update the grid shape.
c
c     irst    = 0 standard update; do all steps described above
c               1 restart - compute the mesh deformation and update
c                 only the grid shape; do not update the grid speeds.
c
c     ismooth = 0 use only transfinite interpolation 
c             > 0 do ismooth smoothing steps using the spring analogy 
c                 and mesh shape preservation near the k=1 and k=kdim
c                 surfaces
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      integer stats
c
      dimension delti(jdim,kdim,3,2)
      dimension deltj(kdim,idim,3,2)
      dimension deltk(jdim,idim,3,2)
      dimension iaesurf(maxbl,maxsegdg)
      dimension ibcinfo(maxbl,maxseg,7,2)
      dimension icsf(maxbl,maxsegdg)
      dimension icsi(maxbl,maxsegdg)
      dimension idefrm(maxbl)
      dimension idfrmseg(maxbl,maxsegdg)
      dimension iskip(maxbl)
      dimension jbcinfo(maxbl,maxseg,7,2)
      dimension jcsf(maxbl,maxsegdg)
      dimension jcsi(maxbl,maxsegdg)
      dimension jskip(maxbl)
      dimension kbcinfo(maxbl,maxseg,7,2)
      dimension kcsf(maxbl,maxsegdg)
      dimension kcsi(maxbl,maxsegdg)
      dimension kskip(maxbl)
      dimension nbci0(maxbl)
      dimension nbcidim(maxbl)
      dimension nbcj0(maxbl)
      dimension nbcjdim(maxbl)
      dimension nbck0(maxbl)
      dimension nbckdim(maxbl)
      dimension nou(nbuf)
      dimension nsegdfrm(maxbl)
      dimension slavept(nslave,3,nmaster,5)
      dimension vel(jdim,kdim,idim,3)
      dimension wk(9*nsurf)
      dimension x(jdim,kdim,idim)
      dimension xnm1(jdim,kdim,idim)
      dimension xnm2(jdim,kdim,idim)
      dimension y(jdim,kdim,idim)
      dimension ynm1(jdim,kdim,idim)
      dimension ynm2(jdim,kdim,idim)
      dimension z(jdim,kdim,idim)
      dimension znm1(jdim,kdim,idim)
      dimension znm2(jdim,kdim,idim)
c
      allocatable :: arci(:,:,:)
      allocatable :: arcj(:,:,:)
      allocatable :: arck(:,:,:)
      allocatable :: dvol(:)
      allocatable :: dx(:,:,:)
      allocatable :: dx1(:,:,:)
      allocatable :: dx2(:,:,:)
      allocatable :: dx3(:,:,:)
      allocatable :: dy(:,:,:)
      allocatable :: dy1(:,:,:)
      allocatable :: dy2(:,:,:)
      allocatable :: dy3(:,:,:)
      allocatable :: dz(:,:,:)
      allocatable :: dz1(:,:,:)
      allocatable :: dz2(:,:,:)
      allocatable :: dz3(:,:,:)
      allocatable :: ibl(:)
      allocatable :: jend(:)
      allocatable :: jsta(:)
      allocatable :: xoi(:,:)
      allocatable :: xoj(:,:)
      allocatable :: xok(:,:)

      common /twod/ i2d
      common /deformz/ beta1,ismooth,negvol
c
c     allocate memory
c
      memuse = 0
c
      allocate( arci(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'arci',memuse,stats)
      allocate( arcj(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'arcj',memuse,stats)
      allocate( arck(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'arck',memuse,stats)
      allocate( dx(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx',memuse,stats)
      allocate( dx1(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx1',memuse,stats)
      allocate( dx2(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx2',memuse,stats)
      allocate( dx3(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx3',memuse,stats)
      allocate( dy(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy',memuse,stats)
      allocate( dy1(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy1',memuse,stats)
      allocate( dy2(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy2',memuse,stats)
      allocate( dy3(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy3',memuse,stats)
      allocate( dz(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz',memuse,stats)
      allocate( dz1(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz1',memuse,stats)
      allocate( dz2(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz2',memuse,stats)
      allocate( dz3(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz3',memuse,stats)
      allocate( ibl(nsurf), stat=stats )
      call umalloc(nsurf,1,'ibl',memuse,stats)
c
c     determine arc lengths for use in TFI blending functions
c
      call arclen(idim,jdim,kdim,arci,arcj,arck,xnm1,ynm1,
     .            znm1,nbl,nou,bou,nbuf,ibufdim,myid)
c
c     initialize deltas
c
      do i=1,idim
         do j=1,jdim
            do k=1,kdim
               dx(j,k,i) = 0.
               dy(j,k,i) = 0.
               dz(j,k,i) = 0.
            end do
         end do
      end do
c
c     determine the master surface point for the 4 corner
c     points of each block subface that is NOT on a solid
c     surface, and calculate the deltas for the corner
c     points based on a decaying function of distance.
c     (the x,y,z and dx,dy,dz of the master surface points
c     are stored in wk); then set deltas interior to such
c     subfaces via TFI 
c
      iskp = iskip(nbl)
      jskp = jskip(nbl)
      kskp = kskip(nbl)
c
c     i=1 subfaces
c
      do nseg=1,nbci0(nbl)
         ii     = 1
         ibctyp = ibcinfo(nbl,nseg,1,ii)
         if (abs(ibctyp).ne.2004 .and.
     .       abs(ibctyp).ne.2014 .and.
     .       abs(ibctyp).ne.1005 .and.
     .       abs(ibctyp).ne.1006 .or.
     .       real(beta1).lt.0.) then
             js = ibcinfo(nbl,nseg,2,ii)
             je = ibcinfo(nbl,nseg,3,ii)
             ks = ibcinfo(nbl,nseg,4,ii)
             ke = ibcinfo(nbl,nseg,5,ii)
             do j=js,je,jskp
                do k=ks,ke,kskp
                   call setcorner(j,k,1,xnm1,ynm1,znm1,dx,dy,
     .                            dz,jdim,kdim,idim,wk,nsurf,iflag,
     .                            ivert,slavept,nslave,nou,bou,
     .                            ibufdim,nbuf,myid,ibl,nmaster,
     .                            iseq)
                end do
             end do
             do j=js,je-jskp,jskp
                do k=ks,ke,kskp
                   call tfiedge(idim,jdim,kdim,dx,dy,dz,1,1,
     .                          j,j+jskp,k,k,arci,arcj,arck,
     .                          nou,bou,nbuf,ibufdim,myid)
                end do
             end do
             do j=js,je,jskp
                do k=ks,ke-kskp,kskp
                   call tfiedge(idim,jdim,kdim,dx,dy,dz,1,1,
     .                          j,j,k,k+kskp,arci,arcj,arck,
     .                          nou,bou,nbuf,ibufdim,myid)
                end do
             end do
c
c            preserve deltas on edges of solid surfaces
c
             if (idefrm(nbl).lt.999 .and. real(beta1).ge.0.) then
                if (js.eq.1 .or. je.eq.jdim .or.
     .              ks.eq.1 .or. ke.eq.kdim) then
                    call bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,
     .                           jcsf,kcsi,kcsf,icsi,icsf,jdim,kdim,
     .                           idim,maxbl,maxsegdg,nsegdfrm)
                end if
             end if
             do j=js,je-jskp,jskp
                do k=ks,ke-kskp,kskp
                   call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                          dx2,dy2,dz2,arci,arcj,arck,1,1,
     .                          j,j+jskp,k,k+kskp,nou,bou,nbuf,
     .                          ibufdim,myid)
                end do
             end do
         end if
      end do
c
c     i=idim subfaces
c
      do nseg=1,nbcidim(nbl)
         ii     = 2
         ibctyp = ibcinfo(nbl,nseg,1,ii)
         if (abs(ibctyp).ne.2004 .and.
     .       abs(ibctyp).ne.2014 .and.
     .       abs(ibctyp).ne.1005 .and.
     .       abs(ibctyp).ne.1006 .or.
     .       real(beta1).lt.0.) then
             js = ibcinfo(nbl,nseg,2,ii)
             je = ibcinfo(nbl,nseg,3,ii)
             ks = ibcinfo(nbl,nseg,4,ii)
             ke = ibcinfo(nbl,nseg,5,ii)
             do j=js,je,jskp
                do k=ks,ke,kskp
                   call setcorner(j,k,idim,xnm1,ynm1,znm1,dx,dy,
     .                            dz,jdim,kdim,idim,wk,nsurf,iflag,
     .                            ivert,slavept,nslave,nou,bou,
     .                            ibufdim,nbuf,myid,ibl,nmaster,
     .                            iseq)
                end do
             end do
             do j=js,je-jskp,jskp
                do k=ks,ke,kskp
                   call tfiedge(idim,jdim,kdim,dx,dy,dz,idim,idim,
     .                          j,j+jskp,k,k,arci,arcj,arck,
     .                          nou,bou,nbuf,ibufdim,myid)
                end do
             end do
             do j=js,je,jskp
                do k=ks,ke-kskp,kskp
                   call tfiedge(idim,jdim,kdim,dx,dy,dz,idim,idim,
     .                          j,j,k,k+kskp,arci,arcj,arck,
     .                          nou,bou,nbuf,ibufdim,myid)
                end do
             end do
c
c            preserve deltas on edges of solid surfaces
c
             if (idefrm(nbl).lt.999 .and. real(beta1).ge.0.) then
                if (js.eq.1 .or. je.eq.jdim .or.
     .              ks.eq.1 .or. ke.eq.kdim) then
                    call bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,
     .                           jcsf,kcsi,kcsf,icsi,icsf,jdim,kdim,
     .                           idim,maxbl,maxsegdg,nsegdfrm)
                end if
             end if
             do j=js,je-jskp,jskp
                do k=ks,ke-kskp,kskp
                   call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                          dx2,dy2,dz2,arci,arcj,arck,idim,idim,
     .                          j,j+jskp,k,k+kskp,nou,bou,nbuf,
     .                          ibufdim,myid)
                end do
             end do
         end if
      end do
c
      if (i2d .eq. 0) then
c
c        j=1 subfaces
c
         do nseg=1,nbcj0(nbl)
            jj     = 1
            jbctyp = jbcinfo(nbl,nseg,1,jj)
            if (abs(jbctyp).ne.2004 .and.
     .          abs(jbctyp).ne.2014 .and.
     .          abs(jbctyp).ne.1005 .and.
     .          abs(jbctyp).ne.1006 .or.
     .          real(beta1).lt.0.) then
                is = jbcinfo(nbl,nseg,2,jj)
                ie = jbcinfo(nbl,nseg,3,jj)
                ks = jbcinfo(nbl,nseg,4,jj)
                ke = jbcinfo(nbl,nseg,5,jj)
                do i=is,ie,iskp
                   do k=ks,ke,kskp
                      call setcorner(1,k,i,xnm1,ynm1,znm1,dx,dy,
     .                               dz,jdim,kdim,idim,wk,nsurf,iflag,
     .                               ivert,slavept,nslave,nou,bou,
     .                               ibufdim,nbuf,myid,ibl,nmaster,
     .                               iseq)
                   end do
                end do
                do i=is,ie-iskp,iskp
                   do k=ks,ke,kskp
                      call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i+iskp,
     .                             1,1,k,k,arci,arcj,arck,
     .                             nou,bou,nbuf,ibufdim,myid)
                   end do
                end do
                do i=is,ie,iskp
                   do k=ks,ke-kskp,kskp
                      call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i,
     .                             1,1,k,k+kskp,arci,arcj,arck,
     .                             nou,bou,nbuf,ibufdim,myid)
                   end do
                end do
c
c               preserve deltas on edges of solid surfaces
c
                if (idefrm(nbl).lt.999 .and. real(beta1).ge.0.) then
                   if (is.eq.1 .or. ie.eq.idim .or.
     .                 ks.eq.1 .or. ke.eq.kdim) then
                       call bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,
     .                              jcsf,kcsi,kcsf,icsi,icsf,jdim,kdim,
     .                              idim,maxbl,maxsegdg,nsegdfrm)
                   end if
                end if
                do i=is,ie-iskp,iskp
                   do k=ks,ke-kskp,kskp
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,i+iskp,
     .                             1,1,k,k+kskp,nou,bou,nbuf,
     .                             ibufdim,myid)
                   end do
                end do
            end if
         end do
c
c        j=jdim subfaces
c
         do nseg=1,nbcjdim(nbl)
            jj     = 2
            jbctyp = jbcinfo(nbl,nseg,1,jj)
            if (abs(jbctyp).ne.2004 .and.
     .          abs(jbctyp).ne.2014 .and.
     .          abs(jbctyp).ne.1005 .and.
     .          abs(jbctyp).ne.1006 .or.
     .          real(beta1).lt.0.) then
                is = jbcinfo(nbl,nseg,2,jj)
                ie = jbcinfo(nbl,nseg,3,jj)
                ks = jbcinfo(nbl,nseg,4,jj)
                ke = jbcinfo(nbl,nseg,5,jj)
                do i=is,ie,iskp
                   do k=ks,ke,kskp
                      call setcorner(jdim,k,i,xnm1,ynm1,znm1,dx,dy,
     .                               dz,jdim,kdim,idim,wk,nsurf,iflag,
     .                               ivert,slavept,nslave,nou,bou,
     .                               ibufdim,nbuf,myid,ibl,nmaster,
     .                               iseq)
                   end do
                end do
                do i=is,ie-iskp,iskp
                   do k=ks,ke,kskp
                      call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i+iskp,
     .                             jdim,jdim,k,k,arci,arcj,arck,
     .                             nou,bou,nbuf,ibufdim,myid)
                   end do
                end do
                do i=is,ie,iskp
                   do k=ks,ke-kskp,kskp
                      call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i,
     .                             jdim,jdim,k,k+kskp,arci,arcj,arck,
     .                             nou,bou,nbuf,ibufdim,myid)
                   end do
                end do
c
c               preserve deltas on edges of solid surfaces
c
                if (idefrm(nbl).lt.999 .and. real(beta1).ge.0.) then
                   if (is.eq.1 .or. ie.eq.idim .or.
     .                 ks.eq.1 .or. ke.eq.kdim) then
                       call bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,
     .                              jcsf,kcsi,kcsf,icsi,icsf,jdim,kdim,
     .                              idim,maxbl,maxsegdg,nsegdfrm)
                   end if
                end if
                do i=is,ie-iskp,iskp
                   do k=ks,ke-kskp,kskp
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,i+iskp,
     .                             jdim,jdim,k,k+kskp,nou,bou,nbuf,
     .                             ibufdim,myid)
                   end do
                end do
            end if
         end do
c
c        k=1 subfaces
c
         do nseg=1,nbck0(nbl)
            kk     = 1
            kbctyp = kbcinfo(nbl,nseg,1,kk)
            if (abs(kbctyp).ne.2004 .and.
     .          abs(kbctyp).ne.2014 .and.
     .          abs(kbctyp).ne.1005 .and.
     .          abs(kbctyp).ne.1006 .or.
     .          real(beta1).lt.0.) then
                is = kbcinfo(nbl,nseg,2,kk)
                ie = kbcinfo(nbl,nseg,3,kk)
                js = kbcinfo(nbl,nseg,4,kk)
                je = kbcinfo(nbl,nseg,5,kk)
                do i=is,ie,iskp
                   do j=js,je,jskp
                      call setcorner(j,1,i,xnm1,ynm1,znm1,dx,dy,
     .                               dz,jdim,kdim,idim,wk,nsurf,iflag,
     .                               ivert,slavept,nslave,nou,bou,
     .                               ibufdim,nbuf,myid,ibl,nmaster,
     .                               iseq)
                   end do
                end do
                do i=is,ie-iskp,iskp
                   do j=js,je,jskp
                      call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i+iskp,
     .                             j,j,1,1,arci,arcj,arck,
     .                             nou,bou,nbuf,ibufdim,myid)
                   end do
                end do
                do i=is,ie,iskp
                   do j=js,je-jskp,jskp
                      call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i,
     .                             j,j+jskp,1,1,arci,arcj,arck,
     .                             nou,bou,nbuf,ibufdim,myid)
                   end do
                end do
c
c               preserve deltas on edges of solid surfaces
c
                if (idefrm(nbl).lt.999 .and. real(beta1).ge.0.) then
                   if (is.eq.1 .or. ie.eq.idim .or.
     .                 js.eq.1 .or. je.eq.jdim) then
                       call bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,
     .                              jcsf,kcsi,kcsf,icsi,icsf,jdim,kdim,
     .                              idim,maxbl,maxsegdg,nsegdfrm)
                   end if
                end if
                do i=is,ie-iskp,iskp
                   do j=js,je-jskp,jskp
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,i+iskp,
     .                             j,j+jskp,1,1,nou,bou,nbuf,
     .                             ibufdim,myid)
                   end do
                end do
            end if
         end do
c
c        k=kdim subfaces
c
         do nseg=1,nbckdim(nbl)
            kk     = 2
            kbctyp = kbcinfo(nbl,nseg,1,kk)
            if (abs(kbctyp).ne.2004 .and.
     .          abs(kbctyp).ne.2014 .and.
     .          abs(kbctyp).ne.1005 .and.
     .          abs(kbctyp).ne.1006 .or.
     .          real(beta1).lt.0.) then
                is = kbcinfo(nbl,nseg,2,kk)
                ie = kbcinfo(nbl,nseg,3,kk)
                js = kbcinfo(nbl,nseg,4,kk)
                je = kbcinfo(nbl,nseg,5,kk)
                do i=is,ie,iskp
                   do j=js,je,jskp
                      call setcorner(j,kdim,i,xnm1,ynm1,znm1,dx,dy,
     .                               dz,jdim,kdim,idim,wk,nsurf,iflag,
     .                               ivert,slavept,nslave,nou,bou,
     .                               ibufdim,nbuf,myid,ibl,nmaster,
     .                               iseq)
                   end do
                end do
                do i=is,ie-iskp,iskp
                   do j=js,je,jskp
                      call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i+iskp,
     .                             j,j,kdim,kdim,arci,arcj,arck,
     .                             nou,bou,nbuf,ibufdim,myid)
                   end do
                end do
                do i=is,ie,iskp
                   do j=js,je-jskp,jskp
                      call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i,
     .                             j,j+jskp,kdim,kdim,arci,arcj,arck,
     .                             nou,bou,nbuf,ibufdim,myid)
                   end do
                end do
c
c               preserve deltas on edges of solid surfaces
c
                if (idefrm(nbl).lt.999 .and. real(beta1).ge.0.) then
                   if (is.eq.1 .or. ie.eq.idim .or.
     .                 js.eq.1 .or. je.eq.jdim) then
                       call bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,
     .                              jcsf,kcsi,kcsf,icsi,icsf,jdim,kdim,
     .                              idim,maxbl,maxsegdg,nsegdfrm)
                   end if
                end if
                do i=is,ie-iskp,iskp
                   do j=js,je-jskp,jskp
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,i+iskp,
     .                             j,j+jskp,kdim,kdim,nou,bou,nbuf,
     .                             ibufdim,myid)
                   end do
                end do
            end if
         end do
c
      end if
c
c     TFI to get deltas throughout the volume
c
      call tfivol(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,dx2,dy2,dz2,
     .            dx3,dy3,dz3,arci,arcj,arck,nou,bou,nbuf,ibufdim,
     .            myid)
c
      if (ismooth.gt.0) then
c
c        This step incorporates the correct translation and rota-
c        tion of the interior mesh near the k=1 and k=kdim boundaries,
c        exponentially decaying to a spring analogy smoothing away
c        from the boundary.
c
c        allocate memory for the smoothing step
c
         allocate( dvol(jdim*idim*kdim), stat=stats )
         call umalloc(jdim*idim*kdim,0,'dvol',memuse,stats)
         allocate( xoi(jdim*kdim,6), stat=stats )
         call umalloc(jdim*kdim*6,0,'xoi',memuse,stats)
         allocate( xoj(kdim*idim,6), stat=stats )
         call umalloc(kdim*idim*6,0,'xoj',memuse,stats)
         allocate( xok(jdim*idim,6), stat=stats )
         call umalloc(jdim*idim*6,0,'xok',memuse,stats)
         allocate( jend(2), stat=stats )
         call umalloc(2,1,'jend',memuse,stats)
         allocate( jsta(2), stat=stats )
         call umalloc(2,1,'jsta',memuse,stats)
c
         call dmsetup0(jdim,kdim,idim,icsi,icsf,jcsi,jcsf,
     .                 kcsi,kcsf,nbl,maxbl,maxsegdg,
     .                 jsta,jend,x,y,z,dx,dy,dz,xoj,xoi,xok,dvol)
c
         call gdsa(time,xmach,dt,dx,dy,dz,x,y,z,xnm1,
     .             ynm1,znm1,dvol,xok,xoi,xoj,irotat,ndeftp1,
     .             iupdat,nbl,nt,idim,jdim,kdim,iaesol,jsta,jend,
     .             ismooth)
c
c        release memory
c
         deallocate(dvol)
         deallocate(xok)
         deallocate(xoi)
         deallocate(xoj)
         deallocate(jend)
         deallocate(jsta)
c
      end if
c
      if (irst.eq.0) then
c
c        compute grid velocites before overwriting current grid
c
         if (real(time).le.real(dt) .or. abs(ita).eq.1) then
c
c           first order backward derivatives
c
            fact = 1./dt
            do i=1,idim
               do k=1,kdim
                  do j=1,jdim
                     vel(j,k,i,1) = fact*((dx(j,k,i)+x(j,k,i))
     .                            - xnm1(j,k,i))
                     vel(j,k,i,2) = fact*((dy(j,k,i)+y(j,k,i))
     .                            - ynm1(j,k,i))
                     vel(j,k,i,3) = fact*((dz(j,k,i)+z(j,k,i))
     .                            - znm1(j,k,i))
                  end do
               end do
            end do
c
         else
c
c           second order backward derivatives
c
            fact = 1./(2.*dt)
            do i=1,idim
               do k=1,kdim
                  do j=1,jdim
                     vel(j,k,i,1) = fact*(3.*(dx(j,k,i)+x(j,k,i))
     .                            - 4*xnm1(j,k,i)+xnm2(j,k,i))
                     vel(j,k,i,2) = fact*(3.*(dy(j,k,i)+y(j,k,i))
     .                            - 4*ynm1(j,k,i)+ynm2(j,k,i))
                     vel(j,k,i,3) = fact*(3.*(dz(j,k,i)+z(j,k,i))
     .                            - 4*znm1(j,k,i)+znm2(j,k,i))
                  end do
               end do
            end do
c
         end if
c
      end if
c
c     for second order case, store current--->old before
c     updating grid
c
      if (abs(ita).gt.1) then
         do i=1,idim
            do k=1,kdim
               do j=1,jdim
                  xnm2(j,k,i) = xnm1(j,k,i)
                  ynm2(j,k,i) = ynm1(j,k,i)
                  znm2(j,k,i) = znm1(j,k,i)
               end do
            end do
         end do
      end if
c
c     add deltas to the current mesh to get the new one
c
      do i=1,idim
         do k=1,kdim
            do j=1,jdim
               x(j,k,i) = x(j,k,i) + dx(j,k,i)
               y(j,k,i) = y(j,k,i) + dy(j,k,i)
               z(j,k,i) = z(j,k,i) + dz(j,k,i)
            end do
         end do
      end do
c
c     release memory
c
      deallocate(dx1)
      deallocate(dy1)
      deallocate(dz1)
      deallocate(dx2)
      deallocate(dy2)
      deallocate(dz2)
      deallocate(dx3)
      deallocate(dy3)
      deallocate(dz3)
      deallocate(dx)
      deallocate(dy)
      deallocate(dz)
      deallocate(arci)
      deallocate(arcj)
      deallocate(arck)
      deallocate(ibl)
c
      return
      end
       subroutine dmsetup(jdm,kdm,idm,nbl,jsta,jend,ipid,imid,rj1jd,
     .                    rj2jd,epsilon) 
c***********************************************************************
c     Purpose:  Sets up information for use in the smoothing step
c
c
c     Bob Bartels NASA Langley Aeroelasticity Branch    (757)864-2813
c       responsible for the following subroutines:
c         gdsa      -      spring analogy smoothing and near surface
c                          orientation preservation.
c         dmsetup0  -      some initialization of xoi,xoj and dx,dy,dz 
c         dmsetup1  -      some initialization of xok and dx,dy,dz 
c         dmsetup   -      some initialization of arrays 
c         bdata     -      computation of block boundary face spacing
c                          arrays 
c         permut_ijkxyz -  permutes indices i and j to better align 
c                          k = 0 and k = kdim surfaces for vectorizing 
c         unperm_ijkxyz -  unpermutes indices i and j after smoothing 
c     
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
       common /grid1/ algepmn,coef(2)   
       common /grid2/ kimx(2),jte1(2),         
     .   jflp(2,2),icsi(2),icso(2),iafmv1(2),jsp1(2),  
     .   jsp2(2),jsp3(2),itecs(2),ispcs(2),itp(2,2)
     .  ,iaxrt(2)
       common /grid4/ time0
       common /sklton/ isklton
       dimension epsilon(kdm,2)
       dimension  ipid(jdm*idm),imid(jdm*idm),rj1jd(jdm*idm),
     .            rj2jd(jdm*idm),jend(2),jsta(2)


       do 30 k = 1,kdm
         epsilon(k,1) = 1. 
         epsilon(k,2) = 1. 
30     continue

       coefa  =  -algepmn
       do 29 n = 1,2
       do 29 k = kimx(n)+1,kdm 
          ro1 = coef(n)*real(k-kimx(n)) 
          arg1= coefa*ro1        
          epsilon(k,n) = exp(-arg1) 
29     continue 
       if(idm.gt.2) then
        do 438 j = 1,jdm*idm  
          ipid(j)= j+jdm           
          imid(j)= j-jdm  
438     continue  
        jdmi1 = jdm*(idm-1)
        do 440 j = 1,jdm 
          imid(j)= j   
          ipid(j+jdmi1)=j+jdmi1  
440     continue 
       else
        do 448 j = 1,jdm     
          imid(j) = j    
          ipid(j) = j+jdm      
448     continue  
        do 449 j = 1+jdm,jdm*idm     
          imid(j) = j-jdm      
          ipid(j) = j      
449     continue  
       end if 

       do 511 i=1,idm   
       jdmi1 = jdm*(i-1)  
       do 511 j = 1,jdm   
           rj1jd(j+jdmi1) = real(jdm-j) 
           rj2jd(j+jdmi1) = real(j-1) 
511    continue 

       return
       end  
       subroutine dmsetup0(jdm,kdm,idm,itp1,itp2,jtp1,jtp2,
     .                 ksta,kend,nbl,maxbl,maxsegdg,
     .                 jsta,jend,x,y,z,dx,dy,dz,xoj,xoi,xok,dvol)

#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      integer stats
c
      dimension dvol(jdm*idm*kdm)
      dimension dx(jdm,kdm,idm)
      dimension dy(jdm,kdm,idm)
      dimension dz(jdm,kdm,idm)
      dimension itp1(maxbl,maxsegdg)
      dimension itp2(maxbl,maxsegdg)
      dimension jend(2)
      dimension jsta(2)
      dimension jtp1(maxbl,maxsegdg)
      dimension jtp2(maxbl,maxsegdg)
      dimension kend(maxbl,maxsegdg)
      dimension ksta(maxbl,maxsegdg)
      dimension x(jdm,kdm,idm)
      dimension xoi(jdm*kdm,6)
      dimension xoj(kdm*idm,6)
      dimension xok(jdm*idm,6)
      dimension y(jdm,kdm,idm)
      dimension z(jdm,kdm,idm)

      allocatable :: dxn(:,:,:)
      allocatable :: dyn(:,:,:)
      allocatable :: dzn(:,:,:)
      allocatable :: epsilon3(:)
      allocatable :: epsilon4(:)

      memuse = 0
c
      allocate( dxn(jdm,kdm,idm), stat=stats )
      call umalloc(jdm*kdm*idm,0,'dxn',memuse,stats)
      allocate( dyn(jdm,kdm,idm), stat=stats )
      call umalloc(jdm*kdm*idm,0,'dyn',memuse,stats)
      allocate( dzn(jdm,kdm,idm), stat=stats )
      call umalloc(jdm*kdm*idm,0,'dzn',memuse,stats)
      allocate( epsilon3(jdm), stat=stats )
      call umalloc(jdm,0,'epsilon3',memuse,stats)
      allocate( epsilon4(idm), stat=stats )
      call umalloc(idm,0,'epsilon4',memuse,stats)
c
        if(idm.gt.2) then
          do j = 1,kdm*jdm
             xoi(j,1)       = x(j,1,1)       + dx(j,1,1)
             xoi(j,2)       = y(j,1,1)       + dy(j,1,1)
             xoi(j,3)       = z(j,1,1)       + dz(j,1,1)
             xoi(j,4)       = x(j,1,idm)    + dx(j,1,idm)
             xoi(j,5)       = y(j,1,idm)    + dy(j,1,idm)
             xoi(j,6)       = z(j,1,idm)    + dz(j,1,idm)
          end do
        end if
        if(jdm.gt.2) then
         do i =1,idm
          do k =1,kdm
           xoj(k+kdm*(i-1),1) = x(1,k,i)       + dx(1,k,i)
           xoj(k+kdm*(i-1),2) = y(1,k,i)       + dy(1,k,i)
           xoj(k+kdm*(i-1),3) = z(1,k,i)       + dz(1,k,i)
           xoj(k+kdm*(i-1),4) = x(jdm,k,i)    + dx(jdm,k,i)
           xoj(k+kdm*(i-1),5) = y(jdm,k,i)    + dy(jdm,k,i)
           xoj(k+kdm*(i-1),6) = z(jdm,k,i)    + dz(jdm,k,i)
          end do
         end do

         do 303 j = 1,jdm-1 
          epsilon3(j)=sqrt((x(j+1,1,1)-x(j,1,1))
     .                  *(x(j+1,1,1)-x(j,1,1))
     .                  +(y(j+1,1,1)-y(j,1,1))
     .                  *(y(j+1,1,1)-y(j,1,1))
     .                  +(z(j+1,1,1)-z(j,1,1))
     .                  *(z(j+1,1,1)-z(j,1,1))) 
303      continue
         reps3         = 0.
         do 304 j = 1,jdm-1 
304      reps3       = reps3 + epsilon3(j)
         reps31      = 0.
         do 305 j = jdm-1,2,-1 
          reps31      = reps31 + epsilon3(j)
          epsilon3(j) = reps31/reps3
305      continue
         epsilon3(1)  = 1. 
         epsilon3(jdm)= 0.


         do 306 j = 1,jdm 
          eps3   = epsilon3(j)  
          omeps3 = 1. - eps3  
          do 306 i = 1,idm
          do 306 k = 1,kdm
           dxn(j,k,i) = eps3*(xoj(k+kdm*(i-1),1)-x(1,k,i))
     .              + omeps3*(xoj(k+kdm*(i-1),4)-x(jdm,k,i))  
           dyn(j,k,i) = eps3*(xoj(k+kdm*(i-1),2)-y(1,k,i))
     .              + omeps3*(xoj(k+kdm*(i-1),5)-y(jdm,k,i))  
           dzn(j,k,i) = eps3*(xoj(k+kdm*(i-1),3)-z(1,k,i))
     .              + omeps3*(xoj(k+kdm*(i-1),6)-z(jdm,k,i))  
 306     continue
       else
         do j = 1,jdm
            do i = 1,idm
               do k = 1,kdm
                  dxn(j,k,i) =0.
                  dyn(j,k,i) =0.
                  dzn(j,k,i) =0.
               end do
            end do
         end do
       end if
       if(idm.gt.2) then
         do 307 i = 1,idm-1 
          epsilon4(i)=sqrt((x(1,1,i+1)-x(1,1,i))
     .                    *(x(1,1,i+1)-x(1,1,i))
     .                    +(y(1,1,i+1)-y(1,1,i))
     .                    *(y(1,1,i+1)-y(1,1,i))
     .                    +(z(1,1,i+1)-z(1,1,i))
     .                    *(z(1,1,i+1)-z(1,1,i))) 
307      continue
         reps4         = 0.
         do 308 i = 1,idm-1 
308      reps4       = reps4 + epsilon4(i)
         reps41      = 0.
         do 309 i = idm-1,2,-1 
           reps41      = reps41 + epsilon4(i)
           epsilon4(i) = reps41/reps4
309      continue
         epsilon4(1)  = 1. 
         epsilon4(idm)= 0.

         do 310 i = 2,idm 
           eps    = epsilon4(i)  
           omeps  = 1. - eps   
           do 310 j = 1,kdm*jdm
            dxn(j,1,i)=dxn(j,1,i)+eps*(xoi(j,1)-x(j,1,1)  -dxn(j,1,1))
     .                        + omeps*(xoi(j,4)-x(j,1,idm)-dxn(j,1,idm))  
            dyn(j,1,i)=dyn(j,1,i)+eps*(xoi(j,2)-y(j,1,1)  -dyn(j,1,1))
     .                        + omeps*(xoi(j,5)-y(j,1,idm)-dyn(j,1,idm))  
            dzn(j,1,i)=dzn(j,1,i)+eps*(xoi(j,3)-z(j,1,1)  -dzn(j,1,1))
     .                        + omeps*(xoi(j,6)-z(j,1,idm)-dzn(j,1,idm))  
310      continue
         do j = 1,kdm*jdm
           dxn(j,1,1) = xoi(j,1)-x(j,1,1)
           dyn(j,1,1) = xoi(j,2)-y(j,1,1)
           dzn(j,1,1) = xoi(j,3)-z(j,1,1)
         end do
       end if

       deallocate(epsilon3)
       deallocate(epsilon4)

       call permut_ijkxyz(jdm,idm,kdm,nbl,dx,dy,dz,dvol)
       call permut_ijkxyz(jdm,idm,kdm,nbl,dxn,dyn,dzn,dvol)
       call permut_ijkxyz(jdm,idm,kdm,nbl,x,y,z,dvol)

       call dmsetup1(jdm,kdm,idm,itp1,itp2,jtp1,jtp2,
     .                 ksta,kend,nbl,maxbl,maxsegdg,
     .                 jsta,jend,x,y,z,dx,dy,dz,dxn,
     .                 dyn,dzn,xok)

       deallocate(dxn)
       deallocate(dyn)
       deallocate(dzn)
       
       return
       end
       subroutine dmsetup1(jdm,kdm,idm,itp1,itp2,jtp1,jtp2,
     .                 ksta,kend,nbl,maxbl,maxsegdg,
     .                 jsta,jend,x,y,z,dx,dy,dz,dxn,
     .                 dyn,dzn,xok)
c***********************************************************************
c     Purpose:  Sets up information for use in the smoothing step
c
c
c     Bob Bartels NASA Langley Aeroelasticity Branch    (757)864-2813
c       responsible for the following subroutines:
c         gdsa      -      spring analogy smoothing and near surface
c                          orientation preservation.
c         dmsetup0  -      some initialization of xoi,xoj and dx,dy,dz 
c         dmsetup1  -      some initialization of xok and dx,dy,dz 
c         dmsetup   -      some initialization of arrays 
c         bdata     -      computation of block boundary face spacing
c                          arrays 
c         permut_ijkxyz -  permutes indices i and j to better align 
c                          k = 0 and k = kdim surfaces for vectorizing 
c         unperm_ijkxyz -  unpermutes indices i and j after smoothing 
c     
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
       common /grid1/ algepmn,coef(2)
       common /grid2/ kimx(2),jte1(2),         
     .   jflp(2,2),icsi(2),icso(2),iafmv1(2),jsp1(2),  
     .   jsp2(2),jsp3(2),itecs(2),ispcs(2),itp(2,2)
     .  ,iaxrt(2)
       common /grid15/ igtyp(2),idmg(2),jdmg(2),kdmg(2)
       common /sklton/ isklton
c
      integer stats
c
       dimension dx(jdm*idm*kdm),dy(jdm*idm*kdm),dz(jdm*idm*kdm)
       dimension dxn(jdm*idm*kdm),dyn(jdm*idm*kdm),dzn(jdm*idm*kdm)
       dimension x(jdm*idm*kdm),y(jdm*idm*kdm),z(jdm*idm*kdm)
       dimension xok(jdm*idm,6)
       dimension itp1(maxbl,maxsegdg),itp2(maxbl,maxsegdg),
     .           jtp1(maxbl,maxsegdg),jtp2(maxbl,maxsegdg),
     .           ksta(maxbl,maxsegdg),kend(maxbl,maxsegdg)
     .          ,jend(2),jsta(2)


      allocatable :: epsilon3(:)

      memuse = 0
c
      allocate( epsilon3(kdm), stat=stats )
      call umalloc(kdm,0,'epsilon3',memuse,stats)
c
         jsta(1)= jtp1(nbl,1)
         jend(1)= jtp2(nbl,1)
         itp(1,1) = itp1(nbl,1)
         itp(2,1) = itp2(nbl,1)

         igtyp(1)=1
         idmg(1) = idm
         jdmg(1) = jdm
         kdmg(1) = kdm

        if(jdmg(1)-jend(1)+1.ne.jsta(1)) then
          jsta(1)= 1
          jend(1)= jdmg(1)
        end if 


        jdid = jdm*idm
        jk2 = jdm*idm*(kdm-1)  
        do 1013 j = 1,jdid 
              xok(j,1) = x(j)     + dx(j) 
              xok(j,2) = y(j)     + dy(j) 
              xok(j,3) = z(j)     + dz(j) 
              xok(j,4) = x(j+jk2) + dx(j+jk2)
              xok(j,5) = y(j+jk2) + dy(j+jk2)
              xok(j,6) = z(j+jk2) + dz(j+jk2)
1010        continue  
1013    continue
       rjdm1 = 1.e+6
       rjdm2 = 1.e+6
       do j = 1,jdm-1 
         rjdm3= ((x(j+1)- x(j)) *(x(j+1)- x(j))
     .          +(y(j+1)- y(j)) *(y(j+1)- y(j))
     .          +(z(j+1)- z(j)) *(z(j+1)- z(j)))
         rjdm4= ((x(j+1+jdm*idm*(kdm-1))- x(j+jdm*idm*(kdm-1))) 
     .          *(x(j+1+jdm*idm*(kdm-1))- x(j+jdm*idm*(kdm-1)))
     .          +(y(j+1+jdm*idm*(kdm-1))- y(j+jdm*idm*(kdm-1)))
     .          *(y(j+1+jdm*idm*(kdm-1))- y(j+jdm*idm*(kdm-1)))
     .          +(z(j+1+jdm*idm*(kdm-1))- z(j+jdm*idm*(kdm-1)))
     .          *(z(j+1+jdm*idm*(kdm-1))- z(j+jdm*idm*(kdm-1))))
         if(real(rjdm1).gt.real(rjdm3)) then
             rjdm1 = rjdm3
             jrjdm1= j
         end if
         if(real(rjdm2).gt.real(rjdm4)) then
             rjdm2 = rjdm4
             jrjdm2= j
         end if
       enddo
       j    = jrjdm1
       rkdm1= ((x(j+jdm*idm)- x(j))
     .        *(x(j+jdm*idm)- x(j))
     .        +(y(j+jdm*idm)- y(j))
     .        *(y(j+jdm*idm)- y(j))
     .        +(z(j+jdm*idm)- z(j))
     .        *(z(j+jdm*idm)- z(j)))
       j    = jrjdm2
       rkdm2= ((x(j+jdm*idm*(kdm-1))-
     .          x(j+jdm*idm*(kdm-2)))
     .        *(x(j+jdm*idm*(kdm-1))-
     .          x(j+jdm*idm*(kdm-2)))
     .        +(y(j+jdm*idm*(kdm-1))-
     .          y(j+jdm*idm*(kdm-2)))
     .        *(y(j+jdm*idm*(kdm-1))-
     .          y(j+jdm*idm*(kdm-2)))
     .        +(z(j+jdm*idm*(kdm-1))-
     .          z(j+jdm*idm*(kdm-2)))
     .        *(z(j+jdm*idm*(kdm-1))-
     .          z(j+jdm*idm*(kdm-2))))
       rkdm = ((x(jsta(1)+jdm*idm*(kdm-1))- x(jsta(1)))
     .        *(x(jsta(1)+jdm*idm*(kdm-1))- x(jsta(1)))
     .        +(y(jsta(1)+jdm*idm*(kdm-1))- y(jsta(1)))
     .        *(y(jsta(1)+jdm*idm*(kdm-1))- y(jsta(1)))
     .        +(z(jsta(1)+jdm*idm*(kdm-1))- z(jsta(1)))
     .        *(z(jsta(1)+jdm*idm*(kdm-1))- z(jsta(1))))


   
       if((real(rkdm1).lt.real(.0004*rkdm)).and.(real(rjdm1).gt.
     .                                      real(2.*rkdm1))) then 
        if(kdm.lt.15) then
          kimx(1) =1
          coef(1) = .016/real(kdm)
        else if(kdm.lt.20) then
          kimx(1) =2 
          coef(1) = .014/real(kdm)
        else if(kdm.lt.30) then
          kimx(1) =3
          coef(1) = .010/real(kdm)
        else
          kimx(1) =3 
          coef(1) = .006/real(kdm)
        end if
       else if(real(rjdm1).lt.real(rkdm1)) then 
        kimx(1) =0  
        coef(1) = .028/real(kdm)
       else 
        kimx(1) =1 
        coef(1) = .018/real(kdm)
       end if
       if((real(rkdm2).lt.real(.0004*rkdm)).and.(real(rjdm2).gt.
     .                                      real(2.*rkdm2))) then 
        if(kdm.lt.15) then
          kimx(2) =1
          coef(2) = .016/real(kdm)
        else if(kdm.lt.20) then
          kimx(2) =2
          coef(2) = .014/real(kdm) 
        else if(kdm.lt.30) then
          kimx(2) =3
          coef(2) = .010/real(kdm)
        else
          kimx(2) =3
          coef(2) = .006/real(kdm) 
        end if
       else if(real(rjdm2).lt.real(rkdm2)) then 
        kimx(2) =0 
        coef(2) = .058/real(kdm)
       else
        kimx(2) =1
        coef(2) = .018/real(kdm)
       end if


        algepmn = -2000. 

       do k = 1,kimx(1)
         epsilon3(k) = 1. 
       enddo
       do k = kdm,kdm-kimx(2),-1 
         epsilon3(k) = 0. 
       enddo
       rkdm = sqrt((x(jsta(1)+jdm*idm*kimx(1))-
     .              x(jsta(1)+jdm*idm*(kdm-kimx(2)-1)))
     .            *(x(jsta(1)+jdm*idm*kimx(1))-
     .              x(jsta(1)+jdm*idm*(kdm-kimx(2)-1)))
     .            +(y(jsta(1)+jdm*idm*kimx(1))-
     .              y(jsta(1)+jdm*idm*(kdm-kimx(2)-1)))
     .            *(y(jsta(1)+jdm*idm*kimx(1))-
     .              y(jsta(1)+jdm*idm*(kdm-kimx(2)-1)))
     .            +(z(jsta(1)+jdm*idm*kimx(1))-
     .              z(jsta(1)+jdm*idm*(kdm-kimx(2)-1)))
     .            *(z(jsta(1)+jdm*idm*kimx(1))-
     .              z(jsta(1)+jdm*idm*(kdm-kimx(2)-1))))
       do 28 k = kimx(1)+1,kdm-kimx(2)       
         epsilon3(k)= sqrt((x(jsta(1)+jdm*idm*(k-1))-
     .                      x(jsta(1)+jdm*idm*(kdm-kimx(2)-1)))
     .                    *(x(jsta(1)+jdm*idm*(k-1))-
     .                      x(jsta(1)+jdm*idm*(kdm-kimx(2)-1)))
     .                    +(y(jsta(1)+jdm*idm*(k-1))-
     .                      y(jsta(1)+jdm*idm*(kdm-kimx(2)-1)))
     .                    *(y(jsta(1)+jdm*idm*(k-1))-
     .                      y(jsta(1)+jdm*idm*(kdm-kimx(2)-1)))
     .                    +(z(jsta(1)+jdm*idm*(k-1))-
     .                      z(jsta(1)+jdm*idm*(kdm-kimx(2)-1)))
     .                    *(z(jsta(1)+jdm*idm*(k-1))-
     .                      z(jsta(1)+jdm*idm*(kdm-kimx(2)-1))))
28     continue
       do 29 k = kimx(1)+1,kdm-kimx(2)    
         epsilon3(k)=epsilon3(k)/rkdm
29     continue

       jdid = jdm*idm  
       jk2  = jdm*idm*(kdm-1)  
       do 300 k = 2,kdm
         eps    = epsilon3(k)
         omeps  = 1.-eps 
         jdidk  = jdid*(k-1)
         do j  = 1,jdm*idm                    
           dxn(j+jdidk)=dxn(j+jdidk)+eps*(xok(j,1)-x(j)    -dxn(j))
     .                            +omeps*(xok(j,4)-x(j+jk2)-dxn(j+jk2))
           dyn(j+jdidk)=dyn(j+jdidk)+eps*(xok(j,2)-y(j)    -dyn(j))
     .                            +omeps*(xok(j,5)-y(j+jk2)-dyn(j+jk2))
           dzn(j+jdidk)=dzn(j+jdidk)+eps*(xok(j,3)-z(j)    -dzn(j))
     .                            +omeps*(xok(j,6)-z(j+jk2)-dzn(j+jk2))
         end do
300    continue 
       do j  = 1,jdm*idm                    
         dxn(j)= xok(j,1)-x(j) 
         dyn(j)= xok(j,2)-y(j)
         dzn(j)= xok(j,3)-z(j)
       end do

      do j = 1,jdm*idm*kdm
         dx(j) = dxn(j)
         dy(j) = dyn(j)
         dz(j) = dzn(j)
      end do
    
      deallocate(epsilon3)

      return
      end 
       subroutine bdata(jdm,kdm,idm,nbl,x,y,z,xj1p,xj1m,xj2p,xj2m
     .                 ,xi1p,xi1m,xi2p,xi2m,xk1p,xk1m,xk2p,xk2m 
     .                 ,xj3p,seta)   
c***********************************************************************
c     Purpose:  Sets up boundary spacing arrays for use in the smoothing step
c
c
c     Bob Bartels NASA Langley Aeroelasticity Branch    (757)864-2813
c       responsible for the following subroutines:
c         gdsa      -      spring analogy smoothing and near surface
c                          orientation preservation.
c         dmsetup0  -      some initialization of xoi,xoj and dx,dy,dz 
c         dmsetup1  -      some initialization of xok and dx,dy,dz 
c         dmsetup   -      some initialization of arrays 
c         bdata     -      computation of block boundary face spacing
c                          arrays 
c         permut_ijkxyz -  permutes indices i and j to better align 
c                          k = 0 and k = kdim surfaces for vectorizing 
c         unperm_ijkxyz -  unpermutes indices i and j after smoothing 
c     
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
       dimension x(jdm*kdm*idm),y(jdm*kdm*idm),z(jdm*kdm*idm) 
       dimension xj1p(jdm*idm),xk1p(kdm), xj2p(jdm*idm)                          
     .          ,xk2p(kdm),xj1m(jdm*idm),xk1m(kdm),
     .           xj2m(jdm*idm), xk2m(kdm),xj3p(jdm*idm,2)  
     .          ,xi1m(jdm*idm),xi1p(jdm*idm),xi2m(jdm*idm),
     .           xi2p(jdm*idm)
       dimension seta(kdm) 

        jdid   = jdm*idm 
        jdidkpo= 1+jdid*(kdm-1)
        jdidkmo=-1+jdid*(kdm-1)
        jdidkd = jdid*(kdm-1)
        do 2 j = 2,idm*jdm-1
           xj1p(j)=(x(j+1)-x(j))*(x(j+1)-x(j))
     1            +(y(j+1)-y(j))*(y(j+1)-y(j))  
     2            +(z(j+1)-z(j))*(z(j+1)-z(j))  
           xj1m(j)=(x(j)-x(j-1))*(x(j)-x(j-1))
     1            +(y(j)-y(j-1))*(y(j)-y(j-1))  
     2            +(z(j)-z(j-1))*(z(j)-z(j-1))  
           xj2p(j)=(x(j+jdidkpo)-x(j+jdidkd))*
     .             (x(j+jdidkpo)-x(j+jdidkd))
     .            +(y(j+jdidkpo)-y(j+jdidkd))*
     .             (y(j+jdidkpo)-y(j+jdidkd))  
     .            +(z(j+jdidkpo)-z(j+jdidkd))*
     .             (z(j+jdidkpo)-z(j+jdidkd))  
           xj2m(j)=(x(j+jdidkd)-x(j+jdidkmo))*
     .             (x(j+jdidkd)-x(j+jdidkmo))
     .            +(y(j+jdidkd)-y(j+jdidkmo))* 
     .             (y(j+jdidkd)-y(j+jdidkmo))  
     .            +(z(j+jdidkd)-z(j+jdidkmo))*
     .             (z(j+jdidkd)-z(j+jdidkmo))  
2       continue
        do 3 j = 2,idm*jdm-1
           xj2m(j) = 1./sqrt(xj2m(j))  
           xj2p(j) = 1./sqrt(xj2p(j))  
           xj1m(j) = 1./sqrt(xj1m(j))  
           xj1p(j) = 1./sqrt(xj1p(j))  
 3      continue 
        nb0 = 0
        nb1 =  1
        nb2 =  2
        do 4 i = 1,idm
         xj1p(nb1+jdm*(i-1)) = xj1p(nb2+jdm*(i-1))  
         xj1m(nb1+jdm*(i-1)) = xj1p(nb2+jdm*(i-1))  
         xj2p(nb1+jdm*(i-1)) = xj1p(nb2+jdm*(i-1))  
         xj2m(nb1+jdm*(i-1)) = xj1p(nb2+jdm*(i-1))  
         xj1p(jdm*i)         = xj1p(jdm*i-1)  
         xj1m(jdm*i)         = xj1p(jdm*i-1)  
         xj2p(jdm*i)         = xj1p(jdm*i-1)  
         xj2m(jdm*i)         = xj1p(jdm*i-1)  
4       continue 

        jdidkd = jdid*(kdm-1)
        do 8 i = 2,idm-1
        do 7 j = 2,jdm-1
        xi1p(j+jdm*(i-1))=(x(j+i*jdm)-x(j+jdm*(i-1)))*
     .                    (x(j+i*jdm)-x(j+jdm*(i-1)))
     .                   +(y(j+i*jdm)-y(j+jdm*(i-1)))*
     .                    (y(j+i*jdm)-y(j+jdm*(i-1)))  
     .                   +(z(j+i*jdm)-z(j+jdm*(i-1)))*
     .                    (z(j+i*jdm)-z(j+jdm*(i-1)))  
        xi1m(j+jdm*(i-1))=(x(j+jdm*(i-1))-x(j+jdm*(i-2)))*
     .                    (x(j+jdm*(i-1))-x(j+jdm*(i-2)))
     .                   +(y(j+jdm*(i-1))-y(j+jdm*(i-2)))*
     .                    (y(j+jdm*(i-1))-y(j+jdm*(i-2)))  
     .                   +(z(j+jdm*(i-1))-z(j+jdm*(i-2)))*
     .                    (z(j+jdm*(i-1))-z(j+jdm*(i-2)))  
        xi2p(j+jdm*(i-1))=(x(j+jdm*i+jdidkd)-
     .                     x(j+jdm*(i-1)+jdidkd))*
     .                    (x(j+jdm*i+jdidkd)-
     .                     x(j+jdm*(i-1)+jdidkd))
     .                   +(y(j+jdm*i+jdidkd)-
     .                     y(j+jdm*(i-1)+jdidkd))*
     .                    (y(j+jdm*i+jdidkd)-
     .                     y(j+jdm*(i-1)+jdidkd))  
     .                   +(z(j+jdm*i+jdidkd)-
     .                     z(j+jdm*(i-1)+jdidkd))*
     .                    (z(j+jdm*i+jdidkd)-
     .                     z(j+jdm*(i-1)+jdidkd))  
        xi2m(j+jdm*(i-1))=(x(j+jdm*(i-1)+jdidkd)-
     .                     x(j+jdm*(i-2)+jdidkd))*
     .                    (x(j+jdm*(i-1)+jdidkd)-
     .                     x(j+jdm*(i-2)+jdidkd))
     .                   +(y(j+jdm*(i-1)+jdidkd)-
     .                     y(j+jdm*(i-2)+jdidkd))* 
     .                    (y(j+jdm*(i-1)+jdidkd)-
     .                     y(j+jdm*(i-2)+jdidkd))  
     .                   +(z(j+jdm*(i-1)+jdidkd)-
     .                     z(j+jdm*(i-2)+jdidkd))*
     .                    (z(j+jdm*(i-1)+jdidkd)-
     .                     z(j+jdm*(i-2)+jdidkd))  
           xi2tp = xi2m(j+jdm*(i-1)) 
           xi2m(j+jdm*(i-1)) = 1./sqrt(xi2tp)   
           xi2tp = xi2p(j+jdm*(i-1)) 
           xi2p(j+jdm*(i-1)) = 1./sqrt(xi2tp)  
           xi2tp = xi1m(j+jdm*(i-1)) 
           xi1m(j+jdm*(i-1)) = 1./sqrt(xi2tp)  
           xi2tp = xi1p(j+jdm*(i-1)) 
           xi1p(j+jdm*(i-1)) = 1./sqrt(xi2tp)  
 7      continue 
        xi2m(jdm*i)=xi2m(jdm*i-1) 
        xi2p(jdm*i)=xi2p(jdm*i-1) 
        xi1m(jdm*i)=xi1m(jdm*i-1) 
        xi1p(jdm*i)=xi1p(jdm*i-1) 
        xi2m(1+jdm*(i-1))=xi2m(2+jdm*(i-1)) 
        xi2p(1+jdm*(i-1))=xi2p(2+jdm*(i-1)) 
        xi1m(1+jdm*(i-1))=xi1m(2+jdm*(i-1)) 
        xi1p(1+jdm*(i-1))=xi1p(2+jdm*(i-1)) 
8       continue  

        jdid = jdm*idm  
          do 14 k = 2,kdm-1 
            xk1p(k)=(x(1+jdid*k)-x(1+jdid*(k-1)))*
     .              (x(1+jdid*k)-x(1+jdid*(k-1)))
     .           +(y(1+jdid*k)-y(1+jdid*(k-1)))*
     .            (y(1+jdid*k)-y(1+jdid*(k-1)))  
     .           +(z(1+jdid*k)-z(1+jdid*(k-1)))*
     .            (z(1+jdid*k)-z(1+jdid*(k-1)))  
            xk1m(k)=(x(1+jdid*(k-1))-x(1+jdid*(k-2)))*
     .            (x(1+jdid*(k-1))-x(1+jdid*(k-2)))
     .           +(y(1+jdid*(k-1))-y(1+jdid*(k-2)))*
     .            (y(1+jdid*(k-1))-y(1+jdid*(k-2)))  
     .           +(z(1+jdid*(k-1))-z(1+jdid*(k-2)))*
     .            (z(1+jdid*(k-1))-z(1+jdid*(k-2)))  
            xk2p(k)=(x(jdm+jdid*k)-x(jdm+jdid*(k-1)))
     .              *(x(jdm+jdid*k)-x(jdm+jdid*(k-1)))
     .           +(y(jdm+jdid*k)-y(jdm+jdid*(k-1)))*
     .            (y(jdm+jdid*k)-y(jdm+jdid*(k-1)))  
     .           +(z(jdm+jdid*k)-z(jdm+jdid*(k-1)))*
     .            (z(jdm+jdid*k)-z(jdm+jdid*(k-1)))  
            xk2m(k)=(x(jdm+jdid*(k-1))-x(jdm+jdid*(k-2)))
     .           *(x(jdm+jdid*(k-1))-x(jdm+jdid*(k-2)))
     .           +(y(jdm+jdid*(k-1))-y(jdm+jdid*(k-2)))*
     .            (y(jdm+jdid*(k-1))-y(jdm+jdid*(k-2)))  
     .           +(z(jdm+jdid*(k-1))-z(jdm+jdid*(k-2)))*
     .            (z(jdm+jdid*(k-1))-z(jdm+jdid*(k-2)))  
14        continue
          do 15 k = 2,kdm-1
           xk2m(k) = 1./sqrt(xk2m(k))  
           xk2p(k) = 1./sqrt(xk2p(k))  
           xk1m(k) = 1./sqrt(xk1m(k))  
           xk1p(k) = 1./sqrt(xk1p(k))  
15        continue
          do 16 k = 2,kdm-1
           seta(k) = 1./xk1m(k) 
16        continue   
          seta(kdm) = 1./xk1p(kdm-1)
 
       jdkdm = jdm*idm*(kdm-2)
       do 20 j = 1,jdm*idm-1
         xk1mt =(x(j+jdid)-x(j))*
     .            (x(j+jdid)-x(j))
     .           +(y(j+jdid)-y(j))*
     .            (y(j+jdid)-y(j))  
     .           +(z(j+jdid)-z(j))*
     .            (z(j+jdid)-z(j))  
         xk2mt =(x(j+jdid+jdkdm)-x(j+jdkdm))*
     .            (x(j+jdid+jdkdm)-x(j+jdkdm))
     .           +(y(j+jdid+jdkdm)-y(j+jdkdm))*
     .            (y(j+jdid+jdkdm)-y(j+jdkdm))  
     .           +(z(j+jdid+jdkdm)-z(j+jdkdm))*
     .            (z(j+jdid+jdkdm)-z(j+jdkdm))  
         xj3p(j,1) = sqrt(xk1mt)/seta(2)
         xj3p(j,2) = sqrt(xk2mt)/seta(kdm)
20     continue
        do 24 i = 1,idm
         xj3p(jdm*i,1)       = xj3p(jdm*i-1,1) 
         xj3p(jdm*i,2)       = xj3p(jdm*i-1,2) 
24      continue 


       return
        end 
       subroutine permut_ijkxyz(jdm,idm,kdm,nbl,x,y,z,dvol)
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
       dimension x(jdm*idm*kdm),y(jdm*idm*kdm),
     .               z(jdm*idm*kdm)
       dimension dvol(jdm*idm*kdm)

c***********************************************************************
c     Purpose:  Permutes the indices i,j,k for more efficient computation 
c               within the smoothing step.
c
c
c     Bob Bartels NASA Langley Aeroelasticity Branch    (757)864-2813
c       responsible for the following subroutines:
c         gdsa      -      spring analogy smoothing and near surface
c                          orientation preservation.
c         dmsetup0  -      some initialization of xoi,xoj and dx,dy,dz 
c         dmsetup1  -      some initialization of xok and dx,dy,dz 
c         dmsetup   -      some initialization of arrays 
c         bdata     -      computation of block boundary face spacing
c                          arrays 
c         permut_ijkxyz -  permutes indices i and j to better align 
c                          k = 0 and k = kdim surfaces for vectorizing 
c         unperm_ijkxyz -  unpermutes indices i and j after smoothing 
c     
c***********************************************************************
c
         jdid = jdm*idm
         jdkd = jdm*kdm 
         do 10 k = 1,kdm 
           jd2 = jdid*(k-1) 
           jk1 = jdm*(k-1)  
           do 10 i = 1,idm 
             ji1 = jdm*(i-1) 
             jdi1= jdkd*(i-1) 
             do 10 j = 1,jdm 
               dvol(j+ji1+jd2) = x(j+jk1+jdi1)
10           continue  
         do 11 j = 1,jdm*kdm*idm  
           x(j)  = dvol(j) 
11       continue  
         do 12 k = 1,kdm 
           jd2 = jdid*(k-1) 
           jk1 = jdm*(k-1)  
           do 12 i = 1,idm 
             ji1 = jdm*(i-1) 
             jdi1= jdkd*(i-1) 
             do 12 j = 1,jdm 
               dvol(j+ji1+jd2) = y(j+jk1+jdi1)
12           continue  
         do 13 j = 1,jdm*kdm*idm  
           y(j)  = dvol(j) 
13       continue  
         do 14 k = 1,kdm 
           jd2 = jdid*(k-1) 
           jk1 = jdm*(k-1)  
           do 14 i = 1,idm 
             ji1 = jdm*(i-1) 
             jdi1= jdkd*(i-1) 
             do 14 j = 1,jdm 
              dvol(j+ji1+jd2) = z(j+jk1+jdi1)
14           continue  
         do 15 j = 1,jdm*kdm*idm  
           z(j)  = dvol(j) 
15       continue  

       return
       end
       subroutine unperm_ijkxyz(jdm,idm,kdm,nbl,x,y,z,dvol)
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
       dimension x(jdm*idm*kdm),y(jdm*idm*kdm),
     .               z(jdm*idm*kdm)
       dimension dvol(jdm*idm*kdm)
c***********************************************************************
c     Purpose:  Reverses the permuted indices i,j,k 
c               after the smoothing step.
c
c
c     Bob Bartels NASA Langley Aeroelasticity Branch    (757)864-2813
c       responsible for the following subroutines:
c         gdsa      -      spring analogy smoothing and near surface
c                          orientation preservation.
c         dmsetup0  -      some initialization of xoi,xoj and dx,dy,dz 
c         dmsetup1  -      some initialization of xok and dx,dy,dz 
c         dmsetup   -      some initialization of arrays 
c         bdata     -      computation of block boundary face spacing
c                          arrays 
c         permut_ijkxyz -  permutes indices i and j to better align 
c                          k = 0 and k = kdim surfaces for vectorizing 
c         unperm_ijkxyz -  unpermutes indices i and j after smoothing 
c     
c***********************************************************************

         jdid = jdm*idm
         jdkd = jdm*kdm 
         do 1010 k = 1,kdm 
           jd2 = jdid*(k-1) 
           jk1 = jdm*(k-1)  
           do 1010 i = 1,idm 
             ji1 = jdm*(i-1) 
             jdi1= jdkd*(i-1) 
             do 1010 j = 1,jdm 
               dvol(j+jk1+jdi1)=x(j+ji1+jd2)                
1010         continue  
         do 1011 j = 1,jdm*kdm*idm  
           x(j)  = dvol(j) 
1011     continue  
         do 1012 k = 1,kdm 
           jd2 = jdid*(k-1) 
           jk1 = jdm*(k-1)  
           do 1012 i = 1,idm 
             ji1 = jdm*(i-1) 
             jdi1= jdkd*(i-1) 
             do 1012 j = 1,jdm 
               dvol(j+jk1+jdi1) = y(j+ji1+jd2)
1012         continue  
         do 1013 j = 1,jdm*kdm*idm  
           y(j)  = dvol(j) 
1013     continue  
         do 1014 k = 1,kdm 
           jd2 = jdid*(k-1) 
           jk1 = jdm*(k-1)  
           do 1014 i = 1,idm 
             ji1 = jdm*(i-1) 
             jdi1= jdkd*(i-1) 
             do 1014 j = 1,jdm 
               dvol(j+jk1+jdi1) = z(j+ji1+jd2)
1014         continue  
         do 1015 j = 1,jdm*kdm*idm  
           z(j)  = dvol(j) 
1015     continue  
       return
       end
      subroutine gdsa(time,xmach,dt,dx,dy,dz,x,y,z,xnm1,ynm1,
     .                znm1,dvol,xok,xoi,xoj,irotat,ndeftp1,
     .                iupdat,nbl,nt,idmt,jdmt,kdmt,iaesol,jsta,jend,
     .                ismooth)
c***********************************************************************
c     Purpose:  Determines increment to grid position due to mesh smoothing
c               caused by movement of the solid surface 
c
c
c     Bob Bartels NASA Langley Aeroelasticity Branch    (757)864-2813
c       responsible for the following subroutines:
c         gdsa      -      spring analogy smoothing and near surface
c                          orientation preservation.
c         dmsetup0  -      some initialization of xoi,xoj and dx,dy,dz 
c         dmsetup1  -      some initialization of xok and dx,dy,dz 
c         dmsetup   -      some initialization of arrays 
c         bdata     -      computation of block boundary face spacing
c                          arrays 
c         permut_ijkxyz -  permutes indices i and j to better align 
c                          k = 0 and k = kdim surfaces for vectorizing 
c         unperm_ijkxyz -  unpermutes indices i and j after smoothing 
c     
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
       common /grid1/ algepmn,coef(2)   
       common /grid2/ kimx(2),jte1(2),         
     .   jflp(2,2),icsi(2),icso(2),iafmv1(2),jsp1(2),  
     .   jsp2(2),jsp3(2),itecs(2),ispcs(2),itp(2,2)
     .  ,iaxrt(2)
       common /grid4/ time0
       common /grid8/ fc3,dxave
       common /grid15/ igtyp(2),idmg(2),jdmg(2),kdmg(2)
       common /sklton/ isklton
       common /zero/ iexp
c
       common /twod/ i2d
c
      integer stats
c
      dimension dvol(jdmt*idmt*kdmt)
      dimension dx(jdmt*kdmt*idmt)
      dimension dy(jdmt*kdmt*idmt)
      dimension dz(jdmt*kdmt*idmt)
      dimension jend(2)
      dimension jsta(2)
      dimension x(jdmt*idmt*kdmt)
      dimension xnm1(jdmt*kdmt*idmt)
      dimension xoi(jdmt*kdmt,6)
      dimension xoj(kdmt*idmt,6)
      dimension xok(jdmt*idmt,6)
      dimension y(jdmt*idmt*kdmt)
      dimension ynm1(jdmt*kdmt*idmt)
      dimension z(jdmt*idmt*kdmt)
      dimension znm1(jdmt*kdmt*idmt)

      allocatable :: dxj(:,:)
      allocatable :: dyj(:,:)
      allocatable :: dzj(:,:)
      allocatable :: epsilon(:,:)
      allocatable :: imid(:)
      allocatable :: ipid(:)
      allocatable :: rj1jd(:)
      allocatable :: rj2jd(:)
      allocatable :: seta(:)
      allocatable :: xi1m(:)
      allocatable :: xi1p(:)
      allocatable :: xi2m(:)
      allocatable :: xi2p(:)
      allocatable :: xid(:,:)
      allocatable :: xim(:,:)
      allocatable :: xip(:,:)
      allocatable :: xj1m(:)
      allocatable :: xj1p(:)
      allocatable :: xj2m(:)
      allocatable :: xj2p(:)
      allocatable :: xj3p(:,:)
      allocatable :: xjd(:,:)
      allocatable :: xjm(:,:)
      allocatable :: xjp(:,:)
      allocatable :: xk1m(:)
      allocatable :: xk1p(:)
      allocatable :: xk2m(:)
      allocatable :: xk2p(:)
      allocatable :: xkm(:,:)
      allocatable :: xkp(:,:)

      memuse = 0
c
      allocate( dxj(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'dxj',memuse,stats)
      allocate( dyj(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'dyj',memuse,stats)
      allocate( dzj(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'dzj',memuse,stats)
      allocate( epsilon(kdmt,2), stat=stats)
      call umalloc(kdmt*2,0,'epsilon',memuse,stats)
      allocate( imid(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,1,'imid',memuse,stats)
      allocate( ipid(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,1,'ipid',memuse,stats)
      allocate( rj1jd(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,0,'rj1jd',memuse,stats)
      allocate( rj2jd(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,0,'rj2jd',memuse,stats)
      allocate( seta(kdmt), stat=stats)
      call umalloc(kdmt,0,'seta',memuse,stats)
      allocate( xi1m(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,0,'xi1m',memuse,stats)
      allocate( xi1p(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,0,'xi1p',memuse,stats)
      allocate( xi2m(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,0,'xi2m',memuse,stats)
      allocate( xi2p(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,0,'xi2p',memuse,stats)
      allocate( xid(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'xid',memuse,stats)
      allocate( xim(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'xim',memuse,stats)
      allocate( xip(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'xip',memuse,stats)
      allocate( xj1m(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,0,'xj1m',memuse,stats)
      allocate( xj1p(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,0,'xj1p',memuse,stats)
      allocate( xj2m(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,0,'xj2m',memuse,stats)
      allocate( xj2p(jdmt*idmt), stat=stats)
      call umalloc(jdmt*idmt,0,'xj2p',memuse,stats)
      allocate( xj3p(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'xj3p',memuse,stats)
      allocate( xjd(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'xjd',memuse,stats)
      allocate( xjm(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'xjm',memuse,stats)
      allocate( xjp(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'xjp',memuse,stats)
      allocate( xk1m(kdmt), stat=stats)
      call umalloc(kdmt,0,'xk1m',memuse,stats)
      allocate( xk1p(kdmt), stat=stats)
      call umalloc(kdmt,0,'xk1p',memuse,stats)
      allocate( xk2m(kdmt), stat=stats)
      call umalloc(kdmt,0,'xk2m',memuse,stats)
      allocate( xk2p(kdmt), stat=stats)
      call umalloc(kdmt,0,'xk2p',memuse,stats)
      allocate( xkm(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'xkm',memuse,stats)
      allocate( xkp(jdmt*idmt,2), stat=stats)
      call umalloc(jdmt*idmt*2,0,'xkp',memuse,stats)
c
      idm = idmg(1)
      jdm = jdmg(1)
      kdm = kdmg(1)

      call dmsetup(jdm,kdm,idm,nbl,jsta,jend,ipid,imid,rj1jd,
     .             rj2jd,epsilon)

      call permut_ijkxyz(jdmt,idmt,kdmt,nbl,xnm1,ynm1,znm1,dvol)

      call bdata(jdm,kdm,idm,nbl,xnm1,ynm1,znm1,xj1p,xj1m,xj2p
     .          ,xj2m,xi1p,xi1m,xi2p,xi2m,xk1p,xk1m,xk2p,xk2m,xj3p
     .          ,seta)   


        fc3 =  600.   

        jdi1 = jdm*idm
        jdi2 = 2*jdm*idm
        jdi3 = 3*jdm*idm
        jdi4 = 4*jdm*idm
        jdi5 = 5*jdm*idm



        eps = 10.**(-iexp+1)
        eps1= 10.*eps
        if(idm.gt.2) then 
         ilow = 2
         ihgh = idm - 1
         jdmi1 = jdm*(ilow-1) 
         jdmi2 = jdm*ihgh    
         do 132 j = 2+jdmi1,jdmi2-1  
          dxip =  (xok(ipid(j),1)
     .                  -xok(j,1))         
          dyip =  (xok(ipid(j),2)
     .                  -xok(j,2)) 
          dzip =  (xok(ipid(j),3)
     .                  -xok(j,3)) 
          dxim =  (xok(imid(j),1)
     .                  -xok(j,1))         
          dyim =  (xok(imid(j),2)
     .                  -xok(j,2)) 
          dzim =  (xok(imid(j),3)
     .                  -xok(j,3)) 
          dxi  =  dxip-dxim
          dyi  =  dyip-dyim
          dzi  =  dzip-dzim
          dyjp =  (xok(j+1,2)-xok(j,2))
          dyjm =  (xok(j-1,2)-xok(j,2))  
          dzjp =  (xok(j+1,3)-xok(j,3))  
          dzjm =  (xok(j-1,3)-xok(j,3)) 
          dxjm =  (xok(j-1,1)-xok(j,1))  
          dxjp =  (xok(j+1,1)-xok(j,1)) 
          dxj1 = dxjp-dxjm  
          dyj1 = dyjp-dyjm
          dzj1 = dzjp-dzjm 
          rjp  = sqrt(dxjp*dxjp+dyjp*dyjp+dzjp*dzjp) 
          rjm  = sqrt(dxjm*dxjm+dyjm*dyjm+dzjm*dzjm)  
          rjpm = rjp+rjm
          rip  = sqrt(dxip*dxip+dyip*dyip+dzip*dzip) 
          rim  = sqrt(dxim*dxim+dyim*dyim+dzim*dzim) 
          ripm = rip+rim
          term1   = (dxi*(dyjp*dzjm-dyjm*dzjp)
     .              +dyi*(dzjp*dxjm-dzjm*dxjp) 
     .              +dzi*(dxjp*dyjm-dxjm*dyjp))-eps1        
          term2   = (ripm*(rjp*rjm-
     .            (dxjp*dxjm+dyjp*dyjm+dzjp*dzjm)+eps))   
          xjd(j,1)= term1/term2
          term3   =-(dxj1*(dyip*dzim-dyim*dzip)
     .             +dyj1*(dzip*dxim-dzim*dxip) 
     .             +dzj1*(dxip*dyim-dxim*dyip))-eps1 
          term4    = (rjpm*(rip*rim-
     .            (dxip*dxim+dyip*dyim+dzip*dzim)+eps)) 
          xid(j,1)= term3/term4
          xjp(j,1) = (dxjp*rjm-dxjm*rjp)/rjpm     
          xjm(j,1) = (dyjp*rjm-dyjm*rjp)/rjpm     
          xkm(j,1) = (dzjp*rjm-dzjm*rjp)/rjpm     
          xkp(j,1) = .5*dxi
          xim(j,1) = .5*dyi
          xip(j,1) = .5*dzi
132      continue  
         do 133 j = 2+jdmi1,jdmi2-1  
          dxip =  (xok(ipid(j),4)
     .                  -xok(j,4))         
          dyip =  (xok(ipid(j),5)
     .                  -xok(j,5)) 
          dzip =  (xok(ipid(j),6)
     .                  -xok(j,6)) 
          dxim =  (xok(imid(j),4)
     .                  -xok(j,4))         
          dyim =  (xok(imid(j),5)
     .                  -xok(j,5)) 
          dzim =  (xok(imid(j),6)
     .                  -xok(j,6)) 
          dxi  =  dxip-dxim
          dyi  =  dyip-dyim
          dzi  =  dzip-dzim
          dyjp =  (xok(j+1,5)-xok(j,5))
          dyjm =  (xok(j-1,5)-xok(j,5))  
          dzjp =  (xok(j+1,6)-xok(j,6))  
          dzjm =  (xok(j-1,6)-xok(j,6)) 
          dxjm =  (xok(j-1,4)-xok(j,4))  
          dxjp =  (xok(j+1,4)-xok(j,4)) 
          dxj1 = dxjp-dxjm  
          dyj1 = dyjp-dyjm
          dzj1 = dzjp-dzjm 
          rjp  = sqrt(dxjp*dxjp+dyjp*dyjp+dzjp*dzjp) 
          rjm  = sqrt(dxjm*dxjm+dyjm*dyjm+dzjm*dzjm)  
          rjpm = rjp+rjm
          rip  = sqrt(dxip*dxip+dyip*dyip+dzip*dzip) 
          rim  = sqrt(dxim*dxim+dyim*dyim+dzim*dzim) 
          ripm = rip+rim
          term1   =-(dxi*(dyjp*dzjm-dyjm*dzjp)
     .              +dyi*(dzjp*dxjm-dzjm*dxjp) 
     .              +dzi*(dxjp*dyjm-dxjm*dyjp))-eps1        
          term2   = (ripm*(rjp*rjm-
     .            (dxjp*dxjm+dyjp*dyjm+dzjp*dzjm)+eps))   
          xjd(j,2)= term1/term2
          term3   = (dxj1*(dyip*dzim-dyim*dzip)
     .             +dyj1*(dzip*dxim-dzim*dxip) 
     .             +dzj1*(dxip*dyim-dxim*dyip))-eps1 
          term4    = (rjpm*(rip*rim-
     .            (dxip*dxim+dyip*dyim+dzip*dzim)+eps)) 
          xid(j,2)= term3/term4
          xjp(j,2) = (dxjp*rjm-dxjm*rjp)/rjpm     
          xjm(j,2) = (dyjp*rjm-dyjm*rjp)/rjpm     
          xkm(j,2) = (dzjp*rjm-dzjm*rjp)/rjpm     
          xkp(j,2) = .5*dxi
          xim(j,2) = .5*dyi
          xip(j,2) = .5*dzi
133      continue  
        else
         ilow = 1
         ihgh = 2 
         jdmi1 = jdm*(ilow-1) 
         jdmi2 = jdm*ihgh    
         do 134 j = 2+jdmi1,jdmi2-1  
          dxip =  (xok(ipid(j),1)
     .            -xok(imid(j),1))         
          dyip =  (xok(ipid(j),2)
     .            -xok(imid(j),2)) 
          dzip =  (xok(ipid(j),3)
     .            -xok(imid(j),3)) 
          dxim = -(xok(ipid(j),1)
     .            -xok(imid(j),1))         
          dyim = -(xok(ipid(j),2)
     .            -xok(imid(j),2)) 
          dzim = -(xok(ipid(j),3)
     .            -xok(imid(j),3)) 
          dxi  =  dxip-dxim
          dyi  =  dyip-dyim
          dzi  =  dzip-dzim
          dyjp =  (xok(j+1,2)-xok(j,2))
          dyjm =  (xok(j-1,2)-xok(j,2))  
          dzjp =  (xok(j+1,3)-xok(j,3))  
          dzjm =  (xok(j-1,3)-xok(j,3)) 
          dxjm =  (xok(j-1,1)-xok(j,1))  
          dxjp =  (xok(j+1,1)-xok(j,1)) 
          dxj1 = dxjp-dxjm  
          dyj1 = dyjp-dyjm
          dzj1 = dzjp-dzjm 
          rjp  = sqrt(dxjp*dxjp+dyjp*dyjp+dzjp*dzjp) 
          rjm  = sqrt(dxjm*dxjm+dyjm*dyjm+dzjm*dzjm)  
          rjpm = rjp+rjm
          rip  = sqrt(dxip*dxip+dyip*dyip+dzip*dzip) 
          rim  = sqrt(dxim*dxim+dyim*dyim+dzim*dzim) 
          ripm = rip+rim
          term1   = (dxi*(dyjp*dzjm-dyjm*dzjp)
     .              +dyi*(dzjp*dxjm-dzjm*dxjp) 
     .              +dzi*(dxjp*dyjm-dxjm*dyjp))-eps1        
          term2   = (ripm*(rjp*rjm-
     .            (dxjp*dxjm+dyjp*dyjm+dzjp*dzjm)+eps))   
          xjd(j,1)= term1/term2
          term3   =-(dxj1*(dyip*dzim-dyim*dzip)
     .              +dyj1*(dzip*dxim-dzim*dxip) 
     .             +dzj1*(dxip*dyim-dxim*dyip))-eps1 
          term4    = (rjpm*(rip*rim-
     .            (dxip*dxim+dyip*dyim+dzip*dzim)+eps)) 
          xid(j,1)= term3/term4
          xjp(j,1) = (dxjp*rjm-dxjm*rjp)/rjpm     
          xjm(j,1) = (dyjp*rjm-dyjm*rjp)/rjpm     
          xkm(j,1) = (dzjp*rjm-dzjm*rjp)/rjpm     
          xkp(j,1) = .5*dxi
          xim(j,1) = .5*dyi
          xip(j,1) = .5*dzi
134      continue  
         do 135 j = 2+jdmi1,jdmi2-1  
          dxip =  (xok(ipid(j),4)
     .            -xok(imid(j),4))         
          dyip =  (xok(ipid(j),5)
     .            -xok(imid(j),5)) 
          dzip =  (xok(ipid(j),6)
     .            -xok(imid(j),6)) 
          dxim = -(xok(ipid(j),4)
     .            -xok(imid(j),4))         
          dyim = -(xok(ipid(j),5)
     .            -xok(imid(j),5)) 
          dzim = -(xok(ipid(j),6)
     .            -xok(imid(j),6)) 
          dxi  =  dxip-dxim
          dyi  =  dyip-dyim
          dzi  =  dzip-dzim
          dyjp =  (xok(j+1,5)-xok(j,5))
          dyjm =  (xok(j-1,5)-xok(j,5))  
          dzjp =  (xok(j+1,6)-xok(j,6))  
          dzjm =  (xok(j-1,6)-xok(j,6)) 
          dxjm =  (xok(j-1,4)-xok(j,4))  
          dxjp =  (xok(j+1,4)-xok(j,4)) 
          dxj1 = dxjp-dxjm  
          dyj1 = dyjp-dyjm
          dzj1 = dzjp-dzjm 
          rjp  = sqrt(dxjp*dxjp+dyjp*dyjp+dzjp*dzjp) 
          rjm  = sqrt(dxjm*dxjm+dyjm*dyjm+dzjm*dzjm)  
          rjpm = rjp+rjm
          rip  = sqrt(dxip*dxip+dyip*dyip+dzip*dzip) 
          rim  = sqrt(dxim*dxim+dyim*dyim+dzim*dzim) 
          ripm = rip+rim
          term1   =-(dxi*(dyjp*dzjm-dyjm*dzjp)
     .              +dyi*(dzjp*dxjm-dzjm*dxjp) 
     .              +dzi*(dxjp*dyjm-dxjm*dyjp))-eps1        
          term2   = (ripm*(rjp*rjm-
     .            (dxjp*dxjm+dyjp*dyjm+dzjp*dzjm)+eps))   
          xjd(j,2)= term1/term2
          term3   = (dxj1*(dyip*dzim-dyim*dzip)
     .             +dyj1*(dzip*dxim-dzim*dxip) 
     .             +dzj1*(dxip*dyim-dxim*dyip))-eps1 
          term4    = (rjpm*(rip*rim-
     .            (dxip*dxim+dyip*dyim+dzip*dzim)+eps)) 
          xid(j,2)= term3/term4
          xjp(j,2) = (dxjp*rjm-dxjm*rjp)/rjpm     
          xjm(j,2) = (dyjp*rjm-dyjm*rjp)/rjpm     
          xkm(j,2) = (dzjp*rjm-dzjm*rjp)/rjpm     
          xkp(j,2) = .5*dxi
          xim(j,2) = .5*dyi
          xip(j,2) = .5*dzi
135      continue  
        end if  
        do 136 j = 2+jdmi1,jdmi2-1 
          xjd(j,1) = exp(fc3*xjd(j,1))
          xid(j,1) = exp(fc3*xid(j,1))
          xjd(j,2) = exp(fc3*xjd(j,2))
          xid(j,2) = exp(fc3*xid(j,2))
136     continue
        arg1 = exp(-12.)
        arg2 = 1. 
        do 138 j = 2+jdmi1,jdmi2-1 
          if(real(xjd(j,1)).lt.real(arg1)) xjd(j,1)=arg1 
          if(real(xjd(j,1)).gt.real(arg2)) xjd(j,1)=arg2 
          if(real(xid(j,1)).lt.real(arg1)) xid(j,1)=arg1 
          if(real(xid(j,1)).gt.real(arg2)) xid(j,1)=arg2 
          if(real(xjd(j,2)).lt.real(arg1)) xjd(j,2)=arg1 
          if(real(xjd(j,2)).gt.real(arg2)) xjd(j,2)=arg2 
          if(real(xid(j,2)).lt.real(arg1)) xid(j,2)=arg1 
          if(real(xid(j,2)).gt.real(arg2)) xid(j,2)=arg2 
138     continue 
        do 141 i=ilow,ihgh  
        jdmi1 = jdm*(i-1)
        xjd(1+jdmi1,1)    = xjd(2+jdmi1,1) 
        xjd(jdm+jdmi1,1)  = xjd(jdm-1+jdmi1,1) 
        xid(1+jdmi1,1)    = xid(2+jdmi1,1) 
        xid(jdm+jdmi1,1)  = xid(jdm-1+jdmi1,1) 
        xjd(1+jdmi1,2)    = xjd(2+jdmi1,2) 
        xjd(jdm+jdmi1,2)  = xjd(jdm-1+jdmi1,2) 
        xid(1+jdmi1,2)    = xid(2+jdmi1,2) 
        xid(jdm+jdmi1,2)  = xid(jdm-1+jdmi1,2) 
        xjp(1+jdmi1,1)   = xjp(2+jdmi1,1)
        xjp(jdm+jdmi1,1) = xjp(jdm-1+jdmi1,1)
        xjm(1+jdmi1,1)   = xjm(2+jdmi1,1)
        xjm(jdm+jdmi1,1) = xjm(jdm-1+jdmi1,1)
        xkm(1+jdmi1,1)   = xkm(2+jdmi1,1)
        xkm(jdm+jdmi1,1) = xkm(jdm-1+jdmi1,1)
        xkp(1+jdmi1,1)   = xkp(2+jdmi1,1)
        xkp(jdm+jdmi1,1) = xkp(jdm-1+jdmi1,1)
        xim(1+jdmi1,1)   = xim(2+jdmi1,1)
        xim(jdm+jdmi1,1) = xim(jdm-1+jdmi1,1)
        xip(1+jdmi1,1)   = xip(2+jdmi1,1)
        xip(jdm+jdmi1,1) = xip(jdm-1+jdmi1,1)
        xjp(1+jdmi1,2)   = xjp(2+jdmi1,2)
        xjp(jdm+jdmi1,2) = xjp(jdm-1+jdmi1,2)
        xjm(1+jdmi1,2)   = xjm(2+jdmi1,2)
        xjm(jdm+jdmi1,2) = xjm(jdm-1+jdmi1,2)
        xkm(1+jdmi1,2)   = xkm(2+jdmi1,2)
        xkm(jdm+jdmi1,2) = xkm(jdm-1+jdmi1,2)
        xkp(1+jdmi1,2)   = xkp(2+jdmi1,2)
        xkp(jdm+jdmi1,2) = xkp(jdm-1+jdmi1,2)
        xim(1+jdmi1,2)   = xim(2+jdmi1,2)
        xim(jdm+jdmi1,2) = xim(jdm-1+jdmi1,2)
        xip(1+jdmi1,2)   = xip(2+jdmi1,2)
        xip(jdm+jdmi1,2) = xip(jdm-1+jdmi1,2)
 141      continue 

        jdmi1 = jdm*(ilow-1) 
        do j = 2+jdmi1,jdmi2-1
          if(real(xjd(j,1)).lt.1.e-2) xjd(j,2) = xjd(j,1)
          if(real(xjd(j,2)).lt.1.e-2) xjd(j,1) = xjd(j,2)
          if(real(xid(j,1)).lt.1.e-2) xid(j,2) = xid(j,1)
          if(real(xid(j,2)).lt.1.e-2) xid(j,1) = xid(j,2)
        enddo
         
         
        jdmi1 = jdm*(ilow-1) 
        jdmi2 = jdm*ihgh    
          do 154 j = 1+jdmi1,jdmi2    
            dzj(j,1) =-(xjp(j,1)*xim(j,1)-xjm(j,1)*xkp(j,1)) 
            dxj(j,1) =-(xjm(j,1)*xip(j,1)-xkm(j,1)*xim(j,1))  
            dzj(j,2) =-(xjp(j,2)*xim(j,2)-xjm(j,2)*xkp(j,2)) 
            dxj(j,2) =-(xjm(j,2)*xip(j,2)-xkm(j,2)*xim(j,2))  
154       continue  
        do 156 j = 1+jdmi1,jdmi2    
          term1 =  sqrt((xkm(j,1)*xkp(j,1)-xjp(j,1)*xip(j,1))  
     .              *(xkm(j,1)*xkp(j,1)-xjp(j,1)*xip(j,1))  
     .                    +dxj(j,1)*dxj(j,1)+dzj(j,1)*dzj(j,1))   
          term2 =  sqrt((xkm(j,2)*xkp(j,2)-xjp(j,2)*xip(j,2))  
     .              *(xkm(j,2)*xkp(j,2)-xjp(j,2)*xip(j,2))  
     .                    +dxj(j,2)*dxj(j,2)+dzj(j,2)*dzj(j,2))   
          dyj(j,1) = 1./(term1+eps) 
          dxj(j,1) = dxj(j,1)*dyj(j,1) 
          dzj(j,1) = dzj(j,1)*dyj(j,1) 
          dyj(j,2) = 1./(term2+eps) 
          dxj(j,2) = dxj(j,2)*dyj(j,2) 
          dzj(j,2) = dzj(j,2)*dyj(j,2) 
156     continue  
          
       if(i2d.eq.0) then            
          do 222 j = 1+jdmi1,jdmi2
           dyj(j,1) =-(xkm(j,1)*xkp(j,1)-xjp(j,1)*xip(j,1))*dyj(j,1)     
           dyj(j,2) =-(xkm(j,2)*xkp(j,2)-xjp(j,2)*xip(j,2))*dyj(j,2)     
 222      continue  
       else 
         do 225 j = 1,jdm*idm
           xid(j,1) = 0.    
           xid(j,2) = 0.    
           dyj(j,1) = 0.
           dyj(j,2) = 0.
 225      continue
       end if 




        if(idm.gt.2) then 
         ilow = 2
         ihgh = idm - 1
         jdmi1 = jdm*(ilow-1) 
         jdmi2 = jdm*ihgh    
         do 252 j = 2+jdmi1,jdmi2-1  
          dxip =  (x(ipid(j)) -x(j))         
          dyip =  (y(ipid(j)) -y(j)) 
          dzip =  (z(ipid(j)) -z(j)) 
          dxim =  (x(imid(j)) -x(j))         
          dyim =  (y(imid(j)) -y(j)) 
          dzim =  (z(imid(j)) -z(j)) 
          dxi  =  dxip-dxim
          dyi  =  dyip-dyim
          dzi  =  dzip-dzim
          dyjp =  (y(j+1)-y(j))
          dyjm =  (y(j-1)-y(j))  
          dzjp =  (z(j+1)-z(j))  
          dzjm =  (z(j-1)-z(j)) 
          dxjm =  (x(j-1)-x(j))  
          dxjp =  (x(j+1)-x(j)) 
          rjp  = sqrt(dxjp*dxjp+dyjp*dyjp+dzjp*dzjp) 
          rjm  = sqrt(dxjm*dxjm+dyjm*dyjm+dzjm*dzjm)  
          rjpm = rjp+rjm
          rip  = sqrt(dxip*dxip+dyip*dyip+dzip*dzip) 
          rim  = sqrt(dxim*dxim+dyim*dyim+dzim*dzim) 
          ripm = rip+rim
          xjp(j,1) = (dxjp*rjm-dxjm*rjp)/rjpm     
          xjm(j,1) = (dyjp*rjm-dyjm*rjp)/rjpm     
          xkm(j,1) = (dzjp*rjm-dzjm*rjp)/rjpm     
          xkp(j,1) = .5*dxi
          xim(j,1) = .5*dyi
          xip(j,1) = .5*dzi
252      continue  
         jdkdm = jdm*idm*(kdm-1)
         do 253 j = 2+jdmi1,jdmi2-1  
          dxip =  (x(jdkdm+ipid(j)) -x(jdkdm+j))         
          dyip =  (y(jdkdm+ipid(j)) -y(jdkdm+j)) 
          dzip =  (z(jdkdm+ipid(j)) -z(jdkdm+j)) 
          dxim =  (x(jdkdm+imid(j)) -x(jdkdm+j))         
          dyim =  (y(jdkdm+imid(j)) -y(jdkdm+j)) 
          dzim =  (z(jdkdm+imid(j)) -z(jdkdm+j)) 
          dxi  =  dxip-dxim
          dyi  =  dyip-dyim
          dzi  =  dzip-dzim
          dyjp =  (y(jdkdm+j+1)-y(jdkdm+j))
          dyjm =  (y(jdkdm+j-1)-y(jdkdm+j))  
          dzjp =  (z(jdkdm+j+1)-z(jdkdm+j))  
          dzjm =  (z(jdkdm+j-1)-z(jdkdm+j)) 
          dxjm =  (x(jdkdm+j-1)-x(jdkdm+j))  
          dxjp =  (x(jdkdm+j+1)-x(jdkdm+j)) 
          rjp  = sqrt(dxjp*dxjp+dyjp*dyjp+dzjp*dzjp) 
          rjm  = sqrt(dxjm*dxjm+dyjm*dyjm+dzjm*dzjm)  
          rjpm = rjp+rjm
          rip  = sqrt(dxip*dxip+dyip*dyip+dzip*dzip) 
          rim  = sqrt(dxim*dxim+dyim*dyim+dzim*dzim) 
          ripm = rip+rim
          xjp(j,2) = (dxjp*rjm-dxjm*rjp)/rjpm     
          xjm(j,2) = (dyjp*rjm-dyjm*rjp)/rjpm     
          xkm(j,2) = (dzjp*rjm-dzjm*rjp)/rjpm     
          xkp(j,2) = .5*dxi
          xim(j,2) = .5*dyi
          xip(j,2) = .5*dzi
253      continue  
        else
         ilow = 1
         ihgh = 2 
         jdmi1 = jdm*(ilow-1) 
         jdmi2 = jdm*ihgh    
         do 254 j = 2+jdmi1,jdmi2-1  
          dxip =  (x(ipid(j)) -x(imid(j)))         
          dyip =  (y(ipid(j)) -y(imid(j))) 
          dzip =  (z(ipid(j)) -z(imid(j))) 
          dxim = -(x(ipid(j)) -x(imid(j)))         
          dyim = -(y(ipid(j)) -y(imid(j))) 
          dzim = -(z(ipid(j)) -z(imid(j))) 
          dxi  =  dxip-dxim
          dyi  =  dyip-dyim
          dzi  =  dzip-dzim
          dyjp =  (y(j+1)-y(j))
          dyjm =  (y(j-1)-y(j))  
          dzjp =  (z(j+1)-z(j))  
          dzjm =  (z(j-1)-z(j)) 
          dxjm =  (x(j-1)-x(j))  
          dxjp =  (x(j+1)-x(j)) 
          rjp  = sqrt(dxjp*dxjp+dyjp*dyjp+dzjp*dzjp) 
          rjm  = sqrt(dxjm*dxjm+dyjm*dyjm+dzjm*dzjm)  
          rjpm = rjp+rjm
          rip  = sqrt(dxip*dxip+dyip*dyip+dzip*dzip) 
          rim  = sqrt(dxim*dxim+dyim*dyim+dzim*dzim) 
          ripm = rip+rim
          xjp(j,1) = (dxjp*rjm-dxjm*rjp)/rjpm     
          xjm(j,1) = (dyjp*rjm-dyjm*rjp)/rjpm     
          xkm(j,1) = (dzjp*rjm-dzjm*rjp)/rjpm     
          xkp(j,1) = .5*dxi
          xim(j,1) = .5*dyi
          xip(j,1) = .5*dzi
254      continue
         jdkdm = jdm*idm*(kdm-1)
         do 255 j = 2+jdmi1,jdmi2-1  
          dxip =  (x(jdkdm+ipid(j)) -x(jdkdm+imid(j)))         
          dyip =  (y(jdkdm+ipid(j)) -y(jdkdm+imid(j))) 
          dzip =  (z(jdkdm+ipid(j)) -z(jdkdm+imid(j))) 
          dxim = -(x(jdkdm+ipid(j)) -x(jdkdm+imid(j)))         
          dyim = -(y(jdkdm+ipid(j)) -y(jdkdm+imid(j))) 
          dzim = -(z(jdkdm+ipid(j)) -z(jdkdm+imid(j))) 
          dxi  =  dxip-dxim
          dyi  =  dyip-dyim
          dzi  =  dzip-dzim
          dyjp =  (y(jdkdm+j+1)-y(jdkdm+j))
          dyjm =  (y(jdkdm+j-1)-y(jdkdm+j))  
          dzjp =  (z(jdkdm+j+1)-z(jdkdm+j))  
          dzjm =  (z(jdkdm+j-1)-z(jdkdm+j)) 
          dxjm =  (x(jdkdm+j-1)-x(jdkdm+j))  
          dxjp =  (x(jdkdm+j+1)-x(jdkdm+j)) 
          rjp  = sqrt(dxjp*dxjp+dyjp*dyjp+dzjp*dzjp) 
          rjm  = sqrt(dxjm*dxjm+dyjm*dyjm+dzjm*dzjm)  
          rjpm = rjp+rjm
          rip  = sqrt(dxip*dxip+dyip*dyip+dzip*dzip) 
          rim  = sqrt(dxim*dxim+dyim*dyim+dzim*dzim) 
          ripm = rip+rim
          xjp(j,2) = (dxjp*rjm-dxjm*rjp)/rjpm     
          xjm(j,2) = (dyjp*rjm-dyjm*rjp)/rjpm     
          xkm(j,2) = (dzjp*rjm-dzjm*rjp)/rjpm     
          xkp(j,2) = .5*dxi
          xim(j,2) = .5*dyi
          xip(j,2) = .5*dzi
255      continue  
        end if  
        do 261 i=ilow,ihgh  
        jdmi1 = jdm*(i-1)
        xjp(1+jdmi1,1)   = xjp(2+jdmi1,1)
        xjp(jdm+jdmi1,1) = xjp(jdm-1+jdmi1,1)
        xjm(1+jdmi1,1)   = xjm(2+jdmi1,1)
        xjm(jdm+jdmi1,1) = xjm(jdm-1+jdmi1,1)
        xkm(1+jdmi1,1)   = xkm(2+jdmi1,1)
        xkm(jdm+jdmi1,1) = xkm(jdm-1+jdmi1,1)
        xkp(1+jdmi1,1)   = xkp(2+jdmi1,1)
        xkp(jdm+jdmi1,1) = xkp(jdm-1+jdmi1,1)
        xim(1+jdmi1,1)   = xim(2+jdmi1,1)
        xim(jdm+jdmi1,1) = xim(jdm-1+jdmi1,1)
        xip(1+jdmi1,1)   = xip(2+jdmi1,1)
        xip(jdm+jdmi1,1) = xip(jdm-1+jdmi1,1)
        xjp(1+jdmi1,2)   = xjp(2+jdmi1,2)
        xjp(jdm+jdmi1,2) = xjp(jdm-1+jdmi1,2)
        xjm(1+jdmi1,2)   = xjm(2+jdmi1,2)
        xjm(jdm+jdmi1,2) = xjm(jdm-1+jdmi1,2)
        xkm(1+jdmi1,2)   = xkm(2+jdmi1,2)
        xkm(jdm+jdmi1,2) = xkm(jdm-1+jdmi1,2)
        xkp(1+jdmi1,2)   = xkp(2+jdmi1,2)
        xkp(jdm+jdmi1,2) = xkp(jdm-1+jdmi1,2)
        xim(1+jdmi1,2)   = xim(2+jdmi1,2)
        xim(jdm+jdmi1,2) = xim(jdm-1+jdmi1,2)
        xip(1+jdmi1,2)   = xip(2+jdmi1,2)
        xip(jdm+jdmi1,2) = xip(jdm-1+jdmi1,2)
 261      continue 
         

        jdmi1 = jdm*(ilow-1) 
        jdmi2 = jdm*ihgh    
        do 274 j = 1+jdmi1,jdmi2  
          dvol(j+jdi2) =-(xjp(j,1)*xim(j,1)-xjm(j,1)*xkp(j,1)) 
          dvol(j)      =-(xjm(j,1)*xip(j,1)-xkm(j,1)*xim(j,1))  
          dvol(j+jdi5) =-(xjp(j,2)*xim(j,2)-xjm(j,2)*xkp(j,2)) 
          dvol(j+jdi3) =-(xjm(j,2)*xip(j,2)-xkm(j,2)*xim(j,2))  
274     continue  
        do 276 j = 1+jdmi1,jdmi2    
          term1 =  sqrt((xkm(j,1)*xkp(j,1)-xjp(j,1)*xip(j,1))  
     .              *(xkm(j,1)*xkp(j,1)-xjp(j,1)*xip(j,1))  
     .            +dvol(j)*dvol(j)+dvol(j+jdi2)*dvol(j+jdi2))   
          term2 =  sqrt((xkm(j,2)*xkp(j,2)-xjp(j,2)*xip(j,2))  
     .              *(xkm(j,2)*xkp(j,2)-xjp(j,2)*xip(j,2))  
     .         +dvol(j+jdi3)*dvol(j+jdi3)+dvol(j+jdi5)*dvol(j+jdi5))   
          dvol(j+jdi1) = 1./(term1+eps) 
          dvol(j)      = dvol(j     )*dvol(j+jdi1) 
          dvol(j+jdi2) = dvol(j+jdi2)*dvol(j+jdi1) 
          dvol(j+jdi4) = 1./(term2+eps) 
          dvol(j+jdi3) = dvol(j+jdi3)*dvol(j+jdi4) 
          dvol(j+jdi5) = dvol(j+jdi5)*dvol(j+jdi4) 
276     continue  
          
       if(i2d.eq.0) then            
          do 292 j = 1+jdmi1,jdmi2
            dvol(j+jdi1) =-(xkm(j,1)*xkp(j,1)-xjp(j,1)*xip(j,1))*
     .                                              dvol(j+jdi1)     
            dvol(j+jdi4) =-(xkm(j,2)*xkp(j,2)-xjp(j,2)*xip(j,2))*
     .                                              dvol(j+jdi4)     
 292      continue  
       else 
         do 295 j = 1,jdm*idm
           dvol(j+jdi1) = 0.
           dvol(j+jdi4) = 0.
 295      continue
       end if 

       do 298 j = 1+jdmi1,jdmi2
         dxj(j,1) = dxj(j,1) - dvol(j)
         dyj(j,1) = dyj(j,1) - dvol(j+jdi1)
         dzj(j,1) = dzj(j,1) - dvol(j+jdi2)
         dxj(j,2) =-dxj(j,2) + dvol(j+jdi3)
         dyj(j,2) =-dyj(j,2) + dvol(j+jdi4)
         dzj(j,2) =-dzj(j,2) + dvol(j+jdi5)
298    continue 

       tol = .5e-9
       zzero = 10.**(-iexp)
       if (real(tol).lt.real(zzero)) then
          tol = 10.**(-iexp+1)
       end if

       ixflag = 0
       iyflag = 0
       izflag = 0
       if(idm.gt.2) then 
         ilow = 2
         ihgh = idm - 1
       else
         ilow = 1
         ihgh = 2 
       end if  
       jdid = jdm*idm  

       rkdok1= sqrt((x(1+jdid*(kdm-1))-x(1))
     .            *(x(1+jdid*(kdm-1))-x(1))+
     .             (y(1+jdid*(kdm-1))-y(1))
     .            *(y(1+jdid*(kdm-1))-y(1))+
     .             (z(1+jdid*(kdm-1))-z(1))
     .            *(z(1+jdid*(kdm-1))-z(1)))
       rkdok2 = rkdok1

       nnt = 1
       kdmot = kdm/2 + 1
       do 400 iter = 1,ismooth                   
        xdel = 0.
        ydel = 0.
        zdel = 0.  
        orjdm = 1./real(jdm-1) 


        xk1mk = xk1m(2)*orjdm  
        xk2mk = xk2m(2)*orjdm   
        xk3mk = xk1p(kdm-1)*orjdm  
        xk4mk = xk2p(kdm-1)*orjdm   
        do 309 j = jdm*(ilow-1)+2,jdm*ihgh-1 
           xkp(j,1) = (rj1jd(j)*xk1mk+
     .          rj2jd(j)*xk2mk)       
           xkp(j,2) = (rj1jd(j)*xk3mk+
     .          rj2jd(j)*xk4mk)       
309     continue 

        rk11  = 0.
        rk21  = 0.



        jdkdm3 = jdid*(kdm-1)
c       do 445 k = 2,kdmot
        do 445 k = 2,kdm-1
           jdk1 = jdid*(k-1)  
           jdk2 = jdid*(k-2)
           jdk  = jdid*k  
           jdkdm1 = jdid*(kdm-k)  
           jdkdm2 = jdid*(kdm-k+1)
           jdkdm  = jdid*(kdm-k-1)  
           rkdm1= real(kdm-1) 
           kdm1 = int(rkdm1)
           rk11   = rk11 + seta(k)
           rk21   = rk21 + seta(kdm-k+2)
           rkdok  = (rkdok1-rk11)/rkdok1
           rk1    = rk11/rkdok1 
           if(real(rkdok).lt.0.) rkdok = 0. 
           if(real(rk1).gt.1.) rk1 = 1. 
           rkdok2 = (rkdok1-rk21)/rkdok1
           rk2    = rk21/rkdok1 
           if(real(rkdok2).lt.0.) rkdok2 = 0. 
           if(real(rk2).gt.1.) rk2 = 1. 
           xk1mk = xk1m(k)  
           xk2mk = xk2m(k)  
           xk1pk = xk1p(k)*orjdm  
           xk2pk = xk2p(k)*orjdm  
           xk3pk = xk1m(kdm-k+1)*orjdm  
           xk4pk = xk2m(kdm-k+1)*orjdm  
           do 310 j = jdm*(ilow-1)+2,jdm*ihgh-1 
           xkm(j,1) = xkp(j,1) 
310        xkm(j,2) = xkp(j,2) 
           do 311 j = jdm*(ilow-1)+2,jdm*ihgh-1 
           xjm(j,1) =(rkdok*xj1m(j)*xjd(j,1)+rk1*xj2m(j)
     .             *xjd(j,2))        
           xjp(j,1) =(rkdok*xj1p(j)*xjd(j,1)+rk1*xj2p(j)          
     .             *xjd(j,2))        
           xkp(j,1) =(rj1jd(j)*xk1pk+rj2jd(j)
     .              *xk2pk)      
           xjm(j,2) =(rkdok2*xj2m(j)*xjd(j,2)+rk2*xj1m(j)
     .             *xjd(j,1))        
           xjp(j,2) =(rkdok2*xj2p(j)*xjd(j,2)+rk2*xj1p(j)          
     .             *xjd(j,1))        
           xkp(j,2) =(rj1jd(j)*xk3pk+rj2jd(j)
     .              *xk4pk)      
311        continue 
           if(idm.gt.2) then
             do 322 j = jdm*(ilow-1)+2,jdm*ihgh-1 
              xim(j,1)=(rkdok*xi1m(j)*xid(j,1)+rk1*xi2m(j)
     .             *xid(j,2))
              xip(j,1)=(rkdok*xi1p(j)*xid(j,1)+rk1*xi2p(j)
     .             *xid(j,2))
              xim(j,2)=(rkdok2*xi2m(j)*xid(j,2)+rk2*xi1m(j)
     .             *xid(j,1))
              xip(j,2)=(rkdok2*xi2p(j)*xid(j,2)+rk2*xi1p(j)
     .             *xid(j,1))
322          continue  
             jdm2 = jdm 
             jdm3 =-jdm 
           else
             do 323 j = 1,jdm*idm
               xim(j,1) = 0.0
               xip(j,1) = 0.0  
               xim(j,2) = 0.0
               xip(j,2) = 0.0  
323          continue 
             jdm2 = 0
             jdm3 = 0 
           end if 
           eps1 = epsilon(k,1)  
           eps2 = 1.-epsilon(k,1)
           eps3 = epsilon(k,2)  
           eps4 = 1.-epsilon(k,2)
           seta1= seta(k)
           seta2= seta(kdm-k+2)
           do 324 j=jdm*(ilow-1)+2,jdm*ihgh-1      
           dvol(j       )=eps2*(xjp(j,1)*dx(j+1+jdk1)+xjm(j,1)*  
     .           dx(j-1+jdk1)+xkp(j,1)*dx(j+jdk)+xkm(j,1)*dx(j+jdk2)         
     .           +xip(j,1)*dx(j+jdm2+jdk1)+xim(j,1)*dx(j+jdm3+jdk1))   
     .          /(xjp(j,1)+xjm(j,1)+xkp(j,1)+xkm(j,1)+xip(j,1)+xim(j,1))           
     .          + eps1*(dx(j+jdk2)+seta1*dxj(j,1)*
     .           xj3p(j,1))    
           dvol(j+jdkdm3)=eps4*(xjp(j,2)*dx(j+1+jdkdm1)+xjm(j,2)*  
     .         dx(j-1+jdkdm1)+xkp(j,2)*dx(j+jdkdm)+xkm(j,2)*dx(j+jdkdm2)         
     .           +xip(j,2)*dx(j+jdm2+jdkdm1)+xim(j,2)*dx(j+jdm3+jdkdm1))   
     .          /(xjp(j,2)+xjm(j,2)+xkp(j,2)+xkm(j,2)+xip(j,2)+xim(j,2))           
     .          + eps3*(dx(j+jdkdm2)+seta2*dxj(j,2)*
     .           xj3p(j,2))    
324        continue 
    
           do 325 i=1,idm 
             dvol(jdm*i       )       = dx(jdm*i+jdk1)  
             dvol(1+jdm*(i-1)       ) = dx(1+jdm*(i-1)+jdk1)  
             dvol(jdm*i+jdkdm3)       = dx(jdm*i+jdkdm1)  
             dvol(1+jdm*(i-1)+jdkdm3) = dx(1+jdm*(i-1)+jdkdm1)  
325        continue  
           if(idm.gt.2) then 
             jdmi2 = jdm*(idm-1) 
             do 326 j  = 1,jdm    
               dvol(j       )      = dx(j+jdk1)
               dvol(j+jdmi2       )= dx(j+jdmi2+jdk1)
               dvol(j+jdkdm3)      = dx(j+jdkdm1)
               dvol(j+jdmi2+jdkdm3)= dx(j+jdmi2+jdkdm1)
326          continue  
           end if 
        do j=1,jdm*idm        
             xdel = xdel+(dvol(j     )-dx(j+jdk1))*
     .                   (dvol(j     )-dx(j+jdk1))
     .                  +(dvol(j+jdkdm3)-dx(j+jdkdm1))*
     .                   (dvol(j+jdkdm3)-dx(j+jdkdm1))
        end do
        if(jdk1.ne.jdkdm1) then
          do  j=1,jdm*idm                      
             dx(j+jdk1  ) = dvol(j       )
             dx(j+jdkdm1) = dvol(j+jdkdm3)
          enddo
        else
          do  j=1,jdm*idm                      
             dx(j+jdk1  ) = .5*(dvol(j       )+dvol(j+jdkdm3)) 
          enddo
        end if
           do 339 j=jdm*(ilow-1)+2,jdm*ihgh-1 
           dvol(j     )=eps2*(xjp(j,1)*dz(j+1+jdk1)+xjm(j,1)*  
     .           dz(j-1+jdk1)+xkp(j,1)*dz(j+jdk)+xkm(j,1)*dz(j+jdk2)         
     .           +xip(j,1)*dz(j+jdm2+jdk1)+xim(j,1)*dz(j+jdm3+jdk1))   
     .          /(xjp(j,1)+xjm(j,1)+xkp(j,1)+xkm(j,1)+xip(j,1)+xim(j,1))           
     .          + eps1*(dz(j+jdk2)+seta1*dzj(j,1)*
     .           xj3p(j,1))    
           dvol(j+jdkdm3)=eps4*(xjp(j,2)*dz(j+1+jdkdm1)+xjm(j,2)*  
     .         dz(j-1+jdkdm1)+xkp(j,2)*dz(j+jdkdm)+xkm(j,2)*dz(j+jdkdm2)         
     .           +xip(j,2)*dz(j+jdm2+jdkdm1)+xim(j,2)*dz(j+jdm3+jdkdm1))   
     .          /(xjp(j,2)+xjm(j,2)+xkp(j,2)+xkm(j,2)+xip(j,2)+xim(j,2))           
     .          + eps3*(dz(j+jdkdm2)+seta2*dzj(j,2)*
     .           xj3p(j,2))    
339        continue  

           do 341 i=1,idm 
             dvol(jdm*i       )       = dz(jdm*i+jdk1)  
             dvol(1+jdm*(i-1)       ) = dz(1+jdm*(i-1)+jdk1)  
             dvol(jdm*i+jdkdm3)       = dz(jdm*i+jdkdm1)  
             dvol(1+jdm*(i-1)+jdkdm3) = dz(1+jdm*(i-1)+jdkdm1)  
341        continue  
           if(idm.gt.2) then 
             jdmi2 = jdm*(idm-1) 
             do 343 j  = 1,jdm    
               dvol(j       )      = dz(j+jdk1)
               dvol(j+jdmi2       )= dz(j+jdmi2+jdk1)
               dvol(j+jdkdm3)      = dz(j+jdkdm1)
               dvol(j+jdmi2+jdkdm3)= dz(j+jdmi2+jdkdm1)
343          continue  
           end if 
               
        do j=1,jdm*idm          
             zdel = zdel+(dvol(j     )-dz(j+jdk1))*
     .                   (dvol(j     )-dz(j+jdk1))
     .                  +(dvol(j+jdkdm3)-dz(j+jdkdm1))*
     .                   (dvol(j+jdkdm3)-dz(j+jdkdm1))
        end do
        if(jdk1.ne.jdkdm1) then
          do  j=1,jdm*idm                      
             dz(j+jdk1  ) = dvol(j       )
             dz(j+jdkdm1) = dvol(j+jdkdm3)
          enddo
        else
          do  j=1,jdm*idm                      
             dz(j+jdk1  ) = .5*(dvol(j       )+dvol(j+jdkdm3)) 
          enddo
        end if
        do 416 j=jdm*(ilow-1)+2,jdm*ihgh-1      
         dvol(j     )=eps2*(xjp(j,1)*dy(j+1+jdk1)+xjm(j,1)*  
     .           dy(j-1+jdk1)+xkp(j,1)*dy(j+jdk)+xkm(j,1)*dy(j+jdk2)         
     .           +xip(j,1)*dy(j+jdm2+jdk1)+xim(j,1)*dy(j+jdm3+jdk1))   
     .          /(xjp(j,1)+xjm(j,1)+xkp(j,1)+xkm(j,1)+xip(j,1)+xim(j,1))           
     .          + eps1*(dy(j+jdk2)+seta1*dyj(j,1)*
     .            xj3p(j,1))    
         dvol(j+jdkdm3)=eps4*(xjp(j,2)*dy(j+1+jdkdm1)+xjm(j,2)*  
     .           dy(j-1+jdkdm1)+xkp(j,2)*dy(j+jdk)+xkm(j,2)*dy(j+jdk2)         
     .           +xip(j,2)*dy(j+jdm2+jdkdm1)+xim(j,2)*dy(j+jdm3+jdkdm1))   
     .          /(xjp(j,2)+xjm(j,2)+xkp(j,2)+xkm(j,2)+xip(j,2)+xim(j,2))           
     .          + eps3*(dy(j+jdk2)+seta2*dyj(j,2)*
     .            xj3p(j,2))    
416     continue 
        do 417 i=1,idm
          dvol(jdm*i       )       = dy(jdm*i+jdk1)
          dvol(1+jdm*(i-1)       ) = dy(1+jdm*(i-1)+jdk1)
          dvol(jdm*i+jdkdm3)       = dy(jdm*i+jdkdm1)
          dvol(1+jdm*(i-1)+jdkdm3) = dy(1+jdm*(i-1)+jdkdm1)
417     continue
        if(idm.gt.2) then
          jdmi2 = jdm*(idm-1)
          do 423 j  = 1,jdm
            dvol(j       )      = dy(j+jdk1)
            dvol(j+jdmi2       )= dy(j+jdmi2+jdk1)
            dvol(j+jdkdm3)      = dy(j+jdkdm1)
            dvol(j+jdmi2+jdkdm3)= dy(j+jdmi2+jdkdm1)
423       continue
        end if
        do 427 j=1,jdm*idm
          ydel = ydel+(dvol(j     )-dy(j+jdk1))*
     .                   (dvol(j     )-dy(j+jdk1))
     .                  +(dvol(j+jdkdm3)-dy(j+jdkdm1))*
     .                   (dvol(j+jdkdm3)-dy(j+jdkdm1))
427     continue
        if(jdk1.ne.jdkdm1) then
          do  j=1,jdm*idm                      
             dy(j+jdk1  ) = dvol(j       )
             dy(j+jdkdm1) = dvol(j+jdkdm3)
          enddo
        else
          do  j=1,jdm*idm                      
             dy(j+jdk1  ) = .5*(dvol(j     )+dvol(j+jdkdm3)) 
          enddo
        end if
445    continue 
       if(real(xdel/(real(jdm)*real(kdm)*real(idm))).lt.real(tol).and.
     1    real(ydel/(real(jdm)*real(kdm)*real(idm))).lt.real(tol).and.
     2    real(zdel/(real(jdm)*real(kdm)*real(idm))).lt.real(tol)) 
     3      goto  401  
1220   format(e12.5,2i7,3(1x,e12.5))   

400    continue
c      write(43,1220) real(time),nbl,ismooth,real(xdel)/(real(jdm)*
c    .            real(kdm)*real(idm)),
c    .               real(ydel)/(real(jdm)*real(kdm)*real(idm)),
c    .               real(zdel)/(real(jdm)*real(kdm)*real(idm))
       goto 402
401    continue
c      write(43,1220) real(time),nbl,iter,real(xdel)/(real(jdm)*
c    .            real(kdm)*real(idm)),
c    .               real(ydel)/(real(jdm)*real(kdm)*real(idm)),
c    .               real(zdel)/(real(jdm)*real(kdm)*real(idm))
402    continue
  
       if(i2d.ne.0) then 
        do 2000 k=1,kdm 
        do 2000 j=1,jdm 
          dx(j+jdm+jdm*idm*(k-1)) = dx(j+jdm*idm*(k-1))  
 2000     dz(j+jdm+jdm*idm*(k-1)) = dz(j+jdm*idm*(k-1))  
       end if 

      deallocate(xj1p)
      deallocate(xj1m)
      deallocate(xj2p)
      deallocate(xj2m)
      deallocate(xi1p)
      deallocate(xi1m)
      deallocate(xi2p)
      deallocate(xi2m)
      deallocate(xk1p)
      deallocate(xk1m)
      deallocate(xk2p)
      deallocate(xk2m)
      deallocate(xj3p)
      deallocate(seta)
      deallocate(ipid)
      deallocate(imid)
      deallocate(rj1jd)
      deallocate(rj2jd)
      deallocate(xip)
      deallocate(xjp)
      deallocate(xkp)
      deallocate(xim)
      deallocate(xjm)
      deallocate(xkm)
      deallocate(xid)
      deallocate(xjd)
      deallocate(dzj)
      deallocate(dxj)
      deallocate(dyj)
      deallocate(epsilon)

      call unperm_ijkxyz(jdmt,idmt,kdmt,nbl,dx,dy,dz,dvol)
      call unperm_ijkxyz(jdmt,idmt,kdmt,nbl,xnm1,ynm1,znm1,dvol)
      call unperm_ijkxyz(jdmt,idmt,kdmt,nbl,x,y,z,dvol)

       return
       end 
