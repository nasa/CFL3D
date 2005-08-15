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
c     isktyp  = 0 use only transfinite interpolation 
c             > 0 do isktyp smoothing steps using the spring analogy 
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
      common /deformz/ beta1,alpha1,isktyp,negvol,meshdef,nsprgit
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
     .       abs(ibctyp).ne.2016 .and.
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
     .       abs(ibctyp).ne.2016 .and.
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
     .          abs(jbctyp).ne.2016 .and.
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
     .          abs(jbctyp).ne.2016 .and.
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
     .          abs(kbctyp).ne.2016 .and.
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
     .          abs(kbctyp).ne.2016 .and.
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
      if (isktyp.gt.0) then
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
     .             isktyp)
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
     .                isktyp)
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
       do 400 iter = 1,isktyp                   
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
c      write(43,1220) real(time),nbl,isktyp,real(xdel)/(real(jdm)*
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
      subroutine deform_el(islavept,nslave,nmaster,u,xst,yst,zst,
     .                     nt,myhost,mycomm,myid,nnodes,
     .                     mblk2nd,nblelst,maxbl,iseqr)
c
c***********************************************************************
c     Purpose: Compute mesh deformation via Elastic Equations
c***********************************************************************
c

      integer stats
c
      dimension islavept(nslave,nmaster,5)
      dimension xst(nslave),yst(nslave),zst(nslave)
      dimension u(3*nslave)
      dimension mblk2nd(maxbl)
      dimension nblelst(maxbl,2)

      allocatable :: sa(:)
      allocatable :: stiffl(:,:)
      allocatable :: ija(:)
      allocatable :: b(:)
      allocatable :: ei(:)
      allocatable :: ripm(:,:)
      allocatable :: ej(:)
      allocatable :: ek(:)
      allocatable :: gij(:)
      allocatable :: gjk(:)
      allocatable :: gik(:)
      allocatable :: xix(:)
      allocatable :: xiy(:)
      allocatable :: xiz(:)
      allocatable :: etax(:)
      allocatable :: etay(:)
      allocatable :: etaz(:)
      allocatable :: zetax(:)
      allocatable :: zetay(:)
      allocatable :: zetaz(:)
      allocatable :: ooj(:)
      allocatable :: volij(:)
      allocatable :: volik(:)
      allocatable :: tinv(:,:)
      allocatable :: cjm1(:,:)
      allocatable :: cjp1(:,:)
      allocatable :: cim1(:,:)
      allocatable :: cip1(:,:)
      allocatable :: ckm1(:,:)
      allocatable :: ckp1(:,:)
      allocatable :: c(:,:)
      allocatable :: t(:,:)
      allocatable :: ul(:,:)
      allocatable :: n1(:)

c
      common /twod/ i2d
      common /zero/ iexp
      common /unst/ time,cfltau,ntstep,ita,iunst
      common /deformz/ beta1,alpha1,isktyp,negvol,meshdef,nsprgit

c
      memuse = 0
c
      allocate( sa(245*nslave+2), stat=stats )
      call umalloc(245*nslave+2,0,'sa',memuse,stats)
      allocate( stiffl(24,24), stat=stats )
      call umalloc(576,0,'stiffl',memuse,stats)
      allocate( ija(245*nslave+2), stat=stats )
      call umalloc(245*nslave+2,1,'ija',memuse,stats)
      allocate( b(3*nslave), stat=stats )
      call umalloc(3*nslave,0,'b',memuse,stats)
      allocate( ripm(7,nslave), stat=stats )
      call umalloc(7*nslave,0,'ripm',memuse,stats)
      allocate( ei(nslave), stat=stats )
      call umalloc(nslave,0,'ei',memuse,stats)
      allocate( ej(nslave), stat=stats )
      call umalloc(nslave,0,'ej',memuse,stats)
      allocate( ek(nslave), stat=stats )
      call umalloc(nslave,0,'ek',memuse,stats)
      allocate( gjk(nslave), stat=stats )
      call umalloc(nslave,0,'gjk',memuse,stats)
      allocate( gij(nslave), stat=stats )
      call umalloc(nslave,0,'gij',memuse,stats)
      allocate( gik(nslave), stat=stats )
      call umalloc(nslave,0,'gik',memuse,stats)
      allocate( xix(nslave), stat=stats )
      call umalloc(nslave,0,'xix',memuse,stats)
      allocate( xiy(nslave), stat=stats )
      call umalloc(nslave,0,'xiy',memuse,stats)
      allocate( xiz(nslave), stat=stats )
      call umalloc(nslave,0,'xiz',memuse,stats)
      allocate( etax(nslave), stat=stats )
      call umalloc(nslave,0,'etax',memuse,stats)
      allocate( etay(nslave), stat=stats )
      call umalloc(nslave,0,'etay',memuse,stats)
      allocate( etaz(nslave), stat=stats )
      call umalloc(nslave,0,'etaz',memuse,stats)
      allocate( zetax(nslave), stat=stats )
      call umalloc(nslave,0,'zetax',memuse,stats)
      allocate( zetay(nslave), stat=stats )
      call umalloc(nslave,0,'zetay',memuse,stats)
      allocate( zetaz(nslave), stat=stats )
      call umalloc(nslave,0,'zetaz',memuse,stats)
      allocate( ooj(nslave), stat=stats )
      call umalloc(nslave,0,'ooj',memuse,stats)
c     allocate( volij(nslave), stat=stats )
c     call umalloc(nslave,0,'volij',memuse,stats)
      allocate( volij(3*nslave), stat=stats )
      call umalloc(3*nslave,0,'volij',memuse,stats)
      allocate( volik(nslave), stat=stats )
      call umalloc(nslave,0,'volik',memuse,stats)
      allocate( tinv(6,6), stat=stats )
      call umalloc(36,0,'tinv',memuse,stats)
      allocate( cjm1(6,6), stat=stats )
      call umalloc(36,0,'cjm1',memuse,stats)
      allocate( cjp1(6,6), stat=stats )
      call umalloc(36,0,'cjp1',memuse,stats)
      allocate( cim1(6,6), stat=stats )
      call umalloc(36,0,'cim1',memuse,stats)
      allocate( cip1(6,6), stat=stats )
      call umalloc(36,0,'cip1',memuse,stats)
      allocate( ckm1(6,6), stat=stats )
      call umalloc(36,0,'ckm1',memuse,stats)
      allocate( ckp1(6,6), stat=stats )
      call umalloc(36,0,'ckp1',memuse,stats)
      allocate( c(6,6), stat=stats )
      call umalloc(36,0,'c',memuse,stats)
      allocate( t(6,6), stat=stats )
      call umalloc(36,0,'t',memuse,stats)
      allocate( ul(6,6), stat=stats )
      call umalloc(36,0,'ul',memuse,stats)
      allocate( n1(20), stat=stats )
      call umalloc(20,1,'n1',memuse,stats)

      eps  = 100.*10.**(-iexp)
c 
c     n       = i,j,k
c     islavept(2) = j - 1
c     islavept(3) = j + 1
c     islavept(4) = k - 1
c     islavept(5) = k + 1
c     islavept(6) = i - 1
c     islavept(7) = i + 1
c     islavept(10) = number of ija and sa array values to reserve for this
c                    point
c     islavept(11) = max number of coincident points at block interfaces
c     islavept(12) = coincident pt 2 (pt n is the first)   
c     islavept(13) = coincident pt 3
c     islavept(14) = coincident pt 4
c     islavept(15) = coincident pt 5
c     islavept(16) = coincident pt 6
c     islavept(17) = coincident pt 7
c     islavept(18) = coincident pt 8
c     islavept(19) = coincident pt 9
c     islavept(20) = coincident pt 10
c     islavept(21) = number of the nearest surface pt
c
      ipt = 3*nslave+2 
      ija(1) = ipt
      do n = 1,nslave 
         ija(3*(n-1)+2)= ija(3*(n-1)+1) + islavept(n,10,iseqr)
         ija(3*(n-1)+3)= ija(3*(n-1)+2) + islavept(n,10,iseqr) 
         ija(3*n    +1)= ija(3*(n-1)+3) + islavept(n,10,iseqr) 
      enddo
      do n = 3*nslave+2,245*nslave+2
        ija(n) = 0
      enddo 
c
c
       if(nt.eq.1.and.meshdef.eq.1) then
          do n = 1,nslave
           iimax = islavept(n,11,iseqr)
           if(iimax.gt.1) then
                do ii2 = 2,iimax
                  n2 = islavept(n,12+ii2-2,iseqr)
                  test  = (xst(n2)-xst(n))*(xst(n2)-xst(n))
     .                   +(yst(n2)-yst(n))*(yst(n2)-yst(n))
     .                   +(zst(n2)-zst(n))*(zst(n2)-zst(n))
                  test = sqrt(test)
                  if(test.gt.100.*eps) then
                   write(11,29102) n,n2,test
29102              format(' WARNING: Macro-elements ',i6,' and ',i6,
     .           ' at 1-1 blocking interface have geometric mismatch'
     .              ,' = ',e16.8)
                  end if
                enddo
           end if
          enddo
       end if

c
      if(isktyp.gt.0) then 
        gini = 10.d0 
        eini = 10.d0 

        mxiter1 = 1      
c
        call elrhs(b,u,islavept,nslave,nmaster,iseqr) 

        do itr1 = 1,mxiter1 

          do n = 1,245*nslave+2
            sa(n) = 0.d0 
          enddo


          if(meshdef.eq.1) then
           if(myid.eq.myhost) then
             write(11,21022) nslave 
21022        format(' Calculating macroelement moduli ',
     .               'for ',i8,' element nodes')
           end if
          end if
          call hookefe(ei,ej,ek,gij,gjk,gik,xst,yst,zst,volij,
     .                 volik,eps,eini,gini,arg1,arg2,islavept,
     .                 nslave,nmaster,nt,nnodes,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,iseqr)

          if(meshdef.eq.1) then
           if(myid.eq.myhost) then
             write(11,21023) 
21023        format(' Calculating macroelement metrics ')
           end if
          end if

          call elmetricsfe(xix,xiy,xiz,etax,etay,etaz,zetax,zetay,
     .                     zetaz,ooj,xst,yst,zst,eps,islavept,nslave,
     .                     nmaster,nnodes,myhost,myid,mycomm,
     .                     mblk2nd,maxbl,iseqr)
 
          if(meshdef.eq.1) then
           if(myid.eq.myhost) then
             write(11,21024) 
21024        format(' Calculating and assembling macroelement stiffness'
     .                ,'  matrix ')
           end if
          end if

          call elglobfe(sa,xst,yst,zst,xix,xiy,xiz,etax,etay,
     .                  etaz,zetax,zetay,zetaz,ei,ej,ek,gij,gjk,
     .                  gik,ooj,eps,stiffl,islavept,ija,nslave,
     .                  nmaster,nnodes,myhost,myid,mycomm,mblk2nd,
     .                  nblelst,maxbl,iseqr)

          if(meshdef.eq.1) then
           if(myid.eq.myhost) then
             write(11,21025) 
21025        format(' Solving the macroelement system ')
           end if
          end if
         
          ncount = 245*nslave+2
          ndim   = 3*nslave
c
          call gaussseidel(ncount,ndim,b,u,sa,ija,1e-6,9*nslave,err
     .                     ,volij,myid,myhost) 
c         call dprec(sa,ija,b,ndim,ncount)

c         call linbcg(ncount,ndim,b,u,sa,ija,1,1e-6,nslave,iter,err
c    .                ,time,eps,myid) 
          
          call coincdef(volij,u,islavept,nslave,nmaster,iseqr,n1)
        enddo
      else

          if(meshdef.eq.1) then
           if(myid.eq.myhost) then
             write(11,21027) 
21027        format(' Solving nodal displacements using exponential'
     .              ,' decay  ')
           end if
          end if
         
          call expdecay(xst,yst,zst,u,volij,ei,ripm,eps,
     .                  islavept,nslave,nmaster,ndim,nt,
     .                  nnodes,myhost,myid,mycomm,mblk2nd,
     .                  maxbl,iseqr)
          call coincdef(volij,u,islavept,nslave,nmaster,iseqr,n1)
      end if

        if(i2d.ne.0) then 
          do n = 1,nslave
            if(islavept(n,6,iseqr).eq.n) then
              nip1     = islavept(n,7,iseqr)
              u(3*(n-1)+1)    = .5*(u(3*(n-1)+1) + u(3*(nip1-1)+1))
              u(3*(nip1-1)+1) = u(3*(n-1)+1)
              u(3*(nip1-1)+2) = 0.
              u(3*(n-1)+2)    = 0. 
              u(3*(n-1)+3)    = .5*(u(3*(n-1)+3) + u(3*(nip1-1)+3))
              u(3*(nip1-1)+3) = u(3*(n-1)+3) 
            end if
          enddo
        end if
c
c     release memory
c
      deallocate(sa)
      deallocate(stiffl)
      deallocate(ija)
      deallocate(b)
      deallocate(ei)
      deallocate(ripm)
      deallocate(ej)
      deallocate(ek)
      deallocate(gij)
      deallocate(gjk)
      deallocate(gik)
      deallocate(xix)
      deallocate(xiy)
      deallocate(xiz)
      deallocate(etax)
      deallocate(etay)
      deallocate(etaz)
      deallocate(zetax)
      deallocate(zetay)
      deallocate(zetaz)
      deallocate(ooj)
      deallocate(volij)
      deallocate(volik)
      deallocate(tinv)
      deallocate(cjm1)
      deallocate(cjp1)
      deallocate(ckm1)
      deallocate(ckp1)
      deallocate(cim1)
      deallocate(cip1)
      deallocate(c)
      deallocate(t)
      deallocate(ul)
      deallocate(n1)


      return
      end
c
      subroutine coincdef(volij,u,islavept,nslave,nmaster,iseqr,n1)
      dimension islavept(nslave,nmaster,5)
      dimension u(3*nslave),volij(3*nslave)
      dimension n1(20) 
c
          volij = 0.
          do n = 1,nslave
           iimax = islavept(n,11,iseqr)
           if(iimax.gt.1) then
              n1(1) = n
              if(islavept(n,8,iseqr).ne.0) then
                ii3 = 0
                do ii2 = 2,iimax
                  n1(ii2) = islavept(n,12+ii2-2,iseqr)
                  if(islavept(n1(ii2),8,iseqr).eq.0) ii3 = ii2
                enddo
                if(ii3.eq.0) then
                  rooiim = 1./real(iimax)
                  do ii2 = 1,iimax
                    do j = 1,3
                      volij(3*(n-1)+j) = volij(3*(n-1)+j) + 
     .                                   u(3*(n1(ii2)-1)+j)*rooiim
                    enddo
                  enddo
                else
                  do j = 1,3
                    volij(3*(n-1)+j) = u(3*(n1(ii3)-1)+j)
                  enddo
                end if
              else
                 do j = 1,3
                   volij(3*(n-1)+j) = u(3*(n-1)+j)
                 enddo
              end if
           end if
          enddo
          do n = 1,nslave
           iimax = islavept(n,11,iseqr)
           if(iimax.gt.1) then
              n1(1) = n
              do ii2 = 2,iimax
                n1(ii2) = islavept(n,12+ii2-2,iseqr)
              enddo
              do ii2 = 1,iimax
                do j = 1,3
                  u(3*(n1(ii2)-1)+j) = volij(3*(n-1)+j)
                enddo
              enddo
           end if
          enddo

      return
      end
c
      SUBROUTINE gaussseidel(nc,n,b,x,sa,ija,tol,itmax,err, 
     .                       xt,myid,myhost) 
      parameter(EPS=1.d-12)  
      INTEGER iter,itmax,itol,n 
      INTEGER j 
      dimension b(n),x(n)
      dimension xt(n)
      dimension ija(nc)
      dimension sa(nc)
      common /deformz/ beta1,alpha1,isktyp,negvol,meshdef,nsprgit

      it1 = 0
      itmx1 =200
      err   = 1.e+10
      do it = 1,itmax
        sum = 0.
        resid = 0.
        do i = 1,n
          ut = b(i)/sa(i)  
          do j = ija(i),ija(i+1)-1
            ut = ut - sa(j)*x(ija(j))/sa(i)
          enddo
          xt(i) = beta1*ut + (1.-beta1)*x(i) 
          sum = sum + (xt(i)-x(i))*(xt(i)-x(i))
          resid1= sqrt((xt(i)-x(i))*(xt(i)-x(i)))
          if(resid1.gt.resid) then
           resid = resid1
           nresmx= i
          end if
          x(i) = xt(i)
        enddo
        errold = err
        err = sqrt(sum)
        nresmx = 1+(nresmx-1)/3
        if(meshdef.eq.1.and.myid.eq.myhost) 
     .               write(197+myid,31029) it,err,nresmx,resid   
31029   format(i8,e16.8,i8,e16.8)
        if(err.lt.tol) goto 2000 
        if(errold.lt.err) then
          it1 = it1 + 1
          if(it1.gt.itmx1) then
             write(197,21022) 
             stop
          end if
21022     format(' Stopping, Gauss-Seidel scheme not converging: '
     .        ,/,' Reduce BETA1 or change Macro-element definition')
        else
          it1 = 0
        end if
      enddo
2000  continue 
      if(myid.eq.myhost) write(197+myid,31029) it,err   
 
      return
      end
      SUBROUTINE linbcg(nc,n,b,x,sa,ija,itol,tol,itmax,iter,err
     .                  ,time,eps,myid) 
      INTEGER iter,itmax,itol,n 
      INTEGER j 
      dimension b(n),x(n)
      dimension ija(nc)
      dimension sa(nc)
      dimension p(n),pp(n),r(n),rr(n),z(n),zz(n)
CU    USES atimes,asolve,snrm 
      iter=0 

      test = 0.
      do 5 j = 1,n
        if(test.lt.abs(b(j))) test = abs(b(j))
5     continue
      if(test.le.eps) return
      call atimes(n,nc,x,r,sa,ija,0) 
      do 11 j=1,n 
        r(j)=b(j)-r(j) 
        rr(j)=r(j) 
11    continue 
      call atimes(n,nc,r,rr,sa,ija,0) 
      znrm=1.d0 
      if(itol.eq.1) then 
        bnrm=snrm(n,b,itol) 
      else if (itol.eq.2) then 
        call asolve(n,nc,b,z,sa,0)
        bnrm=snrm(n,z,itol) 
      else if (itol.eq.3.or.itol.eq.4) then 
        call asolve(n,nc,b,z,sa,0)
        bnrm=snrm(n,z,itol) 
        call asolve(n,nc,r,z,sa,0)
        znrm=snrm(n,z,itol) 
      else 
        write(197+myid,31922)
31922   format('illegal itol in linbcg') 
        stop
      endif 
      if(bnrm.lt.EPS) then
        err = 0.
        znrm= 0.
        goto 200 
      end if
      call asolve(n,nc,r,z,sa,0)
100   if (iter.le.itmax) then 
        iter=iter+1 
        call asolve(n,nc,rr,zz,sa,1)
        bknum=0.d0 
        do 12 j=1,n 
          bknum=bknum+z(j)*rr(j) 
12      continue 
        if(iter.eq.1) then 
          do 13 j=1,n 
            p(j)=z(j) 
            pp(j)=zz(j) 
13        continue 
        else 
          bk=bknum/bkden 
          do 14 j=1,n 
            p(j)=bk*p(j)+z(j) 
            pp(j)=bk*pp(j)+zz(j) 
14        continue 
        endif 
        bkden=bknum 
        call atimes(n,nc,p,z,sa,ija,0) 
        akden=0.d0 
        do 15 j=1,n 
          akden=akden+z(j)*pp(j) 
15      continue 
        ak=bknum/akden 
        call atimes(n,nc,pp,zz,sa,ija,1) 
c
c
        do 16 j=1,n 
          x(j)=x(j)+ak*p(j) 
          r(j)=r(j)-ak*z(j) 
          rr(j)=rr(j)-ak*zz(j) 
16      continue 
c
        
        call asolve(n,nc,r,z,sa,0)
        if(itol.eq.1.) then
          if(bnrm.gt.EPS) then
            err=snrm(n,r,itol)/bnrm 
          else
            err = 0.
          end if
             write(197+myid,5040) iter,log10(err),znrm  
        else if(itol.eq.2) then 
          if(bnrm.gt.EPS) then
            err=snrm(n,z,itol)/bnrm 
          else
            err = 0.
          end if
             write(197+myid,5040) iter,log10(err),znrm  
        else if(itol.eq.3.or.itol.eq.4)then 
          zm1nrm=znrm 
          znrm=snrm(n,z,itol)  
c
          if(abs(zm1nrm-znrm).gt.EPS*znrm) then 
            dxnrm=abs(ak)*snrm(n,p,itol) 
            if(abs(zm1nrm-znrm).gt.EPS) then
              err=znrm/abs(zm1nrm-znrm)*dxnrm 
            else
              err = 0.
            end if
          else
            if(bnrm.gt.EPS) then
              err=znrm/bnrm 
            else
              err = 0.
              znrm= 0. 
            end if
c            write(197+myid,5042) iter,log10(err),znrm  
            goto 100 
          endif 
          xnrm=snrm(n,x,itol) 
c
          if(err.le.0.5d0*xnrm) then 
              if(xnrm.gt.EPS) then
                err=err/xnrm 
              else
                err = 0.
              end if
             write(197+myid,5043) iter,log10(err),znrm  
          else 
            if(bnrm.gt.EPS) then
              err=znrm/bnrm 
            else
              err = 0.
              znrm= 0. 
            end if
c            write(197+myid,5044) iter,log10(err),znrm  
            goto 100 
          endif 
        endif 
      if(err.gt.tol) goto 100 
      endif 
200   continue
5040    format('1',i8,2x,2(1x,e12.5)) 
5041    format('2',i8,2x,6(1x,e12.5)) 
5042    format('3',i8,2x,2(1x,e12.5)) 
5043    format('4',i8,2x,2(1x,e12.5)) 
5044    format('5',i8,2x,2(1x,e12.5)) 

      return 
      END 


C  (C) Copr. 1986-92 Numerical Recipes Software ">u,3. 



      SUBROUTINE atimes(n,nc,x,r,sa,ija,itrnsp) 
      INTEGER n,itrnsp
      dimension ija(nc)
      dimension x(n),r(n),sa(nc) 
CU    USES dsprsax,dsprstx 
      if (itrnsp.eq.0) then 
        call dsprsax(sa,ija,x,r,n,nc) 
      else 
        call dsprstx(sa,ija,x,r,n,nc) 
      endif 
      return 
      END 
C  (C) Copr. 1986-92 Numerical Recipes Software ">u,3.

      SUBROUTINE dsprstx(sa,ija,x,b,n,nc)
      INTEGER n,ija(nc)
      dimension b(n),sa(nc),x(n)
      INTEGER i,j,k
      if (ija(1).ne.n+2) then
        write(197+myid,3111) ija(1),n+2,nc 
3111    format('mismatched vector and matrix in sprstx',3i8) 
        stop
      end if 
      do 11 i=1,n
        b(i)=sa(i)*x(i)
11    continue
      do 13 i=1,n
        do 12 k=ija(i),ija(i+1)-1
          j=ija(k)
          b(j)=b(j)+sa(k)*x(i)
12      continue
13    continue
c
      return
      END
c
c  Bartels preconditioner
c
      SUBROUTINE dprec(sa,ija,b,n,nc)
      INTEGER ija(nc)
      dimension sa(nc),b(n) 
      INTEGER i,k
      if (ija(1).ne.n+2) then
        write(197+myid,3111) ija(1),n+2,nc 
3111    format('mismatched vector and matrix in sprsax',3i8) 
        stop
      end if 
      do 12 i=1,n
        b(i) = b(i)/sa(i)
        do 11 k=ija(i),ija(i+1)-1
          sa(k)=sa(k)/sa(i)
11      continue
        sa(i)=1.
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ">u,3.

      SUBROUTINE dsprsax(sa,ija,x,b,n,nc)
      INTEGER n,ija(nc)
      dimension b(n),sa(nc),x(n)
      INTEGER i,k
      if (ija(1).ne.n+2) then
        write(197+myid,3111) ija(1),n+2,nc 
3111    format('mismatched vector and matrix in sprsax',3i8) 
        stop
      end if 
      do 12 i=1,n
        b(i)=sa(i)*x(i)
        do 11 k=ija(i),ija(i+1)-1
          b(i)=b(i)+sa(k)*x(ija(k))
11      continue
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ">u,3.

      SUBROUTINE asolve(n,nc,b,x,sa,itrnsp) 
      INTEGER n,itrnsp,i 
      dimension x(n),b(n),sa(nc)
      do 11 i=1,n 
        x(i)=b(i)/sa(i) 
11    continue 

      return 
      END 
C  (C) Copr. 1986-92 Numerical Recipes Software ">u,3.

      FUNCTION snrm(n,sx,itol) 
      INTEGER n,itol,i,isamax 
      dimension sx(n)
      if (itol.le.3)then 
        snrm=0. 
        do 11 i=1,n 
          snrm=snrm+sx(i)**2 
11      continue 
        snrm=sqrt(snrm) 
      else 
        isamax=1 
        do 12 i=1,n 
          if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i 
12      continue 
        snrm=abs(sx(isamax)) 
      endif 
      return 
      END 
c
************************************************************************
*
*         SUBROUTINE FOR THE INVERSION OF COMIN USING GAUSS-JORDAN
*         WITH PIVOTING.
*
************************************************************************
      SUBROUTINE INVDET(UL,COMIN,N,DTNRM,DETM)
      DIMENSION  UL(N,N),COMIN(N,N)
      DIMENSION C(N,N),J(550)
      DO 10 K = 1,n
        DO 12 I = 1,n
          C(I,K) = COMIN(I,K)
12      CONTINUE
10    CONTINUE
      PD = 1.D0
      DO 124 L = 1,N
        DD = 0.D0
        DO 123 K = 1,N
123     DD = DD + C(L,K)*C(L,K)
        DD = SQRT(DD)
124   PD = PD*DD
      DETM = 1.D0
      DO 125 L = 1,N
125   J(L+20) = L
      DO 144 L = 1,N
        CC = 0.D0
        M  = L
        DO 135 K = L,N
          IF((ABS(CC)-ABS(C(L,K))).GE.0.D0) GOTO 135
126       M = K
          CC = C(L,K)
135     CONTINUE
127     IF(L.EQ.M) GOTO 138
128     K = J(M+20)
        J(M+20) = J(L+20)
        J(L+20) = K
        DO 137 K = 1,N
          S = C(K,L)
          C(K,L) = C(K,M)
137     C(K,M) = S
138     C(L,L) = 1.D0
        DETM = DETM*CC
        DO 139 M = 1,N
139     C(L,M) = C(L,M)/CC
        DO 142 M = 1,N
          IF(L.EQ.M) GOTO 142
129       CC = C(M,L)
          IF(CC.EQ.0.D0) GOTO 142
130       C(M,L) = 0.D0
          DO 141 K = 1,N
141       C(M,K) = C(M,K) - CC*C(L,K)
142     CONTINUE
144   CONTINUE
      DO 143 L = 1,N
        IF(J(L+20).EQ.L) GOTO 143
131     M = L
132     M = M + 1
        IF(J(M+20).EQ.L) GOTO 133
136     IF(N.GT.M) GOTO 132
133     J(M+20) = J(L+20)
        DO 163 K = 1,N
          CC = C(L,K)
          C(L,K) = C(M,K)
163     C(M,K) = CC
        J(L+20) = L
143   CONTINUE
      DETM = ABS(DETM)
      DTNRM = DETM/PD
      DO 210 k = 1,N
        DO 212 I = 1,N
          UL(I,K) = C(I,K)
212     CONTINUE
210   CONTINUE
      RETURN
      END
      subroutine expdecay(xst,yst,zst,u,ut,r,ripm,eps,
     .                    islavept,nslave,nmaster,ndim,nt,
     .                    nnodes,myhost,myid,mycomm,mblk2nd,
     .                    maxbl,iseqr)
      dimension islavept(nslave,nmaster,5)
      dimension xst(nslave),yst(nslave),zst(nslave)
      dimension ut(ndim),u(ndim)
      dimension r(nslave),ripm(7,nslave)
      dimension mblk2nd(maxbl)
      dimension arg(7),coef(8),dm(7),rsn(7)
      dimension nipm(7),nsf(7),n3(7)
      common /deformz/ beta1,alpha1,isktyp,negvol,meshdef,nsprgit
      common /sgminmax/ rmin,rmax

      relax = 1.0
c
c    Do this at every time step because the closest surface
c    point to a control point will potentially change as a
c    surface moves.
c
        do n1 = 1,nslave
          r(n1) = 1.e6 
          if(islavept(n1,8,iseqr).ne.0) then
           do n2 = 1,nslave
            if(islavept(n2,8,iseqr).eq.0) then
              test  = (xst(n2)-xst(n1))*(xst(n2)-xst(n1))
     .               +(yst(n2)-yst(n1))*(yst(n2)-yst(n1))
     .               +(zst(n2)-zst(n1))*(zst(n2)-zst(n1))
              test  = sqrt(test)
              if(r(n1).gt.test) then
                r(n1)             = test
                islavept(n1,21,iseqr) = n2
              end if
            end if
           enddo
          else
           r(n1) = 0.
           islavept(n1,21,iseqr) = n1
          end if
        enddo
        if(nt.eq.1) then
          rmax = 0. 
          rmin = 1.e+30  
          do n = 1,nslave
            rave =  r(n)
            if(rave.gt.rmax) rmax = rave   
            if(rave.lt.rmin.and.rave.gt.0.) rmin = rave  
          enddo
        end if
        eps = rmin/10000.
        do n = 1,nslave
          nsf(1)  = islavept(n,8,iseqr)
          if(nsf(1).ne.0) then
            nipm(1)= n
            nipm(2)= islavept(n,2,iseqr)
            nipm(3)= islavept(n,3,iseqr)
            nipm(4)= islavept(n,4,iseqr)
            nipm(5)= islavept(n,5,iseqr)
            nipm(6)= islavept(n,6,iseqr)
            nipm(7)= islavept(n,7,iseqr)
            do ii = 1,7
              nsf(ii)  = islavept(nipm(ii),8,iseqr)
              if(n.ne.nipm(ii)) then 
                ripm(ii,n)=(xst(n)-xst(nipm(ii)))*(xst(n)-xst(nipm(ii)))
     .                    +(yst(n)-yst(nipm(ii)))*(yst(n)-yst(nipm(ii)))
     .                    +(zst(n)-zst(nipm(ii)))*(zst(n)-zst(nipm(ii)))
     .                     +eps 
                ripm(ii,n) = sqrt(ripm(ii,n))
              else
                ripm(ii,n) = 1.e+30
              end if  
            enddo
            coef   = 0.
            coef(8)= 1.
            coef(1)= 1.
            do ii = 2,7
              if(n.ne.nipm(ii)) then
                coef(ii)= exp(-100.*ripm(ii,n)/rmax)
                coef(8)= coef(8) + coef(ii)
              end if
            enddo
            n3(1) = islavept(n,21,iseqr)
            rsn(1)= r(n)
            do ii = 2,7
              n3(ii) = islavept(nipm(ii),21,iseqr)
              rsn(ii)= r(nipm(ii))
            enddo
            n2    = n3(1)
            dm(1) = sqrt(u(3*(n2-1)+1)*u(3*(n2-1)+1) + 
     .                   u(3*(n2-1)+2)*u(3*(n2-1)+2) + 
     .                   u(3*(n2-1)+3)*u(3*(n2-1)+3))
            arg(1) = - beta1*(rsn(1)/rmax  - alpha1) 
            if(arg(1).gt.0.) arg(1) = 0.
            coef(1)= coef(1)/coef(8) 
            do ii = 2,7
              if(n.ne.nipm(ii)) then 
                n2    = n3(ii)
                dm(ii) = sqrt(u(3*(n2-1)+1)*u(3*(n2-1)+1) + 
     .                        u(3*(n2-1)+2)*u(3*(n2-1)+2) + 
     .                        u(3*(n2-1)+3)*u(3*(n2-1)+3))
                arg(ii) = - beta1*(rsn(ii)/rmax  - alpha1) 
                if(arg(ii).gt.0.) arg(ii) = 0.
                coef(ii)= coef(ii)/coef(8) 
              else
                coef(ii) = 0.
                arg(ii)  = 0.
              end if
            enddo
            do j = 1,3
              u(3*(n-1)+j) = 0.
              do ii = 1,7
                n2 = n3(ii)
                u(3*(n-1)+j) = u(3*(n-1)+j) + 
     .                      u(3*(n2-1)+j)*coef(ii)*exp(arg(ii))
              enddo
            enddo
          end if
        enddo
       do it = 1,nsprgit   
        do n = 1,nslave
          nipm(2)= islavept(n,2,iseqr)
          nipm(3)= islavept(n,3,iseqr)
          nipm(4)= islavept(n,4,iseqr)
          nipm(5)= islavept(n,5,iseqr)
          nipm(6)= islavept(n,6,iseqr)
          nipm(7)= islavept(n,7,iseqr)
          nsf(1) = islavept(n,8,iseqr)
c         do ii = 2,7
c           nsf(ii) = islavept(nipm(ii),8,iseqr)
c         enddo
          if(nsf(1).ne.0) then
            coef1 = 0. 
            do ii = 2,7
              coef1= coef1 + 1./ripm(ii,n)
            enddo
            do j = 1,3
              ut(3*(n-1)+j) = 0.
              do ii = 2,7
               ut(3*(n-1)+j)= ut(3*(n-1)+j)
     .             + u(3*(nipm(ii)-1)+j)/ripm(ii,n)/coef1   
              enddo
            enddo
          end if
        enddo
        do n = 1,nslave
          nipm(2)= islavept(n,2,iseqr)
          nipm(3)= islavept(n,3,iseqr)
          nipm(4)= islavept(n,4,iseqr)
          nipm(5)= islavept(n,5,iseqr)
          nipm(6)= islavept(n,6,iseqr)
          nipm(7)= islavept(n,7,iseqr)
          nsf(1) = islavept(n,8,iseqr)
c         do ii = 2,7
c           nsf(ii) = islavept(nipm(ii),8,iseqr)
c         enddo
          if(nsf(1).ne.0) then
            u(3*(n-1)+1) = relax*ut(3*(n-1)+1)+(1.-relax)*u(3*(n-1)+1) 
            u(3*(n-1)+2) = relax*ut(3*(n-1)+2)+(1.-relax)*u(3*(n-1)+2) 
            u(3*(n-1)+3) = relax*ut(3*(n-1)+3)+(1.-relax)*u(3*(n-1)+3) 
          end if
        enddo
       enddo

      return
      end

      subroutine hookefe(ei,ej,ek,gij,gjk,gik,xst,yst,zst,volij,
     .                   volik,eps,eini,gini,arg1,arg2,islavept,
     .                   nslave,nmaster,nt,nnodes,myhost,myid,
     .                   mycomm,mblk2nd,maxbl,iseqr)
      dimension islavept(nslave,nmaster,5)
      dimension xst(nslave),yst(nslave),zst(nslave)
      dimension ei(nslave),ej(nslave),ek(nslave),gij(nslave),
     .          gjk(nslave),gik(nslave)
      dimension volik(nslave),volij(3*nslave)
      dimension r(nslave),r1(0:8) 
      dimension ni1(0:8),ni0(0:8),test1(0:8) 
      dimension mblk2nd(maxbl)
      common /sgminmax/ rmin,rmax

      beta2 = 1.000 
c
c    Do this at every time step because the closest surface
c    point to a control point will potentially change as a
c    surface moves.
c
c       if(nt.eq.1) then
          rmax = 0. 
          rmin = 1.e+30  
c       end if  
        do n1 = 1,nslave
          r(n1) = 1.e6 
          ni0(0)    = n1
          ni0(1)    = islavept(n1    ,3,iseqr)
          ni0(2)    = islavept(n1    ,5,iseqr)
          ni0(3)    = islavept(ni0(1),5,iseqr)
          ni0(4)    = islavept(n1    ,7,iseqr)
          ni0(5)    = islavept(ni0(2),7,iseqr)
          ni0(6)    = islavept(ni0(1),7,iseqr)
          ni0(7)    = islavept(ni0(6),5,iseqr)
          if(ni0(1).ne.n1.and.ni0(2).ne.n1.and.ni0(4).ne.n1) then
           xcent = 0.
           ycent = 0.
           zcent = 0.
           do ii = 0,7
             xcent = xcent + xst(ni0(ii))/8. 
             ycent = ycent + yst(ni0(ii))/8. 
             zcent = zcent + zst(ni0(ii))/8. 
           enddo
           do n2 = 1,nslave
            if(islavept(n2,8,iseqr).eq.0) then
              ni1(0)    = n2
              ni1(1)    = islavept(n2    ,3,iseqr)
              ni1(2)    = islavept(n2    ,5,iseqr)
              ni1(3)    = islavept(ni1(1),5,iseqr)
              ni1(4)    = islavept(n2    ,7,iseqr)
              ni1(5)    = islavept(ni1(2),7,iseqr)
              ni1(6)    = islavept(ni1(1),7,iseqr)
              ni1(7)    = islavept(ni1(6),5,iseqr)
              test  = (xst(n2)-xcent)*(xst(n2)-xcent)
     .               +(yst(n2)-ycent)*(yst(n2)-ycent)
     .               +(zst(n2)-zcent)*(zst(n2)-zcent)
              test    = sqrt(test)
              if(r(n1).gt.test) then
                r(n1)                 = test   
                islavept(n1,21,iseqr) = n2
              end if
              if(islavept(ni1(1),8,iseqr).eq.0.and.
     .           islavept(ni1(2),8,iseqr).eq.0) then 
                 xsct = (xst(n2)+xst(ni1(1))+xst(ni1(2))
     .                    +xst(ni1(3)))/4.
                 ysct = (yst(n2)+yst(ni1(1))+yst(ni1(2))
     .                    +yst(ni1(3)))/4.
                 zsct = (zst(n2)+zst(ni1(1))+zst(ni1(2))
     .                    +zst(ni1(3)))/4.
                 test  = (xsct-xcent)*(xsct-xcent)
     .                  +(ysct-ycent)*(ysct-ycent)
     .                  +(zsct-zcent)*(zsct-zcent)
                 test    = sqrt(test)
                 if(r(n1).gt.test) then
                    r(n1)                 = test   
                    islavept(n1,21,iseqr) = n2
                 end if
              end if
              if(islavept(ni1(2),8,iseqr).eq.0.and.
     .           islavept(ni1(4),8,iseqr).eq.0) then 
                 xsct = (xst(n2)+xst(ni1(2))+xst(ni1(4))
     .                    +xst(ni1(5)))/4.
                 ysct = (yst(n2)+yst(ni1(2))+yst(ni1(4))
     .                    +yst(ni1(5)))/4.
                 zsct = (zst(n2)+zst(ni1(2))+zst(ni1(4))
     .                    +zst(ni1(5)))/4.
                 test  = (xsct-xcent)*(xsct-xcent)
     .                  +(ysct-ycent)*(ysct-ycent)
     .                  +(zsct-zcent)*(zsct-zcent)
                 test    = sqrt(test)
                 if(r(n1).gt.test) then
                    r(n1)                 = test   
                    islavept(n1,21,iseqr) = n2
                 end if
              end if
              if(islavept(ni1(1),8,iseqr).eq.0.and.
     .           islavept(ni1(4),8,iseqr).eq.0) then 
                 xsct = (xst(n2)+xst(ni1(1))+xst(ni1(4))
     .                    +xst(ni1(6)))/4.
                 ysct = (yst(n2)+yst(ni1(1))+yst(ni1(4))
     .                    +yst(ni1(6)))/4.
                 zsct = (zst(n2)+zst(ni1(1))+zst(ni1(4))
     .                    +zst(ni1(6)))/4.
                 test  = (xsct-xcent)*(xsct-xcent)
     .                  +(ysct-ycent)*(ysct-ycent)
     .                  +(zsct-zcent)*(zsct-zcent)
                 test    = sqrt(test)
                 if(r(n1).gt.test) then
                    r(n1)                 = test   
                    islavept(n1,21,iseqr) = n2
                 end if
              end if
          
            end if
           enddo
c          if(nt.eq.1) then
             n2    = islavept(n1,21,iseqr)
             if(r(n1).gt.rmax) rmax = r(n1)    
             if(r(n1).lt.rmin.and.r(n1).gt.0.) rmin = r(n1)   
c          end if
          end if 
        enddo

         
 
      e0 = 500000.
      r0 = (rmax*alog(1.-eini/e0)+beta2*rmin)/
     .     (     alog(1.-eini/e0)+beta2     ) 
  
      do n = 1,nslave
       nbl = islavept(n,9,iseqr)
        ni0(1)    = islavept(n     ,3,iseqr)
        ni0(2)    = islavept(n     ,5,iseqr)
        ni0(4)    = islavept(n     ,7,iseqr)
        if(ni0(1).ne.n.and.ni0(2).ne.n.and.ni0(4).ne.n) then
          rave = r(n)
          fact     = 1./(1.-exp(-beta2*(rave-r0)/(rmax-r0))) 
c         fact     = 1.
          ej(n) = eini*fact
          ek(n) = eini*fact
          ei(n) = eini*fact
          gij(n)= gini*fact
          gik(n)= gini*fact
          gjk(n)= gini*fact
        end if
      enddo
      return
      end

      subroutine elmetricsfe(xix,xiy,xiz,etax,etay,etaz,zetax,zetay,
     .                       zetaz,ooj,xs,ys,zs,eps,islavept,nslave,
     .                       nmaster,nnodes,myhost,myid,mycomm,
     .                       mblk2nd,maxbl,iseqr)

      dimension islavept(nslave,nmaster,5)
      dimension xs(nslave),ys(nslave),zs(nslave)
      dimension xix(nslave),xiy(nslave),xiz(nslave),
     .          etax(nslave),etay(nslave),etaz(nslave),
     .          zetax(nslave),zetay(nslave),zetaz(nslave),
     .          ooj(nslave)
      dimension mblk2nd(maxbl)

c
      do n = 1,nslave
       nbl = islavept(n,9,iseqr)
c      if (myid .eq. mblk2nd(nbl)) then 
        njp    = islavept(n,3,iseqr)
        nkp    = islavept(n,5,iseqr)
        njpkp  = islavept(nkp,3,iseqr)
        nip    = islavept(n,7,iseqr)
        nipkp  = islavept(nkp,7,iseqr)
        njpip  = islavept(njp,7,iseqr)
        njpipkp= islavept(njpip,5,iseqr)
        if(njp.ne.n.and.nip.ne.n.and.nkp.ne.n) then
          xxi    =.25*(xs(njp)  -xs(n)  +xs(njpkp)  -xs(nkp)
     .                +xs(njpip)-xs(nip)+xs(njpipkp)-xs(nipkp))
          xeta   =.25*(xs(nkp)  -xs(n)  +xs(njpkp)  -xs(njp)
     .                +xs(nipkp)-xs(nip)+xs(njpipkp)-xs(njpip))
          xzeta  =.25*(xs(nip)  -xs(n)  +xs(nipkp)  -xs(nkp)
     .                +xs(njpip)-xs(njp)+xs(njpipkp)-xs(njpkp))
          yxi    =.25*(ys(njp)  -ys(n)  +ys(njpkp)  -ys(nkp)
     .                +ys(njpip)-ys(nip)+ys(njpipkp)-ys(nipkp))
          yeta   =.25*(ys(nkp)  -ys(n)  +ys(njpkp)  -ys(njp)
     .                +ys(nipkp)-ys(nip)+ys(njpipkp)-ys(njpip))
          yzeta  =.25*(ys(nip)  -ys(n)  +ys(nipkp)  -ys(nkp)
     .                +ys(njpip)-ys(njp)+ys(njpipkp)-ys(njpkp))
          zxi    =.25*(zs(njp)  -zs(n)  +zs(njpkp)  -zs(nkp)
     .                +zs(njpip)-zs(nip)+zs(njpipkp)-zs(nipkp))
          zeta   =.25*(zs(nkp)  -zs(n)  +zs(njpkp)  -zs(njp)
     .                +zs(nipkp)-zs(nip)+zs(njpipkp)-zs(njpip))
          zzeta  =.25*(zs(nip)  -zs(n)  +zs(nipkp)  -zs(nkp)
     .                +zs(njpip)-zs(njp)+zs(njpipkp)-zs(njpkp))
          ooj(n)  = 1./(  xxi*(yeta*zzeta-yzeta*zeta)
     .                  -xeta*(yxi *zzeta-yzeta*zxi )
     .                 +xzeta*(yxi *zeta -yeta *zxi )) 
          xix(n)  =(yeta*zzeta-yzeta*zeta)*ooj(n)
          xiy(n)  =(xzeta*zeta-xeta*zzeta)*ooj(n)
          xiz(n)  =(xeta*yzeta-xzeta*yeta)*ooj(n)
          etax(n) =(yzeta*zxi -yxi *zzeta)*ooj(n)
          etay(n) =(xxi*zzeta -xzeta* zxi)*ooj(n)
          etaz(n) =(xzeta*yxi -xxi *yzeta)*ooj(n)
          zetax(n)=(yxi *zeta -yeta*  zxi)*ooj(n)
          zetay(n)=(xeta*zxi  -xxi * zeta)*ooj(n)
          zetaz(n)=(xxi *yeta -xeta*  yxi)*ooj(n)
        else
          ooj(n)  = 0.
          xix(n)  = 0.
          xiy(n)  = 0.
          xiz(n)  = 0.
          etax(n) = 0.
          etay(n) = 0.
          etaz(n) = 0.
          zetax(n)= 0.
          zetay(n)= 0.
          zetaz(n)= 0.
        end if
c      end if
      enddo
  
      return
      end
c
      subroutine elrhs(b,u,islavept,nslave,nmaster,iseqr) 
      dimension islavept(nslave,nmaster,5)
      dimension b(3*nslave)
      dimension u(3*nslave)
c
c
      do n = 1,nslave
         ii4  = 1
         ni3  = 0 
         if(islavept(n,8,iseqr).eq.0) ii4 = 0
         if(ii4.ne.0) then
           iimax = islavept(n,11,iseqr)
           do ii2 = 2,iimax
             ni3 = islavept(n,12+ii2-2,iseqr)
             if(islavept(ni3,8,iseqr).eq.0) then
               ii4 = 0
               goto 1500
             end if
           enddo
         end if
1500     continue

       if(ii4.eq.0.and.ni3.eq.0) then
         b(3*(n-1)+1) =  u(3*(n-1)+1)
         b(3*(n-1)+2) =  u(3*(n-1)+2)
         b(3*(n-1)+3) =  u(3*(n-1)+3)
       else if(ii4.eq.0.and.ni3.ne.0) then
         b(3*(n-1)+1) =  u(3*(ni3-1)+1)
         b(3*(n-1)+2) =  u(3*(ni3-1)+2)
         b(3*(n-1)+3) =  u(3*(ni3-1)+3)
       else
         b(3*(n-1)+1) =  0.
         b(3*(n-1)+2) =  0.
         b(3*(n-1)+3) =  0.
       end if
      enddo

      return
      end
      subroutine deform_surf(nbl,idim,jdim,kdim,deltj,deltk,delti,
     .                       lw,lw2,icsi,icsf,jcsi,jcsf,kcsi,kcsf,
     .                       maxbl,mseq,time,dt,ita,nou,bou,nbuf,
     .                       ibufdim,myid,idefrm,nbci0,nbcidim,nbcj0,
     .                       nbcjdim,nbck0,nbckdim,ibcinfo,jbcinfo,
     .                       kbcinfo,maxseg,wk,u,nsurf,irst,iflag,
     .                       islavept,nslave,nsegdfrm,idfrmseg,iaesurf,
     .                       maxsegdg,nmaster,iseq)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Compute solid surface deformations for initializing the
c              mesh deformation scheme.
c
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
      dimension lw(65,maxbl),lw2(43,maxbl)
      dimension delti(jdim,kdim,3,2)
      dimension deltj(kdim,idim,3,2)
      dimension deltk(jdim,idim,3,2)
      dimension iaesurf(maxbl,maxsegdg)
      dimension ibcinfo(maxbl,maxseg,7,2)
      dimension icsf(maxbl,maxsegdg)
      dimension icsi(maxbl,maxsegdg)
      dimension idefrm(maxbl)
      dimension idfrmseg(maxbl,maxsegdg)
      dimension jbcinfo(maxbl,maxseg,7,2)
      dimension jcsf(maxbl,maxsegdg)
      dimension jcsi(maxbl,maxsegdg)
      dimension kbcinfo(maxbl,maxseg,7,2)
      dimension kcsf(maxbl,maxsegdg)
      dimension kcsi(maxbl,maxsegdg)
      dimension nbci0(maxbl)
      dimension nbcidim(maxbl)
      dimension nbcj0(maxbl)
      dimension nbcjdim(maxbl)
      dimension nbck0(maxbl)
      dimension nbckdim(maxbl)
      dimension nou(nbuf)
      dimension nsegdfrm(maxbl)
      dimension islavept(nslave,nmaster,5)
      dimension wk(9*nsurf)
      dimension u(3*nslave)
c
      allocatable :: dx(:)
      allocatable :: dy(:)
      allocatable :: dz(:)

      common /twod/ i2d
c
c     allocate memory
c
      memuse = 0
c
      allocate( dx(jdim*kdim*idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx',memuse,stats)
      allocate( dy(jdim*kdim*idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy',memuse,stats)
      allocate( dz(jdim*kdim*idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz',memuse,stats)
c
c     initialize deltas
c
      do i=1,jdim*kdim*idim
        dx(i) = 0.
        dy(i) = 0.
        dz(i) = 0.
      end do
c
c               preserve deltas on edges of solid surfaces
c
       if (idefrm(nbl).lt.999) then
          call bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,
     .                 jcsf,kcsi,kcsf,icsi,icsf,jdim,kdim,
     .                 idim,maxbl,maxsegdg,nsegdfrm)
c
        do n = 1,nslave
          nbl1= islavept(n,9,iseq)
          ll  = islavept(n,1,iseq)
          isrf= islavept(n,8,iseq)
          if(nbl1.eq.nbl.and.isrf.eq.0) then
            u(3*(n-1)+1) = dx(ll+1)
            u(3*(n-1)+2) = dy(ll+1)
            u(3*(n-1)+3) = dz(ll+1)
          end if
        enddo
      end if
c
c     release memory
c
      deallocate(dx)
      deallocate(dy)
      deallocate(dz)


      return
      end
