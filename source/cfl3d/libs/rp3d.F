      subroutine rp3d(x,y,z,jdim,kdim,idim,igrid,ialph,igeom_img)
c***********************************************************************
c     Purpose:  Read grids in plot3d format.
c     ialph     - flag to interpret angle of attack
c               = 0 alpha measured in x-z plane (cfl3d standard)
c               > 0 alpha measured in x-y plane (tlns3d standard)
c     igeom_img - flag for complex grid file
c               = 0 input grid is real valued
c               > 0 input grid is complex valued
c***********************************************************************
c
#ifdef CMPLX
      implicit complex(a-h,o-z)
c
      real xr(jdim,kdim,idim),yr(jdim,kdim,idim),
     .     zr(jdim,kdim,idim)
c
      pointer (ip_xr,xr),(ip_yr,yr),(ip_zr,zr)
#endif
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c
#ifdef CMPLX
c
c     allocate real arrays in which to read grid data
c     (for non-geometric derivatives)
c
      icall  = 0
      memuse = 0
      call umalloc(ip_xr,jdim*kdim*idim,0,icall,memuse)
      call umalloc(ip_yr,jdim*kdim*idim,0,icall,memuse)
      call umalloc(ip_zr,jdim*kdim*idim,0,icall,memuse)
#endif
c
      if (ialph.eq.0) then
#ifdef CMPLX
         if (igeom_img .eq. 0) then
c        
c           for non-geometric derivatives, input grid is real
c
            read(1) (((xr(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .              (((yr(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .              (((zr(j,k,i),i=1,idim),j=1,jdim),k=1,kdim)
            do i=1,idim
               do j=1,jdim
                  do k=1,kdim
                     x(j,k,i) = xr(j,k,i)
                     y(j,k,i) = yr(j,k,i)
                     z(j,k,i) = zr(j,k,i)
                  end do
               end do
            end do
         else
c        
c           for geometric derivatives, input grid is complex
c
            read(1) (((x(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .              (((y(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .              (((z(j,k,i),i=1,idim),j=1,jdim),k=1,kdim)

         end if
#else
         read(1) (((x(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .           (((y(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .           (((z(j,k,i),i=1,idim),j=1,jdim),k=1,kdim)
#endif
      else
#ifdef CMPLX
         if (igeom_img .eq. 0) then
c
c           for non-geometric derivatives, input grid is real
c
            read(1) (((xr(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .              (((zr(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .              (((yr(j,k,i),i=1,idim),j=1,jdim),k=1,kdim)
            do i=1,idim
               do j=1,jdim
                  do k=1,kdim
                     x(j,k,i) = xr(j,k,i)
                     y(j,k,i) = yr(j,k,i)
                     z(j,k,i) = zr(j,k,i)
                  end do
               end do
            end do
         else
c
c           for geometric derivatives, input grid is complex
c
            read(1) (((x(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .              (((z(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .              (((y(j,k,i),i=1,idim),j=1,jdim),k=1,kdim)
         end if
#else
         read(1) (((x(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .           (((z(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
     .           (((y(j,k,i),i=1,idim),j=1,jdim),k=1,kdim)
#endif
         do i=1,idim
            do j=1,jdim
               do k=1,kdim
                  y(j,k,i) = -y(j,k,i)
               end do
            end do
         end do
      end if
#   ifdef CMPLX
c
c     deallocate real arrays in which grid data was read
c
      call ufree(ip_xr)
      call ufree(ip_yr)
      call ufree(ip_zr)
#   endif
c
      return
      end	
