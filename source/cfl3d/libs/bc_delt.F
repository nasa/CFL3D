      subroutine bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,jcsf,
     .                   kcsi,kcsf,icsi,icsf,jdim,kdim,idim,maxbl,
     .                   maxsegdg,nsegdfrm)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Install boundary deltas from forced oscillation or
c     aeroelastic motion into the 3d displacement array
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension dx(jdim,kdim,idim),dy(jdim,kdim,idim),dz(jdim,kdim,idim)
      dimension deltj(kdim,idim,3,2),deltk(jdim,idim,3,2),
     .          delti(jdim,kdim,3,2)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl)
c
      do iseg=1,nsegdfrm(nbl)
c
         if (icsi(nbl,iseg).eq.icsf(nbl,iseg)) then
            i = icsi(nbl,iseg)
            if (icsi(nbl,iseg).eq.1) then
               ii = 1
            else
               ii = 2
            end if
            do j=jcsi(nbl,iseg),jcsf(nbl,iseg)
               do k=kcsi(nbl,iseg),kcsf(nbl,iseg)
                  dx(j,k,i)    = delti(j,k,1,ii)
                  dy(j,k,i)    = delti(j,k,2,ii)
                  dz(j,k,i)    = delti(j,k,3,ii)
               end do
            end do
         end if
c
         if (jcsi(nbl,iseg).eq.jcsf(nbl,iseg)) then
            j = jcsi(nbl,iseg)
            if (jcsi(nbl,iseg).eq.1) then
               jj = 1
            else
               jj = 2
            end if
            do k=kcsi(nbl,iseg),kcsf(nbl,iseg)
               do i=icsi(nbl,iseg),icsf(nbl,iseg)
                  dx(j,k,i)    = deltj(k,i,1,jj)
                  dy(j,k,i)    = deltj(k,i,2,jj)
                  dz(j,k,i)    = deltj(k,i,3,jj)
               end do
            end do
         end if
         if (kcsi(nbl,iseg).eq.kcsf(nbl,iseg)) then
            k = kcsi(nbl,iseg)
            if (kcsi(nbl,iseg).eq.1) then
               kk = 1
            else
               kk = 2
            end if
            do j=jcsi(nbl,iseg),jcsf(nbl,iseg)
               do i=icsi(nbl,iseg),icsf(nbl,iseg)
                  dx(j,k,i)    = deltk(j,i,1,kk)
                  dy(j,k,i)    = deltk(j,i,2,kk)
                  dz(j,k,i)    = deltk(j,i,3,kk)
               end do
            end do
         end if
c
      end do
c
      return
      end
