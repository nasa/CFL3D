      subroutine unld_qc(qc,qtemp,jc,kc,ic,is,ie,js,je,ks,ke,ldim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Install the qtemp array used for message passing into
c     the appropriate section of the qc array. Only the cell-center
c     locations are unloaded, not the fill in locations on the idim,
c     jdim and kdim planes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension qc(jc,kc,ic,ldim),qtemp(je-js+1,ke-ks+1,ie-is+1,ldim)
c
      do l =1,ldim
         jj = js - 1
         do j=1,je-js
            jj = jj + 1
            kk = ks - 1
            do k=1,ke-ks
               kk = kk + 1
               ii = is - 1
               do i=1,ie-is
                  ii = ii + 1
                  qc(jj,kk,ii,l) = qtemp(j,k,i,l)
               end do
            end do
         end do
      end do
      return
      end
