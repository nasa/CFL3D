      subroutine ld_qc(qc,qtemp,jc,kc,ic,is,ie,js,je,ks,ke,ldim)
c***********************************************************************
c     Purpose:  Install the embedded subset of the qc array into qtemp
c     for massage passing.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension qc(jc,kc,ic,ldim),qtemp(je-js+1,ke-ks+1,ie-is+1,ldim)
c
      do l =1,ldim
         jj = 0
         do j=js,je
            jj = jj + 1
            kk = 0
            do k=ks,ke
               kk = kk + 1
               ii = 0
               do i=is,ie
                  ii = ii + 1
                  qtemp(jj,kk,ii,l) = qc(j,k,i,l)
               end do
            end do
         end do
      end do
      return
      end
