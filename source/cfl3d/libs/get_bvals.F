	subroutine get_bvals(jdim1,kdim1,idim1,jmax1,kmax1,
     .                       ldim,q,qbou,mtype,i1)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  To pick out values at the boundary from global values
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
        dimension q(jdim1, kdim1, idim1, ldim)
        dimension qbou(jmax1, kmax1, ldim)
c
        if (mtype. eq.1) then
c
           do j = 1,jmax1
              do k = 1,kmax1
                 do l = 1,ldim
                    qbou(j,k,l) = q(j,k,i1,l)
                 end do
              end do
           end do
c
        else if (mtype. eq. 2) then
c
           do j = 1,jmax1
              do k = 1,kmax1
                 do l = 1,ldim
                    qbou(j,k,l) = q(i1,j,k,l)
                 end do
              end do
           end do             
c
        else if (mtype. eq. 3) then
c
           do j = 1,jmax1
              do k = 1,kmax1
                 do l = 1,ldim
                    qbou(j,k,l) = q(j,i1,k,l)
                 end do
              end do
           end do
c
        end if
c
        return
        end
