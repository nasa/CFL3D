      subroutine avghole(q,blank,jdim,kdim,idim,nbl,ldim,int_updt)
c***********************************************************************
c      Purpose: Replace the solution at any point with blank = 0 with
c               the average of neighboring points. This affects hole,
c               points, orphan points AND fringe points; however, 
c               fringe points are subsequently updated properly by
c               a call to subroutine xupdt. This allows any orhpan 
c               points to be updated, even though no interpolation
c               stencils are available for them. 
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,ldim),blank(jdim,kdim,idim)
c
c
      do l=1,ldim
         do i=1,idim-1
            do k=1,kdim-1
               do j=1,jdim-1
                  if (blank(j,k,i) .eq. 0.) then
                     jp = min(j+1,jdim-1)
                     jm = max(j-1,1)
                     kp = min(k+1,kdim-1)
                     km = max(k-1,1)
                     ip = min(i+1,idim-1)
                     im = max(i-1,1)
                     q(j,k,i,l) = (q(jp, k,ip,l) + q(jm, k,ip,l)
     .                          +  q( j,kp,ip,l) + q( j,km,ip,l)
     .                          +  q(jp, k,im,l) + q(jm, k,im,l)
     .                          +  q( j,kp,im,l) + q( j,km,im,l))*0.125
                  end if
               end do
            end do
         end do
      end do
      return
      end
