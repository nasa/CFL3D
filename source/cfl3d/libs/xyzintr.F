      subroutine xyzintr(jdimf,kdimf,idimf,jdimc,kdimc,idimc,xf,yf,
     .                   zf,xc,yc,zc)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Interpolate coarse grid x,y,z to fine grid for mesh
c              sequencing in deforming mesh cases
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xf(jdimf,kdimf,idimf),yf(jdimf,kdimf,idimf),
     .          zf(jdimf,kdimf,idimf)
      dimension xc(jdimc,kdimc,idimc),yc(jdimc,kdimc,idimc),
     .          zc(jdimc,kdimc,idimc)
c
      common /twod/ i2d
c
      if (i2d.eq.0) then
c
         do jf=1,jdimf,2
            jc = (jf+1)/2
            do kf=1,kdimf,2
               kc = (kf+1)/2
               do if=1,idimf,2
                  ic = (if+1)/2
                  xf(jf,kf,if) = xc(jc,kc,ic)
                  yf(jf,kf,if) = yc(jc,kc,ic)
                  zf(jf,kf,if) = zc(jc,kc,ic)
               end do
            end do
         end do
         do jf=2,jdimf,2
            do kf=1,kdimf,2
               do if=1,idimf,2
                  xf(jf,kf,if) = 0.5*(xf(jf-1,kf,if)+xf(jf+1,kf,if))
                  yf(jf,kf,if) = 0.5*(yf(jf-1,kf,if)+yf(jf+1,kf,if))
                  zf(jf,kf,if) = 0.5*(zf(jf-1,kf,if)+zf(jf+1,kf,if))
               end do
            end do
         end do
         do jf=1,jdimf
            do kf=2,kdimf,2
               do if=1,idimf,2
                  xf(jf,kf,if) = 0.5*(xf(jf,kf-1,if)+xf(jf,kf+1,if))
                  yf(jf,kf,if) = 0.5*(yf(jf,kf-1,if)+yf(jf,kf+1,if))
                  zf(jf,kf,if) = 0.5*(zf(jf,kf-1,if)+zf(jf,kf+1,if))
               end do
            end do
         end do
         do jf=1,jdimf
            do kf=1,kdimf
               do if=2,idimf,2
                  xf(jf,kf,if) = 0.5*(xf(jf,kf,if-1)+xf(jf,kf,if+1))
                  yf(jf,kf,if) = 0.5*(yf(jf,kf,if-1)+yf(jf,kf,if+1))
                  zf(jf,kf,if) = 0.5*(zf(jf,kf,if-1)+zf(jf,kf,if+1))
               end do
            end do
         end do
c
      else
c
         do jf=1,jdimf,2
            jc = (jf+1)/2
            do kf=1,kdimf,2
               kc = (kf+1)/2
               do i=1,2
                  xf(jf,kf,i) = xc(jc,kc,i)
                  yf(jf,kf,i) = yc(jc,kc,i)
                  zf(jf,kf,i) = zc(jc,kc,i)
               end do
            end do
         end do
         do jf=2,jdimf,2
            do kf=1,kdimf,2
               do i=1,2
                  xf(jf,kf,i) = 0.5*(xf(jf-1,kf,i)+xf(jf+1,kf,i))
                  yf(jf,kf,i) = 0.5*(yf(jf-1,kf,i)+yf(jf+1,kf,i))
                  zf(jf,kf,i) = 0.5*(zf(jf-1,kf,i)+zf(jf+1,kf,i))
               end do
            end do
         end do
         do jf=1,jdimf
            do kf=2,kdimf,2
               do i=1,2
                  xf(jf,kf,i) = 0.5*(xf(jf,kf-1,i)+xf(jf,kf+1,i))
                  yf(jf,kf,i) = 0.5*(yf(jf,kf-1,i)+yf(jf,kf+1,i))
                  zf(jf,kf,i) = 0.5*(zf(jf,kf-1,i)+zf(jf,kf+1,i))
               end do
            end do
         end do
c
      end if
c
      return
      end
