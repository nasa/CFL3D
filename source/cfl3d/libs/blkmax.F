      subroutine blkmax(jdim,kdim,idim,res,resmax,jm,km,im)
c*************************************************************************
c     Purpose:  Determine the location of the maximum residual.
c*************************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension res(jdim,kdim,idim-1,5)
      resmax = res(1,1,1,1)
      jm     = 1
      km     = 1
      im     = 1
      do 10 j=1,jdim
      do 10 k=1,kdim
      do 10 i=1,idim-1
      if (real(res(j,k,i,1)).gt.real(resmax)) then
         resmax = res(j,k,i,1)
         jm = j
         km = k
         im = i
      end if
   10 continue
      return
      end
