      subroutine hole(i,npl,jdim,kdim,idim,res,blank)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Zero out the rhs residuals for the blanked points.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension blank(jdim,kdim,idim)
      dimension res(jdim,kdim,idim-1,5) 
c
      nn = jdim*kdim
      do 130 ipl=1,npl
      ii = i+ipl-1
      do 100 l=1,5
cdir$ ivdep
      do 101 n=1,nn 
      res(n,1,ii,l) = res(n,1,ii,l)*blank(n,1,ii)
  101 continue
  100 continue
  130 continue
      return
      end 
