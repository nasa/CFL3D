      subroutine dabciz(k,npl,jdim,kdim,idim,a,b,c,blank)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Use the blank values to modify the coefficient matrices,
c     a,b,c , for the diagonal inversion in the I-direction.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(jdim-1,npl,idim),b(jdim-1,npl,idim), 
     .          c(jdim-1,npl,idim),blank(jdim,kdim,idim)
c
      jdim1 = jdim-1 
c
      do 10 i=1,idim 
      do 10 ipl=1,npl
      kk = k+ipl-1 
cdir$ ivdep
      do 11 j=1,jdim1
      a(j,ipl,i) =  a(j,ipl,i)*blank(j,kk,i)
      c(j,ipl,i) =  c(j,ipl,i)*blank(j,kk,i)
      b(j,ipl,i) = (b(j,ipl,i)*blank(j,kk,i))+(1.0-blank(j,kk,i))
   11 continue
   10 continue
      return 
      end
