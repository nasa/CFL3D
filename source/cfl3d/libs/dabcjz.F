      subroutine dabcjz(i,npl,jdim,kdim,idim,a,b,c,blank)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Use the blank values to modify the coefficient matrices,
c     a,b,c , for the diagonal inversion in the J-direction.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(kdim,npl,jdim),b(kdim,npl,jdim), 
     .          c(kdim,npl,jdim),blank(jdim,kdim,idim)
c
      do 10 j=1,jdim 
      do 10 ipl=1,npl
      ii = i+ipl-1 
cdir$ ivdep
      do 11 k=1,kdim
      a(k,ipl,j) =  a(k,ipl,j)*blank(j,k,ii)
      c(k,ipl,j) =  c(k,ipl,j)*blank(j,k,ii)
      b(k,ipl,j) = (b(k,ipl,j)*blank(j,k,ii))+(1.0-blank(j,k,ii))
   11 continue
   10 continue
      return 
      end
