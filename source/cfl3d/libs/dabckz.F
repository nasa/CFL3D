      subroutine dabckz(i,npl,jdim,kdim,idim,a,b,c,blank)
c***********************************************************************
c     Purpose:  Use the blank values to modify the coefficient matrices,
c     a,b,c , for the diagonal inversion in the K-direction.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(jdim-1,npl,kdim),b(jdim-1,npl,kdim), 
     .          c(jdim-1,npl,kdim),blank(jdim,kdim,idim)
c
      jdim1 = jdim-1 
c
      do 10 k=1,kdim 
      do 10 ipl=1,npl
      ii = i+ipl-1 
cdir$ ivdep
      do 11 j=1,jdim1
      a(j,ipl,k) = a(j,ipl,k)*blank(j,k,ii)
      c(j,ipl,k) = c(j,ipl,k)*blank(j,k,ii)
      b(j,ipl,k) = (b(j,ipl,k)*blank(j,k,ii))+(1.0-blank(j,k,ii))
   11 continue
   10 continue
      return 
      end
