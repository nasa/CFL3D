      subroutine setblk(blank,jdim,kdim,idim,nbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Initialize the blank array.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension blank(jdim,kdim,idim)
      do 10 i=1,idim
      do 10 j=1,jdim
      do 10 k=1,kdim
      blank(j,k,i) = 1.0
   10 continue
      return
      end
