      subroutine fill(jdim,kdim,idim,q,ll)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Fills the edges of the q array for safety using
c     multi-plane vectorization technique.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,ll)
c
c     fill edges of q array for safety
c
      idim1 = idim-1
      kdim1 = kdim-1
      do 3227 l=1,ll
      do 3226 i=1,idim1
cdir$ ivdep
      do 3229 k=1,kdim1
      q(jdim,k,i,l) = q(jdim-1,k,i,l)
 3229 continue
      do 1000 izz=1,jdim
      q(izz,kdim,i,l) = q(izz,kdim-1,i,l)
 1000 continue
 3226 continue
      n = jdim*kdim
cdir$ ivdep
      do 1001 izz=1,n
      q(izz,1,idim,l) = q(izz,1,idim-1,l)
 1001 continue
 3227 continue
      return
      end
