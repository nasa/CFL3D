      subroutine tau2x(jdim,kdim,idim,res,qr,js,ks,is,je,ke,ie,kode)
c***********************************************************************
c     Purpose:  Put the restricted residual from a finer embedded
c     mesh into a coarser mesh.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension qr(je-js+1,ke-ks+1,ie-is,5)
      dimension res(jdim,kdim,idim-1,5)
c
c      put restricted  r from finer embedded mesh into coarser mesh
c
c      jdim,kdim,idim  coarser mesh indices
c      js,ks,is        coarser mesh starting indices
c      je,ke,ie        coarser mesh ending indices
c
      if (kode.ge.2) then
         do 40 n=1,5
         kk = 0
         do 40 k=ks,ke-1
         kk = kk+1
         ii = 0
         do 40 i=is,ie-1
         ii = ii+1
         jj = 0
         do 40 j=js,je-1
         jj = jj+1
         qr(jj,kk,ii,n) =  qr(jj,kk,ii,n) - res(j,k,i,n)  
   40    continue
      end if
c
      do 60 n=1,5
      kk = 0
      do 60 k=ks,ke-1
      kk = kk+1
      ii = 0
      do 60 i=is,ie-1
      ii = ii+1
      jj = 0
      do 60 j=js,je-1
      jj = jj+1
      res(j,k,i,n) = res(j,k,i,n) + qr(jj,kk,ii,n)
   60 continue
      return
      end
