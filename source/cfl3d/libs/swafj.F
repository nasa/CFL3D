      subroutine swafj(i,npl,jdim,kdim,idim,aj,bj,cj,f,nvt,res,iperd,
     .                 gj,hj)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Solve the block 5x5 tridiagonal equations for the 
c     3-factor spatially-split algorithm in the J-direction.
c     Modified for Weiss-Smith preconditioning by J.R. Edwards, NCSU
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension res(jdim,kdim,idim-1,5)
      dimension aj(npl*(kdim-1),jdim,5,5),bj(npl*(kdim-1),jdim,5,5),
     .          cj(npl*(kdim-1),jdim,5,5),gj(npl*(kdim-1),jdim,5,5),
     .          hj(npl*(kdim-1),jdim,5,5), f(npl*(kdim-1),jdim,5)
c
c     j-implicit k-sweep line inversions af
c
      jdim1 = jdim-1
      kdim1 = kdim-1
c
c     load rhs (-residual) into f
c
      kv  = npl*kdim1
      do 1005 ipl=1,npl
      ii  = i+ipl-1
      kv0 = (ipl-1)*kdim1
      do 1005 l=1,5
      do 1005 k=1,kdim1
      k0  = k + kv0
c
      jj  = 1-kv
cdir$ ivdep
      do 7665 jjj=1,jdim1
      jj  = jj+kv
 7665 f(k0+jj-1,1,l) = res(jjj,k,ii,l)
c      call q8vscatp(jdim1,res(1,k,ii,l),kv,jdim1,jdim1,f(k0,1,l))
 1005 continue
c
      n = kv*jdim1
      do 1010 l=1,5
cdir$ ivdep
      do 1000 izz=1,n
      f(izz,1,l) = -f(izz,1,l)
 1000 continue
 1010 continue
c
c     solve matrix equation
c
      il  = 1
      iu  = jdim1
      n   = kv
c
      id1 = npl*(kdim-1)
      if (iperd.eq.1) then
         call bsubp(id1,jdim,aj,bj,cj,f,1,n,il,iu,gj,hj)
      else
         call bsub(id1,jdim,aj,bj,cj,f,1,n,il,iu)
      end if
c
c     update delta q
c
      do 1300 ipl=1,npl
      ii  = i+ipl-1
      kv0 = (ipl-1)*kdim1
      do 1300 l=1,5
      do 1300 k=1,kdim1
      k0  = k + kv0
c
      jj  = 1-kv
cdir$ ivdep
      do 7891 jjj=1,jdim1
      jj  = jj+kv
 7891 res(jjj,k,ii,l) = f(k0+jj-1,1,l)
c      call q8vgathp(jdim1,f(k0,1,l),kv,jdim1,jdim1,res(1,k,ii,l))
 1300 continue
 2000 continue
      return
      end
