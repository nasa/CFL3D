      subroutine gfluxl(i,npl,xkap,idf,jdim,kdim,idim,res,q,qj0,sj,
     .                  dgp,dgm,t,nvtq)
c***********************************************************************
c     Purpose:  Compute the left-hand flux contributions due to the
c     inviscid terms for the J-direction.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension sj(jdim*kdim,idim-1,5)
      dimension q(jdim,kdim,idim,5),qj0(kdim,idim-1,5,4)
      dimension res(jdim,kdim,idim-1,5),t(nvtq,39)
      dimension dgp(jdim*npl*(kdim-1),5,5),dgm(jdim*npl*(kdim-1),5,5)
c
      common /fvfds/ rkap0(3),ifds(3)
c
c     g-flux-p-m jacobians
c
      kdim1 = kdim-1
c
      n     = jdim*kdim1
      n0    = n-1
c
      do 8008 ipl=1,npl
      ii = i+ipl-1
      do 8008 l=1,5
      do 8008 k=1,kdim
      q(jdim,k,ii,l) = qj0(k,ii,l,3)
 8008 continue
c
      do 8000 ipl=1,npl
      jkv1 = 1 + (ipl-1)*jdim*kdim1
      ii   = i+ipl-1
      do 9000 l=1,4
cdir$ ivdep
      do 1000 izz=1,n
      t(izz+jkv1-1,35+l) = sj(izz,ii,l)
 1000 continue
 9000 continue
cdir$ ivdep
      do 1001 izz=1,n
      t(izz+jkv1-1,20) = sj(izz,ii,5)
 1001 continue
 8000 continue
c
      do 40 ipl=1,npl
      ii  = i+ipl-1
      jkv = (ipl-1)*jdim*kdim1
      do 40 l=1,5
cdir$ ivdep
      do 1002 izz=1,n0
      t(izz+1+jkv,20+l) = q(izz,1,ii,l)
 1002 continue
cdir$ ivdep
      do 1003 izz=1,n
      t(izz+jkv,25+l) = q(izz,1,ii,l)
 1003 continue
      do 40 k=1,kdim1
      jk1         = jkv + (k-1)*jdim  +  1
      t(jk1,20+l) = qj0(k,ii,l,1)
      jk2         = jkv + k*jdim
      t(jk2,25+l) = q(jdim,k,ii,l)
   40 continue
c
      n     = npl*n
      jkpro = jdim*kdim1*npl
      if (ifds(2).eq.0) then
         call dfluxpm(t(1,1), t(1,2), t(1,36),t(1,37),t(1,38),t(1,39),
     .                t(1,20),t(1,21),dgp,n,jkpro,nvtq,+1)
      else
         call dfhat(t(1,36),t(1,37),t(1,38),t(1,39),t(1,20),t(1,21),dgp,
     .              n,nvtq,+1)
      end if
c
      if (ifds(2).eq.0) then
         call dfluxpm(t(1,1), t(1,2), t(1,36),t(1,37),t(1,38),t(1,39),
     .                t(1,20),t(1,26),dgm,n,jkpro,nvtq,-1)
      else
         call dfhat(t(1,36),t(1,37),t(1,38),t(1,39),t(1,20),t(1,26),dgm,
     .              n,nvtq,-1)
      end if
      return
      end
