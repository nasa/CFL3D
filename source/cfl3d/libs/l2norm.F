      subroutine l2norm(nbl,ntime,rmsl,irdq,jdim,kdim,idim,res,vol)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute the L2-norm of the residuals or the change in
c     primative variables from one time to the next.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension res(jdim,kdim,idim-1,5)
      dimension vol(jdim,kdim,idim-1)
c
c      l2 norm of residual  / delta q
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
      nplq  = min(idim1,999000/(jdim*kdim))
      npl   = nplq
      rmsl  = 0.e0
      if (irdq.eq.1) go to 1000
c
      do 500 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      do 500 l=1,5
      do 8900 ipl=1,npl
      ii = i+ipl-1
      do 1001 j=1,jdim
      res(j,kdim,ii,l) = 0.
 1001 continue
c
cdir$ ivdep
      do 5671 k=1,kdim1
 5671 res(jdim,k,ii,l) = 0.
 8900 continue
  500 continue
c
      meq = 1
      do 510 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      n = npl*jdim*kdim - jdim -1
      do 510 l=1,meq
      rmsl = rmsl+q8sdot(n,res(1,1,i,l),n,res(1,1,i,l))
  510 continue
c
      return
 1000 continue
      n = jdim*kdim - jdim -1
      do 2000 i=1,idim1
c
      rmsl = 0.e0
      do 610 l=1,5
c
      do 5699 k=1,kdim1
 5699 res(jdim,k,i,l) = 0. 
      rmsl = rmsl+q8sdot(n,res(1,1,i,l),n,res(1,1,i,l))
  610 continue
c
 2000 continue
      return
      end


      subroutine l2normAll(nbl,ntime,rmsl,irdq,jdim,kdim,idim,res,vol)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute the L2-norm of the residuals or the change in
c     primative variables from one time to the next.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension res(jdim,kdim,idim-1,5)
      dimension vol(jdim,kdim,idim-1)
      dimension rmsl(5)
c
c      l2 norm of residual  / delta q
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
      nplq  = min(idim1,999000/(jdim*kdim))
      npl   = nplq
      rmsl  = 0.e0
      if (irdq.eq.1) go to 1000
c
      do 500 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      do 500 l=1,5
      do 8900 ipl=1,npl
      ii = i+ipl-1
      do 1001 j=1,jdim
      res(j,kdim,ii,l) = 0.
 1001 continue
c
cdir$ ivdep
      do 5671 k=1,kdim1
 5671 res(jdim,k,ii,l) = 0.
 8900 continue
  500 continue
c
      meq = 5
      do 510 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      n = npl*jdim*kdim - jdim -1
      do 510 l=1,meq
      rmsl(l) = rmsl(l)+q8sdot(n,res(1,1,i,l),n,res(1,1,i,l))
  510 continue
c
      return
 1000 continue
      n = jdim*kdim - jdim -1
      do 2000 i=1,idim1
c
      rmsl = 0.e0
      do 610 l=1,5
c
      do 5699 k=1,kdim1
 5699 res(jdim,k,i,l) = 0. 
      rmsl(l) = rmsl(l)+q8sdot(n,res(1,1,i,l),n,res(1,1,i,l))
  610 continue
c
 2000 continue
      return
      end
