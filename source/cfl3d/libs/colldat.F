      subroutine colldat(bcdata,mdim,ndim,bcdatac,mdimc,ndimc)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Restrict auxilary boundary condition data arrays to 
c     coarser meshes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension bcdata(mdim,ndim,2,7),bcdatac(mdimc,ndimc,2,7)
c
      do 10 l=1,7
      mm   = 0
      minc = 2
      mp   = 1
      if (mdim.eq.1) then
         minc = 1
         mp = 0
      end if
      do 10 m=1,mdim,minc
      mm   = mm+1
      mm   = min(mm,mdimc)
      nn   = 0
      ninc = 2
      np   = 1
      if (ndim.eq.1) then
         ninc = 1
         np = 0
      end if
      do 10 n=1,ndim,ninc
      nn   = nn+1
      nn   = min(nn,ndimc)
      nnn  = n+np
      nnn  = min(nnn,ndim)
      mmm  = m+mp
      mmm  = min(mmm,mdim)
      bcdatac(mm,nn,1,l) =
     .   .25*(bcdata(m,n,1,l) + bcdata(mmm,n,1,l) +
     .        bcdata(m,nnn,1,l) + bcdata(mmm,nnn,1,l))
      bcdatac(mm,nn,2,l) =
     .   .25*(bcdata(m,n,2,l) + bcdata(mmm,n,2,l) +
     .        bcdata(m,nnn,2,l) + bcdata(mmm,nnn,2,l))
   10 continue
      return 
      end
