      subroutine collxtb(xtt,xttt,mdim,ndim,mm2,nn2,nbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Restrict xtb and atb arrays containing grid boundary 
c     velocity and acceleration to coarser meshes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xtt(mdim,ndim,3,2),xttt(mm2,nn2,3,2)
c
c      restrict xtt  to coarser mesh
c
c      mdim,ndim  finer mesh
c      mm2,nn2    coarser mesh
c
      ninc = 2
      if (ndim.eq.2) ninc = 1
      minc = 2
      if (mdim.eq.2) minc = 1
      do 10 ii=1,2
      do 10 ll=1,3
      nn   = 0
      do 10 n=1,ndim,ninc
      nn   = nn+1
      mm   = 0
      do 10 m=1,mdim,minc
      mm   = mm+1
      xttt(mm,nn,ll,ii) = xtt(m,n,ll,ii)
   10 continue
      return
      end
