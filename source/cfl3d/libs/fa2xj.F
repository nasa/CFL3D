      subroutine fa2xj(jf,kf,if,qjf,js,ks,is,je,ke,ie,
     .                 jbctyp,w,iw,nw,niw,neta,nou,bou,nbuf,ibufdim)
c***********************************************************************
c     Purpose:  Accumulate fluxes in J-direction for use on
c     twice coarser mesh to ensure conservation.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension qjf(kf,if-1,5,4)
      dimension w(1)
      dimension iw(1)
      dimension jbctyp(2)
      common /sklton/ isklton
c
c     js,ks,is - starting indices of coarser grid 
c                defining twice finer embedded grid
c     je,ke,ie - ending indices of coarser grid
c                defining twice finer embedded grid
c
c     jbctyp(1) - mesh indicator for j=0 boundary
c     jbctyp(2) - mesh indicator for j=jdim boundary
c
c     jf,kf,if - grid dimension of twice finer mesh
c     qjf - array containing fine grid fluxes
c
c     w - storage array for coarse grid flux
c     iw - storage array for coarse grid
c          starting and ending indices
c
c     nw - counter for w
c     niw - counter for iw (number of fine
c           grid flux accumulations)
c     neta - number of edges on which to accumulate
c
      jfm1 = jf-1
      kfm1 = kf-1
      ifm1 = if-1
      nsi  = ifm1/(ie-is)
c
      if (jbctyp(1).eq.21) then
c
c     left boundary
c
         mm        = 2
                        iw(niw+1) = js
         iw(niw+2) = ks
         iw(niw+3) = is
                        iw(niw+4) = js
         iw(niw+5) = ke
         iw(niw+6) = ie
         iw(niw+7) = nw+1
         niw       = niw+7
         neta      = neta+1
c
         if (nsi.eq.2) then
         do 100 l=1,5
         do 100 i=1,ifm1,2
         do 100 k=1,kfm1,2
         nw    = nw+1
         w(nw) = qjf(k,i,l,mm)
     .          +qjf(k+1,i,l,mm)
     .          +qjf(k,i+1,l,mm)
     .          +qjf(k+1,i+1,l,mm)
  100    continue
         else
         do 101 l=1,5
         do 101 i=1,ifm1
         do 101 k=1,kfm1,2
         nw    = nw+1
         w(nw) = qjf(k,i,l,mm)
     .          +qjf(k+1,i,l,mm)
  101    continue
         end if
c
         if (isklton.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)  '  summing flux at j=js  niw=',niw,
     .      ' nw=',nw
         end if
c
      end if
c
      if (jbctyp(2).eq.21) then
c
c     right boundary
c
         mm        = 4
                        iw(niw+1) = je
         iw(niw+2) = ks
         iw(niw+3) = is
                        iw(niw+4) = je
         iw(niw+5) = ke
         iw(niw+6) = ie
         iw(niw+7) = nw+1
         niw       = niw+7
         neta      = neta+1
c
         if (nsi.eq.2) then
         do 200 l=1,5
         do 200 i=1,ifm1,2
         do 200 k=1,kfm1,2
         nw    = nw+1
         w(nw) = qjf(k,i,l,mm)
     .          +qjf(k+1,i,l,mm)
     .          +qjf(k,i+1,l,mm)
     .          +qjf(k+1,i+1,l,mm)
  200    continue
         else
         do 201 l=1,5
         do 201 i=1,ifm1
         do 201 k=1,kfm1,2
         nw    = nw+1
         w(nw) = qjf(k,i,l,mm)
     .          +qjf(k+1,i,l,mm)
  201    continue
         end if
         if (isklton.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)  '  summing flux at j=je  niw=',niw,
     .      ' nw=',nw
         end if
c
      end if
      return
      end
