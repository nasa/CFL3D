      subroutine bc2102(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,mdim,ndim,bcdata,
     .                  filname,iuns,nbl,nou,bou,nbuf,ibufdim,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set pressure ratio as a function of time; extrapolate 
c     other quantities
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
#   ifdef CMPLX
      complex lref
#   else
      real lref
#   endif
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension bcdata(mdim,ndim,2,7)
      dimension tursav(jdim,kdim,idim,2),tj0(kdim,idim-1,2,4),
     .          tk0(jdim,idim-1,2,4),ti0(jdim,kdim,2,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
c
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10,tur20
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /unst/ time,cfltau,ntstep,ita,iunst
c
      twopi = 8.0*atan(1.0)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c     this bc makes use of only one plane of data    
c
      ip    = 1
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=2102 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary         set p(time), extrapolate others         bctype 2102
c******************************************************************************
c
      if (nface.eq.3) then
c
      if (isklton.eq.1) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 1 l=1,4
         do 1 ipp=1,2
         do 1 i=ista,iend1
         ii = i-ista+1
         do 1 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   1     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2102:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c 
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) 
     .            real(twopi)*real(rfreqp)/real(lref)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2001)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      do 100 l=1,4
      do 100 i=ista,iend1
      do 100 k=ksta,kend1
      qj0(k,i,l,1) = q(1,k,i,l)
      qj0(k,i,l,2) = qj0(k,i,l,1)
      bcj(k,i,1) = 0.0
  100 continue
      do 102 i=ista,iend1
      ii = i-ista+1
      do 102 k=ksta,kend1
      kk = k-ksta+1
      pratio = bcdata(kk,ii,ip,1)
      deltap = bcdata(kk,ii,ip,2)
      rfreqp = bcdata(kk,ii,ip,3)
      lref   = bcdata(kk,ii,ip,4)
      pset   = (pratio + deltap*sin(twopi*rfreqp/lref*time))*p0
      qj0(k,i,5,1) = pset
      qj0(k,i,5,2) = pset
  102 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,2) = vist3d(1,k,i)
  191   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 101 i=ista,iend1
        do 101 k=ksta,kend1
          tj0(k,i,1,1) = tursav(1,k,i,1)
          tj0(k,i,2,1) = tursav(1,k,i,2)
          tj0(k,i,1,2) = tursav(1,k,i,1)
          tj0(k,i,2,2) = tursav(1,k,i,2)
  101   continue
      end if
      end if
      end if
c
c******************************************************************************
c      j=jdim boundary         set p(time), extrapolate others      bctype 2102
c******************************************************************************
c 
      if (nface.eq.4) then
c
      if (isklton.eq.1) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 2 l=1,4
         do 2 ipp=1,2
         do 2 i=ista,iend1
         ii = i-ista+1
         do 2 k=ksta,kend1
         kk = k-ksta+1
         if (real(bcdata(kk,ii,ipp,l)) .lt. -1.e10) iflg = 1
   2     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2102:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010) 
     .            real(twopi)*real(rfreqp)/real(lref)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2002)ista,iend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      do 200 l=1,4
      do 200 i=ista,iend1
      do 200 k=ksta,kend1
      qj0(k,i,l,3) = q(jdim1,k,i,l)
      qj0(k,i,l,4) = qj0(k,i,l,3)
      bcj(k,i,2) = 0.0
  200 continue
      do 202 i=ista,iend1
      ii = i-ista+1
      do 202 k=ksta,kend1
      kk = k-ksta+1
      pratio = bcdata(kk,ii,ip,1)
      deltap = bcdata(kk,ii,ip,2)
      rfreqp = bcdata(kk,ii,ip,3)
      lref   = bcdata(kk,ii,ip,4)
      pset   = (pratio + deltap*sin(twopi*rfreqp/lref*time))*p0
      qj0(k,i,5,3) = pset
      qj0(k,i,5,4) = pset
  202 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
          vj0(k,i,1,3) = vist3d(jdim1,k,i)
          vj0(k,i,1,4) = vist3d(jdim1,k,i)
  291   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 201 i=ista,iend1
        do 201 k=ksta,kend1
          tj0(k,i,1,3) = tursav(jdim1,k,i,1)
          tj0(k,i,2,3) = tursav(jdim1,k,i,2)
          tj0(k,i,1,4) = tursav(jdim1,k,i,1)
          tj0(k,i,2,4) = tursav(jdim1,k,i,2)
  201   continue
      end if
      end if
      end if
c
c******************************************************************************
c      k=1 boundary         set p(time), extrapolate others         bctype 2102
c******************************************************************************
c
      if (nface.eq.5) then
c
      if (isklton.eq.1) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 3 l=1,4
         do 3 ipp=1,2
         do 3 i=ista,iend1
         ii = i-ista+1
         do 3 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   3     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2102:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010)
     .            real(twopi)*real(rfreqp)/real(lref)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2003)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      do 300 l=1,4
      do 300 i=ista,iend1
      do 300 j=jsta,jend1
      qk0(j,i,l,1) = q(j,1,i,l)
      qk0(j,i,l,2) = qk0(j,i,l,1)
      bck(j,i,1) = 0.0
  300 continue
      do 302 i=ista,iend1
      ii = i-ista+1
      do 302 j=jsta,jend1
      jj = j-jsta+1
      pratio = bcdata(jj,ii,ip,1)
      deltap = bcdata(jj,ii,ip,2)
      rfreqp = bcdata(jj,ii,ip,3)
      lref   = bcdata(jj,ii,ip,4)
      pset   = (pratio + deltap*sin(twopi*rfreqp/lref*time))*p0
      qk0(j,i,5,1) = pset
      qk0(j,i,5,2) = pset
  302 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = vist3d(j,1,i)
  391   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 301 i=ista,iend1
        do 301 j=jsta,jend1
          tk0(j,i,1,1) = tursav(j,1,i,1)
          tk0(j,i,2,1) = tursav(j,1,i,2)
          tk0(j,i,1,2) = tursav(j,1,i,1)
          tk0(j,i,2,2) = tursav(j,1,i,2)
  301   continue
      end if
      end if
      end if
c
c******************************************************************************
c      k=kdim boundary         set p(time), extrapolate others      bctype 2102
c******************************************************************************
c
      if (nface.eq.6) then
c
      if (isklton.eq.1) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 4 l=1,4
         do 4 ipp=1,2
         do 4 i=ista,iend1
         ii = i-ista+1
         do 4 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,ii,ipp,l)) .lt. -1.e10) iflg = 1
   4     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2102:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010)
     .            real(twopi)*real(rfreqp)/real(lref)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2004)ista,iend,jsta,jend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      do 400 l=1,4
      do 400 i=ista,iend1
      do 400 j=jsta,jend1
      qk0(j,i,l,3) = q(j,kdim1,i,l)
      qk0(j,i,l,4) = qk0(j,i,l,3)
      bck(j,i,2) = 0.0
  400 continue
      do 402 i=ista,iend1
      ii = i-ista+1
      do 402 j=jsta,jend1
      jj = j-jsta+1
      pratio = bcdata(jj,ii,ip,1)
      deltap = bcdata(jj,ii,ip,2)
      rfreqp = bcdata(jj,ii,ip,3)
      lref   = bcdata(jj,ii,ip,4)
      pset   = (pratio + deltap*sin(twopi*rfreqp/lref*time))*p0
      qk0(j,i,5,3) = pset
      qk0(j,i,5,4) = pset
  402 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
          vk0(j,i,1,3) = vist3d(j,kdim1,i)
          vk0(j,i,1,4) = vist3d(j,kdim1,i)
  491   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 401 i=ista,iend1
        do 401 j=jsta,jend1
          tk0(j,i,1,3) = tursav(j,kdim1,i,1)
          tk0(j,i,2,3) = tursav(j,kdim1,i,2)
          tk0(j,i,1,4) = tursav(j,kdim1,i,1)
          tk0(j,i,2,4) = tursav(j,kdim1,i,2)
  401   continue
      end if
      end if
      end if
c
c******************************************************************************
c      i=1 boundary         set p(time), extrapolate others         bctype 2102
c******************************************************************************
c
      if (nface.eq.1) then
c
      if (isklton.eq.1) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 5 l=1,4
         do 5 ipp=1,2
         do 5 k=ksta,kend1
         kk = k-ksta+1
         do 5 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   5     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2102:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010)
     .            real(twopi)*real(rfreqp)/real(lref)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2005)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      do 500 l=1,4
      do 500 k=ksta,kend1
      do 500 j=jsta,jend1
      qi0(j,k,l,1) = q(j,k,1,l)
      qi0(j,k,l,2) = qi0(j,k,l,1)
      bci(j,k,1) = 0.0
  500 continue
      do 502 k=ksta,kend1
      kk = k-ksta+1
      do 502 j=jsta,jend1
      jj = j-jsta+1
      pratio = bcdata(jj,kk,ip,1)
      deltap = bcdata(jj,kk,ip,2)
      rfreqp = bcdata(jj,kk,ip,3)
      lref   = bcdata(jj,kk,ip,4)
      pset   = (pratio + deltap*sin(twopi*rfreqp/lref*time))*p0
      qi0(j,k,5,1) = pset
      qi0(j,k,5,2) = pset
  502 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = vist3d(j,k,1)
  591   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 501 k=ksta,kend1
        do 501 j=jsta,jend1
          ti0(j,k,1,1) = tursav(j,k,1,1)
          ti0(j,k,2,1) = tursav(j,k,1,2)
          ti0(j,k,1,2) = tursav(j,k,1,1)
          ti0(j,k,2,2) = tursav(j,k,1,2)
  501   continue
      end if
      end if
      end if
c
c******************************************************************************
c      i=idim boundary         set p(time), extrapolate others      bctype 2102
c******************************************************************************
c
      if (nface.eq.2) then
c
      if (isklton.eq.1) then
c
c        check that the required auxiliary data has been set
         iflg = 0
         do 6 l=1,4
         do 6 ipp=1,2
         do 6 k=ksta,kend1
         kk = k-ksta+1
         do 6 j=jsta,jend1
         jj = j-jsta+1
         if (real(bcdata(jj,kk,ipp,l)) .lt. -1.e10) iflg = 1
   6     continue
         if (iflg.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  Stopping in bc2102:'
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)'  pressure data incorrectly set '
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
         if (filname.eq.'null') then
            pratio = bcdata(1,1,1,1)
            deltap = bcdata(1,1,1,2)
            rfreqp = bcdata(1,1,1,3)
            lref   = bcdata(1,1,1,4)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1008) real(pratio)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1009) real(deltap)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1010)
     .            real(twopi)*real(rfreqp)/real(lref)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2006)jsta,jend,ksta,kend
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2007)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2008)filname
         end if
c
      end if
c
      do 600 l=1,4
      do 600 k=ksta,kend1
      do 600 j=jsta,jend1
      qi0(j,k,l,3) = q(j,k,idim1,l)
      qi0(j,k,l,4) = qi0(j,k,l,3)
      bci(j,k,2) = 0.0
  600 continue
      do 602 k=ksta,kend1
      kk = k-ksta+1
      do 602 j=jsta,jend1
      jj = j-jsta+1
      pratio = bcdata(jj,kk,ip,1)
      deltap = bcdata(jj,kk,ip,2)
      rfreqp = bcdata(jj,kk,ip,3)
      lref   = bcdata(jj,kk,ip,4)
      pset   = (pratio + deltap*sin(twopi*rfreqp/lref*time))*p0
      qi0(j,k,5,3) = pset
      qi0(j,k,5,4) = pset
  602 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
          vi0(j,k,1,3) = vist3d(j,k,idim1)
          vi0(j,k,1,4) = vist3d(j,k,idim1)
  691   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 601 k=ksta,kend1
        do 601 j=jsta,jend1
          ti0(j,k,1,3) = tursav(j,k,idim1,1)
          ti0(j,k,2,3) = tursav(j,k,idim1,2)
          ti0(j,k,1,4) = tursav(j,k,idim1,1)
          ti0(j,k,2,4) = tursav(j,k,idim1,2)
  601   continue
      end if
      end if
      end if
c
 1001 format(' ','  j=   1  time-varying pressure          type 2102',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
 1002 format(' ','  j=jdim  time-varying pressure          type 2102',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
 1003 format(' ','  k=   1  time-varying pressure          type 2102',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
 1004 format(' ','  k=kdim  time-varying pressure          type 2102',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
 1005 format(' ','  i=   1  time-varying pressure          type 2102',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
 1006 format(' ','  i=idim  time-varying pressure          type 2102',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
 1007 format(11x,'P(t)/Pinf  = P1 + eps*sin(k*t)')
 1008 format(11x,'P1         = ',f8.4)
 1009 format(11x,'eps        = ',f8.4)
 1010 format(11x,'k          = ',f8.4)
c
 2001 format(' ','  j=   1  time-varying pressure          type 2102')
 2002 format(' ','  j=jdim  time-varying pressure          type 2102')
 2003 format(' ','  k=   1  time-varying pressure          type 2102')
 2004 format(' ','  k=kdim  time-varying pressure          type 2102')
 2005 format(' ','  i=   1  time-varying pressure          type 2102')
 2006 format(' ','  i=idim  time-varying pressure          type 2102')
 2007 format(11x,'P(t)/Pinf set from file:')
 2008 format('           ',a60)
c
      return 
      end
