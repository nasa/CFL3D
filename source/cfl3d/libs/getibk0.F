      subroutine getibk0(jdim,kdim,idim,nbl,itotb,itoti,maxbl,
     .                   lig,lbg,ibpntsg,iipntsg,nou,bou,ibufdim,nbuf,
     .                   ierrflg,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Read the output from MAGGIE (not the grids); only the
c     data needed to evaluate sizing requirements are stored at this
c     time. (this is a modified version of subroutine getibk)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      real dyint,dzint,dxint,blank
c
      dimension nou(nbuf)
      dimension intpts(4)
      dimension lig(maxbl),lbg(maxbl),ibpntsg(maxbl,4),iipntsg(maxbl)
c
      common /params/ lmaxgr,lmaxbl,lmxseg,lmaxcs,lnplts,lmxbli,lmaxxe,
     .                lnsub1,lintmx,lmxxe,liitot,isum,lncycm,
     .                isum_n,lminnode,isumi,isumi_n,lmptch,
     .                lmsub1,lintmax,libufdim,lnbuf,llbcprd,
     .                llbcemb,llbcrad,lnmds,lmaxaes,lnslave,lmxsegdg,
     .                lnmaster
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
c
      write(66,101)nbl
  101 format(1x,37hreading overlap information for block,1x,i3)
c
      read(21) jchk,kchk,lchk
c
      if (jchk.ne.jdim1 .or. kchk.ne.kdim1 .or. lchk.ne.idim1) then
         write(66,*) '  mismatch in indices.....stopping in getibk'
         write(66,*) '  jdim,kdim,idim=',jdim,kdim,idim
         write(66,*) '  jmax,kmax,lmax=',jchk,kchk,lchk
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      read(21) ibpnts,intpts,iipnts,idum,idum
      write(66,*)'  ibpnts,intpts,iipnts = ',ibpnts,intpts,iipnts
      ibpntsg(nbl,1) = intpts(1)
      ibpntsg(nbl,2) = intpts(2)
      ibpntsg(nbl,3) = intpts(3)
      ibpntsg(nbl,4) = intpts(4)
      iipntsg(nbl)   = iipnts
c
      lsta = lig(nbl)
      lend = lsta+iipnts-1
      read(21)(jji,kki,iii,dyint,dzint,dxint,
     .         l=lsta,lend)
c
      lsta = lbg(nbl)
      lend = lsta+ibpnts-1
c
      read(21)(jjb,kkb,iib,ibc,l=lsta,lend)
c
      read(21)(((blank,j=1,jdim1),k=1,kdim1),i=1,idim1)
c
      if (nbl.lt.maxbl) then
         lig(nbl+1) = lig(nbl)+iipntsg(nbl)
         lbg(nbl+1) = lbg(nbl)+ibpnts
      end if
c
      itoti = itoti + iipntsg(nbl)
      itotb = itotb + ibpnts
c
      itotmx = max(itoti,itotb)
      liitot = max(itotmx,liitot)
c
      return
      end
