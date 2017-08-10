c  ---------------------------------------------------------------------------
c  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
c  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
c  with point-matched, patched, overset, or embedded connectivities. Both
c  multigrid and mesh sequencing are available in time-accurate or
c  steady-state modes.
c
c  Copyright 2001 United States Government as represented by the Administrator
c  of the National Aeronautics and Space Administration. All Rights Reserved.
c 
c  The CFL3D platform is licensed under the Apache License, Version 2.0 
c  (the "License"); you may not use this file except in compliance with the 
c  License. You may obtain a copy of the License at 
c  http://www.apache.org/licenses/LICENSE-2.0. 
c 
c  Unless required by applicable law or agreed to in writing, software 
c  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
c  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
c  License for the specific language governing permissions and limitations 
c  under the License.
c  ---------------------------------------------------------------------------
c
      subroutine forceout(iseq,maxbl,maxgr,maxseg,nblock,iforce,
     .                    igridg,nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                    nbckdim,levelg,ibcinfo,jbcinfo,kbcinfo,
     .                    swett,clt,cdt,cxt,cyt,czt,cmxt,cmyt,cmzt,
     .                    cdpt,cdvt,swetw,clw,cdw,cxw,cyw,czw,cmxw,
     .                    cmyw,cmzw,cdpw,cdvw,ncycmax,myhost,myid,
     .                    mycomm,mblk2nd)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Write out forces and moments on individual blocks, as
c               well as a global force/moment summary
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
#   ifdef DBLE_PRECSN
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_DOUBLE_COMPLEX
#      else
#        define MY_MPI_REAL MPI_DOUBLE_PRECISION
#      endif
#   else
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_COMPLEX
#      else
#        define MY_MPI_REAL MPI_REAL
#      endif
#   endif
      dimension istat(MPI_STATUS_SIZE)
      dimension sfdat(11),sfsdat(11)
c
#endif
c
      dimension levelg(maxbl),igridg(maxbl),nbci0(maxbl),nbcidim(maxbl),
     .          nbcj0(maxbl),nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          ibcinfo(maxbl,maxseg,7,2),jbcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2),iforce(maxbl)
      dimension swett(maxbl),clt(maxbl),cdt(maxbl),cxt(maxbl),
     .          cyt(maxbl),czt(maxbl),cmxt(maxbl),cmyt(maxbl),
     .          cmzt(maxbl),cdpt(maxbl),cdvt(maxbl)
      dimension swetw(ncycmax),clw(ncycmax),cdw(ncycmax),cxw(ncycmax),
     .          cyw(ncycmax),czw(ncycmax),cmxw(ncycmax),cmyw(ncycmax),
     .          cmzw(ncycmax),cdpw(ncycmax),cdvw(ncycmax)
      dimension mblk2nd(maxbl)
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
c
#if defined DIST_MPI
c     initialize summary array for MPI code
      do n=1,11
      sfsdat(n) = 0.
      end do
c
c     set baseline tag values
c
      itag_sfdat = 1
#endif
c     individual block summaries
c
c
      if (myid.eq.myhost) then
         write(11,2311)
 2311    format(/40h***** FORCE AND MOMENT SUMMARIES - FINAL,
     .           32h TIME STEP/MULTIGRID CYCLE *****)
      end if
c
      do 20 n=1,nblock
      ifor  = iforce(n)
      levg = levelg(n)
      if (levg.eq.lglobal .and. ncyc1(iseq).gt.0 .and. ifor.gt.0) then
      if (myid.eq.myhost) then
         ifo   = ifor/100
         jfo   = (ifor - ifo*100)/10
         kfo   = (ifor - ifo*100 - jfo*10)
         write(11,677) n,igridg(n)
         if (ifo .eq.1) then
            write(11,620)
            do 1720 nseg = 1,nbci0(n)
            jsta = ibcinfo(n,nseg,2,1)
            jend = ibcinfo(n,nseg,3,1)
            ksta = ibcinfo(n,nseg,4,1)
            kend = ibcinfo(n,nseg,5,1)
            if(ibcinfo(n,nseg,6,1) .gt. 0)
     .      write(11,720) nseg,jsta,jend,ksta,kend
 1720       continue
         end if
         if (ifo .eq.2) then
            write(11,621)
            do 1721 nseg = 1,nbcidim(n)
            jsta = ibcinfo(n,nseg,2,2)
            jend = ibcinfo(n,nseg,3,2)
            ksta = ibcinfo(n,nseg,4,2)
            kend = ibcinfo(n,nseg,5,2)
            if(ibcinfo(n,nseg,6,2) .gt. 0)
     .      write(11,720) nseg,jsta,jend,ksta,kend
 1721       continue
         end if
         if (ifo .eq.3) then
            write(11,620)
            do 1730 nseg = 1,nbci0(n)
            jsta = ibcinfo(n,nseg,2,1)
            jend = ibcinfo(n,nseg,3,1)
            ksta = ibcinfo(n,nseg,4,1)
            kend = ibcinfo(n,nseg,5,1)
            if(ibcinfo(n,nseg,6,1) .gt. 0)
     .      write(11,720) nseg,jsta,jend,ksta,kend
 1730       continue
            write(11,621)
            do 1731 nseg = 1,nbcidim(n)
            jsta = ibcinfo(n,nseg,2,2)
            jend = ibcinfo(n,nseg,3,2)
            ksta = ibcinfo(n,nseg,4,2)
            kend = ibcinfo(n,nseg,5,2)
            if(ibcinfo(n,nseg,6,2) .gt. 0)
     .      write(11,720) nseg,jsta,jend,ksta,kend
 1731       continue
         end if
         if (jfo .eq.1) then
            write(11,623)
            do 1722 nseg = 1,nbcj0(n)
            ista = jbcinfo(n,nseg,2,1)
            iend = jbcinfo(n,nseg,3,1)
            ksta = jbcinfo(n,nseg,4,1)
            kend = jbcinfo(n,nseg,5,1)
            if(jbcinfo(n,nseg,6,1) .gt. 0)
     .      write(11,721) nseg,ista,iend,ksta,kend
 1722       continue
         end if
         if (jfo .eq.2) then
            write(11,624)
            do 1723 nseg = 1,nbcjdim(n)
            ista = jbcinfo(n,nseg,2,2)
            iend = jbcinfo(n,nseg,3,2)
            ksta = jbcinfo(n,nseg,4,2)
            kend = jbcinfo(n,nseg,5,2)
            if(jbcinfo(n,nseg,6,2) .gt. 0)
     .      write(11,721) nseg,ista,iend,ksta,kend
 1723       continue
         end if
         if (jfo .eq.3) then
            write(11,623)
            do 1732 nseg = 1,nbcj0(n)
            ista = jbcinfo(n,nseg,2,1)
            iend = jbcinfo(n,nseg,3,1)
            ksta = jbcinfo(n,nseg,4,1)
            kend = jbcinfo(n,nseg,5,1)
            if(jbcinfo(n,nseg,6,1) .gt. 0)
     .      write(11,721) nseg,ista,iend,ksta,kend
 1732       continue
            write(11,624)
            do 1733 nseg = 1,nbcjdim(n)
            ista = jbcinfo(n,nseg,2,2)
            iend = jbcinfo(n,nseg,3,2)
            ksta = jbcinfo(n,nseg,4,2)
            kend = jbcinfo(n,nseg,5,2)
            if(jbcinfo(n,nseg,6,2) .gt. 0)
     .      write(11,721) nseg,ista,iend,ksta,kend
 1733       continue
         end if
         if (kfo .eq.1) then
            write(11,626)
            do 1724 nseg = 1,nbck0(n)
            ista = kbcinfo(n,nseg,2,1)
            iend = kbcinfo(n,nseg,3,1)
            jsta = kbcinfo(n,nseg,4,1)
            jend = kbcinfo(n,nseg,5,1)
            if(kbcinfo(n,nseg,6,1) .gt. 0)
     .      write(11,722) nseg,ista,iend,jsta,jend
 1724       continue
         end if
         if (kfo .eq.2) then
            write(11,627)
            do 1725 nseg = 1,nbckdim(n)
            ista = kbcinfo(n,nseg,2,2)
            iend = kbcinfo(n,nseg,3,2)
            jsta = kbcinfo(n,nseg,4,2)
            jend = kbcinfo(n,nseg,5,2)
            if(kbcinfo(n,nseg,6,2) .gt. 0)
     .      write(11,722) nseg,ista,iend,jsta,jend
 1725       continue
         end if
         if (kfo .eq.3) then
            write(11,626)
            do 1734 nseg = 1,nbck0(n)
            ista = kbcinfo(n,nseg,2,1)
            iend = kbcinfo(n,nseg,3,1)
            jsta = kbcinfo(n,nseg,4,1)
            jend = kbcinfo(n,nseg,5,1)
            if(kbcinfo(n,nseg,6,1) .gt. 0)
     .      write(11,722) nseg,ista,iend,jsta,jend
 1734       continue
            write(11,627)
            do 1735 nseg = 1,nbckdim(n)
            ista = kbcinfo(n,nseg,2,2)
            iend = kbcinfo(n,nseg,3,2)
            jsta = kbcinfo(n,nseg,4,2)
            jend = kbcinfo(n,nseg,5,2)
            if(kbcinfo(n,nseg,6,2) .gt. 0)
     .      write(11,722) nseg,ista,iend,jsta,jend
 1735       continue
         end if
c
      end if
#if defined DIST_MPI
c
      if (myid.eq.mblk2nd(n)) then
         sfdat(1) = clt(n)
         sfdat(2) = cdt(n)
         sfdat(3) = cdpt(n)
         sfdat(4) = cdvt(n)
         sfdat(5) = czt(n)
         sfdat(6) = cyt(n)
         sfdat(7) = cxt(n)
         sfdat(8) = swett(n)
         sfdat(9) = cmyt(n)
         sfdat(10) = cmxt(n)
         sfdat(11) = cmzt(n)
         mytag = itag_sfdat + n
         call MPI_Send(sfdat,11,MY_MPI_REAL,myhost,
     .                 mytag,mycomm,ierr)
      end if
c
      if (myid.eq.myhost) then
         nd_srce = mblk2nd(n)
         mytag = itag_sfdat + n
         call MPI_Recv(sfdat,11,MY_MPI_REAL,nd_srce,
     .                 mytag,mycomm,istat,ierr)
         write(11,603)
         write(11,653) (real(sfdat(kk)),kk=1,4)
         write(11,604)
         write(11,653) (real(sfdat(kk)),kk=5,8)
         write(11,605)
         write(11,653) (real(sfdat(kk)),kk=9,11)
         do kk=1,11
         sfsdat(kk) = sfsdat(kk) + sfdat(kk)
         end do
      end if
#else
         write(11,603)
         write(11,653) real(clt(n)),real(cdt(n)),real(cdpt(n)),
     .                 real(cdvt(n))
         write(11,604)
         write(11,653) real(czt(n)),real(cyt(n)),real(cxt(n)),
     .                 real(swett(n))
         write(11,605) 
         write(11,653) real(cmyt(n)),real(cmxt(n)),real(cmzt(n))
#endif
      end if
   20 continue
c
c     summary of global blocks
c
#if defined DIST_MPI
      if (myid.eq.myhost) then
      write(11,777)
      write(11,703)
      write(11,653) (real(sfsdat(kk)),kk=1,4)
      write(11,704)
      write(11,653) (real(sfsdat(kk)),kk=5,8)
      write(11,705)
      write(11,653) (real(sfsdat(kk)),kk=9,11)
      write(11,*)
      end if
#else
      write(11,777)
      write(11,703)
      write(11,653) real(clw(nres)),real(cdw(nres)),
     .              real(cdpw(nres)),real(cdvw(nres))
      write(11,704)
      write(11,653) real(czw(nres)),real(cyw(nres)),
     .              real(cxw(nres)),real(swetw(nres))
      write(11,705)
      write(11,653) real(cmyw(nres)),real(cmxw(nres)),real(cmzw(nres))
      write(11,*)
#endif
c
  677 format(/1x,37hSUMMARY OF FORCES AND MOMENTS - BLOCK,i6,
     .       6h (GRID,i6,1h),/)
  620 format(1x,'forces computed on i=1 surface segment(s):')
  621 format(1x,'forces computed on i=idim surface segment(s):')
  622 format(1x,'forces computed on i=1 and i=idim surface segment(s):')
  623 format(1x,'forces computed on j=1 surface segment(s):')
  624 format(1x,'forces computed on j=jdim surface segment(s):')
  625 format(1x,'forces computed on j=1 and j=jdim surface segment(s):')
  626 format(1x,'forces computed on k=1 surface segment(s):')
  627 format(1x,'forces computed on k=kdim surface segment(s):')
  628 format(1x,'forces computed on k=1 and k=kdim surface segment(s):')
  603 format(/9x,4hCL-b,14x,4hCD-b,14x,5hCDp-b,13x,5hCDv-b)
  604 format(9x,4hCZ-b,14x,4hCY-b,14x,
     .       4hCX-b,10x,11hwetted area)
  605 format(9x,5hCMY-b,13x,5hCMX-b,13x,5hCMZ-b)
  777 format(/1x,49hSUMMARY OF FORCES AND MOMENTS - ALL GLOBAL BLOCKS)
  703 format(/10x,2hCL,16x,2hCD,15x,3hCDp,15x,3hCDv)
  704 format(10x,2hCZ,16x,2hCY,15x,2hCX,12x,11hwetted area)
  705 format(10x,3hCMY,14x,3hCMX,15x,3hCMZ)
  720 format(4x,'segment',i5,'  j=',i5,',',i5,'  k=',i5,',',i5)
  721 format(4x,'segment',i5,'  i=',i5,',',i5,'  k=',i5,',',i5)
  722 format(4x,'segment',i5,'  i=',i5,',',i5,'  j=',i5,',',i5)
  653 format(1x,e18.11,4(1x,e18.11))
c
      return
      end
