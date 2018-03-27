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
      subroutine resp(nbl,ntime,ntt0,jdim,kdim,idim,q,sj,sk,si,vol,
     .                x,y,z,res,wk0,vmuk,vmuj,vmui,wk,nwork,nblstag,
     .                nblendg,nblstat,nblendt,nptsr,nptsrb,nnegb,
     .                bcj,bck,bci,blank,nt,sumn,
     .                negn,rmst,xtbj,xtbk,xtbi,qc0,dqc0,iseq,
     .                qj0,qk0,qi0,maxbl,maxgr,maxseg,rms,rmsb,
     .                rmstb,clw,cdw,cdpw,cdvw,cxw,cyw,
     .                czw,cmxw,cmyw,cmzw,n_clcd,clcd,nblocks_clcd,
     .                blocks_clcd,chdw,swetw,fmdotw,
     .                cfttotw,cftmomw,cftpw,cftvw,rmstr,
     .                nneg,ihstry,ngrid,nblg,iemg,levelg,
     .                iviscg,itrans,irotat,iforce,swett,clt,cdt,
     .                cxt,cyt,czt,cmxt,cmyt,cmzt,cdpt,cdvt,nbci0,
     .                nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,ibcinfo,
     .                jbcinfo,kbcinfo,resmx,imx,jmx,kmx,vormax,
     .                ivmax,jvmax,kvmax,sx,sy,sz,stot,pav,ptav,
     .                tav,ttav,xmav,fmdot,cfxp,cfyp,cfzp,cfdp,
     .                cflp,cftp,cfxv,cfyv,cfzv,cfdv,cflv,cftv,
     .                cfxmom,cfymom,cfzmom,cfdmom,cflmom,cftmom,
     .                cfxtot,cfytot,cfztot,cfdtot,cfltot,cfttot,
     .                ncs,icsinfo,myid,myhost,mycomm,mblk2nd,nou,
     .                bou,nbuf,ibufdim,ncycmax,maxcs,thetay,iadvance,
     .                idefrm,igridg,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute and print residuals and sum the forces.
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
      dimension fcoef(40)
c
#endif
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf), resa(5), rmsba(5)
      dimension wk(nwork),mblk2nd(maxbl),iadvance(maxbl)
      dimension res(jdim*kdim*(idim-1),5),qc0(jdim,kdim,idim-1,5),
     .          dqc0(jdim,kdim,idim-1,5)
      dimension wk0(idim*jdim*22),blank(jdim,kdim,idim)
      dimension vmuk(jdim-1,idim-1,2),vmuj(kdim-1,idim-1,2),
     .          vmui(jdim-1,kdim-1,2)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension x(jdim*kdim*idim),y(jdim*kdim*idim),z(jdim*kdim*idim)
      dimension q(jdim*kdim*idim,5),vol(jdim*kdim*(idim-1))
      dimension si(jdim*kdim*idim,5),sj(jdim*kdim*(idim-1),5),
     .          sk(jdim*kdim*(idim-1),5)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2)
      dimension qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),
     .          qi0(jdim,kdim,5,4)
      integer blocks_clcd
      dimension rms(ncycmax),clw(ncycmax),
     .          cdw(ncycmax),cdpw(ncycmax),cdvw(ncycmax),
     .          cxw(ncycmax),cyw(ncycmax),czw(ncycmax),
     .          cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),
     .          clcd(2,n_clcd,ncycmax), blocks_clcd(2,nblocks_clcd),
     .          chdw(ncycmax),swetw(ncycmax),
     .          fmdotw(ncycmax),cfttotw(ncycmax),
     .          cftmomw(ncycmax),cftpw(ncycmax),cftvw(ncycmax),
     .          rmstr(ncycmax,nummem),nneg(ncycmax,nummem)
      dimension nblg(maxgr),iemg(maxgr),levelg(maxbl),iviscg(maxbl,3),
     .          itrans(maxbl),irotat(maxbl),idefrm(maxbl),iforce(maxbl),
     .          thetay(maxbl)
      dimension swett(maxbl),clt(maxbl),cdt(maxbl),cxt(maxbl),
     .          cyt(maxbl),czt(maxbl),cmxt(maxbl),cmyt(maxbl),
     .          cmzt(maxbl),cdpt(maxbl),cdvt(maxbl),nbci0(maxbl),
     .          nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2),
     .          resmx(maxbl),imx(maxbl),jmx(maxbl),kmx(maxbl),
     .          vormax(maxbl),ivmax(maxbl),jvmax(maxbl),kvmax(maxbl)
      dimension sx(maxcs),sy(maxcs),sz(maxcs),stot(maxcs),
     .          pav(maxcs),ptav(maxcs),tav(maxcs),ttav(maxcs),
     .          xmav(maxcs),fmdot(maxcs),
     .          cfxp(maxcs),cfyp(maxcs),cfzp(maxcs),
     .          cfdp(maxcs),cflp(maxcs),cftp(maxcs),
     .          cfxv(maxcs),cfyv(maxcs),cfzv(maxcs),
     .          cfdv(maxcs),cflv(maxcs),cftv(maxcs),
     .          cfxmom(maxcs),cfymom(maxcs),cfzmom(maxcs),
     .          cfdmom(maxcs),cflmom(maxcs),cftmom(maxcs),
     .          cfxtot(maxcs),cfytot(maxcs),cfztot(maxcs),
     .          cfdtot(maxcs),cfltot(maxcs),cfttot(maxcs),
     .          icsinfo(maxcs,9)
      dimension nnegb(nummem),rmstb(nummem),sumn(nummem),negn(nummem)
c 
      dimension igridg(maxbl)
      common /twod/ i2d
      common /ginfo2/ lq2avg,iskip_blocks,inc_2d(3),inc_coarse(3)
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
      common /account/ iaccnt,ioutsub
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /wrbl/ nwrest
      common /drag/ cdv,cdp
      common /conversion/ radtodeg
      common /subit/ cltsub,cdtsub,cxtsub,cytsub,cztsub,
     .        cmxtsub,cmytsub,cmztsub,cdptsub,cdvtsub,
     .        sxsub,sysub,szsub,stotsub,fmdotsub,cfxpsub,
     .        cfypsub,cfzpsub,cflpsub,cfdpsub,cftpsub,cfxvsub,
     .        cfyvsub,cfzvsub,cflvsub,cfdvsub,cftvsub,cfxmomsub,
     .        cfymomsub,cfzmomsub,cflmomsub,cfdmomsub,cftmomsub,
     .        cfxtotsub,cfytotsub,cfztotsub,cfdtotsub,cfltotsub,
     .        cfttotsub 
      common /igrdtyp/ ip3dgrd,ialph
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
c
c*************************************************************
c     global iteration convergence - monitor density residual,
c     force/moment coefficients (and turbulence residuals,
c     if appropriate), and control surface data
c*************************************************************
c
      nres = ntt
      iflag = 0
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset = maxbl
      itag_rms   = 1
      itag_fcoef = itag_rms + ioffset
c
#endif
      if (real(dt).lt.0. .and.   nt.eq.1) then
         iflag = 1
         resd = 0.
         if (myid.eq.mblk2nd(nbl) .and. iadvance(nbl).ge.0) then
            call l2norm(nbl,icyc,resd,0,jdim,kdim,idim,res,vol) 
            if ((real(resd+resd) .eq. real(resd) .and. 
     .         abs(real(resd)) .gt. 1.e-20) .or. 
     .         .not.(real(resd).lt.abs(real(resd)) .or.
     .         real(resd).ge.0.e0)) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1650) nbl,icyc
 1650          format(' NaN detected after residual evaluation, block',
     .         i4,' cycle',i5)
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
#if defined DIST_MPI
c
         nd_srce = mblk2nd(nbl)
         mytag = itag_rms + nbl
         if (myid.eq.mblk2nd(nbl)) then
            call MPI_Send (resd, 1, MY_MPI_REAL, myhost,
     .                     mytag, mycomm, ierr)
         end if
         if (myid.eq.myhost) then
            call MPI_Recv (resd, 1, MY_MPI_REAL, nd_srce,
     .                     mytag, mycomm, istat, ierr)
         end if
c
#endif
      end if
c
      if (real(dt).gt.0. .and. icyc.eq.ioutsub) then
         iflag = 1
c        l2-norm of density residual EXCLUDING time terms
         resd = 0.
         if (myid.eq.mblk2nd(nbl) .and. iadvance(nbl).ge.0) then
            call l2norm2(nbl,icyc,resd,0,jdim,kdim,idim,res,vol,
     .                   qc0,dqc0,q,blank)
            if ((real(resd+resd) .eq. real(resd) .and.
     .         abs(real(resd)) .gt. 1.e-20) .or. 
     .         .not.(real(resd).lt.abs(real(resd)) .or.
     .         real(resd).ge.0.e0)) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1651) nbl,nt
 1651          format('NaN detected after residual evaluation, block',
     .         i4,' time step',i5)
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
#if defined DIST_MPI
c
         nd_srce = mblk2nd(nbl)
         mytag = itag_rms + nbl
         if (myid.eq.mblk2nd(nbl)) then
            call MPI_Send (resd, 1, MY_MPI_REAL, myhost,
     .                     mytag, mycomm, ierr)
         end if
         if (myid.eq.myhost) then
            call MPI_Recv (resd, 1, MY_MPI_REAL, nd_srce,
     .                     mytag, mycomm, istat, ierr)
         end if
c
#endif
      end if
c
      if (iflag.eq.1) then
c
c        density residual
c
         if (myid.eq.myhost .and. iadvance(nbl).ge.0) then
            rmssum    = 0.
            if(nptsr.gt.0)
     .      rmssum = rms(ntt)*rms(ntt)*float(nptsr)
            t1    = float((jdim-1)*(kdim-1)*(idim-1))
            rmst  = sqrt(resd/t1)
            nptsr = nptsr + (jdim-1)*(kdim-1)*(idim-1)
            rmssum       = rmssum + resd
            rms(ntt)    = sqrt(rmssum/float(nptsr))
         end if
c
c        turbulence residual
c
         if ((iviscg(nbl,1).gt.2 .or. iviscg(nbl,2).gt.2
     .      .or. iviscg(nbl,3).gt.2) .and. iadvance(nbl).ge.0) then
c
#if defined DIST_MPI
            nd_srce = mblk2nd(nbl)
            mytag = itag_fcoef + nbl
c
            if (myid.eq.mblk2nd(nbl)) then
               call MPI_Send (sumn, nummem, MY_MPI_REAL,
     .                        myhost, mytag, mycomm, ierr)
               call MPI_Send (negn, nummem, MPI_INTEGER,
     .                        myhost, mytag, mycomm, ierr)
            end if

            if (myid.eq.myhost) then
               call MPI_Recv (sumn, nummem, MY_MPI_REAL,
     .                        nd_srce, mytag, mycomm, istat, ierr)
               call MPI_Recv (negn, nummem, MPI_INTEGER,
     .                        nd_srce, mytag, mycomm, istat, ierr)
            end if
c
#endif
            if (myid.eq.myhost) then
              do l=1,nummem
               rmstx  = 0.
               if (nptsr.gt.0) then
                  nptsr = nptsr - jdim1*kdim1*idim1
                  rmstx = rmstr(nres,l)*rmstr(nres,l)*float(nptsr)
               end if
               nptsr        = nptsr + jdim1*kdim1*idim1
               rmstx        = rmstx 
     .                      + sumn(l)*sumn(l)*float(jdim1*kdim1*idim1)
               rmstr(nres,l) = sqrt(  rmstx / float( nptsr) )
               nneg(nres,l)  = nneg(nres,l) + negn(l)
              enddo
            end if
         end if
c
c        forces and moments:
c
         if (iforce(nbl).gt.0) then
c
            ifor  = iforce(nbl)
            ifo   = ifor/100
            jfo   = (ifor - ifo*100)/10
            kfo   = (ifor - ifo*100 - jfo*10)
            icall = 0
            iuns = max(itrans(nbl),irotat(nbl),idefrm(nbl))
c
            if (myid.eq.mblk2nd(nbl)) then
               call force(jdim,kdim,idim,x,y,z,sk,sj,si,cl,cd,
     .                    czz,cyy,cxx,cmy,cmx,cmz,chd,swet,icall,
     .                    q(1,2),q(1,3),q(1,4),vmuk,vmuj,vmui,vol,
     .                    ifo,jfo,kfo,bcj,bck,bci,blank,nbl,xtbj,
     .                    xtbk,xtbi,iuns,qj0,qk0,qi0,nbci0,nbcj0,
     .                    nbck0,nbcidim,nbcjdim,nbckdim,ibcinfo,
     .                    jbcinfo,kbcinfo,maxbl,maxseg)
c
               if (ialph.eq.0) then
                  yref    = cref
                  zref    = bref
                  xref    = bref
               else
                  yref    = bref
                  zref    = cref
                  xref    = bref
               end if
               clt(nbl)   = cl/sref
               cmyt(nbl)  = cmy/(sref*yref)
               cmxt(nbl)  = cmx/(sref*xref)
               cmzt(nbl)  = cmz/(sref*zref)
               cdt(nbl)   = cd/sref
               czt(nbl)   = czz/sref
               cyt(nbl)   = cyy/sref
               cxt(nbl)   = cxx/sref
               cdpt(nbl)  = cdp/sref
               cdvt(nbl)  = cdv/sref
               swett(nbl) = swet
            end if
#if defined DIST_MPI
c
            nd_srce = mblk2nd(nbl)
            mytag = itag_fcoef + nbl
c
            if (myid.eq.mblk2nd(nbl)) then
               fcoef(1)  = clt(nbl)
               fcoef(2)  = cdt(nbl)
               fcoef(3)  = cyt(nbl)
               fcoef(4)  = cmyt(nbl)
               fcoef(5)  = cxt(nbl)
               fcoef(6)  = czt(nbl)
               fcoef(7)  = cmxt(nbl)
               fcoef(8)  = cmzt(nbl)
               fcoef(9)  = swett(nbl)
               fcoef(10) = cdpt(nbl)
               fcoef(11) = cdvt(nbl)
               fcoef(12) = chd
               call MPI_Send (fcoef, 12, MY_MPI_REAL,
     .                        myhost, mytag, mycomm, ierr)
            end if
            if (myid.eq.myhost) then
               call MPI_Recv (fcoef, 12, MY_MPI_REAL,
     .                        nd_srce, mytag, mycomm, istat, ierr)
               clt(nbl)   = fcoef(1)
               cdt(nbl)   = fcoef(2)
               cyt(nbl)   = fcoef(3)
               cmyt(nbl)  = fcoef(4)
               cxt(nbl)   = fcoef(5)
               czt(nbl)   = fcoef(6)
               cmxt(nbl)  = fcoef(7)
               cmzt(nbl)  = fcoef(8)
               swett(nbl) = fcoef(9)
               cdpt(nbl)  = fcoef(10)
               cdvt(nbl)  = fcoef(11)
               chd        = fcoef(12)
            end if
c
#endif
c
c           sum contributions of blocks on global level
c
            if (myid.eq.myhost) then
               if (level.eq.lglobal) then
                  chdw(nres)  = chdw(nres)  + chd       
                  swetw(nres) = swetw(nres) + swett(nbl)
                  clw(nres)   = clw(nres)   + clt(nbl)
                  cmyw(nres)  = cmyw(nres)  + cmyt(nbl) 
                  cdw(nres)   = cdw(nres)   + cdt(nbl)
                  cmxw(nres)  = cmxw(nres)  + cmxt(nbl) 
                  cmzw(nres)  = cmzw(nres)  + cmzt(nbl) 
                  czw(nres)   = czw(nres)   + czt(nbl)  
                  cyw(nres)   = cyw(nres)   + cyt(nbl)
                  cxw(nres)   = cxw(nres)   + cxt(nbl)
                  cdpw(nres)  = cdpw(nres)  + cdpt(nbl)
                  cdvw(nres)  = cdvw(nres)  + cdvt(nbl)

                  if (iclcd .eq. 1 .or. iclcd .eq. 2) then
                     do nb = 1, nblocks_clcd
                        if (igridg(nbl) .eq. blocks_clcd(2,nb)) then
                           nn = blocks_clcd(1,nb)
                           clcd(1,nn,nres) = clcd(1,nn,nres) + clt(nbl)
                           clcd(2,nn,nres) = clcd(2,nn,nres) + cdt(nbl)
                        end if
                     end do
                  end if
               end if
            end if
         end if
c
c        control surface data
c
         do ics = 1,ncs
            if (nbl.eq.icsinfo(ics,1) .and. levelg(nbl).ge.lglobal) then
               iuns = max(itrans(nbl),irotat(nbl),idefrm(nbl))
               if (myid.eq.mblk2nd(nbl)) then
                  call csurf(jdim,kdim,idim,x,y,z,sk,sj,si,q,ics,
     .                       q(1,2),q(1,3),q(1,4),vmuk,vmuj,vmui,vol,
     .                       bcj,bck,bci,blank,xtbj,xtbk,xtbi,
     .                       iuns,ncs,icsinfo,sxi,syi,szi,stoti,
     .                       pavi,ptavi,tavi,ttavi,xmavi,fmdoti,
     .                       cfxpi,cfypi,cfzpi,cfdpi,cflpi,cftpi,
     .                       cfxvi,cfyvi,cfzvi,cfdvi,cflvi,cftvi,
     .                       cfxmomi,cfymomi,cfzmomi,cfdmomi,cflmomi,
     .                       cftmomi,cfxtoti,cfytoti,cfztoti,cfdtoti,
     .                       cfltoti,cfttoti,maxcs,qj0,qk0,qi0)
                  sx(ics)     = sxi
                  sy(ics)     = syi
                  sz(ics)     = szi
                  stot(ics)   = stoti
                  pav(ics)    = pavi
                  ptav(ics)   = ptavi
                  tav(ics)    = tavi
                  ttav(ics)   = ttavi
                  xmav(ics)   = xmavi
                  fmdot(ics)  = fmdoti
                  cfxp(ics)   = cfxpi
                  cfyp(ics)   = cfypi
                  cfzp(ics)   = cfzpi
                  cflp(ics)   = cflpi
                  cfdp(ics)   = cfdpi
                  cftp(ics)   = cftpi
                  cfxv(ics)   = cfxvi
                  cfyv(ics)   = cfyvi
                  cfzv(ics)   = cfzvi
                  cflv(ics)   = cflvi 
                  cfdv(ics)   = cfdvi
                  cftv(ics)   = cftvi
                  cfxmom(ics) = cfxmomi
                  cfymom(ics) = cfymomi
                  cfzmom(ics) = cfzmomi
                  cflmom(ics) = cflmomi
                  cfdmom(ics) = cfdmomi
                  cftmom(ics) = cftmomi
                  cfxtot(ics) = cfxtoti
                  cfytot(ics) = cfytoti
                  cfztot(ics) = cfztoti
                  cfdtot(ics) = cfdtoti
                  cfltot(ics) = cfltoti
                  cfttot(ics) = cfttoti
               end if
#if    defined DIST_MPI
c
               nd_srce = mblk2nd(nbl)
               mytag = itag_fcoef + nbl
c
               if (myid.eq.mblk2nd(nbl)) then
                  fcoef(1)  = sx(ics)
                  fcoef(2)  = sy(ics)
                  fcoef(3)  = sz(ics)
                  fcoef(4)  = stot(ics)
                  fcoef(5)  = pav(ics)
                  fcoef(6)  = ptav(ics)
                  fcoef(7)  = tav(ics)
                  fcoef(8)  = ttav(ics)
                  fcoef(9)  = xmav(ics)
                  fcoef(10) = fmdot(ics)
                  fcoef(11) = cfxp(ics)
                  fcoef(12) = cfyp(ics)
                  fcoef(13) = cfzp(ics)
                  fcoef(14) = cflp(ics)
                  fcoef(15) = cfdp(ics)
                  fcoef(16) = cftp(ics)
                  fcoef(17) = cfxv(ics)
                  fcoef(18) = cfyv(ics)
                  fcoef(19) = cfzv(ics)
                  fcoef(20) = cflv(ics)
                  fcoef(21) = cfdv(ics)
                  fcoef(22) = cftv(ics)
                  fcoef(23) = cfxmom(ics)
                  fcoef(24) = cfymom(ics)
                  fcoef(25) = cfzmom(ics)
                  fcoef(26) = cflmom(ics)
                  fcoef(27) = cfdmom(ics)
                  fcoef(28) = cftmom(ics)
                  fcoef(29) = cfxtot(ics)
                  fcoef(30) = cfytot(ics)
                  fcoef(31) = cfztot(ics)
                  fcoef(32) = cfdtot(ics)
                  fcoef(33) = cfltot(ics)
                  fcoef(34) = cfttot(ics)
                  call MPI_Send (fcoef, 34, MY_MPI_REAL,
     .                           myhost, mytag, mycomm, ierr)
               end if
               if (myid.eq.myhost) then
                  call MPI_Recv (fcoef, 34, MY_MPI_REAL,
     .                           nd_srce, mytag, mycomm, istat, ierr)
                  sx(ics)      = fcoef(1)
                  sy(ics)      = fcoef(2)
                  sz(ics)      = fcoef(3)
                  stot(ics)    = fcoef(4)
                  pav(ics)     = fcoef(5)
                  ptav(ics)    = fcoef(6)
                  tav(ics)     = fcoef(7)
                  ttav(ics)    = fcoef(8)
                  xmav(ics)    = fcoef(9)
                  fmdot(ics)   = fcoef(10)
                  cfxp(ics)    = fcoef(11)
                  cfyp(ics)    = fcoef(12)
                  cfzp(ics)    = fcoef(13)
                  cflp(ics)    = fcoef(14)
                  cfdp(ics)    = fcoef(15)
                  cftp(ics)    = fcoef(16)
                  cfxv(ics)    = fcoef(17)
                  cfyv(ics)    = fcoef(18)
                  cfzv(ics)    = fcoef(19)
                  cflv(ics)    = fcoef(20)
                  cfdv(ics)    = fcoef(21)
                  cftv(ics)    = fcoef(22)
                  cfxmom(ics)  = fcoef(23)
                  cfymom(ics)  = fcoef(24)
                  cfzmom(ics)  = fcoef(25)
                  cflmom(ics)  = fcoef(26)
                  cfdmom(ics)  = fcoef(27)
                  cftmom(ics)  = fcoef(28)
                  cfxtot(ics)  = fcoef(29)
                  cfytot(ics)  = fcoef(30)
                  cfztot(ics)  = fcoef(31)
                  cfdtot(ics)  = fcoef(32)
                  cfltot(ics)  = fcoef(33)
                  cfttot(ics)  = fcoef(34)
               end if
c
#endif
            end if
         end do
c
c        sum contributions over all control surfaces with inorm .ne. 0 
c        (global level only)
c
         if (myid.eq.myhost) then
            if (nbl.eq.nblendg .and. levelg(nbl).eq.lglobal) then
               sx(ncs+1)     = 0.e0
               sy(ncs+1)     = 0.e0
               sz(ncs+1)     = 0.e0
               stot(ncs+1)   = 0.e0
               fmdot(ncs+1)  = 0.e0
               cfxp(ncs+1)   = 0.e0
               cfyp(ncs+1)   = 0.e0
               cfzp(ncs+1)   = 0.e0
               cflp(ncs+1)   = 0.e0
               cfdp(ncs+1)   = 0.e0
               cfxv(ncs+1)   = 0.e0
               cfyv(ncs+1)   = 0.e0
               cfzv(ncs+1)   = 0.e0
               cflv(ncs+1)   = 0.e0
               cfdv(ncs+1)   = 0.e0
               cfxmom(ncs+1) = 0.e0
               cfymom(ncs+1) = 0.e0
               cfzmom(ncs+1) = 0.e0
               cflmom(ncs+1) = 0.e0
               cfdmom(ncs+1) = 0.e0
               cfxtot(ncs+1) = 0.e0
               cfytot(ncs+1) = 0.e0
               cfztot(ncs+1) = 0.e0
               cfdtot(ncs+1) = 0.e0
               cfltot(ncs+1) = 0.e0
               do ics=1,ncs
                  fnorm = float(icsinfo(ics,9))
                  nbl1 = icsinfo(ics,1)
                  if (real(fnorm).ne.0.e0 .and.
     .               levelg(nbl1).eq.lglobal) then
                     sx(ncs+1)     = sx(ncs+1)     + sx(ics)
                     sy(ncs+1)     = sy(ncs+1)     + sy(ics)
                     sz(ncs+1)     = sz(ncs+1)     + sz(ics)
                     stot(ncs+1)   = stot(ncs+1)   + stot(ics)
                     fmdot(ncs+1)  = fmdot(ncs+1)  + fmdot(ics)
                     cfxp(ncs+1)   = cfxp(ncs+1)   + cfxp(ics)
                     cfyp(ncs+1)   = cfyp(ncs+1)   + cfyp(ics)
                     cfzp(ncs+1)   = cfzp(ncs+1)   + cfzp(ics)
                     cflp(ncs+1)   = cflp(ncs+1)   + cflp(ics)
                     cfdp(ncs+1)   = cfdp(ncs+1)   + cfdp(ics)
                     cfxv(ncs+1)   = cfxv(ncs+1)   + cfxv(ics)
                     cfyv(ncs+1)   = cfyv(ncs+1)   + cfyv(ics)
                     cfzv(ncs+1)   = cfzv(ncs+1)   + cfzv(ics)
                     cflv(ncs+1)   = cflv(ncs+1)   + cflv(ics)
                     cfdv(ncs+1)   = cfdv(ncs+1)   + cfdv(ics)
                     cfxmom(ncs+1) = cfxmom(ncs+1) + cfxmom(ics)
                     cfymom(ncs+1) = cfymom(ncs+1) + cfymom(ics)
                     cfzmom(ncs+1) = cfzmom(ncs+1) + cfzmom(ics)
                     cflmom(ncs+1) = cflmom(ncs+1) + cflmom(ics)
                     cfdmom(ncs+1) = cfdmom(ncs+1) + cfdmom(ics)
                     cfxtot(ncs+1) = cfxtot(ncs+1) + cfxtot(ics)
                     cfytot(ncs+1) = cfytot(ncs+1) + cfytot(ics)
                     cfztot(ncs+1) = cfztot(ncs+1) + cfztot(ics)
                     cfltot(ncs+1) = cfltot(ncs+1) + cfltot(ics)
                     cfdtot(ncs+1) = cfdtot(ncs+1) + cfdtot(ics)
c
                     cftp(ncs+1)   = sqrt(cfxp(ncs+1)**2 
     .                                  + cfyp(ncs+1)**2 
     .                                  + cfzp(ncs+1)**2)
                     cftv(ncs+1)   = sqrt(cfxv(ncs+1)**2 
     .                                  + cfyv(ncs+1)**2 
     .                                 + cfzv(ncs+1)**2)
                     cftmom(ncs+1) = sqrt(cfxmom(ncs+1)**2 
     .                                  + cfymom(ncs+1)**2 
     .                                  + cfzmom(ncs+1)**2)
                     cfttot(ncs+1) = sqrt(cfxtot(ncs+1)**2 
     .                                  + cfytot(ncs+1)**2
     .                                  + cfztot(ncs+1)**2)
c
c                    store selected totals for convergence history
c
                     fmdotw(nres)  = fmdot(ncs+1)
                     cftmomw(nres) = cftmom(ncs+1)
                     cftpw(nres)   = cftp(ncs+1)
                     cftvw(nres)   = cftv(ncs+1)
                     cfttotw(nres) = cfttot(ncs+1)
                  end if
               end do
            end if
         end if
c
      end if
c
c     print residual and force coefficients
c
      if (myid.eq.myhost) then
         if (real(dt).lt.0. .and. nt.eq.1) then
            iflag = 0
            if (icyc.eq.1) iflag = 1
            if (icyc.eq.2 .and. nbl.eq.nblendg) iflag = 1
            if (icyc.eq.ncyc .and. nbl.eq.nblstat)  iflag = 1
            if (icyc/nwrest*nwrest.eq.icyc .and. nbl.eq.nblstat) 
     .         iflag = 1
            if ((icyc-1)/nwrest*nwrest.eq.(icyc-1) .and.
     .         nbl.eq.nblstag) iflag = 1
            if (iflag.eq.1) then
               write(11,1080)
            end if
            iflag = 0
            if (icyc.eq.1) iflag = 1
            if (icyc.ge.2 .and. nbl.eq.nblendg) iflag = 1
            if (icyc.eq.ncyc .and. nbl.ge.nblstag) iflag = 1
            if (icyc/nwrest*nwrest.eq.icyc .and. nbl.ge.nblstag)
     .         iflag = 1
            if (iflag.gt.0) then
               if (ialph .ne. 0) then
                  write(11,2000) level,nbl,nres,real(rmst),
     .            real(rms(nres)),real(clw(nres)),real(cdw(nres)),
     .            real(czw(nres))
               else
                  write(11,2000) level,nbl,nres,real(rmst),
     .            real(rms(nres)),real(clw(nres)),real(cdw(nres)),
     .            real(cyw(nres))
               end if
            end if
         end if
c
         if (real(dt).gt.0. .and. icyc.eq.ioutsub) then
            iflag = 0
            if (ioutsub.eq.1) then
               if (nt.eq.1) iflag = 1
               if (nt.eq.2 .and. nbl.eq.nblendg) iflag = 1
               if (nt.eq.ntstep .and. nbl.eq.nblstat) iflag = 1
               if (nt/nwrest*nwrest.eq.nt  .and. nbl.eq.nblstat)
     .            iflag = 1
               if (nt.gt.1 .and. (nt-1)/nwrest*nwrest.eq.(nt-1) .and.
     .            nbl.eq.nblstag) iflag = 1
            end if
            if (ioutsub.eq.ncyc) then
               if (nt.eq.1 .and. nbl.eq.nblstat) iflag = 1 
               if (nt.eq.2 .and. nbl.eq.nblendg) iflag = 1
               if (nt.eq.ntstep .and. nbl.eq.nblstat) iflag = 1
               if (nt/nwrest*nwrest.eq.nt  .and. nbl.eq.nblstat)
     .            iflag = 1
               if (nt.gt.1 .and. (nt-1)/nwrest*nwrest.eq.(nt-1) .and.
     .            nbl.eq.nblstag) iflag = 1
            end if
            if (iflag.gt.0) then
               if (iunst .eq. 0) then
                  write(11,2020)
               else
                  write(11,2040)
               end if
            end if
            iflag = 0
            if (ioutsub.eq.1) then
               if (nt.eq.1) iflag = 1
               if (nt.ge.2 .and. nbl.eq.nblendg) iflag = 1
               if (nt.eq.ntstep .and. nbl.ge.nblstag) iflag = 1
               if (nt/nwrest*nwrest.eq.nt .and. nbl.ge.nblstag)
     .            iflag = 1
            end if
            if (ioutsub.eq.ncyc) then
               if (nt.eq.1 .and. nbl.ge.nblstag) iflag = 1
               if (nt.ge.2 .and. nbl.eq.nblendg) iflag = 1
               if (nt.eq.ntstep .and. nbl.ge.nblstag) iflag = 1
               if (nt/nwrest*nwrest.eq.nt .and. nbl.ge.nblstag)
     .            iflag = 1
            end if
            if (iflag.gt.0) then 
               if (iunst .eq. 0) then
                  if (ialph .ne. 0) then
                     write(11,2000) level,nbl,nres,real(rmst),
     .               real(rms(nres)),real(clw(nres)),real(cdw(nres)),
     .               real(czw(nres))
                  else
                     write(11,2000) level,nbl,nres,real(rmst),
     .               real(rms(nres)),real(clw(nres)),real(cdw(nres)),
     .               real(cyw(nres))
                  end if
               else
                  alot = (alpha+thetay(nbl))*radtodeg
                  if (ialph .ne. 0) then
                     write(11,2000) level,nbl,nres,real(rmst),
     .               real(rms(nres)),real(clw(nres)),real(cdw(nres)),
     .               real(czw(nres)),real(time),real(alot)
                  else
                     write(11,2000) level,nbl,nres,real(rmst),
     .               real(rms(nres)),real(clw(nres)),real(cdw(nres)),
     .               real(cyw(nres)),real(time),real(alot)
                  end if
               end if
            end if
         end if
      end if
c
c     print maximum residual/location and, if turbulent,
c     maximum vorticity/location - print whenever restart
c     file is written and at end of computation
c
      iflag = 0
      if (real(dt).le.0. .and. nt.eq.1) then
         if (icyc/nwrest*nwrest.eq.icyc .or. icyc.eq.ncyc) iflag=1
      end if
      if (real(dt).gt.0. .and. icyc.eq.ioutsub) then
         if(nt/nwrest*nwrest.eq.nt .or. nt.eq.ntstep) iflag=1
      end if
c
      if (iflag.gt.0) then
         if (myid.eq.mblk2nd(nbl)) then
            call blkmax(jdim,kdim,idim,res,resmax,jm,km,im)
            resmx(nbl) = resmax
            imx(nbl)   = im
            jmx(nbl)   = jm
            kmx(nbl)   = km
         end if
#if defined DIST_MPI
c
         nd_srce = mblk2nd(nbl)
         mytag = itag_fcoef + nbl
c
         if (myid.eq.mblk2nd(nbl)) then
            fcoef(1) = resmax
            fcoef(2) = float(imx(nbl))
            fcoef(3) = float(jmx(nbl))
            fcoef(4) = float(kmx(nbl))
            fcoef(5) = vormax(nbl)
            fcoef(6) = float(ivmax(nbl))
            fcoef(7) = float(jvmax(nbl))
            fcoef(8) = float(kvmax(nbl))
            call MPI_Send (fcoef, 8, MY_MPI_REAL,
     .                     myhost, mytag, mycomm, ierr)
         end if
         if (myid.eq.myhost) then
            call MPI_Recv (fcoef, 8, MY_MPI_REAL,
     .                     nd_srce, mytag, mycomm, istat, ierr)
            resmx(nbl)  = fcoef(1)
            imx(nbl)    = int(fcoef(2))
            jmx(nbl)    = int(fcoef(3))
            kmx(nbl)    = int(fcoef(4))
            vormax(nbl) = fcoef(5)
            ivmax(nbl)  = int(fcoef(6))
            jvmax(nbl)  = int(fcoef(7))
            kvmax(nbl)  = int(fcoef(8))
         end if
c
#endif
      end if
c
      if (myid.eq.myhost) then
         iflag = 0
         if (real(dt).le.0. .and. nt.eq.1) then
            if (level.eq.lglobal .and. nbl.eq.nblendg) then
               if (icyc/nwrest*nwrest.eq.icyc .or. icyc.eq.ncyc) iflag=1
            end if
         end if
         if (real(dt).gt.0. .and. icyc.eq.ioutsub) then
            if(level.eq.lglobal .and. nbl.eq.nblendg) then
               if(nt/nwrest*nwrest.eq.nt .or. nt.eq.ntstep) iflag=1
            end if
         end if
c
         if (iflag.gt.0) then 
            do 8000 igrid=1,ngrid
            iskip = 0
            if (igrid.eq.1) iskip = 1
            nbl1 = nblg(igrid)
            iem  = iemg(igrid)
            if (iseq.ne.mseq .and. iem.gt.0) go to 8000
            if (iem.eq.0) nbl1 = nbl1 + (mseq-iseq)
            if (iskip.gt.0) then
               write(11,2060)nbl1,igrid,real(resmx(nbl1)),jmx(nbl1),
     .         kmx(nbl1),imx(nbl1)
            else
               write(11,2080)nbl1,igrid,real(resmx(nbl1)),jmx(nbl1),
     .         kmx(nbl1),imx(nbl1)
            end if
 8000       continue
            do 8050 igrid=1,ngrid
            iskip = 0
            if (igrid.eq.1) iskip = 1
            nbl1 = nblg(igrid)
            iem  = iemg(igrid)
            if (iseq.ne.mseq .and. iem.gt.0) go to 8050
            if (iem.eq.0) nbl1 = nbl1 + (mseq-iseq)
            if ((iviscg(nbl1,3).gt.1 .or. iviscg(nbl1,2).gt.1
     .         .or.iviscg(nbl1,1).gt.1) .and. iadvance(nbl).ge.0) then
               if (iskip.gt.0) then
                  write(11,3000)nbl1,igrid,real(vormax(nbl1)),
     .            jvmax(nbl1),kvmax(nbl1),ivmax(nbl1)
               else
                  write(11,3020)nbl1,igrid,real(vormax(nbl1)),
     .            jvmax(nbl1),kvmax(nbl1),ivmax(nbl1)
               end if
            end if
 8050      continue
         end if
      end if
c
c     output unsteady cp data for dynamic mesh cases
c     irite = 0 don't output data
c             1 output data with old headers
c             2 output data with new headers
c
      irite = 0
      iwk1  = 1
      iwk2  = iwk1 + jdim*kdim
      iwk3  = iwk2 + kdim*idim
      if (irite.gt.0 .and. iunst.gt.0) then
         call prntcp(jdim,kdim,idim,wk(iwk1),wk(iwk2),wk(iwk3),q,
     .               nbl,maxbl,maxseg,ibcinfo,jbcinfo,kbcinfo,nbci0,
     .               nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,thetay,
     .               mblk2nd,myid,myhost,mycomm,irite)
      end if
c
c************************************************************
c     sub-iteration convergence - monitor density residual,
c     force/moment coefficients (and turbulence residuals,
c     if appropriate), and control surface data
c************************************************************
c
      if (real(dt).ge.0. .and. ncyc.gt.1) then
c
c        l2-norm of full residual (INCLUDING time terms)
c
         resd = 0.
         if (myid.eq.mblk2nd(nbl) .and. iadvance(nbl).ge.0) then
            call l2norm(nbl,icyc,resd,0,jdim,kdim,idim,res,vol)
         end if
#if defined DIST_MPI
c
         nd_srce = mblk2nd(nbl)
         mytag = itag_rms + nbl
         if (myid.eq.mblk2nd(nbl)) then
            call MPI_Send (resd, 1, MY_MPI_REAL, myhost,
     .                     mytag, mycomm, ierr)
         end if
         if (myid.eq.myhost) then
            call MPI_Recv (resd, 1, MY_MPI_REAL, nd_srce,
     .                     mytag, mycomm, istat, ierr)
         end if
#endif
c
         if (myid.eq.myhost .and. iadvance(nbl).ge.0) then
            rmssum = 0.
            if (nptsrb.gt.0) rmssum = rmsb*rmsb*float(nptsrb)
            t1     = float(1*jdim1*kdim1*idim1)
c           rmstb  = sqrt(resd/t1)
            nptsrb = nptsrb + jdim1*kdim1*idim1
            rmssum = rmssum + resd
            rmsb   = sqrt(  rmssum / float( nptsrb) )
         end if

         if( isubit_r .ne. 0 ) then
c
c        l2-norm of full residual (INCLUDING time terms)
c
            resa = 0.
            if (myid.eq.mblk2nd(nbl) .and. iadvance(nbl).ge.0) then
               call l2normAll(nbl,icyc,resa,0,jdim,kdim,idim,res,vol)
            end if

#if defined DIST_MPI
            nd_srce = mblk2nd(nbl)
            mytag = itag_rms + nbl
            if (myid.eq.mblk2nd(nbl)) then
               call MPI_Send (resa, 5, MY_MPI_REAL, myhost,
     .              mytag, mycomm, ierr)
            end if
            if (myid.eq.myhost) then
               call MPI_Recv (resa, 5, MY_MPI_REAL, nd_srce,
     .              mytag, mycomm, istat, ierr)
            end if
#endif
c
            if (myid.eq.myhost .and. iadvance(nbl).ge.0) then
               do ivar = 1, 5
                  rmssum = 0.
                  if (nptsrb.gt.jdim1*kdim1*idim1) then
                     rmssum = rmsba(ivar)*rmsba(ivar)
     $                    *float(nptsrb-jdim1*kdim1*idim1)
                  end if
                  rmssum = rmssum + resa(ivar)
                  rmsba(ivar) = sqrt(  rmssum / float( nptsrb) )
               end do
            end if

         end if
            
c
c        turbulence residual
c
         if ((iviscg(nbl,1).gt.2 .or. iviscg(nbl,2).gt.2 .or.
     .       iviscg(nbl,3).gt.2) .and. iadvance(nbl).ge.0) then
c
#if defined DIST_MPI
            nd_srce = mblk2nd(nbl)
            mytag = itag_fcoef + nbl
c
            if (myid.eq.mblk2nd(nbl)) then
               call MPI_Send (sumn, nummem, MY_MPI_REAL,
     .                        myhost, mytag, mycomm, ierr)
               call MPI_Send (negn, nummem, MPI_INTEGER,
     .                        myhost, mytag, mycomm, ierr)
            end if

            if (myid.eq.myhost) then
               call MPI_Recv (sumn, nummem, MY_MPI_REAL,
     .                        nd_srce, mytag, mycomm, istat, ierr)
               call MPI_Recv (negn, nummem, MPI_INTEGER,
     .                        nd_srce, mytag, mycomm, istat, ierr)
            end if
#endif
            if (myid.eq.myhost) then
              do l=1,nummem
               rmsts1 = 0.
               if (nptsrb.gt.0) then
                  nptsrb = nptsrb - jdim1*kdim1*idim1
c                 if (nptsrb.gt.0) then
                     rmsts1 = rmstb(l)*rmstb(l)*float(nptsrb)
c                 end if
               end if
               nptsrb = nptsrb + jdim1*kdim1*idim1
               rmsts1 = rmsts1 
     .                + sumn(l)*sumn(l)*float(jdim1*kdim1*idim1)
               rmstb(l) = sqrt(  rmsts1 / float( nptsrb) )
               nnegb(l) = nnegb(l) + negn(l)
              enddo
            end if
         end if
c
c        compute forces and moments (global level only)
c
         if (level.eq.lglobal) then
c
            if (myid.eq.myhost) then
c
c              zero out values if this is the first block on the global level
               if (nbl.eq.nblstag) then
                  cltsub   = 0.
                  cmytsub  = 0.
                  cmxtsub  = 0.
                  cmztsub  = 0. 
                  cdtsub   = 0.
                  cztsub   = 0.
                  cytsub   = 0.
                  cxtsub   = 0.
                  cdptsub  = 0.
                  cdvtsub  = 0.
               end if
            end if
c
            if (iforce(nbl).gt.0) then
               ifor  = iforce(nbl)
               ifo   = ifor/100
               jfo   = (ifor - ifo*100)/10
               kfo   = (ifor - ifo*100 - jfo*10)
               icall = 0
               iuns = max(itrans(nbl),irotat(nbl),idefrm(nbl))
               if (myid.eq.mblk2nd(nbl)) then
                  call force(jdim,kdim,idim,x,y,z,sk,sj,si,cl,cd,czz,
     .                       cyy,cxx,cmy,cmx,cmz,chd,swet,icall,q(1,2),
     .                       q(1,3),q(1,4),vmuk,vmuj,vmui,vol,ifo,jfo,
     .                       kfo,bcj,bck,bci,blank,nbl,xtbj,xtbk,xtbi,
     .                       iuns,qj0,qk0,qi0,nbci0,nbcj0,
     .                       nbck0,nbcidim,nbcjdim,nbckdim,ibcinfo,
     .                       jbcinfo,kbcinfo,maxbl,maxseg)

               end if
#if defined DIST_MPI
c
               nd_srce = mblk2nd(nbl)
               mytag = itag_fcoef + nbl
c
               if (myid.eq.mblk2nd(nbl)) then
c
                  fcoef(1)  = cl
                  fcoef(2)  = cd
                  fcoef(3)  = cyy
                  fcoef(4)  = cmy
                  fcoef(5)  = cxx
                  fcoef(6)  = czz
                  fcoef(7)  = cmx
                  fcoef(8)  = cmz
                  fcoef(9)  = cdp
                  fcoef(10) = cdv
                  call MPI_Send (fcoef, 10, MY_MPI_REAL,
     .                           myhost, mytag, mycomm, ierr)
               end if
               if (myid.eq.myhost) then
                  call MPI_Recv (fcoef, 10, MY_MPI_REAL,
     .                           nd_srce, mytag, mycomm, istat, ierr)
                  cl   = fcoef(1)
                  cd   = fcoef(2)
                  cyy  = fcoef(3)
                  cmy  = fcoef(4)
                  cxx  = fcoef(5)
                  czz  = fcoef(6)
                  cmx  = fcoef(7)
                  cmz  = fcoef(8)
                  cdp  = fcoef(9)
                  cdv  = fcoef(10)
               end if
c
#endif
               if (myid.eq.myhost) then
c
c                 sum contributions of blocks on global level
c
                  if (ialph.eq.0) then
                     yref    = cref
                     zref    = bref
                     xref    = bref
                  else
                     yref    = bref
                     zref    = cref
                     xref    = bref
                  end if
                  cltsub   = cltsub   + cl/sref
                  cmytsub  = cmytsub  + cmy/(sref*yref)
                  cmxtsub  = cmxtsub  + cmx/(sref*xref)
                  cmztsub  = cmztsub  + cmz/(sref*zref)
                  cdtsub   = cdtsub   + cd/sref
                  cztsub   = cztsub   + czz/sref
                  cytsub   = cytsub   + cyy/sref
                  cxtsub   = cxtsub   + cxx/sref
                  cdptsub  = cdptsub  + cdp/sref
                  cdvtsub  = cdvtsub  + cdv/sref
               end if
c
            end if
c
         end if
c
c        control surface data (global level only)
c
         if (level.eq.lglobal) then
            if (myid.eq.myhost) then
c
c              zero out values if this is the first block on the global level
               if (nbl.eq.nblstag) then
                  sxsub     = 0.e0
                  sysub     = 0.e0
                  szsub     = 0.e0
                  stotsub   = 0.e0
                  fmdotsub  = 0.e0
                  cfxpsub   = 0.e0
                  cfypsub   = 0.e0
                  cfzpsub   = 0.e0
                  cflpsub   = 0.e0
                  cfdpsub   = 0.e0
                  cfxvsub   = 0.e0
                  cfyvsub   = 0.e0
                  cfzvsub   = 0.e0
                  cflvsub   = 0.e0
                  cfdvsub   = 0.e0
                  cfxmomsub = 0.e0
                  cfymomsub = 0.e0
                  cfzmomsub = 0.e0
                  cflmomsub = 0.e0
                  cfdmomsub = 0.e0
                  cfxtotsub = 0.e0
                  cfytotsub = 0.e0
                  cfztotsub = 0.e0
                  cfdtotsub = 0.e0
                  cfltotsub = 0.e0
               end if
            end if
c
            do ics = 1,ncs
               fnorm = float(icsinfo(ics,9))
               if (nbl.eq.icsinfo(ics,1) .and.
     .            real(fnorm).ne.0.e0) then
                  iuns = max(itrans(nbl),irotat(nbl),idefrm(nbl))
                  if (myid.eq.mblk2nd(nbl)) then
                     call csurf(jdim,kdim,idim,x,y,z,sk,sj,si,q,ics,
     .                          q(1,2),q(1,3),q(1,4),vmuk,vmuj,vmui,vol,
     .                          bcj,bck,bci,blank,xtbj,xtbk,xtbi,
     .                          iuns,ncs,icsinfo,sxi,syi,szi,stoti,
     .                          pavi,ptavi,tavi,ttavi,xmavi,fmdoti,
     .                          cfxpi,cfypi,cfzpi,cfdpi,cflpi,cftpi,
     .                          cfxvi,cfyvi,cfzvi,cfdvi,cflvi,cftvi,
     .                          cfxmomi,cfymomi,cfzmomi,cfdmomi,cflmomi,
     .                          cftmomi,cfxtoti,cfytoti,cfztoti,cfdtoti,
     .                          cfltoti,cfttoti,maxcs,qj0,qk0,qi0)
                  end if
#if defined DIST_MPI
c
                  nd_srce = mblk2nd(nbl)
                  mytag = itag_fcoef + nbl
c
                  if (myid.eq.mblk2nd(nbl)) then
                     fcoef(1)  = sxi
                     fcoef(2)  = syi
                     fcoef(3)  = szi
                     fcoef(4)  = stoti
                     fcoef(5)  = fmdoti
                     fcoef(6)  = cfxpi
                     fcoef(7)  = cfypi
                     fcoef(8)  = cfzpi
                     fcoef(9)  = cflpi
                     fcoef(10) = cfdpi
                     fcoef(11) = cfxvi
                     fcoef(12) = cfyvi
                     fcoef(13) = cfzvi
                     fcoef(14) = cflvi
                     fcoef(15) = cfdvi
                     fcoef(16) = cfxmomi
                     fcoef(17) = cfymomi
                     fcoef(18) = cfzmomi
                     fcoef(19) = cflmomi
                     fcoef(20) = cfdmomi
                     fcoef(21) = cfxtoti
                     fcoef(22) = cfytoti
                     fcoef(23) = cfztoti
                     fcoef(24) = cfdtoti
                     fcoef(25) = cfltoti
                     call MPI_Send (fcoef, 25, MY_MPI_REAL,
     .                              myhost, mytag, mycomm, ierr)
                  end if
                  if (myid.eq.myhost) then
                     call MPI_Recv (fcoef, 25, MY_MPI_REAL,
     .                              nd_srce, mytag, mycomm, istat, ierr)
                     sxi     = fcoef(1)
                     syi     = fcoef(2)
                     szi     = fcoef(3)
                     stoti   = fcoef(4)
                     fmdoti  = fcoef(5)
                     cfxpi   = fcoef(6)
                     cfypi   = fcoef(7)
                     cfzpi   = fcoef(8)
                     cflpi   = fcoef(9)
                     cfdpi   = fcoef(10)
                     cfxvi   = fcoef(11)
                     cfyvi   = fcoef(12)
                     cfzvi   = fcoef(13)
                     cflvi   = fcoef(14)
                     cfdvi   = fcoef(15)
                     cfxmomi = fcoef(16)
                     cfymomi = fcoef(17)
                     cfzmomi = fcoef(18)
                     cflmomi = fcoef(19)
                     cfdmomi = fcoef(20)
                     cfxtoti = fcoef(21)
                     cfytoti = fcoef(22)
                     cfztoti = fcoef(23)
                     cfdtoti = fcoef(24)
                     cfltoti = fcoef(25)
                  end if
c
#endif
                  if (myid.eq.myhost) then
                     sxsub     = sxsub     + sxi
                     sysub     = sysub     + syi
                     szsub     = szsub     + szi
                     stotsub   = stotsub   + stoti
                     fmdotsub  = fmdotsub  + fmdoti
                     cfxpsub   = cfxpsub   + cfxpi
                     cfypsub   = cfypsub   + cfypi
                     cfzpsub   = cfzpsub   + cfzpi
                     cflpsub   = cflpsub   + cflpi
                     cfdpsub   = cfdpsub   + cfdpi
                     cfxvsub   = cfxvsub   + cfxvi
                     cfyvsub   = cfyvsub   + cfyvi
                     cfzvsub   = cfzvsub   + cfzvi
                     cflvsub   = cflvsub   + cflvi
                     cfdvsub   = cfdvsub   + cfdvi
                     cfxmomsub = cfxmomsub + cfxmomi
                     cfymomsub = cfymomsub + cfymomi
                     cfzmomsub = cfzmomsub + cfzmomi
                     cflmomsub = cflmomsub + cflmomi
                     cfdmomsub = cfdmomsub + cfdmomi
                     cfxtotsub = cfxtotsub + cfxtoti
                     cfytotsub = cfytotsub + cfytoti
                     cfztotsub = cfztotsub + cfztoti
                     cfdtotsub = cfdtotsub + cfdtoti
                     cfltotsub = cfltotsub + cfltoti
c
                     cftpsub   = sqrt(cfxpsub**2 
     .                              + cfypsub**2 
     .                              + cfzpsub**2)
                     cftvsub   = sqrt(cfxvsub**2 
     .                              + cfyvsub**2 
     .                              + cfzvsub**2)
                     cftmomsub = sqrt(cfxmomsub**2 
     .                              + cfymomsub**2 
     .                              + cfzmomsub**2)
                     cfttotsub = sqrt(cfxtotsub**2 
     .                              + cfytotsub**2
     .                              + cfztotsub**2)
                  end if
c
               end if
c
            end do
c
         end if
c
c
c        output subiteration convergence if iwrit3 > 0
c        iwrit3 > 0...density residual and forces (or cont. surf. data) only
c        iwrit3 > 1...density residual and forces (or cont. surf. data) AND 
c                     turbulence residual(s) (field eqn. turb. models only)
c
         iwrit3 = 2
c
         if (myid.eq.myhost) then
            if (ncyc.gt.1 .and. nbl .eq. nblendg) then
               nout = (nt-1)*ncyc+icyc
               if (iwrit3.gt.0) then
                  if (ihstry.eq.0) then
                     if (ialph.eq.0) then
                     if (nout.eq.1) write(23,1000)
                     write(23,1020) nout,log10(real(rmsb)),real(cltsub),
     .               real(cdtsub),real(cytsub),real(cmytsub)
                     else
                     if (nout.eq.1) write(23,1001)
                     write(23,1020) nout,log10(real(rmsb)),real(cltsub),
     .               real(cdtsub),real(cztsub),real(cmztsub)
                     end if
                     if(isubit_r.ne.0) then
                        write(111,'(I6,1x,5(E14.5))') nout,
     $                       log10(real(rmsba(1:5)))
                     end if
                  else if (ihstry.eq.1) then
                     if (nout.eq.1) write(23,1010)
                     write(23,1020) nout,log10(real(rmsb)),
     .               real(fmdotsub),real(cftpsub),real(cftvsub),
     .               real(cftmomsub)
                  else
                     if (nout.eq.1) write(23,1002)
                     write(23,1020) nout,log10(real(rmsb)),real(cltsub),
     .               real(cdtsub),real(cxtsub),real(cytsub),
     .               real(cztsub),real(cmxtsub),real(cmytsub),
     .               real(cmztsub)
                     if(isubit_r.ne.0) then
                        write(111,'(I6,1x,5(E14.5))') nout,
     $                       log10(real(rmsba(1:5)))
                     end if
                 end if
               end if
               if (iwrit3.gt.1) then
                  if (iviscg(nbl,1).gt.2 .or. iviscg(nbl,2).gt.2
     .            .or. iviscg(nbl,3).gt.2) then
                     if (nummem .eq. 2) then
                       if (nout.eq.1) write(24,1040)
                       write(24,1060) nout,log10(real(rmstb(1))),
     .                 log10(real(rmstb(2))),nnegb(1),nnegb(2)
                     else if (nummem .eq. 3) then
                       if (nout.eq.1) write(24,1041)
                       write(24,1061) nout,log10(real(rmstb(1))),
     .                 log10(real(rmstb(2))),log10(real(rmstb(3))),
     .                 nnegb(1),nnegb(2),nnegb(3)
                     else
                       if (nout.eq.1) write(24,1042)
                       write(24,1062) nout,log10(real(rmstb(1))),
     .                 log10(real(rmstb(2))),log10(real(rmstb(3))),
     .                 log10(real(rmstb(4))),log10(real(rmstb(5))),
     .                 log10(real(rmstb(6))),log10(real(rmstb(7))),
     .                 nnegb(1),nnegb(2),nnegb(3),nnegb(4),
     .                 nnegb(5),nnegb(6),nnegb(7)
                     end if
                  end if
               end if
            end if
         end if
c
      end if
c
 1000 format(4x,5hsubit,3x,11hlog(subres),7x,2hcl,12x,2hcd,
     .       12x,2hcy,11x,3hcmy)
 1001 format(4x,5hsubit,3x,11hlog(subres),7x,2hcl,12x,2hcd,
     .       12x,2hcz,11x,3hcmz)
 1002 format(4x,5hsubit,3x,11hlog(subres),7x,2hcl,12x,2hcd,
     .       12x,2hcx,12x,2hcy,12x,2hcz,11x,3hcmx,11x,3hcmy,11x,3hcmz)
 1010 format(4x,5hsubit,3x,11hlog(subres),4x,9hmass_flow,8x,4hcftp,
     .       10x,4hcftv,8x,6hcftmom)
 1020 format(3x,i6,9e14.5)
 1040 format(4x,5hsubit,2x,12hlog(turres1),2x,12hlog(turres2),
     .       2x,5hnneg1,2x,5hnneg2)
 1041 format(4x,5hsubit,2x,12hlog(turres1),2x,12hlog(turres2),
     .       2x,12hlog(turres3),2x,5hnneg1,2x,5hnneg2,2x,5hnneg3)
 1042 format(4x,5hsubit,2x,12hlog(turres1),2x,12hlog(turres2),
     .       2x,12hlog(turres3),2x,12hlog(turres4),2x,12hlog(turres5),
     .       2x,12hlog(turres6),2x,12hlog(turres7),2x,5hnneg1,
     .       2x,5hnneg2,2x,5hnneg3,2x,5hnneg4,2x,5hnneg5,2x,5hnneg6,
     .       2x,5hnneg7)
 1060 format(3x,i6,2e14.5,1x,i6,1x,i6)
 1061 format(3x,i6,3e14.5,1x,i6,1x,i6,1x,i6)
 1062 format(3x,i6,7e14.5,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)
 1080 format(/1x,5hlevel,4x,3hblk,3x,4hiter,3x,8hresidual,3x,
     .       10htotal res.,5x,4hlift,8x,4hdrag,5x,10hside force)
 2000 format(' ',4x,i1,1x,i6,1x,i6,7(1x,e11.4))
 2020 format(/1x,5hlevel,4x,3hblk,2x,5hntime,3x,8hresidual,3x,
     .       10htotal res.,5x,4hlift,8x,4hdrag,5x,10hside force)
 2040 format(/1x,5hlevel,4x,3hblk,2x,5hntime,3x,8hresidual,3x,
     .       10htotal res.,5x,4hlift,8x,4hdrag,5x,10hside force,
     .       5x,4htime,7x,5halpha)
 2060 format(/1x,23hmax residual on block  ,i6,7h (grid ,i6,
     .       5h) is ,e12.5,12h at j,k,i = ,i4,1h,,i4,1h,,i4,1h.)
 2080 format(1x,23hmax residual on block  ,i6,7h (grid ,i6,
     .       5h) is ,e12.5,12h at j,k,i = ,i4,1h,,i4,1h,,i4,1h.)
 3000 format(/1x,23hmax vorticity on block ,i6,7h (grid ,i6,
     .       5h) is ,e12.5,12h at j,k,i = ,i4,1h,,i4,1h,,i4,1h.)
 3020 format(1x,23hmax vorticity on block ,i6,7h (grid ,i6,
     .       5h) is ,e12.5,12h at j,k,i = ,i4,1h,,i4,1h,,i4,1h.)
c
      if (myid.eq.myhost) call my_flush(11)
      if (real(dt).ge.0.) then
         if (myid.eq.myhost) call my_flush(23)
         if (myid.eq.myhost) call my_flush(24)
      end if
c
      return
      end
