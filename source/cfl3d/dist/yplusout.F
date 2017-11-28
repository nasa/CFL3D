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
      subroutine yplusout(iseq,lw,lw2,w,mgwk,wk,nwork,cl,maxbl,maxgr,
     .                    maxseg,nblock,lwdat,levelg,igridg,jdimg,
     .                    kdimg,idimg,nbci0,nbcj0,nbck0,nbcidim,
     .                    nbcjdim,nbckdim,bcfilei,bcfilej,bcfilek,
     .                    itrans,irotat,idefrm,nblg,ibcinfo,jbcinfo,
     .                    kbcinfo,iadvance,iovrlp,myid,myhost,
     .                    mycomm,mblk2nd,idf,jdf,kdf,nou,bou,nbuf,
     .                    ibufdim,bcfiles,mxbcfil,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Driver routine for calculation of yplus at first gridpoint 
c               above solid walls 
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80  bcfiles(mxbcfil)
c
      integer  bcfilei,bcfilej,bcfilek
c
      dimension nou(nbuf)
      dimension w(mgwk),lwdat(maxbl,maxseg,6),lw(65,maxbl),wk(nwork)
      dimension idf(maxbl*6),jdf(maxbl*6),kdf(maxbl*6),lw2(43,maxbl)
      dimension iadvance(maxbl),iovrlp(maxbl)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension mblk2nd(maxbl),nblg(maxgr),levelg(maxbl)
      dimension bcfilei(maxbl,maxseg,2),bcfilej(maxbl,maxseg,2),
     .          bcfilek(maxbl,maxseg,2),igridg(maxbl),itrans(maxbl),
     .          irotat(maxbl),idefrm(maxbl),idimg(maxbl),jdimg(maxbl),
     .          kdimg(maxbl)
c
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /maxiv/ ivmx
      common /twod/ i2d
      common /igrdtyp/ ip3dgrd,ialph
      common /reyue/ reue,tinf,ivisc(3)
      common /ypinfo/ ypsumb,ypsumsqb,ypmaxb,ypminb,dnmaxb,dnminb,ypchk,
     .                nptsb,jypmaxb,kypmaxb,iypmaxb,jypminb,kypminb,
     .                iypminb,jdnmaxb,kdnmaxb,idnmaxb,jdnminb,
     .                kdnminb,idnminb,nypchkb
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
c
c     check for number of points with yplus > ypchk
c
      ypchk = 5.
c
c     ifunc - flag to generate plot3d function files (for 3d cases) with
c             y+, eddy viscosity and normal spacing at first gridpoint. For 2d
c             cases, a single file is output in column format, with
c             headers between different segments.
c             above all viscous walls
c           = 0 do not output function files 
c           > 0 output function files
c
      ifunc = 0
      iunit  = 28
c
      ypsumt   = 0.
      ypsumsqt = 0.
      ypmaxt   = 0.
      iypmaxt  = 0
      jypmaxt  = 0
      ypmint   = 1.e9
      iypmint  = 0
      jypmint  = 0
      dnmaxt   = 0.
      idnmaxt  = 0
      jdnmaxt  = 0
      dnmint   = 1.e9
      idnmint  = 0
      jdnmint  = 0
      nypchkt  = 0
      nptst    = 0
      nbypmaxt = 0
      nbypmint = 0
      nbdnmaxt = 0
      nbdnmint = 0
      iflag    = 0
c
      if (myid.eq.myhost) then
c
      if (ifunc.gt.0) then
c
c        - the surface grid file is written to iunit
c        - the associated y+ file is written to unit+1
c        - the associated dn (normal spacing) file is written to iunit+2
c        - the associated turb. viscosity file is written to iunit+3
c          (fast only alows one scalar variable per function file)
c        - in 2d, only one file is written, to iunit, containing
c          x,z (y if ialph>0),y+,dn and turb. viscosity, in column format.
c      
c
         if (i2d.eq.0) then
            open(unit=iunit,file='surf_xyz.fmt',form='formatted',
     .      status='unknown')
            open(unit=iunit+1,file='surf_y+.fmt',form='formatted',
     .      status='unknown')
            open(unit=iunit+2,file='surf_dn.fmt',form='formatted',
     .      status='unknown')
            open(unit=iunit+3,file='surf_vist.fmt',form='formatted',
     .      status='unknown')
         else
            open(unit=iunit,file='surf_y+_2d.fmt',form='formatted',
     .      status='unknown')
         end if
         nvar   = 1
c
c        count number of surfaces to be output to function file
c
         ncount = 0
         do 600 nbl=1,nblock
         call lead(nbl,lw,lw2,maxbl)
c
c        determine if any surfaces should be output: block should
c        be turbulent *and* have at least one solid surface
c
         ivv1 = 0
         if (ivisc(3).gt.1 .or. ivisc(2).gt.1 .or. ivisc(1).gt.1) ivv1=1
c
         ivv2 = 0
         do ns=1,nbci0(nbl)
            if(abs(ibcinfo(nbl,ns,1,1)) .eq. 2004 .or.
     .         abs(ibcinfo(nbl,ns,1,1)) .eq. 2014 .or.
     .         abs(ibcinfo(nbl,ns,1,1)) .eq. 2024 .or.
     .         abs(ibcinfo(nbl,ns,1,1)) .eq. 2034 .or.
     .         abs(ibcinfo(nbl,ns,1,1)) .eq. 2016) ivv2 = ivv2 + 1
         enddo
         do ns=1,nbcidim(nbl)
            if(abs(ibcinfo(nbl,ns,1,2)) .eq. 2004 .or.
     .         abs(ibcinfo(nbl,ns,1,2)) .eq. 2014 .or.
     .         abs(ibcinfo(nbl,ns,1,2)) .eq. 2024 .or.
     .         abs(ibcinfo(nbl,ns,1,2)) .eq. 2034 .or.
     .         abs(ibcinfo(nbl,ns,1,2)) .eq. 2016) ivv2 = ivv2 + 1
         enddo
         do ns=1,nbcj0(nbl)
            if(abs(jbcinfo(nbl,ns,1,1)) .eq. 2004 .or.
     .         abs(jbcinfo(nbl,ns,1,1)) .eq. 2014 .or.
     .         abs(jbcinfo(nbl,ns,1,1)) .eq. 2024 .or.
     .         abs(jbcinfo(nbl,ns,1,1)) .eq. 2034 .or.
     .         abs(jbcinfo(nbl,ns,1,1)) .eq. 2016) ivv2 = ivv2 + 1
         enddo
         do ns=1,nbcjdim(nbl)
            if(abs(jbcinfo(nbl,ns,1,2)) .eq. 2004 .or.
     .         abs(jbcinfo(nbl,ns,1,2)) .eq. 2014 .or.
     .         abs(jbcinfo(nbl,ns,1,2)) .eq. 2024 .or.
     .         abs(jbcinfo(nbl,ns,1,2)) .eq. 2034 .or.
     .         abs(jbcinfo(nbl,ns,1,2)) .eq. 2016) ivv2 = ivv2 + 1
         enddo
         do ns=1,nbck0(nbl)
            if(abs(kbcinfo(nbl,ns,1,1)) .eq. 2004 .or.
     .         abs(kbcinfo(nbl,ns,1,1)) .eq. 2014 .or.
     .         abs(kbcinfo(nbl,ns,1,1)) .eq. 2024 .or.
     .         abs(kbcinfo(nbl,ns,1,1)) .eq. 2034 .or.
     .         abs(kbcinfo(nbl,ns,1,1)) .eq. 2016) ivv2 = ivv2 + 1
         enddo
         do ns=1,nbckdim(nbl)
            if(abs(kbcinfo(nbl,ns,1,2)) .eq. 2004 .or.
     .         abs(kbcinfo(nbl,ns,1,2)) .eq. 2014 .or.
     .         abs(kbcinfo(nbl,ns,1,2)) .eq. 2024 .or.
     .         abs(kbcinfo(nbl,ns,1,2)) .eq. 2034 .or.
     .         abs(kbcinfo(nbl,ns,1,2)) .eq. 2016) ivv2 = ivv2 + 1
         enddo
c
         ivv = 0
         if (ivv1.gt.0 .and. ivv2.gt.0) ivv = 1
c
         levg = levelg(nbl)
c
         if (levg.ge.lglobal .and. levg.le.levt
     .      .and. ncyc1(iseq).gt.0 .and. ivv.gt.0) then
c
            if (ivisc(3).gt.1) then
               do 500 m=1,2
               if (m.eq.1) then
                  nbc = nbck0(nbl)
               else
                  nbc = nbckdim(nbl)
               end if
               ktest = 0
               do 510 nseg = 1,nbc
               if (abs(kbcinfo(nbl,nseg,1,m)) .eq. 2004 .or.
     .             abs(kbcinfo(nbl,nseg,1,m)) .eq. 2014 .or.
     .             abs(kbcinfo(nbl,nseg,1,m)) .eq. 2024 .or.
     .             abs(kbcinfo(nbl,nseg,1,m)) .eq. 2034 .or.
     .             abs(kbcinfo(nbl,nseg,1,m)) .eq. 2016) ktest = 1
 510           continue
               if (ktest.gt.0) then
                  ncount = ncount+1
                  idf(ncount) = idim
                  if (i2d.eq.1) idf(ncount) = 1 
                  jdf(ncount) = jdim
                  kdf(ncount) = 1
               end if
 500           continue
            end if
            if (ivisc(2).gt.1) then
               do 520 m=1,2
               if (m.eq.1) then
                  nbc = nbcj0(nbl)
               else
                  nbc = nbcjdim(nbl)
               end if
               jtest = 0
               do 530 nseg = 1,nbc
               if (abs(jbcinfo(nbl,nseg,1,m)) .eq. 2004 .or.
     .             abs(jbcinfo(nbl,nseg,1,m)) .eq. 2014 .or.
     .             abs(jbcinfo(nbl,nseg,1,m)) .eq. 2024 .or.
     .             abs(jbcinfo(nbl,nseg,1,m)) .eq. 2034 .or.
     .             abs(jbcinfo(nbl,nseg,1,m)) .eq. 2016) jtest = 1
 530           continue
               if (jtest.gt.0) then
                  ncount = ncount+1
                  idf(ncount) = idim
                  if (i2d.eq.1) idf(ncount) = 1
                  jdf(ncount) = 1
                  kdf(ncount) = kdim
               end if
 520           continue
            end if
            if (ivisc(1).gt.1) then
               do 540 m=1,2
               if (m.eq.1) then
                  nbc = nbci0(nbl)
               else
                  nbc = nbcidim(nbl)
               end if
               itest = 0
               do 550 nseg = 1,nbc
               if (abs(ibcinfo(nbl,nseg,1,m)) .eq. 2004 .or.
     .             abs(ibcinfo(nbl,nseg,1,m)) .eq. 2014 .or.
     .             abs(ibcinfo(nbl,nseg,1,m)) .eq. 2024 .or.
     .             abs(ibcinfo(nbl,nseg,1,m)) .eq. 2034 .or.
     .             abs(ibcinfo(nbl,nseg,1,m)) .eq. 2016) itest = 1
 550           continue
               if (itest.gt.0) then
                  ncount = ncount+1
                  idf(ncount) = 1
                  jdf(ncount) = jdim
                  kdf(ncount) = kdim
               end if
 540           continue
            end if
c
         end if
c
 600     continue 
c
         if (ncount.gt.0) then
            if (i2d.eq.0) then
               write(iunit,*) ncount
               write(iunit,*) (idf(l),jdf(l),kdf(l),l=1,ncount)
               write(iunit+1,*) ncount
               write(iunit+1,*) (idf(l),jdf(l),kdf(l),nvar,l=1,ncount)
               write(iunit+2,*) ncount
               write(iunit+2,*) (idf(l),jdf(l),kdf(l),nvar,l=1,ncount)
               write(iunit+3,*) ncount
               write(iunit+3,*) (idf(l),jdf(l),kdf(l),nvar,l=1,ncount)
            else
               write(iunit,*) 'title = "surface y+ data"'
               if (ialph.eq.0) then
                  write(iunit,*) 'variables = x, z, y+, dn, turb_visc'
               else
                  write(iunit,*) 'variables = x, y, y+, dn, turb_visc'
               end if
            end if
         end if
c
      end if
c
      end if
c
c     calculate y+ data
c
c
      do 10 nbl=1,nblock
c
      if (myid.eq.mblk2nd(nbl).or.myid.eq.myhost) then
c
      call lead(nbl,lw,lw2,maxbl)
c
c     determine if yplus statistics should be output: block should
c     be turbulent *and* have at least one solid surface
c
      ivv1 = 0
      if (ivisc(3).gt.1 .or. ivisc(2).gt.1 .or. ivisc(1).gt.1) ivv1=1 
c
      ivv2 = 0
      do ns=1,nbci0(nbl)
         if(abs(ibcinfo(nbl,ns,1,1)) .eq. 2004 .or.
     .      abs(ibcinfo(nbl,ns,1,1)) .eq. 2014 .or.
     .      abs(ibcinfo(nbl,ns,1,1)) .eq. 2024 .or.
     .      abs(ibcinfo(nbl,ns,1,1)) .eq. 2034 .or.
     .      abs(ibcinfo(nbl,ns,1,1)) .eq. 2016) ivv2 = ivv2 + 1
      enddo
      do ns=1,nbcidim(nbl)
         if(abs(ibcinfo(nbl,ns,1,2)) .eq. 2004 .or.
     .      abs(ibcinfo(nbl,ns,1,2)) .eq. 2014 .or.
     .      abs(ibcinfo(nbl,ns,1,2)) .eq. 2024 .or.
     .      abs(ibcinfo(nbl,ns,1,2)) .eq. 2034 .or.
     .      abs(ibcinfo(nbl,ns,1,2)) .eq. 2016) ivv2 = ivv2 + 1
      enddo
      do ns=1,nbcj0(nbl)
         if(abs(jbcinfo(nbl,ns,1,1)) .eq. 2004 .or.
     .      abs(jbcinfo(nbl,ns,1,1)) .eq. 2014 .or.
     .      abs(jbcinfo(nbl,ns,1,1)) .eq. 2024 .or.
     .      abs(jbcinfo(nbl,ns,1,1)) .eq. 2034 .or.
     .      abs(jbcinfo(nbl,ns,1,1)) .eq. 2016) ivv2 = ivv2 + 1
      enddo
      do ns=1,nbcjdim(nbl)
         if(abs(jbcinfo(nbl,ns,1,2)) .eq. 2004 .or.
     .      abs(jbcinfo(nbl,ns,1,2)) .eq. 2014 .or.
     .      abs(jbcinfo(nbl,ns,1,2)) .eq. 2024 .or.
     .      abs(jbcinfo(nbl,ns,1,2)) .eq. 2034 .or.
     .      abs(jbcinfo(nbl,ns,1,2)) .eq. 2016) ivv2 = ivv2 + 1
      enddo
      do ns=1,nbck0(nbl)
         if(abs(kbcinfo(nbl,ns,1,1)) .eq. 2004 .or.
     .      abs(kbcinfo(nbl,ns,1,1)) .eq. 2014 .or.
     .      abs(kbcinfo(nbl,ns,1,1)) .eq. 2024 .or.
     .      abs(kbcinfo(nbl,ns,1,1)) .eq. 2034 .or.
     .      abs(kbcinfo(nbl,ns,1,1)) .eq. 2016) ivv2 = ivv2 + 1
      enddo
      do ns=1,nbckdim(nbl)
         if(abs(kbcinfo(nbl,ns,1,2)) .eq. 2004 .or.
     .      abs(kbcinfo(nbl,ns,1,2)) .eq. 2014 .or.
     .      abs(kbcinfo(nbl,ns,1,2)) .eq. 2024 .or.
     .      abs(kbcinfo(nbl,ns,1,2)) .eq. 2034 .or.
     .      abs(kbcinfo(nbl,ns,1,2)) .eq. 2016) ivv2 = ivv2 + 1
      enddo
c
      ivv = 0
      if (ivv1.gt.0 .and. ivv2.gt.0) ivv = 1
c
      levg = levelg(nbl)
c
c     output statistics on global and embedded levels
c
      if (levg.ge.lglobal .and. levg.le.levt 
     .    .and. ncyc1(iseq).gt.0 .and. ivv.gt.0) then
c
         if (myid.eq.myhost) then
            write(11,100) nbl,igridg(nbl)
         end if
c
         if (mblk2nd(nbl).eq.myid) then
            if (iadvance(nbl).ge.0) then
               call bc(1,nbl,lw,lw2,w,mgwk,wk,nwork,cl,
     .                 nou,bou,nbuf,ibufdim,maxbl,maxgr,maxseg,itrans,
     .                 irotat,idefrm,igridg,nblg,nbci0,nbcj0,nbck0,
     .                 nbcidim,nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                 bcfilei,bcfilej,bcfilek,lwdat,myid,
     .                 idimg,jdimg,kdimg,bcfiles,mxbcfil,nummem)
            end if
            ldim = 5
            call qface(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .                 w(lbcj),w(lbck),w(lbci),w(lblk),ldim)
            if (ivmx .ge. 2) then
               ldim = 1
               call qface(jdim,kdim,idim,w(lvis),w(lvj0),w(lvk0),
     .                    w(lvi0),w(lbcj),w(lbck),w(lbci),w(lblk),ldim)
            end if
         end if
c
         lypj = 1
         lypk = lypj + kdim*idim*2
         lypi = lypk + jdim*idim*2
         lblj = lypi + jdim*kdim*2
         lblk = lblj + kdim*idim*2
         lbli = lblk + jdim*idim*2
         ldnj = lbli + jdim*kdim*2
         ldnk = ldnj + kdim*idim*2
         ldni = ldnk + jdim*idim*2
         lvtj = ldni + jdim*kdim*2
         lvtk = lvtj + kdim*idim*2
         lvti = lvtk + jdim*idim*2
         call calyplus(jdim,kdim,idim,nbl,w(lq),w(lqi0),w(lqj0),
     .                 w(lqk0),w(lx),w(ly),w(lz),w(lvis),iovrlp(nbl),
     .                 w(lbcj),w(lbck),w(lbci),w(lsj),w(lsk),w(lsi),
     .                 wk(lypj),wk(lypk),wk(lypi),wk(lblj),wk(lblk),
     .                 wk(lbli),wk(ldnj),wk(ldnk),wk(ldni),wk(lvtj),
     .                 wk(lvtk),wk(lvti),ifunc,iunit,
     .                 w(lvj0),w(lvk0),w(lvi0),maxbl,maxseg,nbci0,
     .                 nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,myid,
     .                 myhost,mycomm,mblk2nd,ibcinfo,jbcinfo,
     .                 kbcinfo,w(lvol))
c
c        save yplus statistcs on on global level for final summary
c
         if (levg.eq.lglobal) then
            iflag = 1
            nptst    = nptst    + nptsb
            ypsumt   = ypsumt   + ypsumb
            ypsumsqt = ypsumsqt + ypsumsqb
            nypchkt  = nypchkt  + nypchkb
            if (real(ypmaxb) .gt. real(ypmaxt)) then
               ypmaxt   = ypmaxb
               iypmaxt  = iypmaxb
               jypmaxt  = jypmaxb
               kypmaxt  = kypmaxb
               nbypmaxt = nbl
            end if
            if (real(ypminb) .lt. real(ypmint)) then
               ypmint   = ypminb
               iypmint = iypminb
               jypmint  = jypminb
               kypmint  = kypminb
               nbypmint = nbl
            end if
            if (real(dnmaxb) .gt. real(dnmaxt)) then
               dnmaxt   = dnmaxb
               idnmaxt  = idnmaxb
               jdnmaxt  = jdnmaxb
               kdnmaxt  = kdnmaxb
               nbdnmaxt = nbl
            end if
            if (real(dnminb) .lt. real(dnmint)) then
               dnmint   = dnminb
               idnmint  = idnminb
               jdnmint  = jdnminb
               kdnmint  = kdnminb
               nbdnmint = nbl
            end if
         end if
c
      end if
c
      end if
c
   10 continue
c
      if (myid.eq.myhost) then
         if (iflag.gt.0) then
            ypavt = ypsumt/nptst
            ypstdt = sqrt(ypsumsqt/(nptst-1))
            write(11,101)
            write(11,102)
            write(11,107) real(ypmaxt),iypmaxt,jypmaxt,kypmaxt,
     .                    nbypmaxt,igridg(nbypmaxt)
            write(11,103)
            write(11,107) real(ypmint),iypmint,jypmint,kypmint,
     .                    nbypmint,igridg(nbypmint)
            write(11,104)
            write(11,107) real(dnmaxt),idnmaxt,jdnmaxt,kdnmaxt,
     .                    nbdnmaxt,igridg(nbdnmaxt)
            write(11,105) 
            write(11,107) real(dnmint),idnmint,jdnmint,kdnmint,
     .                    nbdnmint,igridg(nbdnmint)
            write(11,106) int(real(ypchk))
            write(11,108) real(ypavt),real(ypstdt),nypchkt,nptst 
            if (real(ypavt) .gt. 2.5) then
               write(11,114)
            end if
            if (iwarneddy .gt. 0) then
               write(11,113) edvislim
            end if
            if (ifunc.gt.0) then
               if (i2d.eq.0) write(11,109)
               if (i2d.eq.1) then
                  if (ialph.eq.0) then
                     write(11,111)
                  else
                     write(11,112)
                  end if
               end if
            else
               write(11,110)
            end if
         end if
      end if
c
  100 format(/,1x,46hYPLUS STATISTICS (endpts not included) - BLOCK,i6,
     .      6h (GRID,i6,1h))
  101 format(/,1x,38hYPLUS STATISTICS (endpts not included),
     .      20h - ALL GLOBAL BLOCKS)
  102 format(5x,6hY+ MAX,4x,4hILOC,4x,4hJLOC,4x,4hKLOC,
     .       5x,5hBLOCK,6x,4hGRID)
  103 format(5x,6hY+ MIN,4x,4hILOC,4x,4hJLOC,4x,4hKLOC,
     .       5x,5hBLOCK,6x,4hGRID)
  104 format(5x,6hDN MAX,4x,4hILOC,4x,4hJLOC,4x,4hKLOC,
     .       5x,5hBLOCK,6x,4hGRID)
  105 format(5x,6hDN MIN,4x,4hILOC,4x,4hJLOC,4x,4hKLOC,
     .       5x,5hBLOCK,6x,4hGRID)
  106 format(5x,6hY+ AVG,4x,10hY+ STD DEV,7x,5hNY+ >,i2,
     .       3x,4hNPTS)
  107 format(1x,e10.3,4x,i4,4x,i4,4x,i4,4x,i6,4x,i6)
  108 format(1x,e10.3,4x,e10.3,8x,i6,1x,i6)
  109 format(/,1x,'writing function file surf_y+.fmt (y+)',
     .       /,1x,'writing function file surf_dn.fmt (normal spacing)',
     .       /,1x,'writing function file surf_vist.fmt (turb. visc.)',
     .       /,1x,'surf_xyz.fmt is the corresponding surface grid',
     .       /,1x,'(use multi/formatted/blank as FAST read options)')
  110 format(/,1x,'set ifunc = 1 in subroutine yplusout and recompile',
     .       ' if function files of',
     .       /15x,'y+, dn, and turb visc are desired')
  111 format(/,1x,'writing 2d data file surf_y+_2d.fmt',
     .       ' (x,z,y+,dn,turb. visc.)')
  112 format(/,1x,'writing 2d data file surf_y+_2d.fmt',
     .       ' (x,y,y+,dn,turb. visc.)')
  113 format(/,1x,'WARNING: final value of eddy viscosity has been',
     .         1x,'limited at one or more',
     .       /,1x,'points to ',e10.3,'.  Unless this run is far from',
     .         1x,'convergence, increase',
     .       /,1x,'edvislim via keyword input for subsequent runs')
  114 format(/,1x,'WARNING: avg y+ is > 2.5.  It is recommended that',
     .         1x,'you revise the grid',
     .       /,1x,'to have smaller min spacing at walls.  (If you',
     .         1x,'are employing WALL FUNCTIONS,',
     .       /,1x,'large avg y+ values are acceptable, but wall',
     .         1x,'functions are ad hoc',
     .       /,1x,'and not recommended for general use.)')
c
      return
      end
