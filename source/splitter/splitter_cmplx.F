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
      subroutine splitter_cmplx(mbloc,maxseg,npmax,mxbli,intmax,nsub1,
     .                          msegt,msegn,msplt,mtot,nou,bou,nbuf,
     .                          ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose: split grid and input files
c***********************************************************************     
c
      character*120 bou(ibufdim,nbuf)
c
      integer xif1,xif2,xit1,xit2,etf1,etf2,ett1,ett2
      integer stats
c
      dimension nou(nbuf)
c
      allocatable :: bcval(:,:,:,:)
      allocatable :: bcvali(:,:,:,:)
      allocatable :: bcvalj(:,:,:,:)
      allocatable :: bcvalk(:,:,:,:)
      allocatable :: cpchamber(:,:)
      allocatable :: cporous(:,:)
      allocatable :: cq(:,:)
      allocatable :: cradiation(:,:)
      allocatable :: etf1(:,:)
      allocatable :: etf2(:,:)
      allocatable :: ett1(:)
      allocatable :: ett2(:)
      allocatable :: iadvance(:)
      allocatable :: ibcinfo(:,:,:,:)
      allocatable :: ibct(:,:,:)
      allocatable :: ibeg(:)
      allocatable :: ibif1(:)
      allocatable :: ibif2(:)
      allocatable :: iconcat(:)
      allocatable :: idbloc(:)
      allocatable :: idegg(:,:)
      allocatable :: idiagg(:,:)
      allocatable :: idimg(:)
      allocatable :: idmold(:)
      allocatable :: idnext(:)
      allocatable :: idno(:)
      allocatable :: idobl(:)
      allocatable :: idseg(:)
      allocatable :: ieg(:)
      allocatable :: iemg(:)
      allocatable :: ifdsg(:,:)
      allocatable :: iflimg(:,:)
      allocatable :: ifoflg(:,:,:)
      allocatable :: iforce(:)
      allocatable :: ifrom(:,:)
      allocatable :: ifsor(:)
      allocatable :: igridg(:)
      allocatable :: iic0(:)
      allocatable :: iifit(:)
      allocatable :: iiorph(:)
      allocatable :: iitmax(:)
      allocatable :: ijk(:,:)
      allocatable :: il(:)
      allocatable :: ilamhig(:)
      allocatable :: ilamlog(:)
      allocatable :: iln(:)
      allocatable :: ilo(:)
      allocatable :: ilosd(:)
      allocatable :: imap(:,:,:)
      allocatable :: imax(:)
      allocatable :: imin(:)
      allocatable :: inewgg(:)
      allocatable :: intrfc(:)
      allocatable :: iold(:,:)
      allocatable :: iovrlp(:)
      allocatable :: ipatch(:)
      allocatable :: iredundant(:,:)
      allocatable :: isblocin(:)
      allocatable :: isg(:)
      allocatable :: isin(:)
      allocatable :: isin0(:)
      allocatable :: iskp(:)
      allocatable :: isva(:,:,:)
      allocatable :: ito(:)
      allocatable :: itrb(:,:)
      allocatable :: itrb1(:)
      allocatable :: itrb2(:)
      allocatable :: iturbb(:)
      allocatable :: ivisb(:,:)
      allocatable :: iviscg(:,:)
      allocatable :: jbcinfo(:,:,:,:)
      allocatable :: jdimg(:)
      allocatable :: jeg(:)
      allocatable :: jl(:)
      allocatable :: jlamhig(:)
      allocatable :: jlamlog(:)
      allocatable :: jln(:)
      allocatable :: jlo(:)
      allocatable :: jlosd(:)
      allocatable :: jmax(:)
      allocatable :: jmin(:)
      allocatable :: jsg(:)
      allocatable :: jskp(:)
      allocatable :: jtrb1(:)
      allocatable :: jtrb2(:)
      allocatable :: kbcinfo(:,:,:,:)
      allocatable :: kdimg(:)
      allocatable :: keg(:)
      allocatable :: kl(:)
      allocatable :: klamhig(:)
      allocatable :: klamlog(:)
      allocatable :: kln(:)
      allocatable :: klo(:)
      allocatable :: klosd(:)
      allocatable :: kmax(:)
      allocatable :: kmin(:)
      allocatable :: ksg(:)
      allocatable :: kskp(:)
      allocatable :: ktrb1(:)
      allocatable :: ktrb2(:)
      allocatable :: levelg(:)
      allocatable :: limblk(:,:,:)
      allocatable :: llimit(:)
      allocatable :: mit(:,:)
      allocatable :: mmceta(:)
      allocatable :: mmcxie(:)
      allocatable :: nb1(:,:,:)
      allocatable :: nb1s(:)
      allocatable :: nb2(:,:,:)
      allocatable :: nb2s(:)
      allocatable :: nbci0(:)
      allocatable :: nbcidim(:)
      allocatable :: nbcj0(:)
      allocatable :: nbcjdim(:)
      allocatable :: nbck0(:)
      allocatable :: nbckdim(:)
      allocatable :: nbi1(:)
      allocatable :: nbi2(:)
      allocatable :: nbj1(:)
      allocatable :: nbj2(:)
      allocatable :: nbk1(:)
      allocatable :: nbk2(:)
      allocatable :: nblcg(:)
      allocatable :: nblg(:)
      allocatable :: nblkk(:,:)
      allocatable :: nblon(:)
      allocatable :: nbs(:)
      allocatable :: ncgg(:)
      allocatable :: nd11(:)
      allocatable :: nd12(:)
      allocatable :: nd21(:)
      allocatable :: nd22(:)
      allocatable :: ndat(:,:,:)
      allocatable :: ndirin(:)
      allocatable :: ndirin0(:)
      allocatable :: ne1(:,:,:)
      allocatable :: ne1s(:)
      allocatable :: ne2(:,:,:)
      allocatable :: ne2s(:)
      allocatable :: nei1(:)
      allocatable :: nei2(:)
      allocatable :: nej1(:)
      allocatable :: nej2(:)
      allocatable :: nek1(:)
      allocatable :: nek2(:)
      allocatable :: nfb(:)
      allocatable :: ni(:)
      allocatable :: nj(:)
      allocatable :: nk(:)
      allocatable :: npts(:)
      allocatable :: nseg(:)
      allocatable :: nseg1(:)
      allocatable :: nsihi(:)
      allocatable :: nsilo(:)
      allocatable :: nsjhi(:)
      allocatable :: nsjlo(:)
      allocatable :: nskhi(:)
      allocatable :: nsklo(:)
      allocatable :: nsubbl(:)
      allocatable :: nxtseg(:)
      allocatable :: rkap0g(:,:)
      allocatable :: twotref(:,:)
      allocatable :: val(:)
      complex, allocatable :: x(:)
      allocatable :: xdmold(:)
      allocatable :: xdum(:,:)
      allocatable :: xif1(:,:)
      allocatable :: xif2(:,:)
      allocatable :: xit1(:)
      allocatable :: xit2(:)
      allocatable :: xmap(:,:,:)
      complex, allocatable :: y(:)
      complex, allocatable :: z(:)
c
      character*80 gridin,gridout,header
      character*80 imapin,imapout
      character*80 cflinp,cflout,roninp,ronout,sdgridin,
     .             sdgridout,dovrlap,dpatch,dresid
c
      character*80 plt3dg,plt3dq,output,resid,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt
c
      common /unit5/ iunit5
      common /grdinfo/ nbloc
      common /cflfiles/gridin,plt3dg,plt3dq,output,resid,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 gridout,sdgridin,sdgridout,dovrlap,dpatch,dresid
      common /chk1/ ncgmax
      common /info3/ ncg,isnd,ialph,irest,iunst,ntstep,
     .               ita,ihstry,nplot3d,nprint,nwrest,ichk,i2d,
     .               mglev(5),nem(5),mitl(5,5),mtt,ndv,isdform,
     .               ip3dgrd,iplt3dtyp
c
c     allocate memory
c
      memuse = 0
c
      allocate( bcval(mbloc,6,msegn,7), stat=stats )
      call umalloc(mbloc*6*msegn*7,0,'bcval',memuse,stats)
      allocate( bcvali(mbloc,maxseg,7,2), stat=stats )
      call umalloc(mbloc*maxseg*7*2,0,'bcvali',memuse,stats)
      allocate( bcvalj(mbloc,maxseg,7,2), stat=stats )
      call umalloc(mbloc*maxseg*7*2,0,'bcvalj',memuse,stats)
      allocate( bcvalk(mbloc,maxseg,7,2), stat=stats )
      call umalloc(mbloc*maxseg*7*2,0,'bcvalk',memuse,stats)
      allocate( cpchamber(msegn,mbloc), stat=stats )
      call umalloc(msegn*mbloc,0,'cpchamber',memuse,stats)
      allocate( cporous(msegn,mbloc), stat=stats )
      call umalloc(msegn*mbloc,0,'cporous',memuse,stats)
      allocate( cq(msegn,mbloc), stat=stats )
      call umalloc(msegn*mbloc,0,'cq',memuse,stats)
      allocate( cradiation(msegn,mbloc), stat=stats )
      call umalloc(msegn*mbloc,0,'cradiation',memuse,stats)
      allocate( etf1(intmax,nsub1), stat=stats )
      call umalloc(intmax*nsub1,1,'etf1',memuse,stats)
      allocate( etf2(intmax,nsub1), stat=stats )
      call umalloc(intmax*nsub1,1,'etf2',memuse,stats)
      allocate( ett1(intmax), stat=stats )
      call umalloc(intmax,1,'ett1',memuse,stats)
      allocate( ett2(intmax), stat=stats )
      call umalloc(intmax,1,'ett2',memuse,stats)
      allocate( iadvance(mbloc), stat=stats )
      call umalloc(mbloc,1,'iadvance',memuse,stats)
      allocate( ibcinfo(mbloc,maxseg,7,2), stat=stats )
      call umalloc(mbloc*maxseg*7*2,1,'ibcinfo',memuse,stats)
      allocate( ibct(mbloc,6,msegn), stat=stats )
      call umalloc(mbloc*6*msegn,1,'ibct',memuse,stats)
      allocate( ibeg(mbloc), stat=stats )
      call umalloc(mbloc,1,'ibeg',memuse,stats)
      allocate( ibif1(mxbli), stat=stats )
      call umalloc(mxbli,1,'ibif1',memuse,stats)
      allocate( ibif2(mxbli), stat=stats )
      call umalloc(mxbli,1,'ibif2',memuse,stats)
      allocate( iconcat(intmax), stat=stats )
      call umalloc(intmax,1,'iconcat',memuse,stats)
      allocate( idbloc(mtot), stat=stats )
      call umalloc(mtot,1,'idbloc',memuse,stats)
      allocate( idegg(mbloc,3), stat=stats )
      call umalloc(mbloc*3,1,'idegg',memuse,stats)
      allocate( idiagg(mbloc,3), stat=stats )
      call umalloc(mbloc*3,1,'idiagg',memuse,stats)
      allocate( idimg(mbloc), stat=stats )
      call umalloc(mbloc,1,'idimg',memuse,stats)
      allocate( idmold(msegt), stat=stats )
      call umalloc(msegt,1,'idmold',memuse,stats)
      allocate( idnext(mbloc), stat=stats )
      call umalloc(mbloc,1,'idnext',memuse,stats)
      allocate( idno(mbloc), stat=stats )
      call umalloc(mbloc,1,'dno',memuse,stats)
      allocate( idobl(mbloc), stat=stats )
      call umalloc(mbloc,1,'idobl',memuse,stats)
      allocate( idseg(mbloc), stat=stats )
      call umalloc(mbloc,1,'idseg',memuse,stats)
      allocate( ieg(mbloc), stat=stats )
      call umalloc(mbloc,1,'ieg',memuse,stats)
      allocate( iemg(mbloc), stat=stats )
      call umalloc(mbloc,1,'iemg',memuse,stats)
      allocate( ifdsg(mbloc,3), stat=stats )
      call umalloc(mbloc*3,1,'ifdsg',memuse,stats)
      allocate( iflimg(mbloc,3), stat=stats )
      call umalloc(mbloc*3,1,'iflimg',memuse,stats)
      allocate( ifoflg(mbloc,6,msegn), stat=stats )
      call umalloc(mbloc*6*msegn,1,'ifoflg',memuse,stats)
      allocate( iforce(mbloc), stat=stats )
      call umalloc(mbloc,1,'iforce',memuse,stats)
      allocate( ifrom(intmax,nsub1), stat=stats )
      call umalloc(intmax*nsub1,1,'ifrom',memuse,stats)
      allocate( ifsor(mbloc), stat=stats )
      call umalloc(mbloc,1,'ifsor',memuse,stats)
      allocate( igridg(mbloc), stat=stats )
      call umalloc(mbloc,1,'igridg',memuse,stats)
      allocate( iic0(intmax), stat=stats )
      call umalloc(intmax,1,'iic0',memuse,stats)
      allocate( iifit(intmax), stat=stats )
      call umalloc(intmax,1,'iifit',memuse,stats)
      allocate( iiorph(intmax), stat=stats )
      call umalloc(intmax,1,'iiorph',memuse,stats)
      allocate( iitmax(intmax), stat=stats )
      call umalloc(intmax,1,'iitmax',memuse,stats)
      allocate( ijk(6,mbloc), stat=stats )
      call umalloc(6*mbloc,1,'ijk',memuse,stats)
      allocate( il(mbloc), stat=stats )
      call umalloc(mbloc,1,'il',memuse,stats)
      allocate( ilamhig(mbloc), stat=stats )
      call umalloc(mbloc,1,'ilamhig',memuse,stats)
      allocate( ilamlog(mbloc), stat=stats )
      call umalloc(mbloc,1,'ilamlog',memuse,stats)
      allocate( iln(mbloc), stat=stats )
      call umalloc(mbloc,1,'iln',memuse,stats)
      allocate( ilo(mbloc), stat=stats )
      call umalloc(mbloc,1,'ilo',memuse,stats)
      allocate( ilosd(mbloc), stat=stats )
      call umalloc(mbloc,1,'ilosd',memuse,stats)
      allocate( imap(msegt,msegn,mbloc), stat=stats )
      call umalloc(msegt*msegn*mbloc,1,'imap',memuse,stats)
      allocate( imax(mbloc), stat=stats )
      call umalloc(mbloc,1,'imax',memuse,stats)
      allocate( imin(mbloc), stat=stats )
      call umalloc(mbloc,1,'imin',memuse,stats)
      allocate( inewgg(mbloc), stat=stats )
      call umalloc(mbloc,1,'inewgg',memuse,stats)
      allocate( intrfc(mtot), stat=stats )
      call umalloc(mtot,1,'intrfc',memuse,stats)
      allocate( iold(4,mtot), stat=stats )
      call umalloc(4*mtot,1,'iold',memuse,stats)
      allocate( iovrlp(mbloc), stat=stats )
      call umalloc(mbloc,1,'iovrlp',memuse,stats)
      allocate( ipatch(mtot), stat=stats )
      call umalloc(mtot,1,'ipatch',memuse,stats)
      allocate( iredundant(mbloc,msegn), stat=stats )
      call umalloc(mbloc*msegn,1,'iredundant',memuse,stats)
      allocate( isblocin(msplt*mbloc), stat=stats )
      call umalloc(msplt*mbloc,1,'isblocin',memuse,stats)
      allocate( isg(mbloc), stat=stats )
      call umalloc(mbloc,1,'isg',memuse,stats)
      allocate( isin(msplt*mbloc), stat=stats )
      call umalloc(msplt*mbloc,1,'isin',memuse,stats)
      allocate( isin0(msplt), stat=stats )
      call umalloc(msplt,1,'isin0',memuse,stats)
      allocate( iskp(mbloc), stat=stats )
      call umalloc(mbloc,1,'iskp',memuse,stats)
      allocate( isva(2,2,mxbli), stat=stats )
      call umalloc(2*2*mxbli,1,'isva',memuse,stats)
      allocate( ito(intmax), stat=stats )
      call umalloc(intmax,1,'ito',memuse,stats)
      allocate( itrb(7,mbloc), stat=stats )
      call umalloc(7*mbloc,1,'itrb',memuse,stats)
      allocate( itrb1(mbloc), stat=stats )
      call umalloc(mbloc,1,'itrb1',memuse,stats)
      allocate( itrb2(mbloc), stat=stats )
      call umalloc(mbloc,1,'itrb2',memuse,stats)
      allocate( iturbb(mbloc), stat=stats )
      call umalloc(mbloc,1,'iturbb',memuse,stats)
      allocate( ivisb(msegt,mbloc), stat=stats )
      call umalloc(msegt*mbloc,1,'ivisb',memuse,stats)
      allocate( iviscg(mbloc,3), stat=stats )
      call umalloc(mbloc*3,1,'iviscg',memuse,stats)
      allocate( jbcinfo(mbloc,maxseg,7,2), stat=stats )
      call umalloc(mbloc*maxseg*7*2,1,'jbcinfo',memuse,stats)
      allocate( jdimg(mbloc), stat=stats )
      call umalloc(mbloc,1,'jdimg',memuse,stats)
      allocate( jeg(mbloc), stat=stats )
      call umalloc(mbloc,1,'jeg',memuse,stats)
      allocate( jl(mbloc), stat=stats )
      call umalloc(mbloc,1,'jl',memuse,stats)
      allocate( jlamhig(mbloc), stat=stats )
      call umalloc(mbloc,1,'jlamhig',memuse,stats)
      allocate( jlamlog(mbloc), stat=stats )
      call umalloc(mbloc,1,'jlamlog',memuse,stats)
      allocate( jln(mbloc), stat=stats )
      call umalloc(mbloc,1,'jln',memuse,stats)
      allocate( jlo(mbloc), stat=stats )
      call umalloc(mbloc,1,'jlo',memuse,stats)
      allocate( jlosd(mbloc), stat=stats )
      call umalloc(mbloc,1,'jlosd',memuse,stats)
      allocate( jmax(mbloc), stat=stats )
      call umalloc(mbloc,1,'jmax',memuse,stats)
      allocate( jmin(mbloc), stat=stats )
      call umalloc(mbloc,1,'jmin',memuse,stats)
      allocate( jsg(mbloc), stat=stats )
      call umalloc(mbloc,1,'jsg',memuse,stats)
      allocate( jskp(mbloc), stat=stats )
      call umalloc(mbloc,1,'jskp',memuse,stats)
      allocate( jtrb1(mbloc), stat=stats )
      call umalloc(mbloc,1,'jtrb1',memuse,stats)
      allocate( jtrb2(mbloc), stat=stats )
      call umalloc(mbloc,1,'jtrb2',memuse,stats)
      allocate( kbcinfo(mbloc,maxseg,7,2), stat=stats )
      call umalloc(mbloc*maxseg*7*2,1,'kbcinfo',memuse,stats)
      allocate( kdimg(mbloc), stat=stats )
      call umalloc(mbloc,1,'kdimg',memuse,stats)
      allocate( keg(mbloc), stat=stats )
      call umalloc(mbloc,1,'keg',memuse,stats)
      allocate( kl(mbloc), stat=stats )
      call umalloc(mbloc,1,'kl',memuse,stats)
      allocate( klamhig(mbloc), stat=stats )
      call umalloc(mbloc,1,'klamhig',memuse,stats)
      allocate( klamlog(mbloc), stat=stats )
      call umalloc(mbloc,1,'klamlog',memuse,stats)
      allocate( kln(mbloc), stat=stats )
      call umalloc(mbloc,1,'kln',memuse,stats)
      allocate( klo(mbloc), stat=stats )
      call umalloc(mbloc,1,'klo',memuse,stats)
      allocate( klosd(mbloc), stat=stats )
      call umalloc(mbloc,1,'klosd',memuse,stats)
      allocate( kmax(mbloc), stat=stats )
      call umalloc(mbloc,1,'kmax',memuse,stats)
      allocate( kmin(mbloc), stat=stats )
      call umalloc(mbloc,1,'kmin',memuse,stats)
      allocate( ksg(mbloc), stat=stats )
      call umalloc(mbloc,1,'ksg',memuse,stats)
      allocate( kskp(mbloc), stat=stats )
      call umalloc(mbloc,1,'kskp',memuse,stats)
      allocate( ktrb1(mbloc), stat=stats )
      call umalloc(mbloc,1,'ktrb1',memuse,stats)
      allocate( ktrb2(mbloc), stat=stats )
      call umalloc(mbloc,1,'ktrb2',memuse,stats)
      allocate( levelg(mbloc), stat=stats )
      call umalloc(mbloc,1,'levelg',memuse,stats)
      allocate( limblk(2,6,mxbli), stat=stats )
      call umalloc(2*6*mxbli,1,'limblk',memuse,stats)
      allocate( llimit(intmax), stat=stats )
      call umalloc(intmax,1,'llimit',memuse,stats)
      allocate( mit(5,mbloc), stat=stats )
      call umalloc(5*mbloc,1,'mit',memuse,stats)
      allocate( mmceta(intmax), stat=stats )
      call umalloc(intmax,1,'mmceta',memuse,stats)
      allocate( mmcxie(intmax), stat=stats )
      call umalloc(intmax,1,'mmcxie',memuse,stats)
      allocate( nb1(mbloc,6,msegn), stat=stats )
      call umalloc(mbloc*6*msegn,1,'nb1',memuse,stats)
      allocate( nb1s(mbloc), stat=stats )
      call umalloc(mbloc,1,'nb1s',memuse,stats)
      allocate( nb2(mbloc,6,msegn), stat=stats )
      call umalloc(mbloc*6*msegn,1,'nb2',memuse,stats)
      allocate( nb2s(mbloc), stat=stats )
      call umalloc(mbloc,1,'nb2s',memuse,stats)
      allocate( nbci0(mbloc), stat=stats )
      call umalloc(mbloc,1,'nbci0',memuse,stats)
      allocate( nbcidim(mbloc), stat=stats )
      call umalloc(mbloc,1,'nbcidim',memuse,stats)
      allocate( nbcj0(mbloc), stat=stats )
      call umalloc(mbloc,1,'nbcj0',memuse,stats)
      allocate( nbcjdim(mbloc), stat=stats )
      call umalloc(mbloc,1,'nbcjdim',memuse,stats)
      allocate( nbck0(mbloc), stat=stats )
      call umalloc(mbloc,1,'nbck0',memuse,stats)
      allocate( nbckdim(mbloc), stat=stats )
      call umalloc(mbloc,1,'nbckdim',memuse,stats)
      allocate( nbi1(mxbli), stat=stats )
      call umalloc(mxbli,1,'nbi1',memuse,stats)
      allocate( nbi2(mxbli), stat=stats )
      call umalloc(mxbli,1,'nbi2',memuse,stats)
      allocate( nbj1(mxbli), stat=stats )
      call umalloc(mxbli,1,'nbj1',memuse,stats)
      allocate( nbj2(mxbli), stat=stats )
      call umalloc(mxbli,1,'nbj2',memuse,stats)
      allocate( nbk1(mxbli), stat=stats )
      call umalloc(mxbli,1,'nbk1',memuse,stats)
      allocate( nbk2(mxbli), stat=stats )
      call umalloc(mxbli,1,'nbk2',memuse,stats)
      allocate( nblcg(mbloc), stat=stats )
      call umalloc(mbloc,1,'nblcg',memuse,stats)
      allocate( nblg(mbloc), stat=stats )
      call umalloc(mbloc,1,'nblg',memuse,stats)
      allocate( nblkk(2,mxbli), stat=stats )
      call umalloc(2*mxbli,1,'nblkk',memuse,stats)
      allocate( nblon(mxbli), stat=stats )
      call umalloc(mxbli,1,'nblon',memuse,stats)
      allocate( nbs(mbloc), stat=stats )
      call umalloc(mbloc,1,'nbs',memuse,stats)
      allocate( ncgg(mbloc), stat=stats )
      call umalloc(mbloc,1,'ncgg',memuse,stats)
      allocate( nd11(mxbli), stat=stats )
      call umalloc(mxbli,1,'nd11',memuse,stats)
      allocate( nd12(mxbli), stat=stats )
      call umalloc(mxbli,1,'nd12',memuse,stats)
      allocate( nd21(mxbli), stat=stats )
      call umalloc(mxbli,1,'nd21',memuse,stats)
      allocate( nd22(mxbli), stat=stats )
      call umalloc(mxbli,1,'nd22',memuse,stats)
      allocate( ndat(mbloc,6,msegn), stat=stats )
      call umalloc(mbloc*6*msegn,1,'ndat',memuse,stats)
      allocate( ndirin(msplt*mbloc), stat=stats )
      call umalloc(msplt*mbloc,1,'ndirin',memuse,stats)
      allocate( ndirin0(msplt), stat=stats )
      call umalloc(msplt,1,'ndirin0',memuse,stats)
      allocate( ne1(mbloc,6,msegn), stat=stats )
      call umalloc(mbloc*6*msegn,1,'ne1',memuse,stats)
      allocate( ne1s(mbloc), stat=stats )
      call umalloc(mbloc,1,'ne1s',memuse,stats)
      allocate( ne2(mbloc,6,msegn), stat=stats )
      call umalloc(mbloc*6*msegn,1,'ne2',memuse,stats)
      allocate( ne2s(mbloc), stat=stats )
      call umalloc(mbloc,1,'ne2s',memuse,stats)
      allocate( nei1(mxbli), stat=stats )
      call umalloc(mxbli,1,'nei1',memuse,stats)
      allocate( nei2(mxbli), stat=stats )
      call umalloc(mxbli,1,'nei2',memuse,stats)
      allocate( nej1(mxbli), stat=stats )
      call umalloc(mxbli,1,'nej1',memuse,stats)
      allocate( nej2(mxbli), stat=stats )
      call umalloc(mxbli,1,'nej2',memuse,stats)
      allocate( nek1(mxbli), stat=stats )
      call umalloc(mxbli,1,'nek1',memuse,stats)
      allocate( nek2(mxbli), stat=stats )
      call umalloc(mxbli,1,'nek2',memuse,stats)
      allocate( nfb(intmax), stat=stats )
      call umalloc(intmax,1,'nfb',memuse,stats)
      allocate( ni(mbloc), stat=stats )
      call umalloc(mbloc,1,'ni',memuse,stats)
      allocate( nj(mbloc), stat=stats )
      call umalloc(mbloc,1,'nj',memuse,stats)
      allocate( nk(mbloc), stat=stats )
      call umalloc(mbloc,1,'nk',memuse,stats)
      allocate( npts(mbloc), stat=stats )
      call umalloc(mbloc,1,'npts',memuse,stats)
      allocate( nseg(mbloc), stat=stats )
      call umalloc(mbloc,1,'nseg',memuse,stats)
      allocate( nseg1(mbloc), stat=stats )
      call umalloc(mbloc,1,'nseg1',memuse,stats)
      allocate( nsihi(mbloc), stat=stats )
      call umalloc(mbloc,1,'nsihi',memuse,stats)
      allocate( nsilo(mbloc), stat=stats )
      call umalloc(mbloc,1,'nsilo',memuse,stats)
      allocate( nsjhi(mbloc), stat=stats )
      call umalloc(mbloc,1,'nsjhi',memuse,stats)
      allocate( nsjlo(mbloc), stat=stats )
      call umalloc(mbloc,1,'nsjlo',memuse,stats)
      allocate( nskhi(mbloc), stat=stats )
      call umalloc(mbloc,1,'nskhi',memuse,stats)
      allocate( nsklo(mbloc), stat=stats )
      call umalloc(mbloc,1,'nsklo',memuse,stats)
      allocate( nsubbl(mbloc), stat=stats )
      call umalloc(mbloc,1,'nsubbl',memuse,stats)
      allocate( nxtseg(mtot), stat=stats )
      call umalloc(mtot,1,'nxtseg',memuse,stats)
      allocate( rkap0g(mbloc,3), stat=stats )
      call umalloc(mbloc*3,0,'rkap0g',memuse,stats)
      allocate( twotref(msegn,mbloc), stat=stats )
      call umalloc(msegn*mbloc,0,'twotref',memuse,stats)
      allocate( val(mtot), stat=stats )
      call umalloc(mtot,0,'val',memuse,stats)
      allocate( x(npmax), stat=stats )
      call umalloc(npmax,0,'x',memuse,stats)
      allocate( xdmold(msegt), stat=stats )
      call umalloc(msegt,0,'xdmold',memuse,stats)
      allocate( xdum(msegt,mtot), stat=stats )
      call umalloc(msegt*mtot,0,'xdum',memuse,stats)
      allocate( xif1(intmax,nsub1), stat=stats )
      call umalloc(intmax*nsub1,1,'xif1',memuse,stats)
      allocate( xif2(intmax,nsub1), stat=stats )
      call umalloc(intmax*nsub1,1,'xif2',memuse,stats)
      allocate( xit1(intmax), stat=stats )
      call umalloc(intmax,1,'xit1',memuse,stats)
      allocate( xit2(intmax), stat=stats )
      call umalloc(intmax,1,'xit2',memuse,stats)
      allocate( xmap(msegt,msegn,mbloc), stat=stats )
      call umalloc(msegt*msegn*mbloc,0,'xmap',memuse,stats)
      allocate( y(npmax), stat=stats )
      call umalloc(npmax,0,'y',memuse,stats)
      allocate( z(npmax), stat=stats )
      call umalloc(npmax,0,'z',memuse,stats)
c
c     set names for the temporary tlns3d map files
c
      imapin = 'tlnsmap.in'
      imapout= 'tlnsmap.out'
c
c     output banner
c
      write(6,83)
      write(6,83)
      write(6,87)
      write(6,9900)
 9900 format(2(2h *),40h        SPLITTER - CFL3D BLOCK AND INPUT,
     .14h FILE SPLITTER,8x,2(2h *))
      write(6,87)
      write(6,9990)
 9990 format(2(2h *),43h   VERSION 6.7 :  Computational Fluids Lab,,
     .15h Mail Stop 128,,4x,2(2h *),
     ./2(2h *),18x,41hNASA Langley Research Center, Hampton, VA,
     .3x,2(2h *),/2(2h *),18x,33hRelease Date:  February  1, 2017.,
     .11x,2(2h *))
      write(6,87)
      write(6,83)
      write(6,83)
   83 format(35(2h *))
   87 format(2(2h *),62x,2(2h *))
c
#ifdef WKSTN_OFF
c     wkstn_off implies cray (except t3e)
      write(6,12) float(memuse)/1.e6
#else
#   ifdef DBLE_PRECSN
      write(6,12) float(memuse)/1.e6
#   else
      write(6,13) float(memuse)/1.e6
#   endif
#endif
   12 format(/,' memory allocation: ',f12.6,' Mbytes, double precision')
   13 format(/,' memory allocation: ',f12.6,' Mbytes, single precision')
c
c     read preliminary info (file names, etc.)
c
      read(iunit5,*)
      read(iunit5,'(a60)') cflinp
      read(iunit5,'(a60)') roninp
      read(iunit5,'(a60)') gridin
      read(iunit5,'(a60)') sdgridin
      read(iunit5,*)
      read(iunit5,*) icflin,ironin,ibin,isdin
      if(icflin.eq.-5) icflin=5
      read(iunit5,*)
      read(iunit5,'(a60)') cflout
      read(iunit5,'(a60)') ronout
      read(iunit5,'(a60)') gridout
      read(iunit5,'(a60)') sdgridout
      read(iunit5,*)
      read(iunit5,*) icflout,ironout,ibout,isdout
      if(icflout.eq.-5) icflout=5
      read(iunit5,*)
c
      write(*,*)
      write(*,'("input (unsplit) files")')
      write(*,'(''  '',a60)') cflinp
      write(*,'(''  '',a60)') roninp
      write(*,'(''  '',a60)') gridin
      write(*,'(''  '',a60)') sdgridin
      write(*,'("icflver   ironver   igrdfmt    isdfmt")')
      write(*,'(i7,i10,i10,i10)') icflin,ironin,ibin,isdin
      write(*,'("output (split) files")')
c     if any one of the input files is 'null', then corresponding
c     output file must be 'null' too! also, can't split ronnie
c     input files without also splitting cfl3d input file
      if (cflinp.eq.'null') then
         write(*,'("STOPPING: must provide a cfl3d input file!")')
         call termn8(0,-2,ibufdim,nbuf,bou,nou)
      end if
      if (roninp.eq.'null')   ronout    = 'null'
      if (gridin.eq.'null')   gridout   = 'null' 
      if (sdgridin.eq.'null') then
         sdgridout = 'null'
         dovrlap   = 'null'
         dpatch    = 'null'
         dresid    = 'null'
      end if
      write(*,'(''  '',a60)') cflout
      write(*,'(''  '',a60)') ronout
      write(*,'(''  '',a60)') gridout
      write(*,'(''  '',a60)') sdgridout
      write(*,'("icflver   ironver   igrdfmt    isdfmt")')
      write(*,'(i7,i10,i10,i10)') icflout,ironout,ibout,isdout
      write(*,*)
c
c     parallel input version 4 of cfl3d denoted with icflin < 0
c     and or icflout < 0 on input; there is no distinction between
c     parallel and sequential input files for version 6/5
c
      iparin  = 0
      iparout = 0
      if (icflin .lt. 0) then
         iparin = 1
         icflin = abs(icflin)
      end if
      if (icflout .lt. 0) then
         iparout = 1
         icflout = abs(icflout)
      end if
      if (abs(icflin) .gt. 4) then
          icflin = abs(icflin)
      end if
      if (abs(icflout) .gt. 4) then
          icflout = abs(icflout) 
      end if
c
c     open files
c
      open(10,file=cflinp,form="formatted")
      rewind 10
      if (cflout .ne. 'null') then
         open(7,file=cflout,form="formatted")
         rewind 7
      end if
      ipatch0 = 0
      if (roninp .ne. 'null') then
         ipatch0 = 1
         open(15,file=roninp,form="formatted")
         rewind 15
      end if
      if (ronout .ne. 'null') then
         open(9,file=ronout,form="formatted")
         rewind 9
      end if
      open(20,file=imapin,form="formatted")
      rewind 20
      if (gridin .ne. 'null') then
         if(ibin.eq.1) then
c           the following line is needed to use ieee binary file
c           call asnfile (gridin, '-F f77 -N ieee' , IER)
            open(30,file= gridin,form="unformatted")
         else
            open(30,file= gridin,form="formatted")
         end if
         rewind 30
      end if
      if (gridout .ne. 'null') then
         if (ibout.eq.1) then
c           the following line is needed to use ieee binary file
c           call asnfile (gridout, '-F f77 -N ieee' , IER)
            open(40,file=gridout,form="unformatted")
         else
            open(40,file=gridout,form="formatted")
         end if
         rewind 40
      end if
      open(8,file=imapout,form="formatted")
      rewind 8
      isd = 1
      if (sdgridin .ne. 'null') then
         if (isdin .eq. 1) then
            open(50,file=sdgridin,form="unformatted")
         else
            open(50,file=sdgridin,form="formatted")
         end if
         rewind 50 
      end if
      if (sdgridout .ne. 'null') then
         if (isdout .eq. 1) then
            open(60,file=sdgridout,form="unformatted")
         else
            open(60,file=sdgridout,form="formatted")
         end if
         rewind 60
      end if
c
c     convert cfl3d input file (and ronnie input file, if 
c     applicable) to a tlns3d map file
c
      write(*,'("converting unsplit cfl3d input file to ",
     ."tlns3d map file")')
      write(*,'(1x)')
      call cfl3d_to_tlns3d(icflin,ipatch0,ironin,iparin,nnodes,isd,
     .                     mbloc,msegn,msegt,mxbli,il,jl,kl,
     .                     rkap0g,levelg,igridg,iflimg,ifdsg,
     .                     iviscg,jdimg,kdimg,idimg,idiagg,
     .                     nblcg,idegg,jsg,ksg,jeg,keg,
     .                     ieg,mit,ilamlog,ilamhig,jlamlog,
     .                     jlamhig,klamlog,klamhig,
     .                     iredundant,nseg1,nseg,nsilo,nsihi,
     .                     nsjlo,nsjhi,nsklo,nskhi,nb1,ne1,
     .                     nb2,ne2,ibct,bcval,iovrlp,ifoflg,
     .                     ndat,xmap,imap,ivisb,twotref,itrb1,
     .                     itrb2,jtrb1,jtrb2,ktrb1,ktrb2,iturbb,
     .                     ibif1,ibif2,nbi1,nei1,nbj1,nej1,nbk1,
     .                     nek1,nbi2,nei2,nbj2,nej2,nbk2,nek2,
     .                     nd11,nd21,nd12,nd22,ifsor,nb1s,ne1s,
     .                     nb2s,ne2s,nrotat,ntrans,tlref,rlref,ioflag,
     .                     nou,bou,nbuf,ibufdim,cflout)
      rewind 20
c
      nbl = nbloc
c
c     check dimensions
c
      write(*,'("checking dimensions...")')
      write(*,'(1x)')
c
      if (nbl .gt. mbloc) then
         write(*,'("Number of blocks (",i3,") is larger than ",
     .           "the dimensioned value (",i3,")")') nbl,mbloc
         call termn8(0,-3,ibufdim,nbuf,bou,nou)
      end if
c
      iflag = 0
      ibeg(1) = 1
      npts(1) = il(1) * jl(1) * kl(1)
      ntot = ibeg(1) + npts(1) - 1
      do 50 n=2,nbl
         ibeg(n) = ibeg(n-1) + npts(n-1)
         npts(n) = il(n) * jl(n) * kl(n)
         ntot = ibeg(n) + npts(n) - 1
   50 continue
      if(ntot .gt. npmax) then
         write(*,'("number of points ",i8," is larger than ",
     .           "the dimensioned value ",i8)') ntot,npmax
         iflag = 1
      end if
c
      if (iflag .eq. 1) call termn8(0,-4,ibufdim,nbuf,bou,nou)
c
c     read in the input grid file; if input grid file is 'null', get
c     needed block dimensions from unsplit cfl3d input file
c
      if (gridin .ne. 'null') then
c
         if (ip3dgrd .gt. 0) then
            if (ibin.eq.1) then
               read(30) nbl1
               if (nbl1.ne.nbl) then
                  call termn8(0,-1,ibufdim,nbuf,bou,nou)
               end if
               read(30) (ilo(n),jlo(n),klo(n),n=1,nbl)
            else
               read(30,*) nbl1
               if (nbl1.ne.nbl) then
                  call termn8(0,-1,ibufdim,nbuf,bou,nou)
               end if
               read(30,*) (ilo(n),jlo(n),klo(n),n=1,nbl)
            end if
            do n=1,nbl
               if (ilo(n).ne.il(n) .or. jlo(n).ne.jl(n) .or.
     .             klo(n).ne.kl(n)) then
                   call termn8(0,-1,ibufdim,nbuf,bou,nou)
               end if
            end do
         end if
c
         write(*,'("reading grid...")')
         write(*,'("grid: ",a60)') gridin
         write(*,'(1x)')
c
         if (ibin.eq.1) then
            do 100 n=1,nbl
            if (ip3dgrd.eq.0) then
               read(30) jlo(n),klo(n),ilo(n)
               if (ilo(n).ne.il(n) .or. jlo(n).ne.jl(n) .or.
     .             klo(n).ne.kl(n)) then
                   call termn8(0,-1,ibufdim,nbuf,bou,nou)
               end if
            end if
            call ingrd_cmplx(30,ibin,x(ibeg(n)),y(ibeg(n)),z(ibeg(n)),
     .                       ilo(n),jlo(n),klo(n),ip3dgrd)
  100       continue
         else
            do 110 n=1,nbl
            if (ip3dgrd.eq.0) then
               read(30,*) jlo(n),klo(n),ilo(n)
               if (ilo(n).ne.il(n) .or. jlo(n).ne.jl(n) .or.
     .             klo(n).ne.kl(n)) then
                   call termn8(0,-1,ibufdim,nbuf,bou,nou)
               end if
            end if
            call ingrd_cmplx(30,ibin,x(ibeg(n)),y(ibeg(n)),z(ibeg(n)),
     .                 ilo(n),jlo(n),klo(n),ip3dgrd)
  110       continue
         end if
         ntot = ibeg(nbl) + npts(nbl) - 1
c
      else
c
         do n=1,nbl
            ilo(n) = il(n)
            jlo(n) = jl(n)
            klo(n) = kl(n)
         end do
c
      end if
c
      do 200 n=1,nbl
         imin(n) = 0
         imax(n) = 0
         iskp(n) = 0
         jmin(n) = 0
         jmax(n) = 0
         jskp(n) = 0
         kmin(n) = 0
         kmax(n) = 0
         kskp(n) = 0
 200  continue
c
      do 120 n=1,nbl
      write(*,'("block #",i4,": il=",i4,", jl=",i4,", kl=",i4)')
     .                 n,       ilo(n),    jlo(n),    klo(n)
         iln(n) = ilo(n)
         jln(n) = jlo(n)
         kln(n) = klo(n)
         nbs(n) = n
  120 continue
c
c     read in the tlns3d map file
c
      call mapin(20,6,nbl,iln,jln,kln,mbloc,msegt,mtot,
     .           imap,idbloc,ivisb,itrb,val,xdum,iold,
     .           nxtseg,intrfc,ipatch,nsubbl,idobl,nseg,idno,
     .           ijk,idseg,idnext,nou,bou,nbuf,ibufdim)
c
c     read splitting directives and do the splitting
c
      read(iunit5,*) nout
c     shortcut to split all blocks the same as the single block input
      if (nout.lt.0) then
         nout = abs(nout)
         do no = 1,nout
            read(iunit5,*) isblocdum
            read(iunit5,*) ndirin0(no)
            read(iunit5,*) isin0(no)
         end do
         nn = 0
         do nb = 1,nbl
            do no = 1,nout
               nn = nn + 1
               isblocin(nn) = nb
               ndirin(nn)   = ndirin0(no)
               isin(nn)     = isin0(no)
            end do
         end do
         nout = nn
      else
         do nn = 1,nout
            read(iunit5,*) isblocin(nn)
            read(iunit5,*) ndirin(nn)
            read(iunit5,*) isin(nn)
         end do
      end if
      write(*,'(1x)')
      write(*,'("number of splits = ",i3)') nout
      ntotal = nbl + nout
      if (ntotal.gt.mbloc) then
         write(*,'("Number of output blocks (",i3,") is larger than ",
     .             "the dimensioned value (",i3,")")') nbl,mbloc
         call termn8(0,-3,ibufdim,nbuf,bou,nou)
      endif
      write(*,'("  split  block  coord  index")') 
      do 500 no=1,nout
c
      isbloc = isblocin(no)
      ndir   = ndirin(no)
      is     = isin(no)
      if (is.eq.0) then
         if (ndir.eq.1) then
            is = (iln(isbloc)+1) / 2
         else if (ndir.eq.2) then
            is = (jln(isbloc)+1) / 2
         else
            is = (kln(isbloc)+1) / 2
         endif
      endif
      if (ndir.eq.1) then
         write(*,'(i7,i7,"      I",i7)') no,isbloc,is
      else if (ndir.eq.2) then
         write(*,'(i7,i7,"      J",i7)') no,isbloc,is
      else
         write(*,'(i7,i7,"      K",i7)') no,isbloc,is
      endif
c
      call mapsplt(isbloc,ndir,is,mbloc,msegt,mtot,
     .                   imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                   nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                   nseg,idno,ijk,idseg,idnext,xdmold,idmold)
c
  500 continue
c
      call mapout(8,nbloc,ni,nj,nk,mbloc,msegt,mtot,
     .            imap,idbloc,ivisb,itrb,val,xdum,iold,
     .            nxtseg,intrfc,ipatch,nsubbl,idobl,nseg,idno,
     .            ijk,idseg,idnext)
c
      call mapblk(nblk,nbs,imin,imax,jmin,jmax,kmin,kmax,mbloc,
     .            msegt,mtot,imap,idbloc,ivisb,itrb,val,xdum,
     .            iold,nxtseg,intrfc,ipatch,nsubbl,idobl,
     .            nseg,idno,ijk,idseg,idnext)
      do 125 n=1,nblk
         iskp(n) = 1
         jskp(n) = 1
         kskp(n) = 1
         iln(n) = (imax(n)-imin(n))/iskp(n) + 1
         jln(n) = (jmax(n)-jmin(n))/jskp(n) + 1
         kln(n) = (kmax(n)-kmin(n))/kskp(n) + 1
 125  continue
c
c     check multigridability of split-grid basic dimensions
c
      ncgmax = 0
      iflg = 1
      do nnn=1,5
         nn = nnn-1
         do n=1,nblk
            if (iln(n) .gt. 2) then
               idtest = (iln(n)-1)/2**nn + 1
            else
              idtest = 1
            end if
            jdtest = (jln(n)-1)/2**nn + 1
            kdtest = (kln(n)-1)/2**nn + 1
            if ((idtest/2*2 .eq. idtest) .or. 
     .          (jdtest/2*2 .eq. jdtest) .or.
     .          (kdtest/2*2 .eq. kdtest)) then
                iflg = 0
            end if
         end do 
         if (iflg .ne. 0) ncgmax = ncgmax + 1
      end do
c
c     check that coarsened dimensions don't drop below 3 in
c     j or k directions (or i if not a 2d case) - 2 cells
c     in any direction other than i (for 2d) causes problems
c     with 2nd order scheme since 2 ghost cells cannot be
c     filled properly at block interfaces
c
      iflg = 1
      ncgmax1 = 0
      if (ncgmax.gt.0) then
         do nnn=1,ncgmax
            nn = nnn
            do n=1,nblk
               idtest = (iln(n)-1)/2**nn + 1
               jdtest = (jln(n)-1)/2**nn + 1
               kdtest = (kln(n)-1)/2**nn + 1
               if (iln(n) .gt. 2) then
                  if (idtest .le. 2) iflg = 0
               end if
               if (jdtest .le. 2 .or. kdtest .le. 2) iflg = 0 
            end do
            if (iflg .ne. 0) ncgmax1 = ncgmax1 + 1
         end do
      end if
c
      ncgmax = min(ncgmax,ncgmax1)
c
      write(*,'(1x)')
      write(*,'("split-grid basic dimensions are multigridable",
     .          " to ncg = ",i2)') ncgmax
      if (iflg.eq.0) then
         write(*,'("NOTE: ncg was limited to prevent coarsened",
     .             " dimensions less than 3")')
      end if
c
c     print out new grid
c
      if (gridout .ne. 'null') then
         if (ibout.eq.0) then
c
            if (ip3dgrd.gt.0) then   
               write(40,'(i5)') nblk
               write(40,'(3i5)') (iln(n),jln(n),kln(n),n=1,nblk)
            end if
            do 130 n=1,nblk
            if (ip3dgrd.eq.0) then
                write(40,'(3i5)') jln(n),kln(n),iln(n)
            end if 
            isk = iskp(n)
            jsk = jskp(n)
            ksk = kskp(n)
            nn = nbs(n)
            call outgrd_cmplx(40,ibout,x(ibeg(nn)),y(ibeg(nn)),
     .                        z(ibeg(nn)),ilo(nn),jlo(nn),klo(nn),
     .                        imin(n),imax(n),isk,
     .                        jmin(n),jmax(n),jsk,
     .                        kmin(n),kmax(n),ksk,ip3dgrd)
  130       continue
c
         else
c
            if (ip3dgrd.gt.0) then
               write(40) nblk
               write(40) (iln(n),jln(n),kln(n),n=1,nblk)
            end if
            do 140 n=1,nblk
            if (ip3dgrd.eq.0) then
                write(40) jln(n),kln(n),iln(n)
            end if
            isk = iskp(n)
            jsk = jskp(n)
            ksk = kskp(n)
            nn = nbs(n)
            call outgrd_cmplx(40,ibout,x(ibeg(nn)),y(ibeg(nn)),
     .                        z(ibeg(nn)),ilo(nn),jlo(nn),klo(nn),
     .                        imin(n),imax(n),isk,
     .                        jmin(n),jmax(n),jsk,
     .                        kmin(n),kmax(n),ksk,ip3dgrd)
  140       continue
c
         end if
      end if
c
      np = 0
      do 600 n=1,nbl
        np = np + ilo(n)*jlo(n)*klo(n)
 600  continue
      write(*,'(1x)')
      write(*,'("Input  points: ",i8)') np
      np = 0
      do 601 n=1,nblk
        np = np + iln(n)*jln(n)*kln(n)
 601  continue
      write(*,'("Ouput  points: ",i8)') np
      write(*,'(1x)')
c
c     create split cfl3d input file (and split ronnie
c     input file if applicable)
c
      if (cflout .ne. 'null') then
         write(*,'("converting split tlns3d map file to ",
     .   "cfl3d input file")')
c        set number of compute nodes = number of (split) blocks
c        for parallel computations
         nnodes = nblk
         call tlns3d_to_cfl3d(ibout,icflout,ipatch0,ironout,nblk,
     .        iln,jln,kln,iparout,nnodes,isd,ifrom,xif1,xif2,etf1,
     .        etf2,ito,xit1,xit2,ett1,ett2,nfb,iredundant,iconcat,
     .        xmap,msegn,mbloc,iovrlp,nseg,imap,ivisb,twotref,
     .        itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,iturbb,
     .        cq,cporous,cpchamber,cradiation,
     .        llimit,iitmax,mmcxie,mmceta,iifit,iic0,iiorph,
     .        bcvali,bcvalj,bcvalk,nbci0,nbcidim,nbcj0,nbcjdim,
     .        nbck0,nbckdim,ibcinfo,jbcinfo,kbcinfo,nblkk,limblk,
     .        isva,nblon,rkap0g,levelg,igridg,iflimg,ifdsg,iviscg,
     .        jdimg,kdimg,idimg,idiagg,nblcg,idegg,jsg,ksg,isg,
     .        jeg,keg,ieg,mit,ilamlog,ilamhig,jlamlog,jlamhig,
     .        klamlog,klamhig,ncgg,nblg,iemg,inewgg,iadvance,iforce,
     .        intmax,nsub1,msegt,maxseg,mxbli,nrotat,ntrans,tlref,
     .        rlref,ioflag,nou,bou,nbuf,ibufdim)
      end if
c
c     read and split SD grid, if needed
c
      if (sdgridin.eq.'null' .or. sdgridout.eq.'null') isd = 0
c
      if (isd .gt. 0) then
c
         write(*,'(1x)')
         write(*,'("reading SD grid...")')
c
         if (isdin .eq. 0) then
c
            read(50,'(a80)') header
            read(50,*) nblsd,ndv
            read(50,*) (ilosd(n),jlosd(n),klosd(n),n=1,nblsd)
         else
c
            read(50) header
            read(50) nblsd,ndv
            read(50) (ilosd(n),jlosd(n),klosd(n),n=1,nblsd)
         end if
c
         write(*,'("SD grid: ",a60)') sdgridin
         write(*,'(1x)')
         write(*,'("number of design variables = ",i4)') ndv
         do n=1,nblsd
         write(*,'("block #",i4,": il=",i4,", jl=",i4,", kl=",i4)')
     .                    n,       ilosd(n),  jlosd(n),  klosd(n)
         end do
         write(*,'(1x)')
c
c        check consistancy with grid file
c
         if (nblsd .ne. nbl) then
            write(6,*)
     .      'stopping: number of blocks in grid file    = ',nbl
            write(6,*)
     .      'stopping: number of blocks in grid SD file = ',nblsd
            call termn8(0,-5,ibufdim,nbuf,bou,nou)
         end if
         do n=1,nbl
            if (ilo(n) .ne. ilosd(n) .or. jlo(n) .ne. jlosd(n) .or.
     .         klo(n) .ne. klosd(n)) then
               write(6,*)
     .         'stopping: mismatch in grid dimensions and SD ',
     .         'dimensions in block ',nbl
               write(6,*)'grid i,j,k = ',ilo(n),jlo(n),klo(n)
               write(6,*)'SD   i,j,k = ',ilosd(n),jlosd(n),klosd(n)
               call termn8(0,-5,ibufdim,nbuf,bou,nou)
            end if
         end do
c
         if (sdgridout .ne. 'null') then
            if (isdout .eq. 0) then
               write(60,'(a80)') header
               write(60,*) nblk,ndv
               write(60,*) (iln(n),jln(n),kln(n),n=1,nblk)
            else
               write(60) header
               write(60) nblk,ndv
               write(60) (iln(n),jln(n),kln(n),n=1,nblk)
            end if
         end if
c
            do 740 n=1,nblk
c              write(*,'(1x,8i5)') n,imin(n),imax(n),jmin(n),jmax(n),
c     .                            kmin(n),kmax(n),nbs(n)
               isk = iskp(n)
               jsk = jskp(n)
               ksk = kskp(n)
               nn  = nbs(n)
c
c           read grid sensitivity derivatives for (unsplit) block nn
c
               rewind(50)
               if (isdin .eq. 0) then 
                  read(50,*) 
                  read(50,*)
                  read(50,*)
               else
                  read(50)
                  read(50)
                  read(50)
               end if
               do mm=1,nbl
               if (mm.ne.nn) then
                  ilos = ilosd(mm)
                  jlos = jlosd(mm)
                  klos = klosd(mm)
                  do idv=1,ndv
                     if (isdin .eq. 0) then
                        read(50,*) 
     .                  (((xdum0,i=1,ilos),j=1,jlos),k=1,klos),
     .                  (((ydum0,i=1,ilos),j=1,jlos),k=1,klos),
     .                  (((zdum0,i=1,ilos),j=1,jlos),k=1,klos)
                     else
                        read(50)
     .                  (((xdum0,i=1,ilos),j=1,jlos),k=1,klos),
     .                  (((ydum0,i=1,ilos),j=1,jlos),k=1,klos),
     .                  (((zdum0,i=1,ilos),j=1,jlos),k=1,klos)
                     end if
                  end do
               else
                  do idv=1,ndv 
                     if (isdin .eq. 0) then
                        read(50,*)
     .                  (x(i),i=ibeg(nn),ibeg(nn)+npts(nn)-1),
     .                  (y(i),i=ibeg(nn),ibeg(nn)+npts(nn)-1),
     .                  (z(i),i=ibeg(nn),ibeg(nn)+npts(nn)-1)
                     else
                        read(50)
     .                  (x(i),i=ibeg(nn),ibeg(nn)+npts(nn)-1),
     .                  (y(i),i=ibeg(nn),ibeg(nn)+npts(nn)-1),
     .                  (z(i),i=ibeg(nn),ibeg(nn)+npts(nn)-1)
                     end if
                     call outgrd_cmplx(60,isdout,x(ibeg(nn)),
     .                                 y(ibeg(nn)),z(ibeg(nn)),
     .                                 ilo(nn),jlo(nn),klo(nn),
     .                                 imin(n),imax(n),isk,
     .                                 jmin(n),jmax(n),jsk,
     .                                 kmin(n),kmax(n),ksk,1)
                  end do
               end if
            end do
c
  740       continue
c
      end if
c
c     free up memory
c
      deallocate(ifrom)
      deallocate(xif1)
      deallocate(xif2)
      deallocate(etf1)
      deallocate(etf2)
      deallocate(ito)
      deallocate(xit1)
      deallocate(xit2)
      deallocate(ett1)
      deallocate(ett2)
      deallocate(nfb)
      deallocate(iadvance)
      deallocate(iforce)
      deallocate(nblkk)
      deallocate(limblk)
      deallocate(isva)
      deallocate(nblon)
      deallocate(ncgg)
      deallocate(nblg)
      deallocate(iemg)
      deallocate(inewgg)
      deallocate(bcvali)
      deallocate(bcvalj)
      deallocate(bcvalk)
      deallocate(nbci0)
      deallocate(nbcidim)
      deallocate(nbcj0)
      deallocate(nbcjdim)
      deallocate(nbck0)
      deallocate(nbckdim)
      deallocate(ibcinfo)
      deallocate(jbcinfo)
      deallocate(kbcinfo)
      deallocate(llimit)
      deallocate(iitmax)
      deallocate(mmcxie)
      deallocate(mmceta)
      deallocate(iifit)
      deallocate(iic0)
      deallocate(iiorph)
      deallocate(iredundant)
      deallocate(iconcat)
      deallocate(cq)
      deallocate(cporous)
      deallocate(cpchamber)
      deallocate(cradiation)
      deallocate(isblocin)
      deallocate(ndirin)
      deallocate(isin)
      deallocate(ndirin0)
      deallocate(isin0)
      deallocate(nbs)
      deallocate(imin)
      deallocate(jmin)
      deallocate(kmin)
      deallocate(imax)
      deallocate(jmax)
      deallocate(kmax)
      deallocate(iskp)
      deallocate(jskp)
      deallocate(kskp)
      deallocate(ilo)
      deallocate(jlo)
      deallocate(klo)
      deallocate(iln)
      deallocate(jln)
      deallocate(kln)
      deallocate(ibeg)
      deallocate(npts)
      deallocate(x)
      deallocate(y)
      deallocate(z)
      deallocate(ilosd)
      deallocate(jlosd)
      deallocate(klosd)
      deallocate(il)
      deallocate(jl)
      deallocate(kl)
      deallocate(rkap0g)
      deallocate(levelg)
      deallocate(igridg)
      deallocate(iflimg)
      deallocate(ifdsg)
      deallocate(iviscg)
      deallocate(jdimg)
      deallocate(kdimg)
      deallocate(idimg)
      deallocate(idiagg)
      deallocate(nblcg)
      deallocate(idegg)
      deallocate(jsg)
      deallocate(ksg)
      deallocate(isg)
      deallocate(jeg)
      deallocate(keg)
      deallocate(ieg)
      deallocate(mit)
      deallocate(ilamlog)
      deallocate(ilamhig)
      deallocate(jlamlog)
      deallocate(jlamhig)
      deallocate(klamlog)
      deallocate(klamhig)
      deallocate(imap)
      deallocate(idbloc)
      deallocate(ivisb)
      deallocate(itrb)
      deallocate(val)
      deallocate(xdum)
      deallocate(iold)
      deallocate(nxtseg)
      deallocate(intrfc)
      deallocate(ipatch)
      deallocate(ni)
      deallocate(nj)
      deallocate(nk)
      deallocate(xdmold)
      deallocate(idmold)
      deallocate(nseg1)
      deallocate(nsilo)
      deallocate(nsihi)
      deallocate(nsjlo)
      deallocate(nsjhi)
      deallocate(nsklo)
      deallocate(nskhi)
      deallocate(nb1)
      deallocate(ne1)
      deallocate(nb2)
      deallocate(ne2)
      deallocate(ibct)
      deallocate(bcval)
      deallocate(iovrlp)
      deallocate(ifoflg)
      deallocate(ndat)
      deallocate(xmap)
      deallocate(twotref)
      deallocate(itrb1)
      deallocate(itrb2)
      deallocate(jtrb1)
      deallocate(jtrb2)
      deallocate(ktrb1)
      deallocate(ktrb2)
      deallocate(iturbb)
      deallocate(ibif1)
      deallocate(ibif2)
      deallocate(nbi1)
      deallocate(nei1)
      deallocate(nbj1)
      deallocate(nej1)
      deallocate(nbk1)
      deallocate(nek1)
      deallocate(nbi2)
      deallocate(nei2)
      deallocate(nbj2)
      deallocate(nej2)
      deallocate(nbk2)
      deallocate(nek2)
      deallocate(nd11)
      deallocate(nd21)
      deallocate(nd12)
      deallocate(nd22)
      deallocate(ifsor)
      deallocate(nb1s)
      deallocate(ne1s)
      deallocate(nb2s)
      deallocate(ne2s)
      deallocate(nsubbl)
      deallocate(idobl)
      deallocate(nseg)
      deallocate(idno)
      deallocate(ijk)
      deallocate(idseg)
      deallocate(idnext)
c
      write(*,'("SPLITTING/CONVERSION COMPLETED")')
c
      return
      end
