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
      integer ifrom(intmax,nsub1),xif1(intmax,nsub1),
     .        xif2(intmax,nsub1),etf1(intmax,nsub1),
     .        etf2(intmax,nsub1),ito(intmax),xit1(intmax),
     .        xit2(intmax),ett1(intmax),ett2(intmax),nfb(intmax)
c
      dimension nou(nbuf)
      dimension iadvance(mbloc),iforce(mbloc)
      dimension nblkk(2,mxbli),limblk(2,6,mxbli),
     .          isva(2,2,mxbli),nblon(mxbli)
      dimension ncgg(mbloc),nblg(mbloc),iemg(mbloc),
     .          inewgg(mbloc)
      dimension bcvali(mbloc,maxseg,7,2),
     .          bcvalj(mbloc,maxseg,7,2),bcvalk(mbloc,maxseg,7,2),
     .          nbci0(mbloc),nbcidim(mbloc),nbcj0(mbloc),nbcjdim(mbloc),
     .          nbck0(mbloc),nbckdim(mbloc),ibcinfo(mbloc,maxseg,7,2),
     .          jbcinfo(mbloc,maxseg,7,2),kbcinfo(mbloc,maxseg,7,2)
      dimension llimit(intmax),iitmax(intmax),mmcxie(intmax),
     .          mmceta(intmax),iifit(intmax),iic0(intmax),
     .          iiorph(intmax)
      dimension iredundant(mbloc,msegn),iconcat(intmax)
      dimension cq(msegn,mbloc),cporous(msegn,mbloc),
     .          cpchamber(msegn,mbloc),cradiation(msegn,mbloc)
      dimension isblocin(msplt*mbloc),ndirin(msplt*mbloc),
     .          isin(msplt*mbloc),ndirin0(msplt),isin0(msplt)
      dimension  nbs(mbloc)
      dimension imin(mbloc),jmin(mbloc),kmin(mbloc)
      dimension imax(mbloc),jmax(mbloc),kmax(mbloc)
      dimension iskp(mbloc),jskp(mbloc),kskp(mbloc)
      dimension ilo(mbloc),jlo(mbloc),klo(mbloc)
      dimension iln(mbloc),jln(mbloc),kln(mbloc)
      dimension ibeg(mbloc),npts(mbloc)
      complex   x(npmax),y(npmax),z(npmax)
      dimension ilosd(mbloc),jlosd(mbloc),klosd(mbloc)
      dimension il(mbloc),jl(mbloc),kl(mbloc),rkap0g(mbloc,3),
     .          levelg(mbloc),igridg(mbloc),iflimg(mbloc,3),
     .          ifdsg(mbloc,3),iviscg(mbloc,3),jdimg(mbloc),
     .          kdimg(mbloc),idimg(mbloc),idiagg(mbloc,3),
     .          nblcg(mbloc),idegg(mbloc,3),jsg(mbloc),ksg(mbloc),
     .          isg(mbloc),jeg(mbloc),keg(mbloc),ieg(mbloc),
     .          mit(5,mbloc),ilamlog(mbloc),ilamhig(mbloc),
     .          jlamlog(mbloc),jlamhig(mbloc),klamlog(mbloc),
     .          klamhig(mbloc)
      dimension imap(msegt,msegn,mbloc),idbloc(mtot),
     .          ivisb(msegt,mbloc),itrb(7,mbloc),
     .          val(mtot),xdum(msegt,mtot),
     .          iold(4,mtot),nxtseg(mtot),intrfc(mtot),
     .          ipatch(mtot),ni(mbloc),nj(mbloc),nk(mbloc)
      dimension xdmold(msegt),idmold(msegt)
      dimension nseg1(mbloc)
      dimension nsilo(mbloc),nsihi(mbloc),nsjlo(mbloc),
     .          nsjhi(mbloc),nsklo(mbloc),nskhi(mbloc)
      dimension nb1(mbloc,6,msegn),ne1(mbloc,6,msegn),
     .          nb2(mbloc,6,msegn),ne2(mbloc,6,msegn),
     .          ibct(mbloc,6,msegn)
      dimension bcval(mbloc,6,msegn,7)
      dimension iovrlp(mbloc),ifoflg(mbloc,6,msegn),
     .          ndat(mbloc,6,msegn),xmap(msegt,msegn,mbloc)
      dimension twotref(msegn,mbloc)
      dimension itrb1(mbloc),itrb2(mbloc),jtrb1(mbloc),jtrb2(mbloc),
     .          ktrb1(mbloc),ktrb2(mbloc),iturbb(mbloc)
      dimension ibif1(mxbli),ibif2(mxbli)
      dimension nbi1(mxbli),nei1(mxbli),nbj1(mxbli),
     .          nej1(mxbli),nbk1(mxbli),nek1(mxbli),
     .          nbi2(mxbli),nei2(mxbli),nbj2(mxbli),
     .          nej2(mxbli),nbk2(mxbli),nek2(mxbli),
     .          nd11(mxbli),nd21(mxbli),nd12(mxbli),
     .          nd22(mxbli)
      dimension ifsor(mbloc),nb1s(mbloc),ne1s(mbloc),
     .          nb2s(mbloc),ne2s(mbloc)
      dimension nsubbl(mbloc),idobl(mbloc),nseg(mbloc),idno(mbloc),
     .          ijk(6,mbloc),idseg(mbloc),idnext(mbloc)
c
      pointer (ip_ifrom,ifrom),(ip_xif1,xif1),(ip_xif2,xif2),
     .        (ip_etf1,etf1),(ip_etf2,etf2),(ip_ito,ito),(ip_xit1,xit1),
     .        (ip_xit2,xit2),(ip_ett1,ett1),(ip_ett2,ett2),(ip_nfb,nfb)
      pointer (ip_iadvance,iadvance),(ip_iforce,iforce),
     .        (ip_nblkk,nblkk),(ip_limblk,limblk),(ip_isva,isva),
     .        (ip_nblon,nblon),(ip_ncgg,ncgg),(ip_nblg,nblg),
     .        (ip_iemg,iemg),(ip_inewgg,inewgg),(ip_bcvali,bcvali),
     .        (ip_bcvalj,bcvalj),(ip_bcvalk,bcvalk),(ip_nbci0,nbci0),
     .        (ip_nbcidim,nbcidim),(ip_nbcj0,nbcj0),
     .        (ip_nbcjdim,nbcjdim),(ip_nbck0,nbck0),
     .        (ip_nbckdim,nbckdim),(ip_ibcinfo,ibcinfo),
     .        (ip_jbcinfo,jbcinfo),(ip_kbcinfo,kbcinfo)
      pointer (ip_llimit,llimit),(ip_iitmax,iitmax),(ip_mmcxie,mmcxie),
     .        (ip_mmceta,mmceta),(ip_iifit,iifit),(ip_iic0,iic0),
     .        (ip_iiorph,iiorph),(ip_iredundant,iredundant),
     .        (ip_iconcat,iconcat),(ip_cq,cq),(ip_cporous,cporous),
     .        (ip_cpchamber,cpchamber),(ip_cradiation,cradiation),
     .        (ip_isblocin,isblocin),(ip_ndirin,ndirin),(ip_isin,isin),
     .        (ip_ndirin0,ndirin0),(ip_isin0,isin0)
      pointer (ip_nbs,nbs),(ip_imin,imin),(ip_jmin,jmin),(ip_kmin,kmin),
     .        (ip_imax,imax),(ip_jmax,jmax),(ip_kmax,kmax),
     .        (ip_iskp,iskp),(ip_jskp,jskp),(ip_kskp,kskp),(ip_ilo,ilo),
     .        (ip_jlo,jlo),(ip_klo,klo),(ip_iln,iln),(ip_jln,jln),
     .        (ip_kln,kln),(ip_ibeg,ibeg),(ip_npts,npts),(ip_x,x),
     .        (ip_y,y),(ip_z,z),(ip_ilosd,ilosd),(ip_jlosd,jlosd),
     .        (ip_klosd,klosd) 
      pointer (ip_il,il),(ip_jl,jl),(ip_kl,kl),(ip_rkap0g,rkap0g),
     .        (ip_levelg,levelg),(ip_igridg,igridg),(ip_iflimg,iflimg),
     .        (ip_ifdsg,ifdsg),(ip_iviscg,iviscg),(ip_jdimg,jdimg),
     .        (ip_kdimg,kdimg),(ip_idimg,idimg),(ip_idiagg,idiagg),
     .        (ip_nblcg,nblcg),(ip_idegg,idegg),(ip_jsg,jsg),
     .        (ip_ksg,ksg),(ip_isg,isg),(ip_jeg,jeg),(ip_keg,keg),
     .        (ip_ieg,ieg),(ip_mit,mit),(ip_ilamlog,ilamlog),
     .        (ip_ilamhig,ilamhig),(ip_jlamlog,jlamlog),
     .        (ip_jlamhig,jlamhig),(ip_klamlog,klamlog),
     .        (ip_klamhig,klamhig)
      pointer (ip_imap,imap),(ip_idbloc,idbloc),(ip_ivisb,ivisb),
     .        (ip_itrb,itrb),(ip_val,val),(ip_xdum,xdum),
     .        (ip_iold,iold),(ip_nxtseg,nxtseg),(ip_intrfc,intrfc),
     .        (ip_ipatch,ipatch),(ip_ni,ni),(ip_nj,nj),(ip_nk,nk),
     .        (ip_xdmold,xdmold),(ip_idmold,idmold),(ip_nseg1,nseg1)

      pointer (ip_nsilo,nsilo),(ip_nsihi,nsihi),(ip_nsjlo,nsjlo),
     .        (ip_nsjhi,nsjhi),(ip_nsklo,nsklo),(ip_nskhi,nskhi),
     .        (ip_nb1,nb1),(ip_ne1,ne1),(ip_nb2,nb2),(ip_ne2,ne2),
     .        (ip_ibct,ibct),(ip_bcval,bcval),(ip_iovrlp,iovrlp),
     .        (ip_ifoflg,ifoflg),(ip_ndat,ndat),(ip_xmap,xmap)
      pointer (ip_twotref,twotref),(ip_itrb1,itrb1),(ip_itrb2,itrb2),
     .        (ip_jtrb1,jtrb1),(ip_jtrb2,jtrb2),(ip_ktrb1,ktrb1),
     .        (ip_ktrb2,ktrb2),(ip_iturbb,iturbb),(ip_ibif1,ibif1),
     .        (ip_ibif2,ibif2),(ip_nbi1,nbi1),(ip_nei1,nei1),
     .        (ip_nbj1,nbj1),(ip_nej1,nej1),(ip_nbk1,nbk1),
     .        (ip_nek1,nek1),(ip_nbi2,nbi2),(ip_nei2,nei2),
     .        (ip_nbj2,nbj2),(ip_nej2,nej2),(ip_nbk2,nbk2),
     .        (ip_nek2,nek2),(ip_nd11,nd11),(ip_nd21,nd21),
     .        (ip_nd12,nd12),(ip_nd22,nd22),(ip_ifsor,ifsor),
     .        (ip_nb1s,nb1s),(ip_ne1s,ne1s),(ip_nb2s,nb2s),
     .        (ip_ne2s,ne2s),(ip_nsubbl,nsubbl),(ip_idobl,idobl),
     .        (ip_nseg,nseg),(ip_idno,idno),(ip_ijk,ijk),
     .        (ip_idseg,idseg),(ip_idnext,idnext)
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
      icall  = 0
      memuse = 0
c
      call umalloc(ip_ifrom,intmax*nsub1,1,icall,memuse)
      call umalloc(ip_xif1,intmax*nsub1,1,icall,memuse)
      call umalloc(ip_xif2,intmax*nsub1,1,icall,memuse)
      call umalloc(ip_etf1,intmax*nsub1,1,icall,memuse)
      call umalloc(ip_etf2,intmax*nsub1,1,icall,memuse)
      call umalloc(ip_ito,intmax,1,icall,memuse)
      call umalloc(ip_xit1,intmax,1,icall,memuse)
      call umalloc(ip_xit2,intmax,1,icall,memuse)
      call umalloc(ip_ett1,intmax,1,icall,memuse)
      call umalloc(ip_ett2,intmax,1,icall,memuse)
      call umalloc(ip_nfb,intmax,1,icall,memuse)
      call umalloc(ip_iadvance,mbloc,1,icall,memuse)
      call umalloc(ip_iforce,mbloc,1,icall,memuse)
      call umalloc(ip_nblkk,2*mxbli,1,icall,memuse)
      call umalloc(ip_limblk,12*mxbli,1,icall,memuse)
      call umalloc(ip_isva,4*mxbli,1,icall,memuse)
      call umalloc(ip_nblon,mxbli,1,icall,memuse)
      call umalloc(ip_ncgg,mbloc,1,icall,memuse)
      call umalloc(ip_nblg,mbloc,1,icall,memuse)
      call umalloc(ip_iemg,mbloc,1,icall,memuse)
      call umalloc(ip_inewgg,mbloc,1,icall,memuse)
      call umalloc(ip_bcvali,14*mbloc*maxseg,0,icall,memuse)
      call umalloc(ip_bcvalj,14*mbloc*maxseg,0,icall,memuse) 
      call umalloc(ip_bcvalk,14*mbloc*maxseg,0,icall,memuse)
      call umalloc(ip_nbci0,mbloc,1,icall,memuse)
      call umalloc(ip_nbcidim,mbloc,1,icall,memuse)
      call umalloc(ip_nbcj0,mbloc,1,icall,memuse)
      call umalloc(ip_nbcjdim,mbloc,1,icall,memuse)
      call umalloc(ip_nbck0,mbloc,1,icall,memuse)
      call umalloc(ip_nbckdim,mbloc,1,icall,memuse)
      call umalloc(ip_ibcinfo,14*mbloc*maxseg,1,icall,memuse)
      call umalloc(ip_jbcinfo,14*mbloc*maxseg,1,icall,memuse)
      call umalloc(ip_kbcinfo,14*mbloc*maxseg,1,icall,memuse)
      call umalloc(ip_llimit,intmax,1,icall,memuse)
      call umalloc(ip_iitmax,intmax,1,icall,memuse)
      call umalloc(ip_mmcxie,intmax,1,icall,memuse)
      call umalloc(ip_mmceta,intmax,1,icall,memuse)
      call umalloc(ip_iifit,intmax,1,icall,memuse)
      call umalloc(ip_iic0,intmax,1,icall,memuse)
      call umalloc(ip_iiorph,intmax,1,icall,memuse)
      call umalloc(ip_iredundant,mbloc*msegn,1,icall,memuse)
      call umalloc(ip_iconcat,intmax,1,icall,memuse)
      call umalloc(ip_cq,msegn*mbloc,0,icall,memuse)
      call umalloc(ip_cporous,msegn*mbloc,0,icall,memuse)
      call umalloc(ip_cpchamber,msegn*mbloc,0,icall,memuse)
      call umalloc(ip_cradiation,msegn*mbloc,0,icall,memuse)
      call umalloc(ip_isblocin,msplt*mbloc,1,icall,memuse)
      call umalloc(ip_ndirin,msplt*mbloc,1,icall,memuse)
      call umalloc(ip_isin,msplt*mbloc,1,icall,memuse)
      call umalloc(ip_ndirin0,msplt,1,icall,memuse)
      call umalloc(ip_isin0,msplt,1,icall,memuse)
      call umalloc(ip_nbs,mbloc,1,icall,memuse)
      call umalloc(ip_imin,mbloc,1,icall,memuse)
      call umalloc(ip_jmin,mbloc,1,icall,memuse)
      call umalloc(ip_kmin,mbloc,1,icall,memuse)
      call umalloc(ip_imax,mbloc,1,icall,memuse)
      call umalloc(ip_jmax,mbloc,1,icall,memuse)
      call umalloc(ip_kmax,mbloc,1,icall,memuse)
      call umalloc(ip_iskp,mbloc,1,icall,memuse)
      call umalloc(ip_jskp,mbloc,1,icall,memuse)
      call umalloc(ip_kskp,mbloc,1,icall,memuse)
      call umalloc(ip_ilo,mbloc,1,icall,memuse)
      call umalloc(ip_jlo,mbloc,1,icall,memuse)
      call umalloc(ip_klo,mbloc,1,icall,memuse)
      call umalloc(ip_iln,mbloc,1,icall,memuse)
      call umalloc(ip_jln,mbloc,1,icall,memuse)
      call umalloc(ip_kln,mbloc,1,icall,memuse)
      call umalloc(ip_ibeg,mbloc,1,icall,memuse)
      call umalloc(ip_npts,mbloc,1,icall,memuse)
      call umalloc_c(ip_x,npmax,0,icall,memuse)
      call umalloc_c(ip_y,npmax,0,icall,memuse)
      call umalloc_c(ip_z,npmax,0,icall,memuse)
      call umalloc(ip_ilosd,mbloc,1,icall,memuse)
      call umalloc(ip_jlosd,mbloc,1,icall,memuse)
      call umalloc(ip_klosd,mbloc,1,icall,memuse)
      call umalloc(ip_il,mbloc,1,icall,memuse)
      call umalloc(ip_jl,mbloc,1,icall,memuse)
      call umalloc(ip_kl,mbloc,1,icall,memuse)
      call umalloc(ip_rkap0g,3*mbloc,0,icall,memuse)
      call umalloc(ip_levelg,mbloc,1,icall,memuse)
      call umalloc(ip_igridg,mbloc,1,icall,memuse)
      call umalloc(ip_iflimg,3*mbloc,1,icall,memuse)
      call umalloc(ip_ifdsg,3*mbloc,1,icall,memuse)
      call umalloc(ip_iviscg,3*mbloc,1,icall,memuse)
      call umalloc(ip_jdimg,mbloc,1,icall,memuse)
      call umalloc(ip_kdimg,mbloc,1,icall,memuse)
      call umalloc(ip_idimg,mbloc,1,icall,memuse)
      call umalloc(ip_idiagg,3*mbloc,1,icall,memuse)
      call umalloc(ip_nblcg,mbloc,1,icall,memuse)
      call umalloc(ip_idegg,3*mbloc,1,icall,memuse)
      call umalloc(ip_jsg,mbloc,1,icall,memuse)
      call umalloc(ip_ksg,mbloc,1,icall,memuse)
      call umalloc(ip_isg,mbloc,1,icall,memuse)
      call umalloc(ip_jeg,mbloc,1,icall,memuse)
      call umalloc(ip_keg,mbloc,1,icall,memuse)
      call umalloc(ip_ieg,mbloc,1,icall,memuse)
      call umalloc(ip_mit,5*mbloc,1,icall,memuse)
      call umalloc(ip_ilamlog,mbloc,1,icall,memuse)
      call umalloc(ip_ilamhig,mbloc,1,icall,memuse)
      call umalloc(ip_jlamlog,mbloc,1,icall,memuse)
      call umalloc(ip_jlamhig,mbloc,1,icall,memuse)
      call umalloc(ip_klamlog,mbloc,1,icall,memuse)
      call umalloc(ip_klamhig,mbloc,1,icall,memuse)
      call umalloc(ip_imap,msegt*msegn*mbloc,1,icall,memuse)
      call umalloc(ip_idbloc,mtot,1,icall,memuse)
      call umalloc(ip_ivisb,msegt*mtot,1,icall,memuse)
      call umalloc(ip_itrb,7*mbloc,1,icall,memuse)
      call umalloc(ip_val,mtot,0,icall,memuse)
      call umalloc(ip_xdum,msegt*mtot,0,icall,memuse)
      call umalloc(ip_iold,4*mtot,1,icall,memuse)
      call umalloc(ip_nxtseg,mtot,1,icall,memuse)
      call umalloc(ip_intrfc,mtot,1,icall,memuse)
      call umalloc(ip_ipatch,mtot,1,icall,memuse)
      call umalloc(ip_ni,mbloc,1,icall,memuse)
      call umalloc(ip_nj,mbloc,1,icall,memuse)
      call umalloc(ip_nk,mbloc,1,icall,memuse)
      call umalloc(ip_xdmold,msegt,0,icall,memuse)
      call umalloc(ip_idmold,msegt,1,icall,memuse)
      call umalloc(ip_nseg1,mbloc,1,icall,memuse)
      call umalloc(ip_nsilo,mbloc,1,icall,memuse)
      call umalloc(ip_nsihi,mbloc,1,icall,memuse)
      call umalloc(ip_nsjlo,mbloc,1,icall,memuse)
      call umalloc(ip_nsjhi,mbloc,1,icall,memuse)
      call umalloc(ip_nsklo,mbloc,1,icall,memuse)
      call umalloc(ip_nskhi,mbloc,1,icall,memuse)
      call umalloc(ip_nb1,6*mbloc*msegn,1,icall,memuse)
      call umalloc(ip_ne1,6*mbloc*msegn,1,icall,memuse)
      call umalloc(ip_nb2,6*mbloc*msegn,1,icall,memuse)
      call umalloc(ip_ne2,6*mbloc*msegn,1,icall,memuse)
      call umalloc(ip_ibct,6*mbloc*msegn,1,icall,memuse)
      call umalloc(ip_bcval,42*mbloc*msegn,0,icall,memuse)
      call umalloc(ip_iovrlp,mbloc,1,icall,memuse)
      call umalloc(ip_ifoflg,6*mbloc*msegn,1,icall,memuse)
      call umalloc(ip_ndat,6*mbloc*msegn,1,icall,memuse)
      call umalloc(ip_xmap,msegt*msegn*mbloc,0,icall,memuse)
      call umalloc(ip_twotref,msegn*mbloc,0,icall,memuse)
      call umalloc(ip_itrb1,mbloc,1,icall,memuse)
      call umalloc(ip_itrb2,mbloc,1,icall,memuse)
      call umalloc(ip_jtrb1,mbloc,1,icall,memuse)
      call umalloc(ip_jtrb2,mbloc,1,icall,memuse)
      call umalloc(ip_ktrb1,mbloc,1,icall,memuse)
      call umalloc(ip_ktrb2,mbloc,1,icall,memuse)
      call umalloc(ip_iturbb,mbloc,1,icall,memuse)
      call umalloc(ip_ibif1,mxbli,1,icall,memuse)
      call umalloc(ip_ibif2,mxbli,1,icall,memuse)
      call umalloc(ip_nbi1,mxbli,1,icall,memuse)
      call umalloc(ip_nei1,mxbli,1,icall,memuse)
      call umalloc(ip_nbj1,mxbli,1,icall,memuse)
      call umalloc(ip_nej1,mxbli,1,icall,memuse)
      call umalloc(ip_nbk1,mxbli,1,icall,memuse)
      call umalloc(ip_nek1,mxbli,1,icall,memuse)
      call umalloc(ip_nbi2,mxbli,1,icall,memuse)
      call umalloc(ip_nei2,mxbli,1,icall,memuse)
      call umalloc(ip_nbj2,mxbli,1,icall,memuse)
      call umalloc(ip_nej2,mxbli,1,icall,memuse)
      call umalloc(ip_nbk2,mxbli,1,icall,memuse)
      call umalloc(ip_nek2,mxbli,1,icall,memuse)
      call umalloc(ip_nd11,mxbli,1,icall,memuse)
      call umalloc(ip_nd21,mxbli,1,icall,memuse)
      call umalloc(ip_nd12,mxbli,1,icall,memuse)
      call umalloc(ip_nd22,mxbli,1,icall,memuse)
      call umalloc(ip_ifsor,mbloc,1,icall,memuse)
      call umalloc(ip_nb1s,mbloc,1,icall,memuse)
      call umalloc(ip_ne1s,mbloc,1,icall,memuse)
      call umalloc(ip_nb2s,mbloc,1,icall,memuse)
      call umalloc(ip_ne2s,mbloc,1,icall,memuse)
      call umalloc(ip_nsubbl,mbloc,1,icall,memuse)
      call umalloc(ip_idobl,mbloc,1,icall,memuse)
      call umalloc(ip_nseg,mbloc,1,icall,memuse)
      call umalloc(ip_idno,mbloc,1,icall,memuse)
      call umalloc(ip_ijk,6*mbloc,1,icall,memuse)
      call umalloc(ip_idseg,mbloc,1,icall,memuse)
      call umalloc(ip_idnext,mbloc,1,icall,memuse)
c
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
 9990 format(2(2h *),43h   VERSION 6.0 :  Computational Fluids Lab,,
     .15h Mail Stop 128,,4x,2(2h *),
     ./2(2h *),18x,41hNASA Langley Research Center, Hampton, VA,
     .3x,2(2h *),/2(2h *),18x,33hRelease Date:    August 14, 2001.,
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
      call ufree(ip_ifrom)
      call ufree(ip_xif1)
      call ufree(ip_xif2)
      call ufree(ip_etf1)
      call ufree(ip_etf2)
      call ufree(ip_ito)
      call ufree(ip_xit1)
      call ufree(ip_xit2)
      call ufree(ip_ett1)
      call ufree(ip_ett2)
      call ufree(ip_nfb)
      call ufree(ip_iadvance)
      call ufree(ip_iforce)
      call ufree(ip_nblkk)
      call ufree(ip_limblk)
      call ufree(ip_isva)
      call ufree(ip_nblon)
      call ufree(ip_ncgg)
      call ufree(ip_nblg)
      call ufree(ip_iemg)
      call ufree(ip_inewgg)
      call ufree(ip_bcvali)
      call ufree(ip_bcvalj)
      call ufree(ip_bcvalk)
      call ufree(ip_nbci0)
      call ufree(ip_nbcidim)
      call ufree(ip_nbcj0)
      call ufree(ip_nbcjdim)
      call ufree(ip_nbck0)
      call ufree(ip_nbckdim)
      call ufree(ip_ibcinfo)
      call ufree(ip_jbcinfo)
      call ufree(ip_kbcinfo)
      call ufree(ip_llimit)
      call ufree(ip_iitmax)
      call ufree(ip_mmcxie)
      call ufree(ip_mmceta)
      call ufree(ip_iifit)
      call ufree(ip_iic0)
      call ufree(ip_iiorph)
      call ufree(ip_iredundant)
      call ufree(ip_iconcat)
      call ufree(ip_cq)
      call ufree(ip_cporous)
      call ufree(ip_cpchamber)
      call ufree(ip_cradiation)
      call ufree(ip_isblocin)
      call ufree(ip_ndirin)
      call ufree(ip_isin)
      call ufree(ip_ndirin0)
      call ufree(ip_isin0)
      call ufree(ip_nbs)
      call ufree(ip_imin)
      call ufree(ip_jmin)
      call ufree(ip_kmin)
      call ufree(ip_imax)
      call ufree(ip_jmax)
      call ufree(ip_kmax)
      call ufree(ip_iskp)
      call ufree(ip_jskp)
      call ufree(ip_kskp)
      call ufree(ip_ilo)
      call ufree(ip_jlo)
      call ufree(ip_klo)
      call ufree(ip_iln)
      call ufree(ip_jln)
      call ufree(ip_kln)
      call ufree(ip_ibeg)
      call ufree(ip_npts)
      call ufree(ip_x)
      call ufree(ip_y)
      call ufree(ip_z)
      call ufree(ip_ilosd)
      call ufree(ip_jlosd)
      call ufree(ip_klosd)
      call ufree(ip_il)
      call ufree(ip_jl)
      call ufree(ip_kl)
      call ufree(ip_rkap0g)
      call ufree(ip_levelg)
      call ufree(ip_igridg)
      call ufree(ip_iflimg)
      call ufree(ip_ifdsg)
      call ufree(ip_iviscg)
      call ufree(ip_jdimg)
      call ufree(ip_kdimg)
      call ufree(ip_idimg)
      call ufree(ip_idiagg)
      call ufree(ip_nblcg)
      call ufree(ip_idegg)
      call ufree(ip_jsg)
      call ufree(ip_ksg)
      call ufree(ip_isg)
      call ufree(ip_jeg)
      call ufree(ip_keg)
      call ufree(ip_ieg)
      call ufree(ip_mit)
      call ufree(ip_ilamlog)
      call ufree(ip_ilamhig)
      call ufree(ip_jlamlog)
      call ufree(ip_jlamhig)
      call ufree(ip_klamlog)
      call ufree(ip_klamhig)
      call ufree(ip_imap)
      call ufree(ip_idbloc)
      call ufree(ip_ivisb)
      call ufree(ip_itrb)
      call ufree(ip_val)
      call ufree(ip_xdum)
      call ufree(ip_iold)
      call ufree(ip_nxtseg)
      call ufree(ip_intrfc)
      call ufree(ip_ipatch)
      call ufree(ip_ni)
      call ufree(ip_nj)
      call ufree(ip_nk)
      call ufree(ip_xdmold)
      call ufree(ip_idmold)
      call ufree(ip_nseg1)
      call ufree(ip_nsilo)
      call ufree(ip_nsihi)
      call ufree(ip_nsjlo)
      call ufree(ip_nsjhi)
      call ufree(ip_nsklo)
      call ufree(ip_nskhi)
      call ufree(ip_nb1)
      call ufree(ip_ne1)
      call ufree(ip_nb2)
      call ufree(ip_ne2)
      call ufree(ip_ibct)
      call ufree(ip_bcval)
      call ufree(ip_iovrlp)
      call ufree(ip_ifoflg)
      call ufree(ip_ndat)
      call ufree(ip_xmap)
      call ufree(ip_twotref)
      call ufree(ip_itrb1)
      call ufree(ip_itrb2)
      call ufree(ip_jtrb1)
      call ufree(ip_jtrb2)
      call ufree(ip_ktrb1)
      call ufree(ip_ktrb2)
      call ufree(ip_iturbb)
      call ufree(ip_ibif1)
      call ufree(ip_ibif2)
      call ufree(ip_nbi1)
      call ufree(ip_nei1)
      call ufree(ip_nbj1)
      call ufree(ip_nej1)
      call ufree(ip_nbk1)
      call ufree(ip_nek1)
      call ufree(ip_nbi2)
      call ufree(ip_nei2)
      call ufree(ip_nbj2)
      call ufree(ip_nej2)
      call ufree(ip_nbk2)
      call ufree(ip_nek2)
      call ufree(ip_nd11)
      call ufree(ip_nd21)
      call ufree(ip_nd12)
      call ufree(ip_nd22)
      call ufree(ip_ifsor)
      call ufree(ip_nb1s)
      call ufree(ip_ne1s)
      call ufree(ip_nb2s)
      call ufree(ip_ne2s)
      call ufree(ip_nsubbl)
      call ufree(ip_idobl)
      call ufree(ip_nseg)
      call ufree(ip_idno)
      call ufree(ip_ijk)
      call ufree(ip_idseg)
      call ufree(ip_idnext)
c
      write(*,'("SPLITTING/CONVERSION COMPLETED")')
c
      return
      end
