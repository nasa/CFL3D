      subroutine lead(nbl,lw,lw2,maxbl)
c***********************************************************************
c     Purpose:  Install in common blocks the attributes of a block 
c     given a block number.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension lw(65,maxbl),lw2(43,maxbl)
c
      common /degshf/ ideg(3)
      common /fvfds/ rkap0(3),ifds(3)
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /lam/ ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi
      common /reyue/ reue,tinf,ivisc(3)
      common /wallfun/ iwf(3)
c
      jdim      = lw2(1, nbl) 
      kdim      = lw2(2, nbl)
      idim      = lw2(3, nbl)
      nblc      = lw2(4, nbl)
      jj2       = lw2(5, nbl)
      kk2       = lw2(6, nbl)
      ii2       = lw2(7, nbl)
      ilamlo    = lw2(8, nbl)
      ilamhi    = lw2(9, nbl)
      jlamlo    = lw2(10,nbl)
      jlamhi    = lw2(11,nbl)
      klamlo    = lw2(12,nbl)
      klamhi    = lw2(13,nbl)
      ivisc(1)  = lw2(14,nbl)
      ivisc(2)  = lw2(15,nbl)
      ivisc(3)  = lw2(16,nbl)
      ideg(1)   = lw2(17,nbl)
      ideg(2)   = lw2(18,nbl)
      ideg(3)   = lw2(19,nbl)
      iwf(1)    = lw2(20,nbl)
      iwf(2)    = lw2(21,nbl)
      iwf(3)    = lw2(22,nbl)
      idiag(1)  = lw2(23,nbl)
      idiag(2)  = lw2(24,nbl)
      idiag(3)  = lw2(25,nbl)
      iflim(1)  = lw2(26,nbl)
      iflim(2)  = lw2(27,nbl)
      iflim(3)  = lw2(28,nbl)
      ifds(1)   = lw2(29,nbl)
      ifds(2)   = lw2(30,nbl)
      ifds(3)   = lw2(31,nbl)
c     factor dividing rkap0g must be consistant with subroutine pointers!
      rkap0(1)  = lw2(32,nbl)/1e6
      rkap0(2)  = lw2(33,nbl)/1e6
      rkap0(3)  = lw2(34,nbl)/1e6
      js        = lw2(35,nbl)
      ks        = lw2(36,nbl)
      is        = lw2(37,nbl)
      je        = lw2(38,nbl)
      ke        = lw2(39,nbl)
      ie        = lw2(40,nbl)
c
      lq     = lw(1,nbl)
      lqj0   = lw(2,nbl)
      lqk0   = lw(3,nbl)
      lqi0   = lw(4,nbl)
      lsj    = lw(5,nbl)
      lsk    = lw(6,nbl)
      lsi    = lw(7,nbl)
      lvol   = lw(8,nbl)
      ldtj   = lw(9,nbl)
      lx     = lw(10,nbl)
      ly     = lw(11,nbl)
      lz     = lw(12,nbl)
      lvis   = lw(13,nbl)
      lsnk0  = lw(14,nbl)
      lsni0  = lw(15,nbl)
      lq1    = lw(16,nbl)
      lqr    = lw(17,nbl)
      lblk   = lw(18,nbl)
      lxib   = lw(19,nbl)
      lsig   = lw(20,nbl)
      lsqtq  = lw(21,nbl)
      lg     = lw(22,nbl)
      ltj0   = lw(23,nbl)
      ltk0   = lw(24,nbl)
      lti0   = lw(25,nbl)
      lxkb   = lw(26,nbl)
      lnbl   = lw(27,nbl)
      lvj0   = lw(28,nbl)
      lvk0   = lw(29,nbl)
      lvi0   = lw(30,nbl)
      lbcj   = lw(31,nbl)
      lbck   = lw(32,nbl)
      lbci   = lw(33,nbl)
      lqc0   = lw(34,nbl)
      ldqc0  = lw(35,nbl)
      lxtbj  = lw(36,nbl)
      lxtbk  = lw(37,nbl)
      lxtbi  = lw(38,nbl)
      latbj  = lw(39,nbl)
      latbk  = lw(40,nbl)
      latbi  = lw(41,nbl)
      lbcdj  = lw(42,nbl)
      lbcdk  = lw(43,nbl)
      lbcdi  = lw(44,nbl)
      lxib2  = lw(45,nbl)
      lux    = lw(46,nbl)
      lcmuv  = lw(47,nbl)
      lvolj0 = lw(48,nbl)
      lvolk0 = lw(49,nbl)
      lvoli0 = lw(50,nbl)
      lxmdj  = lw(51,nbl)
      lxmdk  = lw(52,nbl)
      lxmdi  = lw(53,nbl)
      lvelg  = lw(54,nbl)
      lxnm2  = lw(55,nbl)
      lynm2  = lw(56,nbl)
      lznm2  = lw(57,nbl)
      ldeltj = lw(58,nbl)
      ldeltk = lw(59,nbl)
      ldelti = lw(60,nbl)
      lxnm1  = lw(61,nbl)
      lynm1  = lw(62,nbl)
      lznm1  = lw(63,nbl)
c
      return
      end
