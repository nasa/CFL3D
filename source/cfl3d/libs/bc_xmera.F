      subroutine bc_xmera(ntime,nbl,lw,lw2,w,mgwk,wk,nwork,maxbl,iitot,
     .                    iviscg,iovrlp,lbg,ibpntsg,qb,iibg,kkbg,
     .                    jjbg,ibcg,nou,bou,nbuf,ibufdim,int_updt)
c***********************************************************************
c      Purpose: Determine boundary data/conditions at edges of overset
c               grids.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension w(mgwk),wk(nwork),lw(65,maxbl),lw2(43,maxbl)
      dimension iviscg(maxbl,3),iovrlp(maxbl),lbg(maxbl),
     .          ibpntsg(maxbl,4),qb(iitot,5,3),iibg(iitot),
     .          kkbg(iitot),jjbg(iitot),ibcg(iitot)
c
      dimension nou(nbuf)
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1
c
c*****************************
c     chimera grid boundaries
c*****************************
c
      if (iovrlp(nbl).eq.1 .and. ntime.gt.0) then
          ldim = 5
          call avghole(w(lq),w(lblk),jdim,kdim,idim,nbl,ldim,
     .    int_updt)
          call xupdt(w(lq),w(lqj0),w(lqk0),w(lqi0),jdim,kdim,idim,nbl,
     .    ldim,w(lbcj),w(lbck),w(lbci),maxbl,iitot,iibg,kkbg,jjbg,
     .    ibcg,lbg,ibpntsg,qb,nou,bou,nbuf,ibufdim,int_updt)
c
c         update turbulence quantities
c
         if (iviscg(nbl,1).ge.2 .or. iviscg(nbl,2).ge.2 .or.
     .       iviscg(nbl,3).ge.2) then
             ldim = 1
             call avghole(w(lvis),w(lblk),jdim,kdim,idim,nbl,ldim,
     .       int_updt)
             call xupdt(w(lvis),w(lvj0),w(lvk0),w(lvi0),jdim,kdim,
     .       idim,nbl,ldim,w(lbcj),w(lbck),w(lbci),maxbl,iitot,iibg,
     .       kkbg,jjbg,ibcg,lbg,ibpntsg,qb,nou,bou,nbuf,ibufdim,
     .       int_updt)
         end if
         if (iviscg(nbl,1).ge.4 .or. iviscg(nbl,2).ge.4 .or.
     .       iviscg(nbl,3).ge.4) then
             ldim = 2
             call avghole(w(lxib),w(lblk),jdim,kdim,idim,nbl,ldim,
     .       int_updt)
             call xupdt(w(lxib),w(ltj0),w(ltk0),w(lti0),jdim,kdim,
     .       idim,nbl,ldim,w(lbcj),w(lbck),w(lbci),maxbl,iitot,iibg,
     .       kkbg,jjbg,ibcg,lbg,ibpntsg,qb,nou,bou,nbuf,ibufdim,
     .       int_updt)
         end if
      end if
      return
      end
