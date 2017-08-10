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
      subroutine bc_xmera(ntime,nbl,lw,lw2,w,mgwk,wk,nwork,maxbl,iitot,
     .                    iviscg,iovrlp,lbg,ibpntsg,qb,iibg,kkbg,
     .                    jjbg,ibcg,nou,bou,nbuf,ibufdim,int_updt,
     .                    nummem)
c
c     $Id$
c
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
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
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
             ldim = nummem
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
