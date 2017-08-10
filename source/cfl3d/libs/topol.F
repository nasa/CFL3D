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
      subroutine topol(jdim,kdim,nsub,jjmax,kkmax,lmax,l,x,y,z,
     .                 xmid,ymid,zmid,xmide,ymide,zmide,limit,xc,yc,zc,
     .                 xie,eta,jimage,kimage,ifit,itmax,igap,iok,lout,
     .                 ic0,itoss0,jto,kto,iself,xif1,xif2,etf1,etf2,
     .                 nou,bou,nbuf,ibufdim,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Search appropriate "from" blocks for current "to" cell
c     center with coordinates xc,yc,zc; determine xie,eta of cell center
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension x(jdim,kdim,nsub),y(jdim,kdim,nsub),z(jdim,kdim,nsub)
      dimension xmid(jdim,kdim,nsub),ymid(jdim,kdim,nsub),
     .          zmid(jdim,kdim,nsub)
      dimension xmide(jdim,kdim,nsub),ymide(jdim,kdim,nsub),
     .          zmide(jdim,kdim,nsub)
      dimension jimage(nsub,jdim,kdim),kimage(nsub,jdim,kdim)
      dimension jjmax(nsub),kkmax(nsub),jfroz(itmax),kfroz(itmax)
      integer   lout(nsub),xif1(nsub),xif2(nsub),etf1(nsub),
     .          etf2(nsub)
c
      common /tol/ epsc,epsc0,epsreen,epscoll
c
      idum1 = 0
      idum2 = 0
      idum3 = 0
      idum4 = 0
      dum1  = 0.
      dum2  = 0.
      dum3  = 0.
c
      do 5 ll=1,lmax
  5   lout(ll) = 0
      iatmpt   = 0
      ichk     = 0
      ifroze   = 0
      ifroz    = 0
      ihuge    = 0
      jp       = xie
      kp       = eta
      lsav     = 1
      xiesav   = xie
      etasav   = eta
      jpsav    = jp
      kpsav    = kp
      jpc      = jp
      kpc      = kp
c
c     to start search, use solution from last cell, if one exists; otherwise
c     start by searching for minimum distance point over all "from" blocks 
c
      if (real(xie).lt.1. .or. real(eta).lt.1.) then
         call dsmin(jdim,kdim,nsub,jjmax,
     .        kkmax,lmax,x,y,z,xc,yc,zc,jp,kp,l,lout,999,
     .        xif1,xif2,etf1,etf2)
      end if
c
  999 continue
c
      iatmpt = iatmpt + 1
c
c     all "from" blocks have been searched if iatmp > lmax
c
      if(iatmpt.gt.lmax) go to 1000
c
      jmax = jjmax(l)
      kmax = kkmax(l)
      js   = xif1(l)
      je   = xif2(l)
      ks   = etf1(l)
      ke   = etf2(l)
c
c     avoid starting outside specified search range
c
      jp = min( jp , je-1)
      kp = min( kp,  ke-1)
      jp = max( jp , js)
      kp = max( kp , ks)
c
      do 5555 intern=1,itmax
c
      jfroz(intern) = jp
      kfroz(intern) = kp
c     call trace(3,intern,idum2,idum3,idum4,dum1,dum2,dum3)
c
c     find local xie, eta via Newton iteraton in current target cell jp,kp
c
      if (itoss0 .eq. 0 ) then
c
c        call general routine which first determines the best direction 
c        for inversion, then call Newton iteration routine
c
         call xe(jdim,kdim,nsub,l,x,y,z,xmid,ymid,zmid,xmide,ymide,
     .           zmide,jp,kp,xc,yc,zc,xie,eta,imiss,ifit,ic0,
     .           nou,bou,nbuf,ibufdim,myid)
      else
c
c        best direction for inversion is known apriori (i.e. interface
c        lies on a constant coordinate surface, or nearly so)
c
         x1 = x(jp,kp,l)
         y1 = y(jp,kp,l)
         z1 = z(jp,kp,l)
         x2 = x(jp+1,kp,l)
         y2 = y(jp+1,kp,l)
         z2 = z(jp+1,kp,l)
         x4 = x(jp,kp+1,l)
         y4 = y(jp,kp+1,l)
         z4 = z(jp,kp+1,l)
         x3 = x(jp+1,kp+1,l)
         y3 = y(jp+1,kp+1,l)
         z3 = z(jp+1,kp+1,l)
         x5 = xmid(jp,kp,l)
         y5 = ymid(jp,kp,l)
         z5 = zmid(jp,kp,l)
         x6 = xmid(jp,kp+1,l)
         y6 = ymid(jp,kp+1,l)
         z6 = zmid(jp,kp+1,l)
         x7 = xmide(jp,kp,l)
         y7 = ymide(jp,kp,l)
         z7 = zmide(jp,kp,l)
         x8 = xmide(jp+1,kp,l)
         y8 = ymide(jp+1,kp,l)
         z8 = zmide(jp+1,kp,l)
c
         if (itoss0.eq.1) then
c           use only y and z equations
            call xe2(y1,y2,y3,y4,y5,y6,y7,y8,yc,z1,z2,z3,z4,z5,
     .               z6,z7,z8,zc,xie,eta,imiss,ifit,
     .               nou,bou,nbuf,ibufdim,myid)
c     call trace(9,l,jp,kp,idum4,dum1,dum2,dum3)
         end if
         if (itoss0.eq.2) then
c           use only x and z equations
            call xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,z1,z2,z3,z4,z5,
     .               z6,z7,z8,zc,xie,eta,imiss,ifit,
     .               nou,bou,nbuf,ibufdim,myid)
c     call trace(10,l,jp,kp,idum4,dum1,dum2,dum3)
         end if
         if (itoss0.eq.3) then
c           use only x and y equations
            call xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,y1,y2,y3,y4,y5,
     .              y6,y7,y8,yc,xie,eta,imiss,ifit,
     .              nou,bou,nbuf,ibufdim,myid)
c     call trace(11,l,jp,kp,idum4,dum1,dum2,dum3)
         end if
      end if
c
c     call trace(4,idum1,idum2,idum3,idum4,xie,eta,dum3)
c
c     check to make sure that search did not find current "to" cell as the
c     target (a possibility for a grid which communicates with itself along
c     a branch cut). If so, move to to a cell on the other side of the
c     branch cut
c
      if (iself.gt.0) then
         if (imiss.eq.0 .and. (jto.eq.jp-1 .and. kto.eq.kp-1)) then
c
c     call trace(51,jto,kto,idum3,idum4,dum1,dum2,dum3)
c
c           determine in which direction the branch cut lies
c
            jpc = jpc - 1
            jpc = max( 1 , jpc )
            if (jpc .ne. jimage(l,jpc,kp)) ibrdir = 1
            jpc = jpc + 1
            jpc = min( jpc, jmax-1 )
            if (jpc .ne. jimage(l,jpc,kp)) ibrdir = 2
            kpc = kpc - 1
            kpc = max( 1 , kpc )
            if (kpc .ne. kimage(l,jp,kpc)) ibrdir = 3
            kpc = kpc + 1
            kpc = min( kpc, kmax-1 )
            if (kpc .ne. kimage(l,jp,kpc)) ibrdir = 4
c
c           set cell index to trigger a branch cut jump
c
            if(ibrdir .eq. 1) jp = jp - 1
            if(ibrdir .eq. 2) jp = jp + 1
            if(ibrdir .eq. 3) kp = kp - 1
            if(ibrdir .eq. 4) kp = kp + 1
c
            jp = min( jp , jmax-1 )
            kp = min( kp , kmax-1 )
            jp = max( 1 , jp )
            kp = max( 1 , kp )
            jp = jimage(l,jp,kp)
            kp = kimage(l,jp,kp)
            imiss = 1
            go to 5555
         end if
      end if
c
c     current target cell correct if imiss = 0
c
      if (imiss.eq.0) go to 5556
c
c     update current guess for target cell based on result of Newton
c     iteration, with max allowable change set by limit
c
      if (real(xie).ge.0) jinc = abs(xie)
      if (real(xie).lt.0) jinc = abs(xie-1)
      if (real(eta).ge.0) kinc = abs(eta)
      if (real(eta).lt.0) kinc = abs(eta-1)
c
      jinc = min( jinc , limit )
      kinc = min( kinc , limit )
c
      if (real(xie).gt.1.0) then
         jp = jp + jinc  
      else if (real(xie).lt.0.) then
         jp = jp - jinc  
      end if
      if (real(eta).gt.1.0) then
         kp = kp + kinc
      else if (real(eta).lt.0.) then
         kp = kp - kinc 
      end if
c
c     keep within bounds of (expanded) "from" block
c
      jp = min( jp , jmax-1 )
      kp = min( kp , kmax-1 )
      jp = max( 1 , jp )
      kp = max( 1 , kp )
c
      jmax = jjmax(l)
      kmax = kkmax(l)
      js   = xif1(l)
      je   = xif2(l)
      ks   = etf1(l)
      ke   = etf2(l)
c
c     avoid cells outside specified search range
c
      jp = min( jp , je-1)
      kp = min( kp,  ke-1)
      jp = max( jp , js)
      kp = max( kp , ks)
c
c     account for any branch cuts 
c
      jpc = jimage(l,jp,kp)
      kpc = kimage(l,jp,kp)
c     call trace(99,jp,kp,jpc,kpc,dum1,dum2,dum3)
      jp = jpc
      kp = kpc
c
c     search routine off track if local xie or eta become huge...
c     try to get back on track via minimum distance search in the
c     current "from" block
c
      huge=1.e+5
      if (abs(real(xie)).gt.real(huge).or.abs(real(eta)).gt.real(huge))
     .   then
         ihuge = ihuge + 1
         if (ihuge.gt.1) go to 1000
c     call trace(41,idum1,idum2,idum3,idum4,dum1,dum2,dum3)
         call dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,z,xc,yc,zc,
     .              jp,kp,l,lout,-999,xif1,xif2,etf1,etf2)
         go to 5555
      end if
c
c     check for frozen convergence: search routine keeps returning
c     to the same point, without 0 < xie,eta < 1 at that point. if frozen, 
c     attempt to break out of cycle by using minimum distance search
c     in the current "from" block
c
      ifroz = 0
      do 77 ii=1,intern
      int = intern-ii+1
      if (jp.eq.jfroz(int).and.kp.eq.kfroz(int)) ifroz  = 1
77    continue
      if (ifroz.eq.1) then
         ifroze = ifroze + 1
         if (ifroze.gt.1) go to 1000
c     call trace(42,jp,kp,l,idum4,dum1,dum2,dum3)
         call dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,z,xc,yc,zc,
     .        jp,kp,l,lout,-999,xif1,xif2,etf1,etf2)
      end if
c
5555  continue
c
1000  continue
c
c     search routine has been unsuccessful in the current "from" block; if
c     all "from" blocks have not been searched, try another; if all "from" 
c     blocks have been searched, use any previously obtained generalized 
c     coordinates which were found in an expanded cell (ichk > 0)  
c     otherwise, exit to try a different basis function and start over
c
      if (iatmpt.lt.lmax) then
c
c        find new starting cell by searching for minimum distance over 
c        all "from" blocks not yet searched
c
         lout(l) = 1
c     call trace(5,l,idum2,idum3,idum4,dum1,dum2,dum3)
          call dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,z,xc,yc,zc,
     .               jp,kp,l,lout,999,xif1,xif2,etf1,etf2)
         ifroze  = 0
         ihuge = 0
         go to 999
      else
         if (ichk.eq.0) then
           iok = 0
           return
         else
           l   = lsav
           jp  = jpsav
           kp  = kpsav
           xie = xiesav
           eta = etasav
         end if
         if (ichk.eq.1) then
c
c          point really does lie in the expanded cell of the "from" block 
c          in which it was first found
c
c     call trace(7,l,idum2,idum3,idum4,dum1,dum2,dum3)
           go to 5557
         end if
         if (ichk.gt.1) then
c
c          point lies in the expanded zone of two or more grids, 
c          possibly due to small gap between blocks on "from" side.
c          flag this case for later use in diagnostics routine
           igap = 1    
           go to 5557
         end if
      end if
5556  continue
c
c     search routine has been successful; however, if there are multiple 
c     blocks on the "from" side of patch surface, check to see if:
c     1) for interfaces flagged as C-0, the cell center coordinates of the
c     "from" cell jp,kp match the coordinates of the "to" cell center to
c     insure the right "from" block has been used;
c     2)the "to" cell center lies in the overlap region formed when the
c     "from" blocks are expanded. in such cases the cell center may have
c     been found in the wrong "from" block; alternatively, a gap may exist
c     between blocks
c
      if (lmax.gt.1) then
c
        if(ic0.gt.0) then
          j = jp
          k = kp
          x1c = 0.25*( x(j,k,l) + x(j+1,k,l)
     .          + x(j+1,k+1,l) + x(j,k+1,l) )
          y1c = 0.25*( y(j,k,l) + y(j+1,k,l)
     .          + y(j+1,k+1,l) + y(j,k+1,l) )
          z1c = 0.25*( z(j,k,l) + z(j+1,k,l)
     .          + z(j+1,k+1,l) + z(j,k+1,l) )
          if(abs(real(x1c-xc)).gt.real(epsc0) .or.
     .       abs(real(y1c-yc)).gt.real(epsc0) .or.
     .       abs(real(z1c-zc)).gt.real(epsc0)) then
          ichk = ichk+1
          if (ichk.eq.1) then
             lsav    = l
             jpsav   = jp
             kpsav   = kp
             xiesav  = xie
             etasav  = eta
          end if
c
c           find alternate block which may also contain the "to"
c           cell center by using a minimum distance search in
c           all "from" blocks not yet searched
c
            lout(l) = 1
c     call trace(6,l,idum2,idum3,idum4,dum1,dum2,dum3)
            call dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,
     .           z,xc,yc,zc,jp,kp,l,lout,999,xif1,xif2,etf1,etf2)
            ifroze = 0
            ihuge = 0
            go to 999
          end if
        end if
c
        if (jp.eq.1.or.jp.eq.jmax-1.or.kp.eq.1.or.kp.eq.kmax-1) then
          ichk = ichk+1
          if (ichk.eq.1) then
             lsav    = l
             jpsav   = jp
             kpsav   = kp
             xiesav  = xie
             etasav  = eta
          end if
c
c         find alternate block which may also contain the "to"
c         cell center by using a minimum distance search in 
c         all "from" blocks not yet searched
c
          lout(l) = 1
c     call trace(6,l,idum2,idum3,idum4,dum1,dum2,dum3)
          call dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,
     .         z,xc,yc,zc,jp,kp,l,lout,999,xif1,xif2,etf1,etf2)
          ifroze = 0
          ihuge = 0
          go to 999
        end if
      end if
c
5557  continue
c
c     completed search for "to" cell center; convert to global xie, eta
c
      iok = 1
      xie = xie + jp
      eta = eta + kp
c
      return
      end
