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
      subroutine diagnos(istop,iout,igap,jdim1,kdim1,msub1,
     .                  msub2,jjmax1,kkmax1,lmax1,x1,y1,z1,x1mid,y1mid,
     .                  z1mid,x1mide,y1mide,z1mide,x2int,y2int,z2int,
     .                  x2fit,y2fit,z2fit,jjmax2,kkmax2,x2,y2,z2,
     .                  xie2,eta2,mblkpt,icheck,intmx,xc,yc,zc,ifit,
     .                  j21,j22,k21,k22,npt,ic0,iorph,xif1,xif2,
     .                  etf1,etf2,itoss0,iself,
     .                  nou,bou,nbuf,ibufdim,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Perform diagnostic checks on interpolation coefficients
c     (generalized coordinates) found via search and inversion routines
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension x1(jdim1,kdim1,msub1),y1(jdim1,kdim1,msub1),
     .          z1(jdim1,kdim1,msub1)
      dimension x1mid(jdim1,kdim1,msub1),y1mid(jdim1,kdim1,msub1),
     .          z1mid(jdim1,kdim1,msub1)
      dimension x1mide(jdim1,kdim1,msub1),y1mide(jdim1,kdim1,msub1),
     .          z1mide(jdim1,kdim1,msub1)
      dimension x2fit(jdim1,kdim1,msub2),y2fit(jdim1,kdim1),
     .          z2fit(jdim1,kdim1)
      dimension x2int(jdim1,kdim1,msub1),y2int(jdim1,kdim1,msub1),
     .          z2int(jdim1,kdim1,msub2)
      dimension x2(jdim1,kdim1,msub2),y2(jdim1,kdim1,msub2),
     .          z2(jdim1,kdim1,msub2)
      dimension xie2(npt),eta2(npt),mblkpt(npt)
      dimension jjmax1(msub1),kkmax1(msub1),jjmax2(msub2),kkmax2(msub2)
      integer   xif1(msub1),xif2(msub1),etf1(msub1),
     .          etf2(msub1)
c
      character*14 titlptchgrd
c
      common /igrdtyp/ ip3dgrd,ialph
      common /tol/ epsc,epsc0,epsreen,epscoll
c
      jmax2 = jjmax2(1) - 1
      kmax2 = kkmax2(1) - 1
c
c     diagnostics for case where search/inversion routines have failed
c
      if(istop.eq.1) then
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),1222)
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),1224)
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),1226)
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),1227)
        go to 552
      end if
c
c     diagnostics for completed generalized coordinate interpolation
c
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),*)'    beginning diagnostic checks'
c
c     ****** first check ******
c
c     check final cell center generalized coordinates to be sure that xie,eta
c     of the "to" grid lies inside the domain of the from grid(s)
c
      neta = 0
      nxie = 0
c
      do 617 k=k21,k22-1
      do 617 j=j21,j22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      l1=mblkpt(ll) 
      if (l1.gt.0) then
      if(real(eta2(ll)).lt.1. .or. real(eta2(ll)).gt.kkmax1(l1)-2) then
        nou(2) = min(nou(2)+1,ibufdim)
        write(bou(nou(2),2),'(''interp. no.: '',i4,'', "to" cell '',
     .             ''center j,k = '',2i4,'' found in "from" block '',
     .             i4)') icheck,j,k,l1
        nou(2) = min(nou(2)+1,ibufdim)
        write(bou(nou(2),2),'(''   with eta = '',e14.7)') 
     .        real(eta2(ll))
        nou(2) = min(nou(2)+1,ibufdim)
        write(bou(nou(2),2),'(''   in this "from" block, legal range'',
     .             '' of eta is 1 to '',i4)') kkmax1(l1)-2
        neta = neta + 1
      end if
      end if
617   continue
c
      do 618 j=j21,j22-1
      do 618 k=k21,k22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      l1=mblkpt(ll) 
      if (l1.gt.0) then
      if(real(xie2(ll)).lt.1. .or. real(xie2(ll)).gt.jjmax1(l1)-2)then
        nou(2) = min(nou(2)+1,ibufdim)
        write(bou(nou(2),2),'(''interp. no.: '',i4,'', "to" cell'',
     .             '' center j,k = '',2i4,'' found in "from" block '',
     .             i4)') icheck,j,k,l1
        nou(2) = min(nou(2)+1,ibufdim)
        write(bou(nou(2),2),'(''   with xie = '',e14.7)') 
     .        real(xie2(ll))
        nou(2) = min(nou(2)+1,ibufdim)
        write(bou(nou(2),2),'(''   in this "from" block, legal range'',
     .             '' of xie is 1 to '',i4)') jjmax1(l1)-2
        nxie = nxie + 1
      end if
      end if
618   continue
c
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),710) nxie
710   format('     ',i4,' points were found outside the legal',
     .' range of xie')
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),711) neta
711   format('     ',i4,' points were found outside the legal',
     .' range of eta')
      if(neta+nxie .gt. 0) then
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),*)'        check fort.9 to make sure',
     .           ' values are not far outside legal range'
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),1222)
          if(igap.gt.0) then
            nou(4) = min(nou(4)+1,ibufdim)
            write(bou(nou(4),4),1225)
          else
            if (nxie.gt.0) then
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),1228)
            end if
            if (neta.gt.0) then
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),1229)
            end if
          end if
          if (iorph.gt.0) then
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),1332)
          end if
        iout = 1
      end if
c
c     ****** second check (grids with orphan points only) ******
c
c     count number of orphans
c
      if(iorph.gt.0) then
        norph = 0
        do 1500 j=j21,j22-1
        do 1500 k=k21,k22-1
        ll = (j22-j21)*(k-k21) + (j-j21+1)
        l1   = mblkpt(ll)
        if(l1.eq.0) norph = norph + 1
1500    continue
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),721)norph
721     format('     ',i4,' points were flagged as orphans')
      end if
c
c     ****** third check (grids that self-connect with a branch cut ) ******
c
c     check for points that get interpolated from themselves
c     (branch cut error)
c
      if(iself.gt.0) then
        ierr = 0
        do 1699 j=j21,j22-1
        do 1699 k=k21,k22-1
        ll = (j22-j21)*(k-k21) + (j-j21+1)
        l1   = mblkpt(ll)
        jc   = int( xie2(ll) )
        kc   = int( eta2(ll) )
        if(jc.eq.j .and. kc.eq.k) then
          ierr = ierr + 1
          nou(2) = min(nou(2)+1,ibufdim)
          write(bou(nou(2),2),'(''branch cut problem at j,k ='',2i4)')
     .    j,k
        end if
1699    continue
        if (ierr.gt.0) then
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),821) ierr
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),822)
        end if
821     format('     ',i4,' points are interpolated from themselves')
822     format(
     .  '          search routine failure due to branch cut - DO NOT ',
     .  'run flow solver')
      end if
c
c     ****** fourth check (C-0 continuous grids only) ******
c
c     for C-0 grids, check to make sure all interpolation coefficients
c     are xie=m+0.5, eta=n+0.5, to within the tolerence set on convergence
c     for the generalized coordinates
c
      if (ic0.gt.0) then
        c0tol  = epsc
        nc00   = (j22-j21)*(k22-k21)
        nc0    = nc00
        errmx  = 0.
        jerrmx = 1
        kerrmx = 1
c
        do 1600 j=j21,j22-1
        do 1600 k=k21,k22-1
        ll = (j22-j21)*(k-k21) + (j-j21+1)
        l1   = mblkpt(ll)
        jc   = int( xie2(ll) )
        kc   = int( eta2(ll) )
        xiec = ccabs(xie2(ll)-jc-.5)
        etac = ccabs(eta2(ll)-kc-.5)
        if (real(xiec).gt.real(errmx)) then
           jerrmx = j
           kerrmx = k
           errmx  = xiec
        end if
        if (real(etac).gt.real(errmx)) then
           jerrmx = j
           kerrmx = k
           errmx  = etac
        end if
        if (real(xiec).gt.real(c0tol).or.real(etac).gt.real(c0tol)) then
           nc0 = nc0 - 1
        end if
1600    continue
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),722)nc00-nc0,real(c0tol)
722     format('     ',i4,' points  ',
     .  'flagged as C-0 have |xie-.5| or |eta-.5| >',e10.3)
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),723) real(errmx),jerrmx,kerrmx
723     format('          the max. error ',e10.3,
     .  ' occurs at j,k = ',i4,',',i4)
        if (nc00-nc0 .ne. 0) then
           iout = 1
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),1222)
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),1330)
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),1224)
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),1331)
         else
           do 1601 j=j21,j22-1
           do 1601 k=k21,k22-1
           ll = (j22-j21)*(k-k21) + (j-j21+1)
           jc   = int( xie2(ll) )
           kc   = int( eta2(ll) )
           xie2(ll) = jc + 0.5
           eta2(ll) = kc + 0.5
1601       continue
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),720)
720        format('          C-0 check acceptable, resetting xie ',
     .     'and eta to exact values (0.5,0.5)')
        end if
      end if
c
c     ****** fifth check ******
c
c     check jacobians of grid in generalized coordinates for anomolies 
c     (negative areas)
c
c     not applicable in 2D cases or cases with orhan points!
c
      jseg = j22 - j21
      kseg = k22 - k21
      if(jseg.gt.1 .and. kseg.gt.1 .and. iorph.le.0)then
c
      nct   = 0
      iflg  = 0
      do 1750 l=1,lmax1
c
c     find j11,k11, the first "to" cell located in "from" block l.
c     the sign of the jacobians for all other "to" cells in block l will
c     be compared to the sign in cell j11,k11. Note: need two consecutive
c     values of xie and eta to be in block l for this check to be appropriate.
c     Also skip test at points that, while adjacent on the "to" side, recieve
c     data from opposite sides of a branch cut on the "from" side - such
c     points are incorrectly flagged as problems by this test
c
      j11 = -1
      k11 = -1
      do 1754 j=j21,j22-2
      do 1755 k=k21,k22-2
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      ll2 = (j22-j21)*(k-k21) + (j+1-j21+1)
      ll3 = (j22-j21)*(k+1-k21) + (j+1-j21+1)
      ll4 = (j22-j21)*(k+1-k21) + (j-j21+1)
      if (mblkpt(ll).eq.l .and. mblkpt(ll2).eq.l .and.
     .   mblkpt(ll3).eq.l .and. mblkpt(ll4).eq.l) then
c        assume neighboring target cells lie on opposite sides of a branch
c        if generalized coordinates differ by more than 1/2 grid dimensions
         jtest = jjmax1(mblkpt(ll))/2
         ktest = kkmax1(mblkpt(ll))/2
         jdif24 = abs(int(xie2(ll2) - xie2(ll4)))
         kdif24 = abs(int(eta2(ll2) - eta2(ll4)))
         jdif31 = abs(int(xie2(ll3) - xie2(ll)))
         kdif31 = abs(int(eta2(ll3) - eta2(ll)))
         ibrnch = 0
         if (jdif24.gt.jtest .or. kdif24.gt.ktest 
     .      .or. jdif31.gt.jtest .or. kdif31.gt.ktest) ibrnch = 1
         if (ibrnch .eq. 0) then
            j11 = j
            k11 = k
            go to 1756
         end if
      end if
1755  continue  
1754  continue
1756  continue
      if(j11.lt.1 .or. k11.lt.1) go to 1750
      iflg = 1
      j     = j11
      k     = k11
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      ll2 = (j22-j21)*(k-k21) + (j+1-j21+1)
      ll3 = (j22-j21)*(k+1-k21) + (j+1-j21+1)
      ll4 = (j22-j21)*(k+1-k21) + (j-j21+1)
      axie  = xie2(ll3) - xie2(ll)
      aeta  = eta2(ll3) - eta2(ll)
      bxie  = xie2(ll4) - xie2(ll2)
      beta  = eta2(ll4) - eta2(ll2)
      bj0   = axie*beta-bxie*aeta
c
      do 1764 j=j21,j22-2
      do 1765 k=k21,k22-2
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      ll2 = (j22-j21)*(k-k21) + (j+1-j21+1)
      ll3 = (j22-j21)*(k+1-k21) + (j+1-j21+1)
      ll4 = (j22-j21)*(k+1-k21) + (j-j21+1)
      if (mblkpt(ll).eq.l .and. mblkpt(ll2).eq.l .and.
     .   mblkpt(ll3).eq.l .and. mblkpt(ll4).eq.l) then
c        assume neighboring target cells lie on opposite sides of a branch
c        if generalized coordinates differ by more than 1/2 grid dimensions
         jtest = jjmax1(mblkpt(ll))/2
         ktest = kkmax1(mblkpt(ll))/2
         jdif24 = abs(int(xie2(ll2) - xie2(ll4)))
         kdif24 = abs(int(eta2(ll2) - eta2(ll4)))
         jdif31 = abs(int(xie2(ll3) - xie2(ll)))
         kdif31 = abs(int(eta2(ll3) - eta2(ll)))
         ibrnch = 0
         if (jdif24.gt.jtest .or. kdif24.gt.ktest
     .      .or. jdif31.gt.jtest .or. kdif31.gt.ktest) ibrnch = 1
         if (ibrnch .eq. 0) then
            axie  = xie2(ll3) - xie2(ll)
            aeta  = eta2(ll3) - eta2(ll)
            bxie  = xie2(ll4) - xie2(ll2)
            beta  = eta2(ll4) - eta2(ll2)
            bj    = axie*beta-bxie*aeta
            if (real(bj)*real(bj0) .le. 0.) then
               nct = nct + 1
               nou(2) = min(nou(2)+1,ibufdim)
               write(bou(nou(2),2),'(''non-unique point at j,k = '',
     .         i5,'','',i5)') j,k
               nou(2) = min(nou(2)+1,ibufdim)
               write(bou(nou(2),2),'(''xie2(j,k),eta2(j,k)         = '',
     .         f11.5,'','',f11.5)') real(xie2(ll)),real(eta2(ll))
               nou(2) = min(nou(2)+1,ibufdim)
               write(bou(nou(2),2),'(''xie2(j+1,k),eta2(j+1,k)     = '',
     .         f11.5,'','',f11.5)') real(xie2(ll2)),real(eta2(ll2))
               nou(2) = min(nou(2)+1,ibufdim)
               write(bou(nou(2),2),'(''xie2(j+1,k+1),eta2(j+1,k+1) = '',
     .         f11.5,'','',f11.5)') real(xie2(ll3)),real(eta2(ll3))
               nou(2) = min(nou(2)+1,ibufdim)
               write(bou(nou(2),2),'(''xie2(j,k+1),eta2(j,k+1)     = '',
     .         f11.5,'','',f11.5)') real(xie2(ll4)),real(eta2(ll4))
            end if
         end if
      end if
1765  continue
1764  continue
c
1750  continue
c
      if(iflg.eq.0) then
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),*)'  caution...no consecutive "to" cells',
     .  ' found in any of the "from" blocks'
      end if
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),724)nct
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),725)
724   format('     ',i4,' non-unique points were found in the',
     .           ' generalized-coordinate ')
725   format(
     .'       mapping between the "to" and "from" grids')
      if(nct .gt. 0) then  
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),*)'        Caution: this may indicate a',
     .' serious problem'
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),1222)
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),1223)
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),1226)
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),1333)
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),1334)
      iout = 1
      end if
      end if
c
c     ****** sixth check ******
c
c     compute cell centers of "to" grid from the interpolation coefficients
c     and compare with a direct calculation from the "to" grid. both calculations
c     are done consistant with the input value of ifit.
c
      nbpt  = 0
      percnt = 0.10
      errmax = 0.
      jerrmx = 1
      kerrmx = 1
c
      do 9113 j=j21,j22-1
      do 9113 k=k21,k22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      l1   = mblkpt(ll)
c
c     skip if orphan point
c
      if(l1.eq.0) go to 9113
c
      jc   = int( xie2(ll) )
      kc   = int( eta2(ll) )
c
c     add 1 to to account for expanded "from" grids
c
      jc   = jc + 1
      kc   = kc + 1
      xiejk = xie2(ll) + 1
      etajk = eta2(ll) + 1
      xiec = float(jc)
      etac = float(kc)
c
      x11 =  x1(jc,kc,l1)      
      y11 =  y1(jc,kc,l1)      
      z11 =  z1(jc,kc,l1)      
      x12 =  x1(jc+1,kc,l1)      
      y12 =  y1(jc+1,kc,l1)      
      z12 =  z1(jc+1,kc,l1)      
      x13 =  x1(jc+1,kc+1,l1)      
      y13 =  y1(jc+1,kc+1,l1)      
      z13 =  z1(jc+1,kc+1,l1)      
      x14 =  x1(jc,kc+1,l1)      
      y14 =  y1(jc,kc+1,l1)      
      z14 =  z1(jc,kc+1,l1)
      x15 =  x1mid(jc,kc,l1)      
      y15 =  y1mid(jc,kc,l1)      
      z15 =  z1mid(jc,kc,l1)      
      x16 =  x1mid(jc,kc+1,l1)      
      y16 =  y1mid(jc,kc+1,l1)      
      z16 =  z1mid(jc,kc+1,l1)      
      x17 =  x1mide(jc,kc,l1)      
      y17 =  y1mide(jc,kc,l1)      
      z17 =  z1mide(jc,kc,l1)      
      x18 =  x1mide(jc+1,kc,l1)      
      y18 =  y1mide(jc+1,kc,l1)      
      z18 =  z1mide(jc+1,kc,l1)
c      
      dx2 = x12 - x11
      dy2 = y12 - y11
      dz2 = z12 - z11
      dx3 = x13 - x11
      dy3 = y13 - y11
      dz3 = z13 - z11
      dx4 = x14 - x11
      dy4 = y14 - y11
      dz4 = z14 - z11
      dx5 = x15 - x11
      dy5 = y15 - y11
      dz5 = z15 - z11
      dx6 = x16 - x11
      dy6 = y16 - y11
      dz6 = z16 - z11
      dx7 = x17 - x11
      dy7 = y17 - y11
      dz7 = z17 - z11
      dx8 = x18 - x11
      dy8 = y18 - y11
      dz8 = z18 - z11
c
      if (ifit.eq.1) then
c     bi-linear fit
      a2 = dx2
      a3 = dx4
      a4 = dx3 - a2 - a3
      a5 = 0.
      a6 = 0.
      a7 = 0.
      a8 = 0.
      b2 = dy2
      b3 = dy4
      b4 = dy3 - b2 - b3
      b5 = 0.
      b6 = 0.
      b7 = 0.
      b8 = 0.
      c2 = dz2
      c3 = dz4
      c4 = dz3 - c2 - c3
      c5 = 0.
      c6 = 0.
      c7 = 0.
      c8 = 0.
      end if
c
      if (ifit.eq.2) then
c     (degenerate) bi-quadradratic fit
c     (quadratic in xie and eta, but without the (xie**2)*(eta**2) term)
      a2  = -dx2 + 4.*dx5
      a3  = -dx4 + 4.*dx7
      a5  = 2.*dx2 - 4.*dx5
      a7  = 2.*dx4 - 4.*dx7
      df1 = dx3 - a2 - a3 - a5 - a7
      df2 = dx6 - .5*a2 - a3 - .25*a5 - a7
      df3 = dx8 - a2 - .5*a3 - a5 - .25*a7
      a4  = -3.*df1 + 4.*df2 + 4.*df3
      a6  = 2.*df1 - 4.*df2
      a8  = 2.*df1 - 4.*df3
      b2  = -dy2 + 4.*dy5
      b3  = -dy4 + 4.*dy7
      b5  = 2.*dy2 - 4.*dy5
      b7  = 2.*dy4 - 4.*dy7
      df1 = dy3 - b2 - b3 - b5 - b7
      df2 = dy6 - .5*b2 - b3 - .25*b5 - b7
      df3 = dy8 - b2 - .5*b3 - b5 - .25*b7
      b4  = -3.*df1 + 4.*df2 + 4.*df3
      b6  = 2.*df1 - 4.*df2
      b8  = 2.*df1 - 4.*df3
      c2  = -dz2 + 4.*dz5
      c3  = -dz4 + 4.*dz7
      c5  = 2.*dz2 - 4.*dz5
      c7  = 2.*dz4 - 4.*dz7
      df1 = dz3 - c2 - c3 - c5 - c7
      df2 = dz6 - .5*c2 - c3 - .25*c5 - c7
      df3 = dz8 - c2 - .5*c3 - c5 - .25*c7
      c4  = -3.*df1 + 4.*df2 + 4.*df3
      c6  = 2.*df1 - 4.*df2
      c8  = 2.*df1 - 4.*df3
      end if
c
c     quadratic fit in xie, linear fit in eta
      if (ifit.eq.3) then
      a3  = dx4
      b3  = dy4
      c3  = dz4
      a2  = -dx2 + 4.*dx5
      b2  = -dy2 + 4.*dy5
      c2  = -dz2 + 4.*dz5
      a5  = 2.*dx2  - 4.*dx5
      b5  = 2.*dy2  - 4.*dy5
      c5  = 2.*dz2  - 4.*dz5
      df1 = x13 - x12 - a3
      df2 = x16 - x15 - a3
      a4  =   -df1 + 4.*df2
      a6  = 2.*df1 - 4.*df2
      df1 = y13 - y12 - b3
      df2 = y16 - y15 - b3
      b4  =   -df1 + 4.*df2
      b6  = 2.*df1 - 4.*df2
      df1 = z13 - z12 - c3
      df2 = z16 - z15 - c3
      c4  =   -df1 + 4.*df2
      c6  = 2.*df1 - 4.*df2
      a7  = 0.
      a8  = 0.
      b7  = 0.
      b8  = 0.
      c7  = 0.
      c8  = 0.
      end if
c 
c     linear fit in xie,quadratic fit in eta
      if (ifit.eq.4) then
      a2  = dx2
      a3  = -dx4 + 4.*dx7
      a7  = 2.*dx4 -4.*dx7
      df1 = dx3 - a2 - a3 - a7
      df2 = dx8 - a2 -.5*a3 - .25*a7
      a4  = -df1 + 4.*df2
      a8  = 2.*df1 - 4.*df2
      a5  = 0.
      a6  = 0.
      b2  = dy2
      b3  = -dy4 + 4.*dy7
      b7  = 2.*dy4 -4.*dy7
      df1 = dy3 - b2 - b3 - b7
      df2 = dy8 - b2 -.5*b3 - .25*b7
      b4  = -df1 + 4.*df2
      b8  = 2.*df1 - 4.*df2
      b5  = 0.
      b6  = 0.
      c2  = dz2
      c3  = -dz4 + 4.*dz7
      c7  = 2.*dz4 -4.*dz7
      df1 = dz3 - c2 - c3 - c7
      df2 = dz8 - c2 -.5*c3 - .25*c7
      c4  = -df1 + 4.*df2
      c8  = 2.*df1 - 4.*df2
      c5  = 0.
      c6  = 0.
      end if
c
      xie = xiejk - xiec
      eta = etajk - etac
c
      x2int(j,k,1) = x11 + a3*eta + eta*( a7*eta + a8*xie*eta )
     .   + xie*( a2 + a4*eta + a5*xie + a6*xie*eta )  
      y2int(j,k,1) = y11 + b3*eta + eta*( b7*eta + b8*xie*eta )
     .   + xie*( b2 + b4*eta + b5*xie + b6*xie*eta )  
      z2int(j,k,1) = z11 + c3*eta + eta*( c7*eta + c8*xie*eta )
     .   + xie*( c2 + c4*eta + c5*xie + c6*xie*eta )  
c
c     compute center of "to" cell directly from the "to" grid, consistent 
c     with ifit
c
      jl = 1
      jr = jmax2
      kl = 1
      kr = kmax2
      kcall = k+1
      call extra(jdim1,kdim1,msub2,1,x2,y2,z2,
     .           j,kcall,jl,jr,x6,y6,z6,icase,ifit)
      call extra(jdim1,kdim1,msub2,1,x2,y2,z2,
     .           j,k,jl,jr,x5,y5,z5,icase,ifit)
      call extrae(jdim1,kdim1,msub2,1,x2,y2,z2,
     .            j,k,kl,kr,x7,y7,z7,icase,ifit)
      jcall = j+1
      call extrae(jdim1,kdim1,msub2,1,x2,y2,z2,
     .            jcall,k,kl,kr,x8,y8,z8,icase,ifit)
c
c     bi-linear
      if (ifit .eq. 1) then
         x2c = 0.25*( x2(j,k,1) + x2(j+1,k,1) 
     .         + x2(j+1,k+1,1) + x2(j,k+1,1) )
         y2c = 0.25*( y2(j,k,1) + y2(j+1,k,1) 
     .         + y2(j+1,k+1,1) + y2(j,k+1,1) )
         z2c = 0.25*( z2(j,k,1) + z2(j+1,k,1) 
     .         + z2(j+1,k+1,1) + z2(j,k+1,1) )
      end if
c     bi-quadratic
      if (ifit .eq. 2) then
         x2c = 0.5* ( x5 + x6 + x7 + x8 )
     .       -0.25*( x2(j,k,1)     + x2(j+1,k,1)
     .       +       x2(j+1,k+1,1) + x2(j,k+1,1) )
         y2c = 0.5* ( y5 + y6 + y7 + y8 )
     .       -0.25*( y2(j,k,1)    + y2(j+1,k,1) 
     .       +      y2(j+1,k+1,1) + y2(j,k+1,1) )
         z2c = 0.5* ( z5 + z6 + z7 + z8 )
     .       -0.25*( z2(j,k,1)     + z2(j+1,k,1) 
     .       +       z2(j+1,k+1,1) + z2(j,k+1,1) )      
      end if
c     quadratic in xie, linear in eta
      if (ifit .eq. 3) then
         x2c = .5*(x5 + x6)
         y2c = .5*(y5 + y6)
         z2c = .5*(z5 + z6)
      end if
c     linear in xie, quadratic in eta
      if (ifit .eq. 4) then
         x2c = .5*(x7 + x8)
         y2c = .5*(y7 + y8)
         z2c = .5*(z7 + z8)
      end if
c
      x2fit(j,k,1) = x2c
      y2fit(j,k) = y2c
      z2fit(j,k) = z2c
c
c     for each "to" cell, compare the "to" cell center locations calculated
c     above via generalized coordinate interpolation with the "to" cell
c     center locations calculated by applying the appropriate ifit to the
c     "to" grid points.  a "percnt" percent difference (scaled with the
c     "to" cell dimensions) is taken as significant.
c
      x00 = x2(j,k,1)
      y00 = y2(j,k,1)
      z00 = z2(j,k,1)
      xscal = 0.
      yscal = 0.
      zscal = 0.
      do 1621 jj=1,2
      do 1621 kk=1,2
      xsc = ccabs(x2(j+jj-1,k+kk-1,1)-x00)
      if(real(xsc).gt.real(xscal))xscal=xsc
      ysc = ccabs(y2(j+jj-1,k+kk-1,1)-y00)
      if(real(ysc).gt.real(yscal))yscal=ysc
      zsc = ccabs(z2(j+jj-1,k+kk-1,1)-z00)
      if(real(zsc).gt.real(zscal))zscal=zsc
1621  continue
      if(real(xscal).eq.0) xscal=-1.
      if(real(yscal).eq.0) yscal=-1.
      if(real(zscal).eq.0) zscal=-1.
c
c     check errors only in transverse directions - not projected direction
c     this eliminates spurious large error messages in nearly planar regions,
c     but therefore also cannot detect large gaps between zones (true error)
c
      xerr = 0
      yerr = 0.
      zerr = 0.
      if (itoss0 .eq. 0) then
         call direct(x15,x16,x17,x18,y15,y16,y17,y18,z15,z16,z17,z18,
     .               a1,a2,a3,itoss,nou,bou,nbuf,ibufdim)
      else 
         itoss = itoss0
      end if
      if(itoss.eq.1) then
        yerr = ccabs(y2int(j,k,1) - y2c)/yscal
        zerr = ccabs(z2int(j,k,1) - z2c)/zscal
      end if
      if(itoss.eq.2) then
        xerr = ccabs(x2int(j,k,1) - x2c)/xscal
        zerr = ccabs(z2int(j,k,1) - z2c)/zscal
      end if
      if(itoss.eq.3) then
        xerr = ccabs(x2int(j,k,1) - x2c)/xscal
        yerr = ccabs(y2int(j,k,1) - y2c)/yscal
      end if
      error = xerr
      if(real(yerr) .gt. real(error)) error = yerr
      if(real(zerr) .gt. real(error)) error = zerr
      if(real(error).gt.real(errmax)) then
        errmax = error
        jerrmx = j
        kerrmx = k
      end if
      if(real(xerr).gt.real(percnt) .or. 
     .   real(yerr).gt.real(percnt) .or.
     .   real(zerr).gt.real(percnt)) then 
        nbpt = nbpt + 1
c
        nou(2) = min(nou(2)+1,ibufdim)
        write(bou(nou(2),2),'('' interpolation no. '',i5,
     .  '' for to cell j,k '',i5,'','',i5)') icheck,j,k
        nou(2) = min(nou(2)+1,ibufdim)
        write(bou(nou(2),2),'(''    xint,yint,zint= '',f11.5,'','',
     .  f11.5,'','',f11.5)') real(x2int(j,k,1)),real(y2int(j,k,1)),
     .  real(z2int(j,k,1))
        nou(2) = min(nou(2)+1,ibufdim)
        write(bou(nou(2),2),'(''    x,y,z= '',f11.5,'','',
     .  f11.5,'','',f11.5)') real(x2c),real(y2c),real(z2c)
        nou(2) = min(nou(2)+1,ibufdim)
        write(bou(nou(2),2),'(''    xscal,yscal,zscal= '',f11.5,'','',
     .  f11.5,'','',f11.5)') real(xscal),real(yscal),real(zscal)
c
      end if
c
9113  continue
c
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),726) nbpt,int(real(percnt)*100)
726   format('     ',i4,' interpolated cell centers differ by',
     .' more than ',i3,' percent')
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),7261) 
7261  format('          from those obtained directly from the grid',
     .' points')
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),727) real(errmax*100.),jerrmx,kerrmx
727   format('          the max. difference ',e10.3,
     .' percent occurs at j,k = ',i4,',',i4)
      if(nbpt .gt. 0) then
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),*) '           note: differences between',
     .' interpolated and directly obtained cell'
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),*) '           centers often (correctly)',
     .' arise when boundaries are being rendered'
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),*) '           coincident. however, if',
     .' differences occur when boundaries are'
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),*) '           not rendered coincident, or',
     .' if the other diagnostic counts are'
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),*) '           non-zero, a problem in the',
     .' input or the grid is likely' 
      end if
c
c     ****** seventh check ******
c
c     check normals on (cell center) grids as calculated 1)directly from
c     the "to" grid, and 2) via generalized coordinate interpolation
c
c
c     not applicable in 2D cases or for cases with orphans!
c
      jseg = j22 - j21
      kseg = k22 - k21
      if(jseg.gt.1 .and. kseg.gt.1 .and. iorph.le.0)then
c
      epsn = .1
      nnorm = 0
      prdsum = 0.
      dnmx = 0.
      jdnmx = 1
      kdnmx = 1
c
      do 7113 j=j21,j22-2
      do 7113 k=k21,k22-2
c
c     normals at "to" cell centers calculated directly from "to" grid points
c
      ax = x2fit(j+1,k+1,1) - x2fit(j,k,1)
      ay = y2fit(j+1,k+1) - y2fit(j,k)
      az = z2fit(j+1,k+1) - z2fit(j,k)
      bx = x2fit(j,k+1,1) - x2fit(j+1,k,1)
      by = y2fit(j,k+1) - y2fit(j+1,k)
      bz = z2fit(j,k+1) - z2fit(j+1,k)
      dnx = ay*bz - az*by
      dny = az*bx - ax*bz
      dnz = ax*by - ay*bx
      d = sqrt(dnx*dnx +dny*dny +dnz*dnz)
      if (real(d) .le. 0.) d=1.
      xnfit = dnx/d
      ynfit = dny/d
      znfit = dnz/d
c
c     normals at "to" cell centers calculated from "to" cell centers
c
      ax = x2int(j+1,k,1) - x2int(j,k,1)
      ay = y2int(j+1,k,1) - y2int(j,k,1)
      az = z2int(j+1,k,1) - z2int(j,k,1)
      bx = x2int(j,k+1,1) - x2int(j,k,1)
      by = y2int(j,k+1,1) - y2int(j,k,1)
      bz = z2int(j,k+1,1) - z2int(j,k,1)
      dnx = ay*bz - az*by
      dny = az*bx - ax*bz
      dnz = ax*by - ay*bx
      d = sqrt(dnx*dnx +dny*dny +dnz*dnz)
      if (real(d) .le. 0) d=1.
      xnint = dnx/d
      ynint = dny/d
      znint = dnz/d
c
c     inner product of the two unit normals; identical unit normals will 
c     give an inner product of 1. 
c
      prod = xnfit*xnint + ynfit*ynint + znfit*znint
      prdsum = prdsum + prod
      if (real(prod) .lt. 1.-real(epsn) .or.
     .    real(prod) .gt. 1.+real(epsn)) then 
        nnorm = nnorm + 1
      end if
      dn = prod - 1.
      if(abs(real(dn)).gt.abs(real(dnmx))) then
        dnmx = dn
        jdnmx = j
        kdnmx = k
      end if
c
7113  continue
c
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),728) nnorm,int(real(epsn)*100)
728   format('     ',i4,' normals at interpolated cell centers', 
     .' differ by more than ',i3,' percent')
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),*) '        from the normals at cell centers',
     .            ' obtained directly from grid points'
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),729) real(dnmx*100.),jdnmx,kdnmx
729   format('       the max. difference ',e10.3,
     .' percent occurs at j,k = ',i4,',',i4)
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),730) real(prdsum)/(j22-j21-1)/(k22-k21-1)
730   format('       the average inner product',
     .' over the interface is ',e10.3)
c
      if(nnorm .gt. 0) then  
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),1222)
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),1223)
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),1226)
      end if
      end if
c
c     ****** eigth check ******
c    
c     Following output for visual check of the interpolation procedure.
c     Output, in order: the "from" grid(s), the "to" grid, and the "to" 
c     cell centers found from generalized coordinate interpolation.  
c
c     This diagnostic can be enabled by setting iifit < 0 in the input file.
c     It is automatically enabled when the diagnostics detect a potential 
c     problem with the interpolation coefficients, or if the search routine 
c     fails. In the latter case, the last grid output is simply the "to" 
c     cell center at which the search failed.  
c     
552   continue
      if(iout .gt. 0) then
c
      iunit=40
c
      if (icheck.gt.99) then
         len1 = 13
         write (titlptchgrd,'("patch_p3d.",i3)') icheck
      else if (icheck.gt.9) then
         len1 = 12
         write (titlptchgrd,'("patch_p3d.",i2)') icheck
      else
         len1 = 11
         write (titlptchgrd,'("patch_p3d.",i1)') icheck
      endif
      do i = len1+1, 14
         titlptchgrd(i:i) = ' '
      end do
      open(iunit,file=titlptchgrd(1:len1),form='formatted',
     .status='unknown')
      rewind(iunit)
      if(iorph.le.0) then
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),791) titlptchgrd(1:len1)
791     format('   plot3d (/mg/for) file for this patch',
     .  ' interface written to',a14)
      else
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),792) titlptchgrd(1:len1)
792     format('   plot3d (/mg/for/blank) file for this patch',
     .  ' interface written to',a14)
      end if
      if(iorph.le.0) then
        write(iunit,*)lmax1+2
      else
        write(iunit,*)lmax1+3
      end if
c
c     iwhole = 0 output original "from" grids
c     iwhole = 1 output expanded "from" grids (useful for examining cases
c                where the search routine has failured near a boundary)
c     
      iwhole = 0
c
      maxx = 0
      if (iwhole.eq.0) then 
         do l=1,lmax1
            if (xif1(l).eq.1) xif1(l) = xif1(l) + 1
            if (etf1(l).eq.1) etf1(l) = etf1(l) + 1
            if (xif2(l).eq.jjmax1(l)) xif2(l) = xif2(l) - 1
            if (etf2(l).eq.kkmax1(l)) etf2(l) = etf2(l) - 1
         end do
      end if
      if(istop.ne.1) then
        if(iorph.le.0) then
          write(iunit,*) ((xif2(l)-xif1(l)+1-maxx),
     .    (etf2(l)-etf1(l)+1-maxx),1,l=1,lmax1),
     .    (j22-j21+1),(k22-k21+1),1,(j22-j21),(k22-k21),1
        else
          write(iunit,*) ((xif2(l)-xif1(l)+1-maxx),
     .    (etf2(l)-etf1(l)+1-maxx),1,l=1,lmax1),
     .    (j22-j21+1),(k22-k21+1),1,(j22-j21),(k22-k21),1,
     .    (j22-j21),(k22-k21),1
        end if
      else
        write(iunit,*) ((xif2(l)-xif1(l)+1-maxx),
     .  (etf2(l)-etf1(l)+1-maxx),1,l=1,lmax1),
     .  (j22-j21+1),(k22-k21+1),1,1,1,1
      end if
      do 9112 l=1,lmax1
      js = xif1(l)
      ks = etf1(l)
      je = xif2(l)
      ke = etf2(l)
      if(iorph.le.0) then
        if(ialph.eq.0) then
          write(iunit,*) ((real(x1(j,k,l)),j=js,je),k=ks,ke),
     .                   ((real(y1(j,k,l)),j=js,je),k=ks,ke),
     .                   ((real(z1(j,k,l)),j=js,je),k=ks,ke)
        else
          write(iunit,*) ((real(x1(j,k,l)),j=js,je),k=ks,ke),
     .                   ((real(z1(j,k,l)),j=js,je),k=ks,ke),
     .                   ((real(-y1(j,k,l)),j=js,je),k=ks,ke)
        end if
      else
        if(ialph.eq.0) then
          write(iunit,*) ((real(x1(j,k,l)),j=js,je),k=ks,ke),
     .                   ((real(y1(j,k,l)),j=js,je),k=ks,ke),
     .                   ((real(z1(j,k,l)),j=js,je),k=ks,ke),
     .                   ((              1,j=js,je),k=ks,ke)
        else
          write(iunit,*) ((real(x1(j,k,l)),j=js,je),k=ks,ke),
     .                   ((real(z1(j,k,l)),j=js,je),k=ks,ke),
     .                   ((real(-y1(j,k,l)),j=js,je),k=ks,ke),
     .                   ((               1,j=js,je),k=ks,ke)
        end if
      end if
9112  continue
      if(iorph.le.0) then
        if(ialph.eq.0) then
          write(iunit,*) ((real(x2(j,k,1)),j=j21,j22),k=k21,k22),
     .                   ((real(y2(j,k,1)),j=j21,j22),k=k21,k22),
     .                   ((real(z2(j,k,1)),j=j21,j22),k=k21,k22)
        else
          write(iunit,*) ((real(x2(j,k,1)),j=j21,j22),k=k21,k22),
     .                   ((real(z2(j,k,1)),j=j21,j22),k=k21,k22),
     .                   ((real(-y2(j,k,1)),j=j21,j22),k=k21,k22)
        end if
      else
        if(ialph.eq.0) then
          write(iunit,*) ((real(x2(j,k,1)),j=j21,j22),k=k21,k22),
     .                   ((real(y2(j,k,1)),j=j21,j22),k=k21,k22),
     .                   ((real(z2(j,k,1)),j=j21,j22),k=k21,k22),
     .                   ((              1,j=j21,j22),k=k21,k22)
        else
          write(iunit,*) ((real(x2(j,k,1)),j=j21,j22),k=k21,k22),
     .                   ((real(z2(j,k,1)),j=j21,j22),k=k21,k22),
     .                   ((real(-y2(j,k,1)),j=j21,j22),k=k21,k22),
     .                   ((               1,j=j21,j22),k=k21,k22)
        end if
      end if
      if(istop.eq.1) then 
        if(ialph.eq.0) then
          write(iunit,*) real(xc),real(yc),real(zc)
        else
          write(iunit,*) real(xc),real(zc),real(-yc)
        end if
      else
      if(iorph.le.0) then
        if(ialph.eq.0) then
          write(iunit,*)
     .          ((real(x2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(y2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(z2int(j,k,1)),j=j21,j22-1),k=k21,k22-1)
        else
          write(iunit,*) 
     .          ((real(x2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(z2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(-y2int(j,k,1)),j=j21,j22-1),k=k21,k22-1)
        end if
      else
        do 1700 j=j21,j22-1
        do 1700 k=k21,k22-1
        ll = (j22-j21)*(k-k21) + (j-j21+1)
        l1   = mblkpt(ll)
        x2fit(j,k,1) = 1.
        if(l1.eq.0) x2fit(j,k,1) = 0.
1700    continue
        if(ialph.eq.0) then
          write(iunit,*)
     .          ((real(x2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(y2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(z2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .           ((int(x2fit(j,k,1)),j=j21,j22-1),k=k21,k22-1)
        else
          write(iunit,*)
     .          ((real(x2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(z2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(-y2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .            ((int(x2fit(j,k,1)),j=j21,j22-1),k=k21,k22-1)
        end if
        do 1710 j=j21,j22-1
        do 1710 k=k21,k22-1
        ll = (j22-j21)*(k-k21) + (j-j21+1)
        l1   = mblkpt(ll)
        x2fit(j,k,1) = 0.
        if(l1.eq.0) then
          x2fit(j,k,1) = 1.
          x2int(j,k,1) = 0.25*( x2(j,k,1) + x2(j+1,k,1)
     .         + x2(j+1,k+1,1) + x2(j,k+1,1) )
          y2int(j,k,1) = 0.25*( y2(j,k,1) + y2(j+1,k,1)
     .         + y2(j+1,k+1,1) + y2(j,k+1,1) )
          z2int(j,k,1) = 0.25*( z2(j,k,1) + z2(j+1,k,1)
     .         + z2(j+1,k+1,1) + z2(j,k+1,1) )

        end if
1710    continue
        if(ialph.eq.0) then
          write(iunit,*)
     .          ((real(x2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(y2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(z2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .           ((int(x2fit(j,k,1)),j=j21,j22-1),k=k21,k22-1)
        else
          write(iunit,*)
     .          ((real(x2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(z2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .          ((real(-y2int(j,k,1)),j=j21,j22-1),k=k21,k22-1),
     .            ((int(x2fit(j,k,1)),j=j21,j22-1),k=k21,k22-1)
        end if
      end if
      end if
      close(iunit)
      if(istop .gt. 0) call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
1222  format(12x,'possible causes:') 
1223  format(14x,'...mmcxie and/or mmceta incorrectly set')
1224  format(14x,'..."from" block input data incorrectly set')
1225  format(14x,'...small gaps between "from" blocks')
1226  format(14x,'..."to" and "from" blocks do not lie on the same',
     .           ' surface')
1227  format(14x,'...factj or factk too small')
1228  format(14x,'...mmcxie incorrectly set')
1229  format(14x,'...mmceta incorrectly set')
1330  format(14x,'...interface not really C-0, check grid and/or',
     .           ' C-0 flag')
1331  format(14x,'...search routine error...uh-oh')
1332  format(14x,'...these may be the orphan points anticipated',
     .           ' by setting iorph > 0')
1333  format(14x,'...branch cut on the "from" side may falsely',
     .           ' trigger this message')
1334  format(14x,'   and the',
     .           ' results may in fact',
     .           ' be OK - verify with patch_p3d.xx file')
c
      return
      end
