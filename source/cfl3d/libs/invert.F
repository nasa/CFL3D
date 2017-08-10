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
      subroutine invert(jdim1,kdim1,msub1,msub2,jjmax1,
     .                  kkmax1,lmax1,x1,y1,z1,x1mid,y1mid,z1mid,
     .                  x1mide,y1mide,z1mide,limit0,jjmax2,kkmax2,
     .                  x2,y2,z2,xie2,eta2,mblkpt,temp,jimage,kimage,
     .                  ifit,itmax,sxie,seta,sxie2,seta2,xie2s,eta2s,
     .                  intmx,icheck,nblkj,nblkk,jmm,kmm,mcxie,mceta,
     .                  lout,j21,j22,k21,k22,npt,ic0,iorph,itoss0,
     .                  ncall,ioutpt,xif1,xif2,etf1,etf2,iself,ifiner,
     .                  xie2f,eta2f,mblkptf,nptf,xi1f,xi2f,et1f,et2f,
     .                  iavg,nou,bou,nbuf,ibufdim,myid,
     .                  mblk2nd,maxbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determine generalized coordinates of cell centers of the
c     "to" grid  in terms of the generalized coordinate system(s) defined
c     on the "from" grid(s)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,avgg,avgq
      character*14 titlptchgrd
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension xie2f(nptf),eta2f(nptf),mblkptf(nptf)
      dimension jimage(msub1,jdim1,kdim1),kimage(msub1,jdim1,kdim1)
      dimension x1(jdim1,kdim1,msub1),y1(jdim1,kdim1,msub1),
     .          z1(jdim1,kdim1,msub1)
      dimension x1mid(jdim1,kdim1,msub1),y1mid(jdim1,kdim1,msub1),
     .          z1mid(jdim1,kdim1,msub1)
      dimension x1mide(jdim1,kdim1,msub1),y1mide(jdim1,kdim1,msub1),
     .          z1mide(jdim1,kdim1,msub1)
      dimension x2(jdim1,kdim1,msub2),y2(jdim1,kdim1,msub2),
     .          z2(jdim1,kdim1,msub2)
      dimension xie2(npt),eta2(npt),temp(jdim1*kdim1),
     .          mblkpt(npt)
      dimension jjmax1(msub1),kkmax1(msub1),jjmax2(msub2),kkmax2(msub2)
      dimension sxie(jdim1,kdim1,msub1),seta(jdim1,kdim1,msub1)
      dimension sxie2(jdim1,kdim1,msub2),seta2(jdim1,kdim1,msub2)
      dimension xie2s(jdim1,kdim1),eta2s(jdim1,kdim1),nblkj(jdim1),
     .          nblkk(kdim1),jmm(kdim1),kmm(jdim1)
      dimension mblk2nd(maxbl)
      integer   lout(msub1),xif1(msub1),xif2(msub1),etf1(msub1),
     .          etf2(msub1)
      integer xi1f,xi2f,et1f,et2f
c
      common /sklt1/isklt1
      common /areas/ ap(3),imaxa
      common /tacos/ iretry
      common /tracer/ itrace
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
c
      l2    = 1
      jmax2 = jjmax2(l2)
      kmax2 = kkmax2(l2)
      istop = 0
      iout  = 0
      igap  = 0
      if(ifit .lt. 0) iout = 1
      ifit  = abs(ifit)
c
c     itrace < 0, do not write search history for current "to" cell
c     itrace = 0, overwrite history from previous "to" cell with current 
c     itrace = 1, retain the search history for ALL cells (may get huge file)
c     trace output found in fort.7
c
      itrace = -1
c
      idum1 = 0
      idum2 = 0
      idum3 = 0
      idum4 = 0
      dum1  = 0.
      dum2  = 0.
      dum3  = 0.
c
      if (isklt1.gt.0) then
      nou(4) = min(nou(4)+1,ibufdim)
      write(bou(nou(4),4),*)'    beginning inversion for xie,',
     .           ' eta of cell centers'
      if (iavg.eq.0 .or. ifiner.eq.0) then
      if (ifit.eq.1) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),*)'      using bi-linear fit'
      end if
      if (ifit.eq.2) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),*)'      using bi-quadratic fit'
      end if
      if (ifit.eq.3) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),*)'      using quadratic fit',
     .                          ' in xie, linear fit in eta' 
      end if
      if (ifit.eq.4) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),*)'      using linear fit in',
     .                          ' xie, quadratic fit in eta'
      end if
      else if (ifiner.gt.0) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),171) ifiner
 171     format('     using averages of finer level:',
     .   ' interpolation number',i4)
      end if
      end if
c
c     call trace(0,icheck,idum2,idum3,idum4,dum1,dum2,dum3)
c
c     compute cell edge midpoints of "from" grid(s) 
c     using quadratic least squares
c
      do 999 l=1,lmax1
      jmax1 = jjmax1(l)
      kmax1 = kkmax1(l)
      do 500 j=1,jmax1-1
      do 500 k=1,kmax1
      jl = 2
      jr = jmax1-2
      call extra(jdim1,kdim1,msub1,l,x1,y1,z1,
     .           j,k,jl,jr,x5,y5,z5,icase,ifit)
      x1mid(j,k,l) = x5
      y1mid(j,k,l) = y5
      z1mid(j,k,l) = z5
  500 continue
      do 600 k=1,kmax1-1
      do 600 j=1,jmax1
      kl = 2
      kr = kmax1-2
      call extrae(jdim1,kdim1,msub1,l,x1,y1,z1,
     .            j,k,kl,kr,x7,y7,z7,icase,ifit)
      x1mide(j,k,l) = x7
      y1mide(j,k,l) = y7
      z1mide(j,k,l) = z7
  600 continue
  999 continue
c
c     don't go through inversion process if coarser level
c     interpolation data are always obtained by averaging 
c     fine level data
c
      if (ifiner.gt.0 .and. iavg.gt.0) go to 552
c
      nblsav= 1
      jl    = 1
      jr    = jmax2-1
      kl    = 1
      kr    = kmax2-1
c
c     loop over all "to" cells
c
      do 1000 j=j21,j22-1
      do 1001 k=k22-1,k21,-1
c
c     compute cell edge midpoints of "to" grid 
c     using quadratic least squares
c
      if (k.eq.k22-1) then
        kcall = k+1
        call extra(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .             j,kcall,jl,jr,x6,y6,z6,icase,ifit)
      end if 
      call extra(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .           j,k,jl,jr,x5,y5,z5,icase,ifit)
      call extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .            j,k,kl,kr,x7,y7,z7,icase,ifit)
      jcall = j+1
      call extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .            jcall,k,kl,kr,x8,y8,z8,icase,ifit)
c
c     compute normalized directed areas/unit normals of "to" cell
c
      if (itoss0 .eq. 0) then
         call direct(x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,
     .               a1,a2,a3,imaxa,nou,bou,nbuf,ibufdim)
         ap(1) = a1
         ap(2) = a2
         ap(3) = a3
      end if
c
      ifits = ifit
      iretry = 0
      icount = 1
      iorfn = 0
17085 continue
c
c     compute center of "to" cell, consistent with ifit
c
c     bi-linear
      if (ifit .eq. 1) then
         xc = 0.25*( x2(j,k,l2) + x2(j+1,k,l2) 
     .         + x2(j+1,k+1,l2) + x2(j,k+1,l2) )
         yc = 0.25*( y2(j,k,l2) + y2(j+1,k,l2) 
     .         + y2(j+1,k+1,l2) + y2(j,k+1,l2) )
         zc = 0.25*( z2(j,k,l2) + z2(j+1,k,l2) 
     .         + z2(j+1,k+1,l2) + z2(j,k+1,l2) )
      end if
c     bi-quadratic
      if (ifit .eq. 2) then
         xc = 0.5* ( x5 + x6 + x7 + x8 )
     .       -0.25*( x2(j,k,l2)     + x2(j+1,k,l2)
     .       +       x2(j+1,k+1,l2) + x2(j,k+1,l2) )
         yc = 0.5* ( y5 + y6 + y7 + y8 )
     .       -0.25*( y2(j,k,l2)    + y2(j+1,k,l2) 
     .       +      y2(j+1,k+1,l2) + y2(j,k+1,l2) )
         zc = 0.5* ( z5 + z6 + z7 + z8 )
     .       -0.25*( z2(j,k,l2)     + z2(j+1,k,l2) 
     .       +       z2(j+1,k+1,l2) + z2(j,k+1,l2) )      
      end if
c     quadratic in xie, linear in eta
      if (ifit .eq. 3) then
         xc = .5*(x5 + x6)
         yc = .5*(y5 + y6)
         zc = .5*(z5 + z6)
      end if
c     linear in xie, quadratic in eta
      if (ifit .eq. 4) then
         xc = .5*(x7 + x8)
         yc = .5*(y7 + y8)
         zc = .5*(z7 + z8)
      end if
c
c     call trace(1,icheck,j,k,ifit,xc,yc,zc)
c
c     set starting point for search routine
c     first time.........use solution at the previous point, or search
c                        for minimum distance point if xie=eta=0.
c     subsequent times...use solution at the same point from the
c                        previous time.
c
c     xiet,etat are generalized coordinates in expanded "from" grid
c     xie2,eta2 are generalized coordinates in original "from" grid
c
      if (iretry.eq.0) then
         if (ncall.eq.1) then
            if (k.lt.k22-1) then
               ll = (j22-j21)*(k+1-k21) + (j-j21+1)
               if (mblkpt(ll).ne.0) then
                  xiet = xie2(ll) + 1.
                  etat = eta2(ll) + 1.
                  l1   = mblkpt(ll) 
               else
                  xiet = 0.
                  etat = 0.
                  l1   = nblsav
                  l1   = max(l1,1)
               end if
            else if (j.eq.j21) then
               xiet = 0.
               etat = 0.
               l1   = 1
            else if (j.gt.j21) then
               ll = (j22-j21)*(k-k21) + (j-1-j21+1)
               if (mblkpt(ll).ne.0) then
                  xiet = xie2(ll) + 1.
                  etat = eta2(ll) + 1.
                  l1   = mblkpt(ll) 
               else
                  xiet = 0.
                  etat = 0.
                  l1   = nblsav
                  l1   = max(l1,1)
               end if
            end if
         else
            ll = (j22-j21)*(k-k21) + (j-j21+1)
            xiet = xie2(ll) + 1.
            etat = eta2(ll) + 1.
            l1   = mblkpt(ll)
            l1   = max(l1,1)
         end if
       else
          ll = (j22-j21)*(k-k21) + (j-j21+1)
          if (mblkpt(ll).ne.0) then
             xiet = xie2(ll) + 1.
             etat = eta2(ll) + 1.
             l1   = mblkpt(ll) 
          else
             xiet = 0.
             etat = 0.
             l1   = nblsav
             l1   = max(l1,1)
          end if
      end if
c
17086 continue
c
c     search over "from" blocks to find xie,eta associated with "to"
c     cell center xc,yc,zc.  
c
      call topol(jdim1,kdim1,msub1,jjmax1,kkmax1,lmax1,l1,x1,y1,z1,
     .           x1mid,y1mid,z1mid,x1mide,y1mide,z1mide,limit0,xc,yc,zc,
     .           xiet,etat,jimage,kimage,ifit,itmax,igap,iok,lout,ic0,
     .           itoss0,j,k,iself,xif1,xif2,etf1,etf2,nou,bou,nbuf,
     .           ibufdim,myid)
c
c     search routine unsuccessful...try an alternative polynomial fit
c
      if (iretry.eq.0 .and. iok.ne.1) then
        if (icount.lt.4) then
           call newfit(ifits,ifit,icount)
           icount = icount + 1
           go to 17085
        else
            if (iorph.le.0) then
               if (ifiner.gt.0) then
                  nou(4) = min(nou(4)+1,ibufdim)
                  write(bou(nou(4),4),172) 
                  nou(4) = min(nou(4)+1,ibufdim)
                  write(bou(nou(4),4),272) j,k
                  nou(4) = min(nou(4)+1,ibufdim)
                  write(bou(nou(4),4),372) ifiner
                  nou(4) = min(nou(4)+1,ibufdim)
                  write(bou(nou(4),4),472)
 172              format('        all attempts to find generalized',
     .            ' coordinates')
 272              format('        of cell center j,k =',i4,
     .            ',',i4,' have been unsuccessful...')
 372              format('        will use  averages of finer level:',
     .            ' interpolation number',i4)
 472              format('        for ALL points on this interface')
                  istop = 1
                  go to 552
c
               else
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
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),'('' program terminated in'',
     .            '' dynamic patching routines - see file '',a60)')
     .            grdmov
                  nou(4) = min(nou(4)+1,ibufdim)
                  write(bou(nou(4),4),173)
                  nou(4) = min(nou(4)+1,ibufdim)
                  write(bou(nou(4),4),273) j,k,titlptchgrd(1:len1)
 173              format('        stopping...all attempts to find',
     .            ' generalized coordinates of')
 273              format(
     .            '        cell center j,k =',i4,',',i4,
     .            ' have been unsuccessful...check',a14)
                  istop = 1
                  iout  = 1
                  go to 552
               end if
c
            else
c
c              mark orphan point as interpolated from block "0" - will
c              not be interpolated in cfl3d - this option used for points 
c              really do lie outside the domain of the from block - for
c              example, an exposed flap surface.
c
               iok = 1
               ll = (j22-j21)*(k-k21) + (j-j21+1)
               nblsav            = mblkpt(ll)
               xie2(ll)          = 1.
               eta2(ll)          = 1.
               mblkpt(ll)        = 0
               x6 = x5
               y6 = y5
               z6 = z5
c     call trace(50,j,k,idum3,idum4,dum1,dum2,dum3)
               go to 1001
            end if
        end if
      end if
c
c     search successful with current polynomial fit
c
      if(iok.eq.1) then
        ll = (j22-j21)*(k-k21) + (j-j21+1)
        xie2(ll)          = xiet-1.
        eta2(ll)          = etat-1.
        mblkpt(ll)        = l1
        ifit1              = ifit
      end if
c
      if(iorph.gt.0) then
c
c       search routine "succesful", but only because it found what might 
c       otherwise be an orphan point in the expanded "from" grid.
c
        ll = (j22-j21)*(k-k21) + (j-j21+1)
        if(real(xie2(ll)).lt.1. .or. real(xie2(ll)).gt.jjmax1(l1)-2 .or.
     .     real(eta2(ll)).lt.1. .or. real(eta2(ll)).gt.kkmax1(l1)-2)then
          iorfn = iorfn +1
          if(iorfn.le.1)then
c
c           redo the search for this point, starting with a minimum 
c           distance search to make sure this point should really
c           be marked as an orhan
c
            iok    = 0
            icount = 1
            iretry = 0
            xiet   = 0
            etat   = 0.
            go to 17086
          else
c
c           mark the point as an orphan
c
            iok = 1
            nblsav            = mblkpt(ll)
            xie2(ll)          = 1.
            eta2(ll)          = 1.
            mblkpt(ll)        = 0
            x6 = x5
            y6 = y5
            z6 = z5
c     call trace(50,j,k,idum3,idum4,dum1,dum2,dum3)
            go to 1001
          end if
        end if
      end if
c
c     if current polynomial fit is not the one which was input, 
c     try again with the input value
c
      if(ifit.ne.ifits) then
        iretry = 1
        ifit = ifits
        go to 17085
      end if
c
c     if second try with input value of ifit not successful, use the ifit
c     value for which the search routine was successful.
c
      if(iretry.eq.1) then
        if(iok.ne.1) then
          if(isklt1.gt.0) then
          nou(4) = min(nou(4)+1,ibufdim)
          write(bou(nou(4),4),*)'         iterations using original',
     .    ' fit not successful at j,k= ',j,k
          if (ifit.eq.1) then
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),*)
     .    '          used bilinear fit instead'
          end if
          if (ifit.eq.2) then
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),*)
     .    '          used biquadratic fit instead'
          end if
          if (ifit.eq.3) then
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),*)
     .    '          used quadratic fit in xie, linear fit in',
     .    ' eta instead'
          end if
          if (ifit.eq.4) then
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),*)
     .    '          used linear fit in xie, quadratic fit in',
     .    ' eta instead'
          end if
          end if
        end if
      end if
c
      x6 = x5
      y6 = y5
      z6 = z5
c
 1001 continue
 1000 continue
c
c     check and, if necessary, correct values of xie and/or eta    
c     near boundaries. 
c
c     ishear > 0: try a shearing correction to make the "to" 
c                 and "from" boundaries coincident; if shearing fails, then 
c                 try an arc length correction.  
c     ishear < 0: use only the arc length correction
c
      ishear = 1
c
      call shear(ishear,istop,iout,igap,jdim1,kdim1,msub1, 
     .           msub2,jjmax1,kkmax1,lmax1,x1,y1,z1,x1mid,y1mid,
     .           z1mid,x1mide,y1mide,z1mide,limit0,jjmax2,kkmax2,
     .           x2,y2,z2,xie2,eta2,mblkpt,temp,jimage,kimage,ifit,
     .           itmax,xc,yc,zc,sxie2,seta2,jcorr,kcorr,intmx,icheck,
     .           nblkj,nblkk,jmm,kmm,mcxie,mceta,lout,j21,j22,k21,k22,
     .           npt,ic0,iorph,itoss0,xif1,xif2,etf1,etf2,iself,ifiner,
     .           nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl)

c
      if (jcorr.ne.0 .or. kcorr.ne.0)
     .   call arc(jdim1,kdim1,msub1,msub2,jjmax1,kkmax1,
     .            lmax1,x1,y1,z1,limit0,jjmax2,kkmax2,x2,y2,z2,xie2,
     .            eta2,mblkpt,ifit,itmax,jcorr,kcorr,sxie,seta,sxie2,
     .            seta2,xie2s,eta2s,intmx,icheck,nblkj,nblkk,jmm,kmm,
     .            j21,j22,k21,k22,npt,xif1,xif2,etf1,etf2,
     .            nou,bou,nbuf,ibufdim,mblk2nd,maxbl)
c
552   continue
c
c     for coarse level interfaces where the standard search
c     algorithim has failed, use average of finer level data
c
      if (ifiner.gt.0 .and. istop.eq.1) then
         call avgint(xie2,eta2,mblkpt,npt,xie2f,eta2f,mblkptf,
     .               nptf,j21,j22,k21,k22,xi1f,xi2f,et1f,et2f)
         istop = 0
      end if
c
c     perform diagnostic checks on xie,eta values found for "to" cell
c     centers
c
      if (ioutpt .gt.0.or. istop.gt.0) then
      call diagnos(istop,iout,igap,jdim1,kdim1,msub1,
     .             msub2,jjmax1,kkmax1,lmax1,x1,y1,z1,x1mid,y1mid,
     .             z1mid,x1mide,y1mide,z1mide,sxie,seta,sxie2,seta2,
     .             xie2s,eta2s,jjmax2,kkmax2,x2,y2,
     .             z2,xie2,eta2,mblkpt,icheck,intmx,xc,yc,zc,ifit,
     .             j21,j22,k21,k22,npt,ic0,iorph,xif1,xif2,
     .             etf1,etf2,itoss0,iself,
     .             nou,bou,nbuf,ibufdim,myid)
      end if
      return
      end
