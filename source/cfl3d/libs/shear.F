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
      subroutine shear(ishear,istop,iout,igap,jdim1,kdim1,msub1,
     .                 msub2,jjmax1,kkmax1,lmax1,x1,y1,z1,x1mid,
     .                 y1mid,z1mid,x1mide,y1mide,z1mide,limit0,jjmax2,
     .                 kkmax2,x2,y2,z2,xie2,eta2,mblkpt,temp,jimage,
     .                 kimage,ifit,itmax,xc,yc,zc,sxie2,seta2,jcorr,
     .                 kcorr,intmx,icheck,nblkj,nblkk,jmm,kmm,mcxie,
     .                 mceta,lout,j21,j22,k21,k22,npt,ic0,iorph,itoss0,
     .                 xif1,xif2,etf1,etf2,iself,ifiner,nou,bou,nbuf,
     .                 ibufdim,myid,mblk2nd,maxbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determine generalized coordinates of cell edge midpoints
c     on xie=0 and eta=0 boundaries, and determine the requisite shearing 
c     correction to the generalized coordinates near xie=0 and/or eta=0 
c     boundaries.
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
     .          mblk2nd(maxbl),mblkpt(npt)
      dimension jjmax1(msub1),kkmax1(msub1),jjmax2(msub2),kkmax2(msub2)
      dimension sxie2(jdim1,kdim1,msub2),seta2(jdim1,kdim1,msub2)
      dimension nblkj(jdim1),nblkk(kdim1),jmm(kdim1),kmm(jdim1)
      integer   lout(msub1),xif1(msub1),xif2(msub1),etf1(msub1),
     .          etf2(msub1)
c
      common /sklt1/isklt1
      common /areas/ ap(3),imaxa
      common /tacos/ iretry
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
c
      l2    = 1
      jmax2 = jjmax2(l2)
      kmax2 = kkmax2(l2)
      jl = 1
      jr = jmax2-1
      kl = 1
      kr = kmax2-1
c
      if (isklt1.gt.0) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),*)'    beginning check of boundary values'
      end if
c
c***************************************************************************
c      correct boundary values near eta=0
c***************************************************************************
c
      kcorr = 0
      ncoin = 0
c
      if (mceta.eq.0) then 
        if (isklt1.gt.0) then
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),*) '      eta=0 boundaries not rendered', 
     .     ' coincident'
        end if
        go to 670
      end if
      if(ishear.ge.0) then
        if (isklt1.gt.0) then
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),*) '      eta=0 boundaries being',
     .     ' rendered coincident via shearing method'
         end if
      else
        if (isklt1.gt.0) then
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),*) '      eta=0 boundaries being',
     .     ' rendered coincident via arc length method'
        end if
      end if
c
c     loop over all "to" cell on eta=0 boundary
c
      do 2000 j=j21,j22-1
c
c     compute edge midpoints of first layer of "to" grid cells
c     along the eta=0 boundary using quadratic least squares
c
      jcall = j
      kcall = k21
         call extra(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .              jcall,kcall,jl,jr,x5,y5,z5,icase,ifit)
      if(j. eq. j21) then
         call extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .              jcall,kcall,kl,kr,x7,y7,z7,icase,ifit)
      end if
      kcall = k21+1
         call extra(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .              jcall,kcall,jl,jr,x6,y6,z6,icase,ifit)
      jcall = j+1
      kcall = k21
         call extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .              jcall,kcall,kl,kr,x8,y8,z8,icase,ifit)

c
c     compute normalized directed areas/unit normals of "to" cell
c
      if (itoss0 .eq. 0) then
         call direct(x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,
     .                     a1,a2,a3,imaxa,nou,bou,nbuf,ibufdim)
         ap(1) = a1
         ap(2) = a2
         ap(3) = a3
      end if
c
      ifits = ifit
      iretry = 0
      icount = 1
18085 continue
c
c     compute midpoint of "to" cell on eta=0 boundary, consistent with ifit 
c
c     bilinear or linear in xie
      if (ifit.eq.1.or.ifit.eq.4) then
         xc = 0.5*(x2(j,k21,l2)+x2(j+1,k21,l2))
         yc = 0.5*(y2(j,k21,l2)+y2(j+1,k21,l2))
         zc = 0.5*(z2(j,k21,l2)+z2(j+1,k21,l2))
      else
c     biquadratic or quadratic in xie
         xc = x5
         yc = y5
         zc = z5
      end if
c
c     call trace(2,icheck,j,kcall,ifit,xc,yc,zc)
c
c     set starting point for search routine...use solution at cell center
c     of current cell
c
c     xiet,etat are generalized coordinates in expanded "from" grid
c     xie2,eta2 are generalized coordinates in original "from" grid
c
      ll = (j-j21+1)
      xiet   = xie2(ll) + 1.
      etat   = eta2(ll) + 1.
      l1     = mblkpt(ll)
c
c     search over "from" blocks to find xie,eta associated with "to"
c     cell midpoint xc,yc,zc on eta=0
c
      call topol(jdim1,kdim1,msub1,jjmax1,kkmax1,lmax1,l1,x1,y1,z1,
     .           x1mid,y1mid,z1mid,x1mide,y1mide,z1mide,limit0,xc,yc,zc,
     .           xiet,etat,jimage,kimage,ifit,itmax,igap,iok,lout,ic0,
     .           itoss0,j,1,iself,xif1,xif2,etf1,etf2,nou,bou,nbuf,
     .           ibufdim,myid)
c
c     search routine unsuccessful...try alternative polynomial fit
c
      if (iretry.eq.0 .and. iok.ne.1) then
        if (icount.lt.4) then
           call newfit(ifits,ifit,icount)
           icount = icount + 1
           go to 18085
        else
c
c         don't shear orphan points
c
          if(iorph.gt.0) then
            if (isklt1.gt.0) then
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),*)' quitting boundary correction',
     .      ' because points have been flagged as orphans'
            end if
            return
          end if
c
          if (ifiner.gt.0) then
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),172) 
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),272) j,1
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),372) ifiner
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),472)
 172         format('        all attempts to find generalized',
     .       ' coordinates')
 272         format('        of cell center j,k =',i4,',',
     .       i4,' have been unsuccessful...')
 372         format(
     .       '        will use  averages of finer level:',
     .       ' interpolation number',i4)
 472         format(
     .       '        for ALL points on this interface')
          else
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
             write(bou(nou(1),1),'('' program terminated in dynamic '',
     .             ''patching routines - see file '',a60)') grdmov
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),*)'         stopping...all attempts',
     .       ' to find generalized coordinates'
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),*)'         of cell midpoint j,k = ',
     .       j,1,' have been unsuccessful...check',titlptchgrd(1:len1)
          end if
          istop = 1
          iout  = 1
          return
        end if
      end if
c
c     search routine successful with current polynomial fit
c
      if(iok.eq.1) then
        etat1 = etat
        xiet1 = xiet
        l11 = l1
        ifit1 = ifit
      end if
c
c     try again with input polynomial fit
c
      if(ifit.ne.ifits) then
        iretry = 1
        ifit = ifits
        go to 18085
      end if
c
c     if second try with input value of ifit not successful, use the ifit 
c     value for which the search routine was successful
c
      if(iretry.eq.1) then
        if(iok.ne.1) then
          etat = etat1
          xiet = xiet1
          l1 = l11
          if (isklt1.gt.0) then
          nou(4) = min(nou(4)+1,ibufdim)
          write(bou(nou(4),4),*)'         iterations using original',
     .    ' fit  not successful at j,k= ',j,k21
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
c     block locations for k=k21 boundary points are stored in nblkj(j)
c
      nblkj(j) = l1 
      ifit     = ifits
      etac     = etat - 2.0
      fact     = 1.
c     
c     For some cases, boundaries eta=0 on "to"
c     side and "from" side are not supposed to be coincident. Check
c     to see if this is the case and, if so, skip over corrections to
c     make the edges coincident. 
c
      mcoin = mceta*(j22-j21)/100
      if (real(etac).gt.1.) then
      ncoin = ncoin+1
      if (ncoin.ge.mcoin) then
      if(isklt1.gt.0) then
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),*)'      quitting check of boundary...'
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),'(''      "to" and "from" eta=0 lines'',
     .  '' apparently not supposed to be coincident at'',
     .  '' this boundary'')')
      end if
      go to 670
      end if
      end if
c
      if (real(etat).lt.2.0) then
c
c     the following correction accounts for the widely varying cell 
c     sizes near the boundary of a "from" grid  after extension
c
        js    = xiet
        dis1  = (x1(js,3,l1)-x1(js,2,l1))**2
     .         +(y1(js,3,l1)-y1(js,2,l1))**2
     .         +(z1(js,3,l1)-z1(js,2,l1))**2
        dis0  = (x1(js,2,l1)-x1(js,1,l1))**2
     .         +(y1(js,2,l1)-y1(js,1,l1))**2
     .         +(z1(js,2,l1)-z1(js,1,l1))**2
        fact1 = sqrt(dis0/dis1)
        js    = js + 1
        dis1  = (x1(js,3,l1)-x1(js,2,l1))**2
     .         +(y1(js,3,l1)-y1(js,2,l1))**2
     .         +(z1(js,3,l1)-z1(js,2,l1))**2
        dis0  = (x1(js,2,l1)-x1(js,1,l1))**2
     .         +(y1(js,2,l1)-y1(js,1,l1))**2
     .         +(z1(js,2,l1)-z1(js,1,l1))**2
        fact2 =  sqrt(dis0/dis1)
        fact  = 0.5*( fact1 + fact2 )
        etac  = etac*fact
      end if
c
      temp(j) = etac
c
      do 1870 k=k21,k22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      if (mblkpt(ll).eq.l1) then
      seta2(j,k,1) = eta2(ll)
      if (real(eta2(ll)).lt.1.) seta2(j,k,1) = (eta2(ll)-1.)*fact+1.
      end if
 1870 continue
c
c     find km, the last k-point on the current j-line in the "to" grid that
c     lies in the same "from" block as the k=1 boundary point
c
      kmm(j) = 0
      do 2110 k=k21,k22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      if (mblkpt(ll).ne.nblkj(j)) go to 2111
      kmm(j) = k
 2110 continue
 2111 continue
      x7 = x8
      y7 = y8
      z7 = z8
 2000 continue
c
c     shearing correction
c
      if(ishear.lt.0) then
        kcorr=1
        go to 670
      end if
c
      beta   = 25.                      
      do 2100 j=j21,j22-1
      kmaxck = kkmax1(nblkj(j))
      km     = kmm(j)
      if (km.le.1) go to 2100
      do 2088 k=k21,km
      fact         = 1. - (float(km-k)-0.5)/float(km-k21)
      red          = exp(-beta*fact*fact)
      seta2(j,k,1) = seta2(j,k,1) - temp(j)*red
c
c     check to see if sheared generalized coordinates will remain in the 
c     legal range...if not try the arc length correction near the surface
c
c
      if (real(seta2(j,k,1)).lt.1. .or. 
     .    real(seta2(j,k,1)).gt.real(float(kmaxck-2))) then
      if(isklt1.gt.0) then
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),*)'        quitting eta shearing',
     .   ' correction at j,k=',j,k
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),*)'        eta2(j,k),temp(j)= ',
     .   real(seta2(j,k,1)),real(temp(j))
         nou(4) = min(nou(4)+1,ibufdim)
         write(bou(nou(4),4),*)'        will try arc length correction',
     .   ' near eta boundary '
      end if
      kcorr = 1
      go to 670
      end if
 2088 continue
c
      do 2166  k=k21,k22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      if (mblkpt(ll).eq.nblkj(j)) eta2(ll) = seta2(j,k,1)
 2166 continue
 2100 continue
  670 continue
c
c***************************************************************************
c      correct boundary values near xie=0
c***************************************************************************
c
      jcorr = 0
      ncoin = 0
c
      if (mcxie.eq.0) then 
        if(isklt1.gt.0) then
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),*) '      xie=0 boundaries not rendered', 
     .     ' coincident'
        end if
        go to 770
      end if
      if(ishear.ge.0) then
        if(isklt1.gt.0) then
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),*) '      xie=0 boundaries being',
     .     ' rendered coincident via shearing method'
        end if
      else
        if(isklt1.gt.0) then
           nou(4) = min(nou(4)+1,ibufdim)
           write(bou(nou(4),4),*) '      xie=0 boundaries being',
     .     ' rendered coincident via arc length method'
        end if
      end if
c
c     loop over all "to" cell on eta=0 boundary
c
      do 4000 k=k21,k22-1
c
c     compute edge midpoints of first layer of "to" grid cells
c     along the xie=0 boundary using quadratic least squares
c
      kcall = k
      jcall = j21
         call extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .              jcall,kcall,kl,kr,x7,y7,z7,icase,ifit)
      if(k. eq. k21) then
         call extra(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .              jcall,kcall,jl,jr,x5,y5,z5,icase,ifit)
      end if
      jcall = j21+1
         call extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .              jcall,kcall,kl,kr,x8,y8,z8,icase,ifit)
      jcall = j21
      kcall = k+1
         call extra(jdim1,kdim1,msub2,l2,x2,y2,z2,
     .              jcall,kcall,jl,jr,x6,y6,z6,icase,ifit)
c
c     compute normalized directed areas/unit normals of "to" cell
c
      if (itoss0 .eq. 0) then
         call direct(x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,
     .                   a1,a2,a3,imaxa,nou,bou,nbuf,ibufdim)
         ap(1) = a1
         ap(2) = a2
         ap(3) = a3
      end if
c
      ifits = ifit
      iretry = 0
      icount = 1
19085 continue
c
c     compute midpoint of "to" cell on xie=0 boundary, consistent with ifit 
c
c     bilinear or linear in eta
      if(ifit.eq.1 .or. ifit.eq.3) then
        xc = 0.5*(x2(j21,k+1,l2)+x2(j21,k,l2))
        yc = 0.5*(y2(j21,k+1,l2)+y2(j21,k,l2))
        zc = 0.5*(z2(j21,k+1,l2)+z2(j21,k,l2))
      else
c     biquadratic or quadratic in eta
        xc = x7
        yc = y7
        zc = z7
      end if
c
c     call trace(2,icheck,jcall,k,ifit,xc,yc,zc)
c
c     set starting point for search routine...use solution at cell center
c     of current cell
c
c     xiet,etat are generalized coordinates in expanded "from" grid
c     xie2,eta2 are generalized coordinates in original "from" grid
c
c
      ll = (j22-j21)*(k-k21) + 1
      xiet   = xie2(ll) + 1.
      etat   = eta2(ll) + 1.
      l1     = mblkpt(ll)
c
c     search over "from" blocks to find xie,eta associated with "to"
c     cell midpoint xc,yc,zc on xie=0
c
      call topol(jdim1,kdim1,msub1,jjmax1,kkmax1,lmax1,l1,x1,y1,z1,
     .           x1mid,y1mid,z1mid,x1mide,y1mide,z1mide,limit0,xc,yc,zc,
     .           xiet,etat,jimage,kimage,ifit,itmax,igap,iok,lout,ic0,
     .           itoss0,1,k,iself,xif1,xif2,etf1,etf2,nou,bou,nbuf,
     .           ibufdim,myid)
c
c     search routine unsuccessful...try alternative polynomial fit
c
      if (iretry.eq.0 .and. iok.ne.1) then
        if (icount.lt.4) then
           call newfit(ifits,ifit,icount)
           icount = icount + 1
           go to 19085
        else
c
c         don't shear orphan points
c
          if(iorph.gt.0) then
            if(isklt1.gt.0) then
               nou(4) = min(nou(4)+1,ibufdim)
               write(bou(nou(4),4),*)' quitting boundary correction',
     .        ' because points have been flagged as orphans'
            end if
            return
          end if
c
          if (ifiner.gt.0) then
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),173) 
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),273) 1,k
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),373) ifiner
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),473)
 173         format('        all attempts to find generalized',
     .       ' coordinates')
 273         format('        of cell center j,k =',i4,',',
     .       i4,' have been unsuccessful...')
 373         format(
     .       '        will use  averages of finer level:',
     .       ' interpolation number',i4)
 473         format(
     .       '        for ALL points on this interface')
          else
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
             write(bou(nou(1),1),'('' program terminated in dynamic '',
     .             ''patching routines - see file '',a60)') grdmov
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),*)'         stopping...all attempts',
     .       ' to find generalized coordinates'
             nou(4) = min(nou(4)+1,ibufdim)
             write(bou(nou(4),4),*)'         of cell center j,k = ',
     .       ' 1,k,have been unsuccessful...check',titlptchgrd(1:len1)
          end if
          istop = 1
          iout  = 1
          return
        end if
      end if
c
c     search routine successful with current polynomial fit
c
      if(iok.eq.1) then
        xiet1 = xiet
        etat1 = etat
        l11 = l1
        ifit1 = ifit
      end if
c
c     try again with input polynomial fit
c
      if(ifit.ne.ifits) then
        iretry = 1
        ifit = ifits
        go to 19085
      end if
c
c     if second try with input value of ifit not successful, use the ifit 
c     value for which the search routine was successful
c
      if(iretry.eq.1) then
        if(iok.ne.1) then
          xiet = xiet1
          etat = etat1
          l1 =l11
          if(isklt1.gt.0) then
          nou(4) = min(nou(4)+1,ibufdim)
          write(bou(nou(4),4),*)'         iterations using original',
     .    ' fit not successful at j,k= ',j21,k
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
c     block locations for j=j21 boundary points stored in nblkk(k)
c
      nblkk(k) = l1
      ifit     = ifits
      xiec     = xiet - 2.0
      fact     = 1.
c     
c     For some cases, boundaries xie=0 on "to"
c     side and "from" side are not supposed to be coincident. Check
c     to see if this is the case and, if so, skip over corrections to
c     make the edges coincident 
c
      mcoin = mcxie*(k22-k21)/100
      if (real(xiec).gt.1.) then
      ncoin = ncoin+1
      if (ncoin.ge.mcoin) then
      if(isklt1.gt.0) then
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),*)'      quitting check of boundary...'
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),'(''      "to" and "from" xie=0 lines'',
     .  '' apparently not supposed to be coincident at'',
     .  '' this boundary'')')
      end if
      go to 770
      end if
      end if
c
      if (real(xiet).lt.2.0) then
c
c     the following correction accounts for the widely varying cell 
c     sizes near the boundary of a "from" grid  after extension
c
        ks    = etat
        dis1  = (x1(3,ks,l1)-x1(2,ks,l1))**2
     .         +(y1(3,ks,l1)-y1(2,ks,l1))**2
     .         +(z1(3,ks,l1)-z1(2,ks,l1))**2
        dis0  = (x1(2,ks,l1)-x1(1,ks,l1))**2
     .         +(y1(2,ks,l1)-y1(1,ks,l1))**2
     .         +(z1(2,ks,l1)-z1(1,ks,l1))**2
        fact1 = sqrt(dis0/dis1)
        ks    = ks + 1
        dis1  = (x1(3,ks,l1)-x1(2,ks,l1))**2
     .         +(y1(3,ks,l1)-y1(2,ks,l1))**2
     .         +(z1(3,ks,l1)-z1(2,ks,l1))**2
        dis0  = (x1(2,ks,l1)-x1(1,ks,l1))**2
     .         +(y1(2,ks,l1)-y1(1,ks,l1))**2
     .         +(z1(2,ks,l1)-z1(1,ks,l1))**2
        fact2 =  sqrt(dis0/dis1)
        fact  = 0.5*( fact1 + fact2 )
        xiec  = xiec*fact
      end if
c
      temp(k) = xiec
c
      do 3870 j=j21,j22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      if (mblkpt(ll).eq.l1) then
         sxie2(j,k,1) = xie2(ll)
         if (real(xie2(ll)).lt.j21) sxie2(j,k,1) = (xie2(ll)-1.)*fact+1.
      end if
 3870 continue
c
c     find jm, the last j-point on the current k-line in the "to" grid that
c     lies in the same "from" block as the j=1 boundary point
c
      jmm(k) = 0
      do 4110 j=j21,j22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      if (mblkpt(ll).ne.nblkk(k)) go to 4111
      jmm(k) = j
 4110 continue
 4111 continue
      x5 = x6
      y5 = y6
      z5 = z6
 4000 continue
c
c     shearing correction
c
      if(ishear.lt.0) then
        jcorr=1
        go to 770
      end if
c
      beta = 25.                      
      do 4100 k=k21,k22-1
      jmaxck = jjmax1(nblkk(k))
      jm     = jmm(k)
      if (jm.le.j21) go to 4100
      do 4088 j=j21,jm
      fact         = 1. - (float(jm-j)-0.5)/float(jm-j21)
      red          = exp(-beta*fact*fact)
      sxie2(j,k,1) = sxie2(j,k,1) - temp(k)*red
c
c     check to see if sheared generalized coordinates will remain in the 
c     legal range...if not try the arc length correction near the surface
c
      if (real(sxie2(j,k,1)).lt.1. .or. 
     .    real(sxie2(j,k,1)).gt.float(jmaxck-2)) then
      if(isklt1.gt.0) then
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),*) '        quitting xie shearing',
     .  ' correction at j,k=',j,k
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),*)'        xie2(j,k),temp(k)= ',
     .  real(sxie2(j,k,1)),real(temp(k))
        nou(4) = min(nou(4)+1,ibufdim)
        write(bou(nou(4),4),*) '        will try arc length',
     .  ' correction near xie boundary'
      end if
      jcorr = 1
      go to 770
      end if
 4088 continue
      do 4166 j=j21,j22-1
      ll = (j22-j21)*(k-k21) + (j-j21+1)
      if (mblkpt(ll).eq.nblkk(k)) xie2(ll) = sxie2(j,k,1)
 4166 continue
 4100 continue
  770 continue
      return
      end
