      program XINTOUT_to_ovrlp
c
c***********************************************************************
c     Purpose: Converts the XINTOUT overset grid interpolation file
c     from PEGSUS to the ovrlp.bin file used by CFL3D
c***********************************************************************
c
      parameter (mxbl=9999)
c
      dimension jdim(mxbl),kdim(mxbl),ldim(mxbl),iisptr(mxbl),
     .          iieptr(mxbl)
c
      character*80 peggrd
c
c     query for input
c
c     ipeg: flag for pegsus version used
c
 10   continue
      write(6,*)
      write(6,'(''enter 4 if pegsus version 4.X was used'')')
      write(6,'(''enter 5 if pegsus version 5.X was used'')')
      read(5,*) ipeg
      if (ipeg.eq.4 .or. ipeg.eq.5) then
      else
          write(6,'(''incorrect response for the PEGSUS version '',
     .    ''- try again'')')
          go to 10
      end if
c
c     need a file with grid data for plot3d diagnostic files
c
      write(6,*)
      if (ipeg .eq. 4) then
         write(6,'(''enter the name of the COMPOUT file created by '',
     .             ''PEGSUS'')')
         write(6,'(''(unless you have renamed it, enter COMPOUT)'')')
      else
         write(6,'(''enter the name of the cell-center plot3d input '',
     .             ''to PEGSUS'')')
      end if
      read(5,'(a80)') peggrd
c
c     i2d > 0: convert to 2d grid data for cfl3d (2 planes/1 cell)
c
      write(6,*)
      write(6,'(''enter 0 to leave overset data as 3d'')')
      write(6,'(''enter 1 to convert overset data to 2d '',
     .          ''(2 i-planes)'')')
      read(5,*) i2d
c
c     allow for reordering jkl <---> klj if necessary
c
      write(6,*)
      write(6,'(''enter 0 to use the maggie/CFL3D-ijk index '',
     .          ''definition'')')
      write(6,'(''enter 1 to use the pegsus/OVERFLOW-jkl index '',
     .          ''definition'')')
      read(5,*) iswap
c 
c     set debugging flags:
c
c     iogrd > 0: write out a cell-center grid, augmented with
c     iblank data: 0 for holes, 1 for field points, and -ig for
c     fringe/outer boundary points that get interpolated from
c     mesh ig (unformatted, written to "source_chk.unf")
c
      iogrd = 1
c
c     ioorph > 0: write out a list of phantom points in plot3d /mg/unf
c     can be useful for quickly finding location of orphan points when
c     used in conjunction with  the plot3d file generated if iogrd = 1
c     format to "orphan_chk.unf"
c
      ioorph = 1
c
c     nppblk is the max Number of orphan Points Per BLocK in the output file
c     (required to meet the 1-D array limits of some plotting packages)
c
      nppblk = 1660
c
c     iverbose > 0: write out all diagnostic info for target points
c     that get removed
c     (formatted, written to "diagnos_info.dat")
c
      iverbose = 0
c
c     iosten > 0: output the final stencil data for each target point
c     (formatted, written to "stencil_info.dat")
c
      iosten = 0
c
c     open files
c
      open(unit=1,file=peggrd,form='unformatted',
     .     status='old')
      open(unit=2,file='XINTOUT',form='unformatted',
     .     status='old')
      open(unit=3,file='ovrlp.bin',form='unformatted',
     .     status='unknown')
      open(unit=7,file='stencil_error.dat',form='formatted',
     .        status='unknown')
      open(unit=8,file='ORPHAN',form='formatted',
     .     status='unknown')
      if (iogrd .gt. 0) then
         open(unit=4,file='source_chk.unf',form='unformatted',
     .        status='unknown')
      end if
      if (ioorph .gt. 0) then
         open(unit=14,file='orphan_chk.unf',form='unformatted',
     .        status='unknown')
      end if
      if (iverbose .gt. 0) then
         open(unit=76,file='diagnos_info.dat',form='formatted',
     .        status='unknown')
      end if
      if (iosten .gt. 0) then
         open(unit=80,file='stencil_info.dat',form='formatted',
     .        status='unknown')
      end if
c
c     obtain required array dimensions
c
      maxbl  = 1
      lmax   = 1
      jmax   = 1
      kmax   = 1
      ibsptr = 1
      ibeptr = 1
      if (ipeg .ne. 4) then
         read(1) ngrid
         read(1) (jdim(ig),kdim(ig),ldim(ig),ig=1,ngrid)
      end if
      igend = 9999
      if (ipeg .eq.5) igend = ngrid
      do ig = 1,igend
         if (ipeg .eq. 4) then
            read(1,end=999) nm,jd,kd,ld
         else
            jd = jdim(ig)
            kd = kdim(ig)
            ld = ldim(ig)
         end if
         jdim(ig) = jd - 1
         kdim(ig) = kd - 1
         ldim(ig) = ld - 1
         write(6,*)
         write(6,'(''COMPOUT grid '',i3,
     .   '': dimensions j x k x l = '',3i7)') ig,jd,kd,ld
         maxbl = max(maxbl,ig)
         lmax  = max(lmax,ld)
         jmax  = max(jmax,jd)
         kmax  = max(kmax,kd)
         if (ipeg .eq. 4) then
            read(1) (xx,ll=1,jd*kd*ld),
     .              (yy,ll=1,jd*kd*ld),
     .              (zz,ll=1,jd*kd*ld)
         end if
         read(2) ibpnts,iipnts,iieptr(ig),iisptr(ig)
         write(6,'(''XINTOUT grid '',i3,
     .   '':        ibpnts,iipnts = '',2i7)') ig,ibpnts,iipnts
         is = iisptr(ig)
         ie = iieptr(ig)
         read(2) (ji,i=is,ie), (ki,i=is,ie),
     .           (li,i=is,ie), (dj,i=is,ie),
     .           (dk,i=is,ie), (dl,i=is,ie)
         iitoti = ie
         if (ig .gt. 1) ibsptr = ibeptr + 1
         is          = ibsptr
         ie          = is + ibpnts - 1
         ibeptr = ie
         read(2) (jb,i=is,ie), (kb,i=is,ie),
     .           (lb,i=is,ie),(ibc,i=is,ie)
         read(2) (((ib,j=1,jd),k=1,kd),l=1,ld)
         iitotb = ie
      end do
 999  continue
c
      ngrid = maxbl
      iitot = max(iitoti,iitotb)
c
      write(6,*)
      write(6,'(''required array sizes: maxbl = '',i6)') maxbl
      write(6,'(''                       lmax = '',i6)') lmax
      write(6,'(''                       jmax = '',i6)') jmax
      write(6,'(''                       kmax = '',i6)') kmax
      write(6,'(''                      iitot = '',i6)') iitot
c
      rewind(1)
      rewind(2)

      call convert(i2d,iitot,maxbl,jmax,kmax,lmax,jdim,kdim,ldim,
     .             iisptr,iieptr,iswap,iogrd,iverbose,iosten,
     .             ioorph,nppblk,ipeg)
c
      write(6,*)
      write(6,'(''conversion completed'')')
c
      stop
      end
c
      subroutine convert(i2d,iitot,maxbl,jmax,kmax,lmax,jdim,kdim,ldim,
     .                   iisptr,iieptr,iswap,iogrd,iverbose,iosten,
     .                   ioorph,nppblk,ipeg)
c
      dimension ldim(maxbl),jdim(maxbl),kdim(maxbl),iisptr(maxbl),
     .          iieptr(maxbl)
c
      dimension x(jmax,kmax,lmax),y(jmax,kmax,lmax),
     .          z(jmax,kmax,lmax)
c
c     for input pegsus data
c
      dimension ibpnts(maxbl),iipnts(maxbl),ji(iitot),ki(iitot),
     .          li(iitot),jb(iitot),kb(iitot),lb(iitot),ibc(iitot),
     .          iblank(jmax,kmax,lmax)
      dimension dj(iitot),dk(iitot),dl(iitot)
      dimension ibeptr(maxbl),ibsptr(maxbl),iiblk(iitot),
     .          ibblk(iitot)
      dimension lborph(iitot),jborph(iitot),kborph(iitot)
      dimension igorph(iitot),norph(maxbl)
      dimension intpts(maxbl,4)
c
      pointer (ip_x,x),(ip_y,y),(ip_z,z)
      pointer (ip_ibpnts,ibpnts),(ip_iipnts,iipnts),
     .        (ip_ji,ji),(ip_ki,ki),(ip_li,li),(ip_jb,jb),
     .        (ip_kb,kb),(ip_lb,lb),(ip_ibc,ibc),
     .        (ip_iblank,iblank),(ip_dj,dj),(ip_dk,dk),
     .        (ip_dl,dl),(ip_ibeptr,ibeptr),(ip_ibsptr,ibsptr),
     .        (ip_iiblk,iiblk),(ip_ibblk,ibblk)
      pointer (ip_lborph,lborph),(ip_jborph,jborph),(ip_kborph,kborph),
     .        (ip_igorph,igorph),(ip_norph,norph),(ip_intpts,intpts)
c
      ngrid   = maxbl
c
c     allocate memory
c
      icall  = 0
      memuse = 0
      call umalloc_r(ip_x,lmax*jmax*kmax,0,icall,memuse)
      call umalloc_r(ip_y,lmax*jmax*kmax,0,icall,memuse)
      call umalloc_r(ip_z,lmax*jmax*kmax,0,icall,memuse)
      call umalloc_r(ip_ibpnts,maxbl,1,icall,memuse)
      call umalloc_r(ip_iipnts,maxbl,1,icall,memuse)
      call umalloc_r(ip_ji,iitot,1,icall,memuse)
      call umalloc_r(ip_ki,iitot,1,icall,memuse)
      call umalloc_r(ip_li,iitot,1,icall,memuse)
      call umalloc_r(ip_jb,iitot,1,icall,memuse)
      call umalloc_r(ip_kb,iitot,1,icall,memuse)
      call umalloc_r(ip_lb,iitot,1,icall,memuse)
      call umalloc_r(ip_ibc,iitot,1,icall,memuse)
      call umalloc_r(ip_iblank,lmax*jmax*kmax,1,icall,memuse)
      call umalloc_r(ip_dj,iitot,0,icall,memuse)
      call umalloc_r(ip_dk,iitot,0,icall,memuse)
      call umalloc_r(ip_dl,iitot,0,icall,memuse)
      call umalloc_r(ip_ibeptr,maxbl,1,icall,memuse)
      call umalloc_r(ip_ibsptr,maxbl,1,icall,memuse)
      call umalloc_r(ip_iiblk,iitot,1,icall,memuse)
      call umalloc_r(ip_ibblk,iitot,1,icall,memuse)
      call umalloc_r(ip_lborph,iitot,1,icall,memuse)
      call umalloc_r(ip_jborph,iitot,1,icall,memuse)
      call umalloc_r(ip_kborph,iitot,1,icall,memuse)
      call umalloc_r(ip_igorph,iitot,1,icall,memuse)
      call umalloc_r(ip_norph,maxbl,1,icall,memuse)
      call umalloc_r(ip_intpts,maxbl*4,1,icall,memuse)
c
c     notation: arrays jdim(n),kdim(n),ldim(n) are taken as cfl3d
c     GRID dimensions for zone n. the corresponding cfl3d CELL CENTER
c     dimensions are thus jdim(n)-1,kdim(n)-1,ldim(n)-1.
c
c     the augmented cell center dimensions used by pegsus are 
c     jdim(n)+1,kdim(n)+1,ldim(n)+1 (cell centers plus one rind
c     layer in each of the 6 block faces)
c
      do ig=1,ngrid
         norph(ig) = 0
         is = iisptr(ig)
         ie = iieptr(ig)
         do i=is,ie
            iiblk(i) = ig
         end do
      end do
c
      ibsptr(1) = 1
      ibeptr(1) = 1
c
      if (ipeg .eq. 5) then
         read(1)
         read(1)
      end if
      do ig = 1,ngrid
c
         if (ipeg .eq. 4) then
            read(1) nm,jd0,kd0,ld0
         else
            jd0 = jdim(ig) + 1
            kd0 = kdim(ig) + 1
            ld0 = ldim(ig) +1
         end if
         read(1) (((x(j,k,l),j=1,jd0),k=1,kd0),l=1,ld0),
     .           (((y(j,k,l),j=1,jd0),k=1,kd0),l=1,ld0),
     .           (((z(j,k,l),j=1,jd0),k=1,kd0),l=1,ld0)
c
         read(2) ibpnts(ig),iipnts(ig),iieptr(ig),iisptr(ig)
c
c        read source points in zone ig (provide data to another zone)
c
         is = iisptr(ig)
         ie = iieptr(ig)
         read(2) (ji(i),i=is,ie), (ki(i),i=is,ie),
     .           (li(i),i=is,ie), (dj(i),i=is,ie),
     .           (dk(i),i=is,ie), (dl(i),i=is,ie)
c
c        read target points in zone ig (receive data from another zone)
c
         if (ig .gt. 1) ibsptr(ig) = ibeptr(ig-1) + 1
         is         = ibsptr(ig)
         ie         = is + ibpnts(ig) - 1
         ibeptr(ig) = ie
         read(2) (jb(i),i=is,ie),(kb(i),i=is,ie),
     .           (lb(i),i=is,ie),(ibc(i),i=is,ie)
         do i=is,ie
            ibblk(i) = ig
         end do
c
         read(2) (((iblank(j,k,l),j=1,jd0),k=1,kd0),l=1,ld0)
c
c        read in orphan points, if any
c
         if (ig .eq. 1) then
            norpht0 = 0
            do n=1,iitot
               read(8,*,end=888) igorph(n),jborph(n),kborph(n),
     .                           lborph(n)
               norpht0 = norpht0 + 1
               norph(igorph(n)) = norph(igorph(n)) + 1
            end do
 888        continue
         end if
c
c        check initial stencils to make sure pegsus hasn't output
c        data it shouldn't have
c
         call chksten0(ig,iisptr,iieptr,ibsptr,ibeptr,jdim,kdim,
     .                ldim,maxbl,jb,kb,lb,ji,ki,li,dj,dk,dl,
     .                iitot)
c
c        in cell center mode, pegsus adds one extra layer of cell
c        centers on each zonal boundary. start by shifting indicies
c        in the pegsus data by -1, so that the indicies in the pegsus
c        data start at 0. then index 1 in the pegsus data corresponds
c        to the first cfl3d cell center. stencil points with indicies
c        0,jdim+1,kdim+1,ldim+1 will be used to access cell-face
c        center boundary data in cfl3d. this also allows backward 
c        compatability with maggie interpolation data.
c
         do i=iisptr(ig),iieptr(ig)
            ji(i) = ji(i) - 1
            ki(i) = ki(i) - 1
            li(i) = li(i) - 1
         end do
         do i=ibsptr(ig),ibeptr(ig)
            jb(i) = jb(i) - 1
            kb(i) = kb(i) - 1
            lb(i) = lb(i) - 1
         end do
         if (ig .eq. 1) then
            if (norpht0 .gt.0) then
               do n=1,norpht0
                  jborph(n) = jborph(n) - 1
                  kborph(n) = kborph(n) - 1
                  lborph(n) = lborph(n) - 1
               end do
            end if
         end if
c
c        remove target cells whose index falls out of the range
c        of cfl3d cell centers - only target cells in the range
c        of unaugmented cell centers get updated by cfl3d
c
         if (ig .eq. 1) then
            ii = 0
         else
            ii = ibeptr(ig-1)
         end if
         do i=ibsptr(ig),ibeptr(ig)
            jd1 = jdim(ibblk(i)) - 1
            kd1 = kdim(ibblk(i)) - 1
            ld1 = ldim(ibblk(i)) - 1
            if (i2d .gt. 0) ld1 = 1
            if ((jb(i).ge.1 .and. jb(i).le.jd1) .and.
     .          (kb(i).ge.1 .and. kb(i).le.kd1) .and.
     .          (lb(i).ge.1 .and. lb(i).le.ld1)) then
                ii        = ii+1
                jb(ii)    = jb(i)
                kb(ii)    = kb(i)
                lb(ii)    = lb(i)
                ibc(ii)   = ibc(i)
                ibblk(ii) = ibblk(i)
            else if (iverbose.gt.0) then
                write(76,*)
                write(76,'(''removing target point in block'',i4)')
     .          ibblk(i)
                write(76,'(''  list index  '',i5)') i
                write(76,'(''  jd1,kd1,ld1 '',3i4)') jd1,kd1,ld1
                write(76,'(''  jb0,kb0,lb0 '',3i4)') jb(i),kb(i),lb(i)
            end if
         end do
         if (ig .gt. 1) ibsptr(ig) = ibeptr(ig-1) + 1
         ibeptr(ig) = ii
         ibpnts(ig) = ibeptr(ig) - ibsptr(ig) + 1
         intpts(ig,1) = ibpnts(ig)
         intpts(ig,2) = 0
         intpts(ig,3) = 0
         intpts(ig,4) = 0
c
c        remove any orphan points that are outside the range of cfl3d
c        cell centers...these points won't get updated in cfl3d anyway
c
         if (ig .eq. 1) then
            norpht = 0
            do igg=1,ngrid
               norph(igg) = 0
            end do
            if (norpht0 .gt.0) then
               do n=1,norpht0
                  jd1 = jdim(igorph(n)) - 1
                  kd1 = kdim(igorph(n)) - 1
                  ld1 = ldim(igorph(n)) - 1
                  if (i2d .gt. 0) ld1 = 1
                  if ((jborph(n).ge.1 .and. jborph(n).le.jd1) .and.
     .                (kborph(n).ge.1 .and. kborph(n).le.kd1) .and.
     .                (lborph(n).ge.1 .and. lborph(n).le.ld1)) then
                      norpht = norpht + 1
                      jborph(norpht)        = jborph(n)
                      kborph(norpht)        = kborph(n)
                      lborph(norpht)        = lborph(n)
                      igorph(norpht)        = igorph(n)
                      norph(igorph(norpht)) = norph(igorph(norpht)) + 1
                  else if (iverbose.gt.0) then
                      write(76,*)
                      write(76,'(''removing orphan point in block'',
     .                i4)')igorph(n)
                      write(76,'(''     orphan list index '',i5)') n
                      write(76,'(''           jd1,kd1,ld1 '',3i4)')
     .                jd1,kd1,ld1
                      write(76,'(''  jborph,kborph,lborph '',3i4)')
     .                jborph(n),kborph(n),lborph(n)
                  end if
               end do
            end if
         end if
c
c        output the updated lists in maggie format
c
         jd = jdim(ig) - 1
         kd = kdim(ig) - 1
         ld = ldim(ig) - 1
         if (i2d.gt.0) ld = 1
         if (iswap .eq. 0) then
            write(3) jd,kd,ld
            write(3) ibpnts(ig),(intpts(ig,ll),ll=1,4),iipnts(ig),
     .               iieptr(ig),iisptr(ig)
            is = iisptr(ig)
            ie = iieptr(ig)
            write(3) (ji(i),ki(i),li(i),dj(i),dk(i),dl(i),i=is,ie)
            is = ibsptr(ig)
            ie = ibeptr(ig)
            write(3) (jb(i),kb(i),lb(i),ibc(i),i=is,ie)
            write(3) (((real(iblank(j,k,l)),j=2,jd+1),k=2,kd+1),
     .               l=2,ld+1)
         else
            write(3) kd,ld,jd
            write(3) ibpnts(ig),(intpts(ig,ll),ll=1,4),iipnts(ig),
     .               iieptr(ig),iisptr(ig)
            is = iisptr(ig)
            ie = iieptr(ig)
            write(3) (ki(i),li(i),ji(i),dk(i),dl(i),dj(i),i=is,ie)
            is = ibsptr(ig)
            ie = ibeptr(ig)
            write(3) (kb(i),lb(i),jb(i),ibc(i),i=is,ie)
            write(3) (((real(iblank(j,k,l)),k=2,kd+1),l=2,ld+1),
     .      j=2,jd+1)
         end if
c
c        output to diagnostic plot3d files if desired
c 
         if (iogrd .gt. 0) then
            if (ig.eq.1) then
               write(4) ngrid
               if (i2d .eq. 0) then
                  write(4) (ldim(n)-1,jdim(n)-1,kdim(n)-1,n=1,ngrid)
               else
                  write(4) (jdim(n)-1,kdim(n)-1,n=1,ngrid)
               end if
            end if
c
c           modify iblank using stencil info - use to visualize the
c           source point of each target point as suggested by James 
c           Hager. note: must shift jb,kb,lb back to original values
c           in the iblank array
c
            do i=ibsptr(ig),ibeptr(ig)
               ii = ibc(i)
               iblank(jb(i)+1,kb(i)+1,lb(i)+1) = -iiblk(ii)
            enddo
c
c           set orphan points to blank=+2 to distinguish them from 
c           field, hole, or interpolated points
c
            if (norpht .gt. 0) then
               do n=1,norpht
                  if (igorph(n) .eq. ig) then
                     iblank(jborph(n)+1,kborph(n)+1,lborph(n)+1) = +2
                  end if
               end do
            end if
c
            lunit = 4
            call wgrid(jmax,kmax,lmax,ld,jd,kd,x,y,z,iblank,
     .                 lunit,i2d)
         end if
c
      end do
c
c     check final stencil lists to make sure nothing has escaped that
c     shouldn't have!
c
      ierr = 0
      call chksten(ngrid,ibsptr,ibeptr,jdim,kdim,ldim,maxbl,
     .             jb,kb,lb,ji,ki,li,dj,dk,dl,ibc,ibblk,iiblk,
     .             iitot,norph,i2d,ierr)
c
c     print orphan point summary
c
      write(6,*)
      write(6,'(''number of target points with valid '',
     .''stencils:  '',i6)')ibeptr(ngrid)
      write(6,'(''number of original orphan points '',
     .''removed:     '',i6)') norpht0 - norpht
      write(6,'(''number of orphan points remaining:'',
     .''            '',i6)') norpht
       write(6,'(''percentage of orphan points:'',
     .''                  '',f6.2)')
     .float(norpht)/ibeptr(ngrid)*100
c
c     output stencils if desired
c
      if (iosten .gt. 0) then
         iunit = 80
         call wstencil(iitot,maxbl,jb,kb,lb,ji,ki,li,dj,dk,dl,
     .                 iiblk,ibsptr,ibeptr,ibc,iunit,norpht,
     .                 jborph,kborph,lborph,igorph,ngrid)
      end if
c
      write(6,*)
      if (iogrd .gt. 0) then
         if (i2d .eq. 0) then
            write(6,'(''plot3d file (mg/unf/iblank) written to'',
     .      '' "source_chk.unf"'')')
         else
            write(6,'(''plot3d file (mg/unf/iblank/2d) written to'',
     .      '' "source_chk.unf"'')')
         end if
      end if
      if (iverbose .gt. 0) then
         write(6,'(''diagnostic info file written to'',
     .   '' "diagnos_info.dat"'')')
      end if
      if (iosten .gt. 0) then
         write(6,'(''stencil info file written to'',
     .   '' "stencil_info.dat"'')')
      end if
      if (ierr .gt. 0) then
         write(6,'(''WARNING: errors detected in final stencil list'',
     .   '' see "stencil_error.dat"'')')
         write(6,'(''   number of errors = '',i6)') ierr
      else
         close(7,status='delete')
      end if
      if (ioorph .gt. 0) then
         if (i2d .eq. 0) then
            write(6,'(''plot3d file (mg/unf) written to'',
     .      '' "orphan_chk.unf"'')')
         else
            write(6,'(''plot3d file (mg/unf/2d) written to'',
     .      '' "orphan_chk.unf"'')')
         end if
      end if
c
c     output orphan point locations in a plot3d file, if desired
c
      if (ioorph .gt. 0 .and. norpht .gt. 0) then
         lunit = 14
         lomax = 0
         nbmax = 0
         do ig=1,ngrid
            lomax = max(lomax,norph(ig))
            if(norph(ig).ne.0) then
               nout1 = norph(ig) /nppblk
               noutr = norph(ig) -nout1*nppblk
               if(noutr.eq.0) then
                  ioutd = 0
               else
                  ioutd = 1
               endif
               do iout=1,nout1+ioutd
                  nbmax = nbmax +1
               end do
            end if
         end do
         call worph(jmax,kmax,lmax,lomax,nbmax,x,y,z,ngrid,norpht,
     .              lborph,jborph,kborph,igorph,norph,lunit,nppblk,
     .              i2d,ipeg)

      end if

c     free memory
c
      ifree = 1
      if (ifree.gt.0) then
         call ufree_r(ip_x)
         call ufree_r(ip_y)
         call ufree_r(ip_z)
         call ufree_r(ip_ibpnts)
         call ufree_r(ip_iipnts)
         call ufree_r(ip_ji)
         call ufree_r(ip_ki)
         call ufree_r(ip_li)
         call ufree_r(ip_jb)
         call ufree_r(ip_kb)
         call ufree_r(ip_lb)
         call ufree_r(ip_ibc)
         call ufree_r(ip_iblank)
         call ufree_r(ip_dj)
         call ufree_r(ip_dk)
         call ufree_r(ip_dl)
         call ufree_r(ip_ibeptr)
         call ufree_r(ip_ibsptr)
         call ufree_r(ip_iiblk)
         call ufree_r(ip_ibblk)
         call ufree_r(ip_lborph)
         call ufree_r(ip_jborph)
         call ufree_r(ip_kborph)
         call ufree_r(ip_igorph)
         call ufree_r(ip_norph)
         call ufree_r(ip_intpts)
      end if
c
      return
      end
c
      subroutine chksten0(ig,iisptr,iieptr,ibsptr,ibeptr,jdim,kdim,
     .                    ldim,maxbl,jb,kb,lb,ji,ki,li,dj,dk,dl,
     .                    iitot)
c
c     checks unmodified pegsus stencils for illegal index bounds
c
      dimension ji(iitot),ki(iitot),li(iitot),jb(iitot),kb(iitot),
     .          lb(iitot),dj(iitot),dk(iitot),dl(iitot)
      dimension ibeptr(maxbl),ibsptr(maxbl),iieptr(maxbl),iisptr(maxbl)
      dimension jdim(maxbl),kdim(maxbl),ldim(maxbl)
c
c     check indicies of target point
c
      jbd1 = jdim(ig) + 2
      kbd1 = kdim(ig) + 2
      lbd1 = ldim(ig) + 2
c
      do i=ibsptr(ig),ibeptr(ig)
         if (jb(i).lt.1 .or. jb(i).gt.jbd1) then
            write(6,'(''target cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') i
            write(6,'(''     jb,   kb,   lb = '',3i4)') jb(i),kb(i),
     .      lb(i)
            write(6,'(''jdim1, kdim1, ldim1 =  '',3i4)') jbd1,kbd1,lbd1
         end if
         if (kb(i).lt.1 .or. kb(i).gt.kbd1) then
            write(6,'(''target cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') i
            write(6,'(''     jb,   kb,   lb = '',3i4)') jb(i),kb(i),
     .      lb(i)
            write(6,'(''jdim1, kdim1, ldim1 =  '',3i4)') jbd1,kbd1,lbd1
         end if
         if (lb(i).lt.1 .or. lb(i).gt.lbd1) then
            write(6,'(''target cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') i
            write(6,'(''     jb,   kb,   lb = '',3i4)') jb(i),kb(i),
     .      lb(i)
            write(6,'(''jdim1, kdim1, ldim1 = '',3i4)') jbd1,kbd1,lbd1
         end if
      end do
c
c     check indicies of source point
c
      jid1 = jdim(ig) + 2
      kid1 = kdim(ig) + 2
      lid1 = ldim(ig) + 2
c
      do i=iisptr(ig),iieptr(ig)
         if (ji(i).lt.1 .or. ji(i).gt.jid1) then
            write(6,'(''source cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') i
            write(6,'(''     ji,   ki,   li = '',3i4)') ji(i),
     .      ki(i),li(i)
            write(6,'(''jdim1, kdim1, ldim1 =  '',3i4)') jid1,kid1,lid1
            write(6,'(''   dj = '',3f8.3)')  dj(i)
         end if
         if (ki(i).lt.1 .or. ki(i).gt.kid1) then
            write(6,'(''source cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') i
            write(6,'(''     ji,   ki,   li = '',3i4)') ji(i),
     .      ki(i),li(i)
            write(6,'(''jdim1, kdim1, ldim1 =  '',3i4)') jid1,kid1,lid1
            write(6,'(''   dk = '',3f8.3)') dk(i)
         end if
         if (li(i).lt.1 .or. li(i).gt.lid1) then
            write(6,'(''source cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') i
            write(6,'(''     ji,   ki,   li = '',3i4)') ji(i),
     .      ki(i),li(i)
            write(6,'(''jdim1, kdim1, ldim1 =  '',3i4)') jid1,kid1,lid1 
            write(6,'(''   dl = '',3f8.3)') dl(i)
         end if
      end do
c
      return
      end
c
      subroutine chksten(ngrid,ibsptr,ibeptr,jdim,kdim,ldim,maxbl,
     .                   jb,kb,lb,ji,ki,li,dj,dk,dl,ibc,ibblk,iiblk,
     .                   iitot,norph,i2d,ierr)
c
c     checks final pegsus stencils for illegal index bounds; also
c     prints stencil summary
c
      dimension ji(iitot),ki(iitot),li(iitot),jb(iitot),kb(iitot),
     .          lb(iitot),ibc(iitot),dj(iitot),dk(iitot),dl(iitot)
      dimension ibeptr(maxbl),ibsptr(maxbl),iiblk(iitot),
     .          ibblk(iitot),norph(maxbl)
      dimension jdim(maxbl),kdim(maxbl),ldim(maxbl)
c
      do ig=1,ngrid
c
c     check indicies of target point
c
      jbd1 = jdim(ig) - 1
      kbd1 = kdim(ig) - 1
      lbd1 = ldim(ig) - 1
c 
      do i=ibsptr(ig),ibeptr(ig)
c
         jbd1 = jdim(ibblk(i)) - 1
         kbd1 = kdim(ibblk(i)) - 1
         lbd1 = ldim(ibblk(i)) - 1
         ii   = ibc(i)
         jid1 = jdim(iiblk(ii)) - 1
         kid1 = kdim(iiblk(ii)) - 1
         lid1 = ldim(iiblk(ii)) - 1
c
c        check indicies of target point
c
         if (jb(i).lt.1 .or. jb(i).gt.jbd1) then
            write(7,'(''target cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') i
            write(7,'(''     jb,   kb,   lb = '',3i4)') jb(i),kb(i),
     .      lb(i)
            write(7,'(''jdim1, kdim1, ldim1 = '',3i4)') jbd1,kbd1,lbd1
            ierr = ierr + 1
         end if
         if (kb(i).lt.1 .or. kb(i).gt.kbd1) then
            write(7,'(''target cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') i
            write(7,'(''     jb,   kb,   lb = '',3i4)') jb(i),kb(i),
     .      lb(i)
            write(7,'(''jdim1, kdim1, ldim1 = '',3i4)') jbd1,kbd1,lbd1
            ierr = ierr + 1
         end if
         if (lb(i).lt.1 .or. lb(i).gt.lbd1) then
            write(7,'(''target cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') i
            write(7,'(''     jb,   kb,   lb = '',3i4)') jb(i),kb(i),
     .      lb(i)
            write(7,'(''jdim1, kdim1, ldim1 = '',3i4)') jbd1,kbd1,lbd1
            ierr = ierr + 1
         end if
c
c        check indicies of source point
c
         if (ji(ii).lt.0 .or. ji(ii).gt.jid1+2) then
            write(7,'(''source cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') ii
            write(7,'(''     ji,   ki,   li = '',3i4)') ji(ii),
     .      ki(ii),li(ii)
            write(7,'(''jdim1, kdim1, ldim1 =  '',3i4)') jid1,kid1,lid1
            write(7,'(''   dj = '',3f8.3)')  dj(ii)
            ierr = ierr + 1
         end if
         if (ki(ii).lt.0 .or. ki(ii).gt.kid1+2) then
            write(7,'(''source cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') ii
            write(7,'(''     ji,   ki,   li = '',3i4)') ji(ii),
     .      ki(ii),li(ii)
            write(7,'(''jdim1, kdim1, ldim1 =  '',3i4)') jid1,kid1,lid1
            write(7,'(''   dk = '',3f8.3)') dk(ii)
            ierr = ierr + 1
         end if
         if (li(ii).lt.0 .or. li(ii).gt.lid1+2) then
            write(7,'(''source cell number '',i6,'' has one or '',
     .      ''more illegal indicies'')') ii
            write(7,'(''     ji,   ki,   li = '',3i4)') ji(ii),
     .      ki(ii),li(ii)
            write(7,'(''jdim1, kdim1, ldim1 =  '',3i4)') jid1,kid1,lid1 
            write(7,'(''   dl = '',3f8.3)') dl(ii)
            ierr = ierr + 1
         end if
c
      end do
c
c     count up target points, boundary points and print summary
c
      ld = ldim(ig)
      jd = jdim(ig)
      kd = kdim(ig)
      nl1    = 0
      nl2    = 0
      nj1    = 0
      nj2    = 0
      nk1    = 0
      nk2    = 0
      nldim1 = 0
      nldim2 = 0
      njdim1 = 0
      njdim2 = 0
      nkdim1 = 0
      nkdim2 = 0
      do i=ibsptr(ig),ibeptr(ig)
         if (i2d .eq. 0) then
            if (lb(i).eq.1) nl1 = nl1 + 1
            if (lb(i).eq.2) nl2 = nl2 + 1
            if (lb(i).eq.ld-1) nldim1 = nldim1 + 1
            if (lb(i).eq.ld-2) nldim2 = nldim2 + 1
         end if
         if (jb(i).eq.1) nj1 = nj1 + 1
         if (jb(i).eq.2) nj2 = nj2 + 1
         if (jb(i).eq.jd-1) njdim1 = njdim1 + 1
         if (jb(i).eq.jd-2) njdim2 = njdim2 + 1
         if (kb(i).eq.1) nk1 = nk1 + 1
         if (kb(i).eq.2) nk2 = nk2 + 1
         if (kb(i).eq.kd-1) nkdim1 = nkdim1 + 1
         if (kb(i).eq.kd-2) nkdim2 = nkdim2 + 1
      end do
      write(6,*)
      write(6,'(''target point summary for zone'',i5)') ig
      write(6,'(''  total number of target points = '',i6)')
     .ibeptr(ig)-ibsptr(ig)+1
      write(6,'(''  number of points at j=1       = '',i6)')nj1
      write(6,'(''  number of points at j=2       = '',i6)')nj2
      write(6,'(''  number of points at j=jdim-1  = '',i6)')njdim1
      write(6,'(''  number of points at j=jdim-2  = '',i6)')njdim2
      write(6,'(''  number of points at k=1       = '',i6)')nk1
      write(6,'(''  number of points at k=2       = '',i6)')nk2
      write(6,'(''  number of points at k=kdim-1  = '',i6)')nkdim1
      write(6,'(''  number of points at k=kdim-2  = '',i6)')nkdim2
      if (i2d .eq. 0) then
         write(6,'(''  number of points at l=1       = '',i6)')nl1
         write(6,'(''  number of points at l=2       = '',i6)')nl2
         write(6,'(''  number of points at l=ldim-1  = '',i6)')nldim1
         write(6,'(''  number of points at l=ldim-2  = '',i6)')nldim2
      end if
      write(6,'(''  number of orphan points       = '',i6)')norph(ig)
c
      end do
c
      return
      end
c
      subroutine wstencil(iitot,maxbl,jb,kb,lb,ji,ki,li,dj,dk,dl,
     .                    iiblk,ibsptr,ibeptr,ibc,iunit,norpht,
     .                    jborph,kborph,lborph,igorph,ngrid)
c
      dimension jb(iitot),kb(iitot),lb(iitot),ibc(iitot)
      dimension ji(iitot),ki(iitot),li(iitot)
      dimension dj(iitot),dk(iitot),dl(iitot),iiblk(iitot)
      dimension ibsptr(maxbl),ibeptr(maxbl)
      dimension jborph(iitot),kborph(iitot),lborph(iitot),
     .          igorph(iitot)
c
      do ig=1,ngrid
         write(iunit,*)
         write(iunit,'(''stencils for block: '',i4)') ig
         do i=ibsptr(ig),ibeptr(ig)
            write(iunit,*)
            write(iunit,'(''for target: '',3i4)') jb(i),kb(i),lb(i)
            ii = ibc(i)
            write(iunit,'(''  pointer ibc(i):'',i5)') ibc(i)
            write(iunit,'(''  source block:'',3i4)') iiblk(ii)
            write(iunit,'(''  source point:'',3i4)') ji(ii),ki(ii),
     .      li(ii)
            write(iunit,'(''  interp coeff:'',3f10.4)') dj(ii),dk(ii),
     .      dl(ii)
         end do
      end do
c
      write(iunit,*)
      write(iunit,'(''total number of orphan points remaining: '',
     .i4)') norpht
      if (norpht .gt. 0) then
         do n=1,norpht
            write(iunit,'(''  block '',i4,'' orphan pt.: '',3i4)')
     .      igorph(n),jborph(n),kborph(n),lborph(n)
         end do
      end if
c
      return
      end
c
      subroutine wgrid(jmax,kmax,lmax,ld,jd,kd,x,y,z,iblank,
     .                 lunit,i2d)
c
c     outputs iblanked grid at cfl3d cell centers. the input grid
c     is the augmented cell center grid used by pegsus.
c
      dimension x(jmax,kmax,lmax),y(jmax,kmax,lmax),
     .          z(jmax,kmax,lmax)
      real*4    xout(jmax,kmax,lmax),yout(jmax,kmax,lmax),
     .          zout(jmax,kmax,lmax)
      dimension iblank(jmax,kmax,lmax)
      pointer (ip_xout,xout),(ip_yout,yout),(ip_zout,zout)
c
      icall  = 0
      memuse = 0
      call umalloc_r(ip_xout,lmax*jmax*kmax,-1,icall,memuse)
      call umalloc_r(ip_yout,lmax*jmax*kmax,-1,icall,memuse)
      call umalloc_r(ip_zout,lmax*jmax*kmax,-1,icall,memuse)
c
      jd1 = jd + 1
      kd1 = kd + 1
      ld1 = ld + 1
      jd2 = jd + 2
      kd2 = kd + 2
      ld2 = ld + 2
c
      do k=1,kd2
         do j=1,jd2
            do l=1,ld2
               xout(j,k,l) = x(j,k,l)
               yout(j,k,l) = y(j,k,l)
               zout(j,k,l) = z(j,k,l)
            end do
         end do
      end do
c
      if (i2d .eq. 0) then
         write(lunit) (((xout(j,k,l),l=2,ld1),j=2,jd1),k=2,kd1),
     .                (((yout(j,k,l),l=2,ld1),j=2,jd1),k=2,kd1),
     .                (((zout(j,k,l),l=2,ld1),j=2,jd1),k=2,kd1),
     .              (((iblank(j,k,l),l=2,ld1),j=2,jd1),k=2,kd1)
      else
         write(lunit) (((xout(j,k,l),l=2,ld1),j=2,jd1),k=2,kd1),
     .                (((zout(j,k,l),l=2,ld1),j=2,jd1),k=2,kd1),
     .              (((iblank(j,k,l),l=2,ld1),j=2,jd1),k=2,kd1)
      end if
c
      call ufree_r(ip_xout)
      call ufree_r(ip_yout)
      call ufree_r(ip_zout)
c
      return
      end
c
      subroutine worph(jmax,kmax,lmax,lomax,nbmax,xin,yin,zin,
     .                 ngrid,norpht,lborph,jborph,kborph,igorph,norph,
     .                 lunit,nppblk,i2d,ipeg)
c
      dimension xin(jmax,kmax,lmax),yin(jmax,kmax,lmax),
     .          zin(jmax,kmax,lmax)
      dimension lborph(*),jborph(*),kborph(*),igorph(*),norph(*)
      dimension idim(nbmax),jdim(nbmax),kdim(nbmax),ldim(nbmax)
      real*4    xout(lomax),yout(lomax),zout(lomax)
      pointer (ip_xout,xout),(ip_yout,yout),(ip_zout,zout)
      pointer (ip_idim,idim),(ip_jdim,jdim),(ip_kdim,kdim),
     .        (ip_ldim,ldim)
c
      icall  = 0
      memuse = 0
      call umalloc_r(ip_xout,lomax,-1,icall,memuse)
      call umalloc_r(ip_yout,lomax,-1,icall,memuse)
      call umalloc_r(ip_zout,lomax,-1,icall,memuse)
      call umalloc_r(ip_idim,nbmax, 1,icall,memuse)
      call umalloc_r(ip_jdim,nbmax, 1,icall,memuse)
      call umalloc_r(ip_kdim,nbmax, 1,icall,memuse)
      call umalloc_r(ip_ldim,nbmax, 1,icall,memuse)
c
c     write header
c
      iblko = 0
      do ig=1,ngrid
         if(norph(ig).ne.0) then
            nout1 = norph(ig)/nppblk
            noutr = norph(ig) - nout1*nppblk
            if(noutr.eq.0) then
               ioutd = 0
            else
               ioutd = 1
            endif
            l2 = 0
            ipart = 0
            do iout=1,nout1+ioutd
               iblko = iblko + 1
               ipart = ipart + 1
               if(nout1+ioutd.eq.1) then
                  write(6,'(3x,''zone'',i3,
     .                      '' contains orphans from grid'',i3,
     .                      )') iblko,ig
               else
                  write(6,'(3x,''zone'',i3,
     .                      '' contains orphans from grid'',i3,
     .                      '' (part'',i2,'')''
     .                      )') iblko,ig,ipart
               endif
               l1 = l2 + 1
               if(iout.le.nout1) then
                  l2 = iout *nppblk
               else
                  l2 = norph(ig)
               endif
               idim(iblko) = l2 - l1 + 1
            end do
         end if
      end do
c
      ione = 1
      write(lunit) iblko
      write(lunit) (idim(l),ione,ione,l=1,iblko)
c
c     write data
c
      rewind(1)
      if (ipeg .eq. 5) then
         read(1)
         read(1) (jdim(ig),kdim(ig),ldim(ig),ig=1,ngrid)
      end if
      do ig=1,ngrid
         if (ipeg .eq. 4) then
            read(1) nm,jd0,kd0,ld0
         else
            jd0 = jdim(ig)
            kd0 = kdim(ig)
            ld0 = ldim(ig)
         end if
         read(1) (((xin(j,k,l),j=1,jd0),k=1,kd0),l=1,ld0),
     .           (((yin(j,k,l),j=1,jd0),k=1,kd0),l=1,ld0),
     .           (((zin(j,k,l),j=1,jd0),k=1,kd0),l=1,ld0)
c
         if(norph(ig).gt.0) then
            lo = 0
            do nor = 1,norpht
               if (igorph(nor) .eq. ig) then
                  lo = lo + 1
                  j  = jborph(nor) + 1
                  k  = kborph(nor) + 1
                  l  = lborph(nor) + 1
                  xout(lo) = xin(j,k,l)
                  yout(lo) = yin(j,k,l)
                  zout(lo) = zin(j,k,l)
               end if
            end do
c
            nout1 = norph(ig)/nppblk
            noutr = norph(ig) - nout1*nppblk
            if(noutr.eq.0) then
               ioutd = 0
            else
               ioutd = 1
            endif
            l2 = 0
            do iout=1,nout1+ioutd
               l1 = l2 + 1
               if(iout.le.nout1) then
                  l2 = iout *nppblk
               else
                  l2 = norph(ig)
               endif
c
               if (i2d .eq. 0) then
                  write(lunit) (xout(lo),lo=l1,l2),
     .                         (yout(lo),lo=l1,l2),
     .                         (zout(lo),lo=l1,l2)
               else
                  write(lunit) (xout(lo),lo=l1,l2),
     .                         (zout(lo),lo=l1,l2)
               end if
            end do
         end if
      end do
c
      call ufree_r(ip_xout)
      call ufree_r(ip_yout)
      call ufree_r(ip_zout)
      call ufree_r(ip_idim)
      call ufree_r(ip_jdim)
      call ufree_r(ip_kdim)
      call ufree_r(ip_ldim)
c
      return
      end
