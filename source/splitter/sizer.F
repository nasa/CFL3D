      subroutine sizer(mbloc,maxseg,npmax,mxbli,intmax,nsub1,msegt,
     .                 msegn,msplt,mtot,mbloc0,msplt0,nou,bou,nbuf,
     .                 ibufdim)
c***********************************************************************
c     Purpose: determine the parameters required for the block-splitting
c              code
c        note: the following parameters needed by the block splitter are
c              set to very high values, and should be fine for any 
c              conceivable case (these are set near the bottom of this
c              routine in case they ever need adjusting)
c     mxbli....maximum number if 1-1 zonal interfaces
c     intmax...maximum number of patched interfaces
c     msegt....maximum number of entries in the tlns3d mapfile
c              (currently 53 entries are actually used)
c     maxseg...maximum number of block segments per block face
c
c     in addition, the two input parameters for this sizing routine,
c     mbloc0 and msplt0 are set to very large values in the calling
c     routine, and should be quite sufficient for any conceivable case
c***********************************************************************
c
      character*80 cflinp,cflout,roninp,ronout,sdgridin,
     .             sdgridout,dovrlap,dpatch,dresid,gridin,gridout
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension isplits(mbloc0),jsplits(mbloc0),ksplits(mbloc0)
      dimension idimg(mbloc0),jdimg(mbloc0),kdimg(mbloc0)
      dimension ndirin0(msplt0*mbloc0),isblocin(msplt0*mbloc0),
     .          ndirin(msplt0*mbloc0)
c
      pointer (ip_isplits,isplits),(ip_jsplits,jsplits),
     .        (ip_ksplits,ksplits),(ip_idimg,idimg),(ip_jdimg,jdimg),
     .        (ip_kdimg,kdimg),(ip_ndirin0,ndirin0),
     .        (ip_isblocin,isblocin),(ip_ndirin,ndirin)
c
      common /unit5/ iunit5
c
      icall  = 0
      memuse = 0
      call umalloc(ip_isplits,mbloc0,1,icall,memuse)
      call umalloc(ip_jsplits,mbloc0,1,icall,memuse)
      call umalloc(ip_ksplits,mbloc0,1,icall,memuse)
      call umalloc(ip_idimg,mbloc0,1,icall,memuse)
      call umalloc(ip_jdimg,mbloc0,1,icall,memuse)
      call umalloc(ip_kdimg,mbloc0,1,icall,memuse)
      call umalloc(ip_ndirin0,msplt0*mbloc0,1,icall,memuse)
      call umalloc(ip_ndirin,msplt0*mbloc0,1,icall,memuse)
      call umalloc(ip_isblocin,msplt0*mbloc0,1,icall,memuse)
c
      read(iunit5,*)
      read(iunit5,'(a60)') cflinp
      read(iunit5,'(a60)') roninp
      read(iunit5,'(a60)') gridin
      read(iunit5,'(a60)') sdgridin
      read(iunit5,*)
      read(iunit5,*) icflin,ironin,ibin,isdin
      icflin = abs(icflin)
      read(iunit5,*)
      read(iunit5,'(a60)') cflout
      read(iunit5,'(a60)') ronout
      read(iunit5,'(a60)') gridout
      read(iunit5,'(a60)') sdgridout
      read(iunit5,*)
      read(iunit5,*) icflout,ironout,ibout,isdout
      read(iunit5,*)
c
c     read the either the cfl3d input file or the (plot3d) grid 
c     file to get the number of blocks and grid dimensions
c
      if (cflinp.eq.'null') then
         if (gridin .ne. 'null') then
            if(ibin.eq.1) then
c              the following line is needed to use ieee binary file
c              call asnfile (gridin, '-F f77 -N ieee' , IER)
               open(30,file= gridin,form='unformatted')
            else
               open(30,file= gridin,form='formatted')
            end if
            rewind 30
         end if
         if(ibin.eq.1) then
            read(30) nbl
            read(30) (idimg(n),jdimg(n),kdimg(n),n=1,nbl)
         else
            read(30,*) nbl
            read(30,*) (idimg(n),jdimg(n),kdimg(n),n=1,nbl)
         end if
      else
         open(unit=10,file=cflinp,form='formatted',status='old')
         nskip = 14
         do n=1,nskip
            read(10,*)
         end do
c
c        check for keyword input
c
         iunit5sav = iunit5
         iunit5    = 10
         call readkey(ititr,myid,ibufdim,nbuf,bou,nou,99,-1)
         iunit5    = iunit5sav
c
         if (ititr.eq.0) then
            nskip = 8
         else 
            nskip = 7
         end if
c
         do n=1,nskip
            read(10,*)
         end do
         read(10,*) nbl
         nbl = abs(nbl)
         read(10,*)
         do ibloc=1,nbl
            read(10,*)
         end do
         read(10,*)
         do ibloc=1,nbl
            read(10,*) idimg(ibloc),jdimg(ibloc),kdimg(ibloc)
         end do
      end if
c
      read(iunit5,*) nout
c     shortcut to split all blocks the same as the single block input
      if (nout.lt.0) then
         nout = abs(nout)
         do no = 1,nout
            read(iunit5,*) isblocdum
            read(iunit5,*) ndirin0(no)
            read(iunit5,*)
         end do
         nn = 0
         do nb = 1,nbl
            do no = 1,nout
               nn = nn + 1
               isblocin(nn) = nb
               ndirin(nn)   = ndirin0(no)
            end do
         end do
         nout = nn
      else
         do nn = 1,nout
            read(iunit5,*) isblocin(nn)
            read(iunit5,*) ndirin(nn)
            read(iunit5,*) 
         end do
      end if
c
c     count the number of splits in each direction of each block
c     and determine the maximum number of splits per block
c
      msplt = 0
      do nn=1,nbl
         isplits(nn) = 0
         jsplits(nn) = 0
         ksplits(nn) = 0
         do no = 1,nout
            if (isblocin(no) .eq. nn) then
               if (ndirin(no) .eq. 1) then
                  isplits(nn) = isplits(nn) + 1
               else if (ndirin(no) .eq. 2) then
                  jsplits(nn) = jsplits(nn) + 1
               else if (ndirin(no) .eq. 3) then
                  ksplits(nn) = ksplits(nn) + 1
               end if
            end if
         end do
      mspltnn = isplits(nn) + jsplits(nn) + ksplits(nn)
      if (mspltnn.gt.msplt) msplt = mspltnn
      end do
c
c     use a min. of msplt=1 to avoid allocating 0 bytes
c
      msplt = max(1,msplt)
c
c     determine the number of grid points after splitting
c
      npmax = 0
      do nn=1,nbl
         npmax = npmax + jdimg(nn)*kdimg(nn)*idimg(nn)
         if (isplits(nn).gt.0) then
            do ii=1,isplits(nn)
               npmax = npmax + jdimg(nn)*kdimg(nn)
            end do
            idimg(nn) = idimg(nn) + isplits(nn)
         end if
         if (jsplits(nn).gt.0) then
            do jj=1,jsplits(nn)
               npmax = npmax + idimg(nn)*kdimg(nn)
            end do
            jdimg(nn) = jdimg(nn) + jsplits(nn)
         end if
         if (ksplits(nn).gt.0) then
            do kk=1,ksplits(nn)
               npmax = npmax + jdimg(nn)*idimg(nn)
            end do
            kdimg(nn) = kdimg(nn) + ksplits(nn)
         end if
      end do
c
c     determine number of blocks after splitting
c
      mbloc = 0
      do nn=1,nbl
         isp1 = isplits(nn)+1
         jsp1 = jsplits(nn)+1
         ksp1 = ksplits(nn)+1
         mbloc = mbloc + isp1*jsp1*ksp1
      end do
c
c     the following parameters should be quite sufficient for
c     any case
c
      mxbli  = 2500
      intmax = 1500
      msegt  = 75
      maxseg = 200
c
c     set the rest of the parameters - these depend on ones already
c     set, but if all the ones previously set are large enough, these
c     will be too
c
      nsub1 = mbloc
      msegn = maxseg*6
      mtot  = mbloc*msegn
c
      rewind(iunit5)
      rewind(10)
      rewind(30)
c
c     free up memory
c
      call ufree(ip_isplits)
      call ufree(ip_jsplits)
      call ufree(ip_ksplits)
      call ufree(ip_idimg)
      call ufree(ip_jdimg)
      call ufree(ip_kdimg)
      call ufree(ip_ndirin0)
      call ufree(ip_ndirin)
      call ufree(ip_isblocin)
c
      return
      end
