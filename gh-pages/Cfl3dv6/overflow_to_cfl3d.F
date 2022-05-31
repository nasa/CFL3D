      program cfl3d_to_overflow
c
c***********************************************************************
c     Purpose: Creates a cfl3d input file from an overflow input file
c***********************************************************************
c
      character*80 overinp,cfl3dinp,gridin,gridout
      character*1 yesno_ref
c
c***********************************************************************
c     query for file names and open files
c***********************************************************************
c
      write(6,*)
      write(6,'(''This program tries to correct the one-off BCs'')')
      write(6,'(''that occur with OVERFLOW... it may not work!'')')
      write(6,'(''Please check CFL3D input file after-the-fact.'')')
      write(6,'(''It may work better to edit the original OVERFLOW'')')
      write(6,'(''input file (visual inspection may be easier).'')')
      write(6,*)
      write(6,'(''Input 0 to always choose odd index when there is'')')
      write(6,'(''a choice, or 1 to always ask:'')')
      read(5,*) ichoose
      write(6,*)
      write(6,'(''input name of overflow input file to read'')')
      read(5,'(a80)') overinp
      write(6,*)
      write(6,'(''input name of cfl3d input file to create'')')
      read(5,'(a80)') cfl3dinp
      write(6,*)
      write(6,'(''input name of overflow grid file to read '',
     .          ''(plot3d; mg,unf)'')')
      read(5,'(a80)') gridin
      write(6,*)
      write(6,'(''input name of cfl3d grid file to create '',
     .          ''(plot3d; mg,unf)'')')
      read(5,'(a80)') gridout
      write(6,*)
      write(6,'(''do you wish to enter sref,cref,bref,xmc,ymc,zmc '',
     .          ''now? (y/n)'')')
      write(6,'(''(if not, ??? will be used as placeholders in the '',
     .          ''cfl3d input file'')')
      read(5,'(a1)') yesno_ref
      sref = -999.
      if (yesno_ref .eq. 'y' .or. yesno_ref .eq. 'Y') then
         write(6,*)
         write(6,'(''enter sref (reference area)'')')
         read(5,*) sref
         write(6,*)
         write(6,'(''enter cref (ref. length, typically chord)'')')
         read(5,*) cref
         write(6,*)
         write(6,'(''enter bref (ref. length, typically semispan)'')')
         read(5,*) bref
         write(6,*)
         write(6,'(''enter xmc  (x coordinate of moment center)'')')
         read(5,*) xmc
         write(6,*)
         write(6,'(''enter ymc  (y coordinate of moment center)'')')
         read(5,*) ymc
         write(6,*)
         write(6,'(''enter zmc  (z coordinate of moment center)'')')
         read(5,*) zmc
      end if
c
      open(unit=10,file=overinp,form='formatted',status='old')
      open(unit=15,file=cfl3dinp,form='formatted',status='unknown')
      open(unit=20,file=gridin,form='unformatted',status='old')
      open(unit=25,file=gridout,form='unformatted',status='unknown')
      open(unit=30,file='gridnames.dat',form='formatted',
     .     status='unknown')
c
c***********************************************************************
c     determine  array sizes
c***********************************************************************
c
c***********************************************************************
c     convert to cfl3d input file
c***********************************************************************
c
      call convert(gridin,gridout,nobc,sref,cref,bref,xmc,ymc,zmc,
     +             ichoose)
c
      write(6,*)
      write(6,'(''DONE'')')
      if (yesno_ref .eq. 'n' .or. yesno_ref .eq. 'N') then
         write(6,*)
         write(6,'(''Note: you must edit the cfl3d input file before '',
     .             ''running cfl3d'')')
         write(6,'(''      and insert values for sref, cref, etc.'')')
      end if
      if (nobc .gt. 0) then
         write(6,*)
         write(6,'(''WARNING: some overflow bcs could not be '',
     .             ''translated into cfl3d bcs'')')
         write(6,'(''         these are flagged as ??? in the cfl3d '',
     .             ''input file'')')
      end if
      write(6,*)
      write(6,'(''The file gridnames.dat contains the names of the '',
     .          ''grids used in'')')
      write(6,'(''the overflow input file'')')
c
      stop
      end
c
      subroutine convert(gridin,gridout,nobc,sref,cref,bref,xmc,ymc,zmc,
     +                   ichoose)
c
      parameter (maxbl=500,maxseg=400,maxsegt=400)
c
      character*80 grid,plt3dg,plt3dq,output,resid,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrtf,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,inpfile
      character*40 string
      character*40 name
      character*80 gridin,gridout
      character*3  qm
      character*80 bcfile(maxseg),bcfile0(maxseg,maxbl)
      character*80 bcfilei(maxbl,maxseg,2),bcfilej(maxbl,maxseg,2),
     .             bcfilek(maxbl,maxseg,2)
      character*80 tlfile(maxseg),tlfile0(maxsegt,maxbl)
      character*10 datahdr(10)
c
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),
     .          nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          jbcinfo(maxbl,maxseg,7,2),ibcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2),jdimg(maxbl),kdimg(maxbl),
     .          idimg(maxbl),igridg(maxbl),levelg(maxbl),iovrlp(maxbl),
     .          bcvali(maxbl,maxseg,7,2),bcvalj(maxbl,maxseg,7,2),
     .          bcvalk(maxbl,maxseg,7,2)
      dimension nblk(2,maxseg),limblk(2,6,maxseg),
     .          isva(2,2,maxseg),nblon(maxseg)
      dimension fmgcyc(5),gvec(3)
      dimension iter0(maxbl)
      dimension irhs0(maxbl),ilhs0(maxbl),ilhsit0(maxbl),idiss0(maxbl),
     .          bimin0(maxbl),smoop0(maxbl),smooc0(maxbl),
     .          smoor0(maxbl)
      dimension itime0(maxbl),irelax0(maxbl),dt0(maxbl),tfoso0(maxbl),
     .          cflmin0(maxbl),cflmax0(maxbl)
      dimension ispec0(maxbl),ispecj0(maxbl),ispeck0(maxbl),
     .          ispecl0(maxbl),smoo0(maxbl),dis20(maxbl),dis2j0(maxbl),
     .          dis2k0(maxbl),dis2l0(maxbl),dis40(maxbl),dis4j0(maxbl),
     .          dis4k0(maxbl),dis4l0(maxbl),epse0(maxbl),epssgs0(maxbl),
     .          delta0(maxbl),fso0(maxbl),vepsl0(maxbl),vepsn0(maxbl)
      dimension ittyp(maxsegt),itdir(maxsegt),
     .          jtls(maxsegt),jtle(maxsegt),
     .          ktls(maxsegt),ktle(maxsegt),
     .          ltls(maxsegt),ltle(maxsegt),
     .          tlpar1(maxsegt)
      dimension ibtyp(maxseg),
     .          ibdir(maxseg),jbcs(maxseg),
     .          jbce(maxseg),kbcs(maxseg),
     .          kbce(maxseg),lbcs(maxseg),
     .          lbce(maxseg),bcpar1(maxseg)
      dimension nturb0(maxbl),
     .          ittyp0(maxsegt,maxbl),itdir0(maxsegt,maxbl),
     .          jtls0(maxsegt,maxbl),jtle0(maxsegt,maxbl),
     .          ktls0(maxsegt,maxbl),ktle0(maxsegt,maxbl),
     .          ltls0(maxsegt,maxbl),ltle0(maxsegt,maxbl),
     .          tlpar10(maxsegt,maxbl),
     .          itert0(maxbl),itlhit0(maxbl),cflt0(maxbl)
      dimension nbc0(maxbl),ibtyp0(maxseg,maxbl),
     .          ibdir0(maxseg,maxbl),jbcs0(maxseg,maxbl),
     .          jbce0(maxseg,maxbl),kbcs0(maxseg,maxbl),
     .          kbce0(maxseg,maxbl),lbcs0(maxseg,maxbl),
     .          lbce0(maxseg,maxbl),bcpar10(maxseg,maxbl)
      dimension iterc0(maxbl),cflc0(maxbl),iupc0(maxbl),dis2c0(maxbl),
     .          dis4c0(maxbl)
      dimension scinf(2),smw(2),alt0(2),alt1(2),
     .          alt2(2),alt3(2),alt4(2),aut0(2),aut1(2),aut2(2),
     .          aut3(2),aut4(2)
      dimension ivisci(maxbl),iviscj(maxbl),ivisck(maxbl),iforce(maxbl)
      dimension ilamlog(maxbl),ilamhig(maxbl),jlamlog(maxbl),
     .          jlamhig(maxbl),klamlog(maxbl),klamhig(maxbl)
      dimension epsssc(3),epsssr(3)
c
      logical chimra,cdisc,incore,restrt,fmg
      logical targcl,multig,multig1,multig10(maxbl)
      logical corsvi,corsvi0(maxbl)
      logical roeavg,roeavg0(maxbl)
      logical visc,viscj,visck,viscl,viscx
      logical visc0(maxbl),viscj0(maxbl),visck0(maxbl),viscl0(maxbl),
     .        viscx0(maxbl)
c
      namelist /global/ chimra,cdisc,incore,nsteps,restrt,nsave,nfomo,
     .                  nqt,nqc,multig,fmg,nglvl,fmgcyc,dtphys,nitnwt,
     .                  ordnwt,fsonwt,walldist
      namelist /floinp/ alpha,beta,fsmach,gaminf,rey,pr,prt,tinf,retinf,
     .                  xkinf,froude,gvec,targcl,clalph,ntarg,cltarg 
      namelist /vargam/ igam,ht1,ht2,scinf,smw,alt0,alt1,alt2,alt3,alt4,
     .                  aut0,aut1,aut2,aut3,aut4
      namelist /grdnam/ name
      namelist /niters/ iter
      namelist /metprm/ irhs,ilhs,ilhsit,idiss,bimin,multig1,smoop,
     .                  smooc,smoor,corsvi
      namelist /timacu/ itime,irelax,dt,tfoso,cflmin,cflmax
      namelist /smoacu/ ispec,ispecj,ispeck,ispecl,smoo,dis2,dis2j,
     .                  dis2k,dis2l,dis4,dis4j,dis4k,dis4l,epse,epssgs,
     .                  delta,fso,vepsl,vepsn,roeavg
      namelist /visinp/ visc,viscj,visck,viscl,viscx,nturb,ittyp,itdir,
     .                  jtls,jtle,ktls,ktle,ltls,ltle,tlpar1,tlfile,
     .                  itert,itlhit,cflt,irc
      namelist /bcinp/  nbc,ibtyp,ibdir,jbcs,jbce,kbcs,kbce,lbcs,lbce,
     .                  bcpar1,bcfile
      namelist /sceinp/ iterc,cflc,iupc,dis2c,dis4c
c
      qm = '???'
c
c***********************************************************************
c     read grid file to get number of grids and block dimensions
c***********************************************************************
c
      read(20) ngrid
      if (ngrid .gt. maxbl) then
         write(6,*) 'stopping...'
         write(6,*) 'increase parameter maxbl to ',ngrid
      end if
      read(20) (jdimg(n),kdimg(n),idimg(n),n=1,ngrid)
c
c***********************************************************************
c     output plot3d grid for cfl3d in "i,j,k" form, rather than "j,k,l"
c     form used by overflow. this allows a 1 to 1 correspondence 
c     in the cfl3d and overflow INPUT files between j, k, and i(=l)
c***********************************************************************
c
c---rtb
      i2d     = 1
      i2dmax = 3
      do n=1,ngrid
         if (idimg(n) .ne. i2dmax) i2d = 0
      end do
      if (i2d .ne. 0) then
         do n=1,ngrid
            idimg(n) = 2
         end do
      end if
c
c---rtb
      write(25) ngrid
      write(25) (idimg(n),jdimg(n),kdimg(n),n=1,ngrid)
      do n=1,ngrid
         jd = jdimg(n)
         kd = kdimg(n)
         id = idimg(n)
c---rtb
         if (i2d .ne. 0) id = idimg(n)+1
c---rtb
         call rwgrid(jd,kd,id,n,i2d)
      end do
c
c***********************************************************************
c     initialize overflow namelist data, then read data from input file
c***********************************************************************
c
      do n=1,ngrid
c
         if (n .eq. 1) then
c
c           $global
c
            chimra    = .false.
            cdisc     = .false.
            incore    = .false.
            nsteps    = 0
            restrt    = .false.
            nsave     = 100
            nfomo     = 10
            nqt       = 0
            nqc       = 0
            multig    = .false.
            fmg       = .false.
            nglvl     = 3
            fmgcyc(1) = 300
            fmgcyc(2) = 300
            fmgcyc(3) = 300
            fmgcyc(4) = 300
            fmgcyc(5) = 300
            dtphys    = 0.
            nitnwt    = 0
            ordnwt    = 0.
            fsonwt    = 2.
            read(10,global)
c
c           $floinp
c
            alpha   = 0.
            beta    = 0.
            fsmach  = 0.
            gaminf  = 1.4
            rey     = 0.
            pr      = .72
            prt     = .9
            tinf    = 518.7
            retinf  = 0.1
            xkinf   = 0.0001
            froude  = 0.
            gvec(1) = 0
            gvec(2) = 0
            gvec(3) = 0
            targcl  = .false.
            clalph  = 0.1
            ntarg   = 200
            cltarg  = 0.5
            read(10,floinp)
c
c           $vargam
c
            igam     = 0
            ht1      = 10.
            ht2      = 10.
            scinf(1) = 1.
            scinf(2) = 0.
            smw(1)   = 1.
            smw(2)   = 1.
            do nn=1,2
               alt0(nn)  = gaminf/(gaminf-1.)
               alt1(nn)  = 0.
               alt2(nn)  = 0.
               alt3(nn)  = 0.
               alt4(nn)  = 0.
               aut0(nn)  = gaminf/(gaminf-1.)
               aut1(nn)  = 0.
               aut2(nn)  = 0.
               aut3(nn)  = 0.
               aut4(nn)  = 0.
            end do
            read(10,vargam)
c
         end if
c
c        $grdnam
c
         name = ' '
         read(10,*)
         read(10,grdnam)
c
c        output grid names used in overflow input file
c
         write(30,'(a40)') name
c
c
c        $niters
c
         iter = 1
         read(10,niters)
         iter0(n) = iter 
c
c        $metprm
c
         irhs    = 0
         ilhs    = 2
         ilhsit  = 1
         idiss   = 3
         bimin   = 1.
         multig1 = multig
         smoop   = 0.
         smooc   = 0.
         smoor   = 0.
         corsvi  = .true.
         read(10,metprm)
         irhs0(n)    = irhs
         ilhs0(n)    = ilhs
         ilhsit0(n)  = ilhsit
         idiss0(n)   = idiss
         bimin0(n)   = bimin
         multig10(n) = multig1
         smoop0(n)   = smoop
         smooc0(n)   = smooc
         smoor0(n)   = smoor
         corsvi0(n)  = corsvi
c
c        $timacu
c
         itime  = 1
         irelax = 0
         dt     = 0.5
         tfoso  = 1.
         cflmin = 0.
         cflmax = 0.
         read(10,timacu)
         itime0(n)  = itime
         irelax0(n) = irelax
         dt0(n)     = dt
         tfoso0(n)  = tfoso
         cflmin0(n) = cflmin
         cflmax0(n) = cflmax
c
c        $smoacu
c
         ispec  = 2
         ispecj = ispec
         ispeck = ispec
         ispecl = ispec
         smoo   = 1.
         dis2   = 2.
         dis2j  = dis2
         dis2k  = dis2
         dis2l  = dis2
         dis4   = 0.04
         dis4j  = dis4
         dis4k  = dis4
         dis4l  = dis4
         epse   = 0.05
         epssgs = 0.02
         delta  = 1.
         fso    = 2.0
         if (irhs .eq. 4) fso = 3.0
         vepsl  = 0.
         vepsn  = 0.
         roeavg = .false.
         read(10,smoacu)
         ispec0(n)  = ispec
         ispecj0(n) = ispecj
         ispeck0(n) = ispeck
         ispecl0(n) = ispecl
         smoo0(n)   = smoo
         dis20(n)   = dis2
         dis2j0(n)  = dis2j
         dis2k0(n)  = dis2k
         dis2l0(n)  = dis2l
         dis40(n)   = dis4
         dis4j0(n)  = dis4j
         dis4k0(n)  = dis4k
         dis4l0(n)  = dis4l
         epse0(n)   = epse
         epssgs0(n) = epssgs
         delta0(n)  = delta
         fso0(n)    = fso
         vepsl0(n)  = vepsl
         vepsn0(n)  = vepsn
         roeavg0(n) = roeavg
c
c        $visinp
c
         visc   = .false.
         viscj  = .false.
         visck  = .false.
         viscl  = .false.
         viscx  = .false.
         nturb  = 0
         do nn=1,maxseg
            ittyp(nn)  = 0
            itdir(nn)  = 0
            jtls(nn)   = 1
            jtle(nn)   = 1
            ktls(nn)   = 1
            ktle(nn)   = 1
            ltls(nn)   = 1
            ltle(nn)   = 1
            tlpar1(nn) = 0.
            tlfile(nn) = ' '
         end do
         itert  = 1
         itlhit = 1
         cflt   = 1.
         read(10,visinp)
c        count nturb since it is an optional input
         nturb = 0
         do nn=1,maxsegt
            if (ittyp(nn) .ne. 0) nturb = nturb + 1
         end do
         if (nturb .gt. maxsegt) then
            write(6,*) 'stopping...'
            write(6,*) 'increase parameter maxsegt to ',nturb
            stop
         end if
         visc0(n)   = visc
         viscj0(n)  = viscj
         visck0(n)  = visck
         viscl0(n)  = viscl
         viscx0(n)  = viscx
         nturb0(n)  = nturb
         do nn=1,nturb
            ittyp0(nn,n)  = ittyp(nn)
            itdir0(nn,n)  = itdir(nn)
            jtls0(nn,n)   = jtls(nn)
            jtle0(nn,n)   = jtle(nn)
            ktls0(nn,n)   = ktls(nn)
            ktle0(nn,n)   = ktle(nn)
            ltls0(nn,n)   = ltls(nn)
            ltle0(nn,n)   = ltle(nn)
            tlpar10(nn,n) = tlpar1(nn)
            tlfile0(nn,n) = tlfile(nn)
         end do
         itert0(n)  = itert
         itlhit0(n) = itlhit
         cflt0(n)   = cflt
c
c        $bcinp
c
         nbc = 0 
         do nn=1,maxsegt
            ibtyp(nn)  = 0
            ibdir(nn)  = 1
            jbcs(nn)   = 1
            jbce(nn)   = 1
            kbcs(nn)   = 1
            kbce(nn)   = 1
            lbcs(nn)   = 1
            lbce(nn)   = 1
            bcpar1(nn) = 0
            bcfile(nn) = ' '
         end do
         read(10,bcinp)
c        count nbc since it is an optional input
         nbc = 0
         do nn=1,maxseg
            if (ibtyp(nn) .ne. 0) nbc = nbc + 1
         end do
         if (nbc .gt. maxsegt) then
            write(6,*) 'stopping...'
            write(6,*) 'increase parameter maxsegt to ',nbc
            stop
         end if
         nbc0(n) = nbc
         do nn=1,maxsegt
            ibtyp0(nn,n)  = ibtyp(nn)
            ibdir0(nn,n)  = ibdir(nn)
            jbcs0(nn,n)   = jbcs(nn)
            jbce0(nn,n)   = jbce(nn)
            kbcs0(nn,n)   = kbcs(nn)
            kbce0(nn,n)   = kbce(nn)
            lbcs0(nn,n)   = lbcs(nn)
            lbce0(nn,n)   = lbce(nn)
            bcpar10(nn,n) = bcpar1(nn)
            bcfile0(nn,n) = bcfile(nn)
         end do
c
c        $sceinp
c
         iterc = 1
         cflc  = 1.
         iupc  = 1
         dis2c = 2.
         dis4c = 0.04
         read(10,sceinp)
         iterc0(n) = iterc
         cflc0(n)  = cflc
         iupc0(n)  = iupc
         dis2c0(n) = dis2c
         dis4c0(n) = dis4c
c
      end do 
c
      do n=1,ngrid
c
c        mark bc's that only span lines or points - these
c        have no meaning in a finite volume code, and will be 
c        eliminated later
c
         do nn=1,nbc0(n)
            lrng = 1
            jrng = 1
            krng = 1
            if (lbce0(nn,n) .eq. lbcs0(nn,n)) lrng = 0
            if (jbce0(nn,n) .eq. jbcs0(nn,n)) jrng = 0
            if (kbce0(nn,n) .eq. kbcs0(nn,n)) krng = 0
            if (lrng+jrng+krng .lt. 2) ibtyp0(nn,n) = -999
         end do
c
c        add matching periodic face for periodic bcs
c        (overflow set only one side; cfl3d needs both)
c
         nnn = nbc0(n)
         do nn=1,nbc0(n)
            if (ibtyp0(nn,n) .eq. 10) then
               nnn = nnn + 1
               ibtyp0(nnn,n)  =  ibtyp0(nn,n)
               ibdir0(nnn,n)  = -ibdir0(nn,n)
               jbcs0(nnn,n)   =  jbcs0(nn,n)
               jbce0(nnn,n)   =  jbce0(nn,n)
               kbcs0(nnn,n)   =  kbcs0(nn,n)
               kbce0(nnn,n)   =  kbce0(nn,n)
               lbcs0(nnn,n)   =  lbcs0(nn,n)
               lbce0(nnn,n)   =  lbce0(nn,n)
               bcpar10(nnn,n) =  bcpar10(nn,n)
               bcfile0(nnn,n) =  bcfile0(nn,n)
            end if
         end do
         nbc0(n) = max(nbc0(n),nnn)
c
c        convert overflow "shortcut" input to actual index ranges
c
         do nn=1,nbc0(n)
            if (jbcs0(nn,n) .lt. 0) then
               jdecr = jbcs0(nn,n) + 1
               jbcs0(nn,n) = jdimg(n) + jdecr
            end if
            if (jbce0(nn,n) .lt. 0) then
               jdecr = jbce0(nn,n) + 1
               jbce0(nn,n) = jdimg(n) + jdecr
            end if
            if (kbcs0(nn,n) .lt. 0) then
               kdecr = kbcs0(nn,n) + 1
               kbcs0(nn,n) = kdimg(n) + kdecr
            end if
            if (kbce0(nn,n) .lt. 0) then
               kdecr = kbce0(nn,n) + 1
               kbce0(nn,n) = kdimg(n) + kdecr
            end if
            if (lbcs0(nn,n) .lt. 0) then
               ldecr = lbcs0(nn,n) + 1
               lbcs0(nn,n) = idimg(n) + ldecr
            end if
            if (lbce0(nn,n) .lt. 0) then
               ldecr = lbce0(nn,n) + 1
               lbce0(nn,n) = idimg(n) + ldecr
            end if
         end do
c
c        try to find BCs that are off by 1
c
         if (nbc0(n) .gt. 1) then
         do nn=1,nbc0(n)
           do mm=1,nbc0(n)
             if (abs(jbce0(nn,n)-jbcs0(mm,n)) .eq. 1) then
               if (ichoose .eq. 1) then
                 write(6,'('' off-by-one detected at n='',i5)') n
                 write(6,'('' jbce0(nn,n)='',i5)') jbce0(nn,n)
                 write(6,'('' jbcs0(mm,n)='',i5)') jbcs0(mm,n)
                 write(6,'('' input 1 to select 1st, 2 for 2nd:'')')
                 read(5,*) iselect
                 if (iselect .eq. 1) then
                   jbcs0(mm,n)=jbce0(nn,n)
                 else if (iselect .eq. 2) then
                   jbce0(nn,n)=jbcs0(mm,n)
                 else
                   write(6,'('' no change made'')')
                 end if
               else
                 if (mod(jbce0(nn,n),2) .eq. 0) then
                   jbce0(nn,n)=jbcs0(mm,n)
                 else
                   jbcs0(mm,n)=jbce0(nn,n)
                 end if
               end if
             end if
c
             if (abs(kbce0(nn,n)-kbcs0(mm,n)) .eq. 1) then
               if (ichoose .eq. 1) then
                 write(6,'('' off-by-one detected at n='',i5)') n
                 write(6,'('' kbce0(nn,n)='',i5)') kbce0(nn,n)
                 write(6,'('' kbcs0(mm,n)='',i5)') kbcs0(mm,n)
                 write(6,'('' input 1 to select 1st, 2 for 2nd:'')')
                 read(5,*) iselect
                 if (iselect .eq. 1) then
                   kbcs0(mm,n)=kbce0(nn,n)
                 else if (iselect .eq. 2) then
                   kbce0(nn,n)=kbcs0(mm,n)
                 else
                   write(6,'('' no change made'')')
                 end if
               else
                 if (mod(kbce0(nn,n),2) .eq. 0) then
                   kbce0(nn,n)=kbcs0(mm,n)
                 else
                   kbcs0(mm,n)=kbce0(nn,n)
                 end if
               end if
             end if
c
             if (abs(lbce0(nn,n)-lbcs0(mm,n)) .eq. 1) then
               if (ichoose .eq. 1) then
                 write(6,'('' off-by-one detected at n='',i5)') n
                 write(6,'('' lbce0(nn,n)='',i5)') lbce0(nn,n)
                 write(6,'('' lbcs0(mm,n)='',i5)') lbcs0(mm,n)
                 write(6,'('' input 1 to select 1st, 2 for 2nd:'')')
                 read(5,*) iselect
                 if (iselect .eq. 1) then
                   lbcs0(mm,n)=lbce0(nn,n)
                 else if (iselect .eq. 2) then
                   lbce0(nn,n)=lbcs0(mm,n)
                 else
                   write(6,'('' no change made'')')
                 end if
               else
                 if (mod(lbce0(nn,n),2) .eq. 0) then
                   lbce0(nn,n)=lbcs0(mm,n)
                 else
                   lbcs0(mm,n)=lbce0(nn,n)
                 end if
               end if
             end if
           enddo
         enddo
         end if
c
c        mark turb ranges that only span lines or points - these
c        have no meaning in a finite volume code, and will be 
c        eliminated later
c
         do nn=1,nturb0(n)
            lrng = 1
            jrng = 1
            krng = 1
            if (ltle0(nn,n) .eq. ltls0(nn,n)) lrng = 0
            if (jtle0(nn,n) .eq. jtls0(nn,n)) jrng = 0
            if (ktle0(nn,n) .eq. ktls0(nn,n)) krng = 0
            if (lrng+jrng+krng .lt. 2) ittyp0(nn,n) = -999
         end do
c
c        convert overflow "shortcut" input to actual index ranges
c
         do nn=1,nturb0(n)
            if (jtls0(nn,n) .lt. 0) then
               jdecr = jtls0(nn,n) + 1
               jtls0(nn,n) = jdimg(n) + jdecr
            end if
            if (jtle0(nn,n) .lt. 0) then
               jdecr = jtle0(nn,n) + 1
               jtle0(nn,n) = jdimg(n) + jdecr
            end if
            if (ktls0(nn,n) .lt. 0) then
               kdecr = ktls0(nn,n) + 1
               ktls0(nn,n) = kdimg(n) + kdecr
            end if
            if (ktle0(nn,n) .lt. 0) then
               kdecr = ktle0(nn,n) + 1
               ktle0(nn,n) = kdimg(n) + kdecr
            end if
            if (ltls0(nn,n) .lt. 0) then
               ldecr = ltls0(nn,n) + 1
               ltls0(nn,n) = idimg(n) + ldecr
            end if
            if (ltle0(nn,n) .lt. 0) then
               ldecr = ltle0(nn,n) + 1
               ltle0(nn,n) = idimg(n) + ldecr
            end if
         end do
c
      end do 
c
c***********************************************************************
c     write out cfl3d input file
c***********************************************************************
c
c     set generic cfl3d file names
c
      plt3dg   = 'plot3dg.bin'
      plt3dq   = 'plot3dq.bin'
      output   = 'cfl3d.out'
      resid    = 'cfl3d.res'
      turbres  = 'cfl3d.turres'
      blomx    = 'cfl3d.blomax'
      output2  = 'cfl3d.out15'
      printout = 'cfl3d.prout'
      restrtf  = 'restart.bin'
      pplunge  = 'cfl3d.out20'
      ovrlap   = 'ovrlp.bin'
      patch    = 'patch.bin'
      grid     = gridout
c
      write(15,'(''input/output files:'')')
      write(15,'(a60)') grid
      write(15,'(a60)') plt3dg
      write(15,'(a60)') plt3dq
      write(15,'(a60)') output
      write(15,'(a60)') resid
      write(15,'(a60)') turbres
      write(15,'(a60)') blomx
      write(15,'(a60)') output2
      write(15,'(a60)') printout
      write(15,'(a60)') pplunge
      write(15,'(a60)') ovrlap
      write(15,'(a60)') patch
      write(15,'(a60)') restrtf
c
      write(15,'(''Title'')')
c
c     set flow parameters
c
      ialph  = 0
      ihstry = 0
      write(15,22)
      write(15,23) fsmach,alpha,beta,rey/1.e6,tinf,ialph,ihstry
   22 format(6x,4hMach,5x,5halpha,6x,4hbeta,6x,4hReUe,3x,7hTinf,dR,
     .5x,5hialph,4x,6hihstry)
   23 format(3f10.5,e10.3,f10.5,i10,i10,i10)
c
c     set reference dimensions (not available ===> output a ?)
c
      write(15,25)
      if (sref .gt. 0) then
         write(15,26)sref,cref,bref,xmc,ymc,zmc
      else
         write(15,'(6a10)') qm,qm,qm,qm,qm,qm
      end if
   25 format(6x,4hsref,6x,4hcref,6x,4hbref,7x,3hxmc,7x,3hymc,
     .7x,3hzmc)
   26 format(f10.3,7f10.4)
c
c     set dt/cfl number and related quantities
c
      irest = 0
      if (restrt .eq. .true.) irest = 1
      if (dtphys .ne. 0) dt1 = dtphys
      dt1 = -5.0
      iflagts = 0
      fmax    = 1.0
      iunst   = 0
      cfltau  = 5.0
      write(15,33)
      write(15,34) dt1,irest,iflagts,fmax,iunst,cfltau
   33 format(8x,2hdt,5x,5hirest,3x,7hiflagts,6x,4hfmax,5x,5hiunst,
     .3x,7hcfl_tau)
   34 format(f10.5,2i10,f10.5,i10,3f10.5)
c
c     set ngrid, etc; ngrid<0 for plot3d-type grids used by overflow
c
      nplot3d = ngrid
      nprint  = 0
      nwrest  = 9999
      ichk    = 0
c---rtb
c     i2d     = 1
c     i2dmax = 3
c---rtb
      ntstep = 1
      if (dtphys .ne. 0.) ntstep = nsteps
      ita = -2
c---rtb
c     do n=1,ngrid
c        if (idimg(n) .ne. i2dmax) i2d = 0
c     end do
c---rtb
      write(15,35)
      write(15,36) -ngrid,nplot3d,nprint,nwrest,ichk,i2d,ntstep,ita
   35 format(5x,5hngrid,3x,7hnplot3d,4x,6hnprint,4x,6hnwrest,
     .       6x,4hichk,7x,3hi2d,3x,7h ntstep,4x,6h   ita)
   36 format(10i10)
c
c     set ncg, etc; iforce set based on solid wall bc's
c
      ncg      = nglvl-1
      iem      = 0
      iadvance = 0
      ivmxi    = 0
      ivmxj    = 0
      ivmxk    = 0
      itrbmdl  = 0
      if (nqt .eq. 100) itrbmdl  = 4   !Baldwin-Barth
      if (nqt .eq. 101) itrbmdl  = 5   !Spalart-Almaras
      if (nqt .eq. 102) itrbmdl  = 5   !Spalart-Almaras
      if (nqt .eq. 202) itrbmdl  = 8   !k-w
      if (nqt .eq. 202) itrbmdl  = 7   !SST
      do n=1,ngrid
         ifmin = 0
         ifmax = 0
         jfmin = 0
         jfmax = 0
         kfmin = 0
         kfmax = 0
         ivisci(n) = 0
         iviscj(n) = 0
         ivisck(n) = 0
         do nn=1,nbc0(n)
            if (ibtyp0(nn,n).ge.1 .and. ibtyp0(nn,n).le.8) then
               if (abs(ibdir0(nn,n)) .eq. 3) then
                  if (ibdir0(nn,n) .eq.  3) ifmin = 1
                  if (ibdir0(nn,n) .eq. -3) ifmax = 2
               end if
               if (abs(ibdir0(nn,n)) .eq. 2) then
                  if (ibdir0(nn,n) .eq.  2) kfmin = 1
                  if (ibdir0(nn,n) .eq. -2) kfmax = 2
               end if
               if (abs(ibdir0(nn,n)) .eq. 1) then
                  if (ibdir0(nn,n) .eq.  1) jfmin = 1
                  if (ibdir0(nn,n) .eq. -1) jfmax = 2
               end if
               iforce(n) = (ifmin+ifmax)*100 + (jfmin+jfmax)*10
     .                   + (kfmin+kfmax)
            end if
            if (ibtyp0(nn,n).ge.5 .and. ibtyp0(nn,n).le.8) then
               if (abs(ibdir0(nn,n)) .eq. 3) ivmxi = itrbmdl
               if (abs(ibdir0(nn,n)) .eq. 2) ivmxk = itrbmdl
               if (abs(ibdir0(nn,n)) .eq. 1) ivmxj = itrbmdl
            end if
         end do
      end do
      do n=1,ngrid
           ivisci(n) = ivmxi
           iviscj(n) = ivmxj
           ivisck(n) = ivmxk
      end do
      write(15,61)
      do n=1,ngrid
         write(15,36) ncg,iem,iadvance,iforce(n),
     .                ivisci(n),iviscj(n),ivisck(n)
      end do
   61 format(7x,3hncg,7x,3hiem,2x,8hiadvance,4x,6hiforce,
     .       2x,8hivisc(i),2x,8hivisc(j),2x,8hivisc(k))
c
c     block dimensions
c
      write(15,73)
      do n=1,ngrid
         idimg2d = idimg(n)
         if (i2d .gt. 0) i2dd = 2
         write(15,78) idimg2d,jdimg(n),kdimg(n)
      end do
   73 format(6x,4hidim,6x,4hjdim,6x,4hkdim)
   78 format(3i10)
c
c     ilamlo, ilamhi, etc
c
c---rtb
c     going to punt on this for now...
c
      do n=1,ngrid
         ilamlog(n) = 0
         ilamhig(n) = 0
         jlamlog(n) = 0
         jlamhig(n) = 0
         klamlog(n) = 0
         klamhig(n) = 0
      end do
      write(15,80)
      do n=1,ngrid
         write(15,36) ilamlog(n),ilamhig(n),jlamlog(n),
     .                jlamhig(n),klamlog(n),klamhig(n)
      end do
   80 format(4x,6hilamlo,4x,6hilamhi,4x,6hjlamlo,
     .       4x,6hjlamhi,4x,6hklamlo,4x,6hklamhi)
c---rtb
c
c     embeded grid stuff...no analogy in overflow
c
      inewg  = 0
      igridc = 0
      is     = 0
      ie     = 0
      js     = 0
      je     = 0
      ks     = 0
      ke     = 0
      write(15,82)
      do n=1,ngrid
         write(15,36) inewg,igridc,is,js,ks,ie,je,ke
      end do
   82 format(5x,5hinewg,4x,6higridc,
     .       8x,2his,8x,2hjs,8x,2hks,8x,2hie,8x,2hje,8x,2hke)
c
c     idiag, etc...use diagonal scheme and iflim=4
c
      idiagi = 1
      idiagj = 1
      idiagk = 1
      iflimi = 4
      iflimj = 4
      iflimk = 4 
      write(15,84)
      do n=1,ngrid
         write(15,36) idiagi,idiagj,idiagk,iflimi,iflimj,iflimk
      end do
   84 format(2x,8hidiag(i),2x,8hidiag(j),2x,8hidiag(k),
     .       2x,8hiflim(i),2x,8hiflim(j),2x,8hiflim(k))
c
c     flux scheme...use Roe with kappa=1/3 reconstruction
c
      ifdsi = 1
      ifdsj = 1
      ifdsk = 1 
      rkapi = 0.3333
      rkapj = 0.3333
      rkapk = 0.3333
      write(15,86)
      do n=1,ngrid
         write(15,88) ifdsi,ifdsj,ifdsk,rkapi,rkapj,rkapk
      end do
   86 format(3x,7hifds(i),3x,7hifds(j),3x,7hifds(k),
     .       2x,8hrkap0(i),2x,8hrkap0(j),2x,8hrkap0(k))
   88 format(3i10,3f10.4)
c
c     process bc data
c
      do n=1,ngrid
         nseg = 0
         do nn=1,maxseg
            do m=1,7
               ibcinfo(n,nn,m,1) = 0
            end do
         end do
         do nn=1,nbc0(n)
            if (ibtyp0(nn,n) .gt. -999) then
            if (ibdir0(nn,n).eq.3) then
            nseg = nseg + 1
            call bctyper(maxbl,maxseg,n,nn,ibtyp0,lbcs0,lbce0,jbcs0,
     .                   jbce0,kbcs0,kbce0,ibdir0,ibcinfo,1,nseg,
     .                   idimg,jdimg,kdimg,bcvali,bcvalj,bcvalk,
     .                   bcfile0,bcpar10,bcfilei,bcfilej,bcfilek)
            end if
            end if
         end do
c        if no segments have been set, it must be overset
         if (nseg .eq. 0) then
            nseg = nseg + 1
            ibcinfo(n,nseg,1,1) = 0
            ibcinfo(n,nseg,2,1) = 1
            ibcinfo(n,nseg,3,1) = jdimg(n)
            ibcinfo(n,nseg,4,1) = 1
            ibcinfo(n,nseg,5,1) = kdimg(n)
         end if
         nbci0(n) = nseg
      end do
      do n=1,ngrid
         nseg = 0
         do nn=1,maxseg
            do m=1,7
               ibcinfo(n,nn,m,2) = 0
            end do
         end do
         do nn=1,nbc0(n)
            if (ibtyp0(nn,n) .gt. -999) then
            if (ibdir0(nn,n).eq.-3) then
            nseg = nseg + 1
            call bctyper(maxbl,maxseg,n,nn,ibtyp0,lbcs0,lbce0,jbcs0,
     .                   jbce0,kbcs0,kbce0,ibdir0,ibcinfo,2,nseg,
     .                   idimg,jdimg,kdimg,bcvali,bcvalj,bcvalk,
     .                   bcfile0,bcpar10,bcfilei,bcfilej,bcfilek)

            end if
            end if
         end do
c        if no segments have been set, it must be overset
         if (nseg .eq. 0) then
            nseg = nseg + 1
            ibcinfo(n,nseg,1,2) = 0
            ibcinfo(n,nseg,2,2) = 1
            ibcinfo(n,nseg,3,2) = jdimg(n)
            ibcinfo(n,nseg,4,2) = 1
            ibcinfo(n,nseg,5,2) = kdimg(n)
         end if
         nbcidim(n) = nseg
      end do
      do n=1,ngrid
         nseg = 0
         do nn=1,maxseg
            do m=1,7
               jbcinfo(n,nn,m,1) = 0
            end do
         end do
         do nn=1,nbc0(n)
            if (ibtyp0(nn,n) .gt. -999) then
            if (ibdir0(nn,n).eq.1) then
            nseg = nseg + 1
            call bctyper(maxbl,maxseg,n,nn,ibtyp0,lbcs0,lbce0,jbcs0,
     .                   jbce0,kbcs0,kbce0,ibdir0,jbcinfo,1,nseg,
     .                   idimg,jdimg,kdimg,bcvali,bcvalj,bcvalk,
     .                   bcfile0,bcpar10,bcfilei,bcfilej,bcfilek)
            end if
            end if
         end do
c        if no segments have been set, it must be overset
         if (nseg .eq. 0) then
            nseg = nseg + 1
            jbcinfo(n,nseg,1,1) = 0
            jbcinfo(n,nseg,2,1) = 1
            jbcinfo(n,nseg,3,1) = idimg(n)
            jbcinfo(n,nseg,4,1) = 1
            jbcinfo(n,nseg,5,1) = kdimg(n)
         end if
         nbcj0(n) = nseg
      end do
      do n=1,ngrid
         nseg = 0
         do nn=1,maxseg
            do m=1,7
               jbcinfo(n,nn,m,2) = 0
            end do
         end do
         do nn=1,nbc0(n)
            if (ibtyp0(nn,n) .gt. -999) then
            if (ibdir0(nn,n).eq.-1) then
            nseg = nseg + 1
            call bctyper(maxbl,maxseg,n,nn,ibtyp0,lbcs0,lbce0,jbcs0,
     .                   jbce0,kbcs0,kbce0,ibdir0,jbcinfo,2,nseg,
     .                   idimg,jdimg,kdimg,bcvali,bcvalj,bcvalk,
     .                   bcfile0,bcpar10,bcfilei,bcfilej,bcfilek)
            end if
            end if
         end do
c        if no segments have been set, it must be overset
         if (nseg .eq. 0) then
            nseg = nseg + 1
            jbcinfo(n,nseg,1,2) = 0
            jbcinfo(n,nseg,2,2) = 1
            jbcinfo(n,nseg,3,2) = idimg(n)
            jbcinfo(n,nseg,4,2) = 1
            jbcinfo(n,nseg,5,2) = kdimg(n)
         end if
         nbcjdim(n) = nseg
      end do
      do n=1,ngrid
         nseg = 0
         do nn=1,maxseg
            do m=1,7
               kbcinfo(n,nn,m,1) = 0
            end do
         end do
         do nn=1,nbc0(n)
            if (ibtyp0(nn,n) .gt. -999) then
            if (ibdir0(nn,n).eq.2) then
            nseg = nseg + 1
            call bctyper(maxbl,maxseg,n,nn,ibtyp0,lbcs0,lbce0,jbcs0,
     .                   jbce0,kbcs0,kbce0,ibdir0,kbcinfo,1,nseg,
     .                   idimg,jdimg,kdimg,bcvali,bcvalj,bcvalk,
     .                   bcfile0,bcpar10,bcfilei,bcfilej,bcfilek)
            end if
            end if
         end do
c        if no segments have been set, it must be overset
         if (nseg .eq. 0) then
            nseg = nseg + 1
            kbcinfo(n,nseg,1,1) = 0
            kbcinfo(n,nseg,2,1) = 1
            kbcinfo(n,nseg,3,1) = idimg(n)
            kbcinfo(n,nseg,4,1) = 1
            kbcinfo(n,nseg,5,1) = jdimg(n)
         end if
         nbck0(n) = nseg
      end do
      do n=1,ngrid
         nseg = 0
         do nn=1,maxseg
            do m=1,7
               kbcinfo(n,nn,m,2) = 0
            end do
         end do
         do nn=1,nbc0(n)
            if (ibtyp0(nn,n) .gt. -999) then
            if (ibdir0(nn,n).eq.-2) then
            nseg = nseg + 1
            call bctyper(maxbl,maxseg,n,nn,ibtyp0,lbcs0,lbce0,jbcs0,
     .                   jbce0,kbcs0,kbce0,ibdir0,kbcinfo,2,nseg,
     .                   idimg,jdimg,kdimg,bcvali,bcvalj,bcvalk,
     .                   bcfile0,bcpar10,bcfilei,bcfilej,bcfilek)
            end if
            end if
         end do
c        if no segments have been set, it must be overset
         if (nseg .eq. 0) then
            nseg = nseg + 1
            kbcinfo(n,nseg,1,2) = 0
            kbcinfo(n,nseg,2,2) = 1
            kbcinfo(n,nseg,3,2) = idimg(n)
            kbcinfo(n,nseg,4,2) = 1
            kbcinfo(n,nseg,5,2) = jdimg(n)
         end if
         nbckdim(n) = nseg
      end do
c
c     extract 1-1 wake-cut info before resetting any segments with
c     "bcinfo" type -1 to the correct 1-1 interface type 0
c
      nbli = 0
      do n=1,ngrid
         ll = 0
         do nseg=1,nbci0(n)
            if (ibcinfo(n,nseg,1,1) .eq. -1) then
               ll   = ll + 1
               if (ll .eq. 1) nbli = nbli + 1
               nblk(ll,nbli) = n
               limblk(ll,1,nbli) = 1
               limblk(ll,2,nbli) = ibcinfo(n,nseg,2,1)
               limblk(ll,3,nbli) = ibcinfo(n,nseg,4,1)
               limblk(ll,4,nbli) = 1
               limblk(ll,5,nbli) = ibcinfo(n,nseg,3,1)
               limblk(ll,6,nbli) = ibcinfo(n,nseg,5,1)
               isva(ll,1,nbli)   = 2
               isva(ll,2,nbli)   = 3
               ibcinfo(n,nseg,1,1) = 0
            end if 
         end do
      end do
      do n=1,ngrid
         ll = 0
         do nseg=1,nbcidim(n)
            if (ibcinfo(n,nseg,1,2) .eq. -1) then
               ll   = ll + 1
               if (ll .eq. 1) nbli = nbli + 1
               nblk(ll,nbli) = n
               limblk(ll,1,nbli) = idimg(n)
               limblk(ll,2,nbli) = ibcinfo(n,nseg,2,2)
               limblk(ll,3,nbli) = ibcinfo(n,nseg,4,2)
               limblk(ll,4,nbli) = idimg(n) 
               limblk(ll,5,nbli) = ibcinfo(n,nseg,3,2)
               limblk(ll,6,nbli) = ibcinfo(n,nseg,5,2)
               isva(ll,1,nbli)   = 2
               isva(ll,2,nbli)   = 3
               ibcinfo(n,nseg,1,2) = 0
            end if 
         end do
      end do
c
      do n=1,ngrid
         ll = 0
         do nseg=1,nbcj0(n)
            if (jbcinfo(n,nseg,1,1) .eq. -1) then
               ll   = ll + 1
               if (ll .eq. 1) nbli = nbli + 1
               nblk(ll,nbli) = n
               limblk(ll,1,nbli) = jbcinfo(n,nseg,2,1)
               limblk(ll,2,nbli) = 1
               limblk(ll,3,nbli) = jbcinfo(n,nseg,4,1)
               limblk(ll,4,nbli) = jbcinfo(n,nseg,3,1)
               limblk(ll,5,nbli) = 1
               limblk(ll,6,nbli) = jbcinfo(n,nseg,5,1)
               isva(ll,1,nbli)   = 1
               isva(ll,2,nbli)   = 3
               jbcinfo(n,nseg,1,1) = 0
            end if 
         end do
      end do
      do n=1,ngrid
         ll = 0
         do nseg=1,nbcjdim(n)
            if (jbcinfo(n,nseg,1,2) .eq. -1) then
               ll   = ll + 1
               if (ll .eq. 1) nbli = nbli + 1
               nblk(ll,nbli) = n
               limblk(ll,1,nbli) = jbcinfo(n,nseg,2,2)
               limblk(ll,2,nbli) = jdimg(n)
               limblk(ll,3,nbli) = jbcinfo(n,nseg,4,2)
               limblk(ll,4,nbli) = jbcinfo(n,nseg,3,2)
               limblk(ll,5,nbli) = jdimg(n) 
               limblk(ll,6,nbli) = jbcinfo(n,nseg,5,2)
               isva(ll,1,nbli)   = 1
               isva(ll,2,nbli)   = 3
               jbcinfo(n,nseg,1,2) = 0
            end if 
         end do
      end do
c
      do n=1,ngrid
         ll = 0
         do nseg=1,nbck0(n)
            if (kbcinfo(n,nseg,1,1) .eq. -1) then
               ll   = ll + 1
               if (ll .eq. 1) nbli = nbli + 1
               nblk(ll,nbli) = n
               limblk(ll,1,nbli) = kbcinfo(n,nseg,2,1)
               limblk(ll,2,nbli) = kbcinfo(n,nseg,4,1)
               limblk(ll,3,nbli) = 1
               limblk(ll,4,nbli) = kbcinfo(n,nseg,3,1)
               limblk(ll,5,nbli) = kbcinfo(n,nseg,5,1)
               limblk(ll,6,nbli) = 1
               isva(ll,1,nbli)   = 1
               isva(ll,2,nbli)   = 2
               kbcinfo(n,nseg,1,1) = 0
            end if 
         end do
      end do
      do n=1,ngrid
         ll = 0
         do nseg=1,nbckdim(n)
            if (kbcinfo(n,nseg,1,2) .eq. -1) then
               ll   = ll + 1
               if (ll .eq. 1) nbli = nbli + 1
               nblk(ll,nbli) = n
               limblk(ll,1,nbli) = kbcinfo(n,nseg,2,2)
               limblk(ll,2,nbli) = kbcinfo(n,nseg,4,2)
               limblk(ll,3,nbli) = kdimg(n)
               limblk(ll,4,nbli) = kbcinfo(n,nseg,3,2)
               limblk(ll,5,nbli) = kbcinfo(n,nseg,5,2)
               limblk(ll,6,nbli) = kdimg(n)
               isva(ll,1,nbli)   = 1
               isva(ll,2,nbli)   = 2
               kbcinfo(n,nseg,1,2) = 0
            end if 
         end do
      end do
c
c     extract periodic BC info before resetting any segments with
c     "bcinfo" type -2 to the correct 1-1 interface type 0
c
      do n=1,ngrid
         do nseg=1,nbci0(n)
            if (ibcinfo(n,nseg,1,1) .eq. -2) then
               nbli = nbli + 1
               nblk(1,nbli) = n
               limblk(1,1,nbli) = 1
               limblk(1,2,nbli) = ibcinfo(n,nseg,2,1)
               limblk(1,3,nbli) = ibcinfo(n,nseg,4,1)
               limblk(1,4,nbli) = 1
               limblk(1,5,nbli) = ibcinfo(n,nseg,3,1)
               limblk(1,6,nbli) = ibcinfo(n,nseg,5,1)
               isva(1,1,nbli)   = 2
               isva(1,2,nbli)   = 3
c
               limblk(2,1,nbli) = idimg(n)
               limblk(2,2,nbli) = ibcinfo(n,nseg,2,1)
               limblk(2,3,nbli) = ibcinfo(n,nseg,4,1)
               limblk(2,4,nbli) = idimg(n)
               limblk(2,5,nbli) = ibcinfo(n,nseg,3,1)
               limblk(2,6,nbli) = ibcinfo(n,nseg,5,1)
               isva(2,1,nbli)   = 2
               isva(2,2,nbli)   = 3
c
               ibcinfo(n,nseg,1,1) = 0
               ibcinfo(n,nseg,1,2) = 0
            end if
         end do
      end do
c
      do n=1,ngrid
         do nseg=1,nbcj0(n)
            if (jbcinfo(n,nseg,1,1) .eq. -2) then
               nbli = nbli + 1
               nblk(1,nbli) = n
               limblk(1,1,nbli) = jbcinfo(n,nseg,2,1)
               limblk(1,2,nbli) = 1
               limblk(1,3,nbli) = jbcinfo(n,nseg,4,1)
               limblk(1,4,nbli) = jbcinfo(n,nseg,3,1)
               limblk(1,5,nbli) = 1
               limblk(1,6,nbli) = jbcinfo(n,nseg,5,1)
               isva(1,1,nbli)   = 1
               isva(1,2,nbli)   = 3
c
               limblk(2,1,nbli) = jbcinfo(n,nseg,2,1)
               limblk(2,2,nbli) = jdimg(n)
               limblk(2,3,nbli) = jbcinfo(n,nseg,4,1)
               limblk(2,4,nbli) = jbcinfo(n,nseg,3,1)
               limblk(2,5,nbli) = jdimg(n)
               limblk(2,6,nbli) = jbcinfo(n,nseg,5,1)
               isva(2,1,nbli)   = 1
               isva(2,2,nbli)   = 3
c
               jbcinfo(n,nseg,1,1) = 0
               jbcinfo(n,nseg,1,2) = 0
            end if
         end do
      end do
c
      do n=1,ngrid
         do nseg=1,nbck0(n)
            if (kbcinfo(n,nseg,1,1) .eq. -2) then
               nbli = nbli + 1
               nblk(1,nbli) = n
               limblk(1,1,nbli) = kbcinfo(n,nseg,2,1)
               limblk(1,2,nbli) = kbcinfo(n,nseg,4,1)
               limblk(1,3,nbli) = 1
               limblk(1,4,nbli) = kbcinfo(n,nseg,3,1)
               limblk(1,5,nbli) = kbcinfo(n,nseg,5,1)
               limblk(1,6,nbli) = 1
               isva(1,1,nbli)   = 1
               isva(1,2,nbli)   = 2
c
               nblk(2,nbli) = n
               limblk(2,1,nbli) = kbcinfo(n,nseg,2,1)
               limblk(2,2,nbli) = kbcinfo(n,nseg,4,1)
               limblk(2,3,nbli) = kdimg(n)
               limblk(2,4,nbli) = kbcinfo(n,nseg,3,1)
               limblk(2,5,nbli) = kbcinfo(n,nseg,5,1)
               limblk(2,6,nbli) = kdimg(n)
               isva(2,1,nbli)   = 1
               isva(2,2,nbli)   = 2
c
               kbcinfo(n,nseg,1,1) = 0
               kbcinfo(n,nseg,1,2) = 0
            end if
         end do
      end do

      if (nbli.gt.maxseg) then
         write(6,*)'stopping...maxseg too small for 1-1 data...'
         write(6,*)'increase maxseg to at least', nbli
      end if
c
      nobc = 0
      write(15,90)
      do n=1,ngrid
         iovrlp(n) = 1
         write(15,36) n,nbci0(n),nbcidim(n),nbcj0(n),
     .                nbcjdim(n),nbck0(n),nbckdim(n),
     .                iovrlp(n)
      end do
      write(15,102)
      do n=1,ngrid
         do nseg=1,nbci0(n)
            ndata = ibcinfo(n,nseg,7,1)
            if (ibcinfo(n,nseg,1,1) .ne. -99) then
               write(15,'(8i10)') n,nseg,ibcinfo(n,nseg,1,1),
     .                        ibcinfo(n,nseg,2,1),
     .                        ibcinfo(n,nseg,3,1),
     .                        ibcinfo(n,nseg,4,1),
     .                        ibcinfo(n,nseg,5,1),ndata
            else
               nobc = nobc + 1
               write(15,'(2i10,a10,5i10)') n,nseg,qm,
     .                        ibcinfo(n,nseg,2,1),
     .                        ibcinfo(n,nseg,3,1),
     .                        ibcinfo(n,nseg,4,1),
     .                        ibcinfo(n,nseg,5,1),ndata
            end if
            if (ndata .gt. 0) then
               call getdhdr(datahdr,ibcinfo(n,nseg,1,1),ndata)
               write(15,1210) (datahdr(mm),mm=1,ndata)
               write(15,1212) (bcvali(n,nseg,mm,1),mm=1,ndata)
            else if (ndata .lt. 0) then
               write(15,'(''bc file name'')')
               write(15,'(a60)') bcfilei(n,nseg,1)
            end if
         end do
      end do
      write(15,103)
      do n=1,ngrid
         do nseg=1,nbcidim(n)
            ndata = ibcinfo(n,nseg,7,2)
            if (ibcinfo(n,nseg,1,2) .ne. -99) then
               write(15,'(8i10)') n,nseg,ibcinfo(n,nseg,1,2),
     .                        ibcinfo(n,nseg,2,2),
     .                        ibcinfo(n,nseg,3,2),
     .                        ibcinfo(n,nseg,4,2),
     .                        ibcinfo(n,nseg,5,2),ndata
            else
               nobc = nobc + 1
               write(15,'(2i10,a10,5i10)') n,nseg,qm,
     .                        ibcinfo(n,nseg,2,2),
     .                        ibcinfo(n,nseg,3,2),
     .                        ibcinfo(n,nseg,4,2),
     .                        ibcinfo(n,nseg,5,2),ndata
            end if
            if (ndata .gt. 0) then
               call getdhdr(datahdr,ibcinfo(n,nseg,1,2),ndata)
               write(15,1210) (datahdr(mm),mm=1,ndata)
               write(15,1212) (bcvali(n,nseg,mm,2),mm=1,ndata)
            else if (ndata .lt. 0) then
               write(15,'(''bc file name'')')
               write(15,'(a60)') bcfilei(n,nseg,2)
            end if
         end do
      end do
      write(15,112)
      do n=1,ngrid
         do nseg=1,nbcj0(n)
            ndata = jbcinfo(n,nseg,7,1)
            if (jbcinfo(n,nseg,1,1) .ne. -99) then
               write(15,'(8i10)') n,nseg,jbcinfo(n,nseg,1,1),
     .                        jbcinfo(n,nseg,2,1),
     .                        jbcinfo(n,nseg,3,1),
     .                        jbcinfo(n,nseg,4,1),
     .                        jbcinfo(n,nseg,5,1),ndata
            else
               nobc = nobc + 1
               write(15,'(2i10,a10,5i10)') n,nseg,qm,
     .                        jbcinfo(n,nseg,2,1),
     .                        jbcinfo(n,nseg,3,1),
     .                        jbcinfo(n,nseg,4,1),
     .                        jbcinfo(n,nseg,5,1),ndata
            end if
            if (ndata .gt. 0) then
               call getdhdr(datahdr,jbcinfo(n,nseg,1,1),ndata)
               write(15,1210) (datahdr(mm),mm=1,ndata)
               write(15,1212) (bcvalj(n,nseg,mm,1),mm=1,ndata)
            else if (ndata .lt. 0) then
               write(15,'(''bc file name'')')
               write(15,'(a60)') bcfilej(n,nseg,1)
            end if
         end do
      end do
      write(15,113)
      do n=1,ngrid
         do nseg=1,nbcjdim(n)
            ndata = jbcinfo(n,nseg,7,2)
            if (jbcinfo(n,nseg,1,2) .ne. -99) then
               write(15,'(8i10)') n,nseg,jbcinfo(n,nseg,1,2),
     .                        jbcinfo(n,nseg,2,2),
     .                        jbcinfo(n,nseg,3,2),
     .                        jbcinfo(n,nseg,4,2),
     .                        jbcinfo(n,nseg,5,2),ndata
            else
               nobc = nobc + 1
               write(15,'(2i10,a10,5i10)') n,nseg,qm,
     .                        jbcinfo(n,nseg,2,2),
     .                        jbcinfo(n,nseg,3,2),
     .                        jbcinfo(n,nseg,4,2),
     .                        jbcinfo(n,nseg,5,2),ndata
            end if
            if (ndata .gt. 0) then
               call getdhdr(datahdr,jbcinfo(n,nseg,1,2),ndata)
               write(15,1210) (datahdr(mm),mm=1,ndata)
               write(15,1212) (bcvalj(n,nseg,mm,2),mm=1,ndata)
            else if (ndata .lt. 0) then
               write(15,'(''bc file name'')')
               write(15,'(a60)') bcfilej(n,nseg,2)
            end if
         end do
      end do
      write(15,122)
      do n=1,ngrid
         do nseg=1,nbck0(n)
            ndata = kbcinfo(n,nseg,7,1)
            if (kbcinfo(n,nseg,1,1) .ne. -99) then
               write(15,'(8i10)') n,nseg,kbcinfo(n,nseg,1,1),
     .                        kbcinfo(n,nseg,2,1),
     .                        kbcinfo(n,nseg,3,1),
     .                        kbcinfo(n,nseg,4,1),
     .                        kbcinfo(n,nseg,5,1),ndata
            else
               nobc = nobc + 1
               write(15,'(2i10,a10,5i10)') n,nseg,qm,
     .                        kbcinfo(n,nseg,2,1),
     .                        kbcinfo(n,nseg,3,1),
     .                        kbcinfo(n,nseg,4,1),
     .                        kbcinfo(n,nseg,5,1),ndata
            end if
            if (ndata .gt. 0) then
               call getdhdr(datahdr,kbcinfo(n,nseg,1,1),ndata)
               write(15,1210) (datahdr(mm),mm=1,ndata)
               write(15,1212) (bcvalk(n,nseg,mm,1),mm=1,ndata)
            else if (ndata .lt. 0) then
               write(15,'(''bc file name'')')
               write(15,'(a60)') bcfilek(n,nseg,1)
            end if
         end do
      end do
      write(15,123)
      do n=1,ngrid
         do nseg=1,nbckdim(n)
            ndata = kbcinfo(n,nseg,7,2)
            if (kbcinfo(n,nseg,1,2) .ne. -99) then
               write(15,'(8i10)') n,nseg,kbcinfo(n,nseg,1,2),
     .                        kbcinfo(n,nseg,2,2),
     .                        kbcinfo(n,nseg,3,2),
     .                        kbcinfo(n,nseg,4,2),
     .                        kbcinfo(n,nseg,5,2),ndata
            else
               nobc = nobc + 1
               write(15,'(2i10,a10,5i10)') n,nseg,qm,
     .                        kbcinfo(n,nseg,2,2),
     .                        kbcinfo(n,nseg,3,2),
     .                        kbcinfo(n,nseg,4,2),
     .                        kbcinfo(n,nseg,5,2),ndata
            end if
            if (ndata .gt. 0) then
               call getdhdr(datahdr,kbcinfo(n,nseg,1,2),ndata)
               write(15,1210) (datahdr(mm),mm=1,ndata)
               write(15,1212) (bcvalk(n,nseg,mm,2),mm=1,ndata)
            else if (ndata .lt. 0) then
               write(15,'(''bc file name'')')
               write(15,'(a60)') bcfilek(n,nseg,2)
            end if
         end do
      end do
   90 format(6x,4hgrid,5x,5hnbci0,3x,7hnbcidim,5x,5hnbcj0,3x,7hnbcjdim,
     .                 5x,5hnbck0,3x,7hnbckdim,4x,6hiovrlp)
  102 format(3hi0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hjsta,
     .             6x,4hjend,6x,4hksta,6x,4hkend,5x,5hndata)
  103 format(5hidim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hjsta,
     .               6x,4hjend,6x,4hksta,6x,4hkend,5x,5hndata)
  112 format(3hj0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .             6x,4hiend,6x,4hksta,6x,4hkend,5x,5hndata)
  113 format(5hjdim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .             6x,4hiend,6x,4hksta,6x,4hkend,5x,5hndata)
  122 format(3hk0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .             6x,4hiend,6x,4hjsta,6x,4hjend,5x,5hndata)
  123 format(5hkdim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .               6x,4hiend,6x,4hjsta,6x,4hjend,5x,5hndata)
 1210 format(10a10)
 1212 format(5f10.5,2e10.3)
 1213 format(17h     bc data file)
c
c     sequencing/multigrid info - cfl3d can't do mesh sequencing
c     on overset grids at this time
c
      mseq   = 1
      mgflag = 1
      iconsf = 1
      mtt    = 0
      ngam   = 2
      write(15,131)
  131 format(6x,4hmseq,4x,6hmgflag,4x,6hiconsf,7x,3hmtt,6x,4hngam)
      write(15,361) mseq,mgflag,iconsf,mtt,ngam
  361 format(5i10)
c
c     smoothing data
c
      issc = 0
      issr = 0
      do nn=1,3
         epsssc(nn) = 0.3333
         epsssr(nn) = 0.3333
      end do
      write(15,132)
  132 format(6x,4hissc,1x,9hepsssc(1),1x,9hepsssc(2),1x,9hepsssc(3),
     .       6x,4hissr,1x,9hepsssr(1),1x,9hepsssr(2),1x,9hepsssr(3))
      write(15,360) issc,(epsssc(n),n=1,3),issr,(epsssr(n),n=1,3)
  360 format(2(i10,3f10.4))
c
c     iteration data
c
      write(15,136)
  136 format(6x,4hncyc,4x,6hmglevg,5x,5hnemgl,5x,5hnitfo)
      do m=1,mseq
         write(15,36) max(nsteps,1),nglvl,0,0
      end do 
c
c     mitl() sequencing data
c
      write(15,140)
  140 format(6x,4hmit1,6x,4hmit2,6x,4hmit3,6x,4hmit4,6x,9hmit5  ...)
      do 142 m=1,mseq
      write(15,36)(1,i=1,nglvl)
  142 continue
c
c     1-1 data
c
      write(15,144)
      write(15,146)
      write(15,37) nbli
      write(15,1521)
      if (nbli.gt.0) then
          mm = 0
          do 1541 n=1,nbli
             mm = mm + 1
             write(15,37)mm,nblk(1,n),
     .       (limblk(1,l,n),l=1,6),(isva(1,ind,n),ind=1,2)
 1541    continue
      end if
      write(15,1521)
      if (nbli.gt.0) then
         mm = 0
         do 1571 n=1,nbli
            mm = mm + 1
            write(15,37)mm,nblk(2,n),
     .      (limblk(2,l,n),l=1,6),(isva(2,ind,n),ind=1,2)
 1571    continue
      end if
c
  144 format(19h 1-1 blocking data:)
  146 format(6x,4hnbli)
 1521 format(2x,6hnumber,4x,4hgrid,4x,4hista,4x,4hjsta,4x,4hksta,
     .       4x,4hiend,4x,4hjend,4x,4hkend,
     .       3x,5hisva1,3x,5hisva2)
   37 format(10i8)
c
c     patched grid data - no counterpart in overflow
c
      ninter = 0
      write(15,158)
  158 format(22h patch interface data:)
      write(15,160)
  160 format(4x,6hninter)
      write(15,36) ninter
c
c     plot3d output data
c
      iplt3dtyp = 0
      write(15,'(''plot3d output:'')')
      write(15,162)
  162 format(2x,4hgrid,1x,5hiptyp,2x,4hista,
     .2x,4hiend,2x,4hiinc,2x,4hjsta,2x,4hjend,2x,
     .4hjinc,2x,4hksta,2x,4hkend,2x,4hkinc)
      do nn=1,abs(nplot3d)
         write(15,76) nn,iplt3dtyp,0,0,0,0,0,0,0,0,0
      end do
   76 format(11i6)
c
c     movie flag (no movie)
c
      movie = 0
      write(15,'('' movie'')')
      write(15,'(i6)') movie

c
c     print out data (none)
c
      write(15,'(''print out:'')')
      write(15,162)
c
c     control surface data (none)
c
      ncs = 0
      write(15,'(''control surfaces:'')')
      write(15,'(''  ncs'')')
      write(15,166)ncs
  166 format(i6)
      write(15,168)
  168 format(2x,4hgrid,2x,4hista,2x,4hiend,2x,4hjsta,2x,
     .       4hjend,2x,4hksta,2x,4hkend,1x,5hiwall,1x,5hinorm)
c
      return
      end
c
      subroutine bctyper(maxbl,maxseg,n,nn,ibtyp0,lbcs0,lbce0,jbcs0,
     .                   jbce0,kbcs0,kbce0,ibdir0,mbcinfo,mm,nseg,
     .                   idimg,jdimg,kdimg,bcvali,bcvalj,bcvalk,
     .                   bcfile0,bcpar10,bcfilei,bcfilej,bcfilek)
c
      character*80 bcfile0(maxseg,maxbl)
      character*80 bcfilei(maxbl,maxseg,2),bcfilej(maxbl,maxseg,2),
     .             bcfilek(maxbl,maxseg,2)
      dimension mbcinfo(maxbl,maxseg,7,2),idimg(maxbl),jdimg(maxbl),
     .          kdimg(maxbl),bcvali(maxbl,maxseg,7,2),
     .          bcvalj(maxbl,maxseg,7,2),bcvalk(maxbl,maxseg,7,2)
      dimension ibtyp0(maxseg,maxbl),
     .          ibdir0(maxseg,maxbl),jbcs0(maxseg,maxbl),
     .          jbce0(maxseg,maxbl),kbcs0(maxseg,maxbl),
     .          kbce0(maxseg,maxbl),lbcs0(maxseg,maxbl),
     .          lbce0(maxseg,maxbl),bcpar10(maxseg,maxbl)
      dimension bcdat(7)
c
      ndata = 0
      do l= 1,7
         bcdat(l) = 0
      end do
c
c     map overflow bc types onto cfl3d bctypes - note: not all
c     overflow bctypes have an equivalent cfl3d bc type
c
      if (ibtyp0(nn,n) .eq. 0) then
         ibeecee = 0                      !overset
      else if (ibtyp0(nn,n) .eq. 1) then
         ibeecee = 1005                   !invisid wall, adiabatic
      else if (ibtyp0(nn,n) .eq. 2) then
         ibeecee = 1006                   !invis., adiab., norm. mom.
      else if (ibtyp0(nn,n) .eq. 5) then  !viscous wall
         ibeecee  = 2004
         ndata    = 2
         bcdat(1) = 0.                       !adiabatic wall
         bcdat(2) = 0.
      else if (ibtyp0(nn,n) .eq. 7) then
         ibeecee  = 2004 
         ndata    = 2
         bcdat(1) = -1.                      !const. temp. wal
         bcdat(1) = 0.
c     The following is now switched to 1-to-1, instead of 2005
c     else if (ibtyp0(nn,n) .eq. 10) then
c        ibeecee  = 2005
c        ndata    = 4
c        bcdat(1) = real(n)
c        bcdat(2) = 0.
c        bcdat(3) = 0.
c        bcdat(4) = 0.
      else if (ibtyp0(nn,n) .eq. 14) then
         ibeecee = 1013                   !generic singular axis
      else if (ibtyp0(nn,n) .eq. 15) then 
         ibeecee = 1013                   !generic singular axis
      else if (ibtyp0(nn,n) .eq. 16) then
         ibeecee = 1013                   !generic singular axis
      else if (ibtyp0(nn,n) .eq. 17) then
         ibeecee = 1001                   !symmetry
      else if (ibtyp0(nn,n) .eq. 31) then
         ibeecee = 1003                   !characteristic
      else if (ibtyp0(nn,n) .eq. 33) then !specified pressure
         ibeecee  = 2002
         ndata    = 1
         bcdat(1) = bcpar10(nn,n)            !press. ratio
      else if (ibtyp0(nn,n) .eq. 40) then
         ibeecee = 1003                   !freesream
      else if (ibtyp0(nn,n) .eq. 42) then !prescribed q
         ibeecee  = 2007
         ndata    = -5                       !data file with q
      else if (ibtyp0(nn,n) .eq. 47) then
         ibeecee = 1003                   !characteristic/freestream
c
c     c-cuts require a temporary non-zero bc flag to distinguish
c     them from overset boundaries so that 1-1 connectivity can be 
c     be extracted later on
c
      else if (ibtyp0(nn,n) .eq. 51) then
         ibeecee = -1    !temp flag for c-cut
      else if (ibtyp0(nn,n) .eq. 52) then
         ibeecee = -1    !temp flag for c-cut
      else if (ibtyp0(nn,n) .eq. 53) then
         ibeecee = -1    !temp flag for c-cut
      else if (ibtyp0(nn,n) .eq. 57) then
         ibeecee = -1    !temp flag for c-cut
      else if (ibtyp0(nn,n) .eq. 58) then
         ibeecee = -1    !temp flag for c-cut
      else if (ibtyp0(nn,n) .eq. 59) then
         ibeecee = -1    !temp flag for c-cut
      else if (ibtyp0(nn,n) .eq. 10) then
         ibeecee = -2    !temp flag for periodic
c
c     flag any bc not set this point to -99. This bc type will be 
c     output as a ? in the cfl3d input file, requiring user 
c     interpretation
c
      else
         ibeecee = -99
      end if
c
c     set cfl3d bc segment arrays
c
c     j-boundaries
      if (ibdir0(nn,n).eq.1 .or. ibdir0(nn,n).eq.-1) then
          mbcinfo(n,nseg,1,mm) = ibeecee
          mbcinfo(n,nseg,4,mm) = kbcs0(nn,n)
          mbcinfo(n,nseg,5,mm) = kbce0(nn,n)
          mbcinfo(n,nseg,2,mm) = lbcs0(nn,n)
          mbcinfo(n,nseg,3,mm) = lbce0(nn,n)
          mbcinfo(n,nseg,7,mm) = ndata
          if (ndata .gt. 0) then
             do l=1,ndata
                bcvalj(n,nseg,l,mm) = bcdat(l)
             end do
          else if (ndata .lt. 0) then
             bcfilej(n,nseg,mm) = bcfile0(nn,n)
          end if
c         must expand c-type bc's to both sides of wake cut
          if (ibeecee .eq. -1) then
             nseg = nseg + 1
             if (ibtyp0(nn,n).eq.52 .or. ibtyp0(nn,n).eq.58) then
c               k is c-direction
                mbcinfo(n,nseg,1,mm) = ibeecee
                mbcinfo(n,nseg,4,mm) = kdimg(n) - kbcs0(nn,n) + 1
                mbcinfo(n,nseg,5,mm) = kdimg(n) - kbce0(nn,n) + 1
                mbcinfo(n,nseg,2,mm) = lbcs0(nn,n)
                mbcinfo(n,nseg,3,mm) = lbce0(nn,n)
             end if
             if (ibtyp0(nn,n).eq.53 .or. ibtyp0(nn,n).eq.59) then
c               l is c-direction
                mbcinfo(n,nseg,1,mm) = ibeecee
                mbcinfo(n,nseg,4,mm) = kbcs0(nn,n)
                mbcinfo(n,nseg,5,mm) = kbce0(nn,n)
                mbcinfo(n,nseg,2,mm) = idimg(n) - lbcs0(nn,n) + 1
                mbcinfo(n,nseg,3,mm) = idimg(n) - lbce0(nn,n) + 1
             end if
             mbcinfo(n,nseg,7,mm) = ndata
             if (ndata .gt. 0) then
                do l=1,ndata
                   bcvalj(n,nseg,l,mm) = bcdat(l)
                end do
             else if (ndata .lt. 0) then
                bcfilej(n,nseg,mm) = bcfile0(nn,n)
             end if
          end if
      end if
c     k-boundaries
      if (ibdir0(nn,n).eq.2 .or. ibdir0(nn,n).eq.-2) then
          mbcinfo(n,nseg,1,mm) = ibeecee
          mbcinfo(n,nseg,4,mm) = jbcs0(nn,n)
          mbcinfo(n,nseg,5,mm) = jbce0(nn,n)
          mbcinfo(n,nseg,2,mm) = lbcs0(nn,n)
          mbcinfo(n,nseg,3,mm) = lbce0(nn,n)
          mbcinfo(n,nseg,7,mm) = ndata
          if (ndata .gt. 0) then
             do l=1,ndata
                bcvalk(n,nseg,l,mm) = bcdat(l)
             end do
          else if (ndata .lt. 0) then
             bcfilek(n,nseg,mm) = bcfile0(nn,n)
          end if
c         must expand c-type bc's to both sides of wake cut
          if (ibeecee .eq. -1) then
             nseg = nseg + 1
             if (ibtyp0(nn,n).eq.51 .or. ibtyp0(nn,n).eq.57) then
c               j is c-direction
                mbcinfo(n,nseg,1,mm) = ibeecee
                mbcinfo(n,nseg,4,mm) = jdimg(n) - jbcs0(nn,n) + 1
                mbcinfo(n,nseg,5,mm) = jdimg(n) - jbce0(nn,n) + 1
                mbcinfo(n,nseg,2,mm) = lbcs0(nn,n)
                mbcinfo(n,nseg,3,mm) = lbce0(nn,n)
             end if
             if (ibtyp0(nn,n).eq.53 .or. ibtyp0(nn,n).eq.59) then
c               l is c-direction
                mbcinfo(n,nseg,1,mm) = ibeecee
                mbcinfo(n,nseg,4,mm) = jbcs0(nn,n)
                mbcinfo(n,nseg,5,mm) = jbce0(nn,n)
                mbcinfo(n,nseg,2,mm) = idimg(n) - lbcs0(nn,n) + 1
                mbcinfo(n,nseg,3,mm) = idimg(n) - lbce0(nn,n) + 1
             end if
             mbcinfo(n,nseg,7,mm) = ndata
             if (ndata .gt. 0) then
                do l=1,ndata
                   bcvalk(n,nseg,l,mm) = bcdat(l)
                end do
             else if (ndata .lt. 0) then
               bcfilek(n,nseg,mm) = bcfile0(nn,n)
             end if
          end if
      end if
c     i-boundaries
      if (ibdir0(nn,n).eq.3 .or. ibdir0(nn,n).eq.-3) then
          mbcinfo(n,nseg,1,mm) = ibeecee
          mbcinfo(n,nseg,2,mm) = jbcs0(nn,n)
          mbcinfo(n,nseg,3,mm) = jbce0(nn,n)
          mbcinfo(n,nseg,4,mm) = kbcs0(nn,n)
          mbcinfo(n,nseg,5,mm) = kbce0(nn,n)
          mbcinfo(n,nseg,7,mm) = ndata
          if (ndata .gt. 0) then
             do l=1,ndata
                bcvali(n,nseg,l,mm) = bcdat(l)
             end do
          else if (ndata .lt. 0) then
             bcfilei(n,nseg,mm) = bcfile0(nn,n)
          end if
          if (ibeecee .eq. -1) then
             nseg = nseg + 1
             if (ibtyp0(nn,n).eq.51 .or. ibtyp0(nn,n).eq.57) then
c               j is c-direction
                mbcinfo(n,nseg,1,mm) = ibeecee
                mbcinfo(n,nseg,4,mm) = kbcs0(nn,n)
                mbcinfo(n,nseg,5,mm) = kbce0(nn,n)
                mbcinfo(n,nseg,2,mm) = jdimg(n) - jbcs0(nn,n) + 1
                mbcinfo(n,nseg,3,mm) = jdimg(n) - jbce0(nn,n) + 1
             end if
             if (ibtyp0(nn,n).eq.52 .or. ibtyp0(nn,n).eq.58) then
c               k is c-direction
                mbcinfo(n,nseg,1,mm) = ibeecee
                mbcinfo(n,nseg,4,mm) = kdimg(n) - kbcs0(nn,n) + 1
                mbcinfo(n,nseg,5,mm) = kdimg(n) - kbce0(nn,n) + 1
                mbcinfo(n,nseg,2,mm) = jbcs0(nn,n)
                mbcinfo(n,nseg,3,mm) = jbce0(nn,n)
             end if
             mbcinfo(n,nseg,7,mm) = ndata
             if (ndata .gt. 0) then
                do l=1,ndata
                   bcvali(n,nseg,l,mm) = bcdat(l)
                end do
             else if (ndata .lt. 0) then
                bcfilei(n,nseg,mm) = bcfile0(nn,n)
             end if
          end if
      end if
c
      return
      end
c
      subroutine getdhdr(datahdr,ibctyp,ndata)
c***********************************************************************
c     Purpose:  To set character data for main output file headers 
c     when the 2000 series bc's are used 
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*10 datahdr(10)
c
      do 1 m=1,10
    1 datahdr(m)  = '          '
c
      if (ibctyp.eq.2002) then
         datahdr(1)  = '    p/pinf'
      end if
c
      if (ibctyp.eq.2003) then
         datahdr(1)  = '      mach'
         datahdr(2)  = '   pt/pinf'
         datahdr(3)  = '   tt/tinf'
         datahdr(4)  = '     alpha'
         datahdr(5)  = '      beta'
         if (abs(ndata).gt.5) then
            datahdr(6)  = '     turb1'
         end if
         if (abs(ndata).gt.6) then
            datahdr(7)  = '     turb2'
         end if
      end if
c
      if (abs(ibctyp).eq.2004) then
         datahdr(1)  = '   tw/tinf'
         datahdr(2)  = '        cq'
      end if
c
      if (ibctyp.eq.2005) then
         datahdr(1)  = '      nblp'
         datahdr(2)  = '     dthtx'
         datahdr(3)  = '     dthty'
         datahdr(4)  = '     dthtz'
      end if
c
      if (ibctyp.eq.2006) then
         datahdr(1)  = '      nblc'
         datahdr(2)  = '    p/pinf'
         datahdr(3)  = ' integ dir'
         datahdr(4)  = ' axial dir'
      end if
c
      if (ibctyp.eq.2007) then
         datahdr(1)  = '       rho'
         datahdr(2)  = '         u'
         datahdr(3)  = '         v'
         datahdr(4)  = '         w'
         datahdr(5)  = '         p'
         if (abs(ndata).gt.5) then
            datahdr(6)  = '     turb1'
         end if
         if (abs(ndata).gt.6) then
            datahdr(7)  = '     turb2'
         end if
      end if
c
      if (ibctyp.eq.2102) then
         datahdr(1) = '    p/pinf'
         datahdr(2) = '   dp/pinf'
         datahdr(3) = '    rfreqp'
         datahdr(4) = '      lref'
      end if
c
      if (ibctyp.eq.2008) then
         datahdr(1)  = '       rho'
         datahdr(2)  = '         u'
         datahdr(3)  = '         v'
         datahdr(4)  = '         w'
         if (abs(ndata).gt.4) then
            datahdr(5)  = '     turb1'
         end if
         if (abs(ndata).gt.5) then
            datahdr(6)  = '     turb2'
         end if
      end if
c
      return
      end
c
      subroutine rwgrid(jd,kd,id,n,i2d)
c
      dimension x(jd,kd,id),y(jd,kd,id),z(jd,kd,id)
c
      pointer (ip_x,x),(ip_y,y),(ip_z,z)
c
c     allocate memory
c
      memuse = 0
      icall  = 0
      call umalloc_r(ip_x,id*jd*kd,0,icall,memuse)
      call umalloc_r(ip_y,id*jd*kd,0,icall,memuse)
      call umalloc_r(ip_z,id*jd*kd,0,icall,memuse)
c
c     read block in as overflow jki and write as cfl3d ijk
c     (note: as long as input grid is unformatted, can 
c     ignore any iblank data the overflow grid may contain)
c
      write(6,*)
      write(6,'(''reading zone '',i4,'' as j k i = '',3i4)') 
     .  n,jd,kd,id
      write(6,'(''writing zone '',i4,'' as i j k = '',3i4)')
     .  n,id,jd,kd
      read(20) (((x(j,k,i),j=1,jd),k=1,kd),i=1,id),
     .         (((y(j,k,i),j=1,jd),k=1,kd),i=1,id),
     .         (((z(j,k,i),j=1,jd),k=1,kd),i=1,id)
c
c---rtb
c     write(25) (((x(j,k,i),i=1,id),j=1,jd),k=1,kd),
c    .          (((y(j,k,i),i=1,id),j=1,jd),k=1,kd),
c    .          (((z(j,k,i),i=1,id),j=1,jd),k=1,kd)
      idout = id
      if (i2d .ne. 0) idout = 2
      write(25) (((x(j,k,i),i=1,idout),j=1,jd),k=1,kd),
     .          (((y(j,k,i),i=1,idout),j=1,jd),k=1,kd),
     .          (((z(j,k,i),i=1,idout),j=1,jd),k=1,kd)
c
c     free memory
c
      ifree = 1
      if (ifree.gt.0) then
         call ufree_r(ip_x)
         call ufree_r(ip_y)
         call ufree_r(ip_z)
      end if
c
      return
      end

      subroutine umalloc_r(ipoint,n_words,intflag,icall,memuse)
c***********************************************************************
c     Purpose:  universal malloc...allocates memory via fortran 77
c     pointer extension.
c
c     NOTE: in this special-purpose version, allocates memory
c     for REAL (not complex real) arrays, as well as integer arrays.
c
c     icall   - integer counter to aid in debugging
c     intflag - flag to set variable type:
c              -1...real*4 variable
c               0...real variable
c               1...integer variable
c     memuse  - running total of memory that has been allocated
c
c     original version provied courtesy Steven J. Massey, AS&M Inc.
c     modified by r. t. biedron, January, 1999.
c***********************************************************************
c
      parameter(len_i=4,len_r4=4)
c
      integer errcode,abort
c
      dimension dum(1)
c
      pointer (ipoint, dum)
c
#ifdef CRAY_TIME
c     cray_time implies cray
      len_r = 8
#else
#  ifdef DBLE_PRECSN
      len_r = 8
#  else 
      len_r = 4
#  endif
#endif
c---rtb
      len_r = 8
c---rtb
c
      icall = icall+1
c
      if (intflag.eq.1) then
        len = len_i
      else if (intflag.eq.0) then
        len = len_r
      else if (intflag.eq.-1) then
        len = len_r4
      end if
c
      n_bytes = n_words * len
c
#ifdef CRAY_TIME
c     cray_time implies cray
      abort = 0
      errcode = 0
      nwords = n_words
      call hpalloc( ipoint, nwords, errcode, abort)
#else
      ipoint = malloc(n_bytes)
#endif
c
      if (ipoint.eq.0) then
        write(6,'(a,i4)') 
     .  'stopping...malloc failed to allocate the memory on call ',
     .  icall
        if (n_bytes.lt.1e10) then 
          write(6,12) n_words, n_bytes
        else if (n_bytes.lt.1e15) then 
          write(6,13) n_words, n_bytes
        else 
          write(6,*) 'requested more than 10^(15) bytes'
        end if
        write(6,14) memuse
        stop
      else
c        write(6,*) 'call',icall,' umalloc: allocated ',
c    .   n_bytes,' bytes'
      end if
c
  12  format('requested ', i10, ' words, (',i10, ' bytes)')
  13  format('requested ', i15, ' words, (',i15, ' bytes)')
  14  format('prior to this call, total allocation was',
     .        i15, ' bytes')
c
      memuse  = memuse + n_bytes
c
      return
      end

      subroutine ufree_r(ipoint)
c***********************************************************************
c     Purpose:  universal free...deallocates memory via fortran 77
c     pointer extension.
c
c     NOTE: in this special-purpose version, frees REAL (not complex 
c     real) allocation, as well as integer allocation.
c
c     original version provied courtesy Steven J. Massey, AS&M Inc.
c***********************************************************************
c
      integer errcode,abort
c
      pointer (ipoint,dum)
c
#ifdef CRAY_TIME
c     cray_time implies cray
      abort   = 0
      errcode = 0
      call hpdeallc(ipoint,errcode,abort )
#else
      call free(ipoint)
#endif
      return
      end
