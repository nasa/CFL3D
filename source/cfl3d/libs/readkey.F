      subroutine readkey(ititr,myid,ibufdim,nbuf,bou,nou,iunit11,
     .                   ierrflg)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Read in any keyword-driven inputs. ititr is returned
c     as 1 if the title line has been read; returned as 0 if unread
c     note: most, if not all, data read via this routine will have to
c     be explicitly passed from the host to the nodes for parallel 
c     processing, via subroutine trnsfr_vals. 
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*210 inpstr
      character*120 bou(ibufdim,nbuf)
c
      real realval(20)
c
      dimension nou(nbuf)
c
      common /des/ ides,cdes
      common /avgdata/ iteravg,xnumavg
      common /unit5/ iunit5
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /zero/ iexp
      common /singular/ atol
      common /complx/ xmach_img,alpha_img,beta_img,reue_img,tinf_img,
     .                geom_img,surf_img,xrotrate_img,yrotrate_img,
     .                zrotrate_img
      common /cgns/ icgns,iccg,ibase,nzones,nsoluse,irind,jrind,krind
      common /precond/ cprec,uref,avn
      common /alphait/ cltarg,resupdt,cltol,dalim,dcldal,alphalast,
     .                 ialphit,nttlast,icycupdt
      common /turbconv/ cflturb,edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb
      common /entfix/ epsa_l,epsa_r
      common /key/ nkey
      common /is_blockbc/ is_blk(5),ie_blk(5),ivolint
      common /elastic_ss/ idef_ss
      common /bin/ ibin,iblnk,iblnkfr,ip3dgrad
      common /memory/ memadd,memaddi
      common /deformz/ beta1,ismooth,negvol
      common /ghost/ irghost,iwghost
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10,tur20
c
      ititr = 0
      nkey  = 0
c
c*********************************************
c     define default keyword-driven input data
c*********************************************
c
c     ratio of specific heats
      gamma  = 1.4
      nkey   = nkey + 1
c
c     Prandtl number
      pr     = 0.72
      nkey   = nkey + 1
c
c     Turbulent Prandtl number
      prt    = 0.9
      nkey   = nkey + 1
c
c     reference temperature for Sutherland's Law     
      cbar   = 198.6
      nkey   = nkey + 1
c
c     tolerance for collapsed metrics (10.**(-iexp) is machine zero)
      atol   = max(1.e-07,10.**(-iexp+1))
      nkey   = nkey + 1
c
c     complex perturbation to the Mach number
      xmach_img = 0.
      nkey      = nkey + 1
c
c     complex perturbation to angle of attack, alpha
      alpha_img = 0.
      nkey      = nkey + 1
c
c     complex perturbation to yaw angle, beta
      beta_img = 0.
      nkey     = nkey + 1
c
c     complex perturbation to the unit Reynolds number, reue
      reue_img = 0.
      nkey     = nkey + 1
c
c     complex perturbation to the freestream temperature, tinf 
      tinf_img = 0.
      nkey     = nkey + 1 
c
c     size of complex perturbation to grid
c     (this should be consitant with the input complex grid)
c
      geom_img  = 0.
      nkey      = nkey + 1
c
c     CGNS file flag (0 is standard cfl3d files, 1 is CGNS)
c
      icgns     = 0
      nkey      = nkey + 1
c
c     relative amount of preconditioning (0.0 is none, 1.0 is
c     fully on:
c
      cprec     = 0.
      nkey      = nkey + 1
c
c     limiting velocity for preconditioning
c
      uref      = xmach
      nkey      = nkey + 1
c
c     multiplicative factor for uref**2
c
      avn      = 1.0
      nkey     = nkey + 1
c
c     specified Cl (ialphit is NOT a keyword, just a flag that is
c     0 unless the keyword cltarg is not the default value
c
      cltarg   = 99999.0
      ialphit  = 0
      nkey     = nkey + 1
c
c     no. of orders of resid. reduction before alpha is
c     next updated
c
      resupdt  = -40.
      nkey     = nkey + 1
c
c     no. of subiterations in turbulence model per outer iteration
c
      nsubturb = 1
      nkey     = nkey + 1
c
c     factor governing cfl number for turbulence model; supercedes
c     the hardwired value of "factor" in the turb. model only
c     if nonzero
c
      cflturb  = 0.
      nkey     = nkey + 1
c
c     limit on ratio of minimum to maximum |eigenvalue| - the
c     parameter for Harten's entropy fix or some variation 
c     thereof - see fhat for exact details. The value on the
c     RHS is input, and that value for the LHS is taken as a 
c     factor of 2 times the RHS value (hence, epsa_l is NOT a keyword)
c
      epsa_r   = 0.
      epsa_l   = 2.*epsa_r
      nkey     = nkey + 1
c
c     number of cycles over which to freeze turbulence model; the
c     default of zero gives the standard unfrozen treatment.
c
      nfreeze  = 0
      nkey     = nkey + 1
c
c     flag for using/not using the exact volume terms on 1-1 block
c     boundaries (1...use the exact volumes)
c
      ivolint  = 1
      nkey     = nkey + 1
c
c     flag for mesh deformation in steady-state mode
c     default of zero does not deform mesh to fit a new
c     surface shape
c
      idef_ss  = 0
      nkey     = nkey + 1
c
c     flag for writing unformated/formatted plot3d files
c     default is unformatted (except on T3E)
c
      ibin     = 1
#   ifdef T3E
      ibin     = 0
#   endif
      nkey     = nkey + 1
c
c     flag for enabling/disabling the writing of the iblank
c     array in the plot3d grid file. default is enabled
c
      iblnk    = 1
      nkey     = nkey + 1
c
c     flag for enabling/disabling the blanking of fringe points
c     in the plot3d grid file for overset grids. default is
c     endabled
c
      iblnkfr  = 1
      nkey     = nkey + 1
c
c     flag to switch from solution output to derivative output
c     in the plot3d "q" file. default is solution output
c
      ip3dgrad = 0
      nkey     = nkey + 1
c
c     size of complex perturbation to surface grid onto which
c     the volume grid is to be deformed
c     (this should be consitant with the input complex grid)
c
      surf_img  = 0.
      nkey      = nkey + 1
c
c     additional user-specified memory to allocate to the
c     work array sizes (memadd to mwork, memeaddi to mworki)
c     used to offset any underestimation by the sizer routine
c
      memadd    = 0
      nkey      = nkey + 1
      memaddi   = 0
      nkey      = nkey + 1
c
c     flag to turn off stops when negative volumes/bad metrics
c     are encountered - use only to debug mesh deformation!
c     default is to stop when negative volumes/bad metrics are
c     found
c
      negvol    = 0
      nkey      = nkey + 1
c
c     eddy viscosity limiter for two eqn. turbulence models: limit
c     eddy viscosity to edvislim times the laminar viscosity
c
      edvislim  = 100000.
      nkey      = nkey + 1
c
c     flag to set whether approximate production term (0) or full
c     production term (1) is used in EASM models 8,9,13,14
c
      iturbprod = 0
      nkey      = nkey + 1
c
c     flag to read ghost-cell data from restart file (1) or not (0)
c     newer version 6 restart files will have ghost-cell data;
c     restart files from beta version 6 do not, nor do version 5
c     restart files. default is to read ghost-cell data
c
      irghost   = 1
      nkey      = nkey + 1
c
c     flag to write ghost-cell data to restart file (1) or not (0)
c     newer version 6 restart files will have ghost-cell data;
c     restart files from beta version 6 do not, nor do version 5
c     restart files. default is to write ghost-cell data
c
      iwghost   = 1
      nkey      = nkey + 1
c
c     limit of alpha change (deg.) per update
c
      dalim     = 5.0
      nkey      = nkey + 1
c
c     no. of cycles between alpha updates (if set > 0, this
c     method takes precedence over resupdt)
c
      icycupdt  = 0
      nkey      = nkey + 1
c
c     user input value for dcl/dalpha
c
      dcldal    = 6.2831853
      nkey      = nkey + 1
c
c     non-inertial reference frame flag (0=inertial; 1=noninertial)
c
      noninflag = 0
      nkey      = nkey + 1
c
c     rotation center x-coordinate for non-inertial reference frame
c
      xcentrot = 0.
      nkey     = nkey + 1
c
c     rotation center y-coordinate for non-inertial reference frame
c
      ycentrot = 0.
      nkey     = nkey + 1
c
c     rotation center z-coordinate for non-inertial reference frame
c
      zcentrot = 0.
      nkey     = nkey + 1
c
c     rotation rate in x-direction for non-inertial reference frame
c
      xrotrate = 0.
      nkey     = nkey + 1
c
c     rotation rate in y-direction for non-inertial reference frame
c
      yrotrate = 0.
      nkey     = nkey + 1
c
c     rotation rate in z-direction for non-inertial reference frame
c
      zrotrate = 0.
      nkey     = nkey + 1
c
c     complex perturbation to the rotation rate in x-direction for
c     non-inertial reference frame
      xrotrate_img = 0.
      nkey      = nkey + 1
c
c     complex perturbation to the rotation rate in y-direction for
c     non-inertial reference frame
      yrotrate_img = 0.
      nkey      = nkey + 1
c
c     complex perturbation to the rotation rate in z-direction for
c     non-inertial reference frame
      zrotrate_img = 0.
      nkey      = nkey + 1
c
c     flag to read (1) or skip reading (0) 2nd order-time turbulence 
c     terms and dt in restart file (need to skip if using an older 
c     time-accurate-with-2nd-order-time restart file which does not 
c     have them) - 1=default
c
      itime2read = 1
      nkey     = nkey + 1
c
c     flag to control time-accuracy of turbulence model - in the old
c     method (0), turb eqns were 1st order regardless of order of
c     accuracy of mean flow eqns; new default (1) is they are the same
c     order as the mean flow eqns, as defined by the parameter "ita"
c
      itaturb = 1
      nkey     = nkey + 1
c
c     flag to perform DES with turbulence model (1) or not (0)
c
      ides = 0
      nkey     = nkey + 1
c
c     constant associated with DES
c
      cdes = 0.65
      nkey     = nkey + 1
c
c     flag to store iteration-averaged q values in a PLOT3D file:
c     0 = no averaging or storage, 1 = start averaging now,
c     2 = continue averaging from previous run
c
      iteravg = 0
      nkey     = nkey + 1
c
c     turbulent quantity freestream levels
c     < 0  = use the default value (different for each turb model)
c     >= 0 = use this number as the specified user input value
c     (defaults=-1)
      tur10 = -1.
      tur20 = -1.
      nkey     = nkey + 2
c
c*********************************************
c     check for keyword-driven inputs
c*********************************************
c
      read(iunit5,1593) inpstr
      if (inpstr(1:1).eq.'>') then
         go to 1000
      else if (inpstr(1:1).eq.'<') then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1903)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      else
         ititr = 1
         read(inpstr,10)(realval(i),i=1,20)
         do i=1,20
            title(i) = realval(i)
         end do
         return
      end if
c
c     keyword-driven inputs
c
 1000 continue
      if (iunit11.gt.0) write(iunit11,1101)
 1100 continue
      read(iunit5,1593,end=1900) inpstr
      if (inpstr(1:1).eq.'<') then
         if (iunit11.gt.0) write(iunit11,1102)
         return
      end if
c
      if (inpstr(1:1).eq.'/' .or. inpstr(1:1).eq.' ') goto 1100
c
c     echo the keyword character strings
c
      npos =   1
      lc2  =   0
      call parser(inpstr,npos,lc1,lc2,lcl,-1)
      lclw = min(lcl,80)
      if (iunit11.gt.0) write(iunit11,1594) inpstr(1:lclw)
c
      npos =   1
      lc2  =   0
      lcl  = 210
      call parser(inpstr,npos,lc1,lc2,lcl,1)
c
      if (inpstr(lc1:lc2).eq.'gamma') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         gamma = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'pr') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         pr = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'prt') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         prt = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'cbar') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cbar = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'atol') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         atol = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'xmach_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         xmach_img = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'alpha_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         alpha_img = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'beta_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         beta_img = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'reue_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         reue_img = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tinf_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tinf_img = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'geom_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         geom_img = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'icgns') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) icgns
#ifdef CGNS
#else
c
c        cannot use cgns input unless the code has been 
c        installed with -cgnsdir=...
c
         if (icgns.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...must install with '',
     .           ''-cgnsdir=... in order to use _CGNS_ (icgns=1)'')')
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
#endif
#if defined CGNS
#   if defined CMPLX
c
c        cannot use BOTH complex and CGNS options
c
         if (icgns.eq.1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...cannot use BOTH '',
     .           ''_CGNS_ (icgns=1) and complex variables'')')
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
#endif
#endif
c
      else if (inpstr(lc1:lc2).eq.'cprec') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cprec = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'uref') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         uref = realval(1)

c
      else if (inpstr(lc1:lc2).eq.'avn') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         avn = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'cltarg') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cltarg = realval(1)
         if (cltarg.ne.99999.0) then
            ialphit = 1
         end if
c
      else if (inpstr(lc1:lc2).eq.'resupdt') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         resupdt = -abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'nsubturb') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) nsubturb
c
      else if (inpstr(lc1:lc2).eq.'cflturb') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cflturb = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'epsa_r') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         epsa_r = realval(1)
         epsa_l = 2.*epsa_r
c
      else if (inpstr(lc1:lc2).eq.'nfreeze') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) nfreeze
c
      else if (inpstr(lc1:lc2).eq.'ivolint') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ivolint
c
      else if (inpstr(lc1:lc2).eq.'idef_ss') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) idef_ss
c
      else if (inpstr(lc1:lc2).eq.'ibin') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ibin
c
      else if (inpstr(lc1:lc2).eq.'iblnk') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iblnk
c
      else if (inpstr(lc1:lc2).eq.'iblnkfr') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iblnkfr
c
      else if (inpstr(lc1:lc2).eq.'ip3dgrad') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ip3dgrad
c
      else if (inpstr(lc1:lc2).eq.'surf_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         surf_img = realval(1)
c
c        must also have idef_ss option enabled as well
c
         if (idef_ss.eq.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...must have idef_ss = 1 '',
     .           ''if surf_img .ne. 0.'')')
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
c
      else if (inpstr(lc1:lc2).eq.'memadd') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) memadd
c
      else if (inpstr(lc1:lc2).eq.'memaddi') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) memaddi
c
      else if (inpstr(lc1:lc2).eq.'negvol') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) negvol
c
      else if (inpstr(lc1:lc2).eq.'edvislim') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         edvislim = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'iturbprod') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iturbprod
c
      else if (inpstr(lc1:lc2).eq.'irghost') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) irghost
c
      else if (inpstr(lc1:lc2).eq.'iwghost') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iwghost
c
      else if (inpstr(lc1:lc2).eq.'dalim') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         dalim = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'icycupdt') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) icycupdt
c
      else if (inpstr(lc1:lc2).eq.'dcldal') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         dcldal = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'noninflag') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) noninflag
c
      else if (inpstr(lc1:lc2).eq.'xcentrot') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         xcentrot = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'ycentrot') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         ycentrot = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'zcentrot') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         zcentrot = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'xrotrate') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         xrotrate = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'yrotrate') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         yrotrate = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'zrotrate') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         zrotrate = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'xrotrate_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         xrotrate_img = realval(1)
#ifdef CMPLX
         xrotrate = cmplx(real(xrotrate),real(xrotrate_img))
#endif
c
      else if (inpstr(lc1:lc2).eq.'yrotrate_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         yrotrate_img = realval(1)
#ifdef CMPLX
         yrotrate = cmplx(real(yrotrate),real(yrotrate_img))
#endif
c
      else if (inpstr(lc1:lc2).eq.'zrotrate_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         zrotrate_img = realval(1)
#ifdef CMPLX
         zrotrate = cmplx(real(zrotrate),real(zrotrate_img))
#endif
c
      else if (inpstr(lc1:lc2).eq.'itime2read') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) itime2read
c
      else if (inpstr(lc1:lc2).eq.'itaturb') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) itaturb
c
      else if (inpstr(lc1:lc2).eq.'ides') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ides
c
      else if (inpstr(lc1:lc2).eq.'cdes') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cdes = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'iteravg') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iteravg
c
      else if (inpstr(lc1:lc2).eq.'tur10') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur10 = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tur20') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur20 = realval(1)
c
      else
c
         if (iunit11.gt.0) then
            write(iunit11,1595)
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
c
      end if
c
      go to 1100
c
 1900 continue
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1901)
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1902)
      call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
c
 1101 format('>',21('-'),' begin keyword-driven input section ',
     .            21('-'),'>')
 1102 format( '<',22('-'),' end keyword-driven input section ',
     .            22('-'),'<')
c
   10 format(20a4)
 1593 format(a210)
 1594 format(a)
 1595 format('*** STOPPING: The keyword above is not supported. ***')
c
 1901 format(3x,'ERROR: There was no ''<''-line in the input file')
 1902 format(10x,'to exit the keyword-driven input section.')
 1903 format(3x,'ERROR: keyword input must start with ',
     .          'a ''>''-line in the input file')
c
      return
      end
