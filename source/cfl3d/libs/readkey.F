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
      character*241 avgq2,avgq2pert,clcds,clcdp,output_dir,filename
      logical isthere
c
      real realval(20)
c
      dimension nou(nbuf)
c
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /des/ cdes,ides,cddes
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
      common /ginfo2/ lq2avg,iskip_blocks,inc_2d(3),inc_coarse(3)
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
      common /alphait/ ialphit,cltarg,rlxalph,dalim,dalpha,icycupdt
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /entfix/ epsa_l,epsa_r
      common /key/ nkey
      common /is_blockbc/ is_blk(5),ie_blk(5),ivolint
      common /elastic_ss/ idef_ss
      common /bin/ ibin,iblnk,iblnkfr,ip3dgrad
      common /memory/ memadd,memaddi
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt 
      common /ghost/ irghost,iwghost
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /gridtrans/ roll_angle
      common /konew/ ikoprod,isstdenom,pklimterm,ibeta8kzeta,i_bsl,
     .        keepambient,re_thetat0,i_wilcox06,i_wilcox06_chiw,
     .        i_turbprod_kterm,i_catris_kw,prod2d3dtrace,
     .        i_compress_correct,isstsf,i_wilcox98,i_wilcox98_chiw,
     .        isst2003
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
      common /fullns/ ifullns
      common /maxiv/ ivmx
      common /curvat/ isarc2d,sarccr3,ieasmcc2d,isstrc,sstrc_crc,
     .        isar,crot,isarc3d
      common /rigidbody/ irigb,irbtrim
      common /rbstmt2/ tmass,yinert,uinfrb,qinfrb,greflrb,gaccel,crefrb,
     .                 xtmref,areat
      common /filenam2/ avgq2,avgq2pert,clcds,clcdp,output_dir
      common /trim/ dmtrmn,dmtrmnm,dlcln,dlclnm,trtol,cmy,cnw,alf0,
     .              alf1,dzdt,thtd0,thtd1,zrg0,zrg1,dtrmsmx,dtrmsmn,
     .              dalfmx,ddtmx,ddtrm0,ddtrm1,itrmt,itrminc,fp(4,4),
     .              tp(4,4),zlfct,epstr,relax,ittrst
      common /lesinfo/ les_model,les_wallscale,cs_smagorinsky,cs_wale,
     .                 cs_vreman
      common /sourceterm/ xdir_only_source
      common /random_input/ randomize
      common /mms/ iexact_trunc,iexact_disc,iexact_ring
      common /sminn/ isminc,ismincforce
      common /initfac/ scal_ic
      common /plot3dtyp/ ifunct
      common /memry/ lowmem_ux
      common /constit/ i_nonlin,c_nonlin,snonlin_lim,i_tauijs,i_qcr2000,
     .                 i_qcr2013,i_qcr2013v
      common /easmlim/ cmulim
      common /easmv/ c10,c11,c2,c3,c4,c5,sigk1,cmuc1,ieasm_type
      common /is_patch/ is_pat(5),ie_pat(5),ipatch1st
      common /reystressmodel/ issglrrw2012,i_sas_rsm,i_yapterm
      common /writestuff/ ifort50write, j_ifort50write, i_ifort50write
      common /sa_options/ i_saneg,i_sanoft2,sa_cw2,sa_cw3,sa_cv1,sa_ct3,
     .                    sa_ct4,sa_cb1,sa_cb2,sa_sigma,sa_karman
      common /lam/ ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi,
     .        i_lam_forcezero
      common /specialtop_kmax1001/ i_specialtop_kmax1001,
     .        a_specialtop_kmax1001,xc_specialtop_kmax1001,
     .        sig_specialtop_kmax1001,vtp_specialtop_kmax1001,
     .        wc_specialtop_kmax1001,fac_specialtop_kmax1001,
     .        cc_specialtop_kmax1001,xerf_specialtop_kmax1001,
     .        sigerf_specialtop_kmax1001
      common /iupdate/ iupdatemean
      common /wallvel_2034/ uub_2034,vvb_2034,wwb_2034,vvb_xaxisrot2034
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
c     angle-of-attack relaxation factor for specified Cl
c
      rlxalph  = 1.
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
      cflturb(1)  = 0.
      cflturb(2)  = 0.
      cflturb(3)  = 0.
      cflturb(4)  = 0.
      cflturb(5)  = 0.
      cflturb(6)  = 0.
      cflturb(7)  = 0.
      nkey     = nkey + 7
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
c     flag to turn on/off time marching flow solver.  When 0 the
c     flow solver is turned on.  When 1, only the mesh deformation
c     and grid metric subroutines are turned on.
c
      meshdef   = 0
      nkey      = nkey + 1
c
c     eddy viscosity limiter for two eqn. turbulence models: limit
c     eddy viscosity to edvislim times the laminar viscosity
c
      edvislim  = 1.e10
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
      dalim     = 0.2
      nkey      = nkey + 1
c
c     no. of cycles between alpha updates (if set > 0, this
c     method takes precedence over resupdt)
c
      icycupdt  = 1
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
c     flag to perform DES with turbulence model
c     0 = no DES, 1 = standard DES, 2 = DDES (TCFD 20:181-195 2006), 3 = DDES + prod
c      term mod based on fd
c
      ides = 0
      nkey     = nkey + 1
c
c     constant associated with DES or DDES
c     cdes=0.65 is standard value for SA-DES; SST-DES may need a different value
c     (AIAA 2001-0879 uses 0.78)
c
      cdes = 0.65
      nkey     = nkey + 1
c
c     additional constant associated with Modified DDES (ides .ge. 3)
c
      cddes = 0.975
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
      tur10(1) = -1.
      tur10(2) = -1.
      tur10(3) = -1.
      tur10(4) = -1.
      tur10(5) = -1.
      tur10(6) = -1.
      tur10(7) = -1.
      nkey     = nkey + 7
c
c     value that epsilon or omega is set to for 2-eqn models if 
c     it goes below tur1cutlev
      if (ivmx .eq. 15) then
        tur1cut = -1.
      else
        tur1cut = 1.e-20
      end if
      nkey     = nkey + 1
c
c     roll angle
      roll_angle = 0.
      nkey       = nkey + 1
c     k-omega/sst/k-epsilon parameter controlling production term
      ikoprod = 0
      nkey       = nkey + 1
c
c     sst parameter controlling denominator of mut term
      isstdenom = 0
      nkey       = nkey + 1
c
c     limiter term on Pk in the two-eqn models
      pklimterm = 20.
      nkey       = nkey + 1
c
c     parameter for use with particular axi cases (for which
c     i2d=0 & idim=2); modifies time step based on CFL number
c     so it does not depend on the i-direction metrics
      iaxi2plane = 0
      nkey       = nkey + 1
c
c     parameter controls whether turbulence model advection terms
c     are 1st or 2nd order upwind on RHS (1=1st, 2=2nd)
c     (note: LHS uses 1st order in both cases)
      iturbord = 1
      nkey       = nkey + 1
c
c     value that k is set to for 2-eqn models if 
c     it goes below tur2cutlev
      tur2cut = 1.e-20
      nkey     = nkey + 1
c
c     lower limit on epsilon or omega for 2-eqn models that,
c     when reached, causes it to be reset to tur1cut
      tur1cutlev = 0.
      nkey     = nkey + 1
c
c     lower limit on k for 2-eqn models that,
c     when reached, causes it to be reset to tur2cut
      tur2cutlev = 0.
      nkey     = nkey + 1
c
c     full Navier-Stokes
      ifullns = 0
      nkey     = nkey + 1
c
c     beta8 parameter in k-enstrophy model
      ibeta8kzeta = 0
      nkey     = nkey + 1
c
c     SARC - Spalart-Allmaras curvature correction in 2D sense only
c     does not account for system rotation
      isarc2d = 0
      nkey     = nkey + 1
c
c     SARC - Spalart-Allmaras curvature correction in 3D
c     does not account for system rotation
      isarc3d = 0
      nkey     = nkey + 1
c
c     parameter cr3 in SARC model
      sarccr3 = 1.0
      nkey     = nkey + 1
c
c     EASMCC - EASM curvature correction in 2D sense only
      ieasmcc2d = 0
      nkey     = nkey + 1
c
c     ipertavg - Average r,u,v,w,p and second moments
      ipertavg = 0 
      nkey     = nkey + 1 

c     icoarsemovie - write movie files of coarse planar and volume data
      icoarsemovie = 0
      nkey     = nkey + 1
c
c     i2dmovie - write 2D movie files of coarse planar and volume data
      i2dmovie = 0
      nkey     = nkey + 1
c
c     iclcd - write cl and cd history for separate blocks using clcd.inp
      iclcd    = 0
      nkey     = nkey + 1
c
c     iskip_blocks - skip factor for 2d coarse movie
      iskip_blocks    = 1
      nkey     = nkey + 1
c
c     cfltauMax - Maximum value for the cfltau during subiterations
      cfltauMax= -1.0
      nkey     = nkey + 1
c
c     cfltau0 - Exponent for cfltau growth with subiteration
      cfltau0= 1.0
      nkey     = nkey + 1

c     parameter to perform trim (1) or not to perform trim (0)
      irbtrim = 0
      nkey     = nkey + 1
c
c     parameter to perform rigid body sim (1) or not to perform (0)
      irigb    = 0
      nkey     = nkey + 1
c
c     Conversion factor from rigid body dynamics length scale
c     to CFD grid length scale.
      greflrb  = 1.
      nkey     = nkey + 1
c
c     Total aircraft mass for trim and rigid body motion
      tmass = 1.
      nkey     = nkey + 1
c
c     Total aircraft y-mass moment of inertia for
c     rigid body motion
      yinert = 1.
      nkey     = nkey + 1
c
c     Acceleration of gravity for rigid/flex trim and motion
      gaccel = 1.
      nkey     = nkey + 1
c
c     Relaxation parameter for the trim algorithm
      relax = 0.5
      nkey     = nkey + 1
c
c     Iteration increment to update trim equations
      itrminc = 5
      nkey     = nkey + 1
c
c     dclda to be used for trim
      dclda    = 6.
      nkey     = nkey + 1
c
c     dcldd to be used for trim
      dcldd    = 1.4
      nkey     = nkey + 1
c
c     dcmda to be used for trim
      dcmda    =-0.2
      nkey     = nkey + 1
c
c     dcmdd to be used for trim
      dcmdd    =-0.88
      nkey     = nkey + 1
c
c     Parameter controling dynamic grid input
      ndgrd    = 0
      nkey     = nkey + 1
c
c     Parameter controling dynamic grid output
      ndwrt    = 0
      nkey     = nkey + 1
c
c     i_bsl used to turn on Menter's BSL model
      i_bsl     = 0
      nkey     = nkey + 1
c
c     keepambient used to maintain ambient turb levels for 2-eq models
      keepambient = 0
      nkey     = nkey + 1
c
c     re_thetat0 used for transition model
c     < 0  = use the default computed value
c     >= 0 = use this number as the specified user input value
c     (default=-1)
      re_thetat0 = -1.0
      nkey     = nkey + 1
c
c     turbintensity_inf_percent is freestrem turb intensity
c     setting this will override any value set for tur20=tur10(2)
c     this has no effect on tur10=tur10(1)
c     default=-1.
      turbintensity_inf_percent = -1.0
      nkey     = nkey + 1
c
c     eddy_visc_inf is freestream mu_t/mu_inf
c     setting this will override any value set for tur10=tur10(1)
c     this has no effect on tur20=tur10(2)
c     default=-1.
      eddy_visc_inf = -1.0
      nkey     = nkey + 1
c
c     cs_smagorinsky is the (non-dynamic) Smagorinsky constant for LES 
c     (ivisc=25 + les_model=1)
c     if set to zero, you get implicit LES (no model)
c     typical values are 0.1-0.2
c     default=0.0
      cs_smagorinsky = 0.0
      nkey     = nkey + 1
c
c     xdir_only_source is imposed source term
c     default=0.0
      xdir_only_source = 0.0
      nkey     = nkey + 1
c
c     randomize adds random perturbation to restart
c     value is max multiple of current value of q that can be added
c     default=0.0
      randomize = 0.0
      nkey     = nkey + 1
c
c     iexact_trunc is used to check truncation error against exact
c     solution (MMS=method of manufactured solution)
c     should be run only 1 iteration
c     set=1 for MS1, 2 for MS2, 4 for MS4
c     default=0
      iexact_trunc = 0
      nkey     = nkey + 1
c
c     iexact_disc is used to check discretization error against exact
c     solution (MMS=method of manufactured solution)
c     must be run to convergence
c     set=1 for MS1, 2 for MS2, 4 for MS4
c     default=0
      iexact_disc = 0
      nkey     = nkey + 1
c
c     iexact_ring is used to overwrite the exact solution on the
c     outer 2 ring layers of the grid (MMS=method of manufactured solution)
c     default=0
      iexact_ring = 0
      nkey     = nkey + 1
c
c     i_wilcox06 is used to change Wilcox88 model (ivisc=6) to Wilcox06
c     default=0
      i_wilcox06 = 0
      nkey     = nkey + 1
c
c     i_wilcox06_chiw is used to turn on/off Wilcox06 k-omega vortex
c     stretching parameter - only has effect if ivmx=6 and i_wilcox06=1
c     default=1 (ON)
      i_wilcox06_chiw = 1
      nkey     = nkey + 1
c
c     i_turbprod_kterm is used to determine whether 2/3*rho*k term
c     gets subtracted from turbulence production in ivisc=6,7 models
c     (only does anything if ivisc=6 or 7, ikoprod=1, and 
c     i_turbprod_kterm=1)
c     default=0 (term not included)
      i_turbprod_kterm = 0
      nkey     = nkey + 1
c
c     i_catris_kw set to 1 alters the k-omega turb diffusion terms to
c     include density (compressibility) effects, as given in Catris et al,
c     Aerosp Sci Technol. 4 (2000) 1-11
c     default=0
      i_catris_kw = 0
      nkey     = nkey + 1
c
c     ismincforce overrides normal smin/initvist restart usage -
c     chief use for this would be if restarting a solution on an altered
c     grid of the same dimensions (unless you force smin to be recomputed,
c     the new computation will have the WRONG min distances on the new
c     grid!). Note that if this keyword is used, the user must also 
c     appropriately specify whether initvist should be called or not
c     (it should be called if turb model is being changed upon restart)
c     ismincforce=0 do not compute smin, do not call initvist
c                =1 compute smin       , call initvist
c                =2 compute smin       , do not call initvist
c                =3 do not compute smin, call initvist
c     default=-1 (do not override normal usage)
      ismincforce = -1
      nkey     = nkey + 1
c
c     prod2d3dtrace forces Sij used in 2SijSij to be traceless; in 
c     production term in ivmx=6,7,10,30,40 when ikoprod=1; and in 
c     Wilcox06 stress-limiter term
c                =0.5 in 2-D sense
c                =0.333333 in 3-D sense
c     default=0. (incompressible Sij used in 2SijSij term)
      prod2d3dtrace = 0.
      nkey     = nkey + 1
c
c     i_compress_correct adds dilatation-dissipation type compressibility
c     correction, currently to ivisc=6 or 7 only
c         =1 Wilcox-type (Turbulence Modeling for CFD, ed 3, p. 258)
c         =2 Zeman-type for boundary layers (AIAA 93-0897)
c         =0 (default) or any other number beside 1 or 2 - no correction
      i_compress_correct = 0
      nkey     = nkey + 1
c
c     les_model determines LES subgrid model to use with ivisc=25
c         =0 no model (default)
c         =1 standard Smagorinsky model
c         =2 WALE model (Flow, Turb, & Combust 62:183-200 1999)
c         =3 Vreman model (Phys Fluids 16(10):3670-3681 2004)
      les_model = 0
      nkey     = nkey + 1
c
c     les_wallscale turns on van Driest type wall scaling of Delta in LES model
c     only used in conjunction with Smagorinsky model les_model=1
c         =0 scaling off (default)
c         =1 scaling on
      les_wallscale = 0
      nkey     = nkey + 1
c
c     cs_wale is the (non-dynamic) WALE constant for LES 
c     (ivisc=25 + les_model=2)
c     if set to zero, you get implicit LES (no model)
c     typical values are 0.45-0.6
c     default=0.0
      cs_wale = 0.0
      nkey     = nkey + 1
c
c     cs_vreman is the (non-dynamic) Vreman constant for LES 
c     (ivisc=25 + les_model=3)
c     if set to zero, you get implicit LES (no model)
c     typical values are 0.025-0.1
c     default=0.0
      cs_vreman = 0.0
      nkey     = nkey + 1
c
c     isstrc =1 turns on curvature correction for ivmx=6 or 7 (AIAA 98-2554)
c     isstrc =2 turns on CC for ivmx=6 or 7 (Smirnov & Menter, ASME
c               JoTurbomachinery, Vol 131, 2009, 041010)
c     default=0
      isstrc = 0
      nkey     = nkey + 1
c
c     sstrc_crc is the constant for isstrc=1
c     default=1.4
      sstrc_crc = 1.4
      nkey     = nkey + 1
c
c     isstsf turns on separation fix (tends to reattach sooner) for ivmx=6 or 7
c     default=0
      isstsf = 0
      nkey     = nkey + 1
c
c     scal_ic is a scaling factor for 2-eqn model BL-type approx ICs
c     default=5.e6 (smaller no. makes profile thinner, 0 uses freestream ICs)
      scal_ic = 5.e6
      nkey     = nkey + 1
c
c     ifunct will output a PLOT3D function file (instead of Q-type file)
c     only when iptype=2 in the plot3d section of the input file.
c     There will be ifunct variables in the file.
c     (Use of this requires hardwire mod in plot3t.F, to specify what
c     variables are to be output.)
c     Note that if ifunct is large, you may run out of memory when writing
c     and need to augment with keyword memadd.
      ifunct = 0
      nkey     = nkey + 1
c
c     lowmem_ux=1 allows reversion to lower memory usage of ux array
c     (i.e., only allocate this jdim*kdim*idim*9 array for velocity 
c     derivatives if it is needed for a particular model)
c     By default, after V6.5, the memory is always allocated and ux
c     is always computed 
      lowmem_ux=0
      nkey     = nkey + 1
c
c     SAR - Spalart-Allmaras rotation correction (Dacles-Mariani,
c     AIAA J 33(9):1561-1568, 1995
      isar = 0
      nkey     = nkey + 1
c
c     constant for SAR - Spalart-Allmaras rotation correction
      crot = 2.0
      nkey     = nkey + 1
c
c     nonlinear constitutive relation for use with linear models
      i_nonlin = 0
      nkey     = nkey + 1
c
c     nonlinear constitutive relation constant
      c_nonlin = 0.3
      nkey     = nkey + 1
c
c     limiter for nonlinear constitutive relation constant
c     active in ffluxv, hfluxv, and gfluxv; acts as lower limit on denom term
c     default=1.e-10
      snonlin_lim = 1.e-10
      nkey     = nkey + 1
c
c     isubit_r - write out the subiteration residual for all variables to subit_res
      isubit_r    = 0
      nkey     = nkey + 1
c
c     i_wilcox98 is used to change Wilcox88 model (ivisc=6) to Wilcox06
c     default=0
      i_wilcox98 = 0
      nkey     = nkey + 1
c
c     i_wilcox98_chiw is used to turn on/off Wilcox98 k-omega vortex
c     stretching parameter - only has effect if ivmx=6 and i_wilcox98=1
c     default=1 (ON)
      i_wilcox98_chiw = 1
      nkey     = nkey + 1
c
c     cmulim limits the abs min computed value of cmu for EASM (ivisc=8,9,13,14);
c     may need higher number (.025?) for supersonic BL flow to avoid kink at
c     BL edge on very fine grids;
c     should never be set higher than .04 or so (.08-.09 is the std value of cmu
c     in the log layer)
c     default=.0005
      cmulim = .0005
      nkey     = nkey + 1
c
c     iaxi2planeturb used to force solve of field turb model eqns in j-k 
c     plane only; only allowed if i2d=0 & idim=2
c     default=0 (include i-dir)
      iaxi2planeturb=0
      nkey     = nkey + 1
c
c     istrongturbdis for strong conservation of diss terms in field turb eqns
c     default=0 (weak conservation)
      istrongturbdis=0
      nkey     = nkey + 1
c
c     ieasm_type used to change pressure-strain models within context of
c     existing k-eps or k-omega (EXPERIMENTAL ONLY - not documented officially)
c     default=0 (linear SSG, as used by std Gatski-Rumsey)
c     1=Wallin-Johansson type
c     2=LRR-QI type
      ieasm_type=0
      nkey     = nkey + 1
c
c     ipatch1st=1 forces patching interpolation to be 1st order (default=0)
      ipatch1st= 0
      nkey     = nkey + 1
c
c     isst2003=1 when ivisc=7 forces SST-2003 model, Menter, F. R., Kuntz, M., 
c     and Langtry, R., "Ten Years of Industrial Experience with the SST 
c     Turbulence Model," Turbulence, Heat and Mass Transfer 4, 
c     ed: K. Hanjalic, Y. Nagano, and M. Tummers, Begell House, Inc., 2003, 
c     pp. 625 - 632. 
      isst2003=0
      nkey     = nkey + 1
c
c     issglrrw2012=1 when using ivisc=72 and want SSG/LRR-RSM-w2012 instead
c     of WilcoxRSM-w2006; 
c     issglrrw2012=2 if want SSG/LRR-RSM-w2012 with F1 (blend) forced 
c     to be 1.0, resulting in essentially WilcoxRSM-w1988 (with generalized
c     gradient diffusion)
c     issglrrw2012=3 if want SSG/LRR-RSM-w2012 with simple diffusion
c     issglrrw2012=4 if want SSG/LRR-RSM-w2012 with F1=1 and simple diffusion
c     issglrrw2012=5 if want SSG/LRR-RSM-w2012 with Wilcox simple diffusion
c     (non-blended)
c     issglrrw2012=6 same as 1, except uses g-eqn (=1/sqrt(omega))
c     instead of omega
      issglrrw2012=0
      nkey     = nkey + 1
c
c   optional printout to fort.50 (specifically for k-line; k=1 must be body)
c   writes "plus" turbulent values; for use with single-block cases only
c   iptype must=2 in the PLOT3D output section of the input file
      ifort50write=0
      nkey     = nkey + 1
      j_ifort50write=1
      nkey     = nkey + 1
      i_ifort50write=1
      nkey     = nkey + 1
c
c   i_sas_rsm adds/subtracts SAS-like term to RSM omega eqn (AIAA-2014-0586)
c   Set to 1 to make it MORE like SAS
c   Set to -1 to keep the RANS character, adding more eddy visc in shear layers
      i_sas_rsm=0
      nkey     = nkey + 1
c
c   iforcev0=1 forces v=0 on update
      iforcev0=0
      nkey     = nkey + 1
c
c   i_saneg=1 uses SA-neg version
      i_saneg=0
      nkey     = nkey + 1
c
c   i_sanoft2=1 uses SA-noft2 version
      i_sanoft2=0
      nkey     = nkey + 1
c
c   i_lam_forcezero=1 forces laminar regions to have zero vist3d (i.e.,
c   will not convect it)
      i_lam_forcezero=0
      nkey     = nkey + 1
c
c   i_specialtop_kmax1001=1 only for particular special BC on 1001
c   specifying kmax wall velocity in conjunction with BC1001:
c   v = a*{fac*(xc-x)+(1-fac)}*exp[cc-({(xc-x)/sig}**2)] + 
c       vtp*{0.5*(1-erf[(x-xerf)/sigerf])}
c   and wc is for specifying crossflow
      i_specialtop_kmax1001=0
      a_specialtop_kmax1001=0.
      xc_specialtop_kmax1001=0.
      sig_specialtop_kmax1001=1.
      vtp_specialtop_kmax1001=0.
      wc_specialtop_kmax1001=0.
      fac_specialtop_kmax1001=1.0
      cc_specialtop_kmax1001=0.5
      xerf_specialtop_kmax1001=100000.
      sigerf_specialtop_kmax1001=1.e-9
      nkey     = nkey + 10
c
c   SA constants
      sa_cw2=0.3
      sa_cw3=2.0
      sa_cv1=7.1
      sa_ct3=1.2
      sa_ct4=0.5
      sa_cb1=0.1355
      sa_cb2=0.622
      sa_sigma=2./3.
      sa_karman=0.41
      nkey=nkey + 9
c
c   iupdatemean = 0 forces no mean flow update (turbulence still updated)
      iupdatemean=1
      nkey     = nkey + 1
c
c   uub_2034,vvb_2034,wwb_2034 are wall velocities used with bc2034
c   vvb_xaxisrot2034 is rotational velocity about x-axis with bc2034
      uub_2034=0.
      vvb_2034=0.
      wwb_2034=0.
      vvb_xaxisrot2034=0.
      nkey     = nkey + 4
c
c   i_tauijs adds turbulence via tauij terms (when set to 1)
c   NOTE: like regular method, 2/3*rho*k term is currently ignored
c   Also: not done accurately at boundaries for multiblock cases:
c   currently one sided; need ux derivs from across boundaries
c   to be truly correct
      i_tauijs=0
      nkey     = nkey + 1
c
c   i_qcr2000 turns on QCR2000 (when set to 1)
c   only works when i_tauijs=1
      i_qcr2000=0
      nkey     = nkey + 1
c
c   i_qcr2013 turns on QCR2013 (limited) (when set to 1)
c   only works when i_tauijs=1
      i_qcr2013=0
      nkey     = nkey + 1
c
c   i_qcr2013v turns on QCR2013-V (when set to 1)
c   only works when i_tauijs=1
      i_qcr2013v=0
      nkey     = nkey + 1
c
c   i_yapterm=1 adds Yap-like term to get rid of backbending in RSM
      i_yapterm=0
      nkey     = nkey + 1
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

c        do any order-dependent things here before exiting

#ifdef CMPLX
         xrotrate = cmplx(real(xrotrate),real(xrotrate_img))
         yrotrate = cmplx(real(yrotrate),real(yrotrate_img))
         zrotrate = cmplx(real(zrotrate),real(zrotrate_img))
#endif

c     Check for the existence of the averaging files
         if( iteravg .eq. 2 ) then
            isthere = .false.
            if (ipertavg .eq. 2 ) then
               filename = output_dir(1:len_trim(output_dir)) //
     $              "/" // "cfl3d_avg_ruvwp.p3d"
            else
               filename = output_dir(1:len_trim(output_dir)) //
     $              "/" // "cfl3d_avgq.p3d"
            end if
            inquire (file=filename,exist=isthere)
            if( isthere ) then
               open(150, file=filename, form='unformatted',
     $              status='old')
               read(150,end=1011,err=1011) idummy
               close(150)
               goto 1013
            else
               goto 1012
            end if
         else
            goto 1013
         end if

 1011    close(150)
 1012    iteravg = 1
         if( ipertavg .eq. 2 ) ipertavg = 1
 1013    continue
         
c     Check for the existence of the clcd.bin file
         if( iclcd .eq. 2 ) then
            isthere = .false.
            filename = output_dir(1:len_trim(output_dir)) //
     $           "/" // "clcd.bin"
            inquire (file=filename,exist=isthere)
            if( isthere ) then
               open(unit=150,file=filename,form='unformatted',
     $              status='old')
               read(150,end=1014,err=1014) idummy1, idummy2, idummy3
               close(150)
               goto 1016
            else
               goto 1015
            end if
         else
            goto 1016
         end if

 1014    close(150)
 1015    iclcd = 1
 1016    continue
         
c        the variables turbintensity_inf_percent and eddy_visc_inf
c        do not go to the rest of the code; instead, if set, they
c        modify tur10(2) and tur10(1), respectively
c        (or tur10(1,2,3) and tur10(7) for ivmx=72)
         if (real(turbintensity_inf_percent) .ge. 0.) then
           if (ivmx .eq. 72) then
             if (iunit11.gt.0) write(iunit11,'('' turbintensity_inf'',
     +         ''_percent always overrides tur10,20,30 for ivmx=72'')')
             zkinf=1.5*xmach**2*(turbintensity_inf_percent/100.)**2
             tur10(1)=-2./3.*zkinf
             tur10(2)=-2./3.*zkinf
             tur10(3)=-2./3.*zkinf
           else
             if (iunit11.gt.0) write(iunit11,'('' turbintensity_inf'',
     +         ''_percent always overrides tur20'')')
             tur10(2)=1.5*xmach**2*(turbintensity_inf_percent/100.)**2
           end if
         end if
         if (real(eddy_visc_inf) .ge. 0.) then
           if (ivmx .eq. 72) then
            if (iunit11.gt.0) write(iunit11,'('' eddy_visc_inf always'',
     +        '' overrides tur70 for ivmx=72'')')
           else
            if (iunit11.gt.0) write(iunit11,'('' eddy_visc_inf always'',
     +        '' overrides tur10'')')
           end if
           if (ivmx .eq. 4) then
             tur10(1)=eddy_visc_inf/0.09
           else if (ivmx.eq.6 .or. ivmx.eq.7 .or. ivmx.eq.30 .or.
     +              ivmx.eq.40) then
             tur10(1)=tur10(2)/eddy_visc_inf
           else if (ivmx.eq.8 .or. ivmx.eq.14) then
             tur10(1)=0.0895*tur10(2)/eddy_visc_inf
           else if (ivmx.eq.12) then
             tur10(1)=0.081*tur10(2)/eddy_visc_inf
           else if (ivmx.eq.10 .or. ivmx.eq.15) then
             tur10(1)=0.09*tur10(2)**2/eddy_visc_inf
           else if (ivmx.eq.9 .or. ivmx.eq.13) then
             tur10(1)=0.0885*tur10(2)**2/eddy_visc_inf
           else if (ivmx.eq.11) then
             tur10(1)=0.081*tur10(2)**2/eddy_visc_inf
           else if (ivmx .eq. 5) then
c            for SA model, need to iterate to find nuwiggle:
             tur10(1)=eddy_visc_inf
             do n=1,100
               errornum=tur10(1)**4-eddy_visc_inf*tur10(1)**3-
     +           357.911*eddy_visc_inf
               errorden=4.*tur10(1)**3-3.*eddy_visc_inf*tur10(1)**2
               error=errornum/errorden
               if (abs(error) .lt. 1.e-6) goto 105
               tur10(1)=tur10(1)-error
             enddo
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''stopping...Newton iteration'',
     +         '' using eddy_visc_inf not converged in readkey'')')
             call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
 105         continue
           else if (ivmx .eq. 72) then
             tur10(7)=-3./2.*tur10(1)/eddy_visc_inf
           end if
         end if
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
      else if (inpstr(lc1:lc2).eq.'rlxalph') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         rlxalph = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'nsubturb') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) nsubturb
c
      else if (inpstr(lc1:lc2).eq.'cflturb') then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''stopping...keyword cflturb no '',
     .       ''longer allowed... use cflturb1, 2, etc instead'')')
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
c
      else if (inpstr(lc1:lc2).eq.'cflturb1') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cflturb(1) = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'cflturb2') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cflturb(2) = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'cflturb3') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cflturb(3) = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'cflturb4') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cflturb(4) = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'cflturb5') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cflturb(5) = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'cflturb6') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cflturb(6) = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'cflturb7') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cflturb(7) = abs(realval(1))
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
      else if (inpstr(lc1:lc2).eq.'meshdef') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) meshdef
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
c
      else if (inpstr(lc1:lc2).eq.'yrotrate_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         yrotrate_img = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'zrotrate_img') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         zrotrate_img = realval(1)
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
      else if (inpstr(lc1:lc2).eq.'cddes') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cddes = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'iteravg') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iteravg
c
      else if (inpstr(lc1:lc2).eq.'tur10') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur10(1) = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tur20') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur10(2) = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tur30') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur10(3) = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tur40') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur10(4) = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tur50') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur10(5) = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tur60') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur10(6) = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tur70') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur10(7) = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tur1cut') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur1cut = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'roll_angle') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         roll_angle = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'ikoprod') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ikoprod
c
      else if (inpstr(lc1:lc2).eq.'isstdenom') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) isstdenom
c
      else if (inpstr(lc1:lc2).eq.'pklimterm') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         pklimterm = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'iaxi2plane') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iaxi2plane
c
      else if (inpstr(lc1:lc2).eq.'iturbord') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iturbord
c
      else if (inpstr(lc1:lc2).eq.'tur2cut') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur2cut = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tur1cutlev') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur1cutlev = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tur2cutlev') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tur2cutlev = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'ifullns') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ifullns
c
      else if (inpstr(lc1:lc2).eq.'ibeta8kzeta') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ibeta8kzeta
c
      else if (inpstr(lc1:lc2).eq.'isarc2d') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) isarc2d
c
      else if (inpstr(lc1:lc2).eq.'isarc3d') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) isarc3d
c
      else if (inpstr(lc1:lc2).eq.'sarccr3') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sarccr3 = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'ieasmcc2d') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ieasmcc2d
c
      else if (inpstr(lc1:lc2).eq.'ipertavg') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ipertavg
         if (ipertavg .ne. 0) iteravg = ipertavg
c
      else if (inpstr(lc1:lc2).eq.'icoarsemovie') then 
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) icoarsemovie
c
      else if (inpstr(lc1:lc2).eq.'i2dmovie') then 
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i2dmovie
c
      else if (inpstr(lc1:lc2).eq.'iclcd') then 
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iclcd
c
      else if (inpstr(lc1:lc2).eq.'iskip_blocks') then 
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iskip_blocks
c
      else if (inpstr(lc1:lc2).eq.'cfltauMax') then 
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) cfltauMax
c
      else if (inpstr(lc1:lc2).eq.'cfltau0') then 
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) cfltau0
c
      else if (inpstr(lc1:lc2).eq.'irbtrim') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) irbtrim
c
      else if (inpstr(lc1:lc2).eq.'irigb') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) irigb
c
      else if (inpstr(lc1:lc2).eq.'greflrb') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         greflrb = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'tmass') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         tmass = realval(1) 
c
      else if (inpstr(lc1:lc2).eq.'yinert') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         yinert = realval(1) 
c
      else if (inpstr(lc1:lc2).eq.'gaccel') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         gaccel = realval(1) 
c
      else if (inpstr(lc1:lc2).eq.'relax') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         relax = realval(1) 
c
      else if (inpstr(lc1:lc2).eq.'itrminc') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) itrminc
c
      else if (inpstr(lc1:lc2).eq.'dclda') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         dclda = realval(1) 
c
      else if (inpstr(lc1:lc2).eq.'dcldd') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         dcldd = realval(1) 
c
      else if (inpstr(lc1:lc2).eq.'dcmda') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         dcmda = realval(1) 
c
      else if (inpstr(lc1:lc2).eq.'dcmdd') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         dcmdd = realval(1) 
c
      else if (inpstr(lc1:lc2).eq.'ndgrd') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ndgrd
c
      else if (inpstr(lc1:lc2).eq.'ndwrt') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ndwrt
c
      else if (inpstr(lc1:lc2).eq.'i_bsl') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_bsl
c
      else if (inpstr(lc1:lc2).eq.'keepambient') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) keepambient
c
      else if (inpstr(lc1:lc2).eq.'re_thetat0') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         re_thetat0 = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'turbintensity_inf_percent') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         turbintensity_inf_percent = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'eddy_visc_inf') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         eddy_visc_inf = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'cs_smagorinsky') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cs_smagorinsky = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'xdir_only_source') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         xdir_only_source = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'randomize') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         randomize = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'iexact_trunc') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iexact_trunc
c
      else if (inpstr(lc1:lc2).eq.'iexact_disc') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iexact_disc
c
      else if (inpstr(lc1:lc2).eq.'iexact_ring') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iexact_ring
c
      else if (inpstr(lc1:lc2).eq.'i_wilcox06') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_wilcox06
c
      else if (inpstr(lc1:lc2).eq.'i_wilcox06_chiw') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_wilcox06_chiw
c
      else if (inpstr(lc1:lc2).eq.'i_turbprod_kterm') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_turbprod_kterm
c
      else if (inpstr(lc1:lc2).eq.'i_catris_kw') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_catris_kw
c
      else if (inpstr(lc1:lc2).eq.'ismincforce') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ismincforce
c
      else if (inpstr(lc1:lc2).eq.'prod2d3dtrace') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         prod2d3dtrace = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'i_compress_correct') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_compress_correct
c
      else if (inpstr(lc1:lc2).eq.'les_model') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) les_model
c
      else if (inpstr(lc1:lc2).eq.'les_wallscale') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) les_wallscale
c
      else if (inpstr(lc1:lc2).eq.'cs_wale') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cs_wale = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'cs_vreman') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cs_vreman = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'isstrc') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) isstrc
c
      else if (inpstr(lc1:lc2).eq.'sstrc_crc') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sstrc_crc = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'isstsf') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) isstsf
c
      else if (inpstr(lc1:lc2).eq.'scal_ic') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         scal_ic = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'ifunct') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ifunct
c
      else if (inpstr(lc1:lc2).eq.'lowmem_ux') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) lowmem_ux
c
      else if (inpstr(lc1:lc2).eq.'isar') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) isar
c
      else if (inpstr(lc1:lc2).eq.'crot') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         crot = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'i_nonlin') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_nonlin
c
      else if (inpstr(lc1:lc2).eq.'c_nonlin') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         c_nonlin = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'snonlin_lim') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         snonlin_lim = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'isubit_r') then 
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) isubit_r
c
      else if (inpstr(lc1:lc2).eq.'i_wilcox98') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_wilcox98
c
      else if (inpstr(lc1:lc2).eq.'i_wilcox98_chiw') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_wilcox98_chiw
c
      else if (inpstr(lc1:lc2).eq.'cmulim') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cmulim = abs(realval(1))
c
      else if (inpstr(lc1:lc2).eq.'iaxi2planeturb') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iaxi2planeturb
c
      else if (inpstr(lc1:lc2).eq.'istrongturbdis') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) istrongturbdis
c
      else if (inpstr(lc1:lc2).eq.'ieasm_type') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ieasm_type
c
      else if (inpstr(lc1:lc2).eq.'ipatch1st') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ipatch1st
c
      else if (inpstr(lc1:lc2).eq.'isst2003') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) isst2003
c
      else if (inpstr(lc1:lc2).eq.'issglrrw2012') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) issglrrw2012
c
      else if (inpstr(lc1:lc2).eq.'ifort50write') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) ifort50write
c
      else if (inpstr(lc1:lc2).eq.'j_ifort50write') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) j_ifort50write
c
      else if (inpstr(lc1:lc2).eq.'i_ifort50write') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_ifort50write
c
      else if (inpstr(lc1:lc2).eq.'i_sas_rsm') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_sas_rsm
c
      else if (inpstr(lc1:lc2).eq.'iforcev0') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iforcev0
c
      else if (inpstr(lc1:lc2).eq.'i_saneg') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_saneg
c
      else if (inpstr(lc1:lc2).eq.'i_sanoft2') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_sanoft2
c
      else if (inpstr(lc1:lc2).eq.'i_lam_forcezero') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_lam_forcezero
c
      else if (inpstr(lc1:lc2).eq.'i_specialtop_kmax1001') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_specialtop_kmax1001
      else if (inpstr(lc1:lc2).eq.'a_specialtop_kmax1001') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         a_specialtop_kmax1001 = realval(1)
      else if (inpstr(lc1:lc2).eq.'xc_specialtop_kmax1001') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         xc_specialtop_kmax1001 = realval(1)
      else if (inpstr(lc1:lc2).eq.'sig_specialtop_kmax1001') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sig_specialtop_kmax1001 = realval(1)
      else if (inpstr(lc1:lc2).eq.'vtp_specialtop_kmax1001') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         vtp_specialtop_kmax1001 = realval(1)
      else if (inpstr(lc1:lc2).eq.'wc_specialtop_kmax1001') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         wc_specialtop_kmax1001 = realval(1)
      else if (inpstr(lc1:lc2).eq.'fac_specialtop_kmax1001') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         fac_specialtop_kmax1001 = realval(1)
      else if (inpstr(lc1:lc2).eq.'cc_specialtop_kmax1001') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         cc_specialtop_kmax1001 = realval(1)
      else if (inpstr(lc1:lc2).eq.'xerf_specialtop_kmax1001') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         xerf_specialtop_kmax1001 = realval(1)
      else if (inpstr(lc1:lc2).eq.'sigerf_specialtop_kmax1001') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sigerf_specialtop_kmax1001 = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'sa_cw2') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sa_cw2 = realval(1)
      else if (inpstr(lc1:lc2).eq.'sa_cw3') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sa_cw3 = realval(1)
      else if (inpstr(lc1:lc2).eq.'sa_cv1') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sa_cv1 = realval(1)
      else if (inpstr(lc1:lc2).eq.'sa_ct3') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sa_ct3 = realval(1)
      else if (inpstr(lc1:lc2).eq.'sa_ct4') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sa_ct4 = realval(1)
      else if (inpstr(lc1:lc2).eq.'sa_cb1') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sa_cb1 = realval(1)
      else if (inpstr(lc1:lc2).eq.'sa_cb2') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sa_cb2 = realval(1)
      else if (inpstr(lc1:lc2).eq.'sa_sigma') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sa_sigma = realval(1)
      else if (inpstr(lc1:lc2).eq.'sa_karman') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         sa_karman = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'iupdatemean') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) iupdatemean
c
      else if (inpstr(lc1:lc2).eq.'uub_2034') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         uub_2034 = realval(1)
      else if (inpstr(lc1:lc2).eq.'vvb_2034') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         vvb_2034 = realval(1)
      else if (inpstr(lc1:lc2).eq.'wwb_2034') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         wwb_2034 = realval(1)
      else if (inpstr(lc1:lc2).eq.'vvb_xaxisrot2034') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) realval(1)
         vvb_xaxisrot2034 = realval(1)
c
      else if (inpstr(lc1:lc2).eq.'i_tauijs') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_tauijs
c
      else if (inpstr(lc1:lc2).eq.'i_qcr2000') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_qcr2000
c
      else if (inpstr(lc1:lc2).eq.'i_qcr2013') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_qcr2013
c
      else if (inpstr(lc1:lc2).eq.'i_qcr2013v') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_qcr2013v
c
      else if (inpstr(lc1:lc2).eq.'i_yapterm') then
         lc2 = lc2 +1
         read(inpstr(lc2:lcl),*) i_yapterm
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
c
c     Set trim data
      tp(1,1) = dclda
      tp(1,2) = dcldd
      tp(2,1) = dcmda
      tp(2,2) = dcmdd
c
      go to 1100
c
 1900 continue

c     Make sure we are running ddes if cddes is set

      if( abs(cddes - 0.975d0) > 0.001d0 ) ides = 3 
      if( real(cddes) > 0.999d0) cddes = 0.999d0
      if( ides == 3 ) then
         write(*,*) 'Running ides with cddes = ', cddes
      end if

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
