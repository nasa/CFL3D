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
      subroutine bc(ntime,nbl,lw,lw2,w,mgwk,wk,nwork,cl,nou,bou,nbuf,
     .              ibufdim,maxbl,maxgr,maxseg,itrans,irotat,idefrm,
     .              igridg,nblg,nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .              nbckdim,ibcinfo,jbcinfo,kbcinfo,bcfilei,bcfilej,
     .              bcfilek,lwdat,myid,idimg,jdimg,kdimg,bcfiles,
     .              mxbcfil,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determine boundary data/conditions at edges of grids.
c***********************************************************************
c
c     maxbl  - maximum number of blocks
c     maxgr  - maximum number of grids
c     mxbli  - maximum number of 1-1 interfaces (including coarser blocks!)
c     maxseg - maximum number of segments on a boundary
c     intmax - maximum number of block interpolations
c     maxxe  - maximum dimension of array used to store generalized
c              coordinate data used in interpolations
c     nsub1  - maximum number of blocks a single patch face may be
c                interpolated from
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*80 filname
      character*80 bcfiles(mxbcfil)
      character*120 bou(ibufdim,nbuf)
c
      integer bcfilei,bcfilej,bcfilek
      dimension nou(nbuf)
      dimension w(mgwk),lw(65,maxbl),lw2(43,maxbl),wk(nwork)
c
      dimension lwdat(maxbl,maxseg,6)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension bcfilei(maxbl,maxseg,2),bcfilej(maxbl,maxseg,2),
     .          bcfilek(maxbl,maxseg,2),igridg(maxbl),itrans(maxbl),
     .          irotat(maxbl),idefrm(maxbl),idimg(maxbl),jdimg(maxbl),
     .          kdimg(maxbl)
      dimension nblg(maxgr)
c
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /sklton/ isklton
      common /chkbc/ isym,jsym,ksym,iwrap,jwrap,kwrap
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
c
c***********************************************************************
c     Standard Boundary Condition Types
c***********************************************************************
c
c      all standard boundary condition types are applicable to any segment of
c      any face of a given block:
c
c      I - Physical Boundary Conditions:
c
c      the 1000 series bc's have ndata=0...no auxiliary data required
c 
c      bc1000 - freestream
c      bc1001 - symmetry plane
c      bc1002 - extrapolation from inside the computational domain
c      bc1003 - inflow/outflow
c      bc1005 - inviscid surface (also 1006 - uses norm mom)
c      bc1007 - axisymmetric about x-axis
c      bc1008 - tunnel inflow - v,w set to zero, u set from constant
c               total enthalpy, rho set from entropy at reference
c               (infinity) condition; pressure is extrapolated
c               from inside the computational domain.
c      bc1011 - singular axis - half plane symmmetry
c      bc1012 - singular axis - full plane - flux specified
c      bc1013 - singular axis - extrapolation (partial plane)
c
c      the 2000 series bc's *require* auxiliary data, ndata>0
c
c      bc2002 - pressure ratio (p/pinf) specified; rho,u,v,w extrapolated
c               from inside the computational domain. set ndata=1, and
c               specify the desired pressure ratio 
c      bc2003 - engine inflow: set ndata =5 and specify (estimated) mach number,
c               total pressure ratio (Ptotal/pinf), total temperature ratio
c               (Ttotal/tinf), and flow directions (alpe and betae, in degrees)
c               at the inlet; these are then used as the external state in a 1D 
c               characteristic boundary condition.
c      bc2004 - viscous surface (set ndata=2)
c               Twtype is specified 
c                            if > 0 sets Tw/Tinf = Twtype,
c                            if = 0 sets adiabatic wall,
c                            if < 0 sets Tw at stagnation temperature;
c               also mass-flow coefficient c_q is specified 
c                            if < 0 suction if min. surface,
c                                   blowing if max. surface 
c                            if = 0 solid wall,
c                            if > 0 blowing if min. surface,
c                                   suction if max surface)
c
c                            (i.e. cq>0 gives mass flow in the direction
c                             of increasing computational coordinate,
c                             while cq<0 gives mass flow in the direction
c                             of decreasing computational coordinate)
c               THIS BC SUPERCEDES OLD BC1004 & BC2004
c      bc2014 - same as bc2004, except it indicates that the region should
c               be forced to be LAMINAR
c      bc2024 - same as bc2004, except it includes a parameter denoting the
c               boundary condition for intermittency
c      bc2034 - same as bc2004, except it includes imposed wall velocity
c               (must have keyword inputs uub_2034, vvb_2034, wwb_2034,
c               or vvb_xaxisrot2034 set)
c      bc2006 - set pressure to satisfy the radial equilibrium equation dp/dr =
c               rho*vtheta**2/r; extrapolate the other flow variables. in this
c               boundary condition, the pressure ratio p/pinf is set at either
c               the bottom or the top of the block face, and then the eqn. for
c               dp/dr is integrated in either increasing or decreasing radial
c               direction, so that the pressure is specified over the entire
c               block face. this boundary condition is useful in turbomachinery
c               applications where there is a significant swirling flow.
c               set ndata=4, and specify grid FROM which the integration of 
c               pressure is continued (ngridc; if integration is not continued 
c               from another block, set ngridc=0), the pressure at either the
c               min or max radius in the current block (p/pinf, used if 
c               ngridc = 0), the direction in which integration is to proceed
c               (int dir = +/- 1, 2, or 3 for +/- i, j, or k directions), and
c               the physical direction along which the radial axis lies
c               (axial dir= 1, 2, or 3 for x, y, or z)
c      bc2007 - set all primative variables, with standard cfl3d normalization:
c               rho/rho_inf, u/a_inf, v/a_inf, w/a_inf, p/(rho_inf*a_inf**2)
c               set ndata=5, and input the normalized primatives in the order
c               given above.
c      bc2008 - specify rho, u, v, and w; extrapolate p from the interior.
c      bc2016 - viscous surface (set ndata=7)
c           1) Twtype:
c              Twtype > 0 (fixed wall temperature Tw/Tinf = Twtype)
c              Twtype = 0 (adiabatic wall)
c              Twtype < 0 (fixed wall temperature = stagnation temp)
c
c           2) cq  = steady component of surface jet
c           3) cqu = unsteady component of surface jet
c
c              Set surface jet boundary conditions via mass-flow
c              coefficient:
c                                                    (rho * u)_jet
c              mass-flow coefficient (cq + cqu)  = ------------------
c                                                  (rho * u)_infinity
c
c           4) sjetx = direction number for surface jet (x-direction)
c           5) sjety = direction number for surface jet (y-direction)
c           6) sjetz = direction number for surface jet (z-direction)
c           7) rfreq = reduced frequency of unsteady jet =
c                      cqu*sin(twopi*rfreq*time)
c              NOTES:
c              1.  input direction numbers do not have to be normalized
c              2.  Set (sjetx,sjety,sjetz) = (0.0,0.0,0.0) for a
c                  wall-normal surface jet.  In this case,
c                  cq & cqu > 0 ==> flow into domain at i,j,k=1,1,1
c                  cq & cqu > 0 ==> flow out of domain at i,j,k=idim,jdim,kdim
c                  cq & cqu < 0 ==> flow out of domain at i,j,k=1,1,1
c                  cq & cqu < 0 ==> flow in to domain at i,j,k=idim,jdim,kdim
c              3.  rfreq = freq*lref/a_inf where freq is frequency in Hertz,
c                  lref is reference length, and a_inf is speed of sound
c              4.  For rfreq < 0: ramp the mass flow coefficient from zero to
c                  a terminal value where
c                  cq  = terminal mass flow coefficient value
c                  cqu = constant rate of change
c      bc2026 - sweeping jet (set ndata=9)
c           1) vmag = jet velocity magnitude/a_ref
c           2) rfreq = reduced frequency of unsteady jet sweeping =
c                      thetajet=sideangj*sin(twopi*rfreq*time)
c           3) sideangj = max angle (deg) that the jet "sweeps" in each direction
c                           (+-sideangj)
c           4) sxa = direction number for surface jet prior to sweeping (x-direction)
c           5) sya = direction number for surface jet prior to sweeping (y-direction)
c           6) sza = direction number for surface jet prior to sweeping (z-direction)
c           7) sxb = x-direction number for vector normal to (sxa,sya,sza), in
c                       plane perpendicular to body surface and pointing downstream
c           8) syb = y-direction number for vector normal to (sxa,sya,sza), in
c                       plane perpendicular to body surface and pointing downstream
c           9) szb = z-direction number for vector normal to (sxa,sya,sza), in
c                       plane perpendicular to body surface and pointing downstream
c              NOTES:
c              1.  input direction numbers do not have to be normalized
c              2.  rfreq = freq*lref/a_inf where freq is frequency in Hertz,
c                  lref is reference length, and a_inf is speed of sound
c      bc2018 - specify T, rhou, rhov, and rhow; extrapolate p from interior.
c      bc2028 - specify freq, rhoumax, rhovmax, and rhowmax; extrapolate rho 
c               and p from interior.
c      bc2038 - same as 2008, except add pseudo-random number to density and
c               velocity components, based on max +-5% of set local level, only 
c               for values for which u-vel component is less than 85% of xmach.
c      bc2009 - nozzle total BCs: set ndata =4 and specify
c               total pressure ratio (Ptotal/pinf), total temperature ratio
c               (Ttotal/tinf), and flow directions (alpe and betae, in degrees)
c               at the inlet; pressure is extrapolated from the interior.
c      bc2010 - nozzle total BCs-more consistent than 2009: set ndata =4 and specify
c               total pressure ratio (Ptotal/pinf), total temperature ratio
c               (Ttotal/tinf), and flow directions (alpe and betae, in degrees)
c               at the inlet; pressure is extrapolated from the interior.
c      bc2019 - nozzle total BCs: set ndata =2 and specify
c               Ptotal/Ptotal_inf and Ttotal/Ttotal_inf
c               (this BC is taken from OVERFLOW).
c      bc2102 - pressure ratio specified as a sinusoidal function of time; 
c               rho,u,v,w extrapolated from inside the computational domain.
c               set ndata=4, and the desired baseline (steady) pressure ratio 
c               (p/pinf), the  amplitude of the pressure oscillation
c               (deltap/pinf), the reduced frequency of the pressure
c               oscillation (rfreqp), and the "grid equivalent" length of the
c               dimensional reference length used to define rfreqp
c      bc2013 - same as bc2102, except also includes a phioff angle (in deg)
c               so that pset=p0(pratio+deltap*sin(2*pi*kr*t + phioff/57.2957795)
c      bc9999 - exact solution
c
c
c      II - Block Interface Boundary Conditions:
c
c      all block interface boundary conditions have bctype = 0
c
c           0 - 1-1 blocking (C-0 continuous),
c           0 - patched grid interpolation (point mismatch, no overlap),
c           0 - chimera grid interpolation (point overlap)
c           0 - embeded mesh (regular refinement of coarse mesh)
c 
c
c      NOTE: block interface type conditions (if set) will supercede
c            physical boundary conditions
c
c***********************************************************************
c       Notes on the use of subroutine BC to set boundary data:
c***********************************************************************
c
c       The boundary condition type definitions should be set on
c       each segment of a block face for any variations from the
c       default settings, which are type(a), as below, for every
c       boundary point. Boundary condition types are store in the
c       arrays BCI(..,1), BCJ(..,1), BCK(..,1) for the I=0, J=0 and
c       K=0 faces, respectively and BCI(..,2), BCJ(..,2), BCK(..,2)
c       or the I=IDIM, J=JDIM and K=KDIM faces, respectively
c
c       The boundary condition type is either:
c
c           (a) cell-center type, BCI/BCJ/BCK = 0.0 (DEFAULT)
c     
c           or
c
c           (b) interface type,   BCI/BCJ/BCK = 1.0
c
c       Type (a) assumes dependent variables corresponding to two ghost-cell
c       locations outside of the boundary are supplied. This type is commonly
c       used across boundary cuts where the solution outside the boundary can
c       be determined (i.e. symmetry planes, wake cuts, inflow/outflow,
c       abutments/interpolations from one block to another, etc.).
c
c       The values to the immediate left of the boundary are
c       put in arrays QI0(..,1)/QJ0(..,1)/QK0(..,1) for the I/J/K directions.
c       The values to the left of those are put in QI0(..,2)/QJ0(..,2)/
c       QK0(..,2). The values to the immediate right of the boundary are
c       put in arrays QI0(..,3)/QJ0(..,3)/QK0(..,3) for the I/J/K directions.
c       The values to the right of those are put in QI0(..,4)/QJ0(..,4)/
c       QK0(..,4).
c
c       Here, left for the boundary in the I-direction, for example,
c       refers to the grid points with I=1 and J,K variable. The right
c       boundary, correspondingly, refers to the grid points with I=IDIM.
c
c       Type (b) assumes that both the dependent variables and the
c       gradient of the dependent variables in computational space
c       at locations corresponding to the cell-interfaces of the boundary
c       are supplied. This type is commonly used at wall boundaries,
c       for example.
c
c       The interface values of the dependent variables at the immediate
c       left of the boundary are put in arrays QI0(..,1)/QJ0(..,1)/QK0(..,1)
c       for the I/J/K directions. The gradients of the dependent variables
c       with respect to the computational coordinate associated with the
c       direction away from the boundary of the grid are put in QI0(..,2)/
c       QJ0(..,2)/QK0(..,2). The interface values of the dependent variables
c       at the immediate right of the boundary are put in arrays QI0(..,3)/
c       QJ0(..,3)/QK0(..,3) for the I/J/K directions. The gradients of the
c       dependent variables with respect to the computational coordinate
c       associated with the direction away from the boundary of the grid are
c       put in QI0(..,4)/QJ0(..,4)/QK0(..,4).
c
c***********************************************************************
c
c     temporarily set isklton = 0 to supress output messages...
c     the messages are now output from bc_info
c
      isklt_sav = isklton
      isklton   = 0
c
      iuns = max(irotat(nbl),itrans(nbl),idefrm(nbl))
c
      if (isklton.gt.0) then
         if (iuns.eq.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),7)nbl,igridg(nbl)
         end if
         if (iuns.gt.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),8)nbl,igridg(nbl)
         end if
      end if
    7 format(30h boundary conditions for block,i6,1x,6h(grid ,i6,1h))
    8 format(30h boundary conditions for block,i6,1x,6h(grid ,i6,1h),
     .       1x,14h- dynamic mesh)
c
      igrid = igridg(nbl)
      jvdim = max(jdim*kdim,jdim*idim,kdim*idim)
      ixwk  = 1
c
c     check available storage in work array
c
      if (isklt_sav.gt.0) then
c        check storage for bc1003
         iflg = 0
         do 301 nseg=1,nbci0(nbl)
         if (ibcinfo(nbl,nseg,1,1).eq.1003) iflg = 1
 301     continue
         do 302 nseg=1,nbcidim(nbl)
         if (ibcinfo(nbl,nseg,1,2).eq.1003) iflg = 1
 302     continue
         do 303 nseg=1,nbcj0(nbl)
         if (jbcinfo(nbl,nseg,1,1).eq.1003) iflg = 1
 303     continue
         do 304 nseg=1,nbcjdim(nbl)
         if (jbcinfo(nbl,nseg,1,2).eq.1003) iflg = 1
 304     continue
         do 305 nseg=1,nbck0(nbl)
         if (kbcinfo(nbl,nseg,1,1).eq.1003) iflg = 1
 305     continue
         do 306 nseg=1,nbckdim(nbl)
         if (kbcinfo(nbl,nseg,1,2).eq.1003) iflg = 1
 306     continue
         if (iflg.gt.0) then  
            jvchk = ixwk + jvdim*23 - 1
            if (jvchk.gt.nwork) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)' stopping in bc...insufficient',
     .         ' wk storage for bctype 1003'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
c        check storage for bc2003
         iflg = 0
         do 311 nseg=1,nbci0(nbl)
         if (ibcinfo(nbl,nseg,1,1).eq.2003) iflg = 1
 311     continue
         do 312 nseg=1,nbcidim(nbl)
         if (ibcinfo(nbl,nseg,1,2).eq.2003) iflg = 1
 312     continue
         do 313 nseg=1,nbcj0(nbl)
         if (jbcinfo(nbl,nseg,1,1).eq.2003) iflg = 1
 313     continue
         do 314 nseg=1,nbcjdim(nbl)
         if (jbcinfo(nbl,nseg,1,2).eq.2003) iflg = 1
 314     continue
         do 315 nseg=1,nbck0(nbl)
         if (kbcinfo(nbl,nseg,1,1).eq.2003) iflg = 1
 315     continue
         do 316 nseg=1,nbckdim(nbl)
         if (kbcinfo(nbl,nseg,1,2).eq.2003) iflg = 1
 316     continue
         if (iflg.gt.0) then   
            jvchk = ixwk + jvdim*25 - 1
            if (jvchk.gt.nwork) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),*)' stopping in bc...insufficient',
     .         ' wk storage for bctype 2003'
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
         end if
c
      end if
c
c*******************************************************************************
c     Initialize density and pressure boundary conditions for comparison
c     later in boundary condition check.
c*******************************************************************************
c
      ibcflg = 0
      if (isklton.gt.0 .and. ntime.eq.1)  then
         istop = 1
         call bcchk(idim,jdim,kdim,w(lq),w(lqi0),w(lqj0),w(lqk0),
     .              w(lblk),ibcflg,nbl,nou,bou,nbuf,ibufdim,myid,istop,
     .              igridg,maxbl)
      end if
c
c*******************************************************************************
c     Initialize boundary condition types to all cell-center type.
c*******************************************************************************
c
      lsta = lbcj
      lend = lbci+jdim*kdim*2-1
      do 1 ll = lsta,lend
      w(ll) = 0.0
    1 continue
c
c*******************************************************************************
c      Update physical boundary conditions.
c*******************************************************************************
c
c********************
c     i=0 boundary
c********************
c
      do 802 nseg=1,nbci0(nbl)
      ista = 1
      iend = 1
      jsta = ibcinfo(nbl,nseg,2,1)
      jend = ibcinfo(nbl,nseg,3,1)
      ksta = ibcinfo(nbl,nseg,4,1)
      kend = ibcinfo(nbl,nseg,5,1)
      ldata= lwdat(nbl,nseg,1)
      mdim = jend-jsta
      ndim = kend-ksta
      filname = bcfiles(bcfilei(nbl,nseg,1))
c
      if (ibcinfo(nbl,nseg,1,1).eq.9999)
     .  call bc9999(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(ly),w(lz),w(lsnk0))
c
      if (ibcinfo(nbl,nseg,1,1).eq.1000)
     .  call bc1000(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c     noninertial boundary conditions
      if (ibcinfo(nbl,nseg,1,1).eq.1000.and.noninflag.gt.0) then
        call bcnonin(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,
     .  iuns,nou,bou,nbuf,ibufdim,
     .  w(lx),w(ly),w(lz),nbl)
      endif
c
      if (ibcinfo(nbl,nseg,1,1).eq.1001)
     .  call bc1001(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(lz))
c
      if (ibcinfo(nbl,nseg,1,1).eq.1002)
     .  call bc1002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1003)
     .  call bc1003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,1,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,w(lx),w(lz),cl,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1005)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,0,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1006)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,1,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1008)
     .  call bc1008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.1011) then
        if(isklt_sav.gt.0) 
     .  call chksym(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,1,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,isym,
     .  jsym,ksym,nou,bou,nbuf,ibufdim,myid)
        call bc1011(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),isym,jsym,ksym,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (ibcinfo(nbl,nseg,1,1).eq.1012) then
        if(isklt_sav.gt.0)
     .  call chkrap(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,1,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,iwrap,
     .  jwrap,kwrap,nou,bou,nbuf,ibufdim,myid)
        call bc1012(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),iwrap,jwrap,kwrap,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (ibcinfo(nbl,nseg,1,1).eq.1013)
     .  call bc1013(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2002)
     .  call bc2002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2102)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,0)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2103)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,1)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2003)
     .  call bc2003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,1,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2009)
     .  call bc2009(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,1,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2010)
     .  call bc2010(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,1,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2019)
     .  call bc2019(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,1,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (abs(ibcinfo(nbl,nseg,1,1)).eq.2004 .or.
     .    abs(ibcinfo(nbl,nseg,1,1)).eq.2014)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,ibcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,0)
c
      if (abs(ibcinfo(nbl,nseg,1,1)).eq.2024)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,ibcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,2)
c
      if (abs(ibcinfo(nbl,nseg,1,1)).eq.2034)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,ibcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,1)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ichk = 1
        ngnew = 0
        if (isklt_sav.gt.0 .and. ntime.eq.1) ichk = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
           if (ngnew .lt. nbl) ichk = 1
        end if
        call bc2006(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  w(lx),w(ly),w(lz),ngnew,jdimg(ngnew),kdimg(ngnew),idimg(ngnew),
     .  w(lw(2,ngnew)),w(lw(3,ngnew)),w(lw(4,ngnew)),ichk,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
      end if
c
      if (ibcinfo(nbl,nseg,1,1).eq.2007)
     .  call bc2007(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2008 .or.
     .    ibcinfo(nbl,nseg,1,1).eq.2018 .or.
     .    ibcinfo(nbl,nseg,1,1).eq.2028 .or.
     .    ibcinfo(nbl,nseg,1,1).eq.2038)
     .  call bc2008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,ibcinfo(nbl,nseg,1,1),nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2016) 
     .  call bc2016(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,ibcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,1).eq.2026)
     .  call bc2026(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,1,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,ibcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
  802 continue
c
c********************
c     i=idim boundary
c********************
c
      do 803 nseg=1,nbcidim(nbl)
      ista = idim
      iend = idim
      jsta = ibcinfo(nbl,nseg,2,2)
      jend = ibcinfo(nbl,nseg,3,2)
      ksta = ibcinfo(nbl,nseg,4,2)
      kend = ibcinfo(nbl,nseg,5,2)
      ldata= lwdat(nbl,nseg,2)
      mdim = jend-jsta
      ndim = kend-ksta
      filname = bcfiles(bcfilei(nbl,nseg,2))
c
      if (ibcinfo(nbl,nseg,1,2).eq.9999)
     .  call bc9999(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(ly),w(lz),w(lsnk0))
c
      if (ibcinfo(nbl,nseg,1,2).eq.1000)
     .  call bc1000(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c     noninertial boundary conditions
      if (ibcinfo(nbl,nseg,1,2).eq.1000.and.noninflag.gt.0) then
        call bcnonin(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,
     .  iuns,nou,bou,nbuf,ibufdim,
     .  w(lx),w(ly),w(lz),nbl)
      endif
c
      if (ibcinfo(nbl,nseg,1,2).eq.1001)
     .  call bc1001(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(lz))
c
      if (ibcinfo(nbl,nseg,1,2).eq.1002)
     .  call bc1002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1003)
     .  call bc1003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,2,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,w(lx),w(lz),cl,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1005)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,0,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1006)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,1,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1008)
     .  call bc1008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.1011) then
        if(isklt_sav.gt.0)
     .  call chksym(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,2,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,isym,
     .  jsym,ksym,nou,bou,nbuf,ibufdim,myid)
        call bc1011(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),isym,jsym,ksym,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (ibcinfo(nbl,nseg,1,2).eq.1012) then
        if(isklt_sav.gt.0)
     .  call chkrap(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,2,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,iwrap,
     .  jwrap,kwrap,nou,bou,nbuf,ibufdim,myid)
        call bc1012(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),iwrap,jwrap,kwrap,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (ibcinfo(nbl,nseg,1,2).eq.1013)
     .  call bc1013(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2002)
     .  call bc2002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2102)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,0)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2103)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,1)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2003)
     .  call bc2003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,2,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2009)
     .  call bc2009(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,2,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2010)
     .  call bc2010(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,2,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2019)
     .  call bc2019(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,2,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (abs(ibcinfo(nbl,nseg,1,2)).eq.2004 .or.
     .    abs(ibcinfo(nbl,nseg,1,2)).eq.2014)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,ibcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,0)
c
      if (abs(ibcinfo(nbl,nseg,1,2)).eq.2024)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,ibcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,2)
c
      if (abs(ibcinfo(nbl,nseg,1,2)).eq.2034)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,ibcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,1)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ichk = 1
        ngnew = 0
        if (isklt_sav.gt.0 .and. ntime.eq.1) ichk = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
           if (ngnew .lt. nbl) ichk = 1
        end if
        call bc2006(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  w(lx),w(ly),w(lz),ngnew,jdimg(ngnew),kdimg(ngnew),idimg(ngnew),
     .  w(lw(2,ngnew)),w(lw(3,ngnew)),w(lw(4,ngnew)),ichk,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
      end if
c
      if (ibcinfo(nbl,nseg,1,2).eq.2007)
     .  call bc2007(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2008 .or.
     .    ibcinfo(nbl,nseg,1,2).eq.2018 .or.
     .    ibcinfo(nbl,nseg,1,2).eq.2028 .or.
     .    ibcinfo(nbl,nseg,1,2).eq.2038)
     .  call bc2008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,ibcinfo(nbl,nseg,1,2),nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2016)
     .  call bc2016(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,ibcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
      if (ibcinfo(nbl,nseg,1,2).eq.2026)
     .  call bc2026(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,2,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,ibcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
  803 continue
c
c********************
c     j=0 boundary
c********************
c
      do 804 nseg=1,nbcj0(nbl)
      ista = jbcinfo(nbl,nseg,2,1)
      iend = jbcinfo(nbl,nseg,3,1)
      jsta = 1
      jend = 1
      ksta = jbcinfo(nbl,nseg,4,1)
      kend = jbcinfo(nbl,nseg,5,1)
      ldata= lwdat(nbl,nseg,3)
      mdim = kend-ksta
      ndim = iend-ista
      filname = bcfiles(bcfilej(nbl,nseg,1))
c
      if (jbcinfo(nbl,nseg,1,1).eq.9999)
     .  call bc9999(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(ly),w(lz),w(lsnk0))
c
      if (jbcinfo(nbl,nseg,1,1).eq.1000)
     .  call bc1000(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c     noninertial boundary conditions
      if (jbcinfo(nbl,nseg,1,1).eq.1000.and.noninflag.gt.0) then
        call bcnonin(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,
     .  iuns,nou,bou,nbuf,ibufdim,
     .  w(lx),w(ly),w(lz),nbl)
      endif
c
      if (jbcinfo(nbl,nseg,1,1).eq.1001)
     .  call bc1001(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(lz))
c
      if (jbcinfo(nbl,nseg,1,1).eq.1002)
     .  call bc1002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1003)
     .  call bc1003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,3,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,w(lx),w(lz),cl,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1005)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,0,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1006)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,1,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1008)
     .  call bc1008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.1011) then
        if(isklt_sav.gt.0)
     .  call chksym(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,3,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,isym,
     .  jsym,ksym,nou,bou,nbuf,ibufdim,myid)
        call bc1011(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),isym,jsym,ksym,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (jbcinfo(nbl,nseg,1,1).eq.1012) then
        if(isklt_sav.gt.0)
     .  call chkrap(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,3,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,iwrap,
     .  jwrap,kwrap,nou,bou,nbuf,ibufdim,myid)
        call bc1012(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),iwrap,jwrap,kwrap,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (jbcinfo(nbl,nseg,1,1).eq.1013)
     .  call bc1013(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2002)
     .  call bc2002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2102)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,0)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2103)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,1)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2003)
     .  call bc2003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,3,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2009)
     .  call bc2009(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,3,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2010)
     .  call bc2010(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,3,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2019)
     .  call bc2019(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,3,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (abs(jbcinfo(nbl,nseg,1,1)).eq.2004 .or.
     .    abs(jbcinfo(nbl,nseg,1,1)).eq.2014)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,jbcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,0)
c
      if (abs(jbcinfo(nbl,nseg,1,1)).eq.2024)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,jbcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,2)
c
      if (abs(jbcinfo(nbl,nseg,1,1)).eq.2034)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,jbcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,1)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ichk = 1
        ngnew = 0
        if (isklt_sav.gt.0 .and. ntime.eq.1) ichk = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
           if (ngnew .lt. nbl) ichk = 1
        end if
        call bc2006(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  w(lx),w(ly),w(lz),ngnew,jdimg(ngnew),kdimg(ngnew),idimg(ngnew),
     .  w(lw(2,ngnew)),w(lw(3,ngnew)),w(lw(4,ngnew)),ichk,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
      end if
c
      if (jbcinfo(nbl,nseg,1,1).eq.2007)
     .  call bc2007(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2008 .or.
     .    jbcinfo(nbl,nseg,1,1).eq.2018 .or.
     .    jbcinfo(nbl,nseg,1,1).eq.2028 .or.
     .    jbcinfo(nbl,nseg,1,1).eq.2038)
     .  call bc2008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,jbcinfo(nbl,nseg,1,1),nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2016)
     .  call bc2016(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,jbcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,1).eq.2026)
     .  call bc2026(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,3,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,jbcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
  804 continue
c
c********************
c     j=jdim boundary
c********************
c
      do 805 nseg=1,nbcjdim(nbl)
      ista = jbcinfo(nbl,nseg,2,2)
      iend = jbcinfo(nbl,nseg,3,2)
      jsta = jdim
      jend = jdim 
      ksta = jbcinfo(nbl,nseg,4,2)
      kend = jbcinfo(nbl,nseg,5,2)
      ldata= lwdat(nbl,nseg,4)
      mdim = kend-ksta
      ndim = iend-ista
      filname = bcfiles(bcfilej(nbl,nseg,2))
c
      if (jbcinfo(nbl,nseg,1,2).eq.9999)
     .  call bc9999(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(ly),w(lz),w(lsnk0))
c
      if (jbcinfo(nbl,nseg,1,2).eq.1000)
     .  call bc1000(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c     noninertial boundary conditions
      if (jbcinfo(nbl,nseg,1,2).eq.1000.and.noninflag.gt.0) then
        call bcnonin(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,
     .  iuns,nou,bou,nbuf,ibufdim,
     .  w(lx),w(ly),w(lz),nbl)
      endif
c
      if (jbcinfo(nbl,nseg,1,2).eq.1001)
     .  call bc1001(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(lz))
c
      if (jbcinfo(nbl,nseg,1,2).eq.1002)
     .  call bc1002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1003)
     .  call bc1003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,4,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,w(lx),w(lz),cl,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1005)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,0,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1006)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,1,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1008)
     .  call bc1008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.1011) then
        if(isklt_sav.gt.0)
     .  call chksym(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,4,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,isym,
     .  jsym,ksym,nou,bou,nbuf,ibufdim,myid)
        call bc1011(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),isym,jsym,ksym,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (jbcinfo(nbl,nseg,1,2).eq.1012) then
        if(isklt_sav.gt.0)
     .  call chkrap(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,4,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,iwrap,
     .  jwrap,kwrap,nou,bou,nbuf,ibufdim,myid)
        call bc1012(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),iwrap,jwrap,kwrap,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (jbcinfo(nbl,nseg,1,2).eq.1013)
     .  call bc1013(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2002)
     .  call bc2002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2102)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,0)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2103)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,1)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2003)
     .  call bc2003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,4,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2009)
     .  call bc2009(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,4,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2010)
     .  call bc2010(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,4,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2019)
     .  call bc2019(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,4,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (abs(jbcinfo(nbl,nseg,1,2)).eq.2004 .or.
     .    abs(jbcinfo(nbl,nseg,1,2)).eq.2014)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,jbcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,0)
c
      if (abs(jbcinfo(nbl,nseg,1,2)).eq.2024)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,jbcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,2)
c
      if (abs(jbcinfo(nbl,nseg,1,2)).eq.2034)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,jbcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,1)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ichk = 1
        ngnew = 0
        if (isklt_sav.gt.0 .and. ntime.eq.1) ichk = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
           if (ngnew .lt. nbl) ichk = 1
        end if
        call bc2006(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  w(lx),w(ly),w(lz),ngnew,jdimg(ngnew),kdimg(ngnew),idimg(ngnew),
     .  w(lw(2,ngnew)),w(lw(3,ngnew)),w(lw(4,ngnew)),ichk,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
      end if
c
      if (jbcinfo(nbl,nseg,1,2).eq.2007)
     .  call bc2007(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2008 .or.
     .    jbcinfo(nbl,nseg,1,2).eq.2018 .or.
     .    jbcinfo(nbl,nseg,1,2).eq.2028 .or.
     .    jbcinfo(nbl,nseg,1,2).eq.2038)
     .  call bc2008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,jbcinfo(nbl,nseg,1,2),nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2016)
     .  call bc2016(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,jbcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
      if (jbcinfo(nbl,nseg,1,2).eq.2026)
     .  call bc2026(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,4,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,jbcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
  805 continue
c
c********************
c     k=0 boundary
c********************
c
      do 806 nseg=1,nbck0(nbl)
      ista = kbcinfo(nbl,nseg,2,1)
      iend = kbcinfo(nbl,nseg,3,1)
      jsta = kbcinfo(nbl,nseg,4,1)
      jend = kbcinfo(nbl,nseg,5,1)
      ksta = 1
      kend = 1
      ldata= lwdat(nbl,nseg,5)
      mdim = jend-jsta
      ndim = iend-ista
      filname = bcfiles(bcfilek(nbl,nseg,1))
c
      if (kbcinfo(nbl,nseg,1,1).eq.9999)
     .  call bc9999(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(ly),w(lz),w(lsnk0))
c
      if (kbcinfo(nbl,nseg,1,1).eq.1000)
     .  call bc1000(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c     noninertial boundary conditions
      if (kbcinfo(nbl,nseg,1,1).eq.1000.and.noninflag.gt.0) then
        call bcnonin(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,
     .  iuns,nou,bou,nbuf,ibufdim,
     .  w(lx),w(ly),w(lz),nbl)
      endif
c
      if (kbcinfo(nbl,nseg,1,1).eq.1001)
     .  call bc1001(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(lz))
c
      if (kbcinfo(nbl,nseg,1,1).eq.1002)
     .  call bc1002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1003)
     .  call bc1003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,5,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,w(lx),w(lz),cl,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1005)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,0,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1006)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,1,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1008)
     .  call bc1008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.1011) then
        if(isklt_sav.gt.0)
     .  call chksym(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,5,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,isym,
     .  jsym,ksym,nou,bou,nbuf,ibufdim,myid)
        call bc1011(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),isym,jsym,ksym,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (kbcinfo(nbl,nseg,1,1).eq.1012) then
        if(isklt_sav.gt.0)
     .  call chkrap(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,5,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,iwrap,
     .  jwrap,kwrap,nou,bou,nbuf,ibufdim,myid)
        call bc1012(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),iwrap,jwrap,kwrap,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (kbcinfo(nbl,nseg,1,1).eq.1013)
     .  call bc1013(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2002)
     .  call bc2002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2102)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,0)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2103)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,1)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2003)
     .  call bc2003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,5,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2009)
     .  call bc2009(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,5,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2010)
     .  call bc2010(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,5,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2019)
     .  call bc2019(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,5,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (abs(kbcinfo(nbl,nseg,1,1)).eq.2004 .or.
     .    abs(kbcinfo(nbl,nseg,1,1)).eq.2014)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,kbcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,0)
c
      if (abs(kbcinfo(nbl,nseg,1,1)).eq.2024)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,kbcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,2)
c
      if (abs(kbcinfo(nbl,nseg,1,1)).eq.2034)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,kbcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,1)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ichk = 1
        ngnew = 0
        if (isklt_sav.gt.0 .and. ntime.eq.1) ichk = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
           if (ngnew .lt. nbl) ichk = 1
        end if
        call bc2006(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  w(lx),w(ly),w(lz),ngnew,jdimg(ngnew),kdimg(ngnew),idimg(ngnew),
     .  w(lw(2,ngnew)),w(lw(3,ngnew)),w(lw(4,ngnew)),ichk,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
      end if
c
      if (kbcinfo(nbl,nseg,1,1).eq.2007)
     .  call bc2007(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2008 .or.
     .    kbcinfo(nbl,nseg,1,1).eq.2018 .or.
     .    kbcinfo(nbl,nseg,1,1).eq.2028 .or.
     .    kbcinfo(nbl,nseg,1,1).eq.2038)
     .  call bc2008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,kbcinfo(nbl,nseg,1,1),nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2016)
     .  call bc2016(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,kbcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,1).eq.2026)
     .  call bc2026(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,5,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,kbcinfo(nbl,nseg,1,1),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
  806 continue
c
c********************
c     k=kdim boundary
c********************
c
      do 807 nseg=1,nbckdim(nbl)
      ista = kbcinfo(nbl,nseg,2,2)
      iend = kbcinfo(nbl,nseg,3,2)
      jsta = kbcinfo(nbl,nseg,4,2)
      jend = kbcinfo(nbl,nseg,5,2)
      ksta = kdim
      kend = kdim
      ldata= lwdat(nbl,nseg,6)
      mdim = jend-jsta
      ndim = iend-ista
      filname = bcfiles(bcfilek(nbl,nseg,2))
c
      if (kbcinfo(nbl,nseg,1,2).eq.9999)
     .  call bc9999(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(ly),w(lz),w(lsnk0))
c
      if (kbcinfo(nbl,nseg,1,2).eq.1000)
     .  call bc1000(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c     noninertial boundary conditions
      if (kbcinfo(nbl,nseg,1,2).eq.1000.and.noninflag.gt.0) then
        call bcnonin(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,
     .  iuns,nou,bou,nbuf,ibufdim,
     .  w(lx),w(ly),w(lz),nbl)
      endif
c
      if (kbcinfo(nbl,nseg,1,2).eq.1001)
     .  call bc1001(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem,w(lx),w(lz))
c
      if (kbcinfo(nbl,nseg,1,2).eq.1002)
     .  call bc1002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1003)
     .  call bc1003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,6,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,w(lx),w(lz),cl,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1005)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,0,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1006)
     .  call bc1005(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,1,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1008)
     .  call bc1008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.1011) then
        if(isklt_sav.gt.0)
     .  call chksym(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,6,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,isym,
     .  jsym,ksym,nou,bou,nbuf,ibufdim,myid)
        call bc1011(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),isym,jsym,ksym,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (kbcinfo(nbl,nseg,1,2).eq.1012) then
        if(isklt_sav.gt.0)
     .  call chkrap(nbl,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .  ibcinfo,jbcinfo,kbcinfo,6,idim,jdim,kdim,
     .  maxbl,maxseg,ista,iend,jsta,jend,ksta,kend,iwrap,
     .  jwrap,kwrap,nou,bou,nbuf,ibufdim,myid)
        call bc1012(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),
     .  w(lti0),w(lvis),w(lvj0),w(lvk0),w(lvi0),iwrap,jwrap,kwrap,iuns,
     .  nou,bou,nbuf,ibufdim,nummem)
      end if
c
      if (kbcinfo(nbl,nseg,1,2).eq.1013)
     .  call bc1013(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),iuns,nou,bou,nbuf,ibufdim,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2002)
     .  call bc2002(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2102)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,0)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2103)
     .  call bc2102(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nbl,nou,bou,nbuf,ibufdim,myid,nummem,1)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2003)
     .  call bc2003(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,wk(ixwk),jvdim,6,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2009)
     .  call bc2009(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,6,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2010)
     .  call bc2010(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,6,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2019)
     .  call bc2019(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  ista,iend,jsta,jend,ksta,kend,6,
     .  w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (abs(kbcinfo(nbl,nseg,1,2)).eq.2004 .or.
     .    abs(kbcinfo(nbl,nseg,1,2)).eq.2014)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,kbcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,0)
c
      if (abs(kbcinfo(nbl,nseg,1,2)).eq.2024)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,kbcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,2)
c
      if (abs(kbcinfo(nbl,nseg,1,2)).eq.2034)
     .  call bc2004(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,kbcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem,1)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2006) then
c       ngnew is the block number from which integration of the
c       radial equilibrium eq. is continued into block nbl (if applicable)
        ngh = igridg(nbl)
        nblnum = int(w(ldata))
        ichk = 1
        ngnew = 0
        if (isklt_sav.gt.0 .and. ntime.eq.1) ichk = 0
        if (nblnum .gt. 0) then
           ngnew = nblg(nblnum) + (nbl - nblg(ngh))
           if (ngnew .lt. nbl) ichk = 1
        end if
        call bc2006(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  w(lx),w(ly),w(lz),ngnew,jdimg(ngnew),kdimg(ngnew),idimg(ngnew),
     .  w(lw(2,ngnew)),w(lw(3,ngnew)),w(lw(4,ngnew)),ichk,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
      end if
c
      if (kbcinfo(nbl,nseg,1,2).eq.2007)
     .  call bc2007(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2008 .or.
     .    kbcinfo(nbl,nseg,1,2).eq.2018 .or.
     .    kbcinfo(nbl,nseg,1,2).eq.2028 .or.
     .    kbcinfo(nbl,nseg,1,2).eq.2038)
     .  call bc2008(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),filname,iuns,
     .  nou,bou,nbuf,ibufdim,myid,kbcinfo(nbl,nseg,1,2),nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2016)
     .  call bc2016(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,kbcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
      if (kbcinfo(nbl,nseg,1,2).eq.2026)
     .  call bc2026(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .  w(lsj),w(lsk),w(lsi),w(lbcj),w(lbck),w(lbci),
     .  w(lxtbj),w(lxtbk),w(lxtbi),w(latbj),w(latbk),w(latbi),
     .  ista,iend,jsta,jend,ksta,kend,6,w(lxib),w(ltj0),w(ltk0),w(lti0),
     .  w(lsnk0),w(lvis),w(lvj0),w(lvk0),w(lvi0),mdim,ndim,w(ldata),
     .  filname,iuns,kbcinfo(nbl,nseg,1,2),w(lsni0),w(lsnk0),w(lsni0),
     .  ntime,w(lnbl),w(lxkb),w(lnbl),nou,bou,nbuf,ibufdim,myid,
     .  nummem)
c
  807 continue
c
      isklton = isklt_sav
c
      return
      end
