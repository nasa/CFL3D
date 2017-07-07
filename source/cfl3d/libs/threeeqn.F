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
      subroutine threeeqn(jdim,kdim,idim,q,sj,sk,si,vol,dtj,
     + x,y,z,vist3d,vor,smin,zksav,turre,damp1,blend,timestp,fnu,bx,bx2,
     + cx,cx2,dx,dx2,fx,fx2,workx,by,by2,cy,cy2,dy,dy2,fy,fy2,worky,
     + bz,bz2,cz,cz2,dz,dz2,fz,fz2,workz,ntime,tj0,tk0,ti0,nbl,
     + qj0,qk0,qi0,vj0,vk0,vi0,blank,iover,sumn1,sumn2,sumn3,
     + negn1,negn2,negn3,
     + ux,rhside,zksav2,v3dtmp,bcj,bck,bci,nbci0,nbcidim,
     + nbcj0,nbcjdim,nbck0,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     + maxbl,maxseg,volj0,volk0,voli0,nou,bou,nbuf,ibufdim,iex,iex2,
     + iex3,bx3,cx3,dx3,fx3,by3,cy3,dy3,fy3,bz3,cz3,dz3,fz3,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute turbulent viscosity distributions using
c     3-equation turbulence models
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /lam/ ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi,
     .        i_lam_forcezero
      common /maxiv/ ivmx
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /twod/ i2d
      common /zero/ iexp
      common /wallfun/ iwf(3)
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /konew/ ikoprod,isstdenom,pklimterm,ibeta8kzeta,i_bsl,
     .        keepambient,re_thetat0,i_wilcox06,i_wilcox06_chiw,
     .        i_turbprod_kterm,i_catris_kw,prod2d3dtrace,
     .        i_compress_correct,isstsf,i_wilcox98,i_wilcox98_chiw,
     .        isst2003
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
c
      dimension q(jdim,kdim,idim,5),sj(jdim,kdim,idim-1,5),
     + sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
     +,dtj(jdim,kdim,idim-1),x(jdim,kdim,idim),y(jdim,kdim,idim),
     + z(jdim,kdim,idim),vist3d(jdim,kdim,idim),
     + vor(jdim-1,kdim-1,idim-1),smin(jdim-1,kdim-1,idim-1)
      dimension damp1(jdim-1,kdim-1,idim-1),
     + fnu(0:jdim,0:kdim,0-iex3:idim+iex3)
     + , blend(jdim-1,kdim-1,idim-1),timestp(jdim-1,kdim-1,idim-1)
      dimension bx(kdim-1,jdim-1),bx2(kdim-1,jdim-1),
     + cx(kdim-1,jdim-1),cx2(kdim-1,jdim-1),dx(kdim-1,jdim-1),
     + dx2(kdim-1,jdim-1),fx(kdim-1,jdim-1),fx2(kdim-1,jdim-1),
     + workx(kdim-1,jdim-1),
     +          by(jdim-1,kdim-1),by2(jdim-1,kdim-1),
     + cy(jdim-1,kdim-1),cy2(jdim-1,kdim-1),dy(jdim-1,kdim-1),
     + dy2(jdim-1,kdim-1),fy(jdim-1,kdim-1),fy2(jdim-1,kdim-1),
     + worky(jdim-1,kdim-1),
     +          bz(kdim-1,idim-1),bz2(kdim-1,idim-1),
     + cz(kdim-1,idim-1),cz2(kdim-1,idim-1),dz(kdim-1,idim-1),
     + dz2(kdim-1,idim-1),fz(kdim-1,idim-1),fz2(kdim-1,idim-1),
     + workz(kdim-1,idim-1)
       dimension bx3(kdim-1,jdim-1),cx3(kdim-1,jdim-1),
     +  dx3(kdim-1,jdim-1),fx3(kdim-1,jdim-1),
     +           by3(jdim-1,kdim-1),cy3(jdim-1,kdim-1),
     +  dy3(jdim-1,kdim-1),fy3(jdim-1,kdim-1),
     +           bz3(kdim-1,idim-1),cz3(kdim-1,idim-1),
     +  dz3(kdim-1,idim-1),fz3(kdim-1,idim-1)
      dimension turre(0-iex:jdim+iex,0-iex:kdim+iex,0-iex2:idim+iex2,3),
     + blank(jdim,kdim,idim),
     + zksav(jdim,kdim,idim,nummem),rhside(jdim-1,kdim-1,idim-1,3)
     + ,tj0(kdim,idim-1,nummem,4),tk0(jdim,idim-1,nummem,4),
     +  ti0(jdim,kdim,nummem,4),
     + qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4)
     +,vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),vi0(jdim,kdim,1,4)
      dimension ux(jdim-1,kdim-1,idim-1,9),
     + zksav2(jdim,kdim,idim,2*nummem),
     + v3dtmp(0:jdim,0:kdim,0-iex3:idim+iex3)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     +          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     +          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension volj0(kdim,idim-1,4),
     +          volk0(jdim,idim-1,4),voli0(jdim,kdim,4)
c
c
c   Variables:
c
c      jdim,kdim,idim - dimensions of this block
c      q - primitive variables (rho, u, v, w, p)
c      sj,sk,si - metric terms (defined on cell FACES)
c      vol - cell volume
c      volj0,volk0,voli0 - ghost-cell volumes
c      dtj - vol/dt
c      x,y,z - grid vertex locations
c      vist3d - turbulent eddy viscosity (nondimensionalized by mulamref)
c      vor - vorticity
c      smin - minimum distance to any wall - set negative if laminar region!
c      zksav - saved values of the 3 turbulent quantities (1=omega or epsilon,
c              2=k, 3=intermittency)
c                omega is nondimensionalized by rhoref*aref**2/mulamref
c                epsilon          "             rhoref*aref**4/mulamref
c                k                "             aref**2
c      zksav2 - auxiliary k and omega (or epsilon) for saving quantities
c               from last time step during time-accurate subiterations
c               (k,j,i,1), (j,k,i,2), and (j,k,i,3) = k and omega (or epsilon) 
c                   and intermittency
c               (k,j,i,4), (j,k,i,5), and (j,k,i,6) = delta k and delta omega 
c                   (or epsilon) and delta intermittency
c      turre - working values of 3 turb quantities in this subroutine
c      damp1 - for SST=CD (cross-derivative) term; for others, sometimes used
c              as storage for linearizations of production terms for addition
c              to LHS diagonal elements
c      blend - blending term for 2-layer models (just SST right now); 
c              =1 otherwise
c      timestp - time step array delta t
c      fnu - laminar viscosity (from Sutherland's law), nondim by mulamref
c      bx,bx2,cx,cx2,dx,dx2,fx,fx2 - sub, diag, superdiag, & RHS in 
c                                    eta direction
c      by,by2,cy,cy2,dy,dy2,fy,fy2 - sub, diag, superdiag, & RHS in 
c                                    xi direction
c      bz,bz2,cz,cz2,dz,dz2,fz,fz2 - sub, diag, superdiag, & RHS in 
c                                    zeta direction
c      workx,worky,workz - work arrays for tridiagonal solvers
c      ntime - time counter
c      tj0,tk0,ti0 - BCs for turbulent quantities (1=omega or epsilon, 2=k)
c      nbl = block number currently working on
c      qj0,qk0,qi0 - BCs for q's
c      vj0,vk0,vi0 - BCs for vist3d
c      blank - iblanking array for overset
c      iover - overset gridding parameter
c      sumn1,sumn2,sunm3 - residual for 2 turbulence equations
c      negn1,negn2,negn3 - number of locations where the solution yields a "negative"
c                    turbulence quantity.  This SHOULD be zero, but it seems
c                    to be OK if there are only a few of these.  When they
c                    go negative, the values are artificially limited to be
c                    > 0.  When these numbers are large, it indicates that
c                    the solution is probably going to blow up.  Check your
c                    grid for excessive grid stretching, or try lowering CFL.
c      ux - 9 components of velocity derivative: ux,uy,uz,vx,vy,vz,wx,wy,wz 
c           at cell centers.  Used for nonlinear models only.
c      rhside - right-hand-side terms for 3 eqns
c      v3dtmp - temporary storage for vist3d (needed for nonlinear models,
c               which require turb viscosity in diffusion term WITHOUT
c               variable cmu)
c      bcj,bck,bci - =0 for cell-center BC, =1 for at-face-center BC
c      nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim - no. of BC segments each face
c      ibcinfo,jbcinfo,kbcinfo - gives i,j,k start and end indices for each
c               BC segment, among other things
c      maxbl,maxseg - dimensions used in nbci0,etc and ibcinfo, etc arrays
c
c     Parameter if need to revert to standard SST w no transition (itrans_on=0):
      itrans_on=1
c     Set some of the transition-related parameters:
c     Note: ivmx=31 is currently NOT usable at all!
      tu_percent=100.*sqrt(2.*tur10(2)/(3.*xmach*xmach))
      if (real(re_thetat0) .lt. 0.) then
        if (ivmx .eq. 30) then
          re_thetat=400.*tu_percent**(-5./8.)
c         re_thetat=32.591*(tu_percent/100.)**(-0.621)
        else if (ivmx .eq. 31) then
          re_thetat=803.73*(tu_percent+0.6067)**(-1.027)
        end if
      else
        re_thetat=re_thetat0
      end if
      if (ivmx .eq. 30) then
        re_thetac=0.85*re_thetat
      else if (ivmx .eq. 31) then
        re_thetac=re_thetat
        flength=163.*alog(real(tu_percent))+3.625
        ce2=50.
      end if
c
      if(isklton .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     Computing turbulent'',
     +    '' viscosity using 3-eqns, block='',i5)') nbl
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     WARNING: ivisc=30 still'',
     .    '' under development... use at your own risk!'')')
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     Freestream tur10,tur20,tur30 = '',
     +     3e19.8)') tur10(1),tur10(2),tur10(3)
        if(iturbord .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     1st order advection on RHS'')')
        else
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     2nd order advection on RHS'')')
        end if
      end if
      if(isklton .gt. 0) then
        if(ivmx .eq. 30) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-omega SST (Menter)'',
     +       '' + transition model (Langtry unpublished version)'')')
        else if (ivmx .eq. 31) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-omega SST (Menter)'',
     +       '' + transition model (JTurbomach 128:413 2006 and'',
     +       '' IJHFF 29:48 2008)'')')
        end if
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''     Tu (percent)='',f10.4)') 
     +    tu_percent
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''     mut_inf/muref='',f10.4,
     +    '', Re_thetat='',f10.4)') rho0*tur10(2)/tur10(1),re_thetat
        if(ikoprod .eq. 1) then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''     strain-based production'',
     +    '' term'')')
        else
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''     approx (vort) production'',
     .      '' term'')')
        end if
        if(isstdenom .eq. 0) then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''     vort in denom of mut term'')')
        else
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''     S in denom of mut term'')')
        end if
c
        if(abs(prod2d3dtrace-0.5) .lt. 0.01) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     Sij used in 2SijSij prod'',
     .         '' term forced to be traceless in 2-D sense'')')
        else if(abs(prod2d3dtrace-0.33333333) .lt. 0.01) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     Sij used in 2SijSij prod'',
     .         '' term forced to be traceless in 3-D sense'')')
        end if
c
        if(keepambient .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     3-eqn ambient turbulence'',
     .     '' levels not allowed to decay'')')
        end if
c
        if(iaxi2planeturb .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     3-eqn model ignoring i-dir'')')
        end if
        if(istrongturbdis .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     strong conserv - diss terms'')')
        end if
c
      end if
c
c Note: (10.**(-iexp) is machine zero)
      xminn=10.**(-iexp+1)
c
c Set number of subiterations to solve turbulence field eqn per iteration
c (usually, 1 is sufficient... but if residual diverges then may need more)
c
      nsubit=nsubturb
c
c Set factors that multiply the N-S CFL number, to determine CFL number
c for the turbulence model (this typically can be 2 - 10... optimum value
c is a function of the case).  If turb model seems to be having trouble
c converging, lowering these factors may be one strategy to try. NOTE: factors
c used only for steady state cases, not time accurate cases.
c They can be overridden with keyword "cflturb" if want all CFLs to be the
c same (not recommended), or by cflturb1, cflturb2, etc to set individual
c ones.
c
      factor1=10.
      factor2=10.
      factor3=100.
c
c Overwrite factors with keyword value "cflturb()" if nonzero
c
      if (real(cflturb(1)).ne.0.) then
         factor1 = cflturb(1)
      end if
      if (real(cflturb(2)).ne.0.) then
         factor2 = cflturb(2)
      end if
      if (real(cflturb(3)).ne.0.) then
         factor3 = cflturb(3)
      end if
c factor2 and 3 are set relative to factor1
      factor2=factor2/factor1
      factor3=factor3/factor1
c
c Timestep for turb model
c
      if (real(dt).lt.0) then
         do i=1,idim-1
         do k=1,kdim-1
         do j=1,jdim-1
           timestp(j,k,i)=factor1*vol(j,k,i)/dtj(j,k,i)
           timestp(j,k,i)=ccmincr(timestp(j,k,i),100.)
         enddo
         enddo
         enddo
      else
c          turbulence model advanced with physical time only
c          (pseudo-time term NOT included, even for tau-TS in mean-
c          flow equations, since multigrid is not used for turb. eq.)
         do i=1,idim-1
         do k=1,kdim-1
         do j=1,jdim-1
           timestp(j,k,i)=dt
           factor2=1.
           factor3=1.
         enddo
         enddo
         enddo
      end if
c
c Set up constants
      vk =.41
      a1 =.31
c Constants for Set 1:
c  cmuc1 (normally 0.09):
      cmuc1=0.09
c  beta1 (constant in omega or epsilon destruction term):
      if(ivmx.eq.30 .or. ivmx.eq.31)                 beta1=0.075
c  sigo1 (constant in omega, epsilon, or zeta diffusion term):
      if(ivmx.eq.30 .or. ivmx.eq.31)                 sigo1=0.5
c  sigk1 (constant in k diffusion term):
      if(ivmx.eq.30 .or. ivmx.eq.31)                 sigk1=0.85
c  alp1 (constant in omega or epsilon production term):
      if(ivmx.eq.30 .or. ivmx.eq.31)                 alp1 =
     +    beta1/cmuc1 - sigo1*vk*vk/sqrt(cmuc1)
c
c  sigg (constant in intermittency diffusion term):
      sigg=1.0
c Coefficient for fnu in diffusion term
      sigkmu = 1.0
c
c Constants for Set 2 (2-layer SST model only):
      cmuc2=0.09
      beta2=0.0828
      sigo2=0.856
      sigk2=1.00
      alp2 =beta2/cmuc2 - sigo2*vk*vk/sqrt(cmuc2)
c
c Set up some other needed parameters
      jd2=(jdim-1)/2
      re=reue/xmach
      c2b=cbar/tinf
      c2bp=c2b+1.0
c
      iwrite=0
c
      if (icyc .le. nfreeze) then
        if (isklton .gt. 0) then
           nss=min(ncyc,nfreeze)
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),*)
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'('' turbulence model is frozen '',
     +     ''for '',i5,'' iterations or subits'')') nss
        end if
        sumn1 = 0.
        sumn2 = 0.
        sumn3 = 0.
        negn1 = 0
        negn2 = 0
        negn3 = 0
        return
      end if
      phi=0.
      if (real(dt) .gt. 0.) then
        if (abs(ita) .eq. 2) then
          phi=0.5
        else
          phi=0.
        end if
c   revert to old way (always 1st order for turb model) if itaturb=0
        if (itaturb .eq. 0) then
          phi=0.
          if (isklton .gt. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''   turb model is 1st'',
     +       '' order in time'')')
          end if
        else
          if (isklton .gt. 0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''   turb model is same'',
     +       '' order in time as mean flow eqns'')')
          end if
        end if
      end if
c
c Get laminar viscosity at cell centers
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            tt=gamma*q(j,k,i,5)/q(j,k,i,1)
            fnu(j,k,i)=c2bp*tt*sqrt(tt)/(c2b+tt)
          enddo
        enddo
      enddo
      do i=1,idim-1
        do k=1,kdim-1
          tt=gamma*qj0(k,i,5,1)/qj0(k,i,1,1)
          fnu(0,k,i)=c2bp*tt*sqrt(tt)/(c2b+tt)
          tt=gamma*qj0(k,i,5,3)/qj0(k,i,1,3)
          fnu(jdim,k,i)=c2bp*tt*sqrt(tt)/(c2b+tt)
        enddo
      enddo
      do i=1,idim-1
        do j=1,jdim-1
          tt=gamma*qk0(j,i,5,1)/qk0(j,i,1,1)
          fnu(j,0,i)=c2bp*tt*sqrt(tt)/(c2b+tt)
          tt=gamma*qk0(j,i,5,3)/qk0(j,i,1,3)
          fnu(j,kdim,i)=c2bp*tt*sqrt(tt)/(c2b+tt)
        enddo
      enddo
      if (i2d .ne. 1 .and. iaxi2planeturb .ne. 1) then
      do j=1,jdim-1
        do k=1,kdim-1
          tt=gamma*qi0(j,k,5,1)/qi0(j,k,1,1)
          fnu(j,k,0)=c2bp*tt*sqrt(tt)/(c2b+tt)
          tt=gamma*qi0(j,k,5,3)/qi0(j,k,1,3)
          fnu(j,k,idim)=c2bp*tt*sqrt(tt)/(c2b+tt)
        enddo
      enddo
      end if
c
c Load appropriate turb viscosity at cell centers (NOTE:  the code still
c uses vi0,vj0,vk0 values for turb viscosity at ghost cells - this means
c that for nonlinear models, anutp and anutm have variable cmu in them  
c at the block edges)
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            v3dtmp(j,k,i)=vist3d(j,k,i)
          enddo
        enddo
      enddo
c   Load appropriate vist3d value into ghost cells
      do i=1,idim-1
      do j=1,jdim-1
        v3dtmp(j,0,i)=bck(j,i,1)*(iwf(3)*v3dtmp(j,1,i) +
     +    (1-iwf(3))*2.*vk0(j,i,1,1)-v3dtmp(j,1,i))+
     +    (1.-bck(j,i,1))*vk0(j,i,1,1)
        v3dtmp(j,kdim,i)=bck(j,i,2)*(iwf(3)*v3dtmp(j,kdim-1,i) +
     +    (1-iwf(3))*2.*vk0(j,i,1,3)-v3dtmp(j,kdim-1,i))+
     +    (1.-bck(j,i,2))*vk0(j,i,1,3)
      enddo
      enddo
      do i=1,idim-1
      do k=1,kdim-1
        v3dtmp(0,k,i)=bcj(k,i,1)*(iwf(2)*v3dtmp(1,k,i) +
     +    (1-iwf(2))*2.*vj0(k,i,1,1)-v3dtmp(1,k,i))+
     +    (1.-bcj(k,i,1))*vj0(k,i,1,1)
        v3dtmp(jdim,k,i)=bcj(k,i,2)*(iwf(2)*v3dtmp(jdim-1,k,i) +
     +    (1-iwf(2))*2.*vj0(k,i,1,3)-v3dtmp(jdim-1,k,i))+
     +    (1.-bcj(k,i,2))*vj0(k,i,1,3)
      enddo
      enddo
      if (i2d .ne. 1 .and. iaxi2planeturb .ne. 1) then
      do j=1,jdim-1
      do k=1,kdim-1
        v3dtmp(j,k,0)=bci(j,k,1)*(iwf(1)*v3dtmp(j,k,1) +
     +    (1-iwf(1))*2.*vi0(j,k,1,1)-v3dtmp(j,k,1))+
     +    (1.-bci(j,k,1))*vi0(j,k,1,1)
        v3dtmp(j,k,idim)=bci(j,k,2)*(iwf(1)*v3dtmp(j,k,idim-1) +
     +    (1-iwf(1))*2.*vi0(j,k,1,3)-v3dtmp(j,k,idim-1))+
     +    (1.-bci(j,k,2))*vi0(j,k,1,3)
      enddo
      enddo
      end if
c If this is 1st global subiteration for time-accurate computation,
c save zksav (at time step n):
c zksav2(j,k,i,1), (j,k,i,2), and (j,k,i,3) are turb quantities
c zksav2(j,k,i,4), (j,k,i,5), and (j,k,i,6) are Delta turb quantities
      if (real(dt) .gt. 0. .and. icyc .eq. 1) then
      if (abs(ita) .eq. 2) then
c     if zksav2 at 1st point is zero, then we know that we do not have
c     2nd order data from the restart; no choice but to set
c     zksav2(j,k,i,4&5&6)=deltaQ=0 for 1st iteration
        if (real(zksav2(1,1,1,1)) .eq. 0.) then
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              zksav2(j,k,i,4)=0.
              zksav2(j,k,i,5)=0.
              zksav2(j,k,i,6)=0.
            enddo
          enddo
        enddo
        else
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              zksav2(j,k,i,4)=zksav(j,k,i,1)-zksav2(j,k,i,1)
              zksav2(j,k,i,5)=zksav(j,k,i,2)-zksav2(j,k,i,2)
              zksav2(j,k,i,6)=zksav(j,k,i,3)-zksav2(j,k,i,3)
            enddo
          enddo
        enddo
        end if
      end if
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              zksav2(j,k,i,1)=zksav(j,k,i,1)
              zksav2(j,k,i,2)=zksav(j,k,i,2)
              zksav2(j,k,i,3)=zksav(j,k,i,3)
            enddo
          enddo
        enddo
      end if
c Get TURRE values
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,i,1)=zksav(j,k,i,1)
            turre(j,k,i,2)=zksav(j,k,i,2)
            turre(j,k,i,3)=zksav(j,k,i,3)
          enddo
        enddo
      enddo
c
c Iterate to solve the equations
      do 500 not=1,nsubit
c
c    set up boundary conditions (they are in ghost cells everywhere)
c      (1) k=0 boundary:
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,0,i,1)=tk0(j,i,1,1)
            turre(j,0,i,2)=tk0(j,i,2,1)
            turre(j,0,i,3)=tk0(j,i,3,1)
          enddo
        enddo
c      (2) k=kdim boundary:
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,kdim,i,1)=tk0(j,i,1,3)
            turre(j,kdim,i,2)=tk0(j,i,2,3)
            turre(j,kdim,i,3)=tk0(j,i,3,3)
          enddo
        enddo
c      (3) j=0 boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(0,k,i,1)=tj0(k,i,1,1)
            turre(0,k,i,2)=tj0(k,i,2,1)
            turre(0,k,i,3)=tj0(k,i,3,1)
          enddo
        enddo
c      (4) j=jdim boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(jdim,k,i,1)=tj0(k,i,1,3)
            turre(jdim,k,i,2)=tj0(k,i,2,3)
            turre(jdim,k,i,3)=tj0(k,i,3,3)
          enddo
        enddo
        if (i2d .ne. 1 .and. iaxi2planeturb .ne. 1) then
c      (5) i=0 boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,0,1)=ti0(j,k,1,1)
            turre(j,k,0,2)=ti0(j,k,2,1)
            turre(j,k,0,3)=ti0(j,k,3,1)
          enddo
        enddo
c      (6) i=idim boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,idim,1)=ti0(j,k,1,3)
            turre(j,k,idim,2)=ti0(j,k,2,3)
            turre(j,k,idim,3)=ti0(j,k,3,3)
          enddo
        enddo
        end if
        if (iturbord .ne. 1) then
c      (1) k=0 boundary:
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,-1,i,1)=tk0(j,i,1,2)
            turre(j,-1,i,2)=tk0(j,i,2,2)
            turre(j,-1,i,3)=tk0(j,i,3,2)
          enddo
        enddo
c      (2) k=kdim boundary:
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,kdim+1,i,1)=tk0(j,i,1,4)
            turre(j,kdim+1,i,2)=tk0(j,i,2,4)
            turre(j,kdim+1,i,3)=tk0(j,i,3,4)
          enddo
        enddo
c      (3) j=0 boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(-1,k,i,1)=tj0(k,i,1,2)
            turre(-1,k,i,2)=tj0(k,i,2,2)
            turre(-1,k,i,3)=tj0(k,i,3,2)
          enddo
        enddo
c      (4) j=jdim boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(jdim+1,k,i,1)=tj0(k,i,1,4)
            turre(jdim+1,k,i,2)=tj0(k,i,2,4)
            turre(jdim+1,k,i,3)=tj0(k,i,3,4)
          enddo
        enddo
        if (i2d .ne. 1 .and. iaxi2planeturb .ne. 1) then
c      (5) i=0 boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,-1,1)=ti0(j,k,1,2)
            turre(j,k,-1,2)=ti0(j,k,2,2)
            turre(j,k,-1,3)=ti0(j,k,3,2)
          enddo
        enddo
c      (6) i=idim boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,idim+1,1)=ti0(j,k,1,4)
            turre(j,k,idim+1,2)=ti0(j,k,2,4)
            turre(j,k,idim+1,3)=ti0(j,k,3,4)
          enddo
        enddo
        end if
        end if
c    **************
c
c Get damp1 = CD = cross derivative term for SST:
        if(ivmx .eq. 30 .or. ivmx.eq.31) then
        do k=1,kdim-1
          do i=1,idim-1
            do j=1,jdim-1
              xa=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4) +
     +           sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4) +
     +           sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i)
              za=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4) +
     +           sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i)
              tt=xa*xa+ya*ya+za*za
              ca=2.*sigo2*tt/(turre(j,k,i,1)*re)
              damp1(j,k,i)=0.25*ca*(turre(j,k+1,i,1)-turre(j,k-1,i,1))*
     +                             (turre(j,k+1,i,2)-turre(j,k-1,i,2))
            enddo
          enddo
        enddo
        do j=1,jdim-1
          do i=1,idim-1
            do k=1,kdim-1
              xa=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4) +
     +           sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4) +
     +           sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i)
              za=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4) +
     +           sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i)
              tt=xa*xa+ya*ya+za*za
              ca=2.*sigo2*tt/(turre(j,k,i,1)*re)
              damp1(j,k,i)=damp1(j,k,i)+
     +                0.25*ca*(turre(j+1,k,i,1)-turre(j-1,k,i,1))*
     +                        (turre(j+1,k,i,2)-turre(j-1,k,i,2))
            enddo
          enddo
        enddo
        if(i2d .ne. 1 .and. iaxi2planeturb .ne. 1) then
          do i=1,idim-1
            do k=1,kdim-1
              do j=1,jdim-1
                xa=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4) +
     +             si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i)
                ya=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4) +
     +             si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i)
                za=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4) +
     +             si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i)
                tt=xa*xa+ya*ya+za*za
                ca=2.*sigo2*tt/(turre(j,k,i,1)*re)
                damp1(j,k,i)=damp1(j,k,i)+
     +                0.25*ca*(turre(j,k,i+1,1)-turre(j,k,i-1,1))*
     +                        (turre(j,k,i+1,2)-turre(j,k,i-1,2))
              enddo
            enddo
          enddo
        end if
        end if
c   get blend = F1 factor
        if (ivmx .eq. 30 .or. ivmx.eq.31) then
        do k=1,kdim-1
          do i=1,idim-1
            do j=1,jdim-1
              arg1=sqrt(turre(j,k,i,2))/(.09*re*turre(j,k,i,1)*
     +          ccabs(smin(j,k,i)))
              arg2=500.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*re*re*
     +          smin(j,k,i)*turre(j,k,i,1))
              arga=ccmax(arg1,arg2)
              temp=ccmaxcr(damp1(j,k,i)*re,1.e-20)
              argb=4.*sigo2*turre(j,k,i,2)/(temp*smin(j,k,i)*
     +             smin(j,k,i))
              arg=ccmin(arga,argb)
              blend(j,k,i)=cctanh(arg*arg*arg*arg)
              if (itrans_on .eq. 1) then
                re_y=q(j,k,i,1)*ccabs(smin(j,k,i))*sqrt(turre(j,k,i,2))/
     +               fnu(j,k,i)*re
                fff3=exp(-((re_y/120.)**8))
                blend(j,k,i)=ccmax(blend(j,k,i),fff3)
              end if
            enddo
          enddo
        enddo
        end if
c
c    F_eta_eta viscous terms
c       Interior points
        do k=2,kdim-2
          kl=k-1
          ku=k+1
          do i=1,idim-1
            do j=1,jdim-1
              dfacep=0.5*(blend(j,k,i)+blend(j,ku,i))
              dfacem=0.5*(blend(j,k,i)+blend(j,kl,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volku=vol(j,ku,i)
              xp=sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              yp=sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              zp=sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              volkl=vol(j,kl,i)
              xm=sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              ym=sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              zm=sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              xa=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i)
              za=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,1)=-byy*turre(j,k-1,i,1)
     +          -cyy*turre(j,k,i,1) -dyy*turre(j,k+1,i,1)
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,2)=-byy*turre(j,k-1,i,2)
     +          -cyy*turre(j,k,i,2) -dyy*turre(j,k+1,i,2)
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,3)=-byy*turre(j,k-1,i,3)
     +          -cyy*turre(j,k,i,3) -dyy*turre(j,k+1,i,3)
            enddo
          enddo
        enddo
c
c       K0 boundary points
          k=1
          kl=1
          ku=min(2,kdim-1)
          do i=1,idim-1
            do j=1,jdim-1
              dfacep=0.5*(blend(j,k,i)+blend(j,ku,i))
              dfacem=0.5*(blend(j,k,i)+blend(j,kl,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volku=vol(j,ku,i)
              xp=sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              yp=sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              zp=sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              volkl=volk0(j,i,1)
              xm=sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              ym=sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              zm=sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              xa=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i)
              za=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,1)=-byy*turre(j,k-1,i,1)
     +          -cyy*turre(j,k,i,1) -dyy*turre(j,k+1,i,1)
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,2)=-byy*turre(j,k-1,i,2)
     +          -cyy*turre(j,k,i,2) -dyy*turre(j,k+1,i,2)
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,3)=-byy*turre(j,k-1,i,3)
     +          -cyy*turre(j,k,i,3) -dyy*turre(j,k+1,i,3)
            enddo
          enddo
c
c       KDIM points
          k=kdim-1
          kl=kdim-2
          ku=kdim-1
          do i=1,idim-1
            do j=1,jdim-1
              dfacep=0.5*(blend(j,k,i)+blend(j,ku,i))
              dfacem=0.5*(blend(j,k,i)+blend(j,kl,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volku=volk0(j,i,3)
              xp=sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              yp=sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              zp=sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              volkl=vol(j,kl,i)
              xm=sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              ym=sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              zm=sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              xa=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i)
              za=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,1)=-byy*turre(j,k-1,i,1)
     +          -cyy*turre(j,k,i,1) -dyy*turre(j,k+1,i,1)
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,2)=-byy*turre(j,k-1,i,2)
     +          -cyy*turre(j,k,i,2) -dyy*turre(j,k+1,i,2)
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,3)=-byy*turre(j,k-1,i,3)
     +          -cyy*turre(j,k,i,3) -dyy*turre(j,k+1,i,3)
            enddo
          enddo
c    Advective terms in eta
        if (iturbord .eq. 1) then
        do k=1,kdim-1
          do i=1,idim-1
            do j=1,jdim-1
              xc=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,1)*sk(j,k,i  ,4))/vol(j,k,i)
              yc=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,2)*sk(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,3)*sk(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,5)*sk(j,k,i  ,4))/vol(j,k,i)
              uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
              rhside(j,k,i,1)=rhside(j,k,i,1)-uu*(app*(turre(j,k,i,1)-
     +          turre(j,k-1,i,1)) + apm*(turre(j,k+1,i,1)-
     +          turre(j,k,i,1)))
              rhside(j,k,i,2)=rhside(j,k,i,2)-uu*(app*(turre(j,k,i,2)-
     +          turre(j,k-1,i,2)) + apm*(turre(j,k+1,i,2)-
     +          turre(j,k,i,2)))
              rhside(j,k,i,3)=rhside(j,k,i,3)-uu*(app*(turre(j,k,i,3)-
     +          turre(j,k-1,i,3)) + apm*(turre(j,k+1,i,3)-
     +          turre(j,k,i,3)))
            enddo
          enddo
        enddo
        else
c       2nd order upwind; LHS remains 1st order everywhere
        do k=1,kdim-1
          do i=1,idim-1
            do j=1,jdim-1
              xc=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,1)*sk(j,k,i  ,4))/vol(j,k,i)
              yc=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,2)*sk(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,3)*sk(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,5)*sk(j,k,i  ,4))/vol(j,k,i)
              uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
             rhside(j,k,i,1)=rhside(j,k,i,1)-0.5*uu*app*turre(j,k-2,i,1)
     +                                    +2.*uu*app*turre(j,k-1,i,1)
     +                                   -1.5*uu*app*turre(j,k,i,1)
     +                                   +1.5*uu*apm*turre(j,k,i,1)
     +                                    -2.*uu*apm*turre(j,k+1,i,1)
     +                                   +0.5*uu*apm*turre(j,k+2,i,1)
             rhside(j,k,i,2)=rhside(j,k,i,2)-0.5*uu*app*turre(j,k-2,i,2)
     +                                    +2.*uu*app*turre(j,k-1,i,2)
     +                                   -1.5*uu*app*turre(j,k,i,2)
     +                                   +1.5*uu*apm*turre(j,k,i,2)
     +                                    -2.*uu*apm*turre(j,k+1,i,2)
     +                                   +0.5*uu*apm*turre(j,k+2,i,2)
             rhside(j,k,i,3)=rhside(j,k,i,3)-0.5*uu*app*turre(j,k-2,i,3)
     +                                    +2.*uu*app*turre(j,k-1,i,3)
     +                                   -1.5*uu*app*turre(j,k,i,3)
     +                                   +1.5*uu*apm*turre(j,k,i,3)
     +                                    -2.*uu*apm*turre(j,k+1,i,3)
     +                                   +0.5*uu*apm*turre(j,k+2,i,3)
            enddo
          enddo
        enddo
        end if
c
c    F_xi_xi viscous terms
c       Interior points
        do j=2,jdim-2
          jl=j-1
          ju=j+1
          do i=1,idim-1
            do k=1,kdim-1
              dfacep=0.5*(blend(j,k,i)+blend(ju,k,i))
              dfacem=0.5*(blend(j,k,i)+blend(jl,k,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volju=vol(ju,k,i)
              xp=sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              yp=sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              zp=sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              voljl=vol(jl,k,i)
              xm=sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              ym=sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              zm=sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              xa=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i)
              za=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,1)=rhside(j,k,i,1)-bxx*turre(j-1,k,i,1)
     +          -cxx*turre(j,k,i,1) -dxx*turre(j+1,k,i,1)
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,2)=rhside(j,k,i,2)-bxx*turre(j-1,k,i,2)
     +          -cxx*turre(j,k,i,2) -dxx*turre(j+1,k,i,2)
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,3)=rhside(j,k,i,3)-bxx*turre(j-1,k,i,3)
     +          -cxx*turre(j,k,i,3) -dxx*turre(j+1,k,i,3)
            enddo
          enddo
        enddo
c
c       J0 boundary points
          j=1
          jl=1
          ju=min(2,jdim-1)
          do i=1,idim-1
            do k=1,kdim-1
              dfacep=0.5*(blend(j,k,i)+blend(ju,k,i))
              dfacem=0.5*(blend(j,k,i)+blend(jl,k,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volju=vol(ju,k,i)
              xp=sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              yp=sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              zp=sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              voljl=volj0(k,i,1)
              xm=sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              ym=sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              zm=sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              xa=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i)
              za=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,1)=rhside(j,k,i,1)-bxx*turre(j-1,k,i,1)
     +          -cxx*turre(j,k,i,1) -dxx*turre(j+1,k,i,1)
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,2)=rhside(j,k,i,2)-bxx*turre(j-1,k,i,2)
     +          -cxx*turre(j,k,i,2) -dxx*turre(j+1,k,i,2)
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,3)=rhside(j,k,i,3)-bxx*turre(j-1,k,i,3)
     +          -cxx*turre(j,k,i,3) -dxx*turre(j+1,k,i,3)
            enddo
          enddo
c
c       JDIM boundary points
          j=jdim-1
          jl=jdim-2
          ju=jdim-1
          do i=1,idim-1
            do k=1,kdim-1
              dfacep=0.5*(blend(j,k,i)+blend(ju,k,i))
              dfacem=0.5*(blend(j,k,i)+blend(jl,k,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volju=volj0(k,i,3)
              xp=sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              yp=sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              zp=sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              voljl=vol(jl,k,i)
              xm=sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              ym=sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              zm=sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              xa=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i)
              za=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,1)=rhside(j,k,i,1)-bxx*turre(j-1,k,i,1)
     +          -cxx*turre(j,k,i,1) -dxx*turre(j+1,k,i,1)
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,2)=rhside(j,k,i,2)-bxx*turre(j-1,k,i,2)
     +          -cxx*turre(j,k,i,2) -dxx*turre(j+1,k,i,2)
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,3)=rhside(j,k,i,3)-bxx*turre(j-1,k,i,3)
     +          -cxx*turre(j,k,i,3) -dxx*turre(j+1,k,i,3)
            enddo
          enddo
c    Advective terms in xi
        if (iturbord .eq. 1) then
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              xc=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,1)*sj(j,k,i  ,4))/vol(j,k,i)
              yc=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,2)*sj(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,3)*sj(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,5)*sj(j,k,i  ,4))/vol(j,k,i)
              uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
              rhside(j,k,i,1)=rhside(j,k,i,1)-uu*(app*(turre(j,k,i,1)-
     +          turre(j-1,k,i,1)) + apm*(turre(j+1,k,i,1)-
     +          turre(j,k,i,1)))
              rhside(j,k,i,2)=rhside(j,k,i,2)-uu*(app*(turre(j,k,i,2)-
     +          turre(j-1,k,i,2)) + apm*(turre(j+1,k,i,2)-
     +          turre(j,k,i,2)))
              rhside(j,k,i,3)=rhside(j,k,i,3)-uu*(app*(turre(j,k,i,3)-
     +          turre(j-1,k,i,3)) + apm*(turre(j+1,k,i,3)-
     +          turre(j,k,i,3)))
            enddo
          enddo
        enddo
        else
c       2nd order upwind; LHS remains 1st order everywhere
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              xc=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,1)*sj(j,k,i  ,4))/vol(j,k,i)
              yc=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,2)*sj(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,3)*sj(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,5)*sj(j,k,i  ,4))/vol(j,k,i)
              uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
             rhside(j,k,i,1)=rhside(j,k,i,1)-0.5*uu*app*turre(j-2,k,i,1)
     +                                    +2.*uu*app*turre(j-1,k,i,1)
     +                                   -1.5*uu*app*turre(j,k,i,1)
     +                                   +1.5*uu*apm*turre(j,k,i,1)
     +                                    -2.*uu*apm*turre(j+1,k,i,1)
     +                                   +0.5*uu*apm*turre(j+2,k,i,1)
             rhside(j,k,i,2)=rhside(j,k,i,2)-0.5*uu*app*turre(j-2,k,i,2)
     +                                    +2.*uu*app*turre(j-1,k,i,2)
     +                                   -1.5*uu*app*turre(j,k,i,2)
     +                                   +1.5*uu*apm*turre(j,k,i,2)
     +                                    -2.*uu*apm*turre(j+1,k,i,2)
     +                                   +0.5*uu*apm*turre(j+2,k,i,2)
             rhside(j,k,i,3)=rhside(j,k,i,3)-0.5*uu*app*turre(j-2,k,i,3)
     +                                    +2.*uu*app*turre(j-1,k,i,3)
     +                                   -1.5*uu*app*turre(j,k,i,3)
     +                                   +1.5*uu*apm*turre(j,k,i,3)
     +                                    -2.*uu*apm*turre(j+1,k,i,3)
     +                                   +0.5*uu*apm*turre(j+2,k,i,3)
            enddo
          enddo
        enddo
        end if
c
c    F_zeta_zeta viscous terms
        if(i2d .ne. 1 .and. iaxi2planeturb .ne. 1) then
c         Interior points
          do i=2,idim-2
            il=i-1
            iu=i+1
            do k=1,kdim-1
              do j=1,jdim-1
                dfacep=0.5*(blend(j,k,i)+blend(j,k,iu))
                dfacem=0.5*(blend(j,k,i)+blend(j,k,il))
                sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
                sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
                sigop=dfacep*sigo1+(1.-dfacep)*sigo2
                sigom=dfacem*sigo1+(1.-dfacem)*sigo2
                voliu=vol(j,k,iu)
                xp=si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                yp=si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                zp=si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                volil=vol(j,k,il)
                xm=si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                ym=si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                zm=si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                xa=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4) + 
     +             si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i)
                ya=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4) + 
     +             si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i)
                za=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4) + 
     +             si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i)
c
                ttpo=xp*xa+yp*ya+zp*za
                ttmo=xm*xa+ym*ya+zm*za
                ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/
     +               vol(j,k,i)
                ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/
     +               vol(j,k,i)
c               choose between weak (o) and strong (n) conservation form
                ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
                ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
                anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1))
                anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,1)=rhside(j,k,i,1)-bzz*turre(j,k,i-1,1)
     +            -czz*turre(j,k,i,1) -dzz*turre(j,k,i+1,1)
                cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,2)=rhside(j,k,i,2)-bzz*turre(j,k,i-1,2)
     +            -czz*turre(j,k,i,2) -dzz*turre(j,k,i+1,2)
                cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,3)=rhside(j,k,i,3)-bzz*turre(j,k,i-1,3)
     +            -czz*turre(j,k,i,3) -dzz*turre(j,k,i+1,3)
              enddo
            enddo
          enddo
c
c         I0 boundary points
            i=1
            il=1
            iu=min(2,idim-1)
            do k=1,kdim-1
              do j=1,jdim-1
                dfacep=0.5*(blend(j,k,i)+blend(j,k,iu))
                dfacem=0.5*(blend(j,k,i)+blend(j,k,il))
                sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
                sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
                sigop=dfacep*sigo1+(1.-dfacep)*sigo2
                sigom=dfacem*sigo1+(1.-dfacem)*sigo2
                voliu=vol(j,k,iu)
                xp=si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                yp=si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                zp=si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                volil=voli0(j,k,1)
                xm=si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                ym=si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                zm=si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                xa=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4) + 
     +             si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i)
                ya=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4) + 
     +             si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i)
                za=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4) + 
     +             si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i)
c
                ttpo=xp*xa+yp*ya+zp*za
                ttmo=xm*xa+ym*ya+zm*za
                ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/
     +               vol(j,k,i)
                ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/
     +               vol(j,k,i)
c               choose between weak (o) and strong (n) conservation form
                ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
                ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
                anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1))
                anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,1)=rhside(j,k,i,1)-bzz*turre(j,k,i-1,1)
     +            -czz*turre(j,k,i,1) -dzz*turre(j,k,i+1,1)
                cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,2)=rhside(j,k,i,2)-bzz*turre(j,k,i-1,2)
     +            -czz*turre(j,k,i,2) -dzz*turre(j,k,i+1,2)
                cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,3)=rhside(j,k,i,3)-bzz*turre(j,k,i-1,3)
     +            -czz*turre(j,k,i,3) -dzz*turre(j,k,i+1,3)
              enddo
            enddo
c
c         IDIM boundary points 
            i=idim-1
            il=idim-2
            iu=idim-1
            do k=1,kdim-1
              do j=1,jdim-1
                dfacep=0.5*(blend(j,k,i)+blend(j,k,iu))
                dfacem=0.5*(blend(j,k,i)+blend(j,k,il))
                sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
                sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
                sigop=dfacep*sigo1+(1.-dfacep)*sigo2
                sigom=dfacem*sigo1+(1.-dfacem)*sigo2
                voliu=voli0(j,k,3)
                xp=si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                yp=si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                zp=si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                volil=vol(j,k,il)
                xm=si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                ym=si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                zm=si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                xa=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4) + 
     +             si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i)
                ya=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4) + 
     +             si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i)
                za=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4) + 
     +             si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i)
c
                ttpo=xp*xa+yp*ya+zp*za
                ttmo=xm*xa+ym*ya+zm*za
                ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/
     +               vol(j,k,i)
                ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/
     +               vol(j,k,i)
c               choose between weak (o) and strong (n) conservation form
                ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
                ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
                anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1))
                anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,1)=rhside(j,k,i,1)-bzz*turre(j,k,i-1,1)
     +            -czz*turre(j,k,i,1) -dzz*turre(j,k,i+1,1)
                cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,2)=rhside(j,k,i,2)-bzz*turre(j,k,i-1,2)
     +            -czz*turre(j,k,i,2) -dzz*turre(j,k,i+1,2)
                cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,3)=rhside(j,k,i,3)-bzz*turre(j,k,i-1,3)
     +            -czz*turre(j,k,i,3) -dzz*turre(j,k,i+1,3)
              enddo
            enddo
c    Advective terms in zeta
          if (iturbord .eq. 1) then
          do i=1,idim-1
            do k=1,kdim-1
              do j=1,jdim-1
                xc=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,1)*si(j,k,i  ,4))/vol(j,k,i)
                yc=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,2)*si(j,k,i  ,4))/vol(j,k,i)
                zc=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,3)*si(j,k,i  ,4))/vol(j,k,i)
                tc=0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,5)*si(j,k,i  ,4))/vol(j,k,i)
                uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
                sgnu=ccsignrc(1.,uu)
                app=0.5*(1.+sgnu)
                apm=0.5*(1.-sgnu)
                rhside(j,k,i,1)=rhside(j,k,i,1)-uu*(app*(turre(j,k,i,1)-
     +           turre(j,k,i-1,1)) + apm*(turre(j,k,i+1,1)-
     +           turre(j,k,i,1)))
                rhside(j,k,i,2)=rhside(j,k,i,2)-uu*(app*(turre(j,k,i,2)-
     +           turre(j,k,i-1,2)) + apm*(turre(j,k,i+1,2)-
     +           turre(j,k,i,2)))
                rhside(j,k,i,3)=rhside(j,k,i,3)-uu*(app*(turre(j,k,i,3)-
     +           turre(j,k,i-1,3)) + apm*(turre(j,k,i+1,3)-
     +           turre(j,k,i,3)))
              enddo
            enddo
          enddo
          else
c       2nd order upwind; LHS remains 1st order everywhere
          do i=1,idim-1
            do k=1,kdim-1
              do j=1,jdim-1
                xc=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,1)*si(j,k,i  ,4))/vol(j,k,i)
                yc=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,2)*si(j,k,i  ,4))/vol(j,k,i)
                zc=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,3)*si(j,k,i  ,4))/vol(j,k,i)
                tc=0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,5)*si(j,k,i  ,4))/vol(j,k,i)
                uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
             rhside(j,k,i,1)=rhside(j,k,i,1)-0.5*uu*app*turre(j,k,i-2,1)
     +                                    +2.*uu*app*turre(j,k,i-1,1)
     +                                   -1.5*uu*app*turre(j,k,i,1)
     +                                   +1.5*uu*apm*turre(j,k,i,1)
     +                                    -2.*uu*apm*turre(j,k,i+1,1)
     +                                   +0.5*uu*apm*turre(j,k,i+2,1)
             rhside(j,k,i,2)=rhside(j,k,i,2)-0.5*uu*app*turre(j,k,i-2,2)
     +                                    +2.*uu*app*turre(j,k,i-1,2)
     +                                   -1.5*uu*app*turre(j,k,i,2)
     +                                   +1.5*uu*apm*turre(j,k,i,2)
     +                                    -2.*uu*apm*turre(j,k,i+1,2)
     +                                   +0.5*uu*apm*turre(j,k,i+2,2)
             rhside(j,k,i,3)=rhside(j,k,i,3)-0.5*uu*app*turre(j,k,i-2,3)
     +                                    +2.*uu*app*turre(j,k,i-1,3)
     +                                   -1.5*uu*app*turre(j,k,i,3)
     +                                   +1.5*uu*apm*turre(j,k,i,3)
     +                                    -2.*uu*apm*turre(j,k,i+1,3)
     +                                   +0.5*uu*apm*turre(j,k,i+2,3)
              enddo
            enddo
          enddo
          end if
        end if
c
c    * * * * * * * * * * * *
c    ADD SOURCE TERMS TO RHS
c    * * * * * * * * * * * *
      if(isklton .gt. 0) then
        if(ilamlo.eq.0 .or. jlamlo.eq.0 .or. klamlo.eq.0) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'('' Block #'',i5,'' in 3-eqn turb model'',
     .   '' has no laminar regions'')') nbl
        else
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'('' Block #'',i5,'' in 3-eqn turb model'',
     .   '' - laminar region is:'')') nbl
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'('' i='',i5,'' to'',i5,'', j='',i5,'' to'',
     .   i5,'', k='',i5,'' to'',i5)') ilamlo,ilamhi,jlamlo,jlamhi,
     .   klamlo,klamhi
        if (i_lam_forcezero .eq. 1) then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''    ...forcing vist3d=0'')')
        end if
        end if
      end if
c  Source terms:
      if (ivmx .eq. 30) then
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
              if ((i.ge.ilamlo .and. i.lt.ilamhi .and.
     .             j.ge.jlamlo .and. j.lt.jlamhi .and.
     .             k.ge.klamlo .and. k.lt.klamhi) .or.
     .             real(smin(j,k,i)) .lt. 0.) then
                cutoff=0.
              else
                cutoff=1.
              end if
              betax=blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2
              cmuc= blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2
              alp=  blend(j,k,i)*alp1 +(1.-blend(j,k,i))*alp2
c            Add to RHS:
c            Determine Sij values:
              s11 = ux(j,k,i,1)
              s22 = ux(j,k,i,5)
              s33 = ux(j,k,i,9)
              s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
              s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
              s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
              tracepart=(s11+s22+s33)*prod2d3dtrace
              s11t=s11-tracepart
              s22t=s22-tracepart
              s33t=s33-tracepart
              xis = s11t*s11t + s22t*s22t + s33t*s33t +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              if (ikoprod .eq. 1) then
              rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp/re*
     +         2.*xis - re*betax*turre(j,k,i,1)**2
     +         +(1.-blend(j,k,i))*damp1(j,k,i)
              pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*2.*xis
              else
              rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp/re*
     +         vor(j,k,i)**2 - re*betax*turre(j,k,i,1)**2
     +         +(1.-blend(j,k,i))*damp1(j,k,i)
              pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)**2
              end if
              dk=re*cmuc*turre(j,k,i,1)*turre(j,k,i,2)
              pk=ccmin(pk,(pklimterm*dk))
              dkfactor=ccmaxcr(turre(j,k,i,3),0.1)
              rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk*turre(j,k,i,3)-
     +         dk*ccmincr(dkfactor,1.0)
              rhside(j,k,i,1)=rhside(j,k,i,1) + keepambient*
     +         re*betax*tur10(1)*tur10(1)
              rhside(j,k,i,2)=rhside(j,k,i,2) + keepambient*
     +         re*cmuc*tur10(1)*tur10(2)
c
              re_nu=q(j,k,i,1)*smin(j,k,i)**2*vor(j,k,i)*re/fnu(j,k,i)
              re_t=q(j,k,i,1)*turre(j,k,i,2)/(fnu(j,k,i)*turre(j,k,i,1))
              fonset1=re_nu/(2.193*re_thetac)
              fonset2=ccmax(fonset1,fonset1**4)
              fonset2=ccmincr(fonset2,2.0)
              ftrnlength=cctanh(re_nu/(2.193*0.35*300.))
              trnlength=(0.014*re_thetac+1.25)*ftrnlength+
     +                  2.*(1.-ftrnlength)
              fonset3=ccmaxcr((1.-re_t/trnlength),0.)
              fonset=ccmaxcr(fonset2-fonset3,0.)
              fturb=exp(-((re_t/2.)**2))
              rhside(j,k,i,3)=rhside(j,k,i,3)+
     +         5.*sqrt(2.*xis)*fonset*(1.0-ccmincr(turre(j,k,i,3),1.0))-
     +         0.1*vor(j,k,i)*fturb*cctanh(turre(j,k,i,3)/0.35)
            enddo
          enddo
        enddo
      elseif (ivmx .eq. 31) then
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
              if ((i.ge.ilamlo .and. i.lt.ilamhi .and.
     .             j.ge.jlamlo .and. j.lt.jlamhi .and.
     .             k.ge.klamlo .and. k.lt.klamhi) .or.
     .             real(smin(j,k,i)) .lt. 0.) then
                cutoff=0.
              else
                cutoff=1.
              end if
              betax=blend(j,k,i)*beta1+(1.-blend(j,k,i))*beta2
              cmuc= blend(j,k,i)*cmuc1+(1.-blend(j,k,i))*cmuc2
              alp=  blend(j,k,i)*alp1 +(1.-blend(j,k,i))*alp2
c            for now, re_thetat_tilde taken = re_thetat
c            (may instead need to solve 4th eqn eventually)
              re_thetat_tilde=re_thetat
              re_nu=q(j,k,i,1)*smin(j,k,i)**2*vor(j,k,i)*re/fnu(j,k,i)
              re_t=q(j,k,i,1)*turre(j,k,i,2)/(fnu(j,k,i)*turre(j,k,i,1))
              uuu=sqrt(q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)**2)
              uuu=ccmaxcr(uuu,1.e-20)
              theta_bl=re_thetat_tilde*fnu(j,k,i)/(q(j,k,i,1)*uuu*re)
              delta_bl=15.*theta_bl/2.
              delta=50.*vor(j,k,i)*ccabs(smin(j,k,i))*delta_bl/uuu
              re_omega=q(j,k,i,1)*turre(j,k,i,1)*smin(j,k,i)**2/
     +          fnu(j,k,i)*re*re
              fwake=exp(-(re_omega/1.e5)**2)
              f_thetat1=fwake*exp(-(ccabs(smin(j,k,i))/delta)**4)
              gamma_lim = ccmincr(turre(j,k,i,3),1.0)
              f_thetat2=1.-((ce2*gamma_lim-1.)/(ce2-1.))**2
              f_thetat=ccmax(f_thetat1,f_thetat2)
              f_thetat=ccmincr(f_thetat,1.0)
              f_reattach=exp(-(re_t/15.)**4)
              fonset1=re_nu/(2.193*re_thetac)
              gamma_sep1=ccmaxcr(fonset1-1.,0.)
              gamma_sep=ccmincr(8.0*gamma_sep1*f_reattach,5.0)
              gamma_sep=gamma_sep*f_thetat
              gamma_eff=ccmax(turre(j,k,i,3),gamma_sep)
c            Determine Sij values:
              s11 = ux(j,k,i,1)
              s22 = ux(j,k,i,5)
              s33 = ux(j,k,i,9)
              s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
              s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
              s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
              tracepart=(s11+s22+s33)*prod2d3dtrace
              s11t=s11-tracepart
              s22t=s22-tracepart
              s33t=s33-tracepart
              xis = s11t*s11t + s22t*s22t + s33t*s33t +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              if (ikoprod .eq. 1) then
              rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp/re*
     +         2.*xis - re*betax*turre(j,k,i,1)**2
     +         +(1.-blend(j,k,i))*damp1(j,k,i)
              pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*2.*xis
              else
              rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp/re*
     +         vor(j,k,i)**2 - re*betax*turre(j,k,i,1)**2
     +         +(1.-blend(j,k,i))*damp1(j,k,i)
              pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)**2
              end if
              dk=re*cmuc*turre(j,k,i,1)*turre(j,k,i,2)
              pk=ccmin(pk,(pklimterm*dk))
              dkfactor=ccmaxcr(gamma_eff,0.1)
              rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk*gamma_eff-
     +         dk*ccmincr(dkfactor,1.0)
              rhside(j,k,i,1)=rhside(j,k,i,1) + keepambient*
     +         re*betax*tur10(1)*tur10(1)
              rhside(j,k,i,2)=rhside(j,k,i,2) + keepambient*
     +         re*cmuc*tur10(1)*tur10(2)
c
              fonset2=ccmax(fonset1,fonset1**4)
              fonset2=ccmincr(fonset2,2.0)
              fonset3=ccmaxcr((1.-(re_t/2.5)**3),0.)
              fonset=ccmaxcr(fonset2-fonset3,0.)
              fturb=exp(-((re_t/4.)**4))
              rhside(j,k,i,3)=rhside(j,k,i,3)+
     +         (flength*2.0*sqrt(2.*xis)*(gamma_lim*fonset)**0.5)*
     +         (1.0-gamma_lim) -
     +         (0.06*vor(j,k,i)*gamma_lim*fturb*
     +         (ce2*gamma_lim-1.0))
            enddo
          enddo
        enddo
      end if
c
c    Implicit F_eta_eta viscous terms.  Do over all i's
        do i=1,idim-1
c         Interior points
          do k=2,kdim-2
            kl=k-1
            ku=k+1
            do j=1,jdim-1
              dfacep=0.5*(blend(j,k,i)+blend(j,ku,i))
              dfacem=0.5*(blend(j,k,i)+blend(j,kl,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volku=vol(j,ku,i)
              xp=sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              yp=sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              zp=sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              volkl=vol(j,kl,i)
              xm=sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              ym=sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              zm=sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              xa=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i)
              za=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              by(j,k)=-cdm
              cy(j,k)= cdp+cdm
              dy(j,k)=-cdp
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              by2(j,k)=-cdm
              cy2(j,k)= cdp+cdm
              dy2(j,k)=-cdp
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              by3(j,k)=-cdm
              cy3(j,k)= cdp+cdm
              dy3(j,k)=-cdp
            enddo
            do j=1,jdim-1
              xc=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,1)*sk(j,k,i  ,4))/vol(j,k,i)
              yc=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,2)*sk(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,3)*sk(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,5)*sk(j,k,i  ,4))/vol(j,k,i)
              uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
              by(j,k)=by(j,k) - uu*app
              by2(j,k)=by2(j,k) - uu*app
              by3(j,k)=by3(j,k) - uu*app
              cy(j,k)=cy(j,k) + uu*(app-apm)
              cy2(j,k)=cy2(j,k) + uu*(app-apm)
              cy3(j,k)=cy3(j,k) + uu*(app-apm)
              dy(j,k)=dy(j,k) + uu*apm
              dy2(j,k)=dy2(j,k) + uu*apm
              dy3(j,k)=dy3(j,k) + uu*apm
            enddo
c      Add part of source terms to diagonal LHS in this sweep (helps increase
c      diagonal dominance)
            do j=1,jdim-1
              if (ivmx .eq. 30 .or. ivmx .eq. 31) then
                cmuc=  blend(j,k,i)*cmuc1 + (1.-blend(j,k,i))*cmuc2
                betax= blend(j,k,i)*beta1 + (1.-blend(j,k,i))*beta2
                cyadd= 2.*re*betax*turre(j,k,i,1) +
     +                 ccabs(damp1(j,k,i))/turre(j,k,i,1)
                cy2add=re*cmuc*turre(j,k,i,1)
c               cy3add=???...
              end if
              cy(j,k)=cy(j,k) + cyadd
              cy2(j,k)=cy2(j,k) + cy2add
c             cy3(j,k)=cy3(j,k) + cy3add
            enddo
            do j=1,jdim-1
              by(j,k)=by(j,k)*timestp(j,k,i)
              by2(j,k)=by2(j,k)*timestp(j,k,i)*factor2
              by3(j,k)=by3(j,k)*timestp(j,k,i)*factor3
              cy(j,k)=cy(j,k)*timestp(j,k,i)+1.0*(1.+phi)
              cy2(j,k)=cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.+phi)
              cy3(j,k)=cy3(j,k)*timestp(j,k,i)*factor3+1.0*(1.+phi)
              dy(j,k)=dy(j,k)*timestp(j,k,i)
              dy2(j,k)=dy2(j,k)*timestp(j,k,i)*factor2
              dy3(j,k)=dy3(j,k)*timestp(j,k,i)*factor3
              fy(j,k)=rhside(j,k,i,1)*timestp(j,k,i)
              fy2(j,k)=rhside(j,k,i,2)*timestp(j,k,i)*factor2
              fy3(j,k)=rhside(j,k,i,3)*timestp(j,k,i)*factor3
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-
     +                 turre(j,k,i,1))+phi*zksav2(j,k,i,4)
                fy2(j,k)=fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-
     +                 turre(j,k,i,2))+phi*zksav2(j,k,i,5)
                fy3(j,k)=fy3(j,k)+(1.+phi)*(zksav2(j,k,i,3)-
     +                 turre(j,k,i,3))+phi*zksav2(j,k,i,6)
              enddo
            end if
          enddo
c
c         K0 boundary points
            k=1
            kl=1
            ku=min(2,kdim-1)
            do j=1,jdim-1
              dfacep=0.5*(blend(j,k,i)+blend(j,ku,i))
              dfacem=0.5*(blend(j,k,i)+blend(j,kl,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volku=vol(j,ku,i)
              xp=sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              yp=sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              zp=sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              volkl=volk0(j,i,1)
              xm=sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              ym=sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              zm=sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              xa=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i)
              za=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              by(j,k)=-cdm
              cy(j,k)= cdp+cdm
              dy(j,k)=-cdp
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              by2(j,k)=-cdm
              cy2(j,k)= cdp+cdm
              dy2(j,k)=-cdp
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              by3(j,k)=-cdm
              cy3(j,k)= cdp+cdm
              dy3(j,k)=-cdp
            enddo
            do j=1,jdim-1
              xc=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,1)*sk(j,k,i  ,4))/vol(j,k,i)
              yc=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,2)*sk(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,3)*sk(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,5)*sk(j,k,i  ,4))/vol(j,k,i)
              uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
              by(j,k)=by(j,k) - uu*app
              by2(j,k)=by2(j,k) - uu*app
              by3(j,k)=by3(j,k) - uu*app
              cy(j,k)=cy(j,k) + uu*(app-apm)
              cy2(j,k)=cy2(j,k) + uu*(app-apm)
              cy3(j,k)=cy3(j,k) + uu*(app-apm)
              dy(j,k)=dy(j,k) + uu*apm
              dy2(j,k)=dy2(j,k) + uu*apm
              dy3(j,k)=dy3(j,k) + uu*apm
            enddo
c      Add part of source terms to diagonal LHS in this sweep (helps increase
c      diagonal dominance)
            do j=1,jdim-1
              if (ivmx .eq. 30 .or. ivmx .eq. 31) then
                cmuc=  blend(j,k,i)*cmuc1 + (1.-blend(j,k,i))*cmuc2
                betax= blend(j,k,i)*beta1 + (1.-blend(j,k,i))*beta2
                cyadd= 2.*re*betax*turre(j,k,i,1) +
     +                 ccabs(damp1(j,k,i))/turre(j,k,i,1)
                cy2add=re*cmuc*turre(j,k,i,1)
c               cy3add=???
              end if
              cy(j,k)=cy(j,k) + cyadd
              cy2(j,k)=cy2(j,k) + cy2add
c             cy3(j,k)=cy3(j,k) + cy3add
            enddo
            do j=1,jdim-1
              by(j,k)=by(j,k)*timestp(j,k,i)
              by2(j,k)=by2(j,k)*timestp(j,k,i)*factor2
              by3(j,k)=by3(j,k)*timestp(j,k,i)*factor3
              cy(j,k)=cy(j,k)*timestp(j,k,i)+1.0*(1.+phi)
              cy2(j,k)=cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.+phi)
              cy3(j,k)=cy3(j,k)*timestp(j,k,i)*factor3+1.0*(1.+phi)
              dy(j,k)=dy(j,k)*timestp(j,k,i)
              dy2(j,k)=dy2(j,k)*timestp(j,k,i)*factor2
              dy3(j,k)=dy3(j,k)*timestp(j,k,i)*factor3
              fy(j,k)=rhside(j,k,i,1)*timestp(j,k,i)
              fy2(j,k)=rhside(j,k,i,2)*timestp(j,k,i)*factor2
              fy3(j,k)=rhside(j,k,i,3)*timestp(j,k,i)*factor3
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-
     +                 turre(j,k,i,1))+phi*zksav2(j,k,i,4)
                fy2(j,k)=fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-
     +                 turre(j,k,i,2))+phi*zksav2(j,k,i,5)
                fy3(j,k)=fy3(j,k)+(1.+phi)*(zksav2(j,k,i,3)-
     +                 turre(j,k,i,3))+phi*zksav2(j,k,i,6)
              enddo
            end if
c
c         KDIM boundary points
            k=kdim-1
            kl=kdim-2
            ku=kdim-1
            do j=1,jdim-1
              dfacep=0.5*(blend(j,k,i)+blend(j,ku,i))
              dfacem=0.5*(blend(j,k,i)+blend(j,kl,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volku=volk0(j,i,3)
              xp=sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              yp=sk(j,k+1,i,2)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              zp=sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)
     +          +volku))
              volkl=vol(j,kl,i)
              xm=sk(j,k,i,1)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              ym=sk(j,k,i,2)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              zm=sk(j,k,i,3)*sk(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +volkl))
              xa=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,1)*sk(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,2)*sk(j,k,i,4))/vol(j,k,i)
              za=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4) + 
     +           sk(j,k,i,3)*sk(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volku)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volkl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              by(j,k)=-cdm
              cy(j,k)= cdp+cdm
              dy(j,k)=-cdp
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              by2(j,k)=-cdm
              cy2(j,k)= cdp+cdm
              dy2(j,k)=-cdp
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              by3(j,k)=-cdm
              cy3(j,k)= cdp+cdm
              dy3(j,k)=-cdp
            enddo
            do j=1,jdim-1
              xc=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,1)*sk(j,k,i  ,4))/vol(j,k,i)
              yc=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,2)*sk(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,3)*sk(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+
     +                sk(j,k,i  ,5)*sk(j,k,i  ,4))/vol(j,k,i)
              uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
              by(j,k)=by(j,k) - uu*app
              by2(j,k)=by2(j,k) - uu*app
              by3(j,k)=by3(j,k) - uu*app
              cy(j,k)=cy(j,k) + uu*(app-apm)
              cy2(j,k)=cy2(j,k) + uu*(app-apm)
              cy3(j,k)=cy3(j,k) + uu*(app-apm)
              dy(j,k)=dy(j,k) + uu*apm
              dy2(j,k)=dy2(j,k) + uu*apm
              dy3(j,k)=dy3(j,k) + uu*apm
            enddo
c      Add part of source terms to diagonal LHS in this sweep (helps increase
c      diagonal dominance)
            do j=1,jdim-1
              if (ivmx .eq. 30 .or. ivmx .eq. 31) then
                cmuc=  blend(j,k,i)*cmuc1 + (1.-blend(j,k,i))*cmuc2
                betax= blend(j,k,i)*beta1 + (1.-blend(j,k,i))*beta2
                cyadd= 2.*re*betax*turre(j,k,i,1) +
     +                 ccabs(damp1(j,k,i))/turre(j,k,i,1)
                cy2add=re*cmuc*turre(j,k,i,1)
c               cy3add=???
              end if
              cy(j,k)=cy(j,k) + cyadd
              cy2(j,k)=cy2(j,k) + cy2add
c             cy3(j,k)=cy3(j,k) + cy3add
            enddo
            do j=1,jdim-1
              by(j,k)=by(j,k)*timestp(j,k,i)
              by2(j,k)=by2(j,k)*timestp(j,k,i)*factor2
              by3(j,k)=by3(j,k)*timestp(j,k,i)*factor3
              cy(j,k)=cy(j,k)*timestp(j,k,i)+1.0*(1.+phi)
              cy2(j,k)=cy2(j,k)*timestp(j,k,i)*factor2+1.0*(1.+phi)
              cy3(j,k)=cy3(j,k)*timestp(j,k,i)*factor3+1.0*(1.+phi)
              dy(j,k)=dy(j,k)*timestp(j,k,i)
              dy2(j,k)=dy2(j,k)*timestp(j,k,i)*factor2
              dy3(j,k)=dy3(j,k)*timestp(j,k,i)*factor3
              fy(j,k)=rhside(j,k,i,1)*timestp(j,k,i)
              fy2(j,k)=rhside(j,k,i,2)*timestp(j,k,i)*factor2
              fy3(j,k)=rhside(j,k,i,3)*timestp(j,k,i)*factor3
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-
     +                 turre(j,k,i,1))+phi*zksav2(j,k,i,4)
                fy2(j,k)=fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-
     +                 turre(j,k,i,2))+phi*zksav2(j,k,i,5)
                fy3(j,k)=fy3(j,k)+(1.+phi)*(zksav2(j,k,i,3)-
     +                 turre(j,k,i,3))+phi*zksav2(j,k,i,6)
              enddo
            end if
          if (iover .eq. 1) then
            do k=1,kdim-1
              do j=1,jdim-1
                fy(j,k)=fy(j,k)*blank(j,k,i)
                by(j,k)=by(j,k)*blank(j,k,i)
                dy(j,k)=dy(j,k)*blank(j,k,i)
                cy(j,k)=cy(j,k)*blank(j,k,i)+(1.-blank(j,k,i))
                fy2(j,k)=fy2(j,k)*blank(j,k,i)
                by2(j,k)=by2(j,k)*blank(j,k,i)
                dy2(j,k)=dy2(j,k)*blank(j,k,i)
                cy2(j,k)=cy2(j,k)*blank(j,k,i)+(1.-blank(j,k,i))
                fy3(j,k)=fy3(j,k)*blank(j,k,i)
                by3(j,k)=by3(j,k)*blank(j,k,i)
                dy3(j,k)=dy3(j,k)*blank(j,k,i)
                cy3(j,k)=cy3(j,k)*blank(j,k,i)+(1.-blank(j,k,i))
              enddo
            enddo
          end if
          call triv(jdim-1,kdim-1,1,jdim-1,1,kdim-1,worky,by,cy,dy,fy)
          call triv(jdim-1,kdim-1,1,jdim-1,1,kdim-1,worky,by2,cy2,
     +     dy2,fy2)
          call triv(jdim-1,kdim-1,1,jdim-1,1,kdim-1,worky,by3,cy3,
     +     dy3,fy3)
          do k=1,kdim-1
            do j=1,jdim-1
              rhside(j,k,i,1)=fy(j,k)
              rhside(j,k,i,2)=fy2(j,k)
              rhside(j,k,i,3)=fy3(j,k)
            enddo
          enddo
        enddo
c
c    Implicit F_xi_xi viscous terms.  Do over all i's
        do i=1,idim-1
c         Interior points
          do j=2,jdim-2
            jl=j-1
            ju=j+1
            do k=1,kdim-1
              dfacep=0.5*(blend(j,k,i)+blend(ju,k,i))
              dfacem=0.5*(blend(j,k,i)+blend(jl,k,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volju=vol(ju,k,i)
              xp=sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              yp=sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              zp=sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              voljl=vol(jl,k,i)
              xm=sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              ym=sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              zm=sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              xa=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i)
              za=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              bx(k,j)=-cdm
              cx(k,j)= cdp+cdm
              dx(k,j)=-cdp
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              bx2(k,j)=-cdm
              cx2(k,j)= cdp+cdm
              dx2(k,j)=-cdp
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              bx3(k,j)=-cdm
              cx3(k,j)= cdp+cdm
              dx3(k,j)=-cdp
            enddo
            do k=1,kdim-1
              xc=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,1)*sj(j,k,i  ,4))/vol(j,k,i)
              yc=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,2)*sj(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,3)*sj(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,5)*sj(j,k,i  ,4))/vol(j,k,i)
              uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
              bx(k,j)=bx(k,j) - uu*app
              bx2(k,j)=bx2(k,j) - uu*app
              bx3(k,j)=bx3(k,j) - uu*app
              cx(k,j)=cx(k,j) + uu*(app-apm)
              cx2(k,j)=cx2(k,j) + uu*(app-apm)
              cx3(k,j)=cx3(k,j) + uu*(app-apm)
              dx(k,j)=dx(k,j) + uu*apm
              dx2(k,j)=dx2(k,j) + uu*apm
              dx3(k,j)=dx3(k,j) + uu*apm
            enddo
            do k=1,kdim-1
              bx(k,j)=bx(k,j)*timestp(j,k,i)
              bx2(k,j)=bx2(k,j)*timestp(j,k,i)*factor2
              bx3(k,j)=bx3(k,j)*timestp(j,k,i)*factor3
              cx(k,j)=cx(k,j)*timestp(j,k,i)+1.0*(1.+phi)
              cx2(k,j)=cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.+phi)
              cx3(k,j)=cx3(k,j)*timestp(j,k,i)*factor3+1.0*(1.+phi)
              dx(k,j)=dx(k,j)*timestp(j,k,i)
              dx2(k,j)=dx2(k,j)*timestp(j,k,i)*factor2
              dx3(k,j)=dx3(k,j)*timestp(j,k,i)*factor3
              fx(k,j)=rhside(j,k,i,1)*(1.+phi)
              fx2(k,j)=rhside(j,k,i,2)*(1.+phi)
              fx3(k,j)=rhside(j,k,i,3)*(1.+phi)
            enddo
          enddo
c
c         J0 boundary points
            j=1
            jl=1
            ju=min(2,jdim-1)
            do k=1,kdim-1
              dfacep=0.5*(blend(j,k,i)+blend(ju,k,i))
              dfacem=0.5*(blend(j,k,i)+blend(jl,k,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volju=vol(ju,k,i)
              xp=sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              yp=sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              zp=sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              voljl=volj0(k,i,1)
              xm=sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              ym=sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              zm=sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              xa=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i)
              za=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              bx(k,j)=-cdm
              cx(k,j)= cdp+cdm
              dx(k,j)=-cdp
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              bx2(k,j)=-cdm
              cx2(k,j)= cdp+cdm
              dx2(k,j)=-cdp
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              bx3(k,j)=-cdm
              cx3(k,j)= cdp+cdm
              dx3(k,j)=-cdp
            enddo
            do k=1,kdim-1
              xc=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,1)*sj(j,k,i  ,4))/vol(j,k,i)
              yc=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,2)*sj(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,3)*sj(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,5)*sj(j,k,i  ,4))/vol(j,k,i)
              uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
              bx(k,j)=bx(k,j) - uu*app
              bx2(k,j)=bx2(k,j) - uu*app
              bx3(k,j)=bx3(k,j) - uu*app
              cx(k,j)=cx(k,j) + uu*(app-apm)
              cx2(k,j)=cx2(k,j) + uu*(app-apm)
              cx3(k,j)=cx3(k,j) + uu*(app-apm)
              dx(k,j)=dx(k,j) + uu*apm
              dx2(k,j)=dx2(k,j) + uu*apm
              dx3(k,j)=dx3(k,j) + uu*apm
            enddo
            do k=1,kdim-1
              bx(k,j)=bx(k,j)*timestp(j,k,i)
              bx2(k,j)=bx2(k,j)*timestp(j,k,i)*factor2
              bx3(k,j)=bx3(k,j)*timestp(j,k,i)*factor3
              cx(k,j)=cx(k,j)*timestp(j,k,i)+1.0*(1.+phi)
              cx2(k,j)=cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.+phi)
              cx3(k,j)=cx3(k,j)*timestp(j,k,i)*factor3+1.0*(1.+phi)
              dx(k,j)=dx(k,j)*timestp(j,k,i)
              dx2(k,j)=dx2(k,j)*timestp(j,k,i)*factor2
              dx3(k,j)=dx3(k,j)*timestp(j,k,i)*factor3
              fx(k,j)=rhside(j,k,i,1)*(1.+phi)
              fx2(k,j)=rhside(j,k,i,2)*(1.+phi)
              fx3(k,j)=rhside(j,k,i,3)*(1.+phi)
            enddo
c
c         JDIM boundary points
            j=jdim-1
            jl=jdim-2
            ju=jdim-1
            do k=1,kdim-1
              dfacep=0.5*(blend(j,k,i)+blend(ju,k,i))
              dfacem=0.5*(blend(j,k,i)+blend(jl,k,i))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
              volju=volj0(k,i,3)
              xp=sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              yp=sj(j+1,k,i,2)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              zp=sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)
     +          +volju))
              voljl=vol(jl,k,i)
              xm=sj(j,k,i,1)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              ym=sj(j,k,i,2)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              zm=sj(j,k,i,3)*sj(j,k,i,4)/(0.5*(vol(j,k,i)
     +          +voljl))
              xa=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,1)*sj(j,k,i,4))/vol(j,k,i)
              ya=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,2)*sj(j,k,i,4))/vol(j,k,i)
              za=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4) + 
     +           sj(j,k,i,3)*sj(j,k,i,4))/vol(j,k,i)
c
              ttpo=xp*xa+yp*ya+zp*za
              ttmo=xm*xa+ym*ya+zm*za
              ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+volju)/
     +             vol(j,k,i)
              ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+voljl)/
     +             vol(j,k,i)
c             choose between weak (o) and strong (n) conservation form
              ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
              ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
              bx(k,j)=-cdm
              cx(k,j)= cdp+cdm
              dx(k,j)=-cdp
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
              bx2(k,j)=-cdm
              cx2(k,j)= cdp+cdm
              dx2(k,j)=-cdp
              cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
              bx3(k,j)=-cdm
              cx3(k,j)= cdp+cdm
              dx3(k,j)=-cdp
            enddo
            do k=1,kdim-1
              xc=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,1)*sj(j,k,i  ,4))/vol(j,k,i)
              yc=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,2)*sj(j,k,i  ,4))/vol(j,k,i)
              zc=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,3)*sj(j,k,i  ,4))/vol(j,k,i)
              tc=0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+
     +                sj(j,k,i  ,5)*sj(j,k,i  ,4))/vol(j,k,i)
              uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
              bx(k,j)=bx(k,j) - uu*app
              bx2(k,j)=bx2(k,j) - uu*app
              bx3(k,j)=bx3(k,j) - uu*app
              cx(k,j)=cx(k,j) + uu*(app-apm)
              cx2(k,j)=cx2(k,j) + uu*(app-apm)
              cx3(k,j)=cx3(k,j) + uu*(app-apm)
              dx(k,j)=dx(k,j) + uu*apm
              dx2(k,j)=dx2(k,j) + uu*apm
              dx3(k,j)=dx3(k,j) + uu*apm
            enddo
            do k=1,kdim-1
              bx(k,j)=bx(k,j)*timestp(j,k,i)
              bx2(k,j)=bx2(k,j)*timestp(j,k,i)*factor2
              bx3(k,j)=bx3(k,j)*timestp(j,k,i)*factor3
              cx(k,j)=cx(k,j)*timestp(j,k,i)+1.0*(1.+phi)
              cx2(k,j)=cx2(k,j)*timestp(j,k,i)*factor2+1.0*(1.+phi)
              cx3(k,j)=cx3(k,j)*timestp(j,k,i)*factor3+1.0*(1.+phi)
              dx(k,j)=dx(k,j)*timestp(j,k,i)
              dx2(k,j)=dx2(k,j)*timestp(j,k,i)*factor2
              dx3(k,j)=dx3(k,j)*timestp(j,k,i)*factor3
              fx(k,j)=rhside(j,k,i,1)*(1.+phi)
              fx2(k,j)=rhside(j,k,i,2)*(1.+phi)
              fx3(k,j)=rhside(j,k,i,3)*(1.+phi)
            enddo
          if (iover .eq. 1) then
            do k=1,kdim-1
              do j=1,jdim-1
                fx(k,j)=fx(k,j)*blank(j,k,i)
                bx(k,j)=bx(k,j)*blank(j,k,i)
                dx(k,j)=dx(k,j)*blank(j,k,i)
                cx(k,j)=cx(k,j)*blank(j,k,i)+(1.-blank(j,k,i))
                fx2(k,j)=fx2(k,j)*blank(j,k,i)
                bx2(k,j)=bx2(k,j)*blank(j,k,i)
                dx2(k,j)=dx2(k,j)*blank(j,k,i)
                cx2(k,j)=cx2(k,j)*blank(j,k,i)+(1.-blank(j,k,i))
                fx3(k,j)=fx3(k,j)*blank(j,k,i)
                bx3(k,j)=bx3(k,j)*blank(j,k,i)
                dx3(k,j)=dx3(k,j)*blank(j,k,i)
                cx3(k,j)=cx3(k,j)*blank(j,k,i)+(1.-blank(j,k,i))
              enddo
            enddo
          end if
          call triv(kdim-1,jdim-1,1,kdim-1,1,jdim-1,workx,bx,cx,dx,fx)
          call triv(kdim-1,jdim-1,1,kdim-1,1,jdim-1,workx,bx2,
     +     cx2,dx2,fx2)
          call triv(kdim-1,jdim-1,1,kdim-1,1,jdim-1,workx,bx3,
     +     cx3,dx3,fx3)
          do j=1,jdim-1
            do k=1,kdim-1
              rhside(j,k,i,1)=fx(k,j)
              rhside(j,k,i,2)=fx2(k,j)
              rhside(j,k,i,3)=fx3(k,j)
            enddo
          enddo
        enddo
c
c    Implicit F_zeta_zeta viscous terms.  Do over all j's
        if(i2d .ne. 1 .and. iaxi2planeturb .ne. 1) then
          do j=1,jdim-1
c           Interior points
            do i=2,idim-2
              il=i-1
              iu=i+1
              do k=1,kdim-1
              dfacep=0.5*(blend(j,k,i)+blend(j,k,iu))
              dfacem=0.5*(blend(j,k,i)+blend(j,k,il))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
                voliu=vol(j,k,iu)
                xp=si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                yp=si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                zp=si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                volil=vol(j,k,il)
                xm=si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                ym=si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                zm=si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                xa=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4) + 
     +             si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i)
                ya=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4) + 
     +             si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i)
                za=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4) + 
     +             si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i)
c
                ttpo=xp*xa+yp*ya+zp*za
                ttmo=xm*xa+ym*ya+zm*za
                ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/
     +               vol(j,k,i)
                ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/
     +               vol(j,k,i)
c               choose between weak (o) and strong (n) conservation form
                ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
                ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
                anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1))
                anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
                bz(k,i)=-cdm
                cz(k,i)= cdp+cdm
                dz(k,i)=-cdp
                cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
                bz2(k,i)=-cdm
                cz2(k,i)= cdp+cdm
                dz2(k,i)=-cdp
                cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
                bz3(k,i)=-cdm
                cz3(k,i)= cdp+cdm
                dz3(k,i)=-cdp
              enddo
              do k=1,kdim-1
                xc=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,1)*si(j,k,i  ,4))/vol(j,k,i)
                yc=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,2)*si(j,k,i  ,4))/vol(j,k,i)
                zc=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,3)*si(j,k,i  ,4))/vol(j,k,i)
                tc=0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,5)*si(j,k,i  ,4))/vol(j,k,i)
                uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
                sgnu=ccsignrc(1.,uu)
                app=0.5*(1.+sgnu)
                apm=0.5*(1.-sgnu)
                bz(k,i)=bz(k,i) - uu*app
                bz2(k,i)=bz2(k,i) - uu*app
                bz3(k,i)=bz3(k,i) - uu*app
                cz(k,i)=cz(k,i) + uu*(app-apm)
                cz2(k,i)=cz2(k,i) + uu*(app-apm)
                cz3(k,i)=cz3(k,i) + uu*(app-apm)
                dz(k,i)=dz(k,i) + uu*apm
                dz2(k,i)=dz2(k,i) + uu*apm
                dz3(k,i)=dz3(k,i) + uu*apm
              enddo
              do k=1,kdim-1
                bz(k,i)=bz(k,i)*timestp(j,k,i)
                bz2(k,i)=bz2(k,i)*timestp(j,k,i)*factor2
                bz3(k,i)=bz3(k,i)*timestp(j,k,i)*factor3
                cz(k,i)=cz(k,i)*timestp(j,k,i)+1.0*(1.+phi)
                cz2(k,i)=cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.+phi)
                cz3(k,i)=cz3(k,i)*timestp(j,k,i)*factor3+1.0*(1.+phi)
                dz(k,i)=dz(k,i)*timestp(j,k,i)
                dz2(k,i)=dz2(k,i)*timestp(j,k,i)*factor2
                dz3(k,i)=dz3(k,i)*timestp(j,k,i)*factor3
                fz(k,i)=rhside(j,k,i,1)*(1.+phi)
                fz2(k,i)=rhside(j,k,i,2)*(1.+phi)
                fz3(k,i)=rhside(j,k,i,3)*(1.+phi)
              enddo
            enddo
c
c           I0 boundary points
              i=1
              il=1
              iu=min(2,idim-1)
              do k=1,kdim-1
              dfacep=0.5*(blend(j,k,i)+blend(j,k,iu))
              dfacem=0.5*(blend(j,k,i)+blend(j,k,il))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
                voliu=vol(j,k,iu)
                xp=si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                yp=si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                zp=si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                volil=voli0(j,k,1)
                xm=si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                ym=si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                zm=si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                xa=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4) + 
     +             si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i)
                ya=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4) + 
     +             si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i)
                za=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4) + 
     +             si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i)
c
                ttpo=xp*xa+yp*ya+zp*za
                ttmo=xm*xa+ym*ya+zm*za
                ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/
     +               vol(j,k,i)
                ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/
     +               vol(j,k,i)
c               choose between weak (o) and strong (n) conservation form
                ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
                ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
                anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1))
                anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
                bz(k,i)=-cdm
                cz(k,i)= cdp+cdm
                dz(k,i)=-cdp
                cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
                bz2(k,i)=-cdm
                cz2(k,i)= cdp+cdm
                dz2(k,i)=-cdp
                cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
                bz3(k,i)=-cdm
                cz3(k,i)= cdp+cdm
                dz3(k,i)=-cdp
              enddo
              do k=1,kdim-1
                xc=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,1)*si(j,k,i  ,4))/vol(j,k,i)
                yc=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,2)*si(j,k,i  ,4))/vol(j,k,i)
                zc=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,3)*si(j,k,i  ,4))/vol(j,k,i)
                tc=0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,5)*si(j,k,i  ,4))/vol(j,k,i)
                uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
                sgnu=ccsignrc(1.,uu)
                app=0.5*(1.+sgnu)
                apm=0.5*(1.-sgnu)
                bz(k,i)=bz(k,i) - uu*app
                bz2(k,i)=bz2(k,i) - uu*app
                bz3(k,i)=bz3(k,i) - uu*app
                cz(k,i)=cz(k,i) + uu*(app-apm)
                cz2(k,i)=cz2(k,i) + uu*(app-apm)
                cz3(k,i)=cz3(k,i) + uu*(app-apm)
                dz(k,i)=dz(k,i) + uu*apm
                dz2(k,i)=dz2(k,i) + uu*apm
                dz3(k,i)=dz3(k,i) + uu*apm
              enddo
              do k=1,kdim-1
                bz(k,i)=bz(k,i)*timestp(j,k,i)
                bz2(k,i)=bz2(k,i)*timestp(j,k,i)*factor2
                bz3(k,i)=bz3(k,i)*timestp(j,k,i)*factor3
                cz(k,i)=cz(k,i)*timestp(j,k,i)+1.0*(1.+phi)
                cz2(k,i)=cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.+phi)
                cz3(k,i)=cz3(k,i)*timestp(j,k,i)*factor3+1.0*(1.+phi)
                dz(k,i)=dz(k,i)*timestp(j,k,i)
                dz2(k,i)=dz2(k,i)*timestp(j,k,i)*factor2
                dz3(k,i)=dz3(k,i)*timestp(j,k,i)*factor3
                fz(k,i)=rhside(j,k,i,1)*(1.+phi)
                fz2(k,i)=rhside(j,k,i,2)*(1.+phi)
                fz3(k,i)=rhside(j,k,i,3)*(1.+phi)
              enddo
c
c           IDIM boundary points
              i=idim-1
              il=idim-2
              iu=idim-1
              do k=1,kdim-1
              dfacep=0.5*(blend(j,k,i)+blend(j,k,iu))
              dfacem=0.5*(blend(j,k,i)+blend(j,k,il))
              sigkp=dfacep*sigk1+(1.-dfacep)*sigk2
              sigkm=dfacem*sigk1+(1.-dfacem)*sigk2
              sigop=dfacep*sigo1+(1.-dfacep)*sigo2
              sigom=dfacem*sigo1+(1.-dfacem)*sigo2
                voliu=voli0(j,k,3)
                xp=si(j,k,i+1,1)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                yp=si(j,k,i+1,2)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                zp=si(j,k,i+1,3)*si(j,k,i+1,4)/(0.5*(vol(j,k,i)
     +            +voliu))
                volil=vol(j,k,il)
                xm=si(j,k,i,1)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                ym=si(j,k,i,2)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                zm=si(j,k,i,3)*si(j,k,i,4)/(0.5*(vol(j,k,i)
     +            +volil))
                xa=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4) + 
     +             si(j,k,i,1)*si(j,k,i,4))/vol(j,k,i)
                ya=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4) + 
     +             si(j,k,i,2)*si(j,k,i,4))/vol(j,k,i)
                za=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4) + 
     +             si(j,k,i,3)*si(j,k,i,4))/vol(j,k,i)
c
                ttpo=xp*xa+yp*ya+zp*za
                ttmo=xm*xa+ym*ya+zm*za
                ttpn=(xp*xp+yp*yp+zp*zp)*0.5*(vol(j,k,i)+voliu)/
     +               vol(j,k,i)
                ttmn=(xm*xm+ym*ym+zm*zm)*0.5*(vol(j,k,i)+volil)/
     +               vol(j,k,i)
c               choose between weak (o) and strong (n) conservation form
                ttp=ttpo*(1-istrongturbdis)+ttpn*istrongturbdis
                ttm=ttmo*(1-istrongturbdis)+ttmn*istrongturbdis
c
                anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1))
                anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+sigop*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(fnum+sigom*anutm)*ttm/(q(j,k,i,1)*re)
                bz(k,i)=-cdm
                cz(k,i)= cdp+cdm
                dz(k,i)=-cdp
                cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(q(j,k,i,1)*re)
                bz2(k,i)=-cdm
                cz2(k,i)= cdp+cdm
                dz2(k,i)=-cdp
                cdp=(sigkmu*fnup+sigg*anutp)*ttp/(q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigg*anutm)*ttm/(q(j,k,i,1)*re)
                bz3(k,i)=-cdm
                cz3(k,i)= cdp+cdm
                dz3(k,i)=-cdp
              enddo
              do k=1,kdim-1
                xc=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,1)*si(j,k,i  ,4))/vol(j,k,i)
                yc=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,2)*si(j,k,i  ,4))/vol(j,k,i)
                zc=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,3)*si(j,k,i  ,4))/vol(j,k,i)
                tc=0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+
     +                  si(j,k,i  ,5)*si(j,k,i  ,4))/vol(j,k,i)
                uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc
                sgnu=ccsignrc(1.,uu)
                app=0.5*(1.+sgnu)
                apm=0.5*(1.-sgnu)
                bz(k,i)=bz(k,i) - uu*app
                bz2(k,i)=bz2(k,i) - uu*app
                bz3(k,i)=bz3(k,i) - uu*app
                cz(k,i)=cz(k,i) + uu*(app-apm)
                cz2(k,i)=cz2(k,i) + uu*(app-apm)
                cz3(k,i)=cz3(k,i) + uu*(app-apm)
                dz(k,i)=dz(k,i) + uu*apm
                dz2(k,i)=dz2(k,i) + uu*apm
                dz3(k,i)=dz3(k,i) + uu*apm
              enddo
              do k=1,kdim-1
                bz(k,i)=bz(k,i)*timestp(j,k,i)
                bz2(k,i)=bz2(k,i)*timestp(j,k,i)*factor2
                bz3(k,i)=bz3(k,i)*timestp(j,k,i)*factor3
                cz(k,i)=cz(k,i)*timestp(j,k,i)+1.0*(1.+phi)
                cz2(k,i)=cz2(k,i)*timestp(j,k,i)*factor2+1.0*(1.+phi)
                cz3(k,i)=cz3(k,i)*timestp(j,k,i)*factor3+1.0*(1.+phi)
                dz(k,i)=dz(k,i)*timestp(j,k,i)
                dz2(k,i)=dz2(k,i)*timestp(j,k,i)*factor2
                dz3(k,i)=dz3(k,i)*timestp(j,k,i)*factor3
                fz(k,i)=rhside(j,k,i,1)*(1.+phi)
                fz2(k,i)=rhside(j,k,i,2)*(1.+phi)
                fz3(k,i)=rhside(j,k,i,3)*(1.+phi)
              enddo
            if (iover .eq. 1) then
              do i=1,idim-1
                do k=1,kdim-1
                  fz(k,i)=fz(k,i)*blank(j,k,i)
                  bz(k,i)=bz(k,i)*blank(j,k,i)
                  dz(k,i)=dz(k,i)*blank(j,k,i)
                  cz(k,i)=cz(k,i)*blank(j,k,i)+(1.-blank(j,k,i))
                  fz2(k,i)=fz2(k,i)*blank(j,k,i)
                  bz2(k,i)=bz2(k,i)*blank(j,k,i)
                  dz2(k,i)=dz2(k,i)*blank(j,k,i)
                  cz2(k,i)=cz2(k,i)*blank(j,k,i)+(1.-blank(j,k,i))
                  fz3(k,i)=fz3(k,i)*blank(j,k,i)
                  bz3(k,i)=bz3(k,i)*blank(j,k,i)
                  dz3(k,i)=dz3(k,i)*blank(j,k,i)
                  cz3(k,i)=cz3(k,i)*blank(j,k,i)+(1.-blank(j,k,i))
                enddo
              enddo
            end if
            call triv(kdim-1,idim-1,1,kdim-1,1,idim-1,workz,bz,
     +                cz,dz,fz)
            call triv(kdim-1,idim-1,1,kdim-1,1,idim-1,workz,bz2,
     +                cz2,dz2,fz2)
            call triv(kdim-1,idim-1,1,kdim-1,1,idim-1,workz,bz3,
     +                cz3,dz3,fz3)
            do i=1,idim-1
              do k=1,kdim-1
                rhside(j,k,i,1)=fz(k,i)
                rhside(j,k,i,2)=fz2(k,i)
                rhside(j,k,i,3)=fz3(k,i)
              enddo
            enddo
          enddo
        end if
c
c    Update TURRE
        sumno=0.
        sumnk=0.
        sumni=0.
        negno=0
        negnk=0
        negni=0
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              sumno=sumno+rhside(j,k,i,1)**2
              if((real(turre(j,k,i,1)+rhside(j,k,i,1))) .le. 
     +          real(tur1cutlev)) then
                negno=negno+1
                if (real(tur1cut).gt.0.) turre(j,k,i,1)=tur1cut
              else
                turre(j,k,i,1)=turre(j,k,i,1)+rhside(j,k,i,1)
              end if
c
              sumnk=sumnk+rhside(j,k,i,2)**2
              if((real(turre(j,k,i,2)+rhside(j,k,i,2))) .le. 
     +          real(tur2cutlev)) then
                negnk=negnk+1
                if (real(tur2cut).gt.0.) turre(j,k,i,2)=tur2cut
              else
                turre(j,k,i,2)=turre(j,k,i,2)+rhside(j,k,i,2)
              end if
c
              sumni=sumni+rhside(j,k,i,3)**2
c             following section not done if transition off (SST only)
              if (itrans_on .eq. 1) then
              if((real(turre(j,k,i,3)+rhside(j,k,i,3))) .le.
     +          0.0) then
                negni=negni+1
                turre(j,k,i,3)=0.0
              else if ((real(turre(j,k,i,3)+rhside(j,k,i,3))) .ge.
     +          100.0) then
                turre(j,k,i,3)=100.0
              else
                turre(j,k,i,3)=turre(j,k,i,3)+rhside(j,k,i,3)
              end if
              end if
            enddo
          enddo
        enddo
        sumno=sqrt(sumno)/float((kdim-1)*(jdim-1)*(idim-1))
        sumnk=sqrt(sumnk)/float((kdim-1)*(jdim-1)*(idim-1))
        sumni=sqrt(sumni)/float((kdim-1)*(jdim-1)*(idim-1))
c
c       if(iwrite .eq. 1) then
c         write(15,'('' icyc, it, logreso, negno, logresk, negnk'',
c    +     '', block'',/,2i5,e15.5,i5,e15.5,2i5)') icyc,not,
c    +     real(cclog10(sumno)),negno,real(cclog10(sumnk)),negnk,nbl
c       end if
c
 500  continue
      sumn1 = sumno
      sumn2 = sumnk
      sumn3 = sumni
      negn1 = negno
      negn2 = negnk
      negn3 = negni
c
c Update VIST3D and save omega and k values
      if (ivmx .eq. 30 .or. ivmx .eq. 31) then
        do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            denom1=2./re*sqrt(turre(j,k,i,2))/(.09*turre(j,k,i,1)*
     +       ccabs(smin(j,k,i)))
            denom2=500.*fnu(j,k,i)/(q(j,k,i,1)*smin(j,k,i)*re*re*
     +       smin(j,k,i)*turre(j,k,i,1))
            arg2=ccmax(denom1,denom2)
            f2=cctanh(arg2*arg2)
            if (isstdenom .eq. 1) then
c            Determine Sij values:
              s11 = ux(j,k,i,1)
              s22 = ux(j,k,i,5)
              s33 = ux(j,k,i,9)
              s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
              s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
              s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
            xxx=sqrt(2.*xis)
            denom=ccmax(a1*turre(j,k,i,1),xxx*f2/re)
            else
            denom=ccmax(a1*turre(j,k,i,1),vor(j,k,i)*f2/re)
            end if
            vist3d(j,k,i)=a1*q(j,k,i,1)*turre(j,k,i,2)/denom
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),edvislim)
            zksav(j,k,i,1)=turre(j,k,i,1)
            zksav(j,k,i,2)=turre(j,k,i,2)
            zksav(j,k,i,3)=turre(j,k,i,3)
          enddo
        enddo
        enddo
      end if
c
      if (i_lam_forcezero .eq. 1) then
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              if ((i.ge.ilamlo .and. i.lt.ilamhi .and.
     .             j.ge.jlamlo .and. j.lt.jlamhi .and.
     .             k.ge.klamlo .and. k.lt.klamhi) .or.
     .             real(smin(j,k,i)) .lt. 0.) then
                vist3d(j,k,i)=0.
              end if
            enddo
          enddo
        enddo
      end if
c
c-----BEGIN HARDWIRED OUTPUT SECTION
c
c  If desired, write out turb information here during last iteration:
c  This currently writes out assuming body is at k=1, and assumes
c  utau lies along the same j & i indices;  also, only i=1 index is written
c
c  NOTE:
c
c     ONLY SET IWRITEAUX=1 IF USING THE SEQUENTIAL BUILD OF THE
c     CODE - OTHERWISE, IF USED IN PARALELL, EACH PROCESSOR WILL
c     OVERWRITE UNITS 91 AND 92
c
c     THIS SECTION IS PRETTY MUCH HARDWIRED FOR 2D, SINGLE BLOCK
c     CASES ANYWAY!
c
      iwriteaux=0
      if(iwriteaux .eq. 1) then
      if(icyc.eq.ncyc1(1) .or. icyc.eq.ncyc1(2) .or. icyc.eq.ncyc1(3)
     +.or. icyc.eq.ncyc1(4) .or. icyc.eq.ncyc1(5)) then
      if(ntime .eq. 1) then
c
      jset=int(.779*float(jdim))
      write(92,'(i5,'' jset,x='',i5,e12.5,''  u+,k+,e+,-uv+,logy+'')')
     + kdim-1,jset,real(0.5*(x(jset,1,1)+x(jset+1,1,1)))
      qset=sqrt((q(jset,1,1,2)-qk0(jset,1,2,1))**2+
     +          (q(jset,1,1,3)-qk0(jset,1,3,1))**2+
     +          (q(jset,1,1,4)-qk0(jset,1,4,1))**2)
      utau=sqrt((fnu(jset,1,1)+vk0(jset,1,1,1))*qset/
     +     (ccabs(smin(jset,1,1))*q(jset,1,1,1)*re))
      do k=1,kdim-1
        ypl=re*q(jset,1,1,1)*utau*ccabs(smin(jset,k,1))/fnu(jset,1,1)
        ypl=cclog10(ypl)
        uplus=sqrt(q(jset,k,1,2)**2+q(jset,k,1,3)**2+q(jset,k,1,4)**2)
     +        /utau
        zkplus=turre(jset,k,1,2)/(utau**2)
c     Note: the following works for k-epsilon model only:
        eplus=turre(jset,k,1,1)*re/(utau*utau*qset/
     +        ccabs(smin(jset,1,1)))
        uvplus=vist3d(jset,k,1)*vor(jset,k,1)/
     +         (re*utau*utau*q(jset,k,1,1))
        write(92,'(5e15.5)') real(uplus),real(zkplus),real(eplus),
     +   real(uvplus),real(ypl)
      enddo
c
      nnumb=2
      write(91,'(i5)') nnumb
      write(91,'(''   uv'')')
      write(91,'(''   y'')')
      jset=68
      write(91,'(i5,'' jset='',i5,''  -uv/uinf**2,y  x='',e15.5
     + )') kdim-1,jset,real(0.5*(x(jset,1,1)+x(jset+1,1,1)))
      do k=1,kdim-1
        write(91,'(2e15.5)') 
     +  real(vist3d(jset,k,1)*vor(jset,k,1)/q(jset,k,1,1)/(reue*xmach)),
     +  real(ccabs(smin(jset,k,1)))
      enddo
      jset=93
      write(91,'(i5,'' jset='',i5,''  -uv/uinf**2,y  x='',e15.5
     + )') kdim-1,jset,real(0.5*(x(jset,1,1)+x(jset+1,1,1)))
      do k=1,kdim-1
        write(91,'(2e15.5)') 
     +  real(vist3d(jset,k,1)*vor(jset,k,1)/q(jset,k,1,1)/(reue*xmach)),
     +  real(ccabs(smin(jset,k,1)))
      enddo
      jset=107
      write(91,'(i5,'' jset='',i5,''  -uv/uinf**2,y  x='',e15.5
     + )') kdim-1,jset,real(0.5*(x(jset,1,1)+x(jset+1,1,1)))
      do k=1,kdim-1
        write(91,'(2e15.5)') 
     +  real(vist3d(jset,k,1)*vor(jset,k,1)/q(jset,k,1,1)/(reue*xmach)),
     +  real(ccabs(smin(jset,k,1)))
      enddo
      jset=113
      write(91,'(i5,'' jset='',i5,''  -uv/uinf**2,y  x='',e15.5
     + )') kdim-1,jset,real(0.5*(x(jset,1,1)+x(jset+1,1,1)))
      do k=1,kdim-1
        write(91,'(2e15.5)') 
     +  real(vist3d(jset,k,1)*vor(jset,k,1)/q(jset,k,1,1)/(reue*xmach)),
     +  real(ccabs(smin(jset,k,1)))
      enddo
      jset=123
      write(91,'(i5,'' jset='',i5,''  -uv/uinf**2,y  x='',e15.5
     + )') kdim-1,jset,real(0.5*(x(jset,1,1)+x(jset+1,1,1)))
      do k=1,kdim-1
        write(91,'(2e15.5)') 
     +  real(vist3d(jset,k,1)*vor(jset,k,1)/q(jset,k,1,1)/(reue*xmach)),
     +  real(ccabs(smin(jset,k,1)))
      enddo
      jset=129
      write(91,'(i5,'' jset='',i5,''  -uv/uinf**2,y  x='',e15.5
     + )') kdim-1,jset,real(0.5*(x(jset,1,1)+x(jset+1,1,1)))
      do k=1,kdim-1
        write(91,'(2e15.5)') 
     +  real(vist3d(jset,k,1)*vor(jset,k,1)/q(jset,k,1,1)/(reue*xmach)),
     +  real(ccabs(smin(jset,k,1)))
      enddo
      jset=134
      write(91,'(i5,'' jset='',i5,''  -uv/uinf**2,y  x='',e15.5
     + )') kdim-1,jset,real(0.5*(x(jset,1,1)+x(jset+1,1,1)))
      do k=1,kdim-1
        write(91,'(2e15.5)') 
     +  real(vist3d(jset,k,1)*vor(jset,k,1)/q(jset,k,1,1)/(reue*xmach)),
     +  real(ccabs(smin(jset,k,1)))
      enddo
      end if
      end if
      end if
c
c-----END HARDWIRED OUTPUT SECTION
c
      return
      end
