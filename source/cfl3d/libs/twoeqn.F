      subroutine twoeqn(jdim,kdim,idim,q,sj,sk,si,vol,dtj,
     + x,y,z,vist3d,vor,smin,zksav,turre,damp1,blend,timestp,fnu,bx,bx2,
     + cx,cx2,dx,dx2,fx,fx2,workx,by,by2,cy,cy2,dy,dy2,fy,fy2,worky,
     + bz,bz2,cz,cz2,dz,dz2,fz,fz2,workz,ntime,tj0,tk0,ti0,nbl,
     + qj0,qk0,qi0,vj0,vk0,vi0,blank,iover,sumn1,sumn2,negn1,negn2,
     + ux,rhside,zksav2,v3dtmp,cmuv,bcj,bck,bci,nbci0,nbcidim,
     + nbcj0,nbcjdim,nbck0,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     + maxbl,maxseg,volj0,volk0,voli0,nou,bou,nbuf,ibufdim,iex,iex2,
     + iex3,
     + dkdx,dkdy,dkdz,dzdx,dzdy,dzdz,dw1dx,dw1dy,dw1dz,
     + dw2dx,dw2dy,dw2dz,dw3dx,dw3dy,dw3dz,drdx,drdy,drdz,
     + dpdx,dpdy,dpdz,
     + dkdj,dkdk,dkdi,dzdj,dzdk,dzdi,dw1dj,dw1dk,dw1di,
     + dw2dj,dw2dk,dw2di,dw3dj,dw3dk,dw3di,drdj,drdk,drdi,
     + dpdj,dpdk,dpdi,
     + w1,w2,w3,srce,vx,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute turbulent viscosity distributions using
c     2-equation turbulence models
c     (reprogrammed by Rumsey for CFL3D Version 5 - 5/96)
c     Alterations to EASM #8,9,13,14 by Rumsey - 5/00)
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
      common /easmv/ c1,c2,c3,c4,c5,sigk1,cmuc1
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /lam/ ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi
      common /maxiv/ ivmx
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /twod/ i2d
      common /zero/ iexp
      common /wallfun/ iwf(3)
      common /turbconv/ cflturb,edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /konew/ ikoprod,isstdenom,pklimterm,ibeta8kzeta,i_bsl,
     .        keepambient,re_thetat0,i_wilcox06,i_wilcox06_chiw,
     .        i_turbprod_kterm,i_catris_kw,prod2d3dtrace,
     .        i_compress_correct
      common /curvat/ isarc2d,sarccr3,ieasmcc2d,isstrc,sstrc_crc
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
      dimension turre(0-iex:jdim+iex,0-iex:kdim+iex,0-iex2:idim+iex2,2),
     + blank(jdim,kdim,idim),
     + zksav(jdim,kdim,idim,nummem),rhside(jdim-1,kdim-1,idim-1,2)
     + ,tj0(kdim,idim-1,nummem,4),tk0(jdim,idim-1,nummem,4),
     +  ti0(jdim,kdim,nummem,4),
     + qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4)
     +,vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),vi0(jdim,kdim,1,4)
      dimension ux(jdim-1,kdim-1,idim-1,9),
     + zksav2(jdim,kdim,idim,2*nummem),
     + v3dtmp(0:jdim,0:kdim,0-iex3:idim+iex3)
      dimension cmuv(jdim-1,kdim-1,idim-1)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     +          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     +          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension volj0(kdim,idim-1,4),
     +          volk0(jdim,idim-1,4),voli0(jdim,kdim,4)
      dimension  dkdx(jdim),dkdy(jdim),dkdz(jdim),
     +           dzdx(jdim),dzdy(jdim),dzdz(jdim),
     +           dw1dx(jdim),dw1dy(jdim),dw1dz(jdim),
     +           dw2dx(jdim),dw2dy(jdim),dw2dz(jdim),
     +           dw3dx(jdim),dw3dy(jdim),dw3dz(jdim),
     +           drdx(jdim),drdy(jdim),drdz(jdim),
     +           dpdx(jdim),dpdy(jdim),dpdz(jdim)
      dimension dkdj(jdim),dkdk(jdim),dkdi(jdim),
     +          dzdj(jdim),dzdk(jdim),dzdi(jdim),
     +          dw1dj(jdim),dw1dk(jdim),dw1di(jdim),
     +          dw2dj(jdim),dw2dk(jdim),dw2di(jdim),
     +          dw3dj(jdim),dw3dk(jdim),dw3di(jdim),
     +          drdj(jdim),drdk(jdim),drdi(jdim),
     +          dpdj(jdim),dpdk(jdim),dpdi(jdim)
      dimension w1(0:jdim,0:kdim,0:idim),
     +          w2(0:jdim,0:kdim,0:idim),
     +          w3(0:jdim,0:kdim,0:idim)
      dimension srce(0:jdim,0:kdim,0:idim,2)
      dimension vx(0:jdim,0:kdim,idim-1,4)
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
c      zksav - saved values of the 2 turbulent quantities (1=omega or epsilon,
c              2=k)
c                omega is nondimensionalized by rhoref*aref**2/mulamref
c                epsilon          "             rhoref*aref**4/mulamref
c                k                "             aref**2
c      zksav2 - auxiliary k and omega (or epsilon) for saving quantities
c               from last time step during time-accurate subiterations
c               (k,j,i,1) and (j,k,i,2) = k and omega (or epsilon)
c               (k,j,i,3) and (j,k,i,4) = delta k and delta omega (or epsilon)
c      turre - working values of 2 turb quantities in this subroutine
c      damp1 - for SST=CD (cross-derivative) term; for others, sometimes used
c              as storage for linearizations of production terms for addition
c              to LHS diagonal elements
c      blend - blending term for 2-layer models (just SST right now); 
c              =0 Wilcox06
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
c      sumn1,sumn2 - residual for 2 turbulence equations
c      negn1,negn2 - number of locations where the solution yields a "negative"
c                    turbulence quantity.  This SHOULD be zero, but it seems
c                    to be OK if there are only a few of these.  When they
c                    go negative, the values are artificially limited to be
c                    > 0.  When these numbers are large, it indicates that
c                    the solution is probably going to blow up.  Check your
c                    grid for excessive grid stretching, or try lowering CFL.
c      ux - 9 components of velocity derivative: ux,uy,uz,vx,vy,vz,wx,wy,wz 
c           at cell centers.  Used for nonlinear models only.
c      rhside - right-hand-side terms for 2 eqns
c      v3dtmp - temporary storage for vist3d (needed for nonlinear models,
c               which require turb viscosity in diffusion term WITHOUT
c               variable cmu)
c      cmuv   - storage for variable cmu in model #13
c      bcj,bck,bci - =0 for cell-center BC, =1 for at-face-center BC
c      nbci0,nbcidim,nbcj0,nbcjdim,nbck0,nbckdim - no. of BC segments each face
c      ibcinfo,jbcinfo,kbcinfo - gives i,j,k start and end indices for each
c               BC segment, among other things
c      maxbl,maxseg - dimensions used in nbci0,etc and ibcinfo, etc arrays
c
      if(isklton .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     Computing turbulent'',
     +'' viscosity using 2-eqns, block='',i5)') nbl
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     Freestream tur10,tur20 = '',
     +     2e19.8)') real(tur10(1)),real(tur10(2))
        if(iturbord .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     1st order advection on RHS'')')
        else
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     2nd order advection on RHS'')')
        end if
      end if
      if(isklton .gt. 0) then
        if(ivmx .eq. 6) then
           if(i_wilcox06 .eq. 1) then
             if(i_wilcox06_chiw .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     k-omega (Wilcox 06),'',
     +       '' with vortex stretching term'')')
             else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     k-omega (Wilcox 06),'',
     +       '' without vortex stretching term'')')
             end if
           else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     k-omega (Wilcox 88)'')')
           end if
           if(ikoprod .eq. 1) then
             if(i_turbprod_kterm .ne. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     +       '' term'')')
             else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     +       '' term with 2/3rho*k subtracted'')')
             end if
           else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     approx (vort) production'',
     .         '' term'')')
           end if
        end if
        if(ivmx .eq. 7 .and. i_bsl .ne. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-omega SST (Menter)'')')
           if(ikoprod .eq. 1) then
             if(i_turbprod_kterm .ne. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     +       '' term'')')
             else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     +       '' term with 2/3rho*k subtracted'')')
             end if
           else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     approx (vort) production'',
     .         '' term'')')
           end if
           if(isstdenom .eq. 0) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     vort in denom of mut term'')')
           else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     S in denom of mut term'')')
           end if
        else if(ivmx .eq. 7 .and. i_bsl .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-omega BSL (Menter)'')')
           if(ikoprod .eq. 1) then
             if(i_turbprod_kterm .ne. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     +       '' term'')')
             else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     +       '' term with 2/3rho*k subtracted'')')
             end if
           else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     approx (vort) production'',
     .         '' term'')')
           end if
        end if
        if(ivmx .eq. 8) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-omega (EASM, '',
     .   ''linear)'')')
           if(iturbprod .eq. 0) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     .         '' term'')')
           else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     exact production term'')')
           end if
           if (ieasmcc2d .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     curvature correction ON'')')
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''       WARNING: curv terms are'',
     +      '' active in 2-D sense only!!!'')')
           end if
        end if
        if(ivmx .eq. 9) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-epsilon (EASM,'',
     .   ''linear)'')')
           if(iturbprod .eq. 0) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     .         '' term'')')
           else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     exact production term'')')
           end if
           if (ieasmcc2d .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     curvature correction ON'')')
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''       WARNING: curv terms are'',
     +      '' active in 2-D sense only!!!'')')
           end if
        end if
        if(ivmx .eq.10) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-epsilon (Abid)'')')
           if(ikoprod .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     +       '' term'')')
           else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     approx (vort) production'',
     .         '' term'')')
           end if
        end if
        if(ivmx .eq.11) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-epsilon (G-S EASM, '',
     .   ''nonlinear)'')')
           if (ieasmcc2d .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     curvature correction ON'')')
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''       WARNING: curv terms are'',
     +      '' active in 2-D sense only!!!'')')
           end if
        end if
        if(ivmx .eq.12) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-omega (G-S EASM, '',
     .   ''nonlinear)'')')
           if (ieasmcc2d .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     curvature correction ON'')')
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''       WARNING: curv terms are'',
     +      '' active in 2-D sense only!!!'')')
           end if
        end if
        if(ivmx .eq.13) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-epsilon (EASM,'',
     .   ''nonlinear)'')')
           if(iturbprod .eq. 0) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     .         '' term'')')
           else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     exact production term'')')
           end if
           if (ieasmcc2d .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     curvature correction ON'')')
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''       WARNING: curv terms are'',
     +      '' active in 2-D sense only!!!'')')
           end if
        end if
        if(ivmx .eq.14) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-omega (EASM,'',
     .   ''nonlinear)'')')
           if(iturbprod .eq. 0) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     .         '' term'')')
           else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     exact production term'')')
           end if
           if (ieasmcc2d .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     curvature correction ON'')')
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''       WARNING: curv terms are'',
     +      '' active in 2-D sense only!!!'')')
           end if
        end if
        if(ivmx .eq.15) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     k-enstrophy (Note: default '',
     .   ''tur1cut=-1 if no keyword)'')')
           if(ikoprod .eq. 2) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     exact production term'')')
           else if (ikoprod .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     strain-based production'',
     .         '' term'')')
           else
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     approx (vort) production'',
     .         '' term'')')
           end if
        end if
c
        if(ikoprod .eq. 1 .and. (ivmx.eq.6.or.ivmx.eq.7.or.
     .    ivmx.eq.10)) then
        if(abs(real(prod2d3dtrace)-0.5) .lt. 0.001) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     Sij used in 2SijSij prod'',
     .         '' term forced to be traceless in 2-D sense'')')
        else if(abs(real(prod2d3dtrace)-0.33333333) .lt. 0.001) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     Sij used in 2SijSij prod'',
     .         '' term forced to be traceless in 3-D sense'')')
        end if
        end if
        if(ivmx.eq.6.and.i_wilcox06.eq.1) then
        if(abs(real(prod2d3dtrace)-0.5) .lt. 0.001) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     Sij used in 2SijSij in W06'',
     .        '' stress-limiter forced to be traceless in 2-D sense'')')
        else if(abs(real(prod2d3dtrace)-0.33333333) .lt. 0.001) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     Sij used in 2SijSij in W06'',
     .        '' stress-limiter forced to be traceless in 3-D sense'')')
        end if
        end if
c
        if(ivmx .eq. 9 .or. ivmx .eq.10 .or. ivmx .eq.11 .or.
     .     ivmx .eq. 13.or. ivmx .eq .15) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     WARNING: k-epsilon and '',
     .     ''k-enstrophy models sometimes fail to go turbulent!'')')
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''         We recommend either'',
     .     '' restarting from a different converged model,'')')
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''         or, if starting from'',
     .     '' scratch, freezing the model to its initialzed'')')
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''         levels (using keyword'',
     .     '' NFREEZE) until the flowfield is converged enough'')')
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''         for turbulence to'',
     .     '' sustain itself'')')
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     ALWAYS check vist3d levels to'',
     .     '' insure turbulence has tripped!'')')
           if (ikoprod .ne. 0) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     Also check that there is'',
     .       '' not excessive turbulence at stagnation regions.'')')
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''         Try ikoprod=0 if'',
     .       '' there is.'')')
           end if
        end if
        if(ivmx .eq. 9 .or. ivmx .eq.10 .or. ivmx .eq.11 .or.
     .     ivmx .eq. 13) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     Also, k-e models are NOT'',
     .     '' RECOMMENDED for'')')
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''         adverse-pressure-gradient'',
     .     '' wall-bounded flows.'')')
           if (iturbord .eq. 2) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     WARNING: 2nd order'',
     .       '' advection (iturbord=2)'')')
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''          can converge poorly'',
     .       '' for k-epsilon models'')')
           end if
        end if
        if(ivmx .eq. 11 .or. ivmx .eq. 12 .or. ivmx .eq. 13 .or.
     .     ivmx .eq. 14) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     WARNING: when using nonlinear'',
     .     '' models, the grid'')')
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''         needs to have min y+ of'',
     .     '' order 1.  If it is much larger,'')')
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''         you may need'',
     .     '' WallFunctions to keep from blowing up!'')')
        end if
        if(ivmx .eq. 6 .or. ivmx .eq. 7 .or. ivmx .eq. 8 .or.
     .     ivmx .eq. 12.or. ivmx .eq.14) then
           if (i_catris_kw .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     Catris compressible '',
     .       ''correction for turb diffusion terms being used'')')
           end if
        end if
        if(ivmx .eq. 6 .or. ivmx .eq. 7) then
           if (i_compress_correct .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     Wilcox-type dilatation-'',
     .       ''dissipation compressibility correction employed'')')
           else if (i_compress_correct .eq. 2) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     Zeman-BL-type dilatation-'',
     .       ''dissipation compressibility correction employed'')')
           end if
           if (isstrc .eq. 1) then
             nou(1) = min(nou(1)+1,ibufdim)
             write(bou(nou(1),1),'(''     SSTRC-type curvature'',
     .       '' correction employed (AIAA 98-2554), sstrc_crc='',f5.2)')
     .       sstrc_crc
           end if
        end if
c
        if(keepambient .eq. 1 .and. ivmx .ne. 15) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     2-eqn ambient turbulence'',
     .     '' levels not allowed to decay'')')
        end if
c
      end if
c
c   Note: (10.**(-iexp) is machine zero)
      xminn=10.**(-iexp+1)
c
      catris_kw=0.0
      if (i_catris_kw .eq. 1) catris_kw=1.0
      if (i_compress_correct .eq. 1) then
c        Wilcox-type (Turb Modeling for CFD, Ed 3, 2006, p. 258)
         turb_mach_0 = 0.25
         xsi_star = 2.0
         gam_comp = 0.0
      else if (i_compress_correct .eq. 2) then
c        Zeman B.L.-type (AIAA 93-0897)
         turb_mach_0 = 0.2
         xsi_star = 0.75
         gam_comp = 0.66
      end if
c
c Set number of subiterations to solve turbulence field eqn per iteration
c (usually, 1 is sufficient... but if residual diverges then may need more)
c
      nsubit=nsubturb
c
c Set factor that multiplies the N-S CFL number, to determine CFL number
c for the turbulence model (this typically can be 2 - 10... optimum value
c is a function of the case).  If turb model seems to be having trouble
c converging, lowering this factor may be one strategy to try. NOTE: factor
c used only for steady state cases, not time accurate cases.
c
      if(ivmx .eq. 6) factor=2.
      if(ivmx .eq. 7) factor=10.
      if(ivmx .eq. 8 .or. ivmx .eq. 12 .or. ivmx .eq. 14) factor=5.
      if(ivmx .eq. 10 .or. ivmx .eq. 11) factor=2.
      if(ivmx .eq. 9 .or. ivmx .eq. 13) factor=5.
      if(ivmx .eq. 15) factor=10.
c
c Overwrite with keyword value if nonzero
c
      if (real(cflturb).ne.0.) then
         factor = cflturb
      end if
c
c Timestep for turb model - damp(j,k,i)
c
      if (real(dt).lt.0) then
         do i=1,idim-1
         do k=1,kdim-1
         do j=1,jdim-1
           timestp(j,k,i)=factor*vol(j,k,i)/dtj(j,k,i)
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
         enddo
         enddo
         enddo
      end if
c
c Set up constants
      vk =.41
      a1 =.31
c The following will effectively override the SST denominator option:
      if (i_bsl .eq. 1) a1 = 1.e8
c Constants for Set 1:
c  cmuc1 (normally 0.09):
      cmuc1=0.09
c  constants for EASM models:
      if(ivmx.eq.8 .or. ivmx.eq.11 .or. ivmx .eq. 12 .or.
     .   ivmx.eq.9 .or. ivmx.eq.13 .or. ivmx .eq. 14) then
c     SSG:
        cmuc1=0.081
        if(ivmx .eq. 9 .or. ivmx .eq. 13) cmuc1=0.0885
        if(ivmx .eq. 8 .or. ivmx .eq. 14) cmuc1=0.0895
        c1=6.8
        c2=0.36
        c3=1.25
        c4=0.40
        c5=1.88
        gg=1./(c1/2.+c5-1.)
      end if
c  Note: for i_wilcox06, beta1, sigo1, sigk1, alp1 not used
c  (only beta2, sigo2, sigk2, alp2 are used)
c  beta1 (constant in omega or epsilon destruction term):
      if(ivmx.eq.6 .or. ivmx.eq.7)                  beta1=0.075
      if(ivmx.eq.12)                                beta1=0.83
      if(ivmx.eq.8 .or. ivmx.eq.14)                 beta1=0.83
      if(ivmx.eq.10 .or. ivmx.eq.11 .or.ivmx.eq.9
     .              .or. ivmx.eq.13)                beta1=1.83
c  sigo1 (constant in omega, epsilon, or zeta diffusion term):
      if(ivmx.eq.6 .or. ivmx.eq.7)                  sigo1=0.5
      if(ivmx.eq.12)                                sigo1=1./2.20
      if(ivmx.eq.10)                                sigo1=0.71429
      if(ivmx.eq.11)                                sigo1=1./1.3
      if(ivmx.eq.15)                                sigo1=1.46
c
c  sigk1 (constant in k diffusion term):
      if(ivmx.eq.6)                                 sigk1=0.5
      if(ivmx.eq.7 .and. i_bsl .ne. 1)              sigk1=0.85
      if(ivmx.eq.7 .and. i_bsl .eq. 1)              sigk1=0.5
      if(ivmx.eq.12)                                sigk1=1./1.4
      if(ivmx.eq.10 .or. ivmx.eq.11 .or.ivmx.eq.9
     .              .or. ivmx.eq.13)                sigk1=1.00
      if(ivmx.eq.8 .or. ivmx.eq.14)                 sigk1=1./1.0
      if(ivmx.eq.15)                                sigk1 = 1.8
c
c  alp1 (constant in omega or epsilon production term):
      if(ivmx.eq.6 .or. ivmx.eq.7)                  alp1 =
     +    beta1/cmuc1 - sigo1*vk*vk/sqrt(cmuc1)
      if(ivmx.eq.10)                                alp1 =1.45
      if(ivmx.eq.11 .or. ivmx.eq.12)                alp1 =
     +    beta1 - vk*vk*sigo1/sqrt(cmuc1)
      if(ivmx.eq.9 .or. ivmx.eq.13)                 alp1 = 1.44
      if(ivmx.eq.8 .or. ivmx.eq.14)                 alp1 = 0.53
      if(ivmx.eq.9 .or. ivmx.eq.13)                 sigo1=
     + sqrt(cmuc1)*(beta1-alp1)/(vk*vk)
      if(ivmx.eq.8 .or. ivmx.eq.14)                 sigo1=
     + sqrt(cmuc1)*(beta1-alp1)/(vk*vk)
c
c Coefficient for fnu in diffusion term
      sigkmu = 1.0
      if(ivmx.eq.15) then
        sigkmu = 1./3.
      endif
c
c Constants for Set 2 (2-layer SST model and Wilcox06 only):
      if (ivmx.eq.6 .and. i_wilcox06.eq.1) then
        cmuc2=0.09
        beta2=0.0708
        sigo2=0.5
        sigk2=0.6
        alp2 =13./25.
        sigd0=0.125
      else
        cmuc2=0.09
        beta2=0.0828
        sigo2=0.856
        sigk2=1.00
        alp2 =beta2/cmuc2 - sigo2*vk*vk/sqrt(cmuc2)
      end if
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
        negn1 = 0
        negn2 = 0
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
      if (i2d .ne. 1) then
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
c that for these particular cases it is using slightly different form
c at the block edges)
      if (ivmx .eq. 12) then
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            v3dtmp(j,k,i)=cmuc1*q(j,k,i,1)*zksav(j,k,i,2)/
     +                    zksav(j,k,i,1)
            v3dtmp(j,k,i)=ccmin(v3dtmp(j,k,i),edvislim)
          enddo
        enddo
      enddo
      else if (ivmx .eq. 11 .or. ivmx .eq. 9 .or. ivmx .eq. 13) then
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            v3dtmp(j,k,i)=cmuc1*q(j,k,i,1)*zksav(j,k,i,2)**2/
     +                    zksav(j,k,i,1)
            v3dtmp(j,k,i)=ccmin(v3dtmp(j,k,i),edvislim)
          enddo
        enddo
      enddo
      else if (ivmx .eq. 6 .and. i_wilcox06 .eq. 1) then
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            v3dtmp(j,k,i)=q(j,k,i,1)*zksav(j,k,i,2)/
     +                    zksav(j,k,i,1)
            v3dtmp(j,k,i)=ccmin(v3dtmp(j,k,i),edvislim)
          enddo
        enddo
      enddo
      else
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            v3dtmp(j,k,i)=vist3d(j,k,i)
          enddo
        enddo
      enddo
      end if
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
      if (i2d .ne. 1) then
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
c zksav2(j,k,i,1) and (j,k,i,2) are turb quantities
c zksav2(j,k,i,3) and (j,k,i,4) are Delta turb quantities
      if (real(dt) .gt. 0. .and. icyc .eq. 1) then
      if (abs(ita) .eq. 2) then
c     if zksav2 at 1st point is zero, then we know that we do not have
c     2nd order data from the restart; no choice but to set
c     zksav2(j,k,i,3&4)=deltaQ=0 for 1st iteration
        if (real(zksav2(1,1,1,1)) .eq. 0.) then
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              zksav2(j,k,i,3)=0.
              zksav2(j,k,i,4)=0.
            enddo
          enddo
        enddo
        else
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              zksav2(j,k,i,3)=zksav(j,k,i,1)-zksav2(j,k,i,1)
              zksav2(j,k,i,4)=zksav(j,k,i,2)-zksav2(j,k,i,2)
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
          enddo
        enddo
c      (2) k=kdim boundary:
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,kdim,i,1)=tk0(j,i,1,3)
            turre(j,kdim,i,2)=tk0(j,i,2,3)
          enddo
        enddo
c      (3) j=0 boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(0,k,i,1)=tj0(k,i,1,1)
            turre(0,k,i,2)=tj0(k,i,2,1)
          enddo
        enddo
c      (4) j=jdim boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(jdim,k,i,1)=tj0(k,i,1,3)
            turre(jdim,k,i,2)=tj0(k,i,2,3)
          enddo
        enddo
        if (i2d .ne. 1) then
c      (5) i=0 boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,0,1)=ti0(j,k,1,1)
            turre(j,k,0,2)=ti0(j,k,2,1)
          enddo
        enddo
c      (6) i=idim boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,idim,1)=ti0(j,k,1,3)
            turre(j,k,idim,2)=ti0(j,k,2,3)
          enddo
        enddo
        end if
        if (iturbord .ne. 1) then
c      (1) k=0 boundary:
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,-1,i,1)=tk0(j,i,1,2)
            turre(j,-1,i,2)=tk0(j,i,2,2)
          enddo
        enddo
c      (2) k=kdim boundary:
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,kdim+1,i,1)=tk0(j,i,1,4)
            turre(j,kdim+1,i,2)=tk0(j,i,2,4)
          enddo
        enddo
c      (3) j=0 boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(-1,k,i,1)=tj0(k,i,1,2)
            turre(-1,k,i,2)=tj0(k,i,2,2)
          enddo
        enddo
c      (4) j=jdim boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(jdim+1,k,i,1)=tj0(k,i,1,4)
            turre(jdim+1,k,i,2)=tj0(k,i,2,4)
          enddo
        enddo
        if (i2d .ne. 1) then
c      (5) i=0 boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,-1,1)=ti0(j,k,1,2)
            turre(j,k,-1,2)=ti0(j,k,2,2)
          enddo
        enddo
c      (6) i=idim boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,idim+1,1)=ti0(j,k,1,4)
            turre(j,k,idim+1,2)=ti0(j,k,2,4)
          enddo
        enddo
        end if
        end if
c    **************
c
c Get damp1 = CD = cross derivative term for SST or EASM k-o models only:
        if(ivmx .eq. 7 .or. ivmx .eq. 8 .or. ivmx .eq. 14) then
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
        if(i2d .ne. 1) then
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
        if (ivmx .eq. 8 .or. ivmx .eq. 14) then
          do i=1,idim-1
            do k=1,kdim-1
              do j=1,jdim-1
                temp=cmuc1*cmuc1*damp1(j,k,i)/(2.*sigo2*re*
     +               turre(j,k,i,1)**2)
                temp=ccmaxcr(temp,0.)
                damp1(j,k,i)=(1.+680.*temp*temp)/(1.+400.*temp*temp)
              enddo
            enddo
          enddo
        end if
        end if
c Get damp1 = CD = cross derivative term for Wilcox06 model only:
        if(ivmx .eq. 6 .and. i_wilcox06 .eq. 1) then
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
              xderiv=0.25*(turre(j,k+1,i,1)-turre(j,k-1,i,1))*
     +                    (turre(j,k+1,i,2)-turre(j,k-1,i,2))
              damp1(j,k,i)=xderiv*sigd0*tt/(turre(j,k,i,1)*re)
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
              xderiv=0.25*(turre(j+1,k,i,1)-turre(j-1,k,i,1))*
     +                    (turre(j+1,k,i,2)-turre(j-1,k,i,2))
              damp1(j,k,i)=damp1(j,k,i)+xderiv*sigd0*tt/
     +                     (turre(j,k,i,1)*re)
            enddo
          enddo
        enddo
        if(i2d .ne. 1) then
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
                xderiv=0.25*(turre(j,k,i+1,1)-turre(j,k,i-1,1))*
     +                      (turre(j,k,i+1,2)-turre(j,k,i-1,2))
                damp1(j,k,i)=damp1(j,k,i)+xderiv*sigd0*tt/
     +                       (turre(j,k,i,1)*re)
              enddo
            enddo
          enddo
        end if
c   for Wilcox06, damp1 only active if it is positive
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              damp1(j,k,i)=ccmaxcr(damp1(j,k,i),0.)
            enddo
          enddo
        enddo
        end if
c   get blend = F1 factor
        if (ivmx .eq. 7) then
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
            enddo
          enddo
        enddo
        else if (ivmx .eq. 6 .and. i_wilcox06 .eq. 1) then
        do k=1,kdim-1
          do i=1,idim-1
            do j=1,jdim-1
              blend(j,k,i)=0.0
            enddo
          enddo
        enddo
        else
        do k=1,kdim-1
          do i=1,idim-1
            do j=1,jdim-1
              blend(j,k,i)=1.0
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              rhop=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(q(j,k+1,i,1)+q(j,k,i,1)))
              rhom=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(q(j,k-1,i,1)+q(j,k,i,1)))
              cdp=(fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,1)=(1.0-catris_kw)*(-byy*turre(j,k-1,i,1)
     +          -cyy*turre(j,k,i,1) -dyy*turre(j,k+1,i,1)) +
     +          catris_kw*(-byy*sqrt(q(j,k-1,i,1))*turre(j,k-1,i,1)
     +          -cyy*sqrt(q(j,k,i,1))*turre(j,k,i,1) 
     +          -dyy*sqrt(q(j,k+1,i,1))*turre(j,k+1,i,1))
              rhop=(1.0-catris_kw) +
     +             catris_kw*.5*(q(j,k+1,i,1)+q(j,k,i,1))
              rhom=(1.0-catris_kw) +
     +             catris_kw*.5*(q(j,k-1,i,1)+q(j,k,i,1))
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,2)=(1.0-catris_kw)*(-byy*turre(j,k-1,i,2)
     +          -cyy*turre(j,k,i,2) -dyy*turre(j,k+1,i,2)) +
     +          catris_kw*(-byy*q(j,k-1,i,1)*turre(j,k-1,i,2)
     +          -cyy*q(j,k,i,1)*turre(j,k,i,2) 
     +          -dyy*q(j,k+1,i,1)*turre(j,k+1,i,2))
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              rhop=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(q(j,k+1,i,1)+q(j,k,i,1)))
              rhom=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(qk0(j,i,1,1)+q(j,k,i,1)))
              cdp=(fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,1)=(1.0-catris_kw)*(-byy*turre(j,k-1,i,1)
     +          -cyy*turre(j,k,i,1) -dyy*turre(j,k+1,i,1)) +
     +          catris_kw*(-byy*sqrt(qk0(j,i,1,1))*turre(j,k-1,i,1)
     +          -cyy*sqrt(q(j,k,i,1))*turre(j,k,i,1)
     +          -dyy*sqrt(q(j,k+1,i,1))*turre(j,k+1,i,1))
              rhop=(1.0-catris_kw) +
     +             catris_kw*.5*(q(j,k+1,i,1)+q(j,k,i,1))
              rhom=(1.0-catris_kw) +
     +             catris_kw*.5*(qk0(j,i,1,1)+q(j,k,i,1))
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,2)=(1.0-catris_kw)*(-byy*turre(j,k-1,i,2)
     +          -cyy*turre(j,k,i,2) -dyy*turre(j,k+1,i,2)) +
     +          catris_kw*(-byy*qk0(j,i,1,1)*turre(j,k-1,i,2)
     +          -cyy*q(j,k,i,1)*turre(j,k,i,2)
     +          -dyy*q(j,k+1,i,1)*turre(j,k+1,i,2))
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k+1,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k-1,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              rhop=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(qk0(j,i,1,3)+q(j,k,i,1)))
              rhom=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(q(j,k-1,i,1)+q(j,k,i,1)))
              cdp=(fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,1)=(1.0-catris_kw)*(-byy*turre(j,k-1,i,1)
     +          -cyy*turre(j,k,i,1) -dyy*turre(j,k+1,i,1)) +
     +          catris_kw*(-byy*sqrt(q(j,k-1,i,1))*turre(j,k-1,i,1)
     +          -cyy*sqrt(q(j,k,i,1))*turre(j,k,i,1)
     +          -dyy*sqrt(qk0(j,i,1,3))*turre(j,k+1,i,1))
              rhop=(1.0-catris_kw) +
     +             catris_kw*.5*(qk0(j,i,1,3)+q(j,k,i,1))
              rhom=(1.0-catris_kw) +
     +             catris_kw*.5*(q(j,k-1,i,1)+q(j,k,i,1))
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              byy=-cdm
              cyy= cdp+cdm
              dyy=-cdp
              rhside(j,k,i,2)=(1.0-catris_kw)*(-byy*turre(j,k-1,i,2)
     +          -cyy*turre(j,k,i,2) -dyy*turre(j,k+1,i,2)) +
     +          catris_kw*(-byy*q(j,k-1,i,1)*turre(j,k-1,i,1)
     +          -cyy*q(j,k,i,1)*turre(j,k,i,1)
     +          -dyy*qk0(j,i,1,3)*turre(j,k+1,i,1))
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              rhop=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(q(j+1,k,i,1)+q(j,k,i,1)))
              rhom=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(q(j-1,k,i,1)+q(j,k,i,1)))
              cdp=(fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,1)=rhside(j,k,i,1) +
     +          (1.0-catris_kw)*(-bxx*turre(j-1,k,i,1)
     +          -cxx*turre(j,k,i,1) -dxx*turre(j+1,k,i,1)) +
     +          catris_kw*(-bxx*sqrt(q(j-1,k,i,1))*turre(j-1,k,i,1)
     +          -cxx*sqrt(q(j,k,i,1))*turre(j,k,i,1) 
     +          -dxx*sqrt(q(j+1,k,i,1))*turre(j+1,k,i,1))
              rhop=(1.0-catris_kw) +
     +             catris_kw*.5*(q(j+1,k,i,1)+q(j,k,i,1))
              rhom=(1.0-catris_kw) +
     +             catris_kw*.5*(q(j-1,k,i,1)+q(j,k,i,1))
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,2)=rhside(j,k,i,2) +
     +          (1.0-catris_kw)*(-bxx*turre(j-1,k,i,2)
     +          -cxx*turre(j,k,i,2) -dxx*turre(j+1,k,i,2)) +
     +          catris_kw*(-bxx*q(j-1,k,i,1)*turre(j-1,k,i,2)
     +          -cxx*q(j,k,i,1)*turre(j,k,i,2) 
     +          -dxx*q(j+1,k,i,1)*turre(j+1,k,i,2))
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              rhop=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(q(j+1,k,i,1)+q(j,k,i,1)))
              rhom=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(qj0(k,i,1,1)+q(j,k,i,1)))
              cdp=(fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,1)=rhside(j,k,i,1) +
     +          (1.0-catris_kw)*(-bxx*turre(j-1,k,i,1)
     +          -cxx*turre(j,k,i,1) -dxx*turre(j+1,k,i,1)) +
     +          catris_kw*(-bxx*sqrt(qj0(k,i,1,1))*turre(j-1,k,i,1)
     +          -cxx*sqrt(q(j,k,i,1))*turre(j,k,i,1)
     +          -dxx*sqrt(q(j+1,k,i,1))*turre(j+1,k,i,1))
              rhop=(1.0-catris_kw) +
     +             catris_kw*.5*(q(j+1,k,i,1)+q(j,k,i,1))
              rhom=(1.0-catris_kw) +
     +             catris_kw*.5*(qj0(k,i,1,1)+q(j,k,i,1))
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,2)=rhside(j,k,i,2) +
     +          (1.0-catris_kw)*(-bxx*turre(j-1,k,i,2)
     +          -cxx*turre(j,k,i,2) -dxx*turre(j+1,k,i,2)) +
     +          catris_kw*(-bxx*qj0(k,i,1,1)*turre(j-1,k,i,1)
     +          -cxx*q(j,k,i,1)*turre(j,k,i,1)
     +          -dxx*q(j+1,k,i,1)*turre(j+1,k,i,1))
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j+1,k,i))
              anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j-1,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              rhop=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(qj0(k,i,1,3)+q(j,k,i,1)))
              rhom=(1.0-catris_kw) +
     +             catris_kw*sqrt(.5*(q(j-1,k,i,1)+q(j,k,i,1)))
              cdp=(fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,1)=rhside(j,k,i,1) +
     +          (1.0-catris_kw)*(-bxx*turre(j-1,k,i,1)
     +          -cxx*turre(j,k,i,1) -dxx*turre(j+1,k,i,1)) +
     +          catris_kw*(-bxx*sqrt(q(j-1,k,i,1))*turre(j-1,k,i,1)
     +          -cxx*sqrt(q(j,k,i,1))*turre(j,k,i,1)
     +          -dxx*sqrt(qj0(k,i,1,3))*turre(j+1,k,i,1))
              rhop=(1.0-catris_kw) +
     +             catris_kw*.5*(qj0(k,i,1,3)+q(j,k,i,1))
              rhom=(1.0-catris_kw) +
     +             catris_kw*.5*(q(j-1,k,i,1)+q(j,k,i,1))
              cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re)
              cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re)
              bxx=-cdm
              cxx= cdp+cdm
              dxx=-cdp
              rhside(j,k,i,2)=rhside(j,k,i,2) +
     +          (1.0-catris_kw)*(-bxx*turre(j-1,k,i,2)
     +          -cxx*turre(j,k,i,2) -dxx*turre(j+1,k,i,2)) +
     +          catris_kw*(-bxx*q(j-1,k,i,1)*turre(j-1,k,i,1)
     +          -cxx*q(j,k,i,1)*turre(j,k,i,1)
     +          -dxx*qj0(k,i,1,3)*turre(j+1,k,i,1))
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
            enddo
          enddo
        enddo
        end if
c
c    F_zeta_zeta viscous terms
        if(i2d .ne. 1) then
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
                anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1))
                anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                rhop=(1.0-catris_kw) +
     +               catris_kw*sqrt(.5*(q(j,k,i+1,1)+q(j,k,i,1)))
                rhom=(1.0-catris_kw) +
     +               catris_kw*sqrt(.5*(q(j,k,i-1,1)+q(j,k,i,1)))
                cdp=(fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re)
                cdm=(fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,1)=rhside(j,k,i,1) +
     +            (1.0-catris_kw)*(-bzz*turre(j,k,i-1,1)
     +            -czz*turre(j,k,i,1) -dzz*turre(j,k,i+1,1)) +
     +            catris_kw*(-bzz*sqrt(q(j,k,i-1,1))*turre(j,k,i-1,1)
     +            -czz*sqrt(q(j,k,i,1))*turre(j,k,i,1)
     +            -dzz*sqrt(q(j,k,i+1,1))*turre(j,k,i+1,1))
                rhop=(1.0-catris_kw) +
     +               catris_kw*.5*(q(j,k,i+1,1)+q(j,k,i,1))
                rhom=(1.0-catris_kw) +
     +               catris_kw*.5*(q(j,k,i-1,1)+q(j,k,i,1))
                cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,2)=rhside(j,k,i,2) +
     +            (1.0-catris_kw)*(-bzz*turre(j,k,i-1,2)
     +            -czz*turre(j,k,i,2) -dzz*turre(j,k,i+1,2)) +
     +            catris_kw*(-bzz*q(j,k,i-1,1)*turre(j,k,i-1,1)
     +            -czz*q(j,k,i,1)*turre(j,k,i,1)
     +            -dzz*q(j,k,i+1,1)*turre(j,k,i+1,1))
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
                anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1))
                anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                rhop=(1.0-catris_kw) +
     +               catris_kw*sqrt(.5*(q(j,k,iu,1)+q(j,k,i,1)))
                rhom=(1.0-catris_kw) +
     +               catris_kw*sqrt(.5*(qi0(j,k,1,1)+q(j,k,i,1)))
                cdp=(fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re)
                cdm=(fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,1)=rhside(j,k,i,1) +
     +            (1.0-catris_kw)*(-bzz*turre(j,k,i-1,1)
     +            -czz*turre(j,k,i,1) -dzz*turre(j,k,i+1,1)) +
     +            catris_kw*(-bzz*sqrt(qi0(j,k,1,1))*turre(j,k,i-1,1)
     +            -czz*sqrt(q(j,k,i,1))*turre(j,k,i,1)
     +            -dzz*sqrt(q(j,k,iu,1))*turre(j,k,i+1,1))
                rhop=(1.0-catris_kw) +
     +               catris_kw*.5*(q(j,k,iu,1)+q(j,k,i,1))
                rhom=(1.0-catris_kw) +
     +               catris_kw*.5*(qi0(j,k,1,1)+q(j,k,i,1))
                cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,2)=rhside(j,k,i,2) +
     +            (1.0-catris_kw)*(-bzz*turre(j,k,i-1,2)
     +            -czz*turre(j,k,i,2) -dzz*turre(j,k,i+1,2)) +
     +            catris_kw*(-bzz*qi0(j,k,1,1)*turre(j,k,i-1,1)
     +            -czz*q(j,k,i,1)*turre(j,k,i,1)
     +            -dzz*q(j,k,iu,1)*turre(j,k,i+1,1))
              enddo
            enddo
c
c         IDIM boundary points 
            i=idim-1
            il=max(1,idim-2)
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
                anutp=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i+1))
                anutm=.5*(v3dtmp(j,k,i)+v3dtmp(j,k,i-1))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                rhop=(1.0-catris_kw) +
     +               catris_kw*sqrt(.5*(qi0(j,k,1,3)+q(j,k,i,1)))
                rhom=(1.0-catris_kw) +
     +               catris_kw*sqrt(.5*(q(j,k,il,1)+q(j,k,i,1)))
                cdp=(fnup+sigop*anutp)*ttp/(rhop*q(j,k,i,1)*re)
                cdm=(fnum+sigom*anutm)*ttm/(rhom*q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,1)=rhside(j,k,i,1) +
     +            (1.0-catris_kw)*(-bzz*turre(j,k,i-1,1)
     +            -czz*turre(j,k,i,1) -dzz*turre(j,k,i+1,1)) +
     +            catris_kw*(-bzz*sqrt(q(j,k,il,1))*turre(j,k,i-1,1)
     +            -czz*sqrt(q(j,k,i,1))*turre(j,k,i,1)
     +            -dzz*sqrt(qi0(j,k,1,3))*turre(j,k,i+1,1))
                rhop=(1.0-catris_kw) +
     +               catris_kw*.5*(qi0(j,k,1,3)+q(j,k,i,1))
                rhom=(1.0-catris_kw) +
     +               catris_kw*.5*(q(j,k,il,1)+q(j,k,i,1))
                cdp=(sigkmu*fnup+sigkp*anutp)*ttp/(rhop*q(j,k,i,1)*re)
                cdm=(sigkmu*fnum+sigkm*anutm)*ttm/(rhom*q(j,k,i,1)*re)
                bzz=-cdm
                czz= cdp+cdm
                dzz=-cdp
                rhside(j,k,i,2)=rhside(j,k,i,2) +
     +            (1.0-catris_kw)*(-bzz*turre(j,k,i-1,2)
     +            -czz*turre(j,k,i,2) -dzz*turre(j,k,i+1,2)) +
     +            catris_kw*(-bzz*q(j,k,il,1)*turre(j,k,i-1,1)
     +            -czz*q(j,k,i,1)*turre(j,k,i,1)
     +            -dzz*qi0(j,k,1,3)*turre(j,k,i+1,1))
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
              enddo
            enddo
          enddo
          end if
        end if
c
c   Curvature terms for EASMCC
        if ((ivmx.eq.8 .or. ivmx.eq.9 .or. ivmx.eq.11 .or. 
     .       ivmx.eq.12.or. ivmx.eq.13.or. ivmx.eq.14) .and. 
     .       ieasmcc2d .eq. 1) then
c   Note EASMCC currently only set up for curvature in x-z plane, with
c   y direction (necessarily) in the i-direction
c   Get vx(3)=DS11/Dt and vx(4)=DS13/Dt:
        call sijrate2d(idim,jdim,kdim,q,qj0,qk0,
     .  bcj,bck,vol,sj,sk,vx)
c   compute curvature term Dalpha/Dt: store in vx(1)
c   modify Sij such that its diagonal terms are traceless in the 2D sense
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
              s11 = ux(j,k,i,1)
     +           -(ux(j,k,i,1)+ux(j,k,i,9))/2.
              s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
c             denom=2.*(s11**2+s13**2)+
c    +              (xmach*reue/2.3e6)**2
              denom=2.*(s11**2+s13**2)
              denom=ccmax(denom,xminn)
              vx(j,k,i,1)=(s11*vx(j,k,i,4)-
     +         s13*vx(j,k,i,3))/denom
            enddo
          enddo
        enddo
        end if
c    * * * * * * * * * * * *
c    ADD SOURCE TERMS TO RHS
c    * * * * * * * * * * * *
      if(isklton .gt. 0) then
        if(ilamlo.eq.0 .or. jlamlo.eq.0 .or. klamlo.eq.0) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'('' Block #'',i5,'' in 2-eqn turb model'',
     .   '' has no laminar regions'')') nbl
        else
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'('' Block #'',i5,'' in 2-eqn turb model'',
     .   '' - laminar region is:'')') nbl
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'('' i='',i5,'' to'',i5,'', j='',i5,'' to'',
     .   i5,'', k='',i5,'' to'',i5)') ilamlo,ilamhi,jlamlo,jlamhi,
     .   klamlo,klamhi
        end if
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'('' NOTE:  This particular model'',
     .   '' <<transitions>> on its own, but there is'')')
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'('' no guarantee that it will transition'',
     .   '' at all.  Check vist3d levels if unsure.'')')
      end if
c  Source terms for k-omega (Wilcox & SST):
      if (ivmx .eq. 6 .or. ivmx .eq. 7) then
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
c            Determine Sij etc if needed
c            (Note if want traceless must set prod2d3dtrace=0.5 (2D)
c            or 0.33333 (3D))
              if ((ivmx .eq. 6 .and. i_wilcox06 .eq. 1 .and.
     +            i_wilcox06_chiw .eq. 1) .or. ikoprod .eq. 1
     +            .or. isstrc .eq. 1) then
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
     +                2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
                w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
                w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
                w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
              end if
c  the following modifies betax via the vortex-stretching parameter
c  (default is currently i_wilcox06_chiw = 1)
              if (ivmx .eq. 6 .and. i_wilcox06 .eq. 1 .and.
     +            i_wilcox06_chiw .eq. 1) then
c  in 2-D with x-z plane, the following reduces to: 
c  -w13*w13*[(s11-trace2d)+(s33-trace2d)] = 0
c  for this term, the "tracepart" is always fixed at 0.5*dui/dxi
                tracepart=0.5*(s11+s22+s33)
                xnum=-(w12*w12+w13*w13)*(s11-tracepart)
     +               -(w12*w12+w23*w23)*(s22-tracepart)
     +               -(w13*w13+w23*w23)*(s33-tracepart) - 2.*w13*w23*s12
     +               +2.*w12*w23*s13 - 2.*w12*w13*s23
                chiw=abs(xnum/((re*beta2*turre(j,k,i,1))**3))
                fbeta=(1.+85.*chiw)/(1.+100.*chiw)
                betax=betax*fbeta
              end if
              if (i_compress_correct .eq. 1) then
                turb_mach2 = 2.*q(j,k,i,1)*turre(j,k,i,2)/
     +                      (gamma*q(j,k,i,5))
                f_turb_mach = ccmaxcr(turb_mach2-turb_mach_0**2,0.0)
                betax=betax-cmuc*xsi_star*f_turb_mach
                cmuc=cmuc*(1.+xsi_star*f_turb_mach)
              else if (i_compress_correct .eq. 2) then
                turb_mach = sqrt(2.*q(j,k,i,1)*turre(j,k,i,2)/
     +                      (gamma*q(j,k,i,5)))
                zeman_lag = ccmaxcr(turb_mach-turb_mach_0,0.0)
                f_turb_mach = 1.0-(exp(-(zeman_lag/gam_comp)**2))
                betax=betax-cmuc*xsi_star*f_turb_mach
                cmuc=cmuc*(1.+xsi_star*f_turb_mach)
              end if
              f4 = 1.0
              if (isstrc .eq. 1) then
                sij=sqrt(2.0*xis) + 1.e-20
                ri=(vor(j,k,i)/sij)*(vor(j,k,i)/sij-1.0)
                f4=1.0/(1.0+sstrc_crc*ri)
c               f4=blend(j,k,i)*f4 + 1.0 - blend(j,k,i)
              endif
              if (ikoprod .eq. 1) then
                rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp/re*
     +           2.*xis - f4*re*betax*turre(j,k,i,1)**2
     +           +(1.-blend(j,k,i))*damp1(j,k,i)
                pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*2.*xis
c   if desired, subtract off 2/3k term
c   (can only do this if ikoprod=1):
                if (i_turbprod_kterm .eq. 1) then
                  rhside(j,k,i,1)=rhside(j,k,i,1)-(cutoff*2./3.*alp/
     +              vist3d(j,k,i)*q(j,k,i,1)*turre(j,k,i,2)*
     +              (s11+s22+s33))
                  pk=pk-(2./3.*turre(j,k,i,2)*(s11+s22+s33))
                end if
              else
                rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp/re*
     +          vor(j,k,i)**2 - f4*re*betax*turre(j,k,i,1)**2
     +          +(1.-blend(j,k,i))*damp1(j,k,i)
                pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)**2
              end if
              dk=re*cmuc*turre(j,k,i,1)*turre(j,k,i,2)
              pk=ccmin(pk,(pklimterm*dk))
              rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk-dk
              rhside(j,k,i,1)=rhside(j,k,i,1) + keepambient*
     +         re*betax*tur10(1)*tur10(1)
              rhside(j,k,i,2)=rhside(j,k,i,2) + keepambient*
     +         re*cmuc*tur10(1)*tur10(2)
            enddo
          enddo
        enddo
c  Source terms for nonlinear EASM k-omega:
      else if (ivmx .eq. 12) then
c      "Constants":
        alpa1=(2.-c4)/2.*gg
        alpa2=(2.-c3)*gg
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
c            Determine Sij and Wij values:
              s11 = ux(j,k,i,1)
              s22 = ux(j,k,i,5)
              s33 = ux(j,k,i,9)
              s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
              s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
              s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
              w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
              w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
              w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
c            Limit omega here:
              omegatemp=ccmax(turre(j,k,i,1),tur10(1))
c     Get "factre" term:  changes numerator of cmu for nonlinear terms only
              eta=(2.-c3)**2*(gg*gg/4.)*xis/(omegatemp*re)**2
              squig=(2.-c4)**2*(gg*gg/4.)*wis/(omegatemp*re)**2
              eta=ccmincr(eta,10.)
              squig=ccmincr(squig,10.)
              factre=(3.*(1.+eta)+.2e-8*(eta*eta*eta +
     +          squig*squig*squig))/
     +               (3.*(1.+eta)+   .2*(eta*eta*eta +
     +          squig*squig*squig))
c            Find tauij values:
              t11 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre/
     +                 omegatemp*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              t11=ccmaxcr(t11,0.)
              t22 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre/
     +                 omegatemp*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              t22=ccmaxcr(t22,0.)
              t33 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre/
     +                 omegatemp*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              t33=ccmaxcr(t33,0.)
              t12 =-2.*vist3d(j,k,i)*s12
     +           -2.*alpa1*vist3d(j,k,i)*factre/omegatemp*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*vist3d(j,k,i)*s13
     +           -2.*alpa1*vist3d(j,k,i)*factre/omegatemp*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*vist3d(j,k,i)*s23
     +           -2.*alpa1*vist3d(j,k,i)*factre/omegatemp*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre/omegatemp*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
c            Calculate production term:
              pk = -(t11*ux(j,k,i,1) + t22*ux(j,k,i,5)
     +              +t33*ux(j,k,i,9) + t12*(ux(j,k,i,2)+ux(j,k,i,4))
     +                               + t13*(ux(j,k,i,3)+ux(j,k,i,7))
     +                               + t23*(ux(j,k,i,6)+ux(j,k,i,8)))/
     +              (q(j,k,i,1)*re)
              pk=ccabs(pk)
c---clr
c   alternate 
c             pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)**2
c---clr
c            Add to RHS:
              rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*pk*
     +         turre(j,k,i,1)/turre(j,k,i,2) -
     +         re*beta1*turre(j,k,i,1)**2
              dk = re*turre(j,k,i,1)*turre(j,k,i,2)
              pk = ccmin(pk,(pklimterm*dk))
              rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk - dk
              rhside(j,k,i,1)=rhside(j,k,i,1) + keepambient*
     +         re*beta1*tur10(1)*tur10(1)
              rhside(j,k,i,2)=rhside(j,k,i,2) + keepambient*
     +         re*tur10(1)*tur10(2)
            enddo
          enddo
        enddo
c  Source terms for k-epsilon (Abid):
      else if (ivmx .eq.10) then
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
              rek=q(j,k,i,1)*sqrt(turre(j,k,i,2))*re*ccabs(smin(j,k,i))/
     +            fnu(j,k,i)
              f2=(1.-exp(-rek/12.))
c            Store quantities to be added to certain implicit LHS terms:
              damp1(j,k,i)=1.5*re*beta1*turre(j,k,i,1)*f2/turre(j,k,i,2)
c            Add to RHS:
              if (ikoprod .eq. 1) then
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
              rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*
     +         turre(j,k,i,1)*
     +         vist3d(j,k,i)*2.*xis/(re*turre(j,k,i,2)*
     +         q(j,k,i,1)) - re*beta1*(turre(j,k,i,1)**2)*f2/
     +         turre(j,k,i,2)
              pk = vist3d(j,k,i)/(q(j,k,i,1)*re)*2.*xis
              else
              rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*
     +         turre(j,k,i,1)*
     +         vist3d(j,k,i)*vor(j,k,i)**2/(re*turre(j,k,i,2)*
     +         q(j,k,i,1)) - re*beta1*(turre(j,k,i,1)**2)*f2/
     +         turre(j,k,i,2)
              pk = vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)**2
              end if
              dk = re*turre(j,k,i,1)
              pk = ccmin(pk,(pklimterm*dk))
              rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk - dk
              rhside(j,k,i,1)=rhside(j,k,i,1) + keepambient*
     +         re*beta1*tur10(1)*tur10(1)/tur10(2)
              rhside(j,k,i,2)=rhside(j,k,i,2) + keepambient*
     +         re*tur10(1)
            enddo
          enddo
        enddo
c  Source terms for nonlinear EASM k-epsilon:
      else if (ivmx .eq.11) then
c      "Constants":
        alpa1=(2.-c4)/2.*gg
        alpa2=(2.-c3)*gg
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
              rek=q(j,k,i,1)*sqrt(turre(j,k,i,2))*re*ccabs(smin(j,k,i))/
     +            fnu(j,k,i)
              f2=(1.-exp(-rek/12.))
c            Store quantities to be added to certain implicit LHS terms:
              damp1(j,k,i)=1.5*re*beta1*turre(j,k,i,1)*f2/turre(j,k,i,2)
c            Determine Sij and Wij values:
              s11 = ux(j,k,i,1)
              s22 = ux(j,k,i,5)
              s33 = ux(j,k,i,9)
              s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
              s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
              s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
              w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
              w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
              w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
c     Get "factre" term:  changes numerator of cmu for nonlinear terms only
              eta=(2.-c3)**2*(gg*gg/4.)*xis*
     +            turre(j,k,i,2)**2/(turre(j,k,i,1)*re)**2
              squig=(2.-c4)**2*(gg*gg/4.)*wis*
     +            turre(j,k,i,2)**2/(turre(j,k,i,1)*re)**2
              eta=ccmincr(eta,10.)
              squig=ccmincr(squig,10.)
              factre=(3.*(1.+eta)+.2e-8*(eta*eta*eta +
     +          squig*squig*squig))/
     +               (3.*(1.+eta)+   .2*(eta*eta*eta +
     +          squig*squig*squig))
c            Find nonlinear tauij values:
              t11 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              t11=ccmaxcr(t11,0.)
              t22 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              t22=ccmaxcr(t22,0.)
              t33 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              t33=ccmaxcr(t33,0.)
              t12 =-2.*vist3d(j,k,i)*s12
     +           -2.*alpa1*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*vist3d(j,k,i)*s13
     +           -2.*alpa1*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*vist3d(j,k,i)*s23
     +           -2.*alpa1*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*factre*
     +                 turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
c            Calculate production term:
              pk = -(t11*ux(j,k,i,1) + t22*ux(j,k,i,5)
     +              +t33*ux(j,k,i,9) + t12*(ux(j,k,i,2)+ux(j,k,i,4))
     +                               + t13*(ux(j,k,i,3)+ux(j,k,i,7))
     +                               + t23*(ux(j,k,i,6)+ux(j,k,i,8)))/
     +              (q(j,k,i,1)*re)
              pk=ccabs(pk)
c---clr
c   alternate
c             pk=vist3d(j,k,i)/(q(j,k,i,1)*re)*vor(j,k,i)**2
c---clr
c            Add to RHS:
              rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*pk*
     +         turre(j,k,i,1)/turre(j,k,i,2) -
     +         re*beta1*(turre(j,k,i,1)**2)*f2/
     +         turre(j,k,i,2)
              dk = re*turre(j,k,i,1)
              pk = ccmin(pk,(pklimterm*dk))
              rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk - dk
              rhside(j,k,i,1)=rhside(j,k,i,1) + keepambient*
     +         re*beta1*tur10(1)*tur10(1)/tur10(2)
              rhside(j,k,i,2)=rhside(j,k,i,2) + keepambient*
     +         re*tur10(1)
            enddo
          enddo
        enddo
c  Source terms for nonlinear EASM k-epsilon:
      else if (ivmx .eq. 9 .or. ivmx .eq. 13) then
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
              rek=q(j,k,i,1)*sqrt(turre(j,k,i,2))*re*ccabs(smin(j,k,i))/
     +            fnu(j,k,i)
              f2=(1.-exp(-rek/10.8))
c            Store quantities to be added to certain implicit LHS terms:
              damp1(j,k,i)=1.5*re*beta1*turre(j,k,i,1)*f2/turre(j,k,i,2)
c            Determine Sij*Sij:
              s11 = ux(j,k,i,1)
              s22 = ux(j,k,i,5)
              s33 = ux(j,k,i,9)
              s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
              s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
              s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
c            Calculate production term:
c              Approximate term (exact for 2-D incompressible flow):
              if (iturbprod .eq. 0) then
              pk=2.*vist3d(j,k,i)*xis/(q(j,k,i,1)*re)
              else
c              Exact term:
              al10 = 0.25*c1 -1.
              al1 = 3.8
c   2-line modification to improve center of WAKE:
              al10=al10+1.8864
              al1=al1-2.
              al2 = 0.5*c2 -2./3.
              al3 = 0.5*c3 -1.
              al4 = 0.5*c4 -1.
              w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
              w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
              w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
              eta1=xis*turre(j,k,i,2)**2/(turre(j,k,i,1)*re)**2
              alpa1 = -al4/(al10-eta1*al1*cmuv(j,k,i))
              alpa2 = -2.*al3/(al10-eta1*cmuv(j,k,i)*al1)
c            Find nonlinear tauij values:
              t11 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              t11=ccmaxcr(t11,0.)
              t22 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              t22=ccmaxcr(t22,0.)
              t33 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/
     +                 turre(j,k,i,1)*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              t33=ccmaxcr(t33,0.)
              t12 =-2.*vist3d(j,k,i)*s12
     +           -2.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*vist3d(j,k,i)*s13
     +           -2.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*vist3d(j,k,i)*s23
     +           -2.*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
              pk = -(t11*ux(j,k,i,1) + t22*ux(j,k,i,5)
     +              +t33*ux(j,k,i,9) + t12*(ux(j,k,i,2)+ux(j,k,i,4))
     +                               + t13*(ux(j,k,i,3)+ux(j,k,i,7))
     +                               + t23*(ux(j,k,i,6)+ux(j,k,i,8)))/
     +              (q(j,k,i,1)*re)
              pk=ccabs(pk)
              end if
c            Add to RHS:
              rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*pk*
     +         turre(j,k,i,1)/turre(j,k,i,2) -
     +         re*beta1*(turre(j,k,i,1)**2)*f2/
     +         turre(j,k,i,2)
              dk = re*turre(j,k,i,1)
              pk = ccmin(pk,(pklimterm*dk))
              rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk - dk
              rhside(j,k,i,1)=rhside(j,k,i,1) + keepambient*
     +         re*beta1*tur10(1)*tur10(1)/tur10(2)
              rhside(j,k,i,2)=rhside(j,k,i,2) + keepambient*
     +         re*tur10(1)
            enddo
          enddo
        enddo
c  Source terms for nonlinear EASM var g k-omega:
      else if (ivmx .eq. 8 .or. ivmx .eq. 14) then
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
c            Determine Sij*Sij:
              s11 = ux(j,k,i,1)
              s22 = ux(j,k,i,5)
              s33 = ux(j,k,i,9)
              s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
              s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
              s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
c            Calculate production term:
c              Approximate term (exact for 2-D incompressible flow):
              if (iturbprod .eq. 0) then
              pk=2.*vist3d(j,k,i)*xis/(q(j,k,i,1)*re)
              else
c              Exact term:
              al10 = 0.25*c1 -1.
              al1 = 3.8
c   2-line modification to improve center of WAKE:
              al10=al10+1.8864
              al1=al1-2.
              al2 = 0.5*c2 -2./3.
              al3 = 0.5*c3 -1.
              al4 = 0.5*c4 -1.
              w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
              w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
              w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
              eta1=xis/(turre(j,k,i,1)*re)**2
              alpa1 = -al4/(al10-eta1*al1*cmuv(j,k,i))
              alpa2 = -2.*al3/(al10-eta1*cmuv(j,k,i)*al1)
c            Find nonlinear tauij values:
              t11 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)/
     +                 turre(j,k,i,1)*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              t11=ccmaxcr(t11,0.)
              t22 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)/
     +                 turre(j,k,i,1)*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              t22=ccmaxcr(t22,0.)
              t33 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
     +             -4.*alpa1*vist3d(j,k,i)/
     +                 turre(j,k,i,1)*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              t33=ccmaxcr(t33,0.)
              t12 =-2.*vist3d(j,k,i)*s12
     +           -2.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*vist3d(j,k,i)*s13
     +           -2.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*vist3d(j,k,i)*s23
     +           -2.*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
              pk = -(t11*ux(j,k,i,1) + t22*ux(j,k,i,5)
     +              +t33*ux(j,k,i,9) + t12*(ux(j,k,i,2)+ux(j,k,i,4))
     +                               + t13*(ux(j,k,i,3)+ux(j,k,i,7))
     +                               + t23*(ux(j,k,i,6)+ux(j,k,i,8)))/
     +              (q(j,k,i,1)*re)
              pk=ccabs(pk)
              end if
              fbeta=1.
              rhside(j,k,i,1)=rhside(j,k,i,1)+cutoff*alp1*pk*
     +         turre(j,k,i,1)/turre(j,k,i,2) -
     +         fbeta*re*beta1*turre(j,k,i,1)**2
              dk = damp1(j,k,i)*re*turre(j,k,i,1)*turre(j,k,i,2)
              pk = ccmin(pk,(pklimterm*dk))
              rhside(j,k,i,2)=rhside(j,k,i,2)+cutoff*pk - dk
              rhside(j,k,i,1)=rhside(j,k,i,1) + keepambient*
     +         re*beta1*tur10(1)*tur10(1)
              rhside(j,k,i,2)=rhside(j,k,i,2) + keepambient*
     +         re*tur10(1)*tur10(2)
            enddo
          enddo
        enddo
c
      else if (ivmx .eq. 15) then
c ---- enstrophy equation additions
c      see AIAA J, Vol 36, No 10, Oct 1998, p. 1825-1833
c
c ---- define constants
c
        cmu = 0.09
        rcmu = 1./cmu
        alpha3 = 0.35  
        beta4 = 0.42
        beta5 = 2.37
        beta6 = 0.10
        beta7 = 1.50
c   beta8=0 (default) often required to prevent numerical oscillations
        if (ibeta8kzeta .ne. 0) then
          beta8 = 2.3
        else
          beta8 = 0.0
        end if
        sigr = 0.07
        c1 = 0.6
        czeta1 = 2.1
        ck = 2.0
        sigp = 0.13
        sigrho = 91.9
        fnuef = 198.6/tinf
        delta = 0.1
c
        wt1 = 1.0
        wt2 = (1.-wt1)
c
        wlim = 1.e-6  
        slim = 1.e-6  
c
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              w1(j,k,i) = ux(j,k,i,8) - ux(j,k,i,6)
              w2(j,k,i) = ux(j,k,i,3) - ux(j,k,i,7)
              w3(j,k,i) = ux(j,k,i,4) - ux(j,k,i,2)
            end do
          end do
        end do
c
        do i=1,idim-1
          ip=1
          im=1
          if(i.eq.1)im=0
          if(i.eq.idim-1)ip=0
          if(idim-1.eq.1) then
            ti=1.
          else
            ti=1./float(ip+im)
          end if
          do k=1,kdim-1
            kp=1
            km=1
            if(k.eq.1)km=0
            if(k.eq.kdim-1)kp=0
            if(kdim-1.eq.1) then
              tk=1.
            else
              tk=1./float(kp+km)
            end if
c
c ---- compute gradients of vorticity and turbulence variables
c
            j=1
            dkdj(j)=(turre(j+1,k,i,2)-turre(j,k,i,2))
            dzdj(j)=(turre(j+1,k,i,1)-turre(j,k,i,1))
            dw1dj(j)=(w1(j+1,k,i)-w1(j,k,i))
            dw2dj(j)=(w2(j+1,k,i)-w2(j,k,i))
            dw3dj(j)=(w3(j+1,k,i)-w3(j,k,i))
            dpdj(j)=(q(j+1,k,i,5)-q(j,k,i,5))
            drdj(j)=(q(j+1,k,i,1)-q(j,k,i,1))
c
            do j=2,jdim-2
              dkdj(j)=0.5*(turre(j+1,k,i,2)-turre(j-1,k,i,2))
              dzdj(j)=0.5*(turre(j+1,k,i,1)-turre(j-1,k,i,1))
              dw1dj(j)=0.5*(w1(j+1,k,i)-w1(j-1,k,i))
              dw2dj(j)=0.5*(w2(j+1,k,i)-w2(j-1,k,i))
              dw3dj(j)=0.5*(w3(j+1,k,i)-w3(j-1,k,i))
              dpdj(j)=0.5*(q(j+1,k,i,5)-q(j-1,k,i,5))
              drdj(j)=0.5*(q(j+1,k,i,1)-q(j-1,k,i,1))
            enddo
c
            j=jdim-1
            dkdj(j)=(turre(j,k,i,2)-turre(j-1,k,i,2))
            dzdj(j)=(turre(j,k,i,1)-turre(j-1,k,i,1))
            dw1dj(j)=(w1(j,k,i)-w1(j-1,k,i))
            dw2dj(j)=(w2(j,k,i)-w2(j-1,k,i))
            dw3dj(j)=(w3(j,k,i)-w3(j-1,k,i))
            dpdj(j)=(q(j,k,i,5)-q(j-1,k,i,5))
            drdj(j)=(q(j,k,i,1)-q(j-1,k,i,1))
c
            do j=1,jdim-1
c
              dkdk(j)=tk*(turre(j,k+kp,i,2)-turre(j,k-km,i,2))
              dkdi(j)=ti*(turre(j,k,i+ip,2)-turre(j,k,i-im,2))
c
              dzdk(j)=tk*(turre(j,k+kp,i,1)-turre(j,k-km,i,1))
              dzdi(j)=ti*(turre(j,k,i+ip,1)-turre(j,k,i-im,1))
c
              dw1dk(j)=tk*(w1(j,k+kp,i)-w1(j,k-km,i))
              dw1di(j)=ti*(w1(j,k,i+ip)-w1(j,k,i-im)) 
c
              dw2dk(j)=tk*(w2(j,k+kp,i)-w2(j,k-km,i))
              dw2di(j)=ti*(w2(j,k,i+ip)-w2(j,k,i-im)) 
c
              dw3dk(j)=tk*(w3(j,k+kp,i)-w3(j,k-km,i))
              dw3di(j)=ti*(w3(j,k,i+ip)-w3(j,k,i-im)) 
c
              dpdk(j)=tk*(q(j,k+kp,i,5)-q(j,k-km,i,5))
              dpdi(j)=ti*(q(j,k,i+ip,5)-q(j,k,i-im,5))
c
              drdk(j)=tk*(q(j,k+kp,i,1)-q(j,k-km,i,1))
              drdi(j)=ti*(q(j,k,i+ip,1)-q(j,k,i-im,1))
c
            enddo
            do j=1,jdim-1
c
              rvol = 1.0/vol(j,k,i)
              sj4=sj(j,k,i,4)
              sj4p=sj(j+1,k,i,4)
              xjp=0.5*(sj(j,k,i,1)*sj4+sj(j+1,k,i,1)*sj4p) 
              yjp=0.5*(sj(j,k,i,2)*sj4+sj(j+1,k,i,2)*sj4p)
              zjp=0.5*(sj(j,k,i,3)*sj4+sj(j+1,k,i,3)*sj4p)
c
              sk4=sk(j,k,i,4)
              sk4p=sk(j,k+1,i,4)
              xkp=0.5*(sk(j,k,i,1)*sk4+sk(j,k+1,i,1)*sk4p)
              ykp=0.5*(sk(j,k,i,2)*sk4+sk(j,k+1,i,2)*sk4p)
              zkp=0.5*(sk(j,k,i,3)*sk4+sk(j,k+1,i,3)*sk4p)
c
              si4=si(j,k,i,4)
              si4p=si(j,k,i+1,4)
              xip=0.5*(si(j,k,i,1)*si4+si(j,k,i+1,1)*si4p)
              yip=0.5*(si(j,k,i,2)*si4+si(j,k,i+1,2)*si4p)
              zip=0.5*(si(j,k,i,3)*si4+si(j,k,i+1,3)*si4p)
c
              dkdx(j)=(xjp*dkdj(j)+xkp*dkdk(j)+xip*dkdi(j))*rvol
              dkdy(j)=(yjp*dkdj(j)+ykp*dkdk(j)+yip*dkdi(j))*rvol
              dkdz(j)=(zjp*dkdj(j)+zkp*dkdk(j)+zip*dkdi(j))*rvol
c
              dzdx(j)=(xjp*dzdj(j)+xkp*dzdk(j)+xip*dzdi(j))*rvol
              dzdy(j)=(yjp*dzdj(j)+ykp*dzdk(j)+yip*dzdi(j))*rvol
              dzdz(j)=(zjp*dzdj(j)+zkp*dzdk(j)+zip*dzdi(j))*rvol
c
              dw1dx(j)=(xjp*dw1dj(j)+xkp*dw1dk(j)+xip*dw1di(j))*rvol
              dw1dy(j)=(yjp*dw1dj(j)+ykp*dw1dk(j)+yip*dw1di(j))*rvol
              dw1dz(j)=(zjp*dw1dj(j)+zkp*dw1dk(j)+zip*dw1di(j))*rvol
c
              dw2dx(j)=(xjp*dw2dj(j)+xkp*dw2dk(j)+xip*dw2di(j))*rvol
              dw2dy(j)=(yjp*dw2dj(j)+ykp*dw2dk(j)+yip*dw2di(j))*rvol
              dw2dz(j)=(zjp*dw2dj(j)+zkp*dw2dk(j)+zip*dw2di(j))*rvol
c
              dw3dx(j)=(xjp*dw3dj(j)+xkp*dw3dk(j)+xip*dw3di(j))*rvol
              dw3dy(j)=(yjp*dw3dj(j)+ykp*dw3dk(j)+yip*dw3di(j))*rvol
              dw3dz(j)=(zjp*dw3dj(j)+zkp*dw3dk(j)+zip*dw3di(j))*rvol
c
              dpdx(j)=(xjp*dpdj(j)+xkp*dpdk(j)+xip*dpdi(j))*rvol
              dpdy(j)=(yjp*dpdj(j)+ykp*dpdk(j)+yip*dpdi(j))*rvol
              dpdz(j)=(zjp*dpdj(j)+zkp*dpdk(j)+zip*dpdi(j))*rvol
c
              drdx(j)=(xjp*drdj(j)+xkp*drdk(j)+xip*drdi(j))*rvol
              drdy(j)=(yjp*drdj(j)+ykp*drdk(j)+yip*drdi(j))*rvol
              drdz(j)=(zjp*drdj(j)+zkp*drdk(j)+zip*drdi(j))*rvol
c
            enddo
            do j=1,jdim-1
c
              fnut=vist3d(j,k,i)
c
c ----- determine strain and vorticity values
c               
              s11 = ux(j,k,i,1)
              s22 = ux(j,k,i,5)
              s33 = ux(j,k,i,9)
              divv = s11 + s22 + s33   
c
              s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
              s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
              s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
              wis = w1(j,k,i)**2 + w2(j,k,i)**2 + w3(j,k,i)**2 
              omemag = sqrt(wis)
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
     +            + slim
c
c ----- compute tau_i,j 
c
              t11 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*fnut*(s11 - 0.33333*divv)
              t11=ccmaxcr(t11,0.)
              t22 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*fnut*(s22 - 0.33333*divv)
              t22=ccmaxcr(t22,0.)
              t33 = 2.*q(j,k,i,1)*turre(j,k,i,2)*re/3.
     +             -2.*fnut*(s33 - 0.33333*divv)
              t33=ccmaxcr(t33,0.)
              t12 =-2.*fnut*s12
              t13 =-2.*fnut*s13
              t23 =-2.*fnut*s23
c
c ---- calculate enstrophy production terms
c
              p3 = 2.0*fnut*alpha3*turre(j,k,i,1)*(xis-slim)
     +           / turre(j,k,i,2) 
     +           + 2.0*q(j,k,i,1)*re*turre(j,k,i,1)*divv/3.
c
              p8 = -2.0*beta8/turre(j,k,i,2)/(xis+0.5*wis)*
     +           ((w1(j,k,i)*t13+w2(j,k,i)*t23+w3(j,k,i)*t33)
     +           *(dkdx(j)*dzdy(j)-dkdy(j)*dzdx(j)) +
     +           (w1(j,k,i)*t11+w2(j,k,i)*t12+w3(j,k,i)*t13)
     +           *(dkdy(j)*dzdz(j)-dkdz(j)*dzdy(j)) +
     +           (w1(j,k,i)*t12+w2(j,k,i)*t22+w3(j,k,i)*t23)
     +           *(dkdz(j)*dzdx(j)-dkdx(j)*dzdz(j)))
c         
              wiwjtauij =   
     +            w1(j,k,i)**2*t11 + w2(j,k,i)**2*t22 + w3(j,k,i)**2*t33 
     +      + 2.*(w1(j,k,i)*w2(j,k,i)*t12 + w1(j,k,i)*w3(j,k,i)*t13
     +      + w2(j,k,i)*w3(j,k,i)*t23)
c
              p4p6 = (beta4*turre(j,k,i,1)/(omemag + wlim)
     +           + 2.0*beta6*fnut*omemag/(fnu(j,k,i)*re*re))
     +             *wiwjtauij/turre(j,k,i,2)
c
              p7 = beta7*re*q(j,k,i,1)*turre(j,k,i,1)*
     +           (w1(j,k,i)**2*s11 + w2(j,k,i)**2*s22 + w3(j,k,i)**2*s33 
     +      + 2.*(w1(j,k,i)*w2(j,k,i)*s12 + w1(j,k,i)*w3(j,k,i)*s13 
     +      + w2(j,k,i)*w3(j,k,i)*s23))
     +          /(wis+wlim)
c
              p5 = fnut/(sigr*re*re)*
     +           (2.*(dw1dx(j)**2+dw2dy(j)**2+dw3dz(j)**2)+
     +           (dw2dx(j)+dw1dy(j))**2+(dw3dx(j)+dw1dz(j))**2+
     +           (dw3dy(j)+dw2dz(j))**2)
c
c ---- calculate enstrophy destruction term (already divided by density)
c
              xnu = fnu(j,k,i)/q(j,k,i,1)
              rt = (turre(j,k,i,2)/xnu)**2/turre(j,k,i,1)
              ft = 1.90*(1. + 
     +        (beta5/1.90-1.)*cctanh(turre(j,k,i,1)/tur10(1)))/beta5
              ft= 1.0
              d5 = re*beta5*ft*turre(j,k,i,1)**1.5/(sqrt(rt)+delta)
c
c ---- calculate compressibility terms for k-equation
c
              ct1 = -fnut*(drdx(j)*dpdx(j) + drdy(j)*dpdy(j) 
     +            + drdz(j)*dpdz(j))/(re*ck*q(j,k,i,1)**3)  
              rtaurho = sqrt(turre(j,k,i,2)*(drdx(j)**2 + drdy(j)**2 +
     +                                       drdz(j)**2))/q(j,k,i,1)
              ct2 = -c1*turre(j,k,i,2)*rtaurho
c
c ---- calculate compressibility terms for enstrophy equation
c
              DpDt = q(j,k,i,2)*dpdx(j) + q(j,k,i,3)*dpdy(j) +
     +               q(j,k,i,4)*dpdz(j)

              termn = (q(j,k,i,1)*turre(j,k,i,2)*omemag)/
     +                (xnu*q(j,k,i,5)*sigp)   
              termd = 1. + sigrho*rtaurho*sqrt(rt/turre(j,k,i,1))/re
              ct3 = ccmaxcr(termn*DpDt/termd,0.)
              ct4 = -czeta1*fnut*omemag*turre(j,k,i,1)*rtaurho
     +             /turre(j,k,i,2)
              ct5 = -2.0*re*q(j,k,i,1)*turre(j,k,i,1)*divv
c                     
c ---- calculate source terms for k-equation
c
              if (ikoprod .eq. 2) then
              pk = -(t11*ux(j,k,i,1) + t22*ux(j,k,i,5)
     +              +t33*ux(j,k,i,9) + t12*(ux(j,k,i,2)+ux(j,k,i,4))
     +                               + t13*(ux(j,k,i,3)+ux(j,k,i,7))
     +                               + t23*(ux(j,k,i,6)+ux(j,k,i,8)))/
     +              (q(j,k,i,1)*re)
              else if (ikoprod .eq. 1) then
              pk=2.*vist3d(j,k,i)*xis/(q(j,k,i,1)*re)
              else
              pk=vist3d(j,k,i)*vor(j,k,i)**2/(q(j,k,i,1)*re)
              end if
c
              rtauk=re*fnu(j,k,i)*turre(j,k,i,1)/(turre(j,k,i,2)*
     +           q(j,k,i,1))
c
              dk = turre(j,k,i,2)*rtauk
c
              pklim = ccmin(pk,pklimterm*dk)
              pk = pklim
c
              srce(j,k,i,2) = pk + ct1 + ct2 - dk 
              rhside(j,k,i,2) = rhside(j,k,i,2)
     +       + srce(j,k,i,2)
c
              srce(j,k,i,1) = 
     +        (p3+p4p6 + p5 + p7 + p8 + ct3 + ct4 + ct5)
     +        /(q(j,k,i,1)*re) - d5
              rhside(j,k,i,1) = rhside(j,k,i,1) 
     +       + srce(j,k,i,1)
c
            enddo
          enddo
        enddo
      end if
c
c    For implicit viscous terms, ignore "catris_kw" - use standard
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
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
              cy(j,k)=cy(j,k) + uu*(app-apm)
              cy2(j,k)=cy2(j,k) + uu*(app-apm)
              dy(j,k)=dy(j,k) + uu*apm
              dy2(j,k)=dy2(j,k) + uu*apm
            enddo
c      Add part of source terms to diagonal LHS in this sweep (helps increase
c      diagonal dominance)
            do j=1,jdim-1
              if (ivmx .eq. 6) then
                cyadd= 2.*re*beta1*turre(j,k,i,1)
                cy2add=re*cmuc1*turre(j,k,i,1)
              else if (ivmx .eq. 7) then
                cmuc=  blend(j,k,i)*cmuc1 + (1.-blend(j,k,i))*cmuc2
                betax= blend(j,k,i)*beta1 + (1.-blend(j,k,i))*beta2
                cyadd= 2.*re*betax*turre(j,k,i,1) +
     +                 ccabs(damp1(j,k,i))/turre(j,k,i,1)
                cy2add=re*cmuc*turre(j,k,i,1)
              else if (ivmx .eq. 12) then
                cyadd= 2.*re*beta1*turre(j,k,i,1)
                cy2add=re*turre(j,k,i,1)
              else if (ivmx .eq. 8 .or. ivmx .eq. 14) then
                cyadd= 2.*re*beta1*turre(j,k,i,1)
                cy2add=damp1(j,k,i)*re*turre(j,k,i,1)
              else if (ivmx .eq.  9 .or. ivmx .eq.10 .or.
     +                 ivmx .eq. 11 .or. ivmx .eq.13) then
                cyadd= damp1(j,k,i)
                cy2add=2.*re*turre(j,k,i,1)/turre(j,k,i,2)
              else if (ivmx .eq.  15) then
c
                term1 = -0.5*(srce(j,k,i,1) 
     c                  - ccabs(srce(j,k,i,1)))
                term2 = -0.5*(srce(j,k,i,2) 
     c                  - ccabs(srce(j,k,i,2)))
                cyadd = term1/turre(j,k,i,1)   
                cy2add = term2/turre(j,k,i,2)   
c
              end if
              cy(j,k)=cy(j,k) + cyadd
              cy2(j,k)=cy2(j,k) + cy2add
            enddo
            do j=1,jdim-1
              by(j,k)=by(j,k)*timestp(j,k,i)
              by2(j,k)=by2(j,k)*timestp(j,k,i)
              cy(j,k)=cy(j,k)*timestp(j,k,i)+1.0*(1.+phi)
              cy2(j,k)=cy2(j,k)*timestp(j,k,i)+1.0*(1.+phi)
              dy(j,k)=dy(j,k)*timestp(j,k,i)
              dy2(j,k)=dy2(j,k)*timestp(j,k,i)
              fy(j,k)=rhside(j,k,i,1)*timestp(j,k,i)
              fy2(j,k)=rhside(j,k,i,2)*timestp(j,k,i)
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-
     +                 turre(j,k,i,1))+phi*zksav2(j,k,i,3)
                fy2(j,k)=fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-
     +                 turre(j,k,i,2))+phi*zksav2(j,k,i,4)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
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
              cy(j,k)=cy(j,k) + uu*(app-apm)
              cy2(j,k)=cy2(j,k) + uu*(app-apm)
              dy(j,k)=dy(j,k) + uu*apm
              dy2(j,k)=dy2(j,k) + uu*apm
            enddo
c      Add part of source terms to diagonal LHS in this sweep (helps increase
c      diagonal dominance)
            do j=1,jdim-1
              if (ivmx .eq. 6) then
                cyadd= 2.*re*beta1*turre(j,k,i,1)
                cy2add=re*cmuc1*turre(j,k,i,1)
              else if (ivmx .eq. 7) then
                cmuc=  blend(j,k,i)*cmuc1 + (1.-blend(j,k,i))*cmuc2
                betax= blend(j,k,i)*beta1 + (1.-blend(j,k,i))*beta2
                cyadd= 2.*re*betax*turre(j,k,i,1) +
     +                 ccabs(damp1(j,k,i))/turre(j,k,i,1)
                cy2add=re*cmuc*turre(j,k,i,1)
              else if (ivmx .eq. 12) then
                cyadd= 2.*re*beta1*turre(j,k,i,1)
                cy2add=re*turre(j,k,i,1)
              else if (ivmx .eq. 8 .or. ivmx .eq. 14) then
                cyadd= 2.*re*beta1*turre(j,k,i,1)
                cy2add=damp1(j,k,i)*re*turre(j,k,i,1)
              else if (ivmx .eq.  9 .or. ivmx .eq.10 .or.
     +                 ivmx .eq. 11 .or. ivmx .eq.13) then
                cyadd= damp1(j,k,i)
                cy2add=2.*re*turre(j,k,i,1)/turre(j,k,i,2)
              else if (ivmx .eq.  15) then
c
                term1 = -0.5*(srce(j,k,i,1) 
     c                  - ccabs(srce(j,k,i,1)))
                term2 = -0.5*(srce(j,k,i,2) 
     c                  - ccabs(srce(j,k,i,2)))
                cyadd = term1/turre(j,k,i,1)   
                cy2add = term2/turre(j,k,i,2)   
              end if
              cy(j,k)=cy(j,k) + cyadd
              cy2(j,k)=cy2(j,k) + cy2add
            enddo
            do j=1,jdim-1
              by(j,k)=by(j,k)*timestp(j,k,i)
              by2(j,k)=by2(j,k)*timestp(j,k,i)
              cy(j,k)=cy(j,k)*timestp(j,k,i)+1.0*(1.+phi)
              cy2(j,k)=cy2(j,k)*timestp(j,k,i)+1.0*(1.+phi)
              dy(j,k)=dy(j,k)*timestp(j,k,i)
              dy2(j,k)=dy2(j,k)*timestp(j,k,i)
              fy(j,k)=rhside(j,k,i,1)*timestp(j,k,i)
              fy2(j,k)=rhside(j,k,i,2)*timestp(j,k,i)
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-
     +                 turre(j,k,i,1))+phi*zksav2(j,k,i,3)
                fy2(j,k)=fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-
     +                 turre(j,k,i,2))+phi*zksav2(j,k,i,4)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
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
              cy(j,k)=cy(j,k) + uu*(app-apm)
              cy2(j,k)=cy2(j,k) + uu*(app-apm)
              dy(j,k)=dy(j,k) + uu*apm
              dy2(j,k)=dy2(j,k) + uu*apm
            enddo
c      Add part of source terms to diagonal LHS in this sweep (helps increase
c      diagonal dominance)
            do j=1,jdim-1
              if (ivmx .eq. 6) then
                cyadd= 2.*re*beta1*turre(j,k,i,1)
                cy2add=re*cmuc1*turre(j,k,i,1)
              else if (ivmx .eq. 7) then
                cmuc=  blend(j,k,i)*cmuc1 + (1.-blend(j,k,i))*cmuc2
                betax= blend(j,k,i)*beta1 + (1.-blend(j,k,i))*beta2
                cyadd= 2.*re*betax*turre(j,k,i,1) +
     +                 ccabs(damp1(j,k,i))/turre(j,k,i,1)
                cy2add=re*cmuc*turre(j,k,i,1)
              else if (ivmx .eq. 12) then
                cyadd= 2.*re*beta1*turre(j,k,i,1)
                cy2add=re*turre(j,k,i,1)
              else if (ivmx .eq. 8 .or. ivmx .eq. 14) then
                cyadd= 2.*re*beta1*turre(j,k,i,1)
                cy2add=damp1(j,k,i)*re*turre(j,k,i,1)
              else if (ivmx .eq.  9 .or. ivmx .eq.10 .or.
     +                 ivmx .eq. 11 .or. ivmx .eq.13) then
                cyadd= damp1(j,k,i)
                cy2add=2.*re*turre(j,k,i,1)/turre(j,k,i,2)
              else if (ivmx .eq.  15) then
c
                term1 = -0.5*(srce(j,k,i,1) 
     c                  - ccabs(srce(j,k,i,1)))
                term2 = -0.5*(srce(j,k,i,2) 
     c                  - ccabs(srce(j,k,i,2)))
                cyadd = term1/turre(j,k,i,1)   
                cy2add = term2/turre(j,k,i,2)   
              end if
              cy(j,k)=cy(j,k) + cyadd
              cy2(j,k)=cy2(j,k) + cy2add
            enddo
            do j=1,jdim-1
              by(j,k)=by(j,k)*timestp(j,k,i)
              by2(j,k)=by2(j,k)*timestp(j,k,i)
              cy(j,k)=cy(j,k)*timestp(j,k,i)+1.0*(1.+phi)
              cy2(j,k)=cy2(j,k)*timestp(j,k,i)+1.0*(1.+phi)
              dy(j,k)=dy(j,k)*timestp(j,k,i)
              dy2(j,k)=dy2(j,k)*timestp(j,k,i)
              fy(j,k)=rhside(j,k,i,1)*timestp(j,k,i)
              fy2(j,k)=rhside(j,k,i,2)*timestp(j,k,i)
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+(1.+phi)*(zksav2(j,k,i,1)-
     +                 turre(j,k,i,1))+phi*zksav2(j,k,i,3)
                fy2(j,k)=fy2(j,k)+(1.+phi)*(zksav2(j,k,i,2)-
     +                 turre(j,k,i,2))+phi*zksav2(j,k,i,4)
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
              enddo
            enddo
          end if
          call triv(jdim-1,kdim-1,1,jdim-1,1,kdim-1,worky,by,cy,dy,fy)
          call triv(jdim-1,kdim-1,1,jdim-1,1,kdim-1,worky,by2,cy2,
     +     dy2,fy2)
          do k=1,kdim-1
            do j=1,jdim-1
              rhside(j,k,i,1)=fy(j,k)
              rhside(j,k,i,2)=fy2(j,k)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
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
              cx(k,j)=cx(k,j) + uu*(app-apm)
              cx2(k,j)=cx2(k,j) + uu*(app-apm)
              dx(k,j)=dx(k,j) + uu*apm
              dx2(k,j)=dx2(k,j) + uu*apm
            enddo
            do k=1,kdim-1
              bx(k,j)=bx(k,j)*timestp(j,k,i)
              bx2(k,j)=bx2(k,j)*timestp(j,k,i)
              cx(k,j)=cx(k,j)*timestp(j,k,i)+1.0*(1.+phi)
              cx2(k,j)=cx2(k,j)*timestp(j,k,i)+1.0*(1.+phi)
              dx(k,j)=dx(k,j)*timestp(j,k,i)
              dx2(k,j)=dx2(k,j)*timestp(j,k,i)
              fx(k,j)=rhside(j,k,i,1)*(1.+phi)
              fx2(k,j)=rhside(j,k,i,2)*(1.+phi)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
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
              cx(k,j)=cx(k,j) + uu*(app-apm)
              cx2(k,j)=cx2(k,j) + uu*(app-apm)
              dx(k,j)=dx(k,j) + uu*apm
              dx2(k,j)=dx2(k,j) + uu*apm
            enddo
            do k=1,kdim-1
              bx(k,j)=bx(k,j)*timestp(j,k,i)
              bx2(k,j)=bx2(k,j)*timestp(j,k,i)
              cx(k,j)=cx(k,j)*timestp(j,k,i)+1.0*(1.+phi)
              cx2(k,j)=cx2(k,j)*timestp(j,k,i)+1.0*(1.+phi)
              dx(k,j)=dx(k,j)*timestp(j,k,i)
              dx2(k,j)=dx2(k,j)*timestp(j,k,i)
              fx(k,j)=rhside(j,k,i,1)*(1.+phi)
              fx2(k,j)=rhside(j,k,i,2)*(1.+phi)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
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
              cx(k,j)=cx(k,j) + uu*(app-apm)
              cx2(k,j)=cx2(k,j) + uu*(app-apm)
              dx(k,j)=dx(k,j) + uu*apm
              dx2(k,j)=dx2(k,j) + uu*apm
            enddo
            do k=1,kdim-1
              bx(k,j)=bx(k,j)*timestp(j,k,i)
              bx2(k,j)=bx2(k,j)*timestp(j,k,i)
              cx(k,j)=cx(k,j)*timestp(j,k,i)+1.0*(1.+phi)
              cx2(k,j)=cx2(k,j)*timestp(j,k,i)+1.0*(1.+phi)
              dx(k,j)=dx(k,j)*timestp(j,k,i)
              dx2(k,j)=dx2(k,j)*timestp(j,k,i)
              fx(k,j)=rhside(j,k,i,1)*(1.+phi)
              fx2(k,j)=rhside(j,k,i,2)*(1.+phi)
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
              enddo
            enddo
          end if
          call triv(kdim-1,jdim-1,1,kdim-1,1,jdim-1,workx,bx,cx,dx,fx)
          call triv(kdim-1,jdim-1,1,kdim-1,1,jdim-1,workx,bx2,
     +     cx2,dx2,fx2)
          do j=1,jdim-1
            do k=1,kdim-1
              rhside(j,k,i,1)=fx(k,j)
              rhside(j,k,i,2)=fx2(k,j)
            enddo
          enddo
        enddo
c
c    Implicit F_zeta_zeta viscous terms.  Do over all j's
        if(i2d .ne. 1) then
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
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
                cz(k,i)=cz(k,i) + uu*(app-apm)
                cz2(k,i)=cz2(k,i) + uu*(app-apm)
                dz(k,i)=dz(k,i) + uu*apm
                dz2(k,i)=dz2(k,i) + uu*apm
              enddo
              do k=1,kdim-1
                bz(k,i)=bz(k,i)*timestp(j,k,i)
                bz2(k,i)=bz2(k,i)*timestp(j,k,i)
                cz(k,i)=cz(k,i)*timestp(j,k,i)+1.0*(1.+phi)
                cz2(k,i)=cz2(k,i)*timestp(j,k,i)+1.0*(1.+phi)
                dz(k,i)=dz(k,i)*timestp(j,k,i)
                dz2(k,i)=dz2(k,i)*timestp(j,k,i)
                fz(k,i)=rhside(j,k,i,1)*(1.+phi)
                fz2(k,i)=rhside(j,k,i,2)*(1.+phi)
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
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
                cz(k,i)=cz(k,i) + uu*(app-apm)
                cz2(k,i)=cz2(k,i) + uu*(app-apm)
                dz(k,i)=dz(k,i) + uu*apm
                dz2(k,i)=dz2(k,i) + uu*apm
              enddo
              do k=1,kdim-1
                bz(k,i)=bz(k,i)*timestp(j,k,i)
                bz2(k,i)=bz2(k,i)*timestp(j,k,i)
                cz(k,i)=cz(k,i)*timestp(j,k,i)+1.0*(1.+phi)
                cz2(k,i)=cz2(k,i)*timestp(j,k,i)+1.0*(1.+phi)
                dz(k,i)=dz(k,i)*timestp(j,k,i)
                dz2(k,i)=dz2(k,i)*timestp(j,k,i)
                fz(k,i)=rhside(j,k,i,1)*(1.+phi)
                fz2(k,i)=rhside(j,k,i,2)*(1.+phi)
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
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
                cz(k,i)=cz(k,i) + uu*(app-apm)
                cz2(k,i)=cz2(k,i) + uu*(app-apm)
                dz(k,i)=dz(k,i) + uu*apm
                dz2(k,i)=dz2(k,i) + uu*apm
              enddo
              do k=1,kdim-1
                bz(k,i)=bz(k,i)*timestp(j,k,i)
                bz2(k,i)=bz2(k,i)*timestp(j,k,i)
                cz(k,i)=cz(k,i)*timestp(j,k,i)+1.0*(1.+phi)
                cz2(k,i)=cz2(k,i)*timestp(j,k,i)+1.0*(1.+phi)
                dz(k,i)=dz(k,i)*timestp(j,k,i)
                dz2(k,i)=dz2(k,i)*timestp(j,k,i)
                fz(k,i)=rhside(j,k,i,1)*(1.+phi)
                fz2(k,i)=rhside(j,k,i,2)*(1.+phi)
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
                enddo
              enddo
            end if
            call triv(kdim-1,idim-1,1,kdim-1,1,idim-1,workz,bz,
     +                cz,dz,fz)
            call triv(kdim-1,idim-1,1,kdim-1,1,idim-1,workz,bz2,
     +                cz2,dz2,fz2)
            do i=1,idim-1
              do k=1,kdim-1
                rhside(j,k,i,1)=fz(k,i)
                rhside(j,k,i,2)=fz2(k,i)
              enddo
            enddo
          enddo
        end if
c
c    Update TURRE
        sumno=0.
        sumnk=0.
        negno=0
        negnk=0
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
              sumnk=sumnk+rhside(j,k,i,2)**2
              if((real(turre(j,k,i,2)+rhside(j,k,i,2))) .le. 
     +          real(tur2cutlev)) then
                negnk=negnk+1
                if (real(tur2cut).gt.0.) turre(j,k,i,2)=tur2cut
              else
                turre(j,k,i,2)=turre(j,k,i,2)+rhside(j,k,i,2)
              end if
            enddo
          enddo
        enddo
        sumno=sqrt(sumno)/float((kdim-1)*(jdim-1)*(idim-1))
        sumnk=sqrt(sumnk)/float((kdim-1)*(jdim-1)*(idim-1))
c
        if(iwrite .eq. 1) then
c         write(15,'('' icyc, it, logreso, negno, logresk, negnk'',
c    +     '', block'',/,2i5,e15.5,i5,e15.5,2i5)') icyc,not,
c    +     real(cclog10(sumno)),negno,real(cclog10(sumnk)),negnk,nbl
        end if
c
 500  continue
      sumn1 = sumno
      sumn2 = sumnk
      negn1 = negno
      negn2 = negnk
c
c Update VIST3D and save omega and k values
      if (ivmx .eq. 6) then
        do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            if (i_wilcox06 .eq. 1) then
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
              xxx=sqrt(2.*xis)
              denom=ccmax(turre(j,k,i,1),xxx*0.875/(re*sqrt(cmuc2)))
              vist3d(j,k,i)=q(j,k,i,1)*turre(j,k,i,2)/denom
            else
              vist3d(j,k,i)=q(j,k,i,1)*turre(j,k,i,2)/turre(j,k,i,1)
            end if
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),edvislim)
            zksav(j,k,i,1)=turre(j,k,i,1)
            zksav(j,k,i,2)=turre(j,k,i,2)
          enddo
        enddo
        enddo
      else if (ivmx .eq. 7) then
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
          enddo
        enddo
        enddo
      else if (ivmx .eq. 12) then
        do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
c          Compute variable cmu coefficient:
            alfa1 = (1.333333-c2)*gg/2.
            alfa2 = (2.-c3)**2*(gg*gg/4.)
            alfa3 = (2.-c4)**2*(gg*gg/4.)
c          Determine Sij and Wij values:
            s11 = ux(j,k,i,1)
            s22 = ux(j,k,i,5)
            s33 = ux(j,k,i,9)
            s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
            s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
            s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
            w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
            w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
            if (ieasmcc2d .eq. 1) then
              w13=w13+(2./(c4-2.))*vx(j,k,i,1)
            end if
            w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
            xis = s11*s11 + s22*s22 + s33*s33 +
     +            2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
            wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
            eta   = alfa2*xis/(turre(j,k,i,1)*re)**2
            squig = alfa3*wis/(turre(j,k,i,1)*re)**2
c          Limit eta and squig:
            eta=ccmincr(eta,10.0)
            squig=ccmincr(squig,10.0)
c  (Note:  in Abid journal paper, with LRR, did not have "0.2" terms:)
            cmu = alfa1*((3.*(1.+eta)+0.2*eta*eta*eta+
     +                                0.2*squig*squig*squig)/
     +            (3.+eta+6.*eta*squig+6.*squig+eta*eta*eta+
     +                                          squig*squig*squig))
c          Limit cmu:
            cmu = ccmincr(cmu,.187)
            cmu = ccmaxcr(cmu,.005)
            vist3d(j,k,i)=cmu*q(j,k,i,1)*turre(j,k,i,2)/
     +       turre(j,k,i,1)
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),edvislim)
            zksav(j,k,i,1)=turre(j,k,i,1)
            zksav(j,k,i,2)=turre(j,k,i,2)
          enddo
        enddo
        enddo
      else if (ivmx .eq.10) then
        do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            ret=q(j,k,i,1)*turre(j,k,i,2)**2/(fnu(j,k,i)*
     +          turre(j,k,i,1))
            rek=q(j,k,i,1)*sqrt(turre(j,k,i,2))*re*ccabs(smin(j,k,i))/
     +          fnu(j,k,i)
            fmu=(1.+4.*(ret**(-.75)))*cctanh(.008*rek)
            fmu=ccmincr(fmu,1.0)
            vist3d(j,k,i)=cmuc1*fmu*q(j,k,i,1)*turre(j,k,i,2)**2/
     +       turre(j,k,i,1)
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),edvislim)
            zksav(j,k,i,1)=turre(j,k,i,1)
            zksav(j,k,i,2)=turre(j,k,i,2)
          enddo
        enddo
        enddo
      else if (ivmx .eq.11) then
        do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
c          Compute variable cmu coefficient:
            alfa1 = (1.333333-c2)*gg/2.
            alfa2 = (2.-c3)**2*(gg*gg/4.)
            alfa3 = (2.-c4)**2*(gg*gg/4.)
c            Determine Sij and Wij values:
            s11 = ux(j,k,i,1)
            s22 = ux(j,k,i,5)
            s33 = ux(j,k,i,9)
            s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
            s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
            s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
            w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
            w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
            if (ieasmcc2d .eq. 1) then
              w13=w13+(2./(c4-2.))*vx(j,k,i,1)
            end if
            w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
            xis = s11*s11 + s22*s22 + s33*s33 +
     +            2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
            wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
            eta   = alfa2*xis*turre(j,k,i,2)**2/
     +              (turre(j,k,i,1)*re)**2
            squig = alfa3*wis*turre(j,k,i,2)**2/
     +              (turre(j,k,i,1)*re)**2
c          Limit eta and squig:
            eta=ccmincr(eta,10.0)
            squig=ccmincr(squig,10.0)
c  (Note:  in Abid journal paper, with LRR, did not have "0.2" terms:)
            cmu = alfa1*((3.*(1.+eta)+0.2*eta*eta*eta+
     +                                0.2*squig*squig*squig)/
     +            (3.+eta+6.*eta*squig+6.*squig+eta*eta*eta+
     +                                          squig*squig*squig))
c          Limit cmu:
            cmu = ccmincr(cmu,.187)
            cmu = ccmaxcr(cmu,.005)
            vist3d(j,k,i)=cmu*q(j,k,i,1)*turre(j,k,i,2)**2/
     +       turre(j,k,i,1)
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),edvislim)
            zksav(j,k,i,1)=turre(j,k,i,1)
            zksav(j,k,i,2)=turre(j,k,i,2)
          enddo
        enddo
        enddo
      else if (ivmx .eq. 9 .or. ivmx .eq. 13) then
        pi = acos(-1.)
        tpi3=2.*pi/3.
        al10 = 0.25*c1 -1.
        al1=3.8
c   2-line modification to improve center of WAKE:
        al10=al10+1.8864
        al1=al1-2.
        al2=0.5*c2-2./3.
        al3=0.5*c3-1.
        al4=0.5*c4-1.
        do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
c          Compute variable cmu coefficient:
            fmu=1.0
c            Determine Sij and Wij values:
            s11 = ux(j,k,i,1)
            s22 = ux(j,k,i,5)
            s33 = ux(j,k,i,9)
            s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
            s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
            s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
            w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
            w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
            if (ieasmcc2d .eq. 1) then
              w13=w13+(2./(c4-2.))*vx(j,k,i,1)
            end if
            w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
            xis = s11*s11 + s22*s22 + s33*s33 +
     +            2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
            wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
            eta1  = xis*turre(j,k,i,2)**2/
     +              (turre(j,k,i,1)*re)**2
            eta2  = wis*turre(j,k,i,2)**2/
     +              (turre(j,k,i,1)*re)**2
c          Limit  eta1 and eta2
            eta1=ccmincr(eta1,1200.0)
            eta2=ccmincr(eta2,1200.0)
            if(real(eta1).le.1.e-6) then
            cmuj=al10*al2/(al10*al10+2.*eta2*al4*al4)
            else
            eta11=eta1*al1
            ap=-2.*al10/eta11
            ar=-al10*al2/eta11**2
            aq=al10*al10+eta11*al2-0.666667*eta1*al3*al3
     +         +2.*eta2*al4*al4
            aq=aq/(eta11*eta11)
            aa=(aq-ap*ap/3.)
            ab=(2.*ap*ap*ap-9.*ap*aq+27.*ar)/27.
            ad=(ab*ab/4.)+(aa*aa*aa)/27.
            if(real(ad).gt.0.) then
            raat=-0.5*ab+sqrt(ad)
            rbbt=-0.5*ab-sqrt(ad)
            raa=(ccabs(raat))**(1./3.)
            raa=ccsign(raa,raat)
            rbb=(ccabs(rbbt))**(1./3.)
            rbb=ccsign(rbb,rbbt)
            cmuj=-ap/3.+raa+rbb
c   choose min of real root AND real part of imaj roots:
            cmub=-ap/3.-.5*raa-.5*rbb
            cmuj=ccmin(cmuj,cmub)
            else
            coss=-ab/2./sqrt(-aa*aa*aa/27.)
            theta=ccacos(coss)
c   choose min of real roots
            cmuj=-ap/3.+2.*sqrt(-aa/3.)*cos(theta/3.)
            cmub=-ap/3.+2.*sqrt(-aa/3.)*cos(tpi3+theta/3.)
            cmuc=-ap/3.+2.*sqrt(-aa/3.)*cos(2.*tpi3+theta/3.)
            cmuj=ccmin(cmuj,cmub)
            cmuj=ccmin(cmuj,cmuc)
            end if
            end if
            cmuj = ccmincr(cmuj,-.0005)
c Uncomment next line for constant cmu:
c           cmuj=-cmuc1
            cmu = -fmu*cmuj
            vist3d(j,k,i)=cmu*q(j,k,i,1)*turre(j,k,i,2)**2/
     +       turre(j,k,i,1)
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),edvislim)
            zksav(j,k,i,1)=turre(j,k,i,1)
            zksav(j,k,i,2)=turre(j,k,i,2)
            cmuv(j,k,i) = cmuj
          enddo
        enddo
        enddo
      else if (ivmx .eq. 8 .or. ivmx .eq. 14) then
        pi = acos(-1.)
        tpi3=2.*pi/3.
        al10 = 0.25*c1 -1.
        al1=3.8
c   2-line modification to improve center of WAKE:
        al10=al10+1.8864
        al1=al1-2.
        al2=0.5*c2-2./3.
        al3=0.5*c3-1.
        al4=0.5*c4-1.
        do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
c            Determine Sij and Wij values:
            s11 = ux(j,k,i,1)
            s22 = ux(j,k,i,5)
            s33 = ux(j,k,i,9)
            s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
            s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
            s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
            w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
            w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
            if (ieasmcc2d .eq. 1) then
              w13=w13+(2./(c4-2.))*vx(j,k,i,1)
            end if
            w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
            xis = s11*s11 + s22*s22 + s33*s33 +
     +            2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
            wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
            eta1  = xis/(turre(j,k,i,1)*re)**2
            eta2  = wis/(turre(j,k,i,1)*re)**2
            if(real(eta1).le.1.e-6) then
            cmuj=al10*al2/(al10*al10+2.*eta2*al4*al4)
            else
            eta11=eta1*al1
            ap=-2.*al10/eta11
            ar=-al10*al2/eta11**2
            aq=al10*al10+eta11*al2-0.666667*eta1*al3*al3
     +         +2.*eta2*al4*al4
            aq=aq/(eta11*eta11)
            aa=(aq-ap*ap/3.)
            ab=(2.*ap*ap*ap-9.*ap*aq+27.*ar)/27.
            ad=(ab*ab/4.)+(aa*aa*aa)/27.
            if(real(ad).gt.0.) then
            raat=-0.5*ab+sqrt(ad)
            rbbt=-0.5*ab-sqrt(ad)
            raa=(ccabs(raat))**(1./3.)
            raa=ccsign(raa,raat)
            rbb=(ccabs(rbbt))**(1./3.)
            rbb=ccsign(rbb,rbbt)
            cmuj=-ap/3.+raa+rbb
c   choose min of real root AND real part of imaj roots:
            cmub=-ap/3.-.5*raa-.5*rbb
            cmuj=ccmin(cmuj,cmub)
            else
            coss=-ab/2./sqrt(-aa*aa*aa/27.)
            theta=ccacos(coss)
c   choose min of real roots
            cmuj=-ap/3.+2.*sqrt(-aa/3.)*cos(theta/3.)
            cmub=-ap/3.+2.*sqrt(-aa/3.)*cos(tpi3+theta/3.)
            cmuc=-ap/3.+2.*sqrt(-aa/3.)*cos(2.*tpi3+theta/3.)
            cmuj=ccmin(cmuj,cmub)
            cmuj=ccmin(cmuj,cmuc)
            end if
            end if
            cmuj = ccmincr(cmuj,-.0005)
c Uncomment next line for constant cmu:
c           cmuj=-cmuc1
            cmu = -cmuj
            vist3d(j,k,i)=cmu*q(j,k,i,1)*turre(j,k,i,2)/
     +       turre(j,k,i,1)
            vist3d(j,k,i)=ccmin(vist3d(j,k,i),edvislim)
            zksav(j,k,i,1)=turre(j,k,i,1)
            zksav(j,k,i,2)=turre(j,k,i,2)
            cmuv(j,k,i) = cmuj
          enddo
        enddo
        enddo
c
      else if (ivmx .eq. 15) then
c
        do i=1,idim-1
        do k=1,kdim-1
        do j=1,jdim-1
          vist3d(j,k,i)=q(j,k,i,1)**2*cmu*turre(j,k,i,2)**2/
     +                   (fnu(j,k,i)*turre(j,k,i,1))
          vist3d(j,k,i)=ccmin(vist3d(j,k,i),edvislim)
          zksav(j,k,i,1)=turre(j,k,i,1)
          zksav(j,k,i,2)=turre(j,k,i,2)
        end do
        end do
        end do
c
      end if
c   for wall functions, force k & omega (or epsilon) in 1st cell off solid wall
c     wall function in k-direction
      if(iwf(3) .eq. 1) then
      do 9181 kk=1,2
      if(kk .eq. 1) then
        k=1
        kstop=nbck0(nbl)
      else
        k=kdim-1
        kstop=nbckdim(nbl)
      end if
      do 9171 kset=1,kstop
        if(abs(kbcinfo(nbl,kset,1,kk)).eq.2004 .or.
     .     abs(kbcinfo(nbl,kset,1,kk)).eq.2024 .or.
     .     abs(kbcinfo(nbl,kset,1,kk)).eq.2016) then
          ibeg=kbcinfo(nbl,kset,2,kk)
          iend=kbcinfo(nbl,kset,3,kk)-1
          jbeg=kbcinfo(nbl,kset,4,kk)
          jend=kbcinfo(nbl,kset,5,kk)-1
          do i=ibeg,iend
            do j=jbeg,jend
              if ((i.ge.ilamlo .and. i.lt.ilamhi .and.
     .             j.ge.jlamlo .and. j.lt.jlamhi .and.
     .             k.ge.klamlo .and. k.lt.klamhi) .or.
     .             real(smin(j,k,i)) .lt. 0.) then
                continue
              else
                uu  = sqrt((q(j,k,i,2)-qk0(j,i,2,kk+kk-1))**2 +
     +                     (q(j,k,i,3)-qk0(j,i,3,kk+kk-1))**2 +
     +                     (q(j,k,i,4)-qk0(j,i,4,kk+kk-1))**2 )
                dudy=uu/ccabs(smin(j,k,i))
                tauw=(fnu(j,k,i)+vk0(j,i,1,kk+kk-1))*dudy
                utau=sqrt(tauw/q(j,k,i,1)/re)
                omegatemp=utau/(sqrt(cmuc1)*vk*ccabs(smin(j,k,i))*re)
                omegatemp=ccmax(omegatemp,tur10(1))
                zktemp=utau*utau/sqrt(cmuc1)
                vist3d(j,k,i)=q(j,k,i,1)*zktemp/omegatemp
                if(ivmx .eq. 6 .or. ivmx .eq. 7) then
                  zksav(j,k,i,1)=omegatemp
                  zksav(j,k,i,2)=zktemp
                else if(ivmx .eq. 8 .or. ivmx .eq. 12 .or.
     +                  ivmx .eq. 14) then
                  zksav(j,k,i,1)=omegatemp*cmuc1
                  zksav(j,k,i,2)=zktemp
                else if(ivmx .eq. 15) then
                  zksav(j,k,i,1)=omegatemp*cmuc1*zktemp*q(j,k,i,1)/
     +              fnu(j,k,i)
                  zksav(j,k,i,2)=zktemp
                else
                  zksav(j,k,i,1)=omegatemp*cmuc1*zktemp
                  zksav(j,k,i,2)=zktemp
                end if
              end if
            enddo
          enddo
        end if
 9171 continue
 9181 continue
      end if
c     wall function in j-direction
      if(iwf(2) .eq. 1) then
      do 9182 jj=1,2
      if(jj .eq. 1) then
        j=1
        jstop=nbcj0(nbl)
      else
        j=jdim-1
        jstop=nbcjdim(nbl)
      end if
      do 9172 jset=1,jstop
        if(abs(jbcinfo(nbl,jset,1,jj)).eq.2004 .or.
     .     abs(jbcinfo(nbl,jset,1,jj)).eq.2024 .or.
     .     abs(jbcinfo(nbl,jset,1,jj)).eq.2016) then
          ibeg=jbcinfo(nbl,jset,2,jj)
          iend=jbcinfo(nbl,jset,3,jj)-1
          kbeg=jbcinfo(nbl,jset,4,jj)
          kend=jbcinfo(nbl,jset,5,jj)-1
          do i=ibeg,iend
            do k=kbeg,kend
              if ((i.ge.ilamlo .and. i.lt.ilamhi .and.
     .             j.ge.jlamlo .and. j.lt.jlamhi .and.
     .             k.ge.klamlo .and. k.lt.klamhi) .or.
     .             real(smin(j,k,i)) .lt. 0.) then
                continue
              else
                uu  = sqrt((q(j,k,i,2)-qj0(k,i,2,jj+jj-1))**2 +
     +                     (q(j,k,i,3)-qj0(k,i,3,jj+jj-1))**2 +
     +                     (q(j,k,i,4)-qj0(k,i,4,jj+jj-1))**2 )
                dudy=uu/ccabs(smin(j,k,i))
                tauw=(fnu(j,k,i)+vj0(k,i,1,jj+jj-1))*dudy
                utau=sqrt(tauw/q(j,k,i,1)/re)
                omegatemp=utau/(sqrt(cmuc1)*vk*ccabs(smin(j,k,i))*re)
                omegatemp=ccmax(omegatemp,tur10(1))
                zktemp=utau*utau/sqrt(cmuc1)
                vist3d(j,k,i)=q(j,k,i,1)*zktemp/omegatemp
                if(ivmx .eq. 6 .or. ivmx .eq. 7) then
                  zksav(j,k,i,1)=omegatemp
                  zksav(j,k,i,2)=zktemp
                else if(ivmx .eq. 8 .or. ivmx .eq. 12 .or.
     +                  ivmx .eq. 14) then
                  zksav(j,k,i,1)=omegatemp*cmuc1
                  zksav(j,k,i,2)=zktemp
                else if(ivmx .eq. 15) then
                  zksav(j,k,i,1)=omegatemp*cmuc1*zktemp*q(j,k,i,1)/
     +              fnu(j,k,i)
                  zksav(j,k,i,2)=zktemp
                else
                  zksav(j,k,i,1)=omegatemp*cmuc1*zktemp
                  zksav(j,k,i,2)=zktemp
                end if
              end if
            enddo
          enddo
        end if
 9172 continue
 9182 continue
      end if
c     wall function in i-direction
      if(i2d .ne. 1 .and. iwf(1) .eq. 1) then
      do 9183 ii=1,2
      if(ii .eq. 1) then
        i=1
        istop=nbci0(nbl)
      else
        i=idim-1
        istop=nbcidim(nbl)
      end if
      do 9173 iset=1,istop
        if(abs(ibcinfo(nbl,iset,1,ii)).eq.2004 .or.
     .     abs(ibcinfo(nbl,iset,1,ii)).eq.2024 .or.
     .     abs(ibcinfo(nbl,iset,1,ii)).eq.2016) then
          jbeg=ibcinfo(nbl,iset,2,ii)
          jend=ibcinfo(nbl,iset,3,ii)-1
          kbeg=ibcinfo(nbl,iset,4,ii)
          kend=ibcinfo(nbl,iset,5,ii)-1
          do j=jbeg,jend
            do k=kbeg,kend
              if ((i.ge.ilamlo .and. i.lt.ilamhi .and.
     .             j.ge.jlamlo .and. j.lt.jlamhi .and.
     .             k.ge.klamlo .and. k.lt.klamhi) .or.
     .             real(smin(j,k,i)) .lt. 0.) then
                continue
              else
                uu  = sqrt((q(j,k,i,2)-qi0(k,i,2,ii+ii-1))**2 +
     +                     (q(j,k,i,3)-qi0(k,i,3,ii+ii-1))**2 +
     +                     (q(j,k,i,4)-qi0(k,i,4,ii+ii-1))**2 )
                dudy=uu/ccabs(smin(j,k,i))
                tauw=(fnu(j,k,i)+vi0(k,i,1,ii+ii-1))*dudy
                utau=sqrt(tauw/q(j,k,i,1)/re)
                omegatemp=utau/(sqrt(cmuc1)*vk*ccabs(smin(j,k,i))*re)
                omegatemp=ccmax(omegatemp,tur10(1))
                zktemp=utau*utau/sqrt(cmuc1)
                vist3d(j,k,i)=q(j,k,i,1)*zktemp/omegatemp
                if(ivmx .eq. 6 .or. ivmx .eq. 7) then
                  zksav(j,k,i,1)=omegatemp
                  zksav(j,k,i,2)=zktemp
                else if(ivmx .eq. 8 .or. ivmx .eq. 12 .or.
     +                  ivmx .eq. 14) then
                  zksav(j,k,i,1)=omegatemp*cmuc1
                  zksav(j,k,i,2)=zktemp
                else if(ivmx .eq. 15) then
                  zksav(j,k,i,1)=omegatemp*cmuc1*zktemp*q(j,k,i,1)/
     +              fnu(j,k,i)
                  zksav(j,k,i,2)=zktemp
                else
                  zksav(j,k,i,1)=omegatemp*cmuc1*zktemp
                  zksav(j,k,i,2)=zktemp
                end if
              end if
            enddo
          enddo
        end if
 9173 continue
 9183 continue
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
      write(92,'(''  utau='',e18.8)') utau
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
