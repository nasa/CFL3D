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
      subroutine spalart(jdim,kdim,idim,q,sj,sk,si,vol,dtj,
     + x,y,z,vist3d,vor,smin,tursav,turre,damp1,timestp,fnu,bx,cx,dx,
     + fx,workx,by,cy,dy,fy,worky,bz,cz,dz,fz,workz,ntime,
     + tj0,tk0,ti0,nbl,qj0,qk0,qi0,blank,iover,sumn1,sumn2,negn1,negn2,
     + tursav2,volj0,volk0,voli0,nou,bou,nbuf,ibufdim,iex,iex2,iex3,
     + bcj,bck,bci,ux,vx,nummem,fd,isarcnum)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute turbulent viscosity distributions using
c     1-equation Spalart turbulence model (ivisc=5)
c     (reprogrammed by Rumsey for CFL3D 1/92)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /des/ cdes,ides,cddes
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /lam/ ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi,
     .        i_lam_forcezero
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /twod/ i2d
      common /zero/ iexp
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /curvat/ isarc2d,sarccr3,ieasmcc2d,isstrc,sstrc_crc,
     .        isar,crot,isarc3d
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /mms/ iexact_trunc,iexact_disc,iexact_ring
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
      common /sa_options/ i_saneg,i_sanoft2,sa_cw2,sa_cw3,sa_cv1,sa_ct3,
     .                    sa_ct4,sa_cb1,sa_cb2,sa_sigma,sa_karman
c

      dimension q(jdim,kdim,idim,5),sj(jdim,kdim,idim-1,5),
     + sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
     +,dtj(jdim,kdim,idim-1),x(jdim,kdim,idim),y(jdim,kdim,idim),
     + z(jdim,kdim,idim),vist3d(jdim,kdim,idim),
     + vor(jdim-1,kdim-1,idim-1)
      dimension damp1(jdim-1,kdim-1,idim-1),
     + fnu(0:jdim,0:kdim,0-iex3:idim+iex3),
     + timestp(jdim-1,kdim-1,idim-1),fd(jdim-1,kdim-1,idim-1)
      dimension bx(kdim-1,jdim-1),cx(kdim-1,jdim-1),dx(kdim-1,jdim-1),
     + fx(kdim-1,jdim-1),workx(kdim-1,jdim-1),
     +          by(jdim-1,kdim-1),cy(jdim-1,kdim-1),dy(jdim-1,kdim-1),
     + fy(jdim-1,kdim-1),worky(jdim-1,kdim-1),
     +          bz(kdim-1,idim-1),cz(kdim-1,idim-1),dz(kdim-1,idim-1),
     + fz(kdim-1,idim-1),workz(kdim-1,idim-1)
      dimension turre(0-iex:jdim+iex,0-iex:kdim+iex,0-iex2:idim+iex2),
     + smin(jdim-1,kdim-1,idim-1),tursav(jdim,kdim,idim,nummem)
     + ,tj0(kdim,idim-1,nummem,4),tk0(jdim,idim-1,nummem,4),
     +  ti0(jdim,kdim,nummem,4),blank(jdim,kdim,idim),
     + qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4)
      dimension tursav2(jdim,kdim,idim,2*nummem)
      dimension volj0(kdim,idim-1,4),
     +          volk0(jdim,idim-1,4),voli0(jdim,kdim,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2),
     + ux(jdim-1,kdim-1,idim-1,9),vx(0:jdim,0:kdim,idim-1,isarcnum)
c
      if(isklton .gt. 0) then
        if (i_saneg .eq. 1) then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''   Using SA-neg form'')')
        end if
        if (i_sanoft2 .eq. 1) then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''   Using SA-noft2 form'')')
        end if
        if (iexact_trunc .ne. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         if (iexact_trunc .eq. 1) then
           write(bou(nou(1),1),'(''   Using exact MMS source in'',
     +      '' Spalart routine to evaluate truncation error, MS1'')')
         else if (iexact_trunc .eq. 2) then
           write(bou(nou(1),1),'(''   Using exact MMS source in'',
     +      '' Spalart routine to evaluate truncation error, MS2'')')
         else if (iexact_trunc .eq. 4) then
           write(bou(nou(1),1),'(''   Using exact MMS source in'',
     +      '' Spalart routine to evaluate truncation error, MS4'')')
         end if
        else
         if (ides .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''   computing turbulent'',
     +       '' viscosity using Spalart DES (cdes='',f7.3,
     +       ''), block ='',i4)') cdes,nbl
         else if (ides .eq. 2) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''   computing turbulent'',
     +       '' viscosity using Spalart DDES (cdes='',f7.3,
     +       ''), block ='',i4)') cdes,nbl
         else if (ides .eq. 3) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''   computing turbulent'',
     +       '' viscosity using Spalart MDDES (cdes='',f7.3,
     +       '',cddes='',f7.3,''), block ='',i4)') cdes,cddes,nbl
         else if (ides .ge. 4 .and. ides .le. 6) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''   computing turbulent'',
     +       '' viscosity using Spalart MuDDES (cdes='',f7.3,
     +       '',cddes='',f7.3,''), block ='',i4)') cdes,cddes,nbl
         else if (ides .ge. 7) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''   computing turbulent'',
     +       '' viscosity using S-A fd on nu (cdes='',f7.3,
     +       '',cddes='',f7.3,''), block ='',i4)') cdes,cddes,nbl
         else
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''   computing turbulent'',
     +       '' viscosity using Spalart, block ='',i4)') nbl
         end if
         if (isarc2d .eq. 1 .or. isarc3d .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      Rotation-curvature'',
     +      '' (SARC) terms are ON; sarccr3='',f7.3)') sarccr3
         end if
         if (isarc2d .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      WARNING: curv terms are'',
     +       '' active in 2-D sense only!!!'')')
         end if
         if (isar .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      Dacles-Mariani rotation'',
     +      '' (SAR) correction is ON, crot='',f7.3)') crot
         end if
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     Freestream tur10 = '',
     +     e19.8)') tur10(1)
         if (iexact_disc .ne. 0) then
           nou(1) = min(nou(1)+1,ibufdim)
           if (iexact_disc .eq. 1) then
             write(bou(nou(1),1),'(''     using exact MMS source in'',
     +     '' Spalart routine to evaluate discretization error, MS1'')')
           else if (iexact_disc .eq. 2) then
             write(bou(nou(1),1),'(''     using exact MMS source in'',
     +     '' Spalart routine to evaluate discretization error, MS2'')')
           else if (iexact_disc .eq. 4) then
             write(bou(nou(1),1),'(''     using exact MMS source in'',
     +     '' Spalart routine to evaluate discretization error, MS4'')')
           end if
         end if
         if(iaxi2planeturb .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     SA model ignoring i-dir'')')
         end if
         if(istrongturbdis .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     strong conserv - diss terms'')')
         end if
        end if
        if(iturbord .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     1st order advection on RHS'')')
        else
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     2nd order advection on RHS'')')
        end if
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   SA constants:'')')
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''     cw2='',f15.5)') sa_cw2
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''     cw3='',f15.5)') sa_cw3
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''     cv1='',f15.5)') sa_cv1
        if (i_sanoft2 .eq. 1) then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''     ct3=0'')')
        else
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''     ct3='',f15.5)') sa_ct3
        end if
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''     ct4='',f15.5)') sa_ct4
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''     cb1='',f15.5)') sa_cb1
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''     cb2='',f15.5)') sa_cb2
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''     sigma='',f15.5)') sa_sigma
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''     akarman='',f15.5)') sa_karman
      end if
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
c   Note: (10.**(-iexp) is machine zero)
      xminn=10.**(-iexp+1)
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
      factor=10.
c
c Overwrite with keyword value if nonzero
c
      if (real(cflturb(1)).ne.0.) then
         factor = cflturb(1)
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
c   *********
c Set up constants
c   defaults: akarman=0.41, cb1=0.1355, sigma=2/3, cb2=0.622, cw2=0.3,
c             cw3=2.0, cv1=7.1, ct3=1.2, ct4=0.5
      akarman=sa_karman
      cb1=sa_cb1
      sigma=sa_sigma
      cb2=sa_cb2
      cw1=cb1/akarman**2+(1.+cb2)/sigma
      cw2=sa_cw2
      cw3=sa_cw3
      cv1=sa_cv1
c     ct1=1.0
c     ct2=2.0
c     ct3=1.1
c     ct4=2.0
      if (i_sanoft2 .eq. 1) then
        ct3=0.0
      else
        ct3=sa_ct3
      end if
      ct4=sa_ct4
c Set up some other needed parameters
      jd2=(jdim-1)/2
      re=reue/xmach
      c2b=cbar/tinf
      c2bp=c2b+1.0
      nnit=nsubit
c
      iwrite=0
c
c Get laminar viscosity (divided by rho) at cell centers
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            tt=gamma*q(j,k,i,5)/q(j,k,i,1)
            fnu(j,k,i)=c2bp*tt*sqrt(tt)/(c2b+tt)/q(j,k,i,1)
          enddo
        enddo
      enddo
      do i=1,idim-1
        do k=1,kdim-1
          tt=gamma*qj0(k,i,5,1)/qj0(k,i,1,1)
          fnu(0,k,i)=c2bp*tt*sqrt(tt)/(c2b+tt)/qj0(k,i,1,1)
          tt=gamma*qj0(k,i,5,3)/qj0(k,i,1,3)
          fnu(jdim,k,i)=c2bp*tt*sqrt(tt)/(c2b+tt)/qj0(k,i,1,3)
        enddo
      enddo
      do i=1,idim-1
        do j=1,jdim-1
          tt=gamma*qk0(j,i,5,1)/qk0(j,i,1,1)
          fnu(j,0,i)=c2bp*tt*sqrt(tt)/(c2b+tt)/qk0(j,i,1,1)
          tt=gamma*qk0(j,i,5,3)/qk0(j,i,1,3)
          fnu(j,kdim,i)=c2bp*tt*sqrt(tt)/(c2b+tt)/qk0(j,i,1,3)
        enddo
      enddo
      if (i2d .ne. 1 .and. iaxi2planeturb .ne. 1) then
      do j=1,jdim-1
        do k=1,kdim-1
          tt=gamma*qi0(j,k,5,1)/qi0(j,k,1,1)
          fnu(j,k,0)=c2bp*tt*sqrt(tt)/(c2b+tt)/qi0(j,k,1,1)
          tt=gamma*qi0(j,k,5,3)/qi0(j,k,1,3)
          fnu(j,k,idim)=c2bp*tt*sqrt(tt)/(c2b+tt)/qi0(j,k,1,3)
        enddo
      enddo
      end if
c If this is 1st global subiteration for time-accurate computation,
c save tursav (at time step n):
c tursav2(j,k,i,1) is turb quantity 
c tursav2(j,k,i,3) is Delta turb quantity
c (in SA, only (j,k,i,1) and (j,k,i,3) are used)
      if (real(dt) .gt. 0. .and. icyc .eq. 1) then
      if (abs(ita) .eq. 2) then
c     if tursav2 at 1st point is zero, then we know that we do not have 
c     2nd order data from the restart; no choice but to set 
c     tursav2(j,k,i,3)=deltaQ=0 for 1st iteration
        if (real(tursav2(1,1,1,1)) .eq. 0.) then
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              tursav2(j,k,i,3)=0.
            enddo
          enddo
        enddo
        else
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              tursav2(j,k,i,3)=tursav(j,k,i,1)-tursav2(j,k,i,1)
            enddo
          enddo
        enddo
        end if
      end if
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              tursav2(j,k,i,1)=tursav(j,k,i,1)
            enddo
          enddo
        enddo
      end if
c Get TURRE values
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,i)=tursav(j,k,i,1)
          enddo
        enddo
      enddo
c Now iterate to solve the equation
      do 500 not=1,nnit
c
c    set up boundary conditions (they are in ghost cells everywhere)
c      (1) k=0 boundary:
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,0,i)=tk0(j,i,1,1)
          enddo
        enddo
c      (2) k=kdim boundary:
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,kdim,i)=tk0(j,i,1,3)
          enddo
        enddo
c      (3) j=0 boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(0,k,i)=tj0(k,i,1,1)
          enddo
        enddo
c      (4) j=jdim boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(jdim,k,i)=tj0(k,i,1,3)
          enddo
        enddo
        if (i2d .ne. 1 .and. iaxi2planeturb .ne. 1) then
c      (5) i=0 boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,0)=ti0(j,k,1,1)
          enddo
        enddo
c      (6) i=idim boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,idim)=ti0(j,k,1,3)
          enddo
        enddo
        end if
        if (iturbord .ne. 1) then
c      (1) k=0 boundary:
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,-1,i)=tk0(j,i,1,2)
          enddo
        enddo
c      (2) k=kdim boundary:
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,kdim+1,i)=tk0(j,i,1,4)
          enddo
        enddo
c      (3) j=0 boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(-1,k,i)=tj0(k,i,1,2)
          enddo
        enddo
c      (4) j=jdim boundary:
        do i=1,idim-1
          do k=1,kdim-1
            turre(jdim+1,k,i)=tj0(k,i,1,4)
          enddo
        enddo
        if (i2d .ne. 1 .and. iaxi2planeturb .ne. 1) then
c      (5) i=0 boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,-1)=ti0(j,k,1,2)
          enddo
        enddo
c      (6) i=idim boundary:
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,idim+1)=ti0(j,k,1,4)
          enddo
        enddo
        end if
        end if
c    **************
c
c    F_eta_eta viscous terms
c       Interior points
        do k=2,kdim-2
          kl=k-1
          ku=k+1
          do i=1,idim-1
            do j=1,jdim-1
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              byy=-ccmaxcr(cdm+cam,0.)
              cyy=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dyy=-ccmaxcr(cdp+cap,0.)
              vist3d(j,k,i)=-byy*turre(j,k-1,i)
     +          -cyy*turre(j,k,i) -dyy*turre(j,k+1,i)
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              byy=-ccmaxcr(cdm+cam,0.)
              cyy=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dyy=-ccmaxcr(cdp+cap,0.)
              vist3d(j,k,i)=-byy*turre(j,k-1,i)
     +          -cyy*turre(j,k,i) -dyy*turre(j,k+1,i)
            enddo
          enddo
c
c       KDIM boundary points
          k=kdim-1
          kl=kdim-2
          ku=kdim-1
          do i=1,idim-1
            do j=1,jdim-1
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              byy=-ccmaxcr(cdm+cam,0.)
              cyy=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dyy=-ccmaxcr(cdp+cap,0.)
              vist3d(j,k,i)=-byy*turre(j,k-1,i)
     +          -cyy*turre(j,k,i) -dyy*turre(j,k+1,i)
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
              vist3d(j,k,i)=vist3d(j,k,i)-uu*(app*(turre(j,k,i)-
     +          turre(j,k-1,i)) + apm*(turre(j,k+1,i)-turre(j,k,i)))
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
              vist3d(j,k,i)=vist3d(j,k,i)-0.5*uu*app*turre(j,k-2,i)
     +                                    +2.*uu*app*turre(j,k-1,i)
     +                                   -1.5*uu*app*turre(j,k,i)
     +                                   +1.5*uu*apm*turre(j,k,i)
     +                                    -2.*uu*apm*turre(j,k+1,i)
     +                                   +0.5*uu*apm*turre(j,k+2,i)
            enddo
          enddo
        enddo
        end if
c
c    F_xi_xi viscous terms
c       interior points
        do j=2,jdim-2
          jl=j-1
          ju=j+1
          do i=1,idim-1
            do k=1,kdim-1
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              bxx=-ccmaxcr(cdm+cam,0.)
              cxx=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dxx=-ccmaxcr(cdp+cap,0.)
              vist3d(j,k,i)=vist3d(j,k,i)-bxx*turre(j-1,k,i)
     +          -cxx*turre(j,k,i) -dxx*turre(j+1,k,i)
            enddo
          enddo
        enddo
c       J0 boundary points
          j=1
          jl=1
          ju=min(2,jdim-1)
          do i=1,idim-1
            do k=1,kdim-1
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              bxx=-ccmaxcr(cdm+cam,0.)
              cxx=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dxx=-ccmaxcr(cdp+cap,0.)
              vist3d(j,k,i)=vist3d(j,k,i)-bxx*turre(j-1,k,i)
     +          -cxx*turre(j,k,i) -dxx*turre(j+1,k,i)
            enddo
          enddo
c       JDIM boundary points 
          j=jdim-1
          jl=jdim-2
          ju=jdim-1
          do i=1,idim-1
            do k=1,kdim-1
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              bxx=-ccmaxcr(cdm+cam,0.)
              cxx=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dxx=-ccmaxcr(cdp+cap,0.)
              vist3d(j,k,i)=vist3d(j,k,i)-bxx*turre(j-1,k,i)
     +          -cxx*turre(j,k,i) -dxx*turre(j+1,k,i)
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
              vist3d(j,k,i)=vist3d(j,k,i)-uu*(app*(turre(j,k,i)-
     +          turre(j-1,k,i)) + apm*(turre(j+1,k,i)-turre(j,k,i)))
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
              vist3d(j,k,i)=vist3d(j,k,i)-0.5*uu*app*turre(j-2,k,i)
     +                                    +2.*uu*app*turre(j-1,k,i)
     +                                   -1.5*uu*app*turre(j,k,i)
     +                                   +1.5*uu*apm*turre(j,k,i)
     +                                    -2.*uu*apm*turre(j+1,k,i)
     +                                   +0.5*uu*apm*turre(j+2,k,i)
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
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                if (i_saneg .eq. 1 .and. real(turre(j,k,i)).lt.0.) then
                 chi_p=anutp/fnup
                 chi_m=anutm/fnum
                 diff_facp=(16.+chi_p*chi_p*chi_p)/
     +                     (16.-chi_p*chi_p*chi_p)
                 diff_facm=(16.+chi_m*chi_m*chi_m)/
     +                     (16.-chi_m*chi_m*chi_m)
                else
                 diff_facp=1.0
                 diff_facm=1.0
                end if
                cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
                bzz=-ccmaxcr(cdm+cam,0.)
                czz=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
                dzz=-ccmaxcr(cdp+cap,0.)
                vist3d(j,k,i)=vist3d(j,k,i)-bzz*turre(j,k,i-1)
     +            -czz*turre(j,k,i) -dzz*turre(j,k,i+1)
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
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                if (i_saneg .eq. 1 .and. real(turre(j,k,i)).lt.0.) then
                 chi_p=anutp/fnup
                 chi_m=anutm/fnum
                 diff_facp=(16.+chi_p*chi_p*chi_p)/
     +                     (16.-chi_p*chi_p*chi_p)
                 diff_facm=(16.+chi_m*chi_m*chi_m)/
     +                     (16.-chi_m*chi_m*chi_m)
                else
                 diff_facp=1.0
                 diff_facm=1.0
                end if
                cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
                bzz=-ccmaxcr(cdm+cam,0.)
                czz=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
                dzz=-ccmaxcr(cdp+cap,0.)
                vist3d(j,k,i)=vist3d(j,k,i)-bzz*turre(j,k,i-1)
     +            -czz*turre(j,k,i) -dzz*turre(j,k,i+1)
              enddo
            enddo
c
c         IDIM boundary points
            i=idim-1
            il=idim-2
            iu=idim-1
            do k=1,kdim-1
              do j=1,jdim-1
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
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                if (i_saneg .eq. 1 .and. real(turre(j,k,i)).lt.0.) then
                 chi_p=anutp/fnup
                 chi_m=anutm/fnum
                 diff_facp=(16.+chi_p*chi_p*chi_p)/
     +                     (16.-chi_p*chi_p*chi_p)
                 diff_facm=(16.+chi_m*chi_m*chi_m)/
     +                     (16.-chi_m*chi_m*chi_m)
                else
                 diff_facp=1.0
                 diff_facm=1.0
                end if
                cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
                bzz=-ccmaxcr(cdm+cam,0.)
                czz=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
                dzz=-ccmaxcr(cdp+cap,0.)
                vist3d(j,k,i)=vist3d(j,k,i)-bzz*turre(j,k,i-1)
     +            -czz*turre(j,k,i) -dzz*turre(j,k,i+1)
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
                vist3d(j,k,i)=vist3d(j,k,i)-uu*(app*(turre(j,k,i)-
     +           turre(j,k,i-1)) + apm*(turre(j,k,i+1)-
     +           turre(j,k,i)))
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
              vist3d(j,k,i)=vist3d(j,k,i)-0.5*uu*app*turre(j,k,i-2)
     +                                    +2.*uu*app*turre(j,k,i-1)
     +                                   -1.5*uu*app*turre(j,k,i)
     +                                   +1.5*uu*apm*turre(j,k,i)
     +                                    -2.*uu*apm*turre(j,k,i+1)
     +                                   +0.5*uu*apm*turre(j,k,i+2)
              enddo
            enddo
          enddo
          end if
        end if
c
c    Add source term to RHS
      if(isklton .gt. 0) then
        if(ilamlo.eq.0 .or. jlamlo.eq.0 .or. klamlo.eq.0) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   block '',i4,'' in S-A turb model'',
     .   '' has no laminar regions'')') nbl
        else
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   block '',i4,'' in S-A turb model -'',
     .   '' laminar region is:'')') nbl
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   i='',i5,'' to'',i5,'', j='',i5,
     .   '' to'',i5,'', k='',i5,'' to'',i5)') ilamlo,ilamhi,jlamlo,
     .   jlamhi,klamlo,klamhi
        if (i_lam_forcezero .eq. 1) then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''    ...forcing vist3d=0'')')
        end if
        end if
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   NOTE:  This particular model '',
     .   ''<<transitions>> on its own, but there is'')')
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   no guarantee that it will '',
     .   ''transition at all.  Check vist3d levels if unsure.'')')
      end if
c     can use damp1 to temporarily store distance variable (different for DES)
c     (later, damp1 is used to store LHS implicit contribution)
        if (ides .eq. 1) then
c       DES
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
              deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
              deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
              deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
              delta = ccmax(deltaj,deltak)
c
c Modification to allow DES to run in 2d
c
              if( i2d .ne. 1 .and. iaxi2planeturb .ne. 1 ) then
                 delta = ccmax(delta,deltai)
              end if
 
              damp1(j,k,i) = ccmin(ccabs(smin(j,k,i)),cdes*delta)
            enddo
          enddo
        enddo
        elseif (ides .ge. 2) then
c       DDES (TCFD 20:181-195, 2006)
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
              deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
              deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
              deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
              delta = ccmax(deltaj,deltak)
c
c Modification to allow DES to run in 2d
c
              if( i2d .ne. 1 .and. iaxi2planeturb .ne. 1 ) then
                 delta = ccmax(delta,deltai)
              end if

              dist = ccabs(smin(j,k,i))
              velterm=ux(j,k,i,1)**2 + ux(j,k,i,2)**2 + ux(j,k,i,3)**2 +
     +                ux(j,k,i,4)**2 + ux(j,k,i,5)**2 + ux(j,k,i,6)**2 +
     +                ux(j,k,i,7)**2 + ux(j,k,i,8)**2 + ux(j,k,i,9)**2
              rd = turre(j,k,i)/(sqrt(velterm)*akarman*akarman*
     +             dist*dist*re)
              fd(j,k,i) = 1.0 - cctanh((8.0*rd)*(8.0*rd)*(8.0*rd))
              if( ides .ne. 7 ) then
                 term = ccmaxrc(0.0,dist-cdes*delta) 
                 damp1(j,k,i) = dist - fd(j,k,i)*term
              else
                 damp1(j,k,i) = ccabs(smin(j,k,i))
              end if
            enddo
          enddo
        enddo
        else
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
              damp1(j,k,i) = ccabs(smin(j,k,i))
            enddo
          enddo
        enddo
        end if
c
c   For SARC model, eventually need to compute rtilde as fn of Dalpha/Dt
        if (isarc2d .eq. 1) then
c   Note 2-D SARC currently only set up for curvature in x-z plane, with
c   y direction (necessarily) in the i-direction
c   It also makes incompressible assumption that S11=-S33
c   Get vx(3)=DS11/Dt and vx(4)=DS13/Dt:
          call sijrate2d(idim,jdim,kdim,q,qj0,qk0,
     .    bcj,bck,vol,sj,sk,vx)
c   compute curvature term, part of Dalpha/Dt: store in vx(1)
c   modify Sij such that its diagonal terms are traceless in the 2D sense
          do i=1,idim-1
            do j=1,jdim-1
              do k=1,kdim-1
                s11 = ux(j,k,i,1)
     +             -(ux(j,k,i,1)+ux(j,k,i,9))/2.
                s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
                vx(j,k,i,1)=s11*vx(j,k,i,4)-s13*vx(j,k,i,3)
              enddo
            enddo
          enddo
        end if
        if (isarc3d .eq. 1) then
c   Newer generalized 3-D SARC
c   Get vx(1)=DS11/Dt, vx(2)=DS12/Dt, vx(3)=DS13/Dt,
c   Get vx(4)=DS22/Dt, vx(5)=DS23/Dt, vx(6)=DS33/Dt:
          call sijrate3d(idim,jdim,kdim,q,ux,vol,si,sj,sk,vx)
        end if
c
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
               if ((i.ge.ilamlo .and. i.lt.ilamhi .and.
     .              j.ge.jlamlo .and. j.lt.jlamhi .and.
     .              k.ge.klamlo .and. k.lt.klamhi) .or.
     .              real(smin(j,k,i)) .lt. 0.) then
                  cutoff=0.
               else if( ides .eq. 3 .and. 
     .          real(fd(j,k,i)) .gt. real(cddes) ) then
                  cutoff = (1.d0 - fd(j,k,i))/(1.d0-cddes)
               else if( ides .eq. 4 .and. 
     .          real(fd(j,k,i)) .gt. real(cddes) ) then
                  cutoff = (1.d0 - fd(j,k,i))/(1.d0-cddes)
               else if( ides .eq. 5 .and. 
     .          real(fd(j,k,i)) .gt. 0.999d0 ) then
                  cutoff = (1.d0 - fd(j,k,i))/(1.d0-0.999d0)
               else
                  cutoff=1.
               end if
              ss=vor(j,k,i)
              if (isar .eq. 1) then
                s11 = ux(j,k,i,1)
                s22 = ux(j,k,i,5)
                s33 = ux(j,k,i,9)
                s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
                s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
                s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
                xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
                xisabs = sqrt(2.*xis)
                ss=vor(j,k,i)+crot*ccminrc(0.,xisabs-vor(j,k,i))
              end if
              chi=turre(j,k,i)/fnu(j,k,i)
              if (i_saneg .ne. 1) then
                fv1=chi**3/(chi**3+cv1**3)
                fv2=1.-(chi/(1.+chi*fv1))
                sst=ss+turre(j,k,i)*fv2/(re*(akarman*damp1(j,k,i))**2)
                sst=ccmax(sst,xminn)
                rr=turre(j,k,i)/(sst*re*(akarman*damp1(j,k,i))**2)
                rr=ccmincr(rr,10.)
                gg=rr+cw2*(rr**6-rr)
                gg=ccmax(gg,xminn)
                fw=gg*((1.+cw3**6)/(gg**6+cw3**6))**(1./6.)
                ft2=ct3*exp(-ct4*chi**2)
              else if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .ge. 0.)
     +          then 
c   Alternative to clipping sst (Allmaras et al., ICCFD7-1902):
                fv1=chi**3/(chi**3+cv1**3)
                fv2=1.-(chi/(1.+chi*fv1))
                sst=ss+turre(j,k,i)*fv2/(re*(akarman*damp1(j,k,i))**2)
                sbar=turre(j,k,i)*fv2/(re*(akarman*damp1(j,k,i))**2)
                if (real(sbar) .ge. -0.7*real(ss)) then
                  sst=ss+sbar
                else
                  sst=ss+(ss*(0.7*0.7*ss+0.9*sbar)/((0.9-2.*0.7)*ss-
     +                sbar))
                end if
                if (sst .eq. 0.) then
                  rr=10.
                else
                  rr=turre(j,k,i)/(sst*re*(akarman*damp1(j,k,i))**2)
                  rr=ccmincr(rr,10.)
                end if
                gg=rr+cw2*(rr**6-rr)
                gg=ccmax(gg,xminn)
                fw=gg*((1.+cw3**6)/(gg**6+cw3**6))**(1./6.)
                ft2=ct3*exp(-ct4*chi**2)
              end if
c   Because we split up sst between term1 and term2, need to limit ss >= 0
c   for case when using isar=1 (for which ss itself may be negative).
c   Normally, with ss=vor, it is always non-negative.
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
                term1=cb1*(1.-ct3)*ss
                term2=cw1
              else
                term1=cb1*(1.-ft2)*ccmaxrc(0.,ss)
                term2=cb1*((1.-ft2)*fv2+ft2)/akarman**2-cw1*fw
c             For SARC, note that fr1 multiplies the production term
c             cb1*Stilde*nuwiggle, as specified in AIAA J 38(5):784,2000,
c             where Stilde=vort+nuwiggle*fv2/(kappa^2*d^2)
c             (in Aerospace Science & Tech. 1(5):297,1997 it
c             multiplies cb1*vort*nuwiggle).
c             Ends up making very little difference.
              if (isarc2d .eq. 1 .or. isarc3d .eq. 1) then
c               Determine Sij and Wij values:
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
                ss=ccmax(ss,xminn)
                xisabs=sqrt(2.*xis)
                rstar=xisabs/ss
                xisabs=ccmax(xisabs,xminn)
                if (isarc2d .eq. 1) then
                  rtilde=-4.*vx(j,k,i,1)*w13/
     +               (0.5*(ss**2+xisabs**2))**2
                else if (isarc3d .eq. 1) then
                  rtilde=2./(0.5*(ss**2+xisabs**2))**2*
     +              ( -w12*vx(j,k,i,2)*(s11-s22)
     +                -w13*vx(j,k,i,3)*(s11-s33)
     +                -w23*vx(j,k,i,5)*(s22-s33)
     +                +s12*(-w12*(vx(j,k,i,4)-vx(j,k,i,1))
     +                      -w13*vx(j,k,i,5)-w23*vx(j,k,i,3))
     +                +s13*(-w13*(vx(j,k,i,6)-vx(j,k,i,1))
     +                      -w12*vx(j,k,i,5)+w23*vx(j,k,i,2))
     +                +s23*(-w23*(vx(j,k,i,6)-vx(j,k,i,4))
     +                      +w12*vx(j,k,i,3)+w13*vx(j,k,i,2)) )
                end if
                fr1=4.*rstar/(1.+rstar)*(1.-sarccr3*
     +             ccatan(12.*rtilde))-1.
                term1=cb1*(fr1-ft2)*ss
                term2=cb1*((fr1-ft2)*fv2+ft2)/akarman**2-cw1*fw
              end if
              end if
              dist2i=1./(re*damp1(j,k,i)**2+1.e-20)
              tt=cutoff*term1*turre(j,k,i)+term2*turre(j,k,i)**2*dist2i
c            Store quantity to be added to certain implicit LHS terms:
              damp1(j,k,i)=2.*term2*turre(j,k,i)*dist2i
              dfv1=(fv1-fv1*fv1)*3./turre(j,k,i)
              dfv2=(fv2-1.)/turre(j,k,i)+((1.-fv2)**2)*(fv1/
     +             turre(j,k,i)+dfv1)
              dft2=-(2.*ct4*turre(j,k,i)/(fnu(j,k,i)**2))*ft2
              drr=rr/turre(j,k,i)-rr*rr*(fv2/turre(j,k,i)+dfv2)
              dgg=(1.-cw2+6.*cw2*(rr**5))*drr
              gg=ccmax(gg,xminn*10.)
c             fix for single prec. via TLNS3D via FUN3D vi Tim Barth
              dfw=((1.+cw3**6)/(gg**6+cw3**6))**(1./6.) -
     .            ((1.+cw3**6)**(1./6.)/((gg**6+cw3**6)**(7./6.)))*gg**6
              dfw=dfw*dgg
c old         dfw=(gg**(-7)/(1.+(cw3**(-6))))*(fw**7)*dgg
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
                continue
              else
                damp1(j,k,i)=damp1(j,k,i)+
     +               dist2i*(turre(j,k,i)**2)*(cb1/(akarman**2)*
     +               (dfv2-ft2*dfv2-fv2*dft2+dft2)-cw1*dfw)
              end if
c            Add to RHS:
              vist3d(j,k,i)=vist3d(j,k,i)+tt
            enddo
          enddo
        enddo
        if (iexact_trunc .ne. 0 .or. iexact_disc .ne. 0) then
c         note vol already included in resid and sign is reverse of mean flow
          call exact_turb_force(jdim,kdim,idim,x,y,z,vol,vist3d,smin,
     +          iexact_trunc,iexact_disc)
          if(icyc.eq.ncyc1(1) .or. icyc.eq.ncyc1(2) .or. 
     +       icyc.eq.ncyc1(3) .or. icyc.eq.ncyc1(4) .or. 
     +       icyc.eq.ncyc1(5)) then
            if(ntime .eq. 1) then
              call exact_turb_res(jdim,kdim,idim,vol,vist3d,
     +                            iexact_trunc,iexact_disc)
            end if
          end if
        end if
c
c      Transition specified by jtran1 and jtran2 (Spalart's method)
c      NOTE:  This method currently is set up to only work in 2-D with
c      jtran1 and jtran2 (lower & upper surface) specified.  It needs
c      a lot of work for 3-D or any non-C-mesh airfoil case with j in
c      streamwise direction.  Also, when this is used, you need to set 
c      tur10=0.1 (rather than 1.341946) in subroutine init.
c         jtran1=jlamlo
c         jtran2=jlamhi
c         jnose=(jtran1+jtran2+1)/2-1
c         do i=1,idim-1
c           do k=1,kdim-1
c             do j=1,jdim-1
c               if(j .le. jnose) then
c                 jtran=jtran1
c               else
c                 jtran=jtran2
c               end if
c               vorttran=vor(jtran,1,i)
c               dxtran=1./sqrt((x(jtran,1,i)-x(jtran+1,1,i))**2+
c    +                         (y(jtran,1,i)-y(jtran+1,1,i))**2+
c    +                         (z(jtran,1,i)-z(jtran+1,1,i))**2)
c               ss=vor(j,k,i)
c               cv2=5.0
c               chi=turre(j,k,i)/fnu(j,k,i)
c               chi=ccmaxcr(chi,0.0001)
c               fv1=chi**3/(chi**3+cv1**3)
c               fv2=1./((1.+chi/cv2)**3)
c               fv3=((1.+chi*fv1)*(1.-fv2))/chi
c             sst=fv3*ss+turre(j,k,i)*fv2/(re*(akarman*smin(j,k,i))**2)
c               sst=ccmax(sst,xminn)
c               rr=turre(j,k,i)/(sst*re*(akarman*smin(j,k,i))**2)
c               rr=ccmincr(rr,10.)
c               gg=rr+cw2*(rr**6-rr)
c               gg=ccmax(gg,xminn)
c               fix for single prec. via TLNS3D via FUN3D via Tim Barth
c               fw=gg*((1.+cw3**6)/(gg**6+cw3**6))**(1./6.)
c old           fw=((gg**(-6)+cw3**(-6))/(1.+cw3**(-6)))**(-1./6.)
c               ft2=ct3*exp(-ct4*chi**2)
c               unorm=sqrt(q(j,k,i,2)**2+q(j,k,i,3)**2+
c    +                     q(j,k,i,4)**2)
c               gt=ccminrc(0.1,unorm/(vorttran+1.e-20)*dxtran)
c               xbody=0.25*(x(jtran,1,i  )+x(jtran+1,1,i  )+
c    +                      x(jtran,1,i+1)+x(jtran+1,1,i+1))
c               ybody=0.25*(y(jtran,1,i  )+y(jtran+1,1,i  )+
c    +                      y(jtran,1,i+1)+y(jtran+1,1,i+1))
c               zbody=0.25*(z(jtran,1,i  )+z(jtran+1,1,i  )+
c    +                      z(jtran,1,i+1)+z(jtran+1,1,i+1))
c               xf=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+x(j  ,k+1,i  )
c    +                   +x(j  ,k  ,i+1)+x(j+1,k+1,i  )+x(j+1,k  ,i+1)
c    +                   +x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
c               yf=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+y(j  ,k+1,i  )
c    +                   +y(j  ,k  ,i+1)+y(j+1,k+1,i  )+y(j+1,k  ,i+1)
c    +                   +y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
c               zf=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+z(j  ,k+1,i  )
c    +                   +z(j  ,k  ,i+1)+z(j+1,k+1,i  )+z(j+1,k  ,i+1)
c    +                   +z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
c               dttran=(xf-xbody)**2+(yf-ybody)**2+(zf-zbody)**2
c               argtur=-ct2*(vorttran/unorm)**2*(smin(j,k,i)**2+
c    +                 gt**2*dttran)
c               ft1=ct1*gt*exp(argtur)
c               term1=cb1*(1.-ft2)*ss*fv3
c               term2=cb1*((1.-ft2)*fv2+ft2)/akarman**2-cw1*fw
c               term3=ft1*unorm**2
c               dist2i=1./(re*smin(j,k,i)**2+1.e-20)
c               tt=term1*turre(j,k,i)+term2*turre(j,k,i)**2*dist2i+
c    +             term3*re
c            Store quantity to be added to certain implicit LHS terms:
c               damp1(j,k,i)=2.*term2*turre(j,k,i)*dist2i
c               dfv1=(fv1-fv1*fv1)*3./turre(j,k,i)
c               dfv2=(fv2-1.)/turre(j,k,i)+((1.-fv2)**2)*(fv1/
c    +               turre(j,k,i)+dfv1)
c               dft2=-(2.*ct4*turre(j,k,i)/(fnu(j,k,i)**2))*ft2
c               dfv3=((1.-fv2)*(fv1+turre(j,k,i)*dfv1)-(fnu(j,k,i)+
c    +               turre(j,k,i)*fv1)*dfv2-fv3)/turre(j,k,i)
c               drr=rr/turre(j,k,i)-rr*rr*(fv2/turre(j,k,i)+dfv2+
c    +               akarman**2*vor(j,k,i)*dfv3/(turre(j,k,i)*dist2i))
c               dgg=(1.-cw2+6.*cw2*(rr**5))*drr
c               gg=ccmax(gg,xminn*10.)
c               fix for single prec. via TLNS3D via FUN3D vi Tim Barth
c               dfw=((1.+cw3**6)/(gg**6+cw3**6))**(1./6.) -
c    .            ((1.+cw3**6)**(1./6.)/((gg**6+cw3**6)**(7./6.)))*gg**6
c               dfw=dfw*dgg
c old           dfw=(gg**(-7)/(1.+(cw3**(-6))))*(fw**7)*dgg
c               damp1(j,k,i)=damp1(j,k,i)+
c    +               dist2i*(turre(j,k,i)**2)*(cb1/(akarman**2)*
c    +               (dfv2-ft2*dfv2-fv2*dft2+dft2)-cw1*dfw)
c            Add to RHS:
c               vist3d(j,k,i)=vist3d(j,k,i)+tt
c             enddo
c           enddo
c         enddo
c
c    Implicit F_eta_eta viscous terms.  Do over all i's
        do i=1,idim-1
c         Interior points
          do k=2,kdim-2
            kl=k-1
            ku=k+1
            do j=1,jdim-1
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              by(j,k)=-ccmaxcr(cdm+cam,0.)
              cy(j,k)=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dy(j,k)=-ccmaxcr(cdp+cap,0.)
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
c      Add part of source terms to diagonal LHS in this sweep
              if(real(damp1(j,k,i)) .lt. 0.) then
                cy(j,k)=cy(j,k) + uu*(app-apm) - damp1(j,k,i)
              else
                cy(j,k)=cy(j,k) + uu*(app-apm)
              end if
              dy(j,k)=dy(j,k) + uu*apm
            enddo
            do j=1,jdim-1
              fact=timestp(j,k,i)
              by(j,k)=by(j,k)*fact
              cy(j,k)=cy(j,k)*fact+1.0*(1.+phi)
              dy(j,k)=dy(j,k)*fact
              fy(j,k)=vist3d(j,k,i)*fact
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+(1.+phi)*(tursav2(j,k,i,1)-turre(j,k,i))
     +                 +phi*tursav2(j,k,i,3)
              enddo
            end if
          enddo
c
c         K0 boundary points
            k=1
            kl=1
            ku=min(2,kdim-1)
            do j=1,jdim-1
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              by(j,k)=-ccmaxcr(cdm+cam,0.)
              cy(j,k)=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dy(j,k)=-ccmaxcr(cdp+cap,0.)
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
c      Add part of source terms to diagonal LHS in this sweep
              if(real(damp1(j,k,i)) .lt. 0.) then
                cy(j,k)=cy(j,k) + uu*(app-apm) - damp1(j,k,i)
              else
                cy(j,k)=cy(j,k) + uu*(app-apm)
              end if
              dy(j,k)=dy(j,k) + uu*apm
            enddo
            do j=1,jdim-1
              fact=timestp(j,k,i)
              by(j,k)=by(j,k)*fact
              cy(j,k)=cy(j,k)*fact+1.0*(1.+phi)
              dy(j,k)=dy(j,k)*fact
              fy(j,k)=vist3d(j,k,i)*fact
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+(1.+phi)*(tursav2(j,k,i,1)-turre(j,k,i))
     +                 +phi*tursav2(j,k,i,3)
              enddo
            end if
c
c         KDIM boundary points
            k=kdim-1
            kl=kdim-2
            ku=kdim-1
            do j=1,jdim-1
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              by(j,k)=-ccmaxcr(cdm+cam,0.)
              cy(j,k)=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dy(j,k)=-ccmaxcr(cdp+cap,0.)
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
c      Add part of source terms to diagonal LHS in this sweep
              if(real(damp1(j,k,i)) .lt. 0.) then
                cy(j,k)=cy(j,k) + uu*(app-apm) - damp1(j,k,i)
              else
                cy(j,k)=cy(j,k) + uu*(app-apm)
              end if
              dy(j,k)=dy(j,k) + uu*apm
            enddo
            do j=1,jdim-1
              fact=timestp(j,k,i)
              by(j,k)=by(j,k)*fact
              cy(j,k)=cy(j,k)*fact+1.0*(1.+phi)
              dy(j,k)=dy(j,k)*fact
              fy(j,k)=vist3d(j,k,i)*fact
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+(1.+phi)*(tursav2(j,k,i,1)-turre(j,k,i))
     +                 +phi*tursav2(j,k,i,3)
              enddo
            end if
          if (iover .eq. 1) then
            do k=1,kdim-1
              do j=1,jdim-1
                fy(j,k)=fy(j,k)*blank(j,k,i)
                by(j,k)=by(j,k)*blank(j,k,i)
                dy(j,k)=dy(j,k)*blank(j,k,i)
                cy(j,k)=cy(j,k)*blank(j,k,i)+(1.-blank(j,k,i))
              enddo
            enddo
          end if
          call triv(jdim-1,kdim-1,1,jdim-1,1,kdim-1,worky,by,cy,dy,fy)
          do k=1,kdim-1
            do j=1,jdim-1
              vist3d(j,k,i)=fy(j,k)
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              bx(k,j)=-ccmaxcr(cdm+cam,0.)
              cx(k,j)=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dx(k,j)=-ccmaxcr(cdp+cap,0.)
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
              cx(k,j)=cx(k,j) + uu*(app-apm)
              dx(k,j)=dx(k,j) + uu*apm
            enddo
            do k=1,kdim-1
              fact=timestp(j,k,i)
              bx(k,j)=bx(k,j)*fact
              cx(k,j)=cx(k,j)*fact+1.0*(1.+phi)
              dx(k,j)=dx(k,j)*fact
              fx(k,j)=vist3d(j,k,i)*(1.+phi)
            enddo
          enddo
c
c         J0 boundary points
            j=1
            jl=1
            ju=min(2,jdim-1)
            do k=1,kdim-1
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              bx(k,j)=-ccmaxcr(cdm+cam,0.)
              cx(k,j)=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dx(k,j)=-ccmaxcr(cdp+cap,0.)
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
              cx(k,j)=cx(k,j) + uu*(app-apm)
              dx(k,j)=dx(k,j) + uu*apm
            enddo
            do k=1,kdim-1
              fact=timestp(j,k,i)
              bx(k,j)=bx(k,j)*fact
              cx(k,j)=cx(k,j)*fact+1.0*(1.+phi)
              dx(k,j)=dx(k,j)*fact
              fx(k,j)=vist3d(j,k,i)*(1.+phi)
            enddo
c
c         JDIM boundary points
            j=jdim-1
            jl=jdim-2
            ju=jdim-1
            do k=1,kdim-1
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
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
               chi_p=anutp/fnup
               chi_m=anutm/fnum
               diff_facp=(16.+chi_p*chi_p*chi_p)/(16.-chi_p*chi_p*chi_p)
               diff_facm=(16.+chi_m*chi_m*chi_m)/(16.-chi_m*chi_m*chi_m)
              else
               diff_facp=1.0
               diff_facm=1.0
              end if
              cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
              bx(k,j)=-ccmaxcr(cdm+cam,0.)
              cx(k,j)=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dx(k,j)=-ccmaxcr(cdp+cap,0.)
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
              cx(k,j)=cx(k,j) + uu*(app-apm)
              dx(k,j)=dx(k,j) + uu*apm
            enddo
            do k=1,kdim-1
              fact=timestp(j,k,i)
              bx(k,j)=bx(k,j)*fact
              cx(k,j)=cx(k,j)*fact+1.0*(1.+phi)
              dx(k,j)=dx(k,j)*fact
              fx(k,j)=vist3d(j,k,i)*(1.+phi)
            enddo
          if (iover .eq. 1) then
            do k=1,kdim-1
              do j=1,jdim-1
                fx(k,j)=fx(k,j)*blank(j,k,i)
                bx(k,j)=bx(k,j)*blank(j,k,i)
                dx(k,j)=dx(k,j)*blank(j,k,i)
                cx(k,j)=cx(k,j)*blank(j,k,i)+(1.-blank(j,k,i))
              enddo
            enddo
          end if
          call triv(kdim-1,jdim-1,1,kdim-1,1,jdim-1,workx,bx,cx,dx,fx)
          do j=1,jdim-1
            do k=1,kdim-1
              vist3d(j,k,i)=fx(k,j)
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
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                if (i_saneg .eq. 1 .and. real(turre(j,k,i)).lt.0.) then
                 chi_p=anutp/fnup
                 chi_m=anutm/fnum
                 diff_facp=(16.+chi_p*chi_p*chi_p)/
     +                     (16.-chi_p*chi_p*chi_p)
                 diff_facm=(16.+chi_m*chi_m*chi_m)/
     +                     (16.-chi_m*chi_m*chi_m)
                else
                 diff_facp=1.0
                 diff_facm=1.0
                end if
                cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
                bz(k,i)=-ccmaxcr(cdm+cam,0.)
                cz(k,i)=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
                dz(k,i)=-ccmaxcr(cdp+cap,0.)
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
                cz(k,i)=cz(k,i) + uu*(app-apm)
                dz(k,i)=dz(k,i) + uu*apm
              enddo
              do k=1,kdim-1
                fact=timestp(j,k,i)
                bz(k,i)=bz(k,i)*fact
                cz(k,i)=cz(k,i)*fact+1.0*(1.+phi)
                dz(k,i)=dz(k,i)*fact
                fz(k,i)=vist3d(j,k,i)*(1.+phi)
              enddo
            enddo
c
c           I0 boundary points
              i=1
              il=1
              iu=min(2,idim-1)
              do k=1,kdim-1
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
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                if (i_saneg .eq. 1 .and. real(turre(j,k,i)).lt.0.) then
                 chi_p=anutp/fnup
                 chi_m=anutm/fnum
                 diff_facp=(16.+chi_p*chi_p*chi_p)/
     +                     (16.-chi_p*chi_p*chi_p)
                 diff_facm=(16.+chi_m*chi_m*chi_m)/
     +                     (16.-chi_m*chi_m*chi_m)
                else
                 diff_facp=1.0
                 diff_facm=1.0
                end if
                cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
                bz(k,i)=-ccmaxcr(cdm+cam,0.)
                cz(k,i)=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
                dz(k,i)=-ccmaxcr(cdp+cap,0.)
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
                cz(k,i)=cz(k,i) + uu*(app-apm)
                dz(k,i)=dz(k,i) + uu*apm
              enddo
              do k=1,kdim-1
                fact=timestp(j,k,i)
                bz(k,i)=bz(k,i)*fact
                cz(k,i)=cz(k,i)*fact+1.0*(1.+phi)
                dz(k,i)=dz(k,i)*fact
                fz(k,i)=vist3d(j,k,i)*(1.+phi)
              enddo
c
c           IDIM boundary points
              i=idim-1
              il=idim-2
              iu=idim-1
              do k=1,kdim-1
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
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                if (i_saneg .eq. 1 .and. real(turre(j,k,i)).lt.0.) then
                 chi_p=anutp/fnup
                 chi_m=anutm/fnum
                 diff_facp=(16.+chi_p*chi_p*chi_p)/
     +                     (16.-chi_p*chi_p*chi_p)
                 diff_facm=(16.+chi_m*chi_m*chi_m)/
     +                     (16.-chi_m*chi_m*chi_m)
                else
                 diff_facp=1.0
                 diff_facm=1.0
                end if
                cdp=(fnup+(diff_facp+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(diff_facm+cb2)*anutm)*ttm/(sigma*re)
                bz(k,i)=-ccmaxcr(cdm+cam,0.)
                cz(k,i)=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
                dz(k,i)=-ccmaxcr(cdp+cap,0.)
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
                cz(k,i)=cz(k,i) + uu*(app-apm)
                dz(k,i)=dz(k,i) + uu*apm
              enddo
              do k=1,kdim-1
                fact=timestp(j,k,i)
                bz(k,i)=bz(k,i)*fact
                cz(k,i)=cz(k,i)*fact+1.0*(1.+phi)
                dz(k,i)=dz(k,i)*fact
                fz(k,i)=vist3d(j,k,i)*(1.+phi)
              enddo
            if (iover .eq. 1) then
              do i=1,idim-1
                do k=1,kdim-1
                  fz(k,i)=fz(k,i)*blank(j,k,i)
                  bz(k,i)=bz(k,i)*blank(j,k,i)
                  dz(k,i)=dz(k,i)*blank(j,k,i)
                  cz(k,i)=cz(k,i)*blank(j,k,i)+(1.-blank(j,k,i))
                enddo
              enddo
            end if
            call triv(kdim-1,idim-1,1,kdim-1,1,idim-1,workz,bz,
     +                cz,dz,fz)
            do i=1,idim-1
              do k=1,kdim-1
                vist3d(j,k,i)=fz(k,i)
              enddo
            enddo
          enddo
        end if
c
c    Update TURRE
        if (iexact_ring .eq. 1) then
          if (isklton.gt.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''   Setting SA res=0 in'',
     +         '' outer ring'')')
          end if
          call zero_resid_ring(jdim,kdim,idim,vist3d,jdim,kdim,idim,
     +      1,2,iexact_trunc,iexact_disc)
        end if
        sumn=0.
        negn=0
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              if(i_saneg .ne. 1 .and.
     +          (real(turre(j,k,i)+vist3d(j,k,i))) .lt. 1.e-12) then
                negn=negn+1
                turre(j,k,i)=1.e-12
                vist3d(j,k,i)=0.
              else if (i_saneg .eq. 1 .and.
     +          (real(turre(j,k,i)+vist3d(j,k,i))) .lt. 0.) then
                negn=negn+1
                turre(j,k,i)=turre(j,k,i)+vist3d(j,k,i)
              else
                turre(j,k,i)=turre(j,k,i)+vist3d(j,k,i)
              end if
              sumn=sumn+vist3d(j,k,i)**2
            enddo
          enddo
        enddo
c  
        sumn=sqrt(sumn)/float((kdim-1)*(jdim-1)*(idim-1))
c
        if(iwrite .eq. 1) then
c         write(15,'('' icyc, it, log res, negn'',
c    +     '', block'',/,2i5,e15.5,2i5)') icyc,not,log10(realsumn)),
c    +     negn,nbl
        end if
c
 500  continue
      sumn1 = sumn
      sumn2 = 1.
      negn1 = negn
      negn2 = 0
c
      if (iexact_trunc .ne. 0) then
c       need to load exact solution again (lost it in vist3d
c       because of memory re-use)
        call exact_turb_q(jdim,kdim,idim,x,y,z,tursav,smin,vist3d,
     +      iexact_trunc,iexact_disc)
      else
c Update VIST3D and save TURRE in TURSAV
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            if (i_saneg .eq. 1 .and. real(turre(j,k,i)) .lt. 0.) then
              vist3d(j,k,i)=0.d0
            else
              chi=turre(j,k,i)/fnu(j,k,i)
              fv1=chi**3/(chi**3+cv1**3)
              vist3d(j,k,i)=fv1*turre(j,k,i)*q(j,k,i,1)
            end if
            if( ides .ge. 4 .and. real(fd(j,k,i)) .gt. real(cddes)) then
               vist3d(j,k,i) = vist3d(j,k,i)
     $              * (1.d0 - fd(j,k,i))/(1.d0-cddes)
            end if
            tursav(j,k,i,1)=turre(j,k,i)
          enddo
        enddo
      enddo
      end if
      if (iexact_ring .eq. 1) then
c       overwrite ring of turbulence values in outer 2 rows of grid
        call exact_turb_q_ring(jdim,kdim,idim,x,y,z,tursav,smin,vist3d,
     +      iexact_trunc,iexact_disc)
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
c
c-----BEGIN HARDWIRED OUTPUT SECTION
c
c  If desired, write out turb index here during last iteration:
c  This currently writes out assuming body is at k=1;
c  also, only i=1 index is written
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
        i=1
        write(92,'(''variables="x","i_t"'')')
        do j=1,jdim-1
          zit=turre(j,1,i)/(0.41*sqrt(fnu(j,1,i)*vor(j,1,i))*
     +        smin(j,1,i)*sqrt(re))
          write(92,'(2e15.5)') 0.5*(x(j,1,i)+x(j+1,1,i)),zit
        enddo
      end if
      end if
      end if
c-----END HARDWIRED OUTPUT SECTION
c
      return
      end
