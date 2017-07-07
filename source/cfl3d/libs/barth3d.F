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
      subroutine barth3d(jdim,kdim,idim,q,sj,sk,si,vol,dtj,
     + x,y,z,vist3d,vor,smin,xjb,tursav,xkb,turre,damp1,damp2,
     + timestp,fnu,
     + bx,cx,dx,fx,workx,by,cy,dy,fy,worky,bz,cz,dz,fz,workz,ntime,
     + tj0,tk0,ti0,nbl,blnum,blank,iover,sumn1,sumn2,negn1,negn2,
     + tursav2,volj0,volk0,voli0,nou,bou,nbuf,ibufdim,iex,iex2,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute turbulent viscosity distributions using
c     1-equation Baldwin-Barth turbulence model (ivisc=4)
c     (reprogrammed by Rumsey for CFL3D 7/91)
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
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /lam/ ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi,
     .        i_lam_forcezero
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /twod/ i2d
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
c
      dimension q(jdim,kdim,idim,5),sj(jdim,kdim,idim-1,5),
     + sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
     +,dtj(jdim,kdim,idim-1),x(jdim,kdim,idim),y(jdim,kdim,idim),
     + z(jdim,kdim,idim),vist3d(jdim,kdim,idim),
     + vor(jdim-1,kdim-1,idim-1)
      dimension damp1(jdim-1,kdim-1,idim-1),damp2(jdim-1,kdim-1,idim-1),
     + fnu(jdim-1,kdim-1,idim-1)
     + ,timestp(jdim-1,kdim-1,idim-1)
      dimension bx(kdim-1,jdim-1),cx(kdim-1,jdim-1),dx(kdim-1,jdim-1),
     + fx(kdim-1,jdim-1),workx(kdim-1,jdim-1),
     +          by(jdim-1,kdim-1),cy(jdim-1,kdim-1),dy(jdim-1,kdim-1),
     + fy(jdim-1,kdim-1),worky(jdim-1,kdim-1),
     +          bz(kdim-1,idim-1),cz(kdim-1,idim-1),dz(kdim-1,idim-1),
     + fz(kdim-1,idim-1),workz(kdim-1,idim-1)
      dimension turre(0-iex:jdim+iex,0-iex:kdim+iex,0-iex2:idim+iex2),
     + smin(jdim-1,kdim-1,idim-1),tursav(jdim,kdim,idim,nummem),
     + tj0(kdim,idim-1,nummem,4),tk0(jdim,idim-1,nummem,4),
     + ti0(jdim,kdim,nummem,4),
     + xjb(jdim-1,kdim-1,idim-1),xkb(jdim-1,kdim-1,idim-1),
     + blnum(jdim-1,kdim-1,idim-1),blank(jdim,kdim,idim)
      dimension tursav2(jdim,kdim,idim,2*nummem)
      dimension volj0(kdim,idim-1,4),
     +          volk0(jdim,idim-1,4),voli0(jdim,kdim,4)
c
      if (isklton .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''   computing turbulent'',
     +   '' viscosity using Baldwin-Barth, block ='',i4)') nbl
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     Freestream tur10 = '',
     +     e19.8)') tur10(1)
        if(iturbord .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     1st order advection on RHS'')')
        else
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     2nd order advection on RHS'')')
        end if
        if(iaxi2planeturb .eq. 1) then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''     BB model ignoring i-dir'')')
        end if
        if(istrongturbdis .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     strong conserv - diss terms'')')
        end if
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
c Determine timestep for turb model from global CFL if steady-state computation,
c                                   from CFL=5 if time-accurate computation
c Set number of subiterations to solve turbulence field eqn per iteration
c (usually, 1 is sufficient... but if residual diverges then may need more)
c
      nsubit=1
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
c Set up constants for the NASA TM-102847 model
      akarman=.41
      cmu=.09
      c1e=1.20
      c2e=2.00
      ratt=c1e/c2e
      sige=akarman**2/((c2e-c1e)*sqrt(cmu))
      sigr=sige
      cr1=1.0
      cr2=.01
      cmax=1.0
      b1=1.0
      aplus1=26.
      aplus2=10.
c Set up some other needed parameters
      jd2=(jdim-1)/2
      re=reue/xmach
      c2b=cbar/tinf
      c2bp=c2b+1.0
c
      iwrite=0
c
      nnit=nsubit
c Get laminar viscosity (divided by rho) at cell centers
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            tt=gamma*q(j,k,i,5)/q(j,k,i,1)
            fnu(j,k,i)=c2bp*tt*sqrt(tt)/(c2b+tt)/q(j,k,i,1)
          enddo
        enddo
      enddo
c Calculate the damping factors at cell centers (vorticity, viscosity,
c and density at wall assumed to be equal to values in cell center above wall)
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
c   Set ypls = big number if nearest body is NOT in current block
c   ***NOTE: This is a crude approximation which should be ok as long
c            as block boundaries are not VERY near bodies.  There is
c            really no EASY way to get correct information from across
c            blocks for ypls
            nblb=int(blnum(j,k,i)+0.1)
            if (nblb .ne. nbl) then
              ypls=1000.
            else
              jbb=int(xjb(j,k,i)+0.1)
              kbb=int(xkb(j,k,i)+0.1)
              ibb=int(tursav(j,k,i,2)+0.1)
              if (jbb .eq. jdim) jbb=jbb-1
              if (kbb .eq. kdim) kbb=kbb-1
              if (ibb .eq. idim) ibb=ibb-1
c             do not allow ibb=0
              ibb=max(ibb,1)
              ssw=vor(jbb,kbb,ibb)
              wnu=fnu(jbb,kbb,ibb)
              ra=sqrt(re*(ssw/wnu + 1.e-10*xmach**2))
              ypls=ra*ccabs(smin(j,k,i))
              ypls=ccmaxcr(ypls,.0001)
            end if
            exp1=exp(-(ypls/aplus1))
            exp2=exp(-(ypls/aplus2))
            damp1(j,k,i)=(1.-exp1)*(1.-exp2)
            ddy=exp1/aplus1*(1.-exp2)+exp2/aplus2*(1.-exp1)
            dd=damp1(j,k,i)
            sdd=sqrt(dd)
            damp2(j,k,i)=ratt+(1.-ratt)*(1./(akarman*ypls)+dd)*
     +                   (sdd+ypls*ddy/sdd)
          enddo
        enddo
      enddo
c If this is 1st global subiteration for time-accurate computation, 
c save tursav (at time step n):
c tursav2(j,k,i,1) is turb quantity
c tursav2(j,k,i,3) is Delta turb quantity
c (in BB, only (j,k,i,1) and (j,k,i,3) are used)
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
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
c      (1) k=0 boundary:
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,0,i)=tk0(j,i,1,1)
          enddo
        enddo
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
c      (2) k=kdim boundary:
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
        do i=1,idim-1
          do j=1,jdim-1
            turre(j,kdim,i)=tk0(j,i,1,3)
          enddo
        enddo
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
c      (3) j=0 boundary:
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
        do i=1,idim-1
          do k=1,kdim-1
            turre(0,k,i)=tj0(k,i,1,1)
          enddo
        enddo
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
c      (4) j=jdim boundary:
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
        do i=1,idim-1
          do k=1,kdim-1
            turre(jdim,k,i)=tj0(k,i,1,3)
          enddo
        enddo
        if (i2d .ne. 1 .and. iaxi2planeturb .ne. 1) then
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
c      (5) i=0 boundary:
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
        do k=1,kdim-1
          do j=1,jdim-1
            turre(j,k,0)=ti0(j,k,1,1)
          enddo
        enddo
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
c      (6) i=idim boundary:
c &*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j,k+1,i)+turre(j,k,i))
              trem=.5*(turre(j,k-1,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/
     +            (sige*re)
              byy=ccmincr(-cdm+cam,0.)
              cyy=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dyy=ccmincr(-cdp+cap,0.)
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j,k+1,i)+turre(j,k,i))
              trem=.5*(turre(j,k-1,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/
     +            (sige*re)
              byy=ccmincr(-cdm+cam,0.)
              cyy=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dyy=ccmincr(-cdp+cap,0.)
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j,k+1,i)+turre(j,k,i))
              trem=.5*(turre(j,k-1,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/
     +            (sige*re)
              byy=ccmincr(-cdm+cam,0.)
              cyy=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dyy=ccmincr(-cdp+cap,0.)
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
              uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
              uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
c       Interior points
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j+1,k,i)+turre(j,k,i))
              trem=.5*(turre(j-1,k,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/
     +            (sige*re)
              bxx=ccmincr(-cdm+cam,0.)
              cxx=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dxx=ccmincr(-cdp+cap,0.)
              vist3d(j,k,i)=vist3d(j,k,i)-bxx*turre(j-1,k,i)
     +          -cxx*turre(j,k,i) -dxx*turre(j+1,k,i)
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j+1,k,i)+turre(j,k,i))
              trem=.5*(turre(j-1,k,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/
     +            (sige*re)
              bxx=ccmincr(-cdm+cam,0.)
              cxx=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dxx=ccmincr(-cdp+cap,0.)
              vist3d(j,k,i)=vist3d(j,k,i)-bxx*turre(j-1,k,i)
     +          -cxx*turre(j,k,i) -dxx*turre(j+1,k,i)
            enddo
          enddo
c
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j+1,k,i)+turre(j,k,i))
              trem=.5*(turre(j-1,k,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/
     +            (sige*re)
              bxx=ccmincr(-cdm+cam,0.)
              cxx=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dxx=ccmincr(-cdp+cap,0.)
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
              uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
              uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
                cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +               (2./sige))/re
                cdp=ttp*cnud
                cdm=ttm*cnud
                trep=.5*(turre(j,k,i+1)+turre(j,k,i))
                trem=.5*(turre(j,k,i-1)+turre(j,k,i))
                cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/
     +              (sige*re)
                cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/
     +              (sige*re)
                bzz=ccmincr(-cdm+cam,0.)
                czz=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
                dzz=ccmincr(-cdp+cap,0.)
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
                cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +               (2./sige))/re
                cdp=ttp*cnud
                cdm=ttm*cnud
                trep=.5*(turre(j,k,i+1)+turre(j,k,i))
                trem=.5*(turre(j,k,i-1)+turre(j,k,i))
                cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/
     +              (sige*re)
                cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/
     +              (sige*re)
                bzz=ccmincr(-cdm+cam,0.)
                czz=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
                dzz=ccmincr(-cdp+cap,0.)
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
                cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +               (2./sige))/re
                cdp=ttp*cnud
                cdm=ttm*cnud
                trep=.5*(turre(j,k,i+1)+turre(j,k,i))
                trem=.5*(turre(j,k,i-1)+turre(j,k,i))
                cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/
     +              (sige*re)
                cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/
     +              (sige*re)
                bzz=ccmincr(-cdm+cam,0.)
                czz=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
                dzz=ccmincr(-cdp+cap,0.)
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
                uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
                uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
c    include check for transition - to CRUDELY simulate trip strips
      if(isklton .gt. 0) then
        if(ilamlo.eq.0 .or. jlamlo.eq.0 .or. klamlo.eq.0) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   block '',i4,'' in B-B turb model'',
     .   '' has no laminar regions'')') nbl
        else
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   block '',i4,'' in B-B turb model -'',
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
        write(bou(nou(1),1),'(''   NOTE:  This particular model'',
     .   '' <<transitions>> on its own, but there is'')')
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   no guarantee that it will'',
     .   '' transition at all.  Check vist3d levels if unsure.'')')
      end if
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
              vist3d(j,k,i)=vist3d(j,k,i)+(c2e*damp2(j,k,i)-c1e)*
     +          sqrt(cmu*damp1(j,k,i))*vor(j,k,i)*cutoff*turre(j,k,i)
            enddo
          enddo
        enddo
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j,k+1,i)+turre(j,k,i))
              trem=.5*(turre(j,k-1,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/
     +            (sige*re)
              by(j,k)=ccmincr(-cdm+cam,0.)
              cy(j,k)=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dy(j,k)=ccmincr(-cdp+cap,0.)
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
              uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
              by(j,k)=by(j,k) - uu*app
              cy(j,k)=cy(j,k) + uu*(app-apm)
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j,k+1,i)+turre(j,k,i))
              trem=.5*(turre(j,k-1,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/
     +            (sige*re)
              by(j,k)=ccmincr(-cdm+cam,0.)
              cy(j,k)=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dy(j,k)=ccmincr(-cdp+cap,0.)
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
              uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
              by(j,k)=by(j,k) - uu*app
              cy(j,k)=cy(j,k) + uu*(app-apm)
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
              volku=vol(j,i,3)
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j,k+1,i)+turre(j,k,i))
              trem=.5*(turre(j,k-1,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,ku,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,kl,i))*ttm/
     +            (sige*re)
              by(j,k)=ccmincr(-cdm+cam,0.)
              cy(j,k)=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dy(j,k)=ccmincr(-cdp+cap,0.)
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
              uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
              sgnu=ccsignrc(1.,uu)
              app=0.5*(1.+sgnu)
              apm=0.5*(1.-sgnu)
              by(j,k)=by(j,k) - uu*app
              cy(j,k)=cy(j,k) + uu*(app-apm)
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j+1,k,i)+turre(j,k,i))
              trem=.5*(turre(j-1,k,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/
     +            (sige*re)
              bx(k,j)=ccmincr(-cdm+cam,0.)
              cx(k,j)=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dx(k,j)=ccmincr(-cdp+cap,0.)
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
              uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j+1,k,i)+turre(j,k,i))
              trem=.5*(turre(j-1,k,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/
     +            (sige*re)
              bx(k,j)=ccmincr(-cdm+cam,0.)
              cx(k,j)=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dx(k,j)=ccmincr(-cdp+cap,0.)
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
              uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
              cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +             (2./sige))/re
              cdp=ttp*cnud
              cdm=ttm*cnud
              trep=.5*(turre(j+1,k,i)+turre(j,k,i))
              trem=.5*(turre(j-1,k,i)+turre(j,k,i))
              cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(ju,k,i))*ttp/
     +            (sige*re)
              cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(jl,k,i))*ttm/
     +            (sige*re)
              bx(k,j)=ccmincr(-cdm+cam,0.)
              cx(k,j)=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
              dx(k,j)=ccmincr(-cdp+cap,0.)
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
              uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
                cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +               (2./sige))/re
                cdp=ttp*cnud
                cdm=ttm*cnud
                trep=.5*(turre(j,k,i+1)+turre(j,k,i))
                trem=.5*(turre(j,k,i-1)+turre(j,k,i))
                cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/
     +              (sige*re)
                cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/
     +              (sige*re)
                bz(k,i)=ccmincr(-cdm+cam,0.)
                cz(k,i)=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
                dz(k,i)=ccmincr(-cdp+cap,0.)
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
                uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
              il=idim-2
              iu=idim-1
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
                cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +               (2./sige))/re
                cdp=ttp*cnud
                cdm=ttm*cnud
                trep=.5*(turre(j,k,i+1)+turre(j,k,i))
                trem=.5*(turre(j,k,i-1)+turre(j,k,i))
                cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/
     +              (sige*re)
                cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/
     +              (sige*re)
                bz(k,i)=ccmincr(-cdm+cam,0.)
                cz(k,i)=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
                dz(k,i)=ccmincr(-cdp+cap,0.)
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
                uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
                cnud=(fnu(j,k,i)+cmu*damp1(j,k,i)*turre(j,k,i)*
     +               (2./sige))/re
                cdp=ttp*cnud
                cdm=ttm*cnud
                trep=.5*(turre(j,k,i+1)+turre(j,k,i))
                trem=.5*(turre(j,k,i-1)+turre(j,k,i))
                cap=cmu*trep*0.5*(damp1(j,k,i)+damp1(j,k,iu))*ttp/
     +              (sige*re)
                cam=cmu*trem*0.5*(damp1(j,k,i)+damp1(j,k,il))*ttm/
     +              (sige*re)
                bz(k,i)=ccmincr(-cdm+cam,0.)
                cz(k,i)=ccmaxcr(cdp-cap,0.) + ccmaxcr(cdm-cam,0.)
                dz(k,i)=ccmincr(-cdp+cap,0.)
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
                uu=b1*(xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4))+tc
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
        sumn=0.
        negn=0
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              if((real(turre(j,k,i)+vist3d(j,k,i))) .lt. 1.e-12) then
                negn=negn+1
                turre(j,k,i)=1.e-12
                vist3d(j,k,i)=0.
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
c    +     '', block'',/,2i5,e15.5,2i5)') icyc,not,log10(sumn),
c    +     negn,nbl
        end if
c
 500  continue
      sumn1 = sumn
      sumn2 = 1.
      negn1 = negn
      negn2 = 0
c
c Update VIST3D
        do i=1,idim-1
          do k=1,kdim-1
            do j=1,jdim-1
              tursav(j,k,i,1)=turre(j,k,i)
              vist3d(j,k,i)=cmu*damp1(j,k,i)*turre(j,k,i)*q(j,k,i,1)
            enddo
          enddo
        enddo
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
      return
      end
