      subroutine spalart(jdim,kdim,idim,q,sj,sk,si,vol,dtj,
     + x,y,z,vist3d,vor,smin,tursav,turre,damp1,timestp,fnu,bx,cx,dx,
     + fx,workx,by,cy,dy,fy,worky,bz,cz,dz,fz,workz,ntime,
     + tj0,tk0,ti0,nbl,qj0,qk0,qi0,blank,iover,sumn1,sumn2,negn1,negn2,
     + tursav2,volj0,volk0,voli0,nou,bou,nbuf,ibufdim)
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
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /lam/ ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /twod/ i2d
      common /zero/ iexp
      common /turbconv/ cflturb,edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy
c
      dimension q(jdim,kdim,idim,5),sj(jdim,kdim,idim-1,5),
     + sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
     +,dtj(jdim,kdim,idim-1),x(jdim,kdim,idim),y(jdim,kdim,idim),
     + z(jdim,kdim,idim),vist3d(jdim,kdim,idim),
     + vor(jdim-1,kdim-1,idim-1)
      dimension damp1(jdim-1,kdim-1,idim-1),fnu(0:jdim,0:kdim,0:idim)
     + ,timestp(jdim-1,kdim-1,idim-1)
      dimension bx(kdim-1,jdim-1),cx(kdim-1,jdim-1),dx(kdim-1,jdim-1),
     + fx(kdim-1,jdim-1),workx(kdim-1,jdim-1),
     +          by(jdim-1,kdim-1),cy(jdim-1,kdim-1),dy(jdim-1,kdim-1),
     + fy(jdim-1,kdim-1),worky(jdim-1,kdim-1),
     +          bz(kdim-1,idim-1),cz(kdim-1,idim-1),dz(kdim-1,idim-1),
     + fz(kdim-1,idim-1),workz(kdim-1,idim-1)
      dimension turre(0:jdim,0:kdim,0:idim),
     + smin(jdim-1,kdim-1,idim-1),tursav(jdim,kdim,idim,2)
     + ,tj0(kdim,idim-1,2,4),tk0(jdim,idim-1,2,4),
     +  ti0(jdim,kdim,2,4),blank(jdim,kdim,idim),
     + qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4)
      dimension tursav2(jdim,kdim,idim,2)
      dimension volj0(kdim,idim-1,4),
     +          volk0(jdim,idim-1,4),voli0(jdim,kdim,4)
c
      if(isklton .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''   computing turbulent'',
     +'' viscosity using Spalart, block ='',i4)') nbl
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
c   *********
c Set up constants
      akarman=.41
      cb1=.1355
      sigma=2./3.
      cb2=.622
      cw1=cb1/akarman**2+(1.+cb2)/sigma
      cw2=.30
      cw3=2.0
      cv1=7.1
      ct1=1.0
      ct2=2.0
c     ct3=1.1
c     ct4=2.0
      ct3=1.2
      ct4=0.5
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
      do j=1,jdim-1
        do k=1,kdim-1
          tt=gamma*qi0(j,k,5,1)/qi0(j,k,1,1)
          fnu(j,k,0)=c2bp*tt*sqrt(tt)/(c2b+tt)/qi0(j,k,1,1)
          tt=gamma*qi0(j,k,5,3)/qi0(j,k,1,3)
          fnu(j,k,idim)=c2bp*tt*sqrt(tt)/(c2b+tt)/qi0(j,k,1,3)
        enddo
      enddo
c If this is 1st global subiteration for time-accurate computation,
c save tursav (at time step n):
      if (real(dt) .gt. 0. .and. icyc .eq. 1) then
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
              byy=-ccmaxcr(cdm+cam,0.)
              cyy=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dyy=-ccmaxcr(cdp+cap,0.)
              vist3d(j,k,i)=-byy*turre(j,k-1,i)
     +          -cyy*turre(j,k,i) -dyy*turre(j,k+1,i)
            enddo
          enddo
c    Advective terms in eta
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
              bxx=-ccmaxcr(cdm+cam,0.)
              cxx=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
              dxx=-ccmaxcr(cdp+cap,0.)
              vist3d(j,k,i)=vist3d(j,k,i)-bxx*turre(j-1,k,i)
     +          -cxx*turre(j,k,i) -dxx*turre(j+1,k,i)
            enddo
          enddo
c    Advective terms in xi
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
c
c    F_zeta_zeta viscous terms
        if(i2d .ne. 1) then
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
                bzz=-ccmaxcr(cdm+cam,0.)
                czz=ccmaxcr(cdp+cap,0.) + ccmaxcr(cdm+cam,0.)
                dzz=-ccmaxcr(cdp+cap,0.)
                vist3d(j,k,i)=vist3d(j,k,i)-bzz*turre(j,k,i-1)
     +            -czz*turre(j,k,i) -dzz*turre(j,k,i+1)
              enddo
            enddo
c    Advective terms in zeta
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
        end if
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   NOTE:  This particular model '',
     .   ''<<transitions>> on its own, but there is'')')
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   no guarantee that it will '',
     .   ''transition at all.  Check vist3d levels if unsure.'')')
      end if
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
              if (i.ge.ilamlo .and. i.lt.ilamhi .and.
     .            j.ge.jlamlo .and. j.lt.jlamhi .and.
     .            k.ge.klamlo .and. k.lt.klamhi) then
                cutoff=0.
              else
                cutoff=1.
              end if
              ss=vor(j,k,i)
              chi=turre(j,k,i)/fnu(j,k,i)
              fv1=chi**3/(chi**3+cv1**3)
              fv2=1.-(chi/(1.+chi*fv1))
              sst=ss+turre(j,k,i)*fv2/(re*(akarman*smin(j,k,i))**2)
              sst=ccmax(sst,xminn)
              rr=turre(j,k,i)/(sst*re*(akarman*smin(j,k,i))**2)
              rr=ccmincr(rr,10.)
              gg=rr+cw2*(rr**6-rr)
              gg=ccmax(gg,xminn)
              fw=((gg**(-6)+cw3**(-6))/(1.+cw3**(-6)))**(-1./6.)
              ft2=ct3*exp(-ct4*chi**2)
              term1=cb1*(1.-ft2)*ss
              term2=cb1*((1.-ft2)*fv2+ft2)/akarman**2-cw1*fw
              dist2i=1./(re*smin(j,k,i)**2+1.e-20)
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
              dfw=(gg**(-7)/(1.+(cw3**(-6))))*(fw**7)*dgg
              damp1(j,k,i)=damp1(j,k,i)+
     +             dist2i*(turre(j,k,i)**2)*(cb1/(akarman**2)*
     +             (dfv2-ft2*dfv2-fv2*dft2+dft2)-cw1*dfw)
c            Add to RHS:
              vist3d(j,k,i)=vist3d(j,k,i)+tt
            enddo
          enddo
        enddo
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
c               fw=((gg**(-6)+cw3**(-6))/(1.+cw3**(-6)))**(-1./6.)
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
c               dfw=(gg**(-7)/(1.+(cw3**(-6))))*(fw**7)*dgg
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
              cy(j,k)=cy(j,k)*fact+1.0
              dy(j,k)=dy(j,k)*fact
              fy(j,k)=vist3d(j,k,i)*fact
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+tursav2(j,k,i,1)-turre(j,k,i)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
              cy(j,k)=cy(j,k)*fact+1.0
              dy(j,k)=dy(j,k)*fact
              fy(j,k)=vist3d(j,k,i)*fact
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+tursav2(j,k,i,1)-turre(j,k,i)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j,k+1,i)+turre(j,k,i))
              anutm=.5*(turre(j,k-1,i)+turre(j,k,i))
              fnup=.5*(fnu(j,k+1,i)+fnu(j,k,i))
              fnum=.5*(fnu(j,k-1,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
              cy(j,k)=cy(j,k)*fact+1.0
              dy(j,k)=dy(j,k)*fact
              fy(j,k)=vist3d(j,k,i)*fact
            enddo
            if (real(dt) .gt. 0.) then
              do j=1,jdim-1
                fy(j,k)=fy(j,k)+tursav2(j,k,i,1)-turre(j,k,i)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
              cx(k,j)=cx(k,j)*fact+1.0
              dx(k,j)=dx(k,j)*fact
              fx(k,j)=vist3d(j,k,i)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
              cx(k,j)=cx(k,j)*fact+1.0
              dx(k,j)=dx(k,j)*fact
              fx(k,j)=vist3d(j,k,i)
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
              ttp=xp*xa+yp*ya+zp*za
              ttm=xm*xa+ym*ya+zm*za
              cnud=-cb2*turre(j,k,i)/(sigma*re)
              cap=ttp*cnud
              cam=ttm*cnud
              anutp=.5*(turre(j+1,k,i)+turre(j,k,i))
              anutm=.5*(turre(j-1,k,i)+turre(j,k,i))
              fnup=.5*(fnu(j+1,k,i)+fnu(j,k,i))
              fnum=.5*(fnu(j-1,k,i)+fnu(j,k,i))
              cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
              cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
              cx(k,j)=cx(k,j)*fact+1.0
              dx(k,j)=dx(k,j)*fact
              fx(k,j)=vist3d(j,k,i)
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
        if(i2d .ne. 1) then
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
                cz(k,i)=cz(k,i)*fact+1.0
                dz(k,i)=dz(k,i)*fact
                fz(k,i)=vist3d(j,k,i)
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
                cz(k,i)=cz(k,i)*fact+1.0
                dz(k,i)=dz(k,i)*fact
                fz(k,i)=vist3d(j,k,i)
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
                ttp=xp*xa+yp*ya+zp*za
                ttm=xm*xa+ym*ya+zm*za
                cnud=-cb2*turre(j,k,i)/(sigma*re)
                cap=ttp*cnud
                cam=ttm*cnud
                anutp=.5*(turre(j,k,i+1)+turre(j,k,i))
                anutm=.5*(turre(j,k,i-1)+turre(j,k,i))
                fnup=.5*(fnu(j,k,i+1)+fnu(j,k,i))
                fnum=.5*(fnu(j,k,i-1)+fnu(j,k,i))
                cdp=(fnup+(1.+cb2)*anutp)*ttp/(sigma*re)
                cdm=(fnum+(1.+cb2)*anutm)*ttm/(sigma*re)
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
                cz(k,i)=cz(k,i)*fact+1.0
                dz(k,i)=dz(k,i)*fact
                fz(k,i)=vist3d(j,k,i)
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
c Update VIST3D and save TURRE in TURSAV
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            chi=turre(j,k,i)/fnu(j,k,i)
            fv1=chi**3/(chi**3+cv1**3)
            vist3d(j,k,i)=fv1*turre(j,k,i)*q(j,k,i,1)
            tursav(j,k,i,1)=turre(j,k,i)
          enddo
        enddo
      enddo
c
      return
      end
