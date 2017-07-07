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
      subroutine lesdiag(myid,jdim,kdim,idim,q,ux,vist3d,vol,si,sj,sk,
     +                   vor,smin,xjb,tursav,xkb,blnum,
     +                   nou,bou,nbuf,ibufdim,nbl,nummem,x,y,z)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Compute diagnostics for LES-type runs.  If a 
c     Smagorinsky-type SGS model is being used, the appropriate eddy
c     viscosity is computed and set here.
c     Output goes to files fort.50x, depending on myid (i.e., if 8
c     processors are used, files fort.501-508 will be written).
c
c     Currently coded: standard Smagorinsky model
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /reyue/ reue,tinf,ivisc(3)
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /sklton/ isklton
      common /lesinfo/ les_model,les_wallscale,cs_smagorinsky,cs_wale,
     .                 cs_vreman
c
      dimension q(jdim,kdim,idim,5)
      dimension ux(jdim-1,kdim-1,idim-1,9)
      dimension vist3d(jdim,kdim,idim)
      dimension sj(jdim,kdim,idim-1,5),
     + sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
      dimension vor(jdim-1,kdim-1,idim-1),smin(jdim-1,kdim-1,idim-1)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c     the method for keeping track of nearest wall location (for
c     getting y+ value associated with a field point) is taken from method
c     used for Baldwin-Barth... it's kind of a memory nightmare:
c     xjb=j index, tursav(2)=i index, xkb=k index, blnum=corresp block number
      dimension xjb(jdim-1,kdim-1,idim-1),xkb(jdim-1,kdim-1,idim-1),
     + tursav(jdim,kdim,idim,nummem),blnum(jdim-1,kdim-1,idim-1)
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
c     les_model controls whether or not a model is used (0=no model)
c     (for les_model=1,2,3, this can also be achieved by setting 
c     corresponding cs=0 - e.g., cs_smagorinsky=0 for les_model=1)
c     les_wallscale controls whether the Delta is scaled by van Driest near walls
c     (Smagorinsky model only)
c
c   Note: (10.**(-iexp) is machine zero)
      xminn=10.**(-iexp+1)
c
      c2b=cbar/tinf
      c2bp=c2b+1.0
      re=reue/xmach
c
      if(isklton .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     Computing LES-type turbulent'',
     +'' viscosity, block='',i5)') nbl
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     WARNING: ivisc=25 still under'',
     .    '' development... use at your own risk!'')')
         if (les_model .eq. 0) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      no subgrid model'')')
         else if (les_model .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      Smag model, Cs='',f10.5)') 
     .      cs_smagorinsky
         else if (les_model .eq. 2) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      WALE model, Cs='',f10.5)')
     .      cs_wale
         else if (les_model .eq. 3) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      Vreman model, Cs='',f10.5)')
     .      cs_vreman
         end if
         if (les_wallscale .eq. 0 .or. les_model .ne. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      no wall damping'')')
         else
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      van Driest wall damping on'')')
         end if
c
         if (cs_smagorinsky .ne. 0. .and. les_model .ne. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     WARNING: cs_smagorinsky has'',
     .      '' been set but not used (les_model .ne. 1)'')')
         end if
         if (cs_wale .ne. 0. .and. les_model .ne. 2) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     WARNING: cs_wale has'',
     .      '' been set but not used (les_model .ne. 2)'')')
         end if
         if (cs_vreman .ne. 0. .and. les_model .ne. 3) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''     WARNING: cs_vreman has'',
     .      '' been set but not used (les_model .ne. 3)'')')
         end if
      end if
c
c   acquire mu_t if employing an SGS model
c
      if (les_model .eq. 0) then
      do i=1,idim-1
        do j=1,jdim-1
          do k=1,kdim-1
            vist3d(j,k,i)=0.
          enddo
        enddo
      enddo
      else
      do i=1,idim-1
        do j=1,jdim-1
          do k=1,kdim-1
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
              ibb=max(ibb,1)
              tt=gamma*q(jbb,kbb,ibb,5)/q(jbb,kbb,ibb,1)
              wnu=c2bp*tt*sqrt(tt)/(c2b+tt)/q(jbb,kbb,ibb,1)
              utau=sqrt(wnu*vor(jbb,kbb,ibb)/(q(jbb,kbb,ibb,1)*re))
              ypls=re*q(jbb,kbb,ibb,1)*utau*ccabs(smin(j,k,i))/wnu
            end if
            deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
            deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
            deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
            if( i2d .eq. 0 ) then
               delta = (deltai*deltaj*deltak)**0.333333
            else
               delta = sqrt(deltaj*deltak)
            end if
            if (les_model .eq. 1) then
c             Standard Smagorinsky model:
c             Scale delta by van Driest function
              if (les_wallscale .eq. 1) then
                vandriest=1.0-exp(-ypls/25.)
                delta=delta*vandriest
              end if
              s11 = ux(j,k,i,1)
              s22 = ux(j,k,i,5)
              s33 = ux(j,k,i,9)
              s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
              s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
              s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +            2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              vist3d(j,k,i)=q(j,k,i,1)*(cs_smagorinsky*delta)**2*
     +            sqrt(2.*xis)*reue/xmach
            else if (les_model .eq. 2) then
c             WALE model (Flow, Turb, & Combust 62:183-200 1999)
              s11 = ux(j,k,i,1)
              s22 = ux(j,k,i,5)
              s33 = ux(j,k,i,9)
              s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
              s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
              s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +            2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              g11=ux(j,k,i,1)*ux(j,k,i,1) + ux(j,k,i,2)*ux(j,k,i,4) +
     +            ux(j,k,i,3)*ux(j,k,i,7)
              g22=ux(j,k,i,4)*ux(j,k,i,2) + ux(j,k,i,5)*ux(j,k,i,5) +
     +            ux(j,k,i,6)*ux(j,k,i,8)
              g33=ux(j,k,i,7)*ux(j,k,i,3) + ux(j,k,i,8)*ux(j,k,i,6) +
     +            ux(j,k,i,9)*ux(j,k,i,9)
              g12=0.5*(ux(j,k,i,1)*ux(j,k,i,2) +
     +                 ux(j,k,i,2)*ux(j,k,i,5) +
     +                 ux(j,k,i,3)*ux(j,k,i,8) +
     +                 ux(j,k,i,4)*ux(j,k,i,1) + 
     +                 ux(j,k,i,5)*ux(j,k,i,4) +
     +                 ux(j,k,i,6)*ux(j,k,i,7))
              g13=0.5*(ux(j,k,i,1)*ux(j,k,i,3) +
     +                 ux(j,k,i,2)*ux(j,k,i,6) +
     +                 ux(j,k,i,3)*ux(j,k,i,9) +
     +                 ux(j,k,i,7)*ux(j,k,i,1) + 
     +                 ux(j,k,i,8)*ux(j,k,i,4) +
     +                 ux(j,k,i,9)*ux(j,k,i,7))
              g23=0.5*(ux(j,k,i,4)*ux(j,k,i,3) +
     +                 ux(j,k,i,5)*ux(j,k,i,6) +
     +                 ux(j,k,i,6)*ux(j,k,i,9) +
     +                 ux(j,k,i,7)*ux(j,k,i,2) + 
     +                 ux(j,k,i,8)*ux(j,k,i,5) +
     +                 ux(j,k,i,9)*ux(j,k,i,8))
              g11=g11-0.3333333*(g11+g22+g33)
              g22=g22-0.3333333*(g11+g22+g33)
              g33=g33-0.3333333*(g11+g22+g33)
              g_gamma=g11*g11 + g22*g22 + g33*g33 +
     +            2.*g12*g12 + 2.*g13*g13 + 2.*g23*g23
              denom=xis**2.5 + g_gamma**1.25
              denom=ccmax(denom,xminn)
              vist3d(j,k,i)=q(j,k,i,1)*(cs_wale*delta)**2*
     +            (g_gamma**1.5)/denom*reue/xmach
            else if (les_model .eq. 2) then
c             Vreman model (Phys Fluids 16(10):3670-3681 2004)
              b11=delta**2*(ux(j,k,i,1)*ux(j,k,i,1) +
     +                      ux(j,k,i,4)*ux(j,k,i,4) +
     +                      ux(j,k,i,7)*ux(j,k,i,7))
              b22=delta**2*(ux(j,k,i,2)*ux(j,k,i,2) +
     +                      ux(j,k,i,5)*ux(j,k,i,5) +
     +                      ux(j,k,i,8)*ux(j,k,i,8))
              b33=delta**2*(ux(j,k,i,3)*ux(j,k,i,3) +
     +                      ux(j,k,i,6)*ux(j,k,i,6) +
     +                      ux(j,k,i,9)*ux(j,k,i,9))
              b12=delta**2*(ux(j,k,i,1)*ux(j,k,i,2) +
     +                      ux(j,k,i,4)*ux(j,k,i,5) +
     +                      ux(j,k,i,7)*ux(j,k,i,8))
              b13=delta**2*(ux(j,k,i,1)*ux(j,k,i,3) +
     +                      ux(j,k,i,4)*ux(j,k,i,6) +
     +                      ux(j,k,i,7)*ux(j,k,i,9))
              b23=delta**2*(ux(j,k,i,2)*ux(j,k,i,3) +
     +                      ux(j,k,i,5)*ux(j,k,i,6) +
     +                      ux(j,k,i,8)*ux(j,k,i,9))
              b_beta=b11*b22-b12*b12+b11*b33-b13*b13+b22*b33-b23*b23
              denom=ux(j,k,i,1)**2 + ux(j,k,i,5)**2 + ux(j,k,i,9)**2 +
     +            2.*ux(j,k,i,2)*ux(j,k,i,4) +
     +            2.*ux(j,k,i,3)*ux(j,k,i,7) +
     +            2.*ux(j,k,i,6)*ux(j,k,i,8)
              denom=ccmax(denom,xminn)
              vist3d(j,k,i)=q(j,k,i,1)*cs_vreman*sqrt(b_beta/denom)*
     +            reue/xmach
            end if
          enddo
        enddo
      enddo
      end if
c
c     if(icyc.eq.ncyc1(1) .or. icyc.eq.ncyc1(2) .or. icyc.eq.ncyc1(3)
c    +.or. icyc.eq.ncyc1(4) .or. icyc.eq.ncyc1(5)) then
c       write(51,*) idim-1,jdim-1,kdim-1
c       write(51,*) (((x(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +              (((y(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +              (((z(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1)
c       write(52,*) idim-1,jdim-1,kdim-1
c       time=0.
c       write(52,*) xmach,alpha,reue,time
c       write(52,*) (((xjb(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +              (((xkb(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +         (((tursav(j,k,i,2),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +         (((tursav(j,k,i,1),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +         (((vist3d(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1)
c     end if
c
c   obtain LES diagnostic info
c
      if (icyc .eq. ncyc) then
      c2b=cbar/tinf
      c2bp=c2b+1.0
      eddytime=1./xmach
      sumvel2=0.
      sumux2=0.
      sumvy2=0.
      sumwz2=0.
      sumux3=0.
      sumvy3=0.
      sumwz3=0.
      sumu2=0.
      sumv2=0.
      sumw2=0.
      sumuv=0.
      sumuw=0.
      sumvw=0.
      sumenergy=0.
      fnu=0.
      sumxis=0.
      sumwis=0.
      rho=0.
      n=0
      do i=1,idim-1
      do j=1,jdim-1
      do k=1,kdim-1
        n=n+1
        sumvel2=sumvel2+q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)**2
        sumux2=sumux2+(ux(j,k,i,1)**2)
        sumvy2=sumvy2+(ux(j,k,i,5)**2)
        sumwz2=sumwz2+(ux(j,k,i,9)**2)
        sumux3=sumux3+(ux(j,k,i,1)**3)
        sumvy3=sumvy3+(ux(j,k,i,5)**3)
        sumwz3=sumwz3+(ux(j,k,i,9)**3)
        sumu2=sumu2+q(j,k,i,2)**2
        sumv2=sumv2+q(j,k,i,3)**2
        sumw2=sumw2+q(j,k,i,4)**2
        sumuv=sumuv+q(j,k,i,2)*q(j,k,i,3)
        sumuw=sumuw+q(j,k,i,2)*q(j,k,i,4)
        sumvw=sumvw+q(j,k,i,3)*q(j,k,i,4)
        energy=q(j,k,i,5)/gm1 + 0.5*q(j,k,i,1)*(q(j,k,i,2)**2 +
     +         q(j,k,i,3)**2 + q(j,k,i,4)**2)
        sumenergy=sumenergy+energy
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
     +        2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
        wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
        tt=gamma*q(j,k,i,5)/q(j,k,i,1)
        fnu=fnu+(c2bp*tt*sqrt(tt)/(c2b+tt))
        sumxis=sumxis+xis
        sumwis=sumwis+wis
        rho=rho+q(j,k,i,1)
      enddo
      enddo
      enddo
      avgux2=sumux2/float(n)
      avgvy2=sumvy2/float(n)
      avgwz2=sumwz2/float(n)
      avgux3=sumux3/float(n)
      avgvy3=sumvy3/float(n)
      avgwz3=sumwz3/float(n)
      avgvel2=sumvel2/float(n)
      avgu2=sumu2/float(n)
      avgv2=sumv2/float(n)
      avgw2=sumw2/float(n)
      avguv=sumuv/float(n)
      avguw=sumuw/float(n)
      avgvw=sumvw/float(n)
      avgfnu=fnu/float(n)
      avgrho=rho/float(n)
      avgwis=sumwis/float(n)
      xtime=time/eddytime
c   for each processor, write out:
c   ntt,icyc,xtime,sumvel2,sumxis,avgux2,avgvy2,avgwz2,
c   avgux3,avgvy3,avgwz3,avgvel2,avgu2,avgv2,avgw2,avgfnu,avgrho,avgwis
      write(500+myid,'(2i5,20e15.5)') ntt,icyc,xtime,sumvel2,sumxis,
     +   avgux2,avgvy2,avgwz2,avgux3,avgvy3,avgwz3,avgvel2,
     +   avgu2,avgv2,avgw2,avgfnu,avgrho,avgwis,sumenergy,
     +   avguv,avguw,avgvw
c
      end if
      return
      end
