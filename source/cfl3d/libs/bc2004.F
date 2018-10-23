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
      subroutine bc2004(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,smin,vist3d,vj0,vk0,vi0,mdim,ndim,bcdata,
     .                  filname,iuns,irelv,snj0,snk0,sni0,ntime,
     .                  snjm,snkm,snim,nou,bou,nbuf,ibufdim,myid,nummem,
     .                  iuse3)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set solid wall BC (viscous wall)
c               This same routine is used for 2004, 2014, 2024, and 2034
c               With 2 pieces of auxiliary data:
c
c           1) Set Twtype
c                  Twtype > 0 ... fixed wall temperature Tw/Tinf = Twtype
c                  Twtype = 0 ... adiabatic wall
c                  Twtype < 0 ... fixed wall temperature = stagnation temp
c
c           2) Set suction/blowing boundary conditions via mass-flow 
c              coefficient
c
c                                                  (rho * u)
c                  mass-flow coefficient c_q = ------------------ 
c                                              (rho * u)_infinity
c
c                  c_q = 0...standard solid wall B.C, no flow thru wall
c                  c_q < 0...suction (mass flow OUT of the zone)
c                  c_q > 0...blowing (mass flow INTO the zone)
c
c           When wall functions are employed (iwf=1, resulting when ivisc
c           parameters are negative on input), code uses PAB3D type
c           (see NASA TP 3480).  ONLY the boundary condition values
c           of the turb viscosity vist3d (=vj0,vk0,vi0) are modified.  No
c           changes are made to vist3d in the interior.
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension bcdata(mdim,ndim,2,12)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          smin(jdim-1,kdim-1,idim-1),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
      dimension snj0(jdim-1,kdim-1,idim-1),snk0(jdim-1,kdim-1,idim-1),
     .          sni0(jdim-1,kdim-1,idim-1)
      dimension snjm(jdim-1,kdim-1,idim-1),snkm(jdim-1,kdim-1,idim-1),
     .          snim(jdim-1,kdim-1,idim-1)
      dimension a0(7),a1(5),a2(3)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /wallfun/ iwf(3)
      common /reystressmodel/ issglrrw2012,i_sas_rsm,i_yapterm
      common /wallvel_2034/ uub_2034,vvb_2034,wwb_2034,vvb_xaxisrot2034
c  
c..  data statement (for use with wall functions)
      data a0/2.354039, 0.1179840, -4.2899192e-04, 2.0404148e-06,
     1       -5.1775775e-09, 6.2687308e-12, -2.916958e-15/
      data a1/5.777191, 6.8756983e-02, -7.1582745e-06, 1.5594904e-09,
     1        -1.4865778e-13/
      data a2/31.08654, 5.0429072e-02, -2.0072314e-8/
c   *********
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c     this bc makes use of only one plane of data    
c
      ip    = 1
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=2004 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary        viscous wall with T & cq specified       bctype 2004
c******************************************************************************
      if (nface.eq.3) then
c
      do 400 i=ista,iend1
      ii = i-ista+1
c
      do 300 k=ksta,kend1
      kk = k-ksta+1
      cq            = bcdata(kk,ii,ip,2)
      pb            = q(1,k,i,5)
      dpb           = q(2,k,i,5)-q(1,k,i,5)
      pb            = pb - dpb/2.0
      if (real(pb).le.0.0)  pb = q(1,k,i,5)
      c2            = gamma*q(1,k,i,5)/q(1,k,i,1)
      if (real(bcdata(kk,ii,ip,1)) .gt. 0.) then
        c2  = bcdata(kk,ii,ip,1)
      else if (real(bcdata(kk,ii,ip,1)) .lt. 0.) then
        c2=1.e0+gm1*0.5e0*xmach*xmach
      else
        if (iuns.gt.0 .and. irelv.gt.0) then
        xm2 = (q(1,k,i,2)-xtbj(k,i,1,1))**2+
     +        (q(1,k,i,3)-xtbj(k,i,2,1))**2+
     +        (q(1,k,i,4)-xtbj(k,i,3,1))**2
        else
        xm2 = q(1,k,i,2)**2+q(1,k,i,3)**2+q(1,k,i,4)**2
        end if
        xm2 = xm2/c2
        c2  = c2*(1.+0.5*gm1*xm2)
      end if
c
c     surface velocities
c
      if (iuse3 .eq. 1) then
c       for bc2034 (imposing velocity)
        uub = uub_2034
        vvb = vvb_2034+vvb_xaxisrot2034*sj(1,k,i,3)
        wwb = wwb_2034-vvb_xaxisrot2034*sj(1,k,i,2)
      else
        uub = 0.
        vvb = 0.
        wwb = 0.
      end if
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
        if (iuse3 .eq. 1) then
c         for bc2034 (imposing velocity)
          uub = uub_2034+xtbj(k,i,1,1)
          vvb = vvb_2034+vvb_xaxisrot2034*sj(1,k,i,3)+xtbj(k,i,2,1)
          wwb = wwb_2034-vvb_xaxisrot2034*sj(1,k,i,2)+xtbj(k,i,3,1)
        else
          uub = xtbj(k,i,1,1)
          vvb = xtbj(k,i,2,1)
          wwb = xtbj(k,i,3,1)
        end if
      end if
c
c     add velocity components due to mass flow
c
      uub = uub + xmach*cq*sj(1,k,i,1)*c2/(gamma*pb)
      vvb = vvb + xmach*cq*sj(1,k,i,2)*c2/(gamma*pb)
      wwb = wwb + xmach*cq*sj(1,k,i,3)*c2/(gamma*pb)
c
      qj0(k,i,1,1) = gamma*pb/c2
      qj0(k,i,2,1) = uub
      qj0(k,i,3,1) = vvb
      qj0(k,i,4,1) = wwb
      qj0(k,i,5,1) = pb
      bcj(k,i,1)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      j2 = min(2,jdim1)
      if (j2.eq.1) f23 = 0.0
c
      z1  =   2.0 +1.5*f23
      z2  =       -0.5*f23   
      z3  = -(2.0 +    f23)
c
      qj0(k,i,1,2) = z1*q(1,k,i,1) + z2*q(j2,k,i,1) + z3*qj0(k,i,1,1)
      qj0(k,i,2,2) = z1*q(1,k,i,2) + z2*q(j2,k,i,2) + z3*qj0(k,i,2,1)
      qj0(k,i,3,2) = z1*q(1,k,i,3) + z2*q(j2,k,i,3) + z3*qj0(k,i,3,1)
      qj0(k,i,4,2) = z1*q(1,k,i,4) + z2*q(j2,k,i,4) + z3*qj0(k,i,4,1)
      qj0(k,i,5,2) = z1*q(1,k,i,5) + z2*q(j2,k,i,5) + z3*qj0(k,i,5,1)
  300 continue
  400 continue
c
      if (ivmx.ge.2) then
      if (level .ge. lglobal .and. ntime .ne. 0) then
        if (iwf(2) .eq. 0) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
c   Store mu_t at wall face center, rather than ghost cell for solid wall BC
          vj0(k,i,1,1) = 0.
          vj0(k,i,1,2) = 0.
  191   continue
        else
c   ***Wall Function begin
        do 1191 i=ista,iend1
        do 1191 k=ksta,kend1
          j=1
          c2b=cbar/tinf
          c2bp=c2b+1.0
          uu  = sqrt((q(j,k,i,2)-qj0(k,i,2,1))**2 +
     +               (q(j,k,i,3)-qj0(k,i,3,1))**2 +
     +               (q(j,k,i,4)-qj0(k,i,4,1))**2 )
          if(ivmx.eq.2) then
            dist=snj0(j,k,i)
          else
            dist=ccabs(smin(j,k,i))
          end if
          tt=gamma*qj0(k,i,5,1)/qj0(k,i,1,1)
          fnuw=c2bp*tt*sqrt(tt)/(c2b+tt)
          rc=q(j,k,i,1)*uu*dist/fnuw*reue/xmach
          if (real(rc) .le. 20.24) then
            xnplus = sqrt(rc)
          else if (real(rc) .le. 435.) then
            xnplus = a0(1)+a0(2)*rc+(a0(3)*rc)*rc+(a0(4)*rc)*rc*rc+
     1           (a0(5)*rc)*rc*rc*rc+(a0(6)*rc)*rc*rc*rc*rc+
     2           (a0(7)*rc)*rc*rc*rc*rc*rc
          else if (real(rc) .le. 4000.) then
            xnplus = a1(1)+a1(2)*rc+(a1(3)*rc)*rc+(a1(4)*rc)*rc*rc+
     1           (a1(5)*rc)*rc*rc*rc
          else if (real(rc) .gt. 4000.) then
            xnplus = a2(1)+a2(2)*rc+a2(3)*rc*rc
          end if
          xnplussav=xnplus
          if(real(xnplus) .gt. 10.) then
c      Newton iteration to solve for nplus, assuming it is in log region:
          do num=1,10
            f=rc/xnplus - 2.44*(log(xnplus)) - 5.2
            dfdn=-rc/(xnplus*xnplus) - 2.44/xnplus
            delta=-f/dfdn
            xnplus=ccabs(xnplus+delta)
            if(abs(real(delta)) .lt. 1.e-3) goto 1190
          enddo
c         revert back to approx series soln if Newton iteration fails
          xnplus=xnplussav
 1190     continue
          end if
          dudy = uu/dist
          xmut = fnuw*(xnplus*xnplus*fnuw/(dist**2)/
     +           (q(j,k,i,1)*dudy)*xmach/reue-1.)
          vj0(k,i,1,1) = xmut
          vj0(k,i,1,2) = 0.
 1191   continue
        end if
c   ***Wall Function end
      end if
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal .and. ntime .ne. 0) then
      if (ivmx.eq.4 .or. ivmx.eq.5) then
        do 101 i=ista,iend1
        do 101 k=ksta,kend1
          tj0(k,i,1,1) = -tursav(1,k,i,1)
          tj0(k,i,2,1) = -tursav(1,k,i,2)
          tj0(k,i,1,2) = 2.*tj0(k,i,1,1)-tursav(1,k,i,1)
          tj0(k,i,2,2) = 2.*tj0(k,i,2,1)-tursav(1,k,i,2)
  101   continue
      end if
      if (ivmx.ge.6) then
        c2b=cbar/tinf
        c2bp=c2b+1.0
        re=reue/xmach
        if (     ivmx.eq.6) then
          beta1=.075
        else if (ivmx.eq.7 .or. ivmx.eq.30 .or. ivmx.eq.40) then
          beta1=.075
        else if (ivmx.eq.8 .or. ivmx.eq.12 .or. ivmx.eq.14) then
          beta1=.83
        else if (ivmx.eq.72) then
          beta1=.0708
        end if
        j=1
c       (Note:  linearly interpolates eps & omega BC to ghost cell.  Previous
c       versions of k-omega simply put 60*mu/(rho*beta*smin**2) at ghost cell 
c       rather than at wall, where it "officially" belongs.  However, this ends
c       up making very little, if any, difference in the k-omega solution.)
        if (ivmx.eq. 9 .or. ivmx.eq.10 .or. ivmx.eq.11 .or.
     .      ivmx.eq.13) then
        do 902 i=ista,iend1
        do 902 k=ksta,kend1
          tt=gamma*qj0(k,i,5,1)/qj0(k,i,1,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=sqrt(tursav(j,k,i,2))/ccabs(smin(j,k,i))
          tj0(k,i,1,1) = 2.*(2.*fnu/(q(j,k,i,1)*re**2)*dkdy**2)-
     +                   tursav(j,k,i,1)
          tj0(k,i,2,1) = -tursav(j,k,i,2)
          tj0(k,i,1,2) = 2.*tj0(k,i,1,1)-tursav(1,k,i,1)
          tj0(k,i,2,2) = 2.*tj0(k,i,2,1)-tursav(1,k,i,2)
 902    continue
        else if(ivmx.eq.15) then
        do i=ista,iend1
        do k=ksta,kend1
          tt=gamma*q(j,k,i,5)/q(j,k,i,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=0.5*(tursav(j,k,i,2)+tursav(j+1,k,i,2))/3.
          tj0(k,i,1,1) = 2.0*dkdy/(re*smin(j,k,i))**2 - tursav(j,k,i,1)
          tj0(k,i,2,1) = -tursav(j,k,i,2)
          tj0(k,i,1,2) = 0.
          tj0(k,i,2,2) = 0.
        enddo
        enddo
        else if(ivmx.eq.16) then
        do i=ista,iend1
        do k=ksta,kend1
          tj0(k,i,1,1) = -tursav(j,k,i,1)
          tj0(k,i,2,1) = -tursav(j,k,i,2)
          tj0(k,i,1,2) = 0.
          tj0(k,i,2,2) = 0.
        enddo
        enddo
#   ifdef CMPLX
c         RSM not working in complex mode
#   else
        else if(ivmx.eq.70) then  ! k-zeta Reynolds Stress model - xxiao
           do i=ista,iend1
              do k=ksta,kend1
                 do iv=1,6
                    tj0(k,i,iv,1) = -tursav(j,k,i,iv)
                    tj0(k,i,iv,2) = huge(tke1)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 tke1 = -0.5*sum(tursav(j,k,i,1:3))
                 tke2 = -0.5*sum(tursav(j+1,k,i,1:3))
                 dkdy=0.5*(tke1+tke2)/3.*0.5
                 tj0(k,i,7,1) = 2.0*dkdy/(re*smin(j,k,i))**2 
     $                - tursav(j,k,i,7)
                 tj0(k,i,7,2) = 0.0
              enddo
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.ne.6) then
           do i=ista,iend1
              do k=ksta,kend1
                 do iv=1,6
                    tj0(k,i,iv,1) = -tursav(j,k,i,iv)
                    tj0(k,i,iv,2) = 2*tj0(k,i,iv,1) - tursav(j,k,i,iv)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 dist=ccabs(smin(j,k,i))
                 tj0(k,i,7,1) = 
     $                   2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +                   tursav(j,k,i,7)
                 tj0(k,i,7,2) = 2*tj0(k,i,7,1) - tursav(j,k,i,7) 
              enddo
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.eq.6) then
           do i=ista,iend1
              do k=ksta,kend1
                 do iv=1,7
                    tj0(k,i,iv,1) = -tursav(j,k,i,iv)
                    tj0(k,i,iv,2) = 2*tj0(k,i,iv,1) - tursav(j,k,i,iv)
                 enddo
              enddo
           enddo
#   endif
        else if(ivmx.eq.25) then
c       not used, but set anyway
        do i=ista,iend1
        do k=ksta,kend1
          tj0(k,i,1,1) = 0.
          tj0(k,i,2,1) = 0.
          tj0(k,i,1,2) = 0.
          tj0(k,i,2,2) = 0.
        enddo
        enddo
        else
c       all other models (ivmx=6,7,8,12,14)
        do 102 i=ista,iend1
        do 102 k=ksta,kend1
          tt=gamma*qj0(k,i,5,1)/qj0(k,i,1,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dist=ccabs(smin(j,k,i))
          tj0(k,i,1,1) = 2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +                   tursav(j,k,i,1)
          tj0(k,i,2,1) = -tursav(j,k,i,2)
          tj0(k,i,1,2) = 2.*tj0(k,i,1,1)-tursav(1,k,i,1)
          tj0(k,i,2,2) = 2.*tj0(k,i,2,1)-tursav(1,k,i,2)
 102    continue
        end if
c       for trans model, intermittency at wall is either zero flux or set by user
        if (ivmx .eq. 30) then
        do i=ista,iend1
        ii = i-ista+1
        do k=ksta,kend1
        kk = k-ksta+1
          if (iuse3 .ne. 2) then
c           bc2004 or 2014 -> zero flux:
            tj0(k,i,3,1) = tursav(j,k,i,3)
          else
c           tj0(k,i,3,1) = 2.*bcdata(kk,ii,ip,3)-tursav(j,k,i,3)
c           bc2024 -> d(gamma)/dn=-(bcdata(3))*rho/mu*Uedge:
            tt=gamma*q(j,k,i,5)/q(j,k,i,1)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            hee=(1.+0.5*gm1*xmach**2)/gm1
            rhoee=(gamma*q(j,k,i,5))**(1./gamma)
            uee2=(gm1*rhoee*hee-gamma*q(j,k,i,5))/(0.5*gm1*rhoee)
            if (real(uee2) .gt. 0.) then
              uee=sqrt(uee2)
            else
              uee=xmach
            end if
            tj0(k,i,3,1) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
     +        uee*bcdata(kk,ii,ip,3)*reue/xmach
          end if
          tj0(k,i,3,2) = 2.*tj0(k,i,3,1)-tursav(1,k,i,3)
        enddo
        enddo
        else if (ivmx .eq. 40) then
        do i=ista,iend1
        ii = i-ista+1
        do k=ksta,kend1
        kk = k-ksta+1
          tj0(k,i,3,1) = tursav(j,k,i,3)
          tj0(k,i,4,1) = tursav(j,k,i,4)
        enddo
        enddo
        end if
      end if
      end if
c
      end if
c
c******************************************************************************
c      j=jdim boundary     viscous wall with T & cq specified       bctype 2004
c******************************************************************************
      if (nface.eq.4) then
c
      do 800 i=ista,iend1
      ii = i-ista+1
c
      do 700 k=ksta,kend1
      kk = k-ksta+1
      cq            = bcdata(kk,ii,ip,2)
      pb            = q(jdim1,k,i,5)
      dpb           = q(jdim1,k,i,5)-q(jdim1-1,k,i,5)
      pb            = pb + dpb/2.0
      if (real(pb).le.0.0)  pb = q(jdim1,k,i,5)
      c2            = gamma*q(jdim1,k,i,5)/q(jdim1,k,i,1)
      if (real(bcdata(kk,ii,ip,1)) .gt. 0.) then
        c2  = bcdata(kk,ii,ip,1)
      else if (real(bcdata(kk,ii,ip,1)) .lt. 0.) then
        c2=1.e0+gm1*0.5e0*xmach*xmach
      else
        if (iuns.gt.0 .and. irelv.gt.0) then
        xm2 = (q(jdim1,k,i,2)-xtbj(k,i,1,2))**2+
     +        (q(jdim1,k,i,3)-xtbj(k,i,2,2))**2+
     +        (q(jdim1,k,i,4)-xtbj(k,i,3,2))**2
        else
        xm2 = q(jdim1,k,i,2)**2+q(jdim1,k,i,3)**2+q(jdim1,k,i,4)**2
        end if
        xm2 = xm2/c2
        c2  = c2*(1.+0.5*gm1*xm2)
      end if
c
c     surface velocities
c
      if (iuse3 .eq. 1) then
c       for bc2034 (imposing velocity)
        uub = uub_2034
        vvb = vvb_2034+vvb_xaxisrot2034*sj(jdim,k,i,3)
        wwb = wwb_2034-vvb_xaxisrot2034*sj(jdim,k,i,2)
      else
        uub = 0.
        vvb = 0.
        wwb = 0.
      end if
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
        if (iuse3 .eq. 1) then
c         for bc2034 (imposing velocity)
          uub = uub_2034+xtbj(k,i,1,2)
          vvb = vvb_2034+vvb_xaxisrot2034*sj(jdim,k,i,3)+xtbj(k,i,2,2)
          wwb = wwb_2034-vvb_xaxisrot2034*sj(jdim,k,i,2)+xtbj(k,i,3,2)
        else
          uub = xtbj(k,i,1,2)
          vvb = xtbj(k,i,2,2)
          wwb = xtbj(k,i,3,2)
        end if
      end if
c
c     add velocity components due to mass flow
c
      uub = uub - xmach*cq*sj(jdim,k,i,1)*c2/(gamma*pb)
      vvb = vvb - xmach*cq*sj(jdim,k,i,2)*c2/(gamma*pb)
      wwb = wwb - xmach*cq*sj(jdim,k,i,3)*c2/(gamma*pb)
c
      qj0(k,i,1,3) = gamma*pb/c2
      qj0(k,i,2,3) = uub
      qj0(k,i,3,3) = vvb
      qj0(k,i,4,3) = wwb
      qj0(k,i,5,3) = pb
      bcj(k,i,2)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      j2 = max(1,jdim-2)
      if (j2.eq.1) f23 = 0.0
c
      z1  =  -2.0 -1.5*f23
      z2  =       +0.5*f23   
      z3  = +(2.0 +    f23)
c
      qj0(k,i,1,4) = z1*q(jdim1,k,i,1)+z2*q(j2,k,i,1)+z3*qj0(k,i,1,3)
      qj0(k,i,2,4) = z1*q(jdim1,k,i,2)+z2*q(j2,k,i,2)+z3*qj0(k,i,2,3)
      qj0(k,i,3,4) = z1*q(jdim1,k,i,3)+z2*q(j2,k,i,3)+z3*qj0(k,i,3,3)
      qj0(k,i,4,4) = z1*q(jdim1,k,i,4)+z2*q(j2,k,i,4)+z3*qj0(k,i,4,3)
      qj0(k,i,5,4) = z1*q(jdim1,k,i,5)+z2*q(j2,k,i,5)+z3*qj0(k,i,5,3)
  700 continue
  800 continue
c
      if (ivmx.ge.2) then
      if (level .ge. lglobal .and. ntime .ne. 0) then
        if (iwf(2) .eq. 0) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
c   Store mu_t at wall face center, rather than ghost cell for solid wall BC
          vj0(k,i,1,3) = 0.
          vj0(k,i,1,4) = 0.
  291   continue
        else
c   ***Wall Function begin
        do 1291 i=ista,iend1
        do 1291 k=ksta,kend1
          j=jdim-1
          c2b=cbar/tinf
          c2bp=c2b+1.0
          uu  = sqrt((q(j,k,i,2)-qj0(k,i,2,3))**2 +
     +               (q(j,k,i,3)-qj0(k,i,3,3))**2 +
     +               (q(j,k,i,4)-qj0(k,i,4,3))**2 )
          if(ivmx.eq.2) then
            dist=snjm(j,k,i)
          else
            dist=ccabs(smin(j,k,i))
          end if
          tt=gamma*qj0(k,i,5,3)/qj0(k,i,1,3)
          fnuw=c2bp*tt*sqrt(tt)/(c2b+tt)
          rc=q(j,k,i,1)*uu*dist/fnuw*reue/xmach
          if (real(rc) .le. 20.24) then
            xnplus = sqrt(rc)
          else if (real(rc) .le. 435.) then
            xnplus = a0(1)+a0(2)*rc+(a0(3)*rc)*rc+(a0(4)*rc)*rc*rc+
     1           (a0(5)*rc)*rc*rc*rc+(a0(6)*rc)*rc*rc*rc*rc+
     2           (a0(7)*rc)*rc*rc*rc*rc*rc
          else if (real(rc) .le. 4000.) then
            xnplus = a1(1)+a1(2)*rc+(a1(3)*rc)*rc+(a1(4)*rc)*rc*rc+
     1           (a1(5)*rc)*rc*rc*rc
          else if (real(rc) .gt. 4000.) then
            xnplus = a2(1)+a2(2)*rc+a2(3)*rc*rc
          end if
          xnplussav=xnplus
          if(real(xnplus) .gt. 10.) then
c      Newton iteration to solve for nplus, assuming it is in log region:
          do num=1,10
            f=rc/xnplus - 2.44*(log(xnplus)) - 5.2
            dfdn=-rc/(xnplus*xnplus) - 2.44/xnplus
            delta=-f/dfdn
            xnplus=ccabs(xnplus+delta)
            if(abs(real(delta)) .lt. 1.e-3) goto 1290
          enddo
c         revert back to approx series soln if Newton iteration fails
          xnplus=xnplussav
 1290     continue
          end if
          dudy = uu/dist
          xmut = fnuw*(xnplus*xnplus*fnuw/(dist**2)/
     +           (q(j,k,i,1)*dudy)*xmach/reue-1.)
          vj0(k,i,1,3) = xmut
          vj0(k,i,1,4) = 0.
 1291   continue
        end if
c   ***Wall Function end
      end if
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal .and. ntime .ne. 0) then
      if (ivmx.eq.4 .or. ivmx.eq.5) then
        do 201 i=ista,iend1
        do 201 k=ksta,kend1
          tj0(k,i,1,3) = -tursav(jdim-1,k,i,1)
          tj0(k,i,2,3) = -tursav(jdim-1,k,i,2)
          tj0(k,i,1,4) = 2.*tj0(k,i,1,3)-tursav(jdim-1,k,i,1)
          tj0(k,i,2,4) = 2.*tj0(k,i,2,3)-tursav(jdim-1,k,i,2)
  201   continue
      end if
      if (ivmx.ge.6) then
        c2b=cbar/tinf
        c2bp=c2b+1.0
        re=reue/xmach
        if (     ivmx.eq.6) then
          beta1=.075
        else if (ivmx.eq.7 .or. ivmx.eq.30 .or. ivmx.eq.40) then
          beta1=.075
        else if (ivmx.eq.8 .or. ivmx.eq.12 .or. ivmx.eq.14) then
          beta1=.83
        else if (ivmx.eq.72) then
          beta1=.0708
        end if
        j=jdim-1
c       (Note:  linearly interpolates eps & omega BC to ghost cell.  Previous
c       versions of k-omega simply put 60*mu/(rho*beta*smin**2) at ghost cell 
c       rather than at wall, where it "officially" belongs.  However, this ends
c       up making very little, if any, difference in the k-omega solution.)
        if (ivmx.eq.9 .or. ivmx.eq.10 .or. ivmx.eq.11 .or.
     .      ivmx.eq.13) then
        do 812 i=ista,iend1
        do 812 k=ksta,kend1
          tt=gamma*qj0(k,i,5,3)/qj0(k,i,1,3)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=sqrt(tursav(j,k,i,2))/ccabs(smin(j,k,i))
          tj0(k,i,1,3) = 2.*(2.*fnu/(q(j,k,i,1)*re**2)*dkdy**2)-
     +                   tursav(j,k,i,1)
          tj0(k,i,2,3) = -tursav(j,k,i,2)
          tj0(k,i,1,4) = 2.*tj0(k,i,1,3)-tursav(jdim-1,k,i,1)
          tj0(k,i,2,4) = 2.*tj0(k,i,2,3)-tursav(jdim-1,k,i,2)
 812    continue
        else if(ivmx.eq.15) then 
        do i=ista,iend1
        do k=ksta,kend1
          tt=gamma*q(j,k,i,5)/q(j,k,i,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=0.5*(tursav(j,k,i,2)+tursav(j-1,k,i,2))/3.
          tj0(k,i,1,3) = 2.0*dkdy/(re*smin(j,k,i))**2 - tursav(j,k,i,1)
          tj0(k,i,2,3) = -tursav(j,k,i,2)
          tj0(k,i,1,4) = 0.
          tj0(k,i,2,4) = 0.
        enddo
        enddo
        else if(ivmx.eq.16) then
        do i=ista,iend1
        do k=ksta,kend1
          tj0(k,i,1,3) = -tursav(j,k,i,1)
          tj0(k,i,2,3) = -tursav(j,k,i,2)
          tj0(k,i,1,4) = 0.
          tj0(k,i,2,4) = 0.
        enddo
        enddo
#   ifdef CMPLX
c         RSM not working in complex mode
#   else
        else if(ivmx.eq.70) then  ! k-zeta Reynolds Stress model - xxiao
           do i=ista,iend1
              do k=ksta,kend1
                 do iv=1,6
                    tj0(k,i,iv,3) = -tursav(j,k,i,iv)
                    tj0(k,i,iv,4) = huge(tke1)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 tke1 = -0.5*sum(tursav(j,k,i,1:3))
                 tke2 = -0.5*sum(tursav(j-1,k,i,1:3))
                 dkdy=0.5*(tke1+tke2)/3.*0.5
                 tj0(k,i,7,3) = 2.0*dkdy/(re*smin(j,k,i))**2 
     $                - tursav(j,k,i,7)
                 tj0(k,i,7,4) = 0.0
              enddo
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.ne.6) then
           do i=ista,iend1
              do k=ksta,kend1
                 do iv=1,6
                    tj0(k,i,iv,3) = -tursav(j,k,i,iv)
                    tj0(k,i,iv,4) = 2*tj0(k,i,iv,3) - tursav(j,k,i,iv)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 dist=ccabs(smin(j,k,i))
                 tj0(k,i,7,3) = 
     $                  2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +                  tursav(j,k,i,7)
                 tj0(k,i,7,4) = 2*tj0(k,i,7,3) - tursav(j,k,i,7)
              enddo
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.eq.6) then
           do i=ista,iend1
              do k=ksta,kend1
                 do iv=1,7
                    tj0(k,i,iv,3) = -tursav(j,k,i,iv)
                    tj0(k,i,iv,4) = 2*tj0(k,i,iv,3) - tursav(j,k,i,iv)
                 enddo
              enddo
           enddo
#   endif
        else if(ivmx.eq.25) then
c       not used, but set anyway
        do i=ista,iend1
        do k=ksta,kend1
          tj0(k,i,1,3) = 0.
          tj0(k,i,2,3) = 0.
          tj0(k,i,1,4) = 0.
          tj0(k,i,2,4) = 0.
        enddo
        enddo
        else
c       all other models (ivmx=6,7,8,12,14)
        do 202 i=ista,iend1
        do 202 k=ksta,kend1
          tt=gamma*qj0(k,i,5,3)/qj0(k,i,1,3)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dist=ccabs(smin(j,k,i))
          tj0(k,i,1,3) = 2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +                   tursav(j,k,i,1)
          tj0(k,i,2,3) = -tursav(j,k,i,2)
          tj0(k,i,1,4) = 2.*tj0(k,i,1,3)-tursav(jdim-1,k,i,1)
          tj0(k,i,2,4) = 2.*tj0(k,i,2,3)-tursav(jdim-1,k,i,2)
 202    continue
        end if
c       for trans model, intermittency at wall is either zero flux or set by user
        if (ivmx .eq. 30) then
        do i=ista,iend1
        ii = i-ista+1
        do k=ksta,kend1
        kk = k-ksta+1
          if (iuse3 .ne. 2) then
c           bc2004 or 2014 -> zero flux:
            tj0(k,i,3,3) = tursav(j,k,i,3)
          else
c           tj0(k,i,3,3) = 2.*bcdata(kk,ii,ip,3)-tursav(j,k,i,3)
c           bc2024 -> d(gamma)/dn=-(bcdata(3))*rho/mu*Uedge:
            tt=gamma*q(j,k,i,5)/q(j,k,i,1)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            hee=(1.+0.5*gm1*xmach**2)/gm1
            rhoee=(gamma*q(j,k,i,5))**(1./gamma)
            uee2=(gm1*rhoee*hee-gamma*q(j,k,i,5))/(0.5*gm1*rhoee)
            if (real(uee2) .gt. 0.) then
              uee=sqrt(uee2)
            else
              uee=xmach
            end if
            tj0(k,i,3,3) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
     +        uee*bcdata(kk,ii,ip,3)*reue/xmach
          end if
          tj0(k,i,3,4) = 2.*tj0(k,i,3,3)-tursav(jdim-1,k,i,3)
        enddo
        enddo
        else if (ivmx .eq. 40) then
        do i=ista,iend1
        ii = i-ista+1
        do k=ksta,kend1
        kk = k-ksta+1
          tj0(k,i,3,3) = tursav(j,k,i,3)
          tj0(k,i,4,3) = tursav(j,k,i,4)
        enddo
        enddo
        end if
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=1 boundary        viscous wall with T & cq specified       bctype 2004
c******************************************************************************
      if (nface.eq.5) then
c
      do 1200 i=ista,iend1
      ii = i-ista+1
c
      do 1100 j=jsta,jend1
      jj = j-jsta+1
      cq            = bcdata(jj,ii,ip,2) 
      pb            = q(j,1,i,5)
      dpb           = q(j,2,i,5)-q(j,1,i,5)
      pb            = pb - dpb/2.0
      if (real(pb).le.0.0)  pb = q(j,1,i,5)
      c2            = gamma*q(j,1,i,5)/q(j,1,i,1)
      if (real(bcdata(jj,ii,ip,1)) .gt. 0.) then
        c2  = bcdata(jj,ii,ip,1)
      else if (real(bcdata(jj,ii,ip,1)) .lt. 0.) then
        c2=1.e0+gm1*0.5e0*xmach*xmach
      else
        if (iuns.gt.0 .and. irelv.gt.0) then
        xm2 = (q(j,1,i,2)-xtbk(j,i,1,1))**2+
     +        (q(j,1,i,3)-xtbk(j,i,2,1))**2+
     +        (q(j,1,i,4)-xtbk(j,i,3,1))**2
        else
        xm2 = q(j,1,i,2)**2+q(j,1,i,3)**2+q(j,1,i,4)**2
        end if
        xm2 = xm2/c2
        c2  = c2*(1.+0.5*gm1*xm2)
      end if
c
c     surface velocities
c
      if (iuse3 .eq. 1) then
c       for bc2034 (imposing velocity)
        uub = uub_2034
        vvb = vvb_2034+vvb_xaxisrot2034*sk(j,1,i,3)
        wwb = wwb_2034-vvb_xaxisrot2034*sk(j,1,i,2)
      else
        uub = 0.
        vvb = 0.
        wwb = 0.
      end if
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
        if (iuse3 .eq. 1) then
c         for bc2034 (imposing velocity)
          uub = uub_2034+xtbk(j,i,1,1)
          vvb = vvb_2034+vvb_xaxisrot2034*sk(j,1,i,3)+xtbk(j,i,2,1)
          wwb = wwb_2034-vvb_xaxisrot2034*sk(j,1,i,2)+xtbk(j,i,3,1)
        else
          uub = xtbk(j,i,1,1)
          vvb = xtbk(j,i,2,1)
          wwb = xtbk(j,i,3,1)
        end if
      end if
c
c     add velocity components due to mass flow
c
      uub = uub + xmach*cq*sk(j,1,i,1)*c2/(gamma*pb)
      vvb = vvb + xmach*cq*sk(j,1,i,2)*c2/(gamma*pb)
      wwb = wwb + xmach*cq*sk(j,1,i,3)*c2/(gamma*pb)
c
      qk0(j,i,1,1) = gamma*pb/c2
      qk0(j,i,2,1) = uub
      qk0(j,i,3,1) = vvb
      qk0(j,i,4,1) = wwb
      qk0(j,i,5,1) = pb
      bck(j,i,1)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      k2 = min(2,kdim1)
      if (k2.eq.1) f23 = 0.0
c
      z1  =   2.0 +1.5*f23
      z2  =       -0.5*f23   
      z3  = -(2.0 +    f23)
c
      qk0(j,i,1,2) = z1*q(j,1,i,1) + z2*q(j,k2,i,1) + z3*qk0(j,i,1,1)
      qk0(j,i,2,2) = z1*q(j,1,i,2) + z2*q(j,k2,i,2) + z3*qk0(j,i,2,1)
      qk0(j,i,3,2) = z1*q(j,1,i,3) + z2*q(j,k2,i,3) + z3*qk0(j,i,3,1)
      qk0(j,i,4,2) = z1*q(j,1,i,4) + z2*q(j,k2,i,4) + z3*qk0(j,i,4,1)
      qk0(j,i,5,2) = z1*q(j,1,i,5) + z2*q(j,k2,i,5) + z3*qk0(j,i,5,1)
 1100 continue
 1200 continue
c
      if (ivmx.ge.2) then
      if (level .ge. lglobal .and. ntime .ne. 0) then
        if (iwf(3) .eq. 0) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
c   Store mu_t at wall face center, rather than ghost cell for solid wall BC
          vk0(j,i,1,1) = 0.
          vk0(j,i,1,2) = 0.
  391   continue
        else
c   ***Wall Function begin
        do 1391 i=ista,iend1
        do 1391 j=jsta,jend1
          k=1
          c2b=cbar/tinf
          c2bp=c2b+1.0
          uu  = sqrt((q(j,k,i,2)-qk0(j,i,2,1))**2 +
     +               (q(j,k,i,3)-qk0(j,i,3,1))**2 +
     +               (q(j,k,i,4)-qk0(j,i,4,1))**2 )
          if(ivmx.eq.2) then
            dist=snk0(j,k,i)
          else
            dist=ccabs(smin(j,k,i))
          end if
          tt=gamma*qk0(j,i,5,1)/qk0(j,i,1,1)
          fnuw=c2bp*tt*sqrt(tt)/(c2b+tt)
          rc=q(j,k,i,1)*uu*dist/fnuw*reue/xmach
          if (real(rc) .le. 20.24) then
            xnplus = sqrt(rc)
          else if (real(rc) .le. 435.) then
            xnplus = a0(1)+a0(2)*rc+(a0(3)*rc)*rc+(a0(4)*rc)*rc*rc+
     1           (a0(5)*rc)*rc*rc*rc+(a0(6)*rc)*rc*rc*rc*rc+
     2           (a0(7)*rc)*rc*rc*rc*rc*rc
          else if (real(rc) .le. 4000.) then
            xnplus = a1(1)+a1(2)*rc+(a1(3)*rc)*rc+(a1(4)*rc)*rc*rc+
     1           (a1(5)*rc)*rc*rc*rc
          else if (real(rc) .gt. 4000.) then
            xnplus = a2(1)+a2(2)*rc+a2(3)*rc*rc
          end if
          xnplussav=xnplus
          if(real(xnplus) .gt. 10.) then
c      Newton iteration to solve for nplus, assuming it is in log region:
          do num=1,10
            f=rc/xnplus - 2.44*(log(xnplus)) - 5.2
            dfdn=-rc/(xnplus*xnplus) - 2.44/xnplus
            delta=-f/dfdn
            xnplus=ccabs(xnplus+delta)
            if(abs(real(delta)) .lt. 1.e-3) goto 1390
          enddo
c         revert back to approx series soln if Newton iteration fails
          xnplus=xnplussav
 1390     continue
          end if
          dudy = uu/dist
          xmut = fnuw*(xnplus*xnplus*fnuw/(dist**2)/
     +           (q(j,k,i,1)*dudy)*xmach/reue-1.)
          vk0(j,i,1,1) = xmut
          vk0(j,i,1,2) = 0.
 1391   continue
        end if
c   ***Wall Function end
      end if
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal .and. ntime .ne. 0) then
      if (ivmx.eq.4 .or. ivmx.eq.5) then
        do 301 i=ista,iend1
        do 301 j=jsta,jend1
          tk0(j,i,1,1) = -tursav(j,1,i,1)
          tk0(j,i,2,1) = -tursav(j,1,i,2)
          tk0(j,i,1,2) = 2.*tk0(j,i,1,1)-tursav(j,1,i,1)
          tk0(j,i,2,2) = 2.*tk0(j,i,2,1)-tursav(j,1,i,2)
  301   continue
      end if
      if (ivmx.ge.6) then
        c2b=cbar/tinf
        c2bp=c2b+1.0
        re=reue/xmach
        if (     ivmx.eq.6) then
          beta1=.075
        else if (ivmx.eq.7 .or. ivmx.eq.30 .or. ivmx.eq.40) then
          beta1=.075
        else if (ivmx.eq.8 .or. ivmx.eq.12 .or. ivmx.eq.14) then
          beta1=.83
        else if (ivmx.eq.72) then
          beta1=.0708
        end if
        k=1
c       (Note:  linearly interpolates eps & omega BC to ghost cell.  Previous
c       versions of k-omega simply put 60*mu/(rho*beta*smin**2) at ghost cell 
c       rather than at wall, where it "officially" belongs.  However, this ends
c       up making very little, if any, difference in the k-omega solution.)
        if (ivmx.eq.9 .or. ivmx.eq.10 .or. ivmx.eq.11 .or.
     .      ivmx.eq.13) then
        do 732 i=ista,iend1
        do 732 j=jsta,jend1
          tt=gamma*qk0(j,i,5,1)/qk0(j,i,1,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=sqrt(tursav(j,k,i,2))/ccabs(smin(j,k,i))
          tk0(j,i,1,1) = 2.*(2.*fnu/(q(j,k,i,1)*re**2)*dkdy**2)-
     +                   tursav(j,k,i,1)
          tk0(j,i,2,1) = -tursav(j,k,i,2)
          tk0(j,i,1,2) = 2.*tk0(j,i,1,1)-tursav(j,1,i,1)
          tk0(j,i,2,2) = 2.*tk0(j,i,2,1)-tursav(j,1,i,2)
 732    continue
        else if(ivmx.eq.15) then
        do i=ista,iend1
        do j=jsta,jend1
          tt=gamma*q(j,k,i,5)/q(j,k,i,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=0.5*(tursav(j,k,i,2)+tursav(j,k+1,i,2))/3.
          tk0(j,i,1,1) = 2.0*dkdy/(re*smin(j,k,i))**2 - tursav(j,k,i,1)
          tk0(j,i,2,1) = -tursav(j,k,i,2)
          tk0(j,i,1,2) = 0.
          tk0(j,i,2,2) = 0.
        enddo
        enddo
        else if(ivmx.eq.16) then
        do i=ista,iend1
        do j=jsta,jend1
          tk0(j,i,1,1) = -tursav(j,k,i,1)
          tk0(j,i,2,1) = -tursav(j,k,i,2)
          tk0(j,i,1,2) = 0.
          tk0(j,i,2,2) = 0.
        enddo
        enddo
#   ifdef CMPLX
c         RSM not working in complex mode
#   else
        else if(ivmx.eq.70) then  ! k-zeta Reynolds Stress model - xxiao
           do i=ista,iend1
              do j=jsta,jend1
                 do iv=1,6
                    tk0(j,i,iv,1) = -tursav(j,k,i,iv)
                    tk0(j,i,iv,2) = huge(tke1)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 tke1 = -0.5*sum(tursav(j,k,i,1:3))
                 tke2 = -0.5*sum(tursav(j,k+1,i,1:3))
                 dkdy=0.5*(tke1+tke2)/3.*0.5
                 tk0(j,i,7,1) = 2.0*dkdy/(re*smin(j,k,i))**2 
     $                - tursav(j,k,i,7)
                 tk0(j,i,7,2) = 0.0
              enddo
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.ne.6) then
           do i=ista,iend1
              do j=jsta,jend1
                 do iv=1,6
                    tk0(j,i,iv,1) = -tursav(j,k,i,iv)
                    tk0(j,i,iv,2) = 2*tk0(j,i,iv,1) - tursav(j,k,i,iv)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 dist=ccabs(smin(j,k,i))
                 tk0(j,i,7,1) = 
     $               2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +               tursav(j,k,i,7)
                 tk0(j,i,7,2) = 2*tk0(j,i,7,1) - tursav(j,k,i,7)
              enddo 
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.eq.6) then
           do i=ista,iend1
              do j=jsta,jend1
                 do iv=1,7
                    tk0(j,i,iv,1) = -tursav(j,k,i,iv)
                    tk0(j,i,iv,2) = 2*tk0(j,i,iv,1) - tursav(j,k,i,iv)
                 enddo
              enddo
           enddo
#   endif
        else if(ivmx.eq.25) then
c       not used, but set anyway
        do i=ista,iend1
        do j=jsta,jend1
          tk0(j,i,1,1) = 0.
          tk0(j,i,2,1) = 0.
          tk0(j,i,1,2) = 0.
          tk0(j,i,2,2) = 0.
        enddo
        enddo
        else
c       all other models (ivmx=6,7,8,12,14)
        do 302 i=ista,iend1
        do 302 j=jsta,jend1
          tt=gamma*qk0(j,i,5,1)/qk0(j,i,1,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dist=ccabs(smin(j,k,i))
          tk0(j,i,1,1) = 2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +                   tursav(j,k,i,1)
          tk0(j,i,2,1) = -tursav(j,k,i,2)
          tk0(j,i,1,2) = 2.*tk0(j,i,1,1)-tursav(j,1,i,1)
          tk0(j,i,2,2) = 2.*tk0(j,i,2,1)-tursav(j,1,i,2)
 302    continue
        end if
c       for trans model, intermittency at wall is either zero flux or set by user
        if (ivmx .eq. 30) then
        do i=ista,iend1
        ii = i-ista+1
        do j=jsta,jend1
        jj = j-jsta+1
          if (iuse3 .ne. 2) then
c           bc2004 or 2014 -> zero flux:
            tk0(j,i,3,1) = tursav(j,k,i,3)
          else
c           tk0(j,i,3,1) = 2.*bcdata(jj,ii,ip,3)-tursav(j,k,i,3)
c           bc2024 -> d(gamma)/dn=-(bcdata(3))*rho/mu*Uedge:
            tt=gamma*q(j,k,i,5)/q(j,k,i,1)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            hee=(1.+0.5*gm1*xmach**2)/gm1
            rhoee=(gamma*q(j,k,i,5))**(1./gamma)
            uee2=(gm1*rhoee*hee-gamma*q(j,k,i,5))/(0.5*gm1*rhoee)
            if (real(uee2) .gt. 0.) then
              uee=sqrt(uee2)
            else
              uee=xmach
            end if
            tk0(j,i,3,1) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
     +        uee*bcdata(jj,ii,ip,3)*reue/xmach
          end if
          tk0(j,i,3,2) = 2.*tk0(j,i,3,1)-tursav(j,1,i,3)
        enddo
        enddo
        else if (ivmx .eq. 40) then
        do i=ista,iend1
        ii = i-ista+1
        do j=jsta,jend1
        jj = j-jsta+1
          tk0(j,i,3,1) = tursav(j,k,i,3)
          tk0(j,i,4,1) = tursav(j,k,i,4)
        enddo
        enddo
        end if
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=kdim boundary     viscous wall with T & cq specified       bctype 2004
c******************************************************************************
      if (nface.eq.6) then
c 
      do 1600 i=ista,iend1
      ii = i-ista+1
c
      do 1500 j=jsta,jend1
      jj = j-jsta+1
      cq            = bcdata(jj,ii,ip,2) 
      pb            = q(j,kdim1,i,5)
      dpb           = q(j,kdim1,i,5)-q(j,kdim1-1,i,5)
      pb            = pb + dpb/2.0
      if (real(pb).le.0.0)  pb = q(j,kdim1,i,5)
      c2            = gamma*q(j,kdim1,i,5)/q(j,kdim1,i,1)
      if (real(bcdata(jj,ii,ip,1)) .gt. 0.) then
        c2  = bcdata(jj,ii,ip,1)
      else if (real(bcdata(jj,ii,ip,1)) .lt. 0.) then
        c2=1.e0+gm1*0.5e0*xmach*xmach
      else
        if (iuns.gt.0 .and. irelv.gt.0) then
        xm2 = (q(j,kdim1,i,2)-xtbk(j,i,1,2))**2+
     +        (q(j,kdim1,i,3)-xtbk(j,i,2,2))**2+
     +        (q(j,kdim1,i,4)-xtbk(j,i,3,2))**2
        else
        xm2 = q(j,kdim1,i,2)**2+q(j,kdim1,i,3)**2+q(j,kdim1,i,4)**2
        end if
        xm2 = xm2/c2
        c2  = c2*(1.+0.5*gm1*xm2)
      end if
c
c     surface velocities
c
      if (iuse3 .eq. 1) then
c       for bc2034 (imposing velocity)
        uub = uub_2034
        vvb = vvb_2034+vvb_xaxisrot2034*sk(j,kdim,i,3)
        wwb = wwb_2034-vvb_xaxisrot2034*sk(j,kdim,i,2)
      else
        uub = 0.
        vvb = 0.
        wwb = 0.
      end if
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
        if (iuse3 .eq. 1) then
c         for bc2034 (imposing velocity)
          uub = uub_2034+xtbk(j,i,1,2)
          vvb = vvb_2034+vvb_xaxisrot2034*sk(j,kdim,i,3)+xtbk(j,i,2,2)
          wwb = wwb_2034-vvb_xaxisrot2034*sk(j,kdim,i,2)+xtbk(j,i,3,2)
        else
          uub = xtbk(j,i,1,2)
          vvb = xtbk(j,i,2,2)
          wwb = xtbk(j,i,3,2)
        end if
      end if
c
c     add velocity components due to mass flow
c
      uub = uub - xmach*cq*sk(j,kdim,i,1)*c2/(gamma*pb)
      vvb = vvb - xmach*cq*sk(j,kdim,i,2)*c2/(gamma*pb)
      wwb = wwb - xmach*cq*sk(j,kdim,i,3)*c2/(gamma*pb)
c
      qk0(j,i,1,3) = gamma*pb/c2
      qk0(j,i,2,3) = uub
      qk0(j,i,3,3) = vvb
      qk0(j,i,4,3) = wwb
      qk0(j,i,5,3) = pb
      bck(j,i,2)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      k2 = max(1,kdim-2)
      if (k2.eq.1) f23 = 0.0
c
      z1  =  -2.0 -1.5*f23
      z2  =       +0.5*f23   
      z3  = +(2.0 +    f23)
c
      qk0(j,i,1,4) = z1*q(j,kdim1,i,1)+z2*q(j,k2,i,1)+z3*qk0(j,i,1,3)
      qk0(j,i,2,4) = z1*q(j,kdim1,i,2)+z2*q(j,k2,i,2)+z3*qk0(j,i,2,3)
      qk0(j,i,3,4) = z1*q(j,kdim1,i,3)+z2*q(j,k2,i,3)+z3*qk0(j,i,3,3)
      qk0(j,i,4,4) = z1*q(j,kdim1,i,4)+z2*q(j,k2,i,4)+z3*qk0(j,i,4,3)
      qk0(j,i,5,4) = z1*q(j,kdim1,i,5)+z2*q(j,k2,i,5)+z3*qk0(j,i,5,3)
 1500 continue
 1600 continue
c
      if (ivmx.ge.2) then
      if (level .ge. lglobal .and. ntime .ne. 0) then
        if (iwf(3) .eq. 0) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
c   Store mu_t at wall face center, rather than ghost cell for solid wall BC
          vk0(j,i,1,3) = 0.
          vk0(j,i,1,4) = 0.
  491   continue
        else
c   ***Wall Function begin
        do 1491 i=ista,iend1
        do 1491 j=jsta,jend1
          k=kdim-1
          c2b=cbar/tinf
          c2bp=c2b+1.0
          uu  = sqrt((q(j,k,i,2)-qk0(j,i,2,3))**2 +
     +               (q(j,k,i,3)-qk0(j,i,3,3))**2 +
     +               (q(j,k,i,4)-qk0(j,i,4,3))**2 )
          if(ivmx.eq.2) then
            dist=snkm(j,k,i)
          else
            dist=ccabs(smin(j,k,i))
          end if
          tt=gamma*qk0(j,i,5,3)/qk0(j,i,1,3)
          fnuw=c2bp*tt*sqrt(tt)/(c2b+tt)
          rc=q(j,k,i,1)*uu*dist/fnuw*reue/xmach
          if (real(rc) .le. 20.24) then
            xnplus = sqrt(rc)
          else if (real(rc) .le. 435.) then
            xnplus = a0(1)+a0(2)*rc+(a0(3)*rc)*rc+(a0(4)*rc)*rc*rc+
     1           (a0(5)*rc)*rc*rc*rc+(a0(6)*rc)*rc*rc*rc*rc+
     2           (a0(7)*rc)*rc*rc*rc*rc*rc
          else if (real(rc) .le. 4000.) then
            xnplus = a1(1)+a1(2)*rc+(a1(3)*rc)*rc+(a1(4)*rc)*rc*rc+
     1           (a1(5)*rc)*rc*rc*rc
          else if (real(rc) .gt. 4000.) then
            xnplus = a2(1)+a2(2)*rc+a2(3)*rc*rc
          end if
          xnplussav=xnplus
          if(real(xnplus) .gt. 10.) then
c      Newton iteration to solve for nplus, assuming it is in log region:
          do num=1,10
            f=rc/xnplus - 2.44*(log(xnplus)) - 5.2
            dfdn=-rc/(xnplus*xnplus) - 2.44/xnplus
            delta=-f/dfdn
            xnplus=ccabs(xnplus+delta)
            if(abs(real(delta)) .lt. 1.e-3) goto 1490
          enddo
c         revert back to approx series soln if Newton iteration fails
          xnplus=xnplussav
 1490     continue
          end if
          dudy = uu/dist
          xmut = fnuw*(xnplus*xnplus*fnuw/(dist**2)/
     +           (q(j,k,i,1)*dudy)*xmach/reue-1.)
          vk0(j,i,1,3) = xmut
          vk0(j,i,1,4) = 0.
 1491   continue
        end if
c   ***Wall Function end
      end if
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal .and. ntime .ne. 0) then
      if (ivmx.eq.4 .or. ivmx.eq.5) then
        do 401 i=ista,iend1
        do 401 j=jsta,jend1
          tk0(j,i,1,3) = -tursav(j,kdim-1,i,1)
          tk0(j,i,2,3) = -tursav(j,kdim-1,i,2)
          tk0(j,i,1,4) = 2.*tk0(j,i,1,3)-tursav(j,kdim-1,i,1)
          tk0(j,i,2,4) = 2.*tk0(j,i,2,3)-tursav(j,kdim-1,i,2)
  401   continue
      end if
      if (ivmx.ge.6) then
        c2b=cbar/tinf
        c2bp=c2b+1.0
        re=reue/xmach
        if (     ivmx.eq.6) then
          beta1=.075
        else if (ivmx.eq.7 .or. ivmx.eq.30 .or. ivmx.eq.40) then
          beta1=.075
        else if (ivmx.eq.8 .or. ivmx.eq.12 .or. ivmx.eq.14) then
          beta1=.83
        else if (ivmx.eq.72) then
          beta1=.0708
        end if
        k=kdim-1
c       (Note:  linearly interpolates eps & omega BC to ghost cell.  Previous
c       versions of k-omega simply put 60*mu/(rho*beta*smin**2) at ghost cell 
c       rather than at wall, where it "officially" belongs.  However, this ends
c       up making very little, if any, difference in the k-omega solution.)
        if (ivmx.eq.9 .or. ivmx.eq.10 .or. ivmx.eq.11 .or.
     .      ivmx.eq.13) then
        do 492 i=ista,iend1
        do 492 j=jsta,jend1
          tt=gamma*qk0(j,i,5,3)/qk0(j,i,1,3)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=sqrt(tursav(j,k,i,2))/ccabs(smin(j,k,i))
          tk0(j,i,1,3) = 2.*(2.*fnu/(q(j,k,i,1)*re**2)*dkdy**2)-
     +                   tursav(j,k,i,1)
          tk0(j,i,2,3) = -tursav(j,k,i,2)
          tk0(j,i,1,4) = 2.*tk0(j,i,1,3)-tursav(j,kdim-1,i,1)
          tk0(j,i,2,4) = 2.*tk0(j,i,2,3)-tursav(j,kdim-1,i,2)
 492    continue
        else if(ivmx.eq.15) then
        do i=ista,iend1
        do j=jsta,jend1
          tt=gamma*q(j,k,i,5)/q(j,k,i,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=0.5*(tursav(j,k,i,2)+tursav(j,k-1,i,2))/3.
          tk0(j,i,1,3) = 2.0*dkdy/(re*smin(j,k,i))**2 - tursav(j,k,i,1)
          tk0(j,i,2,3) = -tursav(j,k,i,2)
          tk0(j,i,1,4) = 0.
          tk0(j,i,2,4) = 0.
        enddo
        enddo
        else if(ivmx.eq.16) then
        do i=ista,iend1
        do j=jsta,jend1
          tk0(j,i,1,3) = -tursav(j,k,i,1)
          tk0(j,i,2,3) = -tursav(j,k,i,2)
          tk0(j,i,1,4) = 0.
          tk0(j,i,2,4) = 0.
        enddo
        enddo
#   ifdef CMPLX
c         RSM not working in complex mode
#   else
        else if(ivmx.eq.70) then  ! k-zeta Reynolds Stress model - xxiao
           do i=ista,iend1
              do j=jsta,jend1
                 do iv=1,6
                    tk0(j,i,iv,3) = -tursav(j,k,i,iv)
                    tk0(j,i,iv,4) = huge(tke1)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 tke1 = -0.5*sum(tursav(j,k,i,1:3))
                 tke2 = -0.5*sum(tursav(j,k-1,i,1:3))
                 dkdy=0.5*(tke1+tke2)/3.*0.5
                 tk0(j,i,7,3) = 2.0*dkdy/(re*smin(j,k,i))**2 
     $                - tursav(j,k,i,7) 
                 tk0(j,i,7,4) = 0.0
              enddo
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.ne.6) then
           do i=ista,iend1
              do j=jsta,jend1
                 do iv=1,6
                    tk0(j,i,iv,3) = -tursav(j,k,i,iv)
                    tk0(j,i,iv,4) = 2*tk0(j,i,iv,3)- tursav(j,k,i,iv)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 dist=ccabs(smin(j,k,i))
                 tk0(j,i,7,3) = 
     $                 2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +                 tursav(j,k,i,7)
                 tk0(j,i,7,4) = 2*tk0(j,i,7,3)- tursav(j,k,i,7)
              enddo
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.eq.6) then
           do i=ista,iend1
              do j=jsta,jend1
                 do iv=1,7
                    tk0(j,i,iv,3) = -tursav(j,k,i,iv)
                    tk0(j,i,iv,4) = 2*tk0(j,i,iv,3)- tursav(j,k,i,iv)
                 enddo
              enddo
           enddo
#   endif
        else if(ivmx.eq.25) then
c       not used, but set anyway
        do i=ista,iend1
        do j=jsta,jend1
          tk0(j,i,1,3) = 0.
          tk0(j,i,2,3) = 0.
          tk0(j,i,1,4) = 0.
          tk0(j,i,2,4) = 0.
        enddo
        enddo
        else
c       all other models (ivmx=6,7,8,12,14)
        do 402 i=ista,iend1
        do 402 j=jsta,jend1
          tt=gamma*qk0(j,i,5,3)/qk0(j,i,1,3)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dist=ccabs(smin(j,k,i))
          tk0(j,i,1,3) = 2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +                   tursav(j,k,i,1)
          tk0(j,i,2,3) = -tursav(j,k,i,2)
          tk0(j,i,1,4) = 2.*tk0(j,i,1,3)-tursav(j,kdim-1,i,1)
          tk0(j,i,2,4) = 2.*tk0(j,i,2,3)-tursav(j,kdim-1,i,2)
 402    continue
        end if
c       for trans model, intermittency at wall is either zero flux or set by user
        if (ivmx .eq. 30) then
        do i=ista,iend1
        ii = i-ista+1
        do j=jsta,jend1
        jj = j-jsta+1
          if (iuse3 .ne. 2) then
c           bc2004 or 2014 -> zero flux:
            tk0(j,i,3,3) = tursav(j,k,i,3)
          else
c           tk0(j,i,3,3) = 2.*bcdata(jj,ii,ip,3)-tursav(j,k,i,3)
c           bc2024 -> d(gamma)/dn=-(bcdata(3))*rho/mu*Uedge:
            tt=gamma*q(j,k,i,5)/q(j,k,i,1)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            hee=(1.+0.5*gm1*xmach**2)/gm1
            rhoee=(gamma*q(j,k,i,5))**(1./gamma)
            uee2=(gm1*rhoee*hee-gamma*q(j,k,i,5))/(0.5*gm1*rhoee)
            if (real(uee2) .gt. 0.) then
              uee=sqrt(uee2)
            else
              uee=xmach
            end if
            tk0(j,i,3,3) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
     +        uee*bcdata(jj,ii,ip,3)*reue/xmach
          end if
          tk0(j,i,3,4) = 2.*tk0(j,i,3,3)-tursav(j,kdim-1,i,3)
        enddo
        enddo
        else if (ivmx .eq. 40) then
        do i=ista,iend1
        ii = i-ista+1
        do j=jsta,jend1
        jj = j-jsta+1
          tk0(j,i,3,3) = tursav(j,k,i,3)
          tk0(j,i,4,3) = tursav(j,k,i,4)
        enddo
        enddo
        end if
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=1 boundary        viscous wall with T & cq specified       bctype 2004
c******************************************************************************
      if (nface.eq.1) then
c 
      do 2000 k=ksta,kend1
      kk = k-ksta+1
c
      do 1900 j=jsta,jend1
      jj = j-jsta+1
      cq            = bcdata(jj,kk,ip,2) 
      pb            = q(j,k,1,5)
      dpb           = q(j,k,2,5)-q(j,k,1,5)
      pb            = pb - dpb/2.0
      if (real(pb).le.0.0)  pb = q(j,k,1,5)
      c2            = gamma*q(j,k,1,5)/q(j,k,1,1)
      if (real(bcdata(jj,kk,ip,1)) .gt. 0.) then
        c2  = bcdata(jj,kk,ip,1)
      else if (real(bcdata(jj,kk,ip,1)) .lt. 0.) then
        c2=1.e0+gm1*0.5e0*xmach*xmach
      else
        if (iuns.gt.0 .and. irelv.gt.0) then
        xm2 = (q(j,k,1,2)-xtbi(j,k,1,1))**2+
     +        (q(j,k,1,3)-xtbi(j,k,2,1))**2+
     +        (q(j,k,1,4)-xtbi(j,k,3,1))**2
        else
        xm2 = q(j,k,1,2)**2+q(j,k,1,3)**2+q(j,k,1,4)**2
        end if
        xm2 = xm2/c2
        c2  = c2*(1.+0.5*gm1*xm2)
      end if
c
c     surface velocities
c
      if (iuse3 .eq. 1) then
c       for bc2034 (imposing velocity)
        uub = uub_2034
        vvb = vvb_2034+vvb_xaxisrot2034*si(j,k,1,3)
        wwb = wwb_2034-vvb_xaxisrot2034*si(j,k,1,2)
      else
        uub = 0.
        vvb = 0.
        wwb = 0.
      end if
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
        if (iuse3 .eq. 1) then
c         for bc2034 (imposing velocity)
          uub = uub_2034+xtbi(j,k,1,1)
          vvb = vvb_2034+vvb_xaxisrot2034*si(j,k,1,3)+xtbi(j,k,2,1)
          wwb = wwb_2034-vvb_xaxisrot2034*si(j,k,1,2)+xtbi(j,k,3,1)
        else
          uub = xtbi(j,k,1,1)
          vvb = xtbi(j,k,2,1)
          wwb = xtbi(j,k,3,1)
        end if
      end if
c
c     add velocity components due to mass flow
c
      uub = uub + xmach*cq*si(j,k,1,1)*c2/(gamma*pb)
      vvb = vvb + xmach*cq*si(j,k,1,2)*c2/(gamma*pb)
      wwb = wwb + xmach*cq*si(j,k,1,3)*c2/(gamma*pb)
c
      qi0(j,k,1,1) = gamma*pb/c2
      qi0(j,k,2,1) = uub
      qi0(j,k,3,1) = vvb
      qi0(j,k,4,1) = wwb
      qi0(j,k,5,1) = pb
      bci(j,k,1)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      i2 = min(2,idim1)
      if (i2.eq.1) f23 = 0.0
c
      z1 =   2.0 +1.5*f23
      z2 =       -0.5*f23   
      z3 = -(2.0 +    f23)
c
      qi0(j,k,1,2) = z1*q(j,k,1,1) + z2*q(j,k,i2,1) + z3*qi0(j,k,1,1)
      qi0(j,k,2,2) = z1*q(j,k,1,2) + z2*q(j,k,i2,2) + z3*qi0(j,k,2,1)
      qi0(j,k,3,2) = z1*q(j,k,1,3) + z2*q(j,k,i2,3) + z3*qi0(j,k,3,1)
      qi0(j,k,4,2) = z1*q(j,k,1,4) + z2*q(j,k,i2,4) + z3*qi0(j,k,4,1)
      qi0(j,k,5,2) = z1*q(j,k,1,5) + z2*q(j,k,i2,5) + z3*qi0(j,k,5,1)
 1900 continue
 2000 continue
c
      if (ivmx.ge.2) then
      if (level .ge. lglobal .and. ntime .ne. 0) then
        if (iwf(1) .eq. 0) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
c   Store mu_t at wall face center, rather than ghost cell for solid wall BC
          vi0(j,k,1,1) = 0.
          vi0(j,k,1,2) = 0.
  591   continue
        else
c   ***Wall Function begin
        do 1591 k=ksta,kend1
        do 1591 j=jsta,jend1
          i=1
          c2b=cbar/tinf
          c2bp=c2b+1.0
          uu  = sqrt((q(j,k,i,2)-qi0(j,k,2,1))**2 +
     +               (q(j,k,i,3)-qi0(j,k,3,1))**2 +
     +               (q(j,k,i,4)-qi0(j,k,4,1))**2 )
          if(ivmx.eq.2) then
            dist=sni0(j,k,i)
          else
            dist=ccabs(smin(j,k,i))
          end if
          tt=gamma*qi0(j,k,5,1)/qi0(j,k,1,1)
          fnuw=c2bp*tt*sqrt(tt)/(c2b+tt)
          rc=q(j,k,i,1)*uu*dist/fnuw*reue/xmach
          if (real(rc) .le. 20.24) then
            xnplus = sqrt(rc)
          else if (real(rc) .le. 435.) then
            xnplus = a0(1)+a0(2)*rc+(a0(3)*rc)*rc+(a0(4)*rc)*rc*rc+
     1           (a0(5)*rc)*rc*rc*rc+(a0(6)*rc)*rc*rc*rc*rc+
     2           (a0(7)*rc)*rc*rc*rc*rc*rc
          else if (real(rc) .le. 4000.) then
            xnplus = a1(1)+a1(2)*rc+(a1(3)*rc)*rc+(a1(4)*rc)*rc*rc+
     1           (a1(5)*rc)*rc*rc*rc
          else if (real(rc) .gt. 4000.) then
            xnplus = a2(1)+a2(2)*rc+a2(3)*rc*rc
          end if
          xnplussav=xnplus
          if(real(xnplus) .gt. 10.) then
c      Newton iteration to solve for nplus, assuming it is in log region:
          do num=1,10
            f=rc/xnplus - 2.44*(log(xnplus)) - 5.2
            dfdn=-rc/(xnplus*xnplus) - 2.44/xnplus
            delta=-f/dfdn
            xnplus=ccabs(xnplus+delta)
            if(abs(real(delta)) .lt. 1.e-3) goto 1590
          enddo
c         revert back to approx series soln if Newton iteration fails
          xnplus=xnplussav
 1590     continue
          end if
          dudy = uu/dist
          xmut = fnuw*(xnplus*xnplus*fnuw/(dist**2)/
     +           (q(j,k,i,1)*dudy)*xmach/reue-1.)
          vi0(j,k,1,1) = xmut
          vi0(j,k,1,2) = 0.
 1591   continue
        end if
c   ***Wall Function end
      end if
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal .and. ntime .ne. 0) then
      if (ivmx.eq.4 .or. ivmx.eq.5) then
        do 501 k=ksta,kend1
        do 501 j=jsta,jend1
          ti0(j,k,1,1) = -tursav(j,k,1,1)
          ti0(j,k,2,1) = -tursav(j,k,1,2)
          ti0(j,k,1,2) = 2.*ti0(j,k,1,1)-tursav(j,k,1,1)
          ti0(j,k,2,2) = 2.*ti0(j,k,2,1)-tursav(j,k,1,2)
  501   continue
      end if
      if (ivmx.ge.6) then
        c2b=cbar/tinf
        c2bp=c2b+1.0
        re=reue/xmach
        if (     ivmx.eq.6) then
          beta1=.075
        else if (ivmx.eq.7 .or. ivmx.eq.30 .or. ivmx.eq.40) then
          beta1=.075
        else if (ivmx.eq.8 .or. ivmx.eq.12 .or. ivmx.eq.14) then
          beta1=.83
        else if (ivmx.eq.72) then
          beta1=.0708
        end if
        i=1
c       (Note:  linearly interpolates eps & omega BC to ghost cell.  Previous
c       versions of k-omega simply put 60*mu/(rho*beta*smin**2) at ghost cell 
c       rather than at wall, where it "officially" belongs.  However, this ends
c       up making very little, if any, difference in the k-omega solution.)
        if (ivmx.eq.9 .or. ivmx.eq.10 .or. ivmx.eq.11 .or.
     .      ivmx.eq.13) then
        do 652 k=ksta,kend1
        do 652 j=jsta,jend1
          tt=gamma*qi0(j,k,5,1)/qi0(j,k,1,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=sqrt(tursav(j,k,i,2))/ccabs(smin(j,k,i))
          ti0(j,k,1,1) = 2.*(2.*fnu/(q(j,k,i,1)*re**2)*dkdy**2)-
     +                   tursav(j,k,i,1)
          ti0(j,k,2,1) = -tursav(j,k,i,2)
          ti0(j,k,1,2) = 2.*ti0(j,k,1,1)-tursav(j,k,1,1)
          ti0(j,k,2,2) = 2.*ti0(j,k,2,1)-tursav(j,k,1,2)
 652    continue
        else if(ivmx.eq.15) then
        do k=ksta,kend1
        do j=jsta,jend1
          tt=gamma*q(j,k,i,5)/q(j,k,i,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=0.5*(tursav(j,k,i,2)+tursav(j,k,i+1,2))/3.
          ti0(j,k,1,1) = 2.0*dkdy/(re*smin(j,k,i))**2 - tursav(j,k,i,1)
          ti0(j,k,2,1) = -tursav(j,k,i,2)
          ti0(j,k,1,2) = 0.
          ti0(j,k,2,2) = 0.
        enddo
        enddo
        else if(ivmx.eq.16) then
        do k=ksta,kend1
        do j=jsta,jend1
          ti0(j,k,1,1) = -tursav(j,k,i,1)
          ti0(j,k,2,1) = -tursav(j,k,i,2)
          ti0(j,k,1,2) = 0.
          ti0(j,k,2,2) = 0.
        enddo
        enddo
#   ifdef CMPLX
c         RSM not working in complex mode
#   else
        else if(ivmx.eq.70) then  ! k-zeta Reynolds Stress model - xxiao
           do k=ksta,kend1
              do j=jsta,jend1
                 do iv=1,6
                    ti0(j,k,iv,1) = -tursav(j,k,i,iv)
                    ti0(j,k,iv,2) = huge(tke1)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 tke1 = -0.5*sum(tursav(j,k,i,1:3))
                 tke2 = -0.5*sum(tursav(j,k,i+1,1:3))
                 dkdy=0.5*(tke1+tke2)/3.*0.5
                 ti0(j,k,7,1) = 2.0*dkdy/(re*smin(j,k,i))**2 
     $                - tursav(j,k,i,7)
                 ti0(j,k,7,2) = 0.0
              enddo
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.ne.6) then
           do k=ksta,kend1
              do j=jsta,jend1
                 do iv=1,6
                    ti0(j,k,iv,1) = -tursav(j,k,i,iv)
                    ti0(j,k,iv,2) = 2.*ti0(j,k,iv,1)-tursav(j,k,i,iv)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 dist=ccabs(smin(j,k,i))
                 ti0(j,k,7,1) = 
     $               2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +                   tursav(j,k,i,7)
                 ti0(j,k,7,2) = 2.*ti0(j,k,7,1)-tursav(j,k,i,7)
              enddo
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.eq.6) then
           do k=ksta,kend1
              do j=jsta,jend1
                 do iv=1,7
                    ti0(j,k,iv,1) = -tursav(j,k,i,iv)
                    ti0(j,k,iv,2) = 2.*ti0(j,k,iv,1)-tursav(j,k,i,iv)
                 enddo
              enddo
           enddo
#   endif
        else if(ivmx.eq.25) then
c       not used, but set anyway
        do k=ksta,kend1
        do j=jsta,jend1
          ti0(j,k,1,1) = 0.
          ti0(j,k,2,1) = 0.
          ti0(j,k,1,2) = 0.
          ti0(j,k,2,2) = 0.
        enddo
        enddo
        else
c       all other models (ivmx=6,7,8,12,14)
        do 502 k=ksta,kend1
        do 502 j=jsta,jend1
          tt=gamma*qi0(j,k,5,1)/qi0(j,k,1,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dist=ccabs(smin(j,k,i))
          ti0(j,k,1,1) = 2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +                   tursav(j,k,i,1)
          ti0(j,k,2,1) = -tursav(j,k,i,2)
          ti0(j,k,1,2) = 2.*ti0(j,k,1,1)-tursav(j,k,1,1)
          ti0(j,k,2,2) = 2.*ti0(j,k,2,1)-tursav(j,k,1,2)
 502    continue
        end if
c       for trans model, intermittency at wall is either zero flux or set by user
        if (ivmx .eq. 30) then
        do k=ksta,kend1
        kk = k-ksta+1
        do j=jsta,jend1
        jj = j-jsta+1
          if (iuse3 .ne. 2) then
c           bc2004 or 2014 -> zero flux:
            ti0(j,k,3,1) = tursav(j,k,i,3)
          else
c           ti0(j,k,3,1) = 2.*bcdata(jj,kk,ip,3)-tursav(j,k,i,3)
c           bc2024 -> d(gamma)/dn=-(bcdata(3))*rho/mu*Uedge:
            tt=gamma*q(j,k,i,5)/q(j,k,i,1)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            hee=(1.+0.5*gm1*xmach**2)/gm1
            rhoee=(gamma*q(j,k,i,5))**(1./gamma)
            uee2=(gm1*rhoee*hee-gamma*q(j,k,i,5))/(0.5*gm1*rhoee)
            if (real(uee2) .gt. 0.) then
              uee=sqrt(uee2)
            else
              uee=xmach
            end if
            ti0(j,k,3,1) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
     +        uee*bcdata(jj,kk,ip,3)*reue/xmach
          end if
          ti0(j,k,3,2) = 2.*ti0(j,k,3,1)-tursav(j,k,1,3)
        enddo
        enddo
        else if (ivmx .eq. 40) then
        do k=ksta,kend1
        kk = k-ksta+1
        do j=jsta,jend1
        jj = j-jsta+1
          ti0(j,k,3,1) = tursav(j,k,i,3)
          ti0(j,k,4,1) = tursav(j,k,i,4)
        enddo
        enddo
        end if
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=idim boundary     viscous wall with T & cq specified       bctype 2004
c******************************************************************************
      if (nface.eq.2) then
c
      do 2400 k=ksta,kend1
      kk = k-ksta+1
c
      do 2300 j=jsta,jend1
      jj = j-jsta+1
      cq            = bcdata(jj,kk,ip,2)
      pb            = q(j,k,idim1,5)
      dpb           = q(j,k,idim1,5)-q(j,k,idim1-1,5)
      pb            = pb + dpb/2.0
      if (real(pb).le.0.0)  pb = q(j,k,idim1,5)
      c2            = gamma*q(j,k,idim1,5)/q(j,k,idim1,1)
      if (real(bcdata(jj,kk,ip,1)) .gt. 0.) then
        c2  = bcdata(jj,kk,ip,1)
      else if (real(bcdata(jj,kk,ip,1)) .lt. 0.) then
        c2=1.e0+gm1*0.5e0*xmach*xmach
      else
        if (iuns.gt.0 .and. irelv.gt.0) then
        xm2 = (q(j,k,idim1,2)-xtbi(j,k,1,2))**2+
     +        (q(j,k,idim1,3)-xtbi(j,k,2,2))**2+
     +        (q(j,k,idim1,4)-xtbi(j,k,3,2))**2
        else
        xm2 = q(j,k,idim1,2)**2+q(j,k,idim1,3)**2+q(j,k,idim1,4)**2
        end if
        xm2 = xm2/c2
        c2  = c2*(1.+0.5*gm1*xm2)
      end if
c
c     surface velocities
c
      if (iuse3 .eq. 1) then
c       for bc2034 (imposing velocity)
        uub = uub_2034
        vvb = vvb_2034+vvb_xaxisrot2034*si(j,k,idim,3)
        wwb = wwb_2034-vvb_xaxisrot2034*si(j,k,idim,2)
      else
        uub = 0.
        vvb = 0.
        wwb = 0.
      end if
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
        if (iuse3 .eq. 1) then
c         for bc2034 (imposing velocity)
          uub = uub_2034+xtbi(j,k,1,2)
          vvb = vvb_2034+vvb_xaxisrot2034*si(j,k,idim,3)+xtbi(j,k,2,2)
          wwb = wwb_2034-vvb_xaxisrot2034*si(j,k,idim,2)+xtbi(j,k,3,2)
        else
          uub = xtbi(j,k,1,2)
          vvb = xtbi(j,k,2,2)
          wwb = xtbi(j,k,3,2)
        end if
      end if
c
c     add velocity components due to mass flow
c
      uub = uub - xmach*cq*si(j,k,idim,1)*c2/(gamma*pb)
      vvb = vvb - xmach*cq*si(j,k,idim,2)*c2/(gamma*pb)
      wwb = wwb - xmach*cq*si(j,k,idim,3)*c2/(gamma*pb)
c
      qi0(j,k,1,3) = gamma*pb/c2
      qi0(j,k,2,3) = uub
      qi0(j,k,3,3) = vvb
      qi0(j,k,4,3) = wwb
      qi0(j,k,5,3) = pb
      bci(j,k,2)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      i2 = max(1,idim-2)
      if (i2.eq.1) f23 = 0.0
c
      z1 =  -2.0 -1.5*f23
      z2 =       +0.5*f23   
      z3 = +(2.0 +    f23)
c
      qi0(j,k,1,4) = z1*q(j,k,idim1,1)+z2*q(j,k,i2,1)+z3*qi0(j,k,1,3)
      qi0(j,k,2,4) = z1*q(j,k,idim1,2)+z2*q(j,k,i2,2)+z3*qi0(j,k,2,3)
      qi0(j,k,3,4) = z1*q(j,k,idim1,3)+z2*q(j,k,i2,3)+z3*qi0(j,k,3,3)
      qi0(j,k,4,4) = z1*q(j,k,idim1,4)+z2*q(j,k,i2,4)+z3*qi0(j,k,4,3)
      qi0(j,k,5,4) = z1*q(j,k,idim1,5)+z2*q(j,k,i2,5)+z3*qi0(j,k,5,3)
 2300 continue
 2400 continue
c
      if (ivmx.ge.2) then
      if (level .ge. lglobal .and. ntime .ne. 0) then
        if (iwf(1) .eq. 0) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
c   Store mu_t at wall face center, rather than ghost cell for solid wall BC
          vi0(j,k,1,3) = 0.
          vi0(j,k,1,4) = 0.
  691   continue
        else
c   ***Wall Function begin
        do 1691 k=ksta,kend1
        do 1691 j=jsta,jend1
          i=idim-1
          c2b=cbar/tinf
          c2bp=c2b+1.0
          uu  = sqrt((q(j,k,i,2)-qi0(j,k,2,3))**2 +
     +               (q(j,k,i,3)-qi0(j,k,3,3))**2 +
     +               (q(j,k,i,4)-qi0(j,k,4,3))**2 )
          if(ivmx.eq.2) then
            dist=snim(j,k,i)
          else
            dist=ccabs(smin(j,k,i))
          end if
          tt=gamma*qi0(j,k,5,3)/qi0(j,k,1,3)
          fnuw=c2bp*tt*sqrt(tt)/(c2b+tt)
          rc=q(j,k,i,1)*uu*dist/fnuw*reue/xmach
          if (real(rc) .le. 20.24) then
            xnplus = sqrt(rc)
          else if (real(rc) .le. 435.) then
            xnplus = a0(1)+a0(2)*rc+(a0(3)*rc)*rc+(a0(4)*rc)*rc*rc+
     1           (a0(5)*rc)*rc*rc*rc+(a0(6)*rc)*rc*rc*rc*rc+
     2           (a0(7)*rc)*rc*rc*rc*rc*rc
          else if (real(rc) .le. 4000.) then
            xnplus = a1(1)+a1(2)*rc+(a1(3)*rc)*rc+(a1(4)*rc)*rc*rc+
     1           (a1(5)*rc)*rc*rc*rc
          else if (real(rc) .gt. 4000.) then
            xnplus = a2(1)+a2(2)*rc+a2(3)*rc*rc
          end if
          xnplussav=xnplus
          if(real(xnplus) .gt. 10.) then
c      Newton iteration to solve for nplus, assuming it is in log region:
          do num=1,10
            f=rc/xnplus - 2.44*(log(xnplus)) - 5.2
            dfdn=-rc/(xnplus*xnplus) - 2.44/xnplus
            delta=-f/dfdn
            xnplus=ccabs(xnplus+delta)
            if(abs(real(delta)) .lt. 1.e-3) goto 1690
          enddo
c         revert back to approx series soln if Newton iteration fails
          xnplus=xnplussav
 1690     continue
          end if
          dudy = uu/dist
          xmut = fnuw*(xnplus*xnplus*fnuw/(dist**2)/
     +           (q(j,k,i,1)*dudy)*xmach/reue-1.)
          vi0(j,k,1,3) = xmut
          vi0(j,k,1,4) = 0.
 1691   continue
        end if
c   ***Wall Function end
      end if
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal .and. ntime .ne. 0) then
      if (ivmx.eq.4 .or. ivmx.eq.5) then
        do 601 k=ksta,kend1
        do 601 j=jsta,jend1
          ti0(j,k,1,3) = -tursav(j,k,idim-1,1)
          ti0(j,k,2,3) = -tursav(j,k,idim-1,2)
          ti0(j,k,1,4) = 2.*ti0(j,k,1,3)-tursav(j,k,idim-1,1)
          ti0(j,k,2,4) = 2.*ti0(j,k,2,3)-tursav(j,k,idim-1,2)
  601   continue
      end if
      if (ivmx.ge.6) then
        c2b=cbar/tinf
        c2bp=c2b+1.0
        re=reue/xmach
        if (     ivmx.eq.6) then
          beta1=.075
        else if (ivmx.eq.7 .or. ivmx.eq.30 .or. ivmx.eq.40) then
          beta1=.075
        else if (ivmx.eq.8 .or. ivmx.eq.12 .or. ivmx.eq.14) then
          beta1=.83
        else if (ivmx.eq.72) then
          beta1=.0708
        end if
        i=idim-1
c       (Note:  linearly interpolates eps & omega BC to ghost cell.  Previous
c       versions of k-omega simply put 60*mu/(rho*beta*smin**2) at ghost cell 
c       rather than at wall, where it "officially" belongs.  However, this ends
c       up making very little, if any, difference in the k-omega solution.)
        if (ivmx.eq.9 .or. ivmx.eq.10 .or. ivmx.eq.11 .or.
     .      ivmx.eq.13) then
        do 632 k=ksta,kend1
        do 632 j=jsta,jend1
          tt=gamma*qi0(j,k,5,3)/qi0(j,k,1,3)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=sqrt(tursav(j,k,i,2))/ccabs(smin(j,k,i))
          ti0(j,k,1,3) = 2.*(2.*fnu/(q(j,k,i,1)*re**2)*dkdy**2)-
     +                   tursav(j,k,i,1)
          ti0(j,k,2,3) = -tursav(j,k,i,2)
          ti0(j,k,1,4) = 2.*ti0(j,k,1,3)-tursav(j,k,idim-1,1)
          ti0(j,k,2,4) = 2.*ti0(j,k,2,3)-tursav(j,k,idim-1,2)
 632    continue
        else if(ivmx.eq.15) then
        do k=ksta,kend1
        do j=jsta,jend1
          tt=gamma*q(j,k,i,5)/q(j,k,i,1)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dkdy=0.5*(tursav(j,k,i,2)+tursav(j,k,i-1,2))/3.
          ti0(j,k,1,3) = 2.0*dkdy/(re*smin(j,k,i))**2 - tursav(j,k,i,1)
          ti0(j,k,2,3) = -tursav(j,k,i,2)
          ti0(j,k,1,4) = 0.
          ti0(j,k,2,4) = 0.
        enddo
        enddo
        else if(ivmx.eq.16) then
        do k=ksta,kend1
        do j=jsta,jend1
          ti0(j,k,1,3) = -tursav(j,k,i,1)
          ti0(j,k,2,3) = -tursav(j,k,i,2)
          ti0(j,k,1,4) = 0.
          ti0(j,k,2,4) = 0.
        enddo
        enddo
#   ifdef CMPLX
c         RSM not working in complex mode
#   else
        else if(ivmx.eq.70) then  ! k-zeta Reynolds Stress model - xxiao
           do k=ksta,kend1
              do j=jsta,jend1
                 do iv=1,6
                    ti0(j,k,iv,3) = -tursav(j,k,i,iv)
                    ti0(j,k,iv,4) = huge(tke1)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 tke1 = -0.5*sum(tursav(j,k,i,1:3))
                 tke2 = -0.5*sum(tursav(j,k,i-1,1:3))
                 dkdy=0.5*(tke1+tke2)/3.*0.5
                 ti0(j,k,7,3) = 2.0*dkdy/(re*smin(j,k,i))**2 
     $                - tursav(j,k,i,7) 
                 ti0(j,k,7,4) = 0.0
              enddo
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.ne.6) then
           do k=ksta,kend1
              do j=jsta,jend1
                 do iv=1,6
                    ti0(j,k,iv,3) = -tursav(j,k,i,iv)
                    ti0(j,k,iv,4) = 2.*ti0(j,k,iv,3) - tursav(j,k,i,iv)
                 enddo
                 tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                 fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
                 dist=ccabs(smin(j,k,i))
                 ti0(j,k,7,3) = 
     $                2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +                 tursav(j,k,i,7)
                 ti0(j,k,7,4) = 2.*ti0(j,k,7,3) - tursav(j,k,i,7)
              enddo
           enddo
        else if(ivmx.eq.72 .and. issglrrw2012.eq.6) then
           do k=ksta,kend1
              do j=jsta,jend1
                 do iv=1,7
                    ti0(j,k,iv,3) = -tursav(j,k,i,iv)
                    ti0(j,k,iv,4) = 2.*ti0(j,k,iv,3) - tursav(j,k,i,iv)
                 enddo
              enddo
           enddo
#   endif
        else if(ivmx.eq.25) then
c       not used, but set anyway
        do k=ksta,kend1
        do j=jsta,jend1
          ti0(j,k,1,3) = 0.
          ti0(j,k,2,3) = 0.
          ti0(j,k,1,4) = 0.
          ti0(j,k,2,4) = 0.
        enddo
        enddo
        else
c       all other models (ivmx=6,7,8,12,14)
        do 602 k=ksta,kend1
        do 602 j=jsta,jend1
          tt=gamma*qi0(j,k,5,3)/qi0(j,k,1,3)
          fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
          dist=ccabs(smin(j,k,i))
          ti0(j,k,1,3) = 2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +                   tursav(j,k,i,1)
          ti0(j,k,2,3) = -tursav(j,k,i,2)
          ti0(j,k,1,4) = 2.*ti0(j,k,1,3)-tursav(j,k,idim-1,1)
          ti0(j,k,2,4) = 2.*ti0(j,k,2,3)-tursav(j,k,idim-1,2)
 602    continue
        end if
c       for trans model, intermittency at wall is either zero flux or set by user
        if (ivmx .eq. 30) then
        do k=ksta,kend1
        kk = k-ksta+1
        do j=jsta,jend1
        jj = j-jsta+1
          if (iuse3 .ne. 2) then
c           bc2004 or 2014 -> zero flux:
            ti0(j,k,3,3) = tursav(j,k,i,3)
          else
c           ti0(j,k,3,3) = 2.*bcdata(jj,kk,ip,3)-tursav(j,k,i,3)
c           bc2024 -> d(gamma)/dn=-(bcdata(3))*rho/mu*Uedge:
            tt=gamma*q(j,k,i,5)/q(j,k,i,1)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            hee=(1.+0.5*gm1*xmach**2)/gm1
            rhoee=(gamma*q(j,k,i,5))**(1./gamma)
            uee2=(gm1*rhoee*hee-gamma*q(j,k,i,5))/(0.5*gm1*rhoee)
            if (real(uee2) .gt. 0.) then
              uee=sqrt(uee2)
            else
              uee=xmach
            end if
            ti0(j,k,3,3) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
     +        uee*bcdata(jj,kk,ip,3)*reue/xmach
          end if
          ti0(j,k,3,4) = 2.*ti0(j,k,3,3)-tursav(j,k,idim-1,3)
        enddo
        enddo
        else if (ivmx .eq. 40) then
        do k=ksta,kend1
        kk = k-ksta+1
        do j=jsta,jend1
        jj = j-jsta+1
          ti0(j,k,3,3) = tursav(j,k,i,3)
          ti0(j,k,4,3) = tursav(j,k,i,4)
        enddo
        enddo
        end if
      end if
      end if
      end if
c
      return
      end
