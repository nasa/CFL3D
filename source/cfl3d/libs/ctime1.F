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
      subroutine ctime1(nbl,jdim,kdim,idim,q,vol,sj,sk,si,dtj,t,delt,
     .                 vist3d,itur,dtmin,iout,ntime,nou,bou,nbuf,
     .                 ibufdim,idef)
c
c     $Id$
c
c************************************************************************
c     Purpose:  Calculate the time step for an input fixed Courant number
c     or calculate the Courant number based on an input value of delta t.
c
c     formerly called ctime...changed to ctime1 to prevent warning
c     messages about conflicts with intrinsic ctime in SGI libraries
c
c     Modified for Weiss-Smith preconditioning by J.R. Edwards, NCSU
c       cprec = 0 ---> original code used
c             > 0 ---> modified code used
c************************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension vol(jdim,kdim,idim-1),vist3d(jdim,kdim,idim)
      dimension sj(jdim*kdim,idim-1,5),sk(jdim*kdim,idim-1,5),
     .          si(jdim*kdim,idim,5)
      dimension q(jdim,kdim,idim,5)
      dimension t(jdim*kdim,17),dtj(jdim,kdim,idim-1)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /reyue/ reue,tinf,ivisc(3)
      common /twod/ i2d
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /precond/ cprec,uref,avn
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
c
      dttemp=dt
c
c     set flags for contributions from pseudo- and physical- time steps
c
      if (real(dt).gt.0.0) then
         if (ita.lt.0) then
            if( real(cfltauMax) > real(cfltau) .and. ncyc > 2 ) then
               dt = -(cfltau +
     $              (cfltauMax-cfltau)*((icyc-1.)/(ncyc-1.))**cfltau0)
            else
               dt=-cfltau
            end if
            fact1 = 1.0
            fact2 = 1.0
         else
            fact1 = 0.0
            fact2 = 1.0
         end if
         if (abs(ita) .eq. 1) then
           tfact=0.e0
         else
           tfact=0.5e0
         end if
         tfacp1=tfact+1.e0
      else
         fact1 = 1.0
         fact2 = 0.0
         tfacp1=1.e0
      end if
c
      if (iout.gt.0 .and. real(dttemp).lt.0.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),7) nbl
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),17)
      end if
    7 format(1x,42hcomputing time step distribution for block,i5)
   17 format(3x,25hsummary of time step data)
c
      if (iout.gt.0 .and. real(dttemp).gt.0.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),8) nbl
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),18)
      end if
    8 format(1x,43hcomputing CFL number distribution for block,i5)
   18 format(3x,26hsummary of CFL number data)
c
      vterm = ccmaxrc( 1.3333 , gamma/pr )*xmach/reue
      term  = 1.
      if (i2d.eq.1 .or. iaxi2plane.eq.1) term = 0.
      delt  = 0.0
      if (iunst.gt.0.and.real(dttemp).gt.0.0) delt = dttemp
      cfl1  = 1.e0/ccabs(dt)
      dt2   = dt*dt
      n     = jdim*kdim-jdim-1
c
      if (iout.gt.0) then
         if (real(dttemp).lt.0.) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1241)
         else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),1247)
         end if
      end if
 1241 format(3x,1hI,6x,3hcfl,10x,5hdtrms,9x,5hdtmin,9x,5hdtmax)
 1247 format(3x,1hI,7x,2hdt,9x,6hcflrms,8x,6hcflmin,8x,6hcflmax)
c
      do 9000 i=1,idim1
c
      n = jdim*kdim - jdim - 1
c
      if (itur.eq.0) then
         do 200 izz = 1,n
         t(izz,1)    = 0.
  200    continue
      else
         do 350 k=1,kdim-1
         izz           = jdim*(k-1)
         do 300 j=1,jdim-1
         t(izz+j,2)    = gamma*q(j,k,i,5)/q(j,k,i,1)
  300    continue
         call xmukin(jdim-1,t(izz+1,2),t(izz+1,1),tinf)
         t(izz+jdim,1) = t(izz+jdim-1,1)
  350    continue
      end if
c
      if (itur.ge.2) then
         do 450 k=1,kdim-1
         izz           = jdim*(k-1) 
         do 400 j=1,jdim-1
         t(izz+j,1)    = t(izz+j,1) + vist3d(j,k,i)
  400    continue
  450    continue
      end if
c
      if (real(cprec) .eq. 0) then
cdir$ ivdep
         do 1000 izz=1,n
         vrho = t(izz,1)*vterm/( vol(izz,1,i)*q(izz,1,i,1) )
c
         t11 = sj(izz,i,1)+sj(izz+1,i,1)
         t12 = sj(izz,i,2)+sj(izz+1,i,2)
         t13 = sj(izz,i,3)+sj(izz+1,i,3)
         t14 = sj(izz,i,4)+sj(izz+1,i,4)
         t1  = (t11*q(izz,1,i,2)+t12*q(izz,1,i,3)
     .       +  t13*q(izz,1,i,4))*0.5e0
c        add cell face speed for moving grids
         t20 = sj(izz,i,5)+sj(izz+1,i,5)
         t1  = t1+t20*0.5e0
         t1  = ccabs(t1) + t14*vrho
c
         t11 = sk(izz,i,1)+sk(izz+jdim,i,1)
         t12 = sk(izz,i,2)+sk(izz+jdim,i,2)
         t13 = sk(izz,i,3)+sk(izz+jdim,i,3)
         t15 = sk(izz,i,4)+sk(izz+jdim,i,4)
         t2  = (t11*q(izz,1,i,2)+t12*q(izz,1,i,3)
     .       +  t13*q(izz,1,i,4))*0.5e0
c        add cell face speed for moving grids
         t20 = sk(izz,i,5)+sk(izz+jdim,i,5)
         t2  = t2+t20*0.5e0
         t2  = ccabs(t2) + t15*vrho
c
         t11 = si(izz,i,1)+si(izz,i+1,1)
         t12 = si(izz,i,2)+si(izz,i+1,2)
         t13 = si(izz,i,3)+si(izz,i+1,3)
         t16 = si(izz,i,4)+si(izz,i+1,4)
         t3  = (t11*q(izz,1,i,2)+t12*q(izz,1,i,3)
     .       +  t13*q(izz,1,i,4))*0.5e0
c        add cell face speed for moving grids
         t20 = si(izz,i,5)+si(izz,i+1,5)
         t3  = t3+t20*0.5e0
         t3  = ccabs(t3) + t16*vrho
c
         t11 = gamma*q(izz,1,i,5)/q(izz,1,i,1)
         t11 = sqrt(t11)
c
         dtj(izz,1,i) = (t1+t11)*t14+(t2+t11)*t15+(t3+t11)*t16*term
         dtj(izz,1,i) = 0.5e0*cfl1*dtj(izz,1,i)
         t(izz,10)    = vol(izz,1,i)/dtj(izz,1,i)*ccabs(dttemp/dt)
 1000 continue
      else
cdir$ ivdep
         do 10001 izz=1,n
         vrho = t(izz,1)*vterm/( vol(izz,1,i)*q(izz,1,i,1) )
         c2 = gamma*q(izz,1,i,5)/q(izz,1,i,1)
         c = sqrt(c2)
         vmag1 =  q(izz,1,i,2)**2 + q(izz,1,i,3)**2 + q(izz,1,i,4)**2
         vel2 = ccmax(vmag1,avn*uref**2)
         vel = sqrt(ccmin(c2,vel2))
         vel = cprec*vel + (1.-cprec)*c
         xm2 = (vel/c)**2
c
         t11 = sj(izz,i,1)+sj(izz+1,i,1)
         t12 = sj(izz,i,2)+sj(izz+1,i,2)
         t13 = sj(izz,i,3)+sj(izz+1,i,3)
         t14 = sj(izz,i,4)+sj(izz+1,i,4)
         t1  = (t11*q(izz,1,i,2)+t12*q(izz,1,i,3)
     .       +  t13*q(izz,1,i,4))*0.5e0
         xmave = t1/c
         tt1j = 0.5*(1.+xm2)
         tt2j = 0.5*sqrt(xmave**2*(1.-xm2)**2 + 4.0*xm2)
c        add cell face speed for moving grids
         t20 = sj(izz,i,5)+sj(izz+1,i,5)
         t1  = tt1j*t1+t20*0.5e0
         t1  = ccabs(t1) + t14*vrho
c
         t11 = sk(izz,i,1)+sk(izz+jdim,i,1)
         t12 = sk(izz,i,2)+sk(izz+jdim,i,2)
         t13 = sk(izz,i,3)+sk(izz+jdim,i,3)
         t15 = sk(izz,i,4)+sk(izz+jdim,i,4)
         t2  = (t11*q(izz,1,i,2)+t12*q(izz,1,i,3)
     .       +  t13*q(izz,1,i,4))*0.5e0
         xmave = t2/c
         tt1k = 0.5*(1.+xm2)
         tt2k = 0.5*sqrt(xmave**2*(1.-xm2)**2 + 4.0*xm2)
c        add cell face speed for moving grids
         t20 = sk(izz,i,5)+sk(izz+jdim,i,5)
         t2  = tt1k*t2+t20*0.5e0
         t2  = ccabs(t2) + t15*vrho
c
         t11 = si(izz,i,1)+si(izz,i+1,1)
         t12 = si(izz,i,2)+si(izz,i+1,2)
         t13 = si(izz,i,3)+si(izz,i+1,3)
         t16 = si(izz,i,4)+si(izz,i+1,4)
         t3  = (t11*q(izz,1,i,2)+t12*q(izz,1,i,3)
     .       +  t13*q(izz,1,i,4))*0.5e0
         xmave = t3/c
         tt1i = 0.5*(1.+xm2)
         tt2i = 0.5*sqrt(xmave**2*(1.-xm2)**2 + 4.0*xm2)
c        add cell face speed for moving grids
         t20 = si(izz,i,5)+si(izz,i+1,5)
         t3  = tt1i*t3+t20*0.5e0
         t3  = ccabs(t3) + t16*vrho
c
         c1j = c*tt2j
         c1k = c*tt2k
         c1i = c*tt2i
c
         dtj(izz,1,i) = (t1+c1j)*t14+(t2+c1k)*t15+(t3+c1i)*t16*term
         dtj(izz,1,i) = 0.5e0*cfl1*dtj(izz,1,i)
         t(izz,10)    = vol(izz,1,i)/dtj(izz,1,i)*ccabs(dttemp/dt)
10001    continue
      end if
      if (real(dttemp).gt.0.e0) then
         dt2 = dttemp*dttemp
cdir$ ivdep
         do 1001 izz=1,n
         t(izz,10)    = dt2/t(izz,10)
c   code can only do 1st order temporal for pseudo-time term (with 
c   subiterations).  Therefore, pseudo-time term is divided here by
c   tfacp1.  In several subroutines in af3f, dtj is multiplied by
c   tfacp1.
         dtj(izz,1,i) = fact1*dtj(izz,1,i)/tfacp1 + 
     .                  fact2*vol(izz,1,i)/dttemp
 1001    continue
         dt2 = dt*dt
c
c        geometric conservation law terms for deforming grids
c
         if (idef.gt.0) then
            do 1002 izz=1,n
            t(izz,17) = sj(izz+1,i,5)*sj(izz+1,i,4)
     .                - sj(izz,i,5)*sj(izz,i,4)
     .                + sk(izz+jdim,i,5)*sk(izz+jdim,i,4)
     .                - sk(izz,i,5)*sk(izz,i,4)
     .                + si(izz,i+1,5)*si(izz,i+1,4)
     .                - si(izz,i,5)*si(izz,i,4)
1002        continue
            do 1005 izz=jdim,n,jdim
            t(izz,17) = sj(izz,i,5)*sj(izz,i,4)
     .                - sj(izz-1,i,5)*sj(izz-1,i,4)
     .                + sk(izz+jdim,i,5)*sk(izz+jdim,i,4)
     .                - sk(izz,i,5)*sk(izz,i,4)
     .                + si(izz,i+1,5)*si(izz,i+1,4)
     .                - si(izz,i,5)*si(izz,i,4)
1005        continue
            do 1006 izz=1,n
            dtj(izz,1,i)= dtj(izz,1,i)+t(izz,17)/tfacp1
1006        continue
         end if
      end if
c
      do 2000 kk=1,kdim1
      jk       = jdim*kk
      t(jk,10) = 0.e0
 2000 continue
      dtrms    = q8sdot(n,t(1,10),n,t(1,10))
      dtrms    = sqrt(dtrms/float(jdim1*kdim1))
      dtmax    = q8smax(n,t(1,10))
      do 3000 kk=1,kdim1
      jk       = jdim*kk
      t(jk,10) = dtmax
 3000 continue
      dtmin    = q8smin(n,t(1,10))
      dtpr     = ccabs(dttemp)
      if (iout.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1201) i,real(dtpr),real(dtrms),
     .                             real(dtmin),real(dtmax)
      end if
 1201 format(i4,1x,e12.5,3(2x,e12.5))
c
cdir$ ivdep
      do 1003 izz=1,jdim+1
      dtj(izz+jdim-1,kdim1,i) = 1.0e0
 1003 continue
 9000 continue
      if (dt.ne.dttemp) then
        dt=dttemp
      end if
c
      return
      end
