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
      subroutine init(nbl,jdim,kdim,idim,q,qj0,qk0,qi0,tj0,tk0,ti0,
     .                vol,volj0,volk0,voli0,nummem,x,y,z,
     .                nou,bou,nbuf,ibufdim,iflagprnt)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Initialize the initial conditions on a mesh to be
c     freestream. Also initialize boundary volume arrays.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension tj0(kdim,idim-1,nummem,4),tk0(jdim,idim-1,nummem,4),
     .          ti0(jdim,kdim  ,nummem,4)
      dimension volj0(kdim,idim-1,4),volk0(jdim,idim-1,4),
     .          voli0(jdim,kdim  ,4),vol(jdim,kdim,idim-1)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c
      character*120 bou(ibufdim,nbuf)
      dimension nou(nbuf)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /reyue/ reue,tinf,ivisc(3)
      common /maxiv/ ivmx
      common /mms/ iexact_trunc,iexact_disc,iexact_ring
      common /reystressmodel/ issglrrw2012,i_sas_rsm,i_yapterm
c
c     write(15,904)nbl
  904 format(29h initial conditions for block,i3)
c
      rho0   = 1.e0
      c0     = 1.e0
      p0     = rho0*c0*c0/gamma
c   The wind axis system follows NASA SP-3070 (1972), with the exception that
c   positive beta is in the opposite direction
      u0     = xmach*cos(alpha)*cos(beta)
      w0     = xmach*sin(alpha)*cos(beta)
      v0     = -xmach*sin(beta)
      ei0    = p0/((gamma-1.e0)*rho0)
      et0    = rho0*(ei0+.5e0*(u0*u0+v0*v0+w0*w0))
      h0     = (et0+p0)/rho0
      pt0    = p0*(1.e0+.5e0*gm1*xmach*xmach)**(gamma/gm1)
      rhot0  = rho0*(1.e0+.5e0*gm1*xmach*xmach)**(1.e0/gm1)
      qiv(1) = rho0
      qiv(2) = u0
      qiv(3) = v0
      qiv(4) = w0
      qiv(5) = p0
c     write(15,499) idim,jdim,kdim
  499 format(/1x,20hinit idim,jdim,kdim=,3i10)
c     write(15,500) real(rho0),real(c0),real(p0),real(u0),real(v0),
c    .real(w0),real(ei0),real(et0),real(h0),real(pt0),real(qiv)
  500 format(1x,37h rho,c,p,u,v,w,ei,e,h,pt(0),qiv(1-5)=/(1x,5e12.5))
      a1    = u0
      a2    = v0
      a3    = w0
      jkpro = jdim*kdim
      idim1 = idim-1
      npl   = 999000/jkpro
      nplq  = npl
c
      do 1010 i=1,idim,nplq
      if (i+npl-1.gt.idim) npl = idim-i+1
      nv = npl*jkpro
cdir$ ivdep
      do 1000 izz=1,nv
      q(izz,1,i,1) = rho0
      q(izz,1,i,2) = a1
      q(izz,1,i,3) = a2
      q(izz,1,i,4) = a3
      q(izz,1,i,5) = p0
 1000 continue
 1010 continue
c
      npl  = 999000/kdim
      nplq = npl
      do 1020 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      nv   = npl*kdim
      do 10 m=1,4
cdir$ ivdep
      do 1001 izz=1,nv
      qj0(izz,i,1,m) = rho0
      qj0(izz,i,2,m) = a1
      qj0(izz,i,3,m) = a2
      qj0(izz,i,4,m) = a3
      qj0(izz,i,5,m) = p0
 1001 continue
   10 continue
 1020 continue
c
      npl  = 999000/kdim
      nplq = npl
      do 1030 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      nv   = jdim*npl
      do 20 m=1,4
cdir$ ivdep
      do 1002 izz=1,nv
      qk0(izz,i,1,m) = rho0
      qk0(izz,i,2,m) = a1
      qk0(izz,i,3,m) = a2
      qk0(izz,i,4,m) = a3
      qk0(izz,i,5,m) = p0
 1002 continue
   20 continue
 1030 continue
c
      do 30 m=1,4
cdir$ ivdep
      do 1003 izz=1,jkpro
      qi0(izz,1,1,m) = rho0
      qi0(izz,1,2,m) = a1
      qi0(izz,1,3,m) = a2
      qi0(izz,1,4,m) = a3
      qi0(izz,1,5,m) = p0
 1003 continue
   30 continue
c   For field eqn turbulence models:
      if (ivisc(1).ge.4 .or. ivisc(2).ge.4 .or. ivisc(3).ge.4) then
        if (ivisc(1).eq.4 .or. ivisc(2).eq.4 .or. ivisc(3).eq.4) then
          if (real(tur10(1)) .lt. 0.) tur10(1)=0.1
          if (real(tur10(2)) .lt. 0.) tur10(2)=0.
        else if (ivisc(1).eq.5 .or. ivisc(2).eq.5 .or.
     .           ivisc(3).eq.5) then
          if (real(tur10(1)) .lt. 0.) tur10(1)=1.341946
          if (real(tur10(2)) .lt. 0.) tur10(2)=0.
        else if (ivisc(1).eq.11 .or. ivisc(2).eq.11 .or.
     .           ivisc(3).eq.11 .or.
     .           ivisc(1).eq.10 .or. ivisc(2).eq.10 .or.
     .           ivisc(3).eq.10 .or.
     .           ivisc(1).eq. 9 .or. ivisc(2).eq. 9 .or.
     .           ivisc(3).eq. 9 .or. ivisc(1).eq.13 .or.
     .           ivisc(2).eq.13 .or. ivisc(3).eq.13 .or.
     .           ivisc(1).eq.15 .or. ivisc(2).eq.15 .or.
     .           ivisc(3).eq.15) then
          if (real(tur10(1)) .lt. 0.) tur10(1)=1.e-17
          if (real(tur10(2)) .lt. 0.) tur10(2)=1.e-9
        else if (ivisc(1).eq. 8 .or. ivisc(2).eq. 8 .or.
     .           ivisc(3).eq. 8 .or.
     .           ivisc(1).eq.12 .or. ivisc(2).eq.12 .or.
     .           ivisc(3).eq.12 .or.
     .           ivisc(1).eq.14 .or. ivisc(2).eq.14 .or.
     .           ivisc(3).eq.14) then
          if (real(tur10(1)) .lt. 0.) tur10(1)=9.e-8
          if (real(tur10(2)) .lt. 0.) tur10(2)=9.e-9
        else if (ivisc(1).eq.16 .or. ivisc(2).eq.16 .or.
     .           ivisc(3).eq.16) then
          if (real(tur10(1)) .lt. 0.) tur10(1)=1.5589e-6
          if (real(tur10(2)) .lt. 0.) tur10(2)=9.e-9
        else if (ivisc(1).eq.25 .or. ivisc(2).eq.25 .or.
     .           ivisc(3).eq.25) then
c         tur10 and tur20 not used for model 25, but we
c         set them anyway for safety:
          if (real(tur10(1)) .lt. 0.) tur10(1)=1.e-6
          if (real(tur10(2)) .lt. 0.) tur10(2)=9.e-9
        else if (ivisc(1).eq. 6 .or. ivisc(2).eq. 6 .or.
     .           ivisc(3).eq. 6 .or.
     .           ivisc(1).eq. 7 .or. ivisc(2).eq. 7 .or.
     .           ivisc(3).eq. 7) then
          if (real(tur10(1)) .lt. 0.) tur10(1)=1.e-6
          if (real(tur10(2)) .lt. 0.) tur10(2)=9.e-9
        else if (ivisc(1).eq.30 .or. ivisc(2).eq.30 .or.
     .           ivisc(3).eq.30) then
          if (real(tur10(1)) .lt. 0.) tur10(1)=1.e-6
          if (real(tur10(2)) .lt. 0.) tur10(2)=9.e-9
          if (real(tur10(3)) .lt. 0.) tur10(3)=1.0
        else if (ivisc(1).eq.40 .or. ivisc(2).eq.40 .or.
     .           ivisc(3).eq.40) then
          if (real(tur10(1)) .lt. 0.) tur10(1)=1.e-6
          if (real(tur10(2)) .lt. 0.) tur10(2)=9.e-9
          if (real(tur10(3)) .lt. 0.) tur10(3)=1.0
          if (real(tur10(4)) .lt. 0.) then
            turbintensity_inf_percent=100./xmach*sqrt(tur10(2)/1.5)
            if (real(turbintensity_inf_percent) .le. 1.3) then
             tur10(4)=1173.51-589.428*turbintensity_inf_percent+
     .          0.2196/(turbintensity_inf_percent**2)
            else
             tur10(4)=331.5*(turbintensity_inf_percent-0.5658)**(-0.671)
            end if
          end if
        else if (ivisc(1).eq.72 .or. ivisc(2).eq.72 .or.
     .           ivisc(3).eq.72) then
          if (real(tur10(1)) .eq. -1.) tur10(1)=-2./3.*9.e-9
          if (real(tur10(2)) .eq. -1.) tur10(2)=-2./3.*9.e-9
          if (real(tur10(3)) .eq. -1.) tur10(3)=-2./3.*9.e-9
          if (real(tur10(4)) .eq. -1.) tur10(4)=0.
          if (real(tur10(5)) .eq. -1.) tur10(5)=0.
          if (real(tur10(6)) .eq. -1.) tur10(6)=0.
          if (issglrrw2012 .eq. 6) then
            if (real(tur10(7)) .eq. -1.) tur10(7)=1.0/sqrt(1.e-6)
          else
            if (real(tur10(7)) .eq. -1.) tur10(7)=1.e-6
          end if
c    Write out some diagnostic info:
c    CFL3D has default tur10 values; can be overridden via
c    turbintensity_inf_percent and/or eddy_visc_inf (see readkey.F)
          if (nbl .eq. 1 .and. iflagprnt .eq. 0) then
          tke=-(tur10(1)+tur10(2)+tur10(3))/2.
          tu=sqrt(tke/1.5)/xmach
          xnu_ratio=tke/tur10(7)
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' For Stress-Omega RSM'',
     .      '', ivisc=72 (all blocks):'')')
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' initial Tu ='',es12.5)') tu
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' initial xnu_ratio ='',es12.5)') 
     .      xnu_ratio
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' initial tke ='',es12.5)') tke
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' initial omega ='',es12.5)') tur10(7)
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' tur10(1) ='',es12.5)') tur10(1)
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' tur10(2) ='',es12.5)') tur10(2)
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' tur10(3) ='',es12.5)') tur10(3)
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' tur10(4) ='',es12.5)') tur10(4)
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' tur10(5) ='',es12.5)') tur10(5)
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' tur10(6) ='',es12.5)') tur10(6)
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' tur10(7) ='',es12.5)') tur10(7)
          end if
        end if
        do 4003 nn=1,nummem
        do 4001 m=1,4
          do 4000 i=1,idim-1
            do 3999 j=1,jdim
              tk0(j,i,nn,m)=tur10(nn)
 3999       continue
            do 3998 k=1,kdim
              tj0(k,i,nn,m)=tur10(nn)
 3998       continue
 4000     continue
          do 4002 k=1,kdim
            do 4002 j=1,jdim
              ti0(j,k,nn,m)=tur10(nn)
 4002     continue
 4001   continue
 4003   continue
c
      end if
c
c     boundary volumes (default to interior cells)
c
      if (ivmx.gt.0) then
         do m=1,4
            if (m.eq.1) kk = 1
            if (m.eq.2) kk = min(2,kdim-1)
            if (m.eq.3) kk = kdim-1
            if (m.eq.4) kk = max(1,kdim-2)
            do i=1,idim-1
               do j=1,jdim
                 volk0(j,i,m) = vol(j,kk,i)
               end do
            end do
         end do
         do m=1,4
            if (m.eq.1) jj = 1
            if (m.eq.2) jj = min(2,jdim-1)
            if (m.eq.3) jj = jdim-1
            if (m.eq.4) jj = max(1,jdim-2)
c           the sgi f77 compiler with -O3 won't do the following 
c           assigmnment correctly if the k-loop is innermost
            do k=1,kdim
               do i=1,idim-1
                 volj0(k,i,m) = vol(jj,k,i)
               end do
            end do
         end do
         do m=1,4
            if (m.eq.1) ii = 1
            if (m.eq.2) ii = min(2,idim-1)
            if (m.eq.3) ii = idim-1
            if (m.eq.4) ii = max(1,idim-2)
            do k=1,kdim
               do j=1,jdim
                 voli0(j,k,m) = vol(j,k,ii)
               end do
            end do
         end do
      end if
c
c   Overwrite with exact soln if doing MMS
      if (iexact_trunc .ne. 0 .or. iexact_disc .ne. 0) then
        call exact_flow_q(jdim,kdim,idim,x,y,z,q,iexact_trunc,
     +                    iexact_disc)
      end if
c
      return
      end
