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
      subroutine bc2010(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,ista,iend,jsta,jend,ksta,kend,nface,
     .                  tursav,tj0,tk0,ti0,vist3d,vj0,vk0,vi0,
     .                  mdim,ndim,bcdata,filname,iuns,
     .                  nou,bou,nbuf,ibufdim,myid,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Set inflow boundary conditions (typically for nozzle,
c              duct or engine flows), given total pressure ratio, 
c              total temperature and flow angle. The pressure is 
c              extrapolated (zeroth order) from the interior of 
c              the domain, and the remaining variables are determined
c              from the extrapolated pressure and the input data.
c              Also can input up to two turbulence quantities.
c
c     pte   = total to free stream static pressure ratio at 
c             inlet (Ptotal/pinf)
c     tte   = total to free stram static temperature ratio at 
c             inlet (Ttotal/tinf)
c     alpe  = alpha angle of inlet flow, deg
c     betae = beta angle of inlet flow, deg
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
      dimension bcdata(mdim,ndim,2,12)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
c
      common /maxiv/ ivmx
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /conversion/ radtodeg
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
c            * standard boundary condition bctype=2010 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary             nozzle total BCs                   type 2010
c******************************************************************************
c
      if (nface.eq.3) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,5)) .gt. -1.e10) itrflg1 = 1
c
      j=1
      do 300 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(kend-ksta)+1
c
      do 800 k=ksta,kend1
      kk = k-ksta+1
c
      pte   = bcdata(kk,ii,ip,1)
      tte   = bcdata(kk,ii,ip,2)
      alpe  = bcdata(kk,ii,ip,3)
      betae = bcdata(kk,ii,ip,4)
c
      pte=pte*p0
      alpe=alpe/radtodeg
      betae=betae/radtodeg
c
c   The following method is used by FUN3D:
      unormi = q(j,k,i,2)*sj(1,k,i,1) +
     +         q(j,k,i,3)*sj(1,k,i,2) +
     +         q(j,k,i,4)*sj(1,k,i,3)
      umag_sq = q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)*2
      umag = sqrt(umag_sq)
      ai = sqrt(gamma*q(j,k,i,5)/q(j,k,i,1))
      rminus = abs(unormi) - (2.*ai)/gm1
c     Note: xmach_local should be < 1 (not checking here)
      xmach_local = umag/ai
      hti = (q(j,k,i,5)/q(j,k,i,1))*(gamma/gm1)+0.5*umag_sq
      qa = 1.+2./gm1
      qb = 2.*rminus
      qc = 0.5*gm1*rminus**2 - gm1*hti
      a_radical = qb**2 - 4.*qa*qc
      if (real(a_radical) .gt. 0.) then
        a_plus = -qb/(2.*qa)+sqrt(a_radical)/(2.*qa)
        a_minus = -qb/(2.*qa)-sqrt(a_radical)/(2.*qa)
      else if (real(a_radical) .lt. 0.) then
        a_plus = ai
        a_minus = ai
      else
        a_plus = -qb/(2.*qa)
        a_minus = -qb/(2.*qa)
      end if
      if (real(a_plus) .gt. real(a_minus)) then
        a_new = a_plus
      else
        a_new = a_minus
      end if
      u_new = ((2.*a_new)/gm1) + rminus
      xmach_new = u_new/a_new
      pressure = pte*(1.+0.5*gm1*xmach_new**2)**(-gamma/(gamma-1.))
      temperature = tte*(pressure/pte)**(gm1/gamma)
      qj0(k,i,1,1) = gamma*pressure/temperature
      qj0(k,i,2,1) = u_new*cos(alpe)*cos(betae)
      qj0(k,i,3,1) = -u_new*sin(betae)
      qj0(k,i,4,1) = u_new*sin(alpe)*cos(betae)
      qj0(k,i,5,1) = pressure
      qj0(k,i,1,2) = qj0(k,i,1,1)
      qj0(k,i,2,2) = qj0(k,i,2,1)
      qj0(k,i,3,2) = qj0(k,i,3,1)
      qj0(k,i,4,2) = qj0(k,i,4,1)
      qj0(k,i,5,2) = qj0(k,i,5,1)
c
      bcj(k,i,1)   = 0.0
c
  800 continue
  300 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,2) = 0.0
  191   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 101 i=ista,iend1
        ii = i-ista+1
        do 101 k=ksta,kend1
          kk=k-ksta+1
          ubar=-(qj0(k,i,2,1)*sj(1,k,i,1)+qj0(k,i,3,1)*sj(1,k,i,2)+
     +           qj0(k,i,4,1)*sj(1,k,i,3))
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,4+l)
          if (real(ubar) .lt. 0.) then
             tj0(k,i,l,1) = t1
             tj0(k,i,l,2) = t1
          else
             tj0(k,i,l,1) = tursav(1,k,i,l)
             tj0(k,i,l,2) = tj0(k,i,l,1)
          end if
  101   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      j=jdim boundary          nozzle total BCs                   type 2010
c******************************************************************************
c
      if (nface.eq.4) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,5)) .gt. -1.e10) itrflg1 = 1
c
      j=jdim1
      do 310 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(kend-ksta)+1
c
      do 810 k=ksta,kend1
      kk = k-ksta+1
c
      pte   = bcdata(kk,ii,ip,1)
      tte   = bcdata(kk,ii,ip,2)
      alpe  = bcdata(kk,ii,ip,3)
      betae = bcdata(kk,ii,ip,4)
c
      pte=pte*p0
      alpe=alpe/radtodeg
      betae=betae/radtodeg
c
c   The following method is used by FUN3D:
      unormi = q(j,k,i,2)*sj(jdim,k,i,1) +
     +         q(j,k,i,3)*sj(jdim,k,i,2) +
     +         q(j,k,i,4)*sj(jdim,k,i,3)
      umag_sq = q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)*2
      umag = sqrt(umag_sq)
      ai = sqrt(gamma*q(j,k,i,5)/q(j,k,i,1))
      rminus = abs(unormi) - (2.*ai)/gm1
c     Note: xmach_local should be < 1 (not checking here)
      xmach_local = umag/ai
      hti = (q(j,k,i,5)/q(j,k,i,1))*(gamma/gm1)+0.5*umag_sq
      qa = 1.+2./gm1
      qb = 2.*rminus
      qc = 0.5*gm1*rminus**2 - gm1*hti
      a_radical = qb**2 - 4.*qa*qc
      if (real(a_radical) .gt. 0.) then
        a_plus = -qb/(2.*qa)+sqrt(a_radical)/(2.*qa)
        a_minus = -qb/(2.*qa)-sqrt(a_radical)/(2.*qa)
      else if (real(a_radical) .lt. 0.) then
        a_plus = ai
        a_minus = ai
      else
        a_plus = -qb/(2.*qa)
        a_minus = -qb/(2.*qa)
      end if
      if (real(a_plus) .gt. real(a_minus)) then
        a_new = a_plus
      else
        a_new = a_minus
      end if
      u_new = ((2.*a_new)/gm1) + rminus
      xmach_new = u_new/a_new
      pressure = pte*(1.+0.5*gm1*xmach_new**2)**(-gamma/(gamma-1.))
      temperature = tte*(pressure/pte)**(gm1/gamma)
      qj0(k,i,1,3) = gamma*pressure/temperature
      qj0(k,i,2,3) = u_new*cos(alpe)*cos(betae)
      qj0(k,i,3,3) = -u_new*sin(betae)
      qj0(k,i,4,3) = u_new*sin(alpe)*cos(betae)
      qj0(k,i,5,3) = pressure
      qj0(k,i,1,4) = qj0(k,i,1,3)
      qj0(k,i,2,4) = qj0(k,i,2,3)
      qj0(k,i,3,4) = qj0(k,i,3,3)
      qj0(k,i,4,4) = qj0(k,i,4,3)
      qj0(k,i,5,4) = qj0(k,i,5,3)
c
      bcj(k,i,2)   = 0.0
c
  810 continue
  310 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
          vj0(k,i,1,3) = vist3d(jdim-1,k,i)
          vj0(k,i,1,4) = 0.0
  291   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 201 i=ista,iend1
        ii = i-ista+1
        do 201 k=ksta,kend1
          kk = k-ksta+1
          ubar=qj0(k,i,2,3)*sj(jdim,k,i,1)+qj0(k,i,3,3)*sj(jdim,k,i,2)+
     +         qj0(k,i,4,3)*sj(jdim,k,i,3)
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,4+l)
          if (real(ubar) .lt. 0.) then
             tj0(k,i,l,3) = t1
             tj0(k,i,l,4) = t1
          else
             tj0(k,i,l,3) = tursav(jdim-1,k,i,l)
             tj0(k,i,l,4) = tj0(k,i,l,3)
          end if
  201   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=1 boundary             nozzle total BCs                   type 2010
c******************************************************************************
c
      if (nface.eq.5) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,5)) .gt. -1.e10) itrflg1 = 1
c
      k=1
      do 320 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(jend-jsta)+1
      do 820 j=jsta,jend1
      jj = j-jsta+1
c
      pte   = bcdata(jj,ii,ip,1)
      tte   = bcdata(jj,ii,ip,2)
      alpe  = bcdata(jj,ii,ip,3)
      betae = bcdata(jj,ii,ip,4)
c
      pte=pte*p0
      alpe=alpe/radtodeg
      betae=betae/radtodeg
c
c   The following method is used by FUN3D:
      unormi = q(j,k,i,2)*sk(j,1,i,1) +
     +         q(j,k,i,3)*sk(j,1,i,2) +
     +         q(j,k,i,4)*sk(j,1,i,3)
      umag_sq = q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)*2
      umag = sqrt(umag_sq)
      ai = sqrt(gamma*q(j,k,i,5)/q(j,k,i,1))
      rminus = abs(unormi) - (2.*ai)/gm1
c     Note: xmach_local should be < 1 (not checking here)
      xmach_local = umag/ai
      hti = (q(j,k,i,5)/q(j,k,i,1))*(gamma/gm1)+0.5*umag_sq
      qa = 1.+2./gm1
      qb = 2.*rminus
      qc = 0.5*gm1*rminus**2 - gm1*hti
      a_radical = qb**2 - 4.*qa*qc
      if (real(a_radical) .gt. 0.) then
        a_plus = -qb/(2.*qa)+sqrt(a_radical)/(2.*qa)
        a_minus = -qb/(2.*qa)-sqrt(a_radical)/(2.*qa)
      else if (real(a_radical) .lt. 0.) then
        a_plus = ai
        a_minus = ai
      else
        a_plus = -qb/(2.*qa)
        a_minus = -qb/(2.*qa)
      end if
      if (real(a_plus) .gt. real(a_minus)) then
        a_new = a_plus
      else
        a_new = a_minus
      end if
      u_new = ((2.*a_new)/gm1) + rminus
      xmach_new = u_new/a_new
      pressure = pte*(1.+0.5*gm1*xmach_new**2)**(-gamma/(gamma-1.))
      temperature = tte*(pressure/pte)**(gm1/gamma)
      qk0(j,i,1,1) = gamma*pressure/temperature
      qk0(j,i,2,1) = u_new*cos(alpe)*cos(betae)
      qk0(j,i,3,1) = -u_new*sin(betae)
      qk0(j,i,4,1) = u_new*sin(alpe)*cos(betae)
      qk0(j,i,5,1) = pressure
      qk0(j,i,1,2) = qk0(j,i,1,1)
      qk0(j,i,2,2) = qk0(j,i,2,1)
      qk0(j,i,3,2) = qk0(j,i,3,1)
      qk0(j,i,4,2) = qk0(j,i,4,1)
      qk0(j,i,5,2) = qk0(j,i,5,1)
c
      bck(j,i,1)   = 0.0
c
  820 continue
  320 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = 0.0
  391   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 301 i=ista,iend1
        ii = i-ista+1
        do 301 j=jsta,jend1
          jj = j-jsta+1
          ubar=-(qk0(j,i,2,1)*sk(j,1,i,1)+qk0(j,i,3,1)*sk(j,1,i,2)+
     +           qk0(j,i,4,1)*sk(j,1,i,3))
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,4+l)
          if (real(ubar) .lt. 0.) then
             tk0(j,i,l,1) = t1
             tk0(j,i,l,2) = t1
          else
             tk0(j,i,l,1) = tursav(j,1,i,l)
             tk0(j,i,l,2) = tk0(j,i,l,1)
          end if
  301   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=kdim boundary          nozzle total BCs                   type 2010
c******************************************************************************
c
      if (nface.eq.6) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,5)) .gt. -1.e10) itrflg1 = 1
c
      k=kdim1
      do 330 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(jend-jsta)+1
c
      do 830 j=jsta,jend1
      jj = j-jsta+1
c
      pte   = bcdata(jj,ii,ip,1)
      tte   = bcdata(jj,ii,ip,2)
      alpe  = bcdata(jj,ii,ip,3)
      betae = bcdata(jj,ii,ip,4)
c
      pte=pte*p0
      alpe=alpe/radtodeg
      betae=betae/radtodeg
c
c
c   The following method is used by FUN3D:
      unormi = q(j,k,i,2)*sk(j,kdim,i,1) +
     +         q(j,k,i,3)*sk(j,kdim,i,2) +
     +         q(j,k,i,4)*sk(j,kdim,i,3)
      umag_sq = q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)*2
      umag = sqrt(umag_sq)
      ai = sqrt(gamma*q(j,k,i,5)/q(j,k,i,1))
      rminus = abs(unormi) - (2.*ai)/gm1
c     Note: xmach_local should be < 1 (not checking here)
      xmach_local = umag/ai
      hti = (q(j,k,i,5)/q(j,k,i,1))*(gamma/gm1)+0.5*umag_sq
      qa = 1.+2./gm1
      qb = 2.*rminus
      qc = 0.5*gm1*rminus**2 - gm1*hti
      a_radical = qb**2 - 4.*qa*qc
      if (real(a_radical) .gt. 0.) then
        a_plus = -qb/(2.*qa)+sqrt(a_radical)/(2.*qa)
        a_minus = -qb/(2.*qa)-sqrt(a_radical)/(2.*qa)
      else if (real(a_radical) .lt. 0.) then
        a_plus = ai
        a_minus = ai
      else
        a_plus = -qb/(2.*qa)
        a_minus = -qb/(2.*qa)
      end if
      if (real(a_plus) .gt. real(a_minus)) then
        a_new = a_plus
      else
        a_new = a_minus
      end if
      u_new = ((2.*a_new)/gm1) + rminus
      xmach_new = u_new/a_new
      pressure = pte*(1.+0.5*gm1*xmach_new**2)**(-gamma/(gamma-1.))
      temperature = tte*(pressure/pte)**(gm1/gamma)
      qk0(j,i,1,3) = gamma*pressure/temperature
      qk0(j,i,2,3) = u_new*cos(alpe)*cos(betae)
      qk0(j,i,3,3) = -u_new*sin(betae)
      qk0(j,i,4,3) = u_new*sin(alpe)*cos(betae)
      qk0(j,i,5,3) = pressure
      qk0(j,i,1,4) = qk0(j,i,1,3)
      qk0(j,i,2,4) = qk0(j,i,2,3)
      qk0(j,i,3,4) = qk0(j,i,3,3)
      qk0(j,i,4,4) = qk0(j,i,4,3)
      qk0(j,i,5,4) = qk0(j,i,5,3)
c
      bck(j,i,2)   = 0.0
c
  830 continue
  330 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
          vk0(j,i,1,3) = vist3d(j,kdim-1,i)
          vk0(j,i,1,4) = 0.0
  491   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 401 i=ista,iend1
        ii = i-ista+1
        do 401 j=jsta,jend1
          jj = j-jsta+1
          ubar=qk0(j,i,2,3)*sk(j,kdim,i,1)+qk0(j,i,3,3)*sk(j,kdim,i,2)+
     +         qk0(j,i,4,3)*sk(j,kdim,i,3)
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,4+l)
          if (real(ubar) .lt. 0.) then
             tk0(j,i,l,3) = t1
             tk0(j,i,l,4) = t1
          else
             tk0(j,i,l,3) = tursav(j,kdim-1,i,l)
             tk0(j,i,l,4) = tk0(j,i,l,3)
          end if
  401   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=1 boundary             nozzle total BCs                   type 2010
c******************************************************************************
c
      if (nface.eq.1) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,5)) .gt. -1.e10) itrflg1 = 1
c
      i=1
      do 340 k=ksta,kend1
      kk = k-ksta+1
      js = (k-ksta)*(jend-jsta)+1
c
      do 840 j=jsta,jend1
      jj = j-jsta+1
c
      pte   = bcdata(jj,kk,ip,1)
      tte   = bcdata(jj,kk,ip,2)
      alpe  = bcdata(jj,kk,ip,3)
      betae = bcdata(jj,kk,ip,4)
c
      pte=pte*p0
      alpe=alpe/radtodeg
      betae=betae/radtodeg
c
c   The following method is used by FUN3D:
      unormi = q(j,k,i,2)*si(j,k,1,1) +
     +         q(j,k,i,3)*si(j,k,1,2) +
     +         q(j,k,i,4)*si(j,k,1,3)
      umag_sq = q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)*2
      umag = sqrt(umag_sq)
      ai = sqrt(gamma*q(j,k,i,5)/q(j,k,i,1))
      rminus = abs(unormi) - (2.*ai)/gm1
c     Note: xmach_local should be < 1 (not checking here)
      xmach_local = umag/ai
      hti = (q(j,k,i,5)/q(j,k,i,1))*(gamma/gm1)+0.5*umag_sq
      qa = 1.+2./gm1
      qb = 2.*rminus
      qc = 0.5*gm1*rminus**2 - gm1*hti
      a_radical = qb**2 - 4.*qa*qc
      if (real(a_radical) .gt. 0.) then
        a_plus = -qb/(2.*qa)+sqrt(a_radical)/(2.*qa)
        a_minus = -qb/(2.*qa)-sqrt(a_radical)/(2.*qa)
      else if (real(a_radical) .lt. 0.) then
        a_plus = ai
        a_minus = ai
      else
        a_plus = -qb/(2.*qa)
        a_minus = -qb/(2.*qa)
      end if
      if (real(a_plus) .gt. real(a_minus)) then
        a_new = a_plus
      else
        a_new = a_minus
      end if
      u_new = ((2.*a_new)/gm1) + rminus
      xmach_new = u_new/a_new
      pressure = pte*(1.+0.5*gm1*xmach_new**2)**(-gamma/(gamma-1.))
      temperature = tte*(pressure/pte)**(gm1/gamma)
      qi0(j,k,1,1) = gamma*pressure/temperature
      qi0(j,k,2,1) = u_new*cos(alpe)*cos(betae)
      qi0(j,k,3,1) = -u_new*sin(betae)
      qi0(j,k,4,1) = u_new*sin(alpe)*cos(betae)
      qi0(j,k,5,1) = pressure
      qi0(j,k,1,2) = qi0(j,k,1,1)
      qi0(j,k,2,2) = qi0(j,k,2,1)
      qi0(j,k,3,2) = qi0(j,k,3,1)
      qi0(j,k,4,2) = qi0(j,k,4,1)
      qi0(j,k,5,2) = qi0(j,k,5,1)
c
      bci(j,k,1)   = 0.0
c
  840 continue
  340 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = 0.0
  591   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 501 k=ksta,kend1
        kk = k-ksta+1
        do 501 j=jsta,jend1
          jj = j-jsta+1
          ubar=-(qi0(j,k,2,1)*si(j,k,1,1)+qi0(j,k,3,1)*si(j,k,1,2)+
     +           qi0(j,k,4,1)*si(j,k,1,3))
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,4+l)
          if (real(ubar) .lt. 0.) then
             ti0(j,k,l,1) = t1
             ti0(j,k,l,2) = t1
          else
             ti0(j,k,l,1) = tursav(j,k,1,l)
             ti0(j,k,l,2) = ti0(j,k,l,1)
          end if
  501   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=idim boundary          nozzle total BCs                   type 2010
c******************************************************************************
c
      if (nface.eq.2) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,5)) .gt. -1.e10) itrflg1 = 1
c
      i=idim1
      do 350 k=ksta,kend1
      kk = k-ksta+1
      js = (k-ksta)*(jend-jsta)+1
      do 850 j=jsta,jend1
      jj = j-jsta+1
c
      pte   = bcdata(jj,kk,ip,1)
      tte   = bcdata(jj,kk,ip,2)
      alpe  = bcdata(jj,kk,ip,3)
      betae = bcdata(jj,kk,ip,4)
c
      pte=pte*p0
      alpe=alpe/radtodeg
      betae=betae/radtodeg
c
c   The following method is used by FUN3D:
      unormi = q(j,k,i,2)*si(j,k,idim,1) +
     +         q(j,k,i,3)*si(j,k,idim,2) +
     +         q(j,k,i,4)*si(j,k,idim,3)
      umag_sq = q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)*2
      umag = sqrt(umag_sq)
      ai = sqrt(gamma*q(j,k,i,5)/q(j,k,i,1))
      rminus = abs(unormi) - (2.*ai)/gm1
c     Note: xmach_local should be < 1 (not checking here)
      xmach_local = umag/ai
      hti = (q(j,k,i,5)/q(j,k,i,1))*(gamma/gm1)+0.5*umag_sq
      qa = 1.+2./gm1
      qb = 2.*rminus
      qc = 0.5*gm1*rminus**2 - gm1*hti
      a_radical = qb**2 - 4.*qa*qc
      if (real(a_radical) .gt. 0.) then
        a_plus = -qb/(2.*qa)+sqrt(a_radical)/(2.*qa)
        a_minus = -qb/(2.*qa)-sqrt(a_radical)/(2.*qa)
      else if (real(a_radical) .lt. 0.) then
        a_plus = ai
        a_minus = ai
      else
        a_plus = -qb/(2.*qa)
        a_minus = -qb/(2.*qa)
      end if
      if (real(a_plus) .gt. real(a_minus)) then
        a_new = a_plus
      else
        a_new = a_minus
      end if
      u_new = ((2.*a_new)/gm1) + rminus
      xmach_new = u_new/a_new
      pressure = pte*(1.+0.5*gm1*xmach_new**2)**(-gamma/(gamma-1.))
      temperature = tte*(pressure/pte)**(gm1/gamma)
      qi0(j,k,1,3) = gamma*pressure/temperature
      qi0(j,k,2,3) = u_new*cos(alpe)*cos(betae)
      qi0(j,k,3,3) = -u_new*sin(betae)
      qi0(j,k,4,3) = u_new*sin(alpe)*cos(betae)
      qi0(j,k,5,3) = pressure
      qi0(j,k,1,4) = qi0(j,k,1,3)
      qi0(j,k,2,4) = qi0(j,k,2,3)
      qi0(j,k,3,4) = qi0(j,k,3,3)
      qi0(j,k,4,4) = qi0(j,k,4,3)
      qi0(j,k,5,4) = qi0(j,k,5,3)
c
      bci(j,k,2)   = 0.0
c
  850 continue
  350 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
          vi0(j,k,1,3) = vist3d(j,k,idim-1)
          vi0(j,k,1,4) = 0.0
  691   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 601 k=ksta,kend1
        kk = k-ksta+1
        do 601 j=jsta,jend1
          jj = j-jsta+1
          ubar=qi0(j,k,2,3)*si(j,k,idim,1)+qi0(j,k,3,3)*si(j,k,idim,2)+
     +         qi0(j,k,4,3)*si(j,k,idim,3)
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,4+l)
          if (real(ubar) .lt. 0.) then
             ti0(j,k,l,3) = t1
             ti0(j,k,l,4) = t1
          else
             ti0(j,k,l,3) = tursav(j,k,idim-1,l)
             ti0(j,k,l,4) = ti0(j,k,l,3)
          end if
  601   continue
        enddo
      end if
      end if
c
      end if
c
      return
      end
