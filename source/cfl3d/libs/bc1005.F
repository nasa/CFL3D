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
      subroutine bc1005(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,iuns,inormmom,
     .                  nou,bou,nbuf,ibufdim,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set inviscid surface boundary conditions
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /zero/iexp
c
c   Note: (10.**(-iexp) is machine zero)
      xminn=10.**(-iexp+1)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c            * * * * * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=1005 or 1006 *
c            * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary             inviscid surface                    bctype 1005
c******************************************************************************
      if (nface.eq.3) then
c 
      do 300 i=ista,iend1
c
      do 300 k=ksta,kend1
      pb           = q(1,k,i,5)
      qj0(k,i,1,1) = q(1,k,i,1)
c 
      contra       = q(1,k,i,2)*sj(1,k,i,1)
     .              +q(1,k,i,3)*sj(1,k,i,2) 
     .              +q(1,k,i,4)*sj(1,k,i,3)+sj(1,k,i,5)
c
      qj0(k,i,2,1) = q(1,k,i,2) - contra*sj(1,k,i,1)
      qj0(k,i,3,1) = q(1,k,i,3) - contra*sj(1,k,i,2)
      qj0(k,i,4,1) = q(1,k,i,4) - contra*sj(1,k,i,3)
c
      if (inormmom .eq. 1) then
        sixa=si(1,k,i,1)*si(1,k,i,4) + si(1,k,i+1,1)*si(1,k,i+1,4)
        siya=si(1,k,i,2)*si(1,k,i,4) + si(1,k,i+1,2)*si(1,k,i+1,4)
        siza=si(1,k,i,3)*si(1,k,i,4) + si(1,k,i+1,3)*si(1,k,i+1,4)
c
        skxa=sk(1,k,i,1)*sk(1,k,i,4) + sk(1,k+1,i,1)*sk(1,k+1,i,4)
        skya=sk(1,k,i,2)*sk(1,k,i,4) + sk(1,k+1,i,2)*sk(1,k+1,i,4)
        skza=sk(1,k,i,3)*sk(1,k,i,4) + sk(1,k+1,i,3)*sk(1,k+1,i,4)
c
        sjxa=2.*sj(1,k,i,1)*sj(1,k,i,4)
        sjya=2.*sj(1,k,i,2)*sj(1,k,i,4)
        sjza=2.*sj(1,k,i,3)*sj(1,k,i,4)
c
        ip=min(i+1,iend1)
        im=max(i-1,ista)
        factor=float(ip-im)
        factor=ccmaxcr(factor,1.)
        rxi = (sj(1,k,ip,1) - sj(1,k,im,1))/factor
        ryi = (sj(1,k,ip,2) - sj(1,k,im,2))/factor
        rzi = (sj(1,k,ip,3) - sj(1,k,im,3))/factor
        pi  = (q(1,k,ip,5)  - q(1,k,im,5)) /factor
c
        kp=min(k+1,kend1)
        km=max(k-1,ksta)
        factor=float(kp-km)
        factor=ccmaxcr(factor,1.)
        rxk = (sj(1,kp,i,1) - sj(1,km,i,1))/factor
        ryk = (sj(1,kp,i,2) - sj(1,km,i,2))/factor
        rzk = (sj(1,kp,i,3) - sj(1,km,i,3))/factor
        pk  = (q(1,kp,i,5)  - q(1,km,i,5)) /factor
c
        sii=sj(1,k,i,1)*sixa + sj(1,k,i,2)*siya + sj(1,k,i,3)*siza
        sjj=sj(1,k,i,1)*sjxa + sj(1,k,i,2)*sjya + sj(1,k,i,3)*sjza
        skk=sj(1,k,i,1)*skxa + sj(1,k,i,2)*skya + sj(1,k,i,3)*skza
c
        qi = qj0(k,i,2,1)*sixa + qj0(k,i,3,1)*siya + qj0(k,i,4,1)*siza
        qk = qj0(k,i,2,1)*skxa + qj0(k,i,3,1)*skya + qj0(k,i,4,1)*skza
c
        dp=  ( (qi*(qj0(k,i,2,1)*rxi + qj0(k,i,3,1)*ryi +
     +              qj0(k,i,4,1)*rzi)
     +         +qk*(qj0(k,i,2,1)*rxk + qj0(k,i,3,1)*ryk +
     +              qj0(k,i,4,1)*rzk))*qj0(k,i,1,1)
     +         -sii*pi - skk*pk) / sjj
        jp2= min(2,jdim1)
        pb = 1.125*pb -.125*q(jp2,k,i,5)- dp*.375
        pb = ccmax(pb,xminn)
c       now correct density for entropy and velocity for enthalpy
        sent         = q(1,k,i,5)/q(1,k,i,1)**gamma
        hent         = ( gamma*q(1,k,i,5)/q(1,k,i,1) )/gm1 +
     .      0.5*( q(1,k,i,2)**2 + q(1,k,i,3)**2 + q(1,k,i,4)**2 )
        rhob         = ( pb/sent )**(1./gamma)
        term2        = rhob/qj0(k,i,1,1)
        vmag         = 2.*( hent - (gamma*pb/rhob)/gm1 )
        term1        = sqrt( ccabs(vmag) /
     +         (qj0(k,i,2,1)**2+qj0(k,i,3,1)**2+qj0(k,i,4,1)**2))
        uvelb = qj0(k,i,2,1)*term1
        vvelb = qj0(k,i,3,1)*term1
        wvelb = qj0(k,i,4,1)*term1
c
        qj0(k,i,1,1) = rhob
        qj0(k,i,2,1) = uvelb
        qj0(k,i,3,1) = vvelb
        qj0(k,i,4,1) = wvelb
      end if
c
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
c
  300 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,2) = vist3d(1,k,i)
  191   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 101 i=ista,iend1
        do 101 k=ksta,kend1
          tj0(k,i,l,1) = tursav(1,k,i,l)
          tj0(k,i,l,2) = tursav(1,k,i,l)
  101   continue
        enddo
      end if
      end if
c
      end if
c 
c******************************************************************************
c      j=jdim boundary          inviscid surface                    bctype 1005
c******************************************************************************
      if (nface.eq.4) then
c 
      do 600 i=ista,iend1
c
      do 600 k=ksta,kend1
      pb           = q(jdim1,k,i,5)
      qj0(k,i,1,3) = q(jdim1,k,i,1)
c 
      contra       = q(jdim1,k,i,2)*sj(jdim,k,i,1)
     .              +q(jdim1,k,i,3)*sj(jdim,k,i,2) 
     .              +q(jdim1,k,i,4)*sj(jdim,k,i,3)+sj(jdim,k,i,5)
c
      qj0(k,i,2,3) = q(jdim1,k,i,2) - contra*sj(jdim,k,i,1)
      qj0(k,i,3,3) = q(jdim1,k,i,3) - contra*sj(jdim,k,i,2)
      qj0(k,i,4,3) = q(jdim1,k,i,4) - contra*sj(jdim,k,i,3)
c
      if (inormmom .eq. 1) then
        sixa=si(jdim1,k,i,1)*si(jdim1,k,i,4) +
     +       si(jdim1,k,i+1,1)*si(jdim1,k,i+1,4)
        siya=si(jdim1,k,i,2)*si(jdim1,k,i,4) +
     +       si(jdim1,k,i+1,2)*si(jdim1,k,i+1,4)
        siza=si(jdim1,k,i,3)*si(jdim1,k,i,4) +
     +       si(jdim1,k,i+1,3)*si(jdim1,k,i+1,4)
c
        skxa=sk(jdim1,k,i,1)*sk(jdim1,k,i,4) +
     +       sk(jdim1,k+1,i,1)*sk(jdim1,k+1,i,4)
        skya=sk(jdim1,k,i,2)*sk(jdim1,k,i,4) +
     +       sk(jdim1,k+1,i,2)*sk(jdim1,k+1,i,4)
        skza=sk(jdim1,k,i,3)*sk(jdim1,k,i,4) +
     +       sk(jdim1,k+1,i,3)*sk(jdim1,k+1,i,4)
c
        sjxa=2.*sj(jdim,k,i,1)*sj(jdim,k,i,4)
        sjya=2.*sj(jdim,k,i,2)*sj(jdim,k,i,4)
        sjza=2.*sj(jdim,k,i,3)*sj(jdim,k,i,4)
c
        ip=min(i+1,iend1)
        im=max(i-1,ista)
        factor=float(ip-im)
        factor=ccmaxcr(factor,1.)
        rxi = (sj(jdim,k,ip,1) - sj(jdim,k,im,1))/factor
        ryi = (sj(jdim,k,ip,2) - sj(jdim,k,im,2))/factor
        rzi = (sj(jdim,k,ip,3) - sj(jdim,k,im,3))/factor
        pi  = (q(jdim1,k,ip,5)  - q(jdim1,k,im,5)) /factor
c
        kp=min(k+1,kend1)
        km=max(k-1,ksta)
        factor=float(kp-km)
        factor=ccmaxcr(factor,1.)
        rxk = (sj(jdim,kp,i,1) - sj(jdim,km,i,1))/factor
        ryk = (sj(jdim,kp,i,2) - sj(jdim,km,i,2))/factor
        rzk = (sj(jdim,kp,i,3) - sj(jdim,km,i,3))/factor
        pk  = (q(jdim1,kp,i,5)  - q(jdim1,km,i,5)) /factor
c
        sii=sj(jdim,k,i,1)*sixa + sj(jdim,k,i,2)*siya +
     +      sj(jdim,k,i,3)*siza
        sjj=sj(jdim,k,i,1)*sjxa + sj(jdim,k,i,2)*sjya +
     +      sj(jdim,k,i,3)*sjza
        skk=sj(jdim,k,i,1)*skxa + sj(jdim,k,i,2)*skya +
     +      sj(jdim,k,i,3)*skza
c
        qi = qj0(k,i,2,3)*sixa + qj0(k,i,3,3)*siya + qj0(k,i,4,3)*siza
        qk = qj0(k,i,2,3)*skxa + qj0(k,i,3,3)*skya + qj0(k,i,4,3)*skza
c
        dp= -( (qi*(qj0(k,i,2,3)*rxi + qj0(k,i,3,3)*ryi +
     +              qj0(k,i,4,3)*rzi)
     +         +qk*(qj0(k,i,2,3)*rxk + qj0(k,i,3,3)*ryk +
     +              qj0(k,i,4,3)*rzk))*qj0(k,i,1,3)
     +         -sii*pi - skk*pk) / sjj
        jm2= max(jdim-2,1)
        pb = 1.125*pb -.125*q(jm2,k,i,5)- dp*.375
        pb = ccmax(pb,xminn)
c       now correct density for entropy and velocity for enthalpy
        sent         = q(jdim1,k,i,5)/q(jdim1,k,i,1)**gamma
        hent         = ( gamma*q(jdim1,k,i,5)/q(jdim1,k,i,1) )/gm1 +
     .      0.5*( q(jdim1,k,i,2)**2 + q(jdim1,k,i,3)**2 +
     +            q(jdim1,k,i,4)**2 )
        rhob         = ( pb/sent )**(1./gamma)
        term2        = rhob/qj0(k,i,1,3)
        vmag         = 2.*( hent - (gamma*pb/rhob)/gm1 )
        term1        = sqrt( ccabs(vmag) /
     +         (qj0(k,i,2,3)**2+qj0(k,i,3,3)**2+qj0(k,i,4,3)**2))
        uvelb = qj0(k,i,2,3)*term1
        vvelb = qj0(k,i,3,3)*term1
        wvelb = qj0(k,i,4,3)*term1
c
        qj0(k,i,1,3) = rhob
        qj0(k,i,2,3) = uvelb
        qj0(k,i,3,3) = vvelb
        qj0(k,i,4,3) = wvelb
      end if
c
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
c
  600 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
          vj0(k,i,1,3) = vist3d(jdim1,k,i)
          vj0(k,i,1,4) = vist3d(jdim1,k,i)
  291   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 201 i=ista,iend1
        do 201 k=ksta,kend1
          tj0(k,i,l,3) = tursav(jdim1,k,i,l)
          tj0(k,i,l,4) = tursav(jdim1,k,i,l)
  201   continue
        enddo
      end if
      end if
c
      end if
c 
c******************************************************************************
c      k=1 boundary             inviscid surface                    bctype 1005
c******************************************************************************
      if (nface.eq.5) then
c
      do 900 i=ista,iend1
c
      do 900 j=jsta,jend1
      pb           = q(j,1,i,5)
      qk0(j,i,1,1) = q(j,1,i,1)
c 
      contra       = q(j,1,i,2)*sk(j,1,i,1)
     .              +q(j,1,i,3)*sk(j,1,i,2) 
     .              +q(j,1,i,4)*sk(j,1,i,3)+sk(j,1,i,5)
c
      qk0(j,i,2,1) = q(j,1,i,2) - contra*sk(j,1,i,1)
      qk0(j,i,3,1) = q(j,1,i,3) - contra*sk(j,1,i,2)
      qk0(j,i,4,1) = q(j,1,i,4) - contra*sk(j,1,i,3)
c
      if (inormmom .eq. 1) then
        sixa=si(j,1,i,1)*si(j,1,i,4) + si(j,1,i+1,1)*si(j,1,i+1,4)
        siya=si(j,1,i,2)*si(j,1,i,4) + si(j,1,i+1,2)*si(j,1,i+1,4)
        siza=si(j,1,i,3)*si(j,1,i,4) + si(j,1,i+1,3)*si(j,1,i+1,4)
c
        sjxa=sj(j,1,i,1)*sj(j,1,i,4) + sj(j+1,1,i,1)*sj(j+1,1,i,4)
        sjya=sj(j,1,i,2)*sj(j,1,i,4) + sj(j+1,1,i,2)*sj(j+1,1,i,4)
        sjza=sj(j,1,i,3)*sj(j,1,i,4) + sj(j+1,1,i,3)*sj(j+1,1,i,4)
c
        skxa=2.*sk(j,1,i,1)*sk(j,1,i,4)
        skya=2.*sk(j,1,i,2)*sk(j,1,i,4)
        skza=2.*sk(j,1,i,3)*sk(j,1,i,4)
c
        ip=min(i+1,iend1)
        im=max(i-1,ista)
        factor=float(ip-im)
        factor=ccmaxcr(factor,1.)
        rxi = (sk(j,1,ip,1) - sk(j,1,im,1))/factor
        ryi = (sk(j,1,ip,2) - sk(j,1,im,2))/factor
        rzi = (sk(j,1,ip,3) - sk(j,1,im,3))/factor
        pi  = (q(j,1,ip,5)  - q(j,1,im,5)) /factor
c
        jp=min(j+1,jend1)
        jm=max(j-1,jsta)
        factor=float(jp-jm)
        factor=ccmaxcr(factor,1.)
        rxj = (sk(jp,1,i,1) - sk(jm,1,i,1))/factor
        ryj = (sk(jp,1,i,2) - sk(jm,1,i,2))/factor
        rzj = (sk(jp,1,i,3) - sk(jm,1,i,3))/factor
        pj  = (q(jp,1,i,5)  - q(jm,1,i,5)) /factor
c
        sii=sk(j,1,i,1)*sixa + sk(j,1,i,2)*siya + sk(j,1,i,3)*siza
        sjj=sk(j,1,i,1)*sjxa + sk(j,1,i,2)*sjya + sk(j,1,i,3)*sjza
        skk=sk(j,1,i,1)*skxa + sk(j,1,i,2)*skya + sk(j,1,i,3)*skza
c
        qi = qk0(j,i,2,1)*sixa + qk0(j,i,3,1)*siya + qk0(j,i,4,1)*siza
        qj = qk0(j,i,2,1)*sjxa + qk0(j,i,3,1)*sjya + qk0(j,i,4,1)*sjza
c
        dp=  ( ( qi*(qk0(j,i,2,1)*rxi + qk0(j,i,3,1)*ryi +
     +               qk0(j,i,4,1)*rzi)
     +          +qj*(qk0(j,i,2,1)*rxj + qk0(j,i,3,1)*ryj +
     +               qk0(j,i,4,1)*rzj))*qk0(j,i,1,1)
     +         -sii*pi - sjj*pj ) / skk
        kp2= min(2,kdim1)
        pb = 1.125*pb -.125*q(j,kp2,i,5)- dp*.375
        pb = ccmax(pb,xminn)
c       now correct density for entropy and velocity for enthalpy
        sent         = q(j,1,i,5)/q(j,1,i,1)**gamma
        hent         = ( gamma*q(j,1,i,5)/q(j,1,i,1) )/gm1 +
     .      0.5*( q(j,1,i,2)**2 + q(j,1,i,3)**2 + q(j,1,i,4)**2 )
        rhob         = ( pb/sent )**(1./gamma)
        term2        = rhob/qk0(j,i,1,1)
        vmag         = 2.*( hent - (gamma*pb/rhob)/gm1 )
        term1        = sqrt( ccabs(vmag) /
     +         (qk0(j,i,2,1)**2+qk0(j,i,3,1)**2+qk0(j,i,4,1)**2))
        uvelb = qk0(j,i,2,1)*term1
        vvelb = qk0(j,i,3,1)*term1
        wvelb = qk0(j,i,4,1)*term1
c
        qk0(j,i,1,1) = rhob
        qk0(j,i,2,1) = uvelb
        qk0(j,i,3,1) = vvelb
        qk0(j,i,4,1) = wvelb
      end if
c
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
c
  900 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = vist3d(j,1,i)
  391   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 301 i=ista,iend1
        do 301 j=jsta,jend1
          tk0(j,i,l,1) = tursav(j,1,i,l)
          tk0(j,i,l,2) = tursav(j,1,i,l)
  301   continue
        enddo
      end if
      end if
c
      end if
c 
c******************************************************************************
c      k=kdim boundary          inviscid surface                    bctype 1005
c******************************************************************************
      if (nface.eq.6) then
c 
      do 1200 i=ista,iend1
c
      do 1200 j=jsta,jend1
      pb           = q(j,kdim1,i,5)
      qk0(j,i,1,3) = q(j,kdim1,i,1)
c 
      contra       = q(j,kdim1,i,2)*sk(j,kdim,i,1)
     .              +q(j,kdim1,i,3)*sk(j,kdim,i,2) 
     .              +q(j,kdim1,i,4)*sk(j,kdim,i,3)+sk(j,kdim,i,5)
c
      qk0(j,i,2,3) = q(j,kdim1,i,2) - contra*sk(j,kdim,i,1)
      qk0(j,i,3,3) = q(j,kdim1,i,3) - contra*sk(j,kdim,i,2)
      qk0(j,i,4,3) = q(j,kdim1,i,4) - contra*sk(j,kdim,i,3)
c
      if (inormmom .eq. 1) then
        sixa=si(j,kdim1,i,1)*si(j,kdim1,i,4) +
     +       si(j,kdim1,i+1,1)*si(j,kdim1,i+1,4)
        siya=si(j,kdim1,i,2)*si(j,kdim1,i,4) +
     +       si(j,kdim1,i+1,2)*si(j,kdim1,i+1,4)
        siza=si(j,kdim1,i,3)*si(j,kdim1,i,4) +
     +       si(j,kdim1,i+1,3)*si(j,kdim1,i+1,4)
c
        sjxa=sj(j,kdim1,i,1)*sj(j,kdim1,i,4) +
     +       sj(j+1,kdim1,i,1)*sj(j+1,kdim1,i,4)
        sjya=sj(j,kdim1,i,2)*sj(j,kdim1,i,4) +
     +       sj(j+1,kdim1,i,2)*sj(j+1,kdim1,i,4)
        sjza=sj(j,kdim1,i,3)*sj(j,kdim1,i,4) +
     +       sj(j+1,kdim1,i,3)*sj(j+1,kdim1,i,4)
c
        skxa=2.*sk(j,kdim,i,1)*sk(j,kdim,i,4)
        skya=2.*sk(j,kdim,i,2)*sk(j,kdim,i,4)
        skza=2.*sk(j,kdim,i,3)*sk(j,kdim,i,4)
c
        ip=min(i+1,iend1)
        im=max(i-1,ista)
        factor=float(ip-im)
        factor=ccmaxcr(factor,1.)
        rxi = (sk(j,kdim,ip,1) - sk(j,kdim,im,1))/factor
        ryi = (sk(j,kdim,ip,2) - sk(j,kdim,im,2))/factor
        rzi = (sk(j,kdim,ip,3) - sk(j,kdim,im,3))/factor
        pi  = (q(j,kdim1,ip,5)  - q(j,kdim1,im,5)) /factor
c
        jp=min(j+1,jend1)
        jm=max(j-1,jsta)
        factor=float(jp-jm)
        factor=ccmaxcr(factor,1.)
        rxj = (sk(jp,kdim,i,1) - sk(jm,kdim,i,1))/factor
        ryj = (sk(jp,kdim,i,2) - sk(jm,kdim,i,2))/factor
        rzj = (sk(jp,kdim,i,3) - sk(jm,kdim,i,3))/factor
        pj  = (q(jp,kdim1,i,5)  - q(jm,kdim1,i,5)) /factor
c
        sii=sk(j,kdim,i,1)*sixa + sk(j,kdim,i,2)*siya +
     +      sk(j,kdim,i,3)*siza
        sjj=sk(j,kdim,i,1)*sjxa + sk(j,kdim,i,2)*sjya +
     +      sk(j,kdim,i,3)*sjza
        skk=sk(j,kdim,i,1)*skxa + sk(j,kdim,i,2)*skya +
     +      sk(j,kdim,i,3)*skza
c
        qi = qk0(j,i,2,3)*sixa + qk0(j,i,3,3)*siya + qk0(j,i,4,3)*siza
        qj = qk0(j,i,2,3)*sjxa + qk0(j,i,3,3)*sjya + qk0(j,i,4,3)*sjza
c
        dp= -( ( qi*(qk0(j,i,2,3)*rxi + qk0(j,i,3,3)*ryi +
     +               qk0(j,i,4,3)*rzi)
     +          +qj*(qk0(j,i,2,3)*rxj + qk0(j,i,3,3)*ryj +
     +               qk0(j,i,4,3)*rzj))*qk0(j,i,1,3)
     +         -sii*pi - sjj*pj ) / skk
        km2= max(kdim-2,1)
        pb = 1.125*pb -.125*q(j,km2,i,5)- dp*.375
        pb = ccmax(pb,xminn)
c       now correct density for entropy and velocity for enthalpy
        sent         = q(j,kdim1,i,5)/q(j,kdim1,i,1)**gamma
        hent         = ( gamma*q(j,kdim1,i,5)/q(j,kdim1,i,1) )/gm1 +
     .      0.5*( q(j,kdim1,i,2)**2 + q(j,kdim1,i,3)**2 +
     +            q(j,kdim1,i,4)**2 )
        rhob         = ( pb/sent )**(1./gamma)
        term2        = rhob/qk0(j,i,1,3)
        vmag         = 2.*( hent - (gamma*pb/rhob)/gm1 )
        term1        = sqrt( ccabs(vmag) /
     +         (qk0(j,i,2,3)**2+qk0(j,i,3,3)**2+qk0(j,i,4,3)**2))
        uvelb = qk0(j,i,2,3)*term1
        vvelb = qk0(j,i,3,3)*term1
        wvelb = qk0(j,i,4,3)*term1
c
        qk0(j,i,1,3) = rhob
        qk0(j,i,2,3) = uvelb
        qk0(j,i,3,3) = vvelb
        qk0(j,i,4,3) = wvelb
      end if
c
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
c
 1200 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
          vk0(j,i,1,3) = vist3d(j,kdim1,i)
          vk0(j,i,1,4) = vist3d(j,kdim1,i)
  491   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 401 i=ista,iend1
        do 401 j=jsta,jend1
          tk0(j,i,l,3) = tursav(j,kdim1,i,l)
          tk0(j,i,l,4) = tursav(j,kdim1,i,l)
  401   continue
        enddo
      end if
      end if
c
      end if
c 
c******************************************************************************
c      i=1 boundary             inviscid surface                    bctype 1005
c******************************************************************************
      if (nface.eq.1) then
c 
      do 1500 k=ksta,kend1
c
      do 1500 j=jsta,jend1
      pb           = q(j,k,1,5)
      qi0(j,k,1,1) = q(j,k,1,1)
c 
      contra       = q(j,k,1,2)*si(j,k,1,1)
     .              +q(j,k,1,3)*si(j,k,1,2) 
     .              +q(j,k,1,4)*si(j,k,1,3)+si(j,k,1,5)
c
      qi0(j,k,2,1) = q(j,k,1,2) - contra*si(j,k,1,1)
      qi0(j,k,3,1) = q(j,k,1,3) - contra*si(j,k,1,2)
      qi0(j,k,4,1) = q(j,k,1,4) - contra*si(j,k,1,3)
c
      if (inormmom .eq. 1) then
        skxa=sk(j,k,1,1)*sk(j,k,1,4) + sk(j,k+1,1,1)*sk(j,k+1,1,4)
        skya=sk(j,k,1,2)*sk(j,k,1,4) + sk(j,k+1,1,2)*sk(j,k+1,1,4)
        skza=sk(j,k,1,3)*sk(j,k,1,4) + sk(j,k+1,1,3)*sk(j,k+1,1,4)
c
        sjxa=sj(j,k,1,1)*sj(j,k,1,4) + sj(j+1,k,1,1)*sj(j+1,k,1,4)
        sjya=sj(j,k,1,2)*sj(j,k,1,4) + sj(j+1,k,1,2)*sj(j+1,k,1,4)
        sjza=sj(j,k,1,3)*sj(j,k,1,4) + sj(j+1,k,1,3)*sj(j+1,k,1,4)
c
        sixa=2.*si(j,k,1,1)*si(j,k,1,4)
        siya=2.*si(j,k,1,2)*si(j,k,1,4)
        siza=2.*si(j,k,1,3)*si(j,k,1,4)
c
        kp=min(k+1,kend1)
        km=max(k-1,ksta)
        factor=float(kp-km)
        factor=ccmaxcr(factor,1.)
        rxk = (si(j,kp,1,1) - si(j,km,1,1))/factor
        ryk = (si(j,kp,1,2) - si(j,km,1,2))/factor
        rzk = (si(j,kp,1,3) - si(j,km,1,3))/factor
        pk  = (q(j,kp,1,5)  - q(j,km,1,5)) /factor
c
        jp=min(j+1,jend1)
        jm=max(j-1,jsta)
        factor=float(jp-jm)
        factor=ccmaxcr(factor,1.)
        rxj = (si(jp,k,1,1) - si(jm,k,1,1))/factor
        ryj = (si(jp,k,1,2) - si(jm,k,1,2))/factor
        rzj = (si(jp,k,1,3) - si(jm,k,1,3))/factor
        pj  = (q(jp,k,1,5)  - q(jm,k,1,5)) /factor
c
        sii=si(j,k,1,1)*sixa + si(j,k,1,2)*siya + si(j,k,1,3)*siza
        sjj=si(j,k,1,1)*sjxa + si(j,k,1,2)*sjya + si(j,k,1,3)*sjza
        skk=si(j,k,1,1)*skxa + si(j,k,1,2)*skya + si(j,k,1,3)*skza
c
        qk = qi0(j,k,2,1)*skxa + qi0(j,k,3,1)*skya + qi0(j,k,4,1)*skza
        qj = qi0(j,k,2,1)*sjxa + qi0(j,k,3,1)*sjya + qi0(j,k,4,1)*sjza
c
        dp=  ( (qk*(qi0(j,k,2,1)*rxk + qi0(j,k,3,1)*ryk +
     +              qi0(j,k,4,1)*rzk)
     +         +qj*(qi0(j,k,2,1)*rxj + qi0(j,k,3,1)*ryj +
     +              qi0(j,k,4,1)*rzj))*qi0(j,k,1,1)
     +         -skk*pk - sjj*pj) / sii
        ip2= min(2,idim1)
        pb = 1.125*pb -.125*q(j,k,ip2,5)- dp*.375
        pb = ccmax(pb,xminn)
c       now correct density for entropy and velocity for enthalpy
        sent         = q(j,k,1,5)/q(j,k,1,1)**gamma
        hent         = ( gamma*q(j,k,1,5)/q(j,k,1,1) )/gm1 +
     .      0.5*( q(j,k,1,2)**2 + q(j,k,1,3)**2 + q(j,k,1,4)**2 )
        rhob         = ( pb/sent )**(1./gamma)
        term2        = rhob/qi0(j,k,1,1)
        vmag         = 2.*( hent - (gamma*pb/rhob)/gm1 )
        term1        = sqrt( ccabs(vmag) /
     +         (qi0(j,k,2,1)**2+qi0(j,k,3,1)**2+qi0(j,k,4,1)**2))
        uvelb = qi0(j,k,2,1)*term1
        vvelb = qi0(j,k,3,1)*term1
        wvelb = qi0(j,k,4,1)*term1
c
        qi0(j,k,1,1) = rhob
        qi0(j,k,2,1) = uvelb
        qi0(j,k,3,1) = vvelb
        qi0(j,k,4,1) = wvelb
      end if
c
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
      z1  =   2.0 +1.5*f23
      z2  =       -0.5*f23   
      z3  = -(2.0 +    f23)
c
      qi0(j,k,1,2) = z1*q(j,k,1,1) + z2*q(j,k,i2,1) + z3*qi0(j,k,1,1)
      qi0(j,k,2,2) = z1*q(j,k,1,2) + z2*q(j,k,i2,2) + z3*qi0(j,k,2,1)
      qi0(j,k,3,2) = z1*q(j,k,1,3) + z2*q(j,k,i2,3) + z3*qi0(j,k,3,1)
      qi0(j,k,4,2) = z1*q(j,k,1,4) + z2*q(j,k,i2,4) + z3*qi0(j,k,4,1)
      qi0(j,k,5,2) = z1*q(j,k,1,5) + z2*q(j,k,i2,5) + z3*qi0(j,k,5,1)
 1500 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = vist3d(j,k,1)
  591   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 501 k=ksta,kend1
        do 501 j=jsta,jend1
          ti0(j,k,l,1) = tursav(j,k,1,l)
          ti0(j,k,l,2) = tursav(j,k,1,l)
  501   continue
        enddo
      end if
      end if
c
      end if
c 
c******************************************************************************
c      i=idim boundary          inviscid surface                    bctype 1005
c******************************************************************************
      if (nface.eq.2) then
c 
      do 1800 k=ksta,kend1
c
      do 1800 j=jsta,jend1
      pb           = q(j,k,idim1,5)
      qi0(j,k,1,3) = q(j,k,idim1,1)
c
      contra       = q(j,k,idim1,2)*si(j,k,idim,1)
     .              +q(j,k,idim1,3)*si(j,k,idim,2) 
     .              +q(j,k,idim1,4)*si(j,k,idim,3)+si(j,k,idim,5)
c
      qi0(j,k,2,3) = q(j,k,idim1,2) - contra*si(j,k,idim,1)
      qi0(j,k,3,3) = q(j,k,idim1,3) - contra*si(j,k,idim,2)
      qi0(j,k,4,3) = q(j,k,idim1,4) - contra*si(j,k,idim,3)
c
      if (inormmom .eq. 1) then
        skxa=sk(j,k,idim1,1)*sk(j,k,idim1,4) +
     +       sk(j,k+1,idim1,1)*sk(j,k+1,idim1,4)
        skya=sk(j,k,idim1,2)*sk(j,k,idim1,4) +
     +       sk(j,k+1,idim1,2)*sk(j,k+1,idim1,4)
        skza=sk(j,k,idim1,3)*sk(j,k,idim1,4) +
     +       sk(j,k+1,idim1,3)*sk(j,k+1,idim1,4)
c
        sjxa=sj(j,k,idim1,1)*sj(j,k,idim1,4) +
     +       sj(j+1,k,idim1,1)*sj(j+1,k,idim1,4)
        sjya=sj(j,k,idim1,2)*sj(j,k,idim1,4) +
     +       sj(j+1,k,idim1,2)*sj(j+1,k,idim1,4)
        sjza=sj(j,k,idim1,3)*sj(j,k,idim1,4) +
     +       sj(j+1,k,idim1,3)*sj(j+1,k,idim1,4)
c
        sixa=2.*si(j,k,idim,1)*si(j,k,idim,4)
        siya=2.*si(j,k,idim,2)*si(j,k,idim,4)
        siza=2.*si(j,k,idim,3)*si(j,k,idim,4)
c
        kp=min(k+1,kend1)
        km=max(k-1,ksta)
        factor=float(kp-km)
        factor=ccmaxcr(factor,1.)
        rxk = (si(j,kp,idim,1) - si(j,km,idim,1))/factor
        ryk = (si(j,kp,idim,2) - si(j,km,idim,2))/factor
        rzk = (si(j,kp,idim,3) - si(j,km,idim,3))/factor
        pk  = (q(j,kp,idim1,5)  - q(j,km,idim1,5)) /factor
c
        jp=min(j+1,jend1)
        jm=max(j-1,jsta)
        factor=float(jp-jm)
        factor=ccmaxcr(factor,1.)
        rxj = (si(jp,k,idim,1) - si(jm,k,idim,1))/factor
        ryj = (si(jp,k,idim,2) - si(jm,k,idim,2))/factor
        rzj = (si(jp,k,idim,3) - si(jm,k,idim,3))/factor
        pj  = (q(jp,k,idim1,5)  - q(jm,k,idim1,5)) /factor
c
        sii=si(j,k,idim,1)*sixa + si(j,k,idim,2)*siya +
     +      si(j,k,idim,3)*siza
        sjj=si(j,k,idim,1)*sjxa + si(j,k,idim,2)*sjya +
     +      si(j,k,idim,3)*sjza
        skk=si(j,k,idim,1)*skxa + si(j,k,idim,2)*skya +
     +      si(j,k,idim,3)*skza
c
        qk = qi0(j,k,2,3)*skxa + qi0(j,k,3,3)*skya + qi0(j,k,4,3)*skza
        qj = qi0(j,k,2,3)*sjxa + qi0(j,k,3,3)*sjya + qi0(j,k,4,3)*sjza
c
        dp= -( (qk*(qi0(j,k,2,3)*rxk + qi0(j,k,3,3)*ryk +
     +              qi0(j,k,4,3)*rzk)
     +         +qj*(qi0(j,k,2,3)*rxj + qi0(j,k,3,3)*ryj +
     +              qi0(j,k,4,3)*rzj))*qi0(j,k,1,3)
     +         -skk*pk - sjj*pj) / sii
        im2= max(idim-2,1)
        pb = 1.125*pb -.125*q(j,k,im2,5)- dp*.375
        pb = ccmax(pb,xminn)
c       now correct density for entropy and velocity for enthalpy
        sent         = q(j,k,idim1,5)/q(j,k,idim1,1)**gamma
        hent         = ( gamma*q(j,k,idim1,5)/q(j,k,idim1,1) )/gm1 +
     .      0.5*( q(j,k,idim1,2)**2 + q(j,k,idim1,3)**2 +
     +      q(j,k,idim1,4)**2 )
        rhob         = ( pb/sent )**(1./gamma)
        term2        = rhob/qi0(j,k,1,3)
        vmag         = 2.*( hent - (gamma*pb/rhob)/gm1 )
        term1        = sqrt( ccabs(vmag) /
     +         (qi0(j,k,2,3)**2+qi0(j,k,3,3)**2+qi0(j,k,4,3)**2))
        uvelb = qi0(j,k,2,3)*term1
        vvelb = qi0(j,k,3,3)*term1
        wvelb = qi0(j,k,4,3)*term1
c
        qi0(j,k,1,3) = rhob
        qi0(j,k,2,3) = uvelb
        qi0(j,k,3,3) = vvelb
        qi0(j,k,4,3) = wvelb
      end if
c
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
      z1  =  -2.0 -1.5*f23
      z2  =       +0.5*f23   
      z3  = +(2.0 +    f23)
c
      qi0(j,k,1,4) = z1*q(j,k,idim1,1)+z2*q(j,k,i2,1)+z3*qi0(j,k,1,3)
      qi0(j,k,2,4) = z1*q(j,k,idim1,2)+z2*q(j,k,i2,2)+z3*qi0(j,k,2,3)
      qi0(j,k,3,4) = z1*q(j,k,idim1,3)+z2*q(j,k,i2,3)+z3*qi0(j,k,3,3)
      qi0(j,k,4,4) = z1*q(j,k,idim1,4)+z2*q(j,k,i2,4)+z3*qi0(j,k,4,3)
      qi0(j,k,5,4) = z1*q(j,k,idim1,5)+z2*q(j,k,i2,5)+z3*qi0(j,k,5,3)
 1800 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
          vi0(j,k,1,3) = vist3d(j,k,idim1)
          vi0(j,k,1,4) = vist3d(j,k,idim1)
  691   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 601 k=ksta,kend1
        do 601 j=jsta,jend1
          ti0(j,k,l,3) = tursav(j,k,idim1,l)
          ti0(j,k,l,4) = tursav(j,k,idim1,l)
  601   continue
        enddo
      end if
      end if
      end if
c
      return
      end
