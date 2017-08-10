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
      subroutine l2norm2(nbl,ntime,rmsl,irdq,jdim,kdim,idim,res,vol,
     .                   qc0,dqc0,q,blank)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute the L2-norm of the residuals, after subtracting 
c     out the contribution of the unsteady terms that where added in 
c     subroutine resadd.  
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension res(jdim,kdim,idim-1,5),q(jdim,kdim,idim,5),
     .          qc0(jdim,kdim,idim-1,5),dqc0(jdim,kdim,idim-1,5),
     .          vol(jdim,kdim,idim-1),blank(jdim,kdim,idim)
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
c
c     l2 norm of (residual - terms due to qc0 - terms due to dqc0)
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
      nplq  = min(idim1,999000/(jdim*kdim))
      npl   = nplq
      nt    = jdim*kdim
      rmsl  = 0.e0
      if (abs(ita) .eq. 1) then
        tfact=0.e0
      else
        tfact=0.5e0/dt
      end if
      tfacp1=tfact+1.e0/dt
c
      do 500 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      do 500 l=1,5
      do 600 ipl=1,npl
      ii = i+ipl-1
      do 700 j=1,jdim
      res(j,kdim,ii,l)  = 0.
      qc0(j,kdim,ii,l)  = q(j,kdim,ii,l)
      dqc0(j,kdim,ii,l) = 0.
  700 continue
c
cdir$ ivdep
      do 800 k=1,kdim1
      qc0(jdim,k,ii,l)  = q(jdim,k,ii,l)
      dqc0(jdim,k,ii,l) = 0.
  800 res(jdim,k,ii,l)  = 0.
  600 continue
  500 continue
c
      if (real(dt).gt.0) then
         factdqc0 = min(abs(ita)-1.,1.)
         if (real(factdqc0).lt.0.) factdqc0 = 0.
         factqc0 = 0.
         if (ncyc.gt.1) factqc0 = 1.
      else
         factqc0 = 0.
         factdqc0 = 0.
      end if
c
c     only density residual is monitored, so only modified 
c     density residual is needed
c
      do 900 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      n = npl*nt - jdim -1
cdir$ ivdep
      do 1000 izz=1,n
      resminus = res(izz,1,i,1)
     .         + factdqc0*tfact*vol(izz,1,i)*dqc0(izz,1,i,1)
     .         - factqc0*tfacp1*vol(izz,1,i)
     .         *(q(izz,1,i,1)-qc0(izz,1,i,1))
      rmsl     = rmsl+resminus*resminus*blank(izz,1,i)
 1000 continue
  900 continue
c
      return
      end
