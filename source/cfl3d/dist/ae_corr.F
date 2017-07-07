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
      subroutine ae_corr(stm,stmi,xs,xxn,gforcn,gforcs,gforcnm,
     .                   gf0,lw,lw2,w,mgwk,maxbl,maxseg,
     .                   aesrfdat,nmds,maxaes,nt,mblk2nd,iseqr,
     .                   levelg,iadvance,nblock,icsi,icsf,jcsi,jcsf,
     .                   kcsi,kcsf,myid,nsegdfrm,idfrmseg,iaesurf,
     .                   perturb,aehist,ncycmax,maxsegdg,myhost,mycomm)
c
c     $Id$
c
c***********************************************************************
c     Purpose: update the aeroelastic generalized forces and update
c              the modal displacements and velocities via a corrector
c              step of the aeroelastic equations of motion.
c
c        Reference: Cunningham, H.J., Batina, J.T., and Bennett, R.M,
c                  "Modern Wing Flutter Analysis by Computational Fluid
c                   Dynamics Methods," J. Aircraft, Vol. 25, No. 10,
c                   October 1988, pp. 962-968.
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
#   ifdef DBLE_PRECSN
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_DOUBLE_COMPLEX
#      else
#        define MY_MPI_REAL MPI_DOUBLE_PRECISION
#      endif
#   else
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_COMPLEX
#      else
#        define MY_MPI_REAL MPI_REAL
#      endif
#   endif
      dimension istat(MPI_STATUS_SIZE)
#endif
c
      dimension w(mgwk),lw(65,maxbl),lw2(43,maxbl),levelg(maxbl)
      dimension gforcn(2*nmds,maxaes),gforcnm(2*nmds,maxaes),
     .          gforcs(2*nmds,maxaes),stm(2*nmds,2*nmds,maxaes),
     .          stmi(2*nmds,2*nmds,maxaes),xs(2*nmds,maxaes),
     .          xxn(2*nmds,maxaes),gf0(2*nmds,maxaes),
     .          aesrfdat(5,maxaes),fcoef(3),mblk2nd(maxbl),
     .          iadvance(maxbl),perturb(nmds,maxaes,4)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension idfrmseg(maxbl,maxsegdg),iaesurf(maxbl,maxsegdg),
     .          nsegdfrm(maxbl),aehist(ncycmax,3,nmds,maxaes)
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /elastic/ ndefrm,naesrf
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset = maxbl
      itag_frc   = 1
c
#endif
c
c     compute aeroelastic generalized forces
c
      do iaes=1,naesrf
         iskyhk = aesrfdat(1,iaes)
         grefl  = aesrfdat(2,iaes)
         uinf   = aesrfdat(3,iaes)
         qinf   = aesrfdat(4,iaes)
         nmodes = aesrfdat(5,iaes)
         do n=2,2*nmodes,2
            n2 = n/2
            cx = 0.
            cy = 0.
            cz = 0.
            do nbl=1,nblock
               cxb = 0.
               cyb = 0.
               czb = 0.
               if (iadvance(nbl).ge.0     .and.
     .             levelg(nbl).ge.lglobal .and.
     .             levelg(nbl).le.levelt(iseqr)) then
                  if (myid.eq.mblk2nd(nbl)) then
                     call lead(nbl,lw,lw2,maxbl)
                     call genforce(jdim,kdim,idim,w(lsk),w(lsj),
     .                             w(lsi),czb,cyb,cxb,w(lbcj),
     .                             w(lbck),w(lbci),w(lblk),nbl,
     .                             w(lqj0),w(lqk0),w(lqi0),maxbl,
     .                             maxseg,n2,w(lxmdj),w(lxmdk),
     .                             w(lxmdi),aesrfdat,nmds,maxaes,
     .                             maxsegdg,nsegdfrm,jcsi,jcsf,
     .                             kcsi,kcsf,icsi,icsf,idfrmseg,
     .                             iaes,iaesurf)
                  end if
#if defined DIST_MPI
c
                  nd_srce = mblk2nd(nbl)
                  mytag   = itag_frc + nbl
c
                  if (myid.eq.mblk2nd(nbl)) then
                     fcoef(1)  = cxb
                     fcoef(2)  = cyb
                     fcoef(3)  = czb
                     call MPI_Send(fcoef,3,MY_MPI_REAL,
     .                             myhost,mytag,mycomm,ierr)
                  end if
                  if (myid.eq.myhost) then
                     call MPI_Recv(fcoef,3,MY_MPI_REAL,
     .                             nd_srce,mytag,mycomm,istat,ierr)
                     cxb = fcoef(1)
                     cyb = fcoef(2)
                     czb = fcoef(3)
                  end if
#endif
                  if (myid.eq.myhost) then
                     cx  = cx + cxb
                     cy  = cy + cyb
                     cz  = cz + czb
                  end if
c
               end if
c
            end do
#if defined DIST_MPI
            if (myid.eq.myhost) then
               fcoef(1)  = cx
               fcoef(2)  = cy
               fcoef(3)  = cz
            end if
            call MPI_Bcast(fcoef,3,MY_MPI_REAL,myhost,
     .                  mycomm,ierr)
            cx = fcoef(1)
            cy = fcoef(2)
            cz = fcoef(3) 
#endif
            if (iskyhk.eq.1 .and. nt.eq.1) then
               gf0(n,iaes) = -qinf*grefl*grefl*(cx+cy+cz)+gf0(n,iaes)
            end if
            gforcs(n-1,iaes) = 0.
            gforcs(n,iaes)   = qinf*grefl*grefl*(cx+cy+cz)-gf0(n,iaes)
         end do
      end do
c
c     modal displacement and velocity correction
c
      do iaes=1,naesrf
         nmodes = aesrfdat(5,iaes)
         do n=1,2*nmodes
c           don't update if the modal time variation is specified
            moddfl = perturb((n+1)/2,iaes,1)
            if (moddfl .eq. 0) then
               xs(n,iaes) = 0.
               do j=1,2*nmodes
                  xs(n,iaes) = xs(n,iaes) + stm(n,j,iaes)*xxn(j,iaes) 
     .             + .5*stmi(n,j,iaes)*(gforcs(j,iaes)+gforcn(j,iaes))
               end do
            end if
         end do
      end do
c
c     update variables (set current values to old for next time step)
c
      do iaes=1,naesrf
         nmodes = aesrfdat(5,iaes)
         do n=1,2*nmodes
            gforcnm(n,iaes) = gforcn(n,iaes)
            gforcn(n,iaes)  = gforcs(n,iaes)
            xxn(n,iaes)     = xs(n,iaes)
         end do
      end do
c
c     store off the variables into the time-history array
c
      if (myid.eq.myhost) then
         do iaes=1,naesrf
            nmodes = aesrfdat(5,iaes)
            do n=1,nmodes
               aehist(ntt,1,n,iaes) = xs(2*n-1,iaes)
               aehist(ntt,2,n,iaes) = xs(2*n,iaes)
               aehist(ntt,3,n,iaes) = gforcs(2*n,iaes)
            end do
         end do
      end if
c
      return
      end
