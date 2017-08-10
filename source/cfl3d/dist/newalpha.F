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
      subroutine newalpha(ncycmax,rms,clw,myid,myhost,mycomm)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  determine a new estimate for alpha that will yield
c     the specified value of Cl
c
c     original coding courtesy Neal Frink, cfl3d mods by R. Biedron
c     more mods by C. Rumsey, 3/01
c
c     recoded by S. Allmaras (5/03) to implmement TLNS3D Cl matching
c     algorithm; at each iteration perform:
c     (1) relaxation of alpha based on error in Cl
c     (2) instantaneous velocity rotation for all points in flowfield;
c         the rotation angle is the change in alpha
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
c
#endif
      dimension rms(ncycmax),clw(ncycmax)
c
      common /alphait/ ialphit,cltarg,rlxalph,dalim,dalpha,icycupdt
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /conversion/ radtodeg
c
c     cltarg...target cl
c     rlxalph..relaxation parameter for updating alpha (default = 1)
c     dalim....limit on magnitude of alpha change (default = 0.2 deg)
c     dalpha...alpha change
c     ioalph...if > 0, output alpha history to cfl3d.alpha
      ioalph = 1
c
c     update alpha based on error in cl; 
      if (myid.eq.myhost) then
         rlx    = rlxalph / real(radtodeg)
         dalpha = -rlx*(clw(ntt) - cltarg)
         if (abs(real(dalpha)) .gt. real(dalim)/real(radtodeg)) then
           sgn    = ccsignrc(1., dalpha)
           dalpha = sgn*(dalim/radtodeg)
         end if
         alpha     = alpha + dalpha
c
c        output to cfl3d.alpha
         if (ioalph .gt. 0) then
            write(27,'(3x,i6,4e16.8)') ntt,log10(real(rms(ntt))),
     .                                 real(clw(ntt)),
     .                                 real(alpha)*real(radtodeg),
     .                                 real(dalpha)*real(radtodeg)
         end if
      end if
c
#if defined DIST_MPI
      call mpi_Bcast(alpha,1,MY_MPI_REAL,myhost,mycomm,ierr)
      call mpi_Bcast(dalpha,1,MY_MPI_REAL,myhost,mycomm,ierr)
#endif
c     reset freestream for new alpha
      call init_mast
c
      return
      end
