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
      subroutine csout(iseq,maxbl,maxcs,igridg,levelg,ncs,sx,sy,sz,
     .                 stot,pav,ptav,tav,ttav,xmav,fmdot,cfxp,cfyp,
     .                 cfzp,cfdp,cflp,cftp,cfxv,cfyv,cfzv,cfdv,cflv,
     .                 cftv,cfxmom,cfymom,cfzmom,cfdmom,cflmom,cftmom,
     .                 cfxtot,cfytot,cfztot,cfdtot,cfltot,cfttot,
     .                 icsinfo)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Write out control surface mass flow and momemtum/forces
c               Original coding by R. Cedar, GE Aircraft Engines
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension levelg(maxbl),igridg(maxbl)
      dimension sx(maxcs),sy(maxcs),sz(maxcs),stot(maxcs),
     .          pav(maxcs),ptav(maxcs),tav(maxcs),ttav(maxcs),
     .          xmav(maxcs),fmdot(maxcs),
     .          cfxp(maxcs),cfyp(maxcs),cfzp(maxcs),
     .          cfdp(maxcs),cflp(maxcs),cftp(maxcs),
     .          cfxv(maxcs),cfyv(maxcs),cfzv(maxcs),
     .          cfdv(maxcs),cflv(maxcs),cftv(maxcs),
     .          cfxmom(maxcs),cfymom(maxcs),cfzmom(maxcs),
     .          cfdmom(maxcs),cflmom(maxcs),cftmom(maxcs),
     .          cfxtot(maxcs),cfytot(maxcs),cfztot(maxcs),
     .          cfdtot(maxcs),cfltot(maxcs),cfttot(maxcs),
     .          icsinfo(maxcs,9)
c
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
c
c..   print results for each surface on the current global/embedded level
c
      ics1 = 0
      iprttot = 0
      do ics=1,ncs
         nbl1 = icsinfo(ics,1)
         if (levelg(nbl1) .ge. lglobal .and.
     .       levelg(nbl1) .le. levt) then
            ics1 = ics1 + 1
            if (icsinfo(ics,9).ne.0) iprttot = 1
            write(17,1000) ics1
            igrid = igridg(icsinfo(ics,1))
            write(17,1100) igrid,
     .                     icsinfo(ics,1),icsinfo(ics,2),
     .                     icsinfo(ics,3),icsinfo(ics,4),
     .                     icsinfo(ics,5),icsinfo(ics,6),
     .                     icsinfo(ics,7),icsinfo(ics,8),
     .                     icsinfo(ics,9)
            write(17,1113) real(sx(ics)), real(sy(ics)),
     .                     real(sz(ics)),real(stot(ics))
            if (icsinfo(ics,8).eq.0) then
               write(17,1114)
               write(17,1115) real(pav(ics)),real(ptav(ics)),
     .                        real(tav(ics)),real(ttav(ics)),
     .                        real(xmav(ics))
               write(17,1116) real(fmdot(ics))
            end if
            write(17,1200)
            write(17,1300) real(cfxp(ics)),real(cfyp(ics)),
     .                     real(cfzp(ics)),real(cftp(ics)),
     .                     real(cflp(ics)),real(cfdp(ics))
            if (icsinfo(ics,8).eq.1)
     .      write(17,1400) real(cfxv(ics)),real(cfyv(ics)),
     .                     real(cfzv(ics)),real(cftv(ics)),
     .                     real(cflv(ics)),real(cfdv(ics))
            if (icsinfo(ics,8).eq.0)
     .      write(17,1500) real(cfxmom(ics)),real(cfymom(ics)),
     .                     real(cfzmom(ics)),real(cftmom(ics)),
     .                     real(cflmom(ics)),real(cfdmom(ics))
            write(17,1600) real(cfxtot(ics)),real(cfytot(ics)),
     .                     real(cfztot(ics)),real(cfttot(ics)),
     .                     real(cfltot(ics)),real(cfdtot(ics))
         end if
      end do
c
c..   print out totals
c
      if (iprttot.gt.0) then
         write(17,1700)
         write(17,1113) real(sx(ncs+1)),real(sy(ncs+1)),
     .                  real(sz(ncs+1)),real(stot(ncs+1))
         write(17,1116) real(fmdot(ncs+1))
         write(17,1200)
         write(17,1300) real(cfxp(ncs+1)),real(cfyp(ncs+1)),
     .                  real(cfzp(ncs+1)),real(cftp(ncs+1)),
     .                  real(cflp(ncs+1)),real(cfdp(ncs+1))
         write(17,1400) real(cfxv(ncs+1)),real(cfyv(ncs+1)),
     .                  real(cfzv(ncs+1)),real(cftv(ncs+1)),
     .                  real(cflv(ncs+1)),real(cfdv(ncs+1))
         write(17,1500) real(cfxmom(ncs+1)),real(cfymom(ncs+1)),
     .                  real(cfzmom(ncs+1)),real(cftmom(ncs+1)),
     .                  real(cflmom(ncs+1)),real(cfdmom(ncs+1))
         write(17,1600) real(cfxtot(ncs+1)),real(cfytot(ncs+1)),
     .                  real(cfztot(ncs+1)),real(cfttot(ncs+1)),
     .                  real(cfltot(ncs+1)),real(cfdtot(ncs+1))
      end if
c
c..   formats
c
 1000 format(/,/,' Control Surface',i3,/,
     .          ' ===================',/)
 1100 format(' Grid',i3,' (Block',i3,')  i =',i3,',',i3,'  j =',i3,',',
     .                  i3,'  k =',i3,',',i3,'  iwall = ',i3,
     .                  '  Normal = ',i3,/)
 1113 format(' x-area = ',e12.4,' y-area = ',e12.4,' z-area = ',e12.4,
     .       ' total-area = ',e12.4)
 1114 format(/,' Mass averaged properties')
 1115 format(' P/Pinf      = ',e12.4,'         Pt/Pinf = ',e12.4,/
     .       ' T/Tinf      = ',e12.4,'         Tt/Tinf = ',e12.4,/
     .       ' Mach number = ',e12.4)
 1116 format(/,' Mass flow / (rhoinf*vinf*(L_R)**2) = ',e14.5)
 1200 format(/,22x,'x-force',7x,'y-force',7x,'z-force',
     .         3x,'resultant-force',1x,'lift-force',4x,'drag-force')
 1300 format(' Pressure force  ',6e14.5)
 1400 format(' Viscous force   ',6e14.5)
 1500 format(' Thrust force    ',6e14.5)
 1600 format(' Total force     ',6e14.5,/)            
 1700 format(/,/,' Totals for all sufaces (global level) with normal',
     .          ' .ne. 0',/,
     .          ' =========================================',
     .          '==============='/)
c
      return
      end
