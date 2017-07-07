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
      subroutine trace(icall,idum1,idum2,idum3,idum4,dum1,dum2,dum3)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Writes the search routine history for the current "to"
c     cell to unit 7.  
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /tracer/ itrace
c
c     itrace < 0, do not write history of current "to" cell
c     itrace = 0, overwrite history from previous "to" cell with current 
c     itrace = 1, retain the search history for ALL cells
c
      iunit = 7
c
      if(itrace.lt.0) return
c
      if(icall.eq.0) then
      write(iunit,100)idum1
  100 format(' ',2x,'*** Search History During Generalized',
     .       ' Coordinate Interpolation Number ',i3,' ***',/)
      end if
c
      if(icall.eq.1 .or. icall.eq.2) then
        if(itrace.eq.0) then
          rewind(iunit)
          write(iunit,100)idum1
        end if
      end if
c
      if(icall.eq.1) then
      write(iunit,101)idum2,idum3,real(dum1),real(dum2),real(dum3)
  101 format(' ',2x,'iterations for "to" grid:  cell j,k = ',i3,
     .',',i3,'  with center xc,yc,zc = ',e11.4,',',e11.4,',',e11.4)
      write(iunit,111)idum4
 111  format(' ',3x,'using ifit = ',i3)
      end if
c
      if(icall.eq.2) then
      write(iunit,102)idum2,idum3,real(dum1),real(dum2),real(dum3)
  102 format(' ',2x,'iterations for "to" grid:  cell j,k=',i3,
     .',',i3,'  with edge mid-point xc,yc,zc=',e11.4,',',
     .e11.4,',',e11.4)
      write(iunit,111)idum4
      end if
c
      if(icall.eq.3) then
      write(iunit,103) idum1
  103 format(' ',4x,'intern=',i3)
      end if
c
      if(icall.eq.4) then
      write(iunit,104)real(dum1),real(dum2)
  104 format(' ',8x,'xie,eta=',e10.3,',',e10.3,' (local values)')
      end if
c
      if(icall.eq.5) then
      write(iunit,105) idum1
105   format(' ','  could not find "to" cell center in "from"',
     . ' block',i3,/,'   will try another block')
      end if
c
      if(icall.eq.6) then
      write(iunit,106)idum1
  106 format(' ',4x,'checking to see if this point really lies'
     .       ,' in block ',i3)
      end if
c
      if(icall.eq.7) then
      write(iunit,107) idum1
  107 format(' ','    point actually lies in block ',i3)
      write(iunit,1071)
 1071 format(' ','    future searches to be done in this block')
      end if
c
      if(icall.eq.8) then
      write(iunit,108)
  108 format(' ','    incompatable orientation of "from" cell',
     .       ' and "to" cell')
      write(iunit,1081) 
 1081 format(' ','      ...will try in a nearby cell')
      end if
c
      if(icall.eq.9) then
      write(iunit,109)idum1,idum2,idum3
  109 format(' ',6x,'using y-z to invert in "from" grid',i3,
     .       '  cell j,k=',i3,',',i3)
      end if
c
      if(icall.eq.10) then
      write(iunit,1010)idum1,idum2,idum3
 1010 format(' ',6x,'using x-z to invert in "from" grid',i3,
     .       '  cell j,k=',i3,',',i3)
      end if
c
      if(icall.eq.11) then
      write(iunit,1011)idum1,idum2,idum3
 1011 format(' ',6x,'using x-y to invert in "from" grid',i3,
     .       '  cell j,k=',i3,',',i3)
      end if
c
      if(icall.eq.41) then
      write(iunit,1041)
 1041 format(' ',4x,'search routine off track...using minimum',
     .       ' distance search to get back on track')
      end if
c
      if(icall.eq.42) then
      write(iunit,1042)
1042  format(' ',4x,'frozen convergence...attempting to break',
     .       ' cycle via minimum distance search')
      end if
c
      if(icall.eq.50) then
      write(iunit,1050)idum1,idum2
 1050 format(' ','marking point j,k=',i3,i3,
     .      ' as an orphan ')
      end if
c
      if(icall.eq.99) then
      if (idum1.ne.idum3.or.idum2.ne.idum4)
     .write(iunit,1099)idum1,idum2,idum3,idum4
 1099 format(' ','branch cut modification: before, jp,kp=',i3,i3,
     .      ' after, jpc,kpc=',i3,i3)
      end if
      return
      end
