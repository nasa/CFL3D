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
      subroutine cellvol(jdim,kdim,idim,x,y,z,sj,sk,si,vol,t,nou,bou,
     .                   nbuf,ibufdim,myid,mblk2nd,maxbl,nbl,iflagv1,
     .                   iflagv,imin,imax,jmin,jmax,kmin,kmax)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Calculate the cell volumes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),
     .          z(jdim,kdim,idim),vol(jdim,kdim,idim-1)
      dimension t(jdim*kdim,15),mblk2nd(maxbl)
      dimension sj(jdim*kdim,idim-1,5),sk(jdim*kdim,idim-1,5),
     .          si(jdim*kdim,idim,5)
c
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt 
      common /zero/ iexp
c
      eps  =  10.**(-iexp+1)
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
c     directed area and magnitude (i=constant face)
c
      i = 1
      n = jdim*kdim1-1
c
cdir$ ivdep
      do 1001 izz=1,n
      t(izz,7) = si(izz,i,1)*si(izz,i,4)
      t(izz,8) = si(izz,i,2)*si(izz,i,4)
      t(izz,9) = si(izz,i,3)*si(izz,i,4)
c
c      average point in i=constant face
c
      t(izz,13) = x(izz,1,i)+x(izz+1,1,i)+x(izz,2,i)+x(izz+1,2,i)
      t(izz,14) = y(izz,1,i)+y(izz+1,1,i)+y(izz,2,i)+y(izz+1,2,i)
      t(izz,15) = z(izz,1,i)+z(izz+1,1,i)+z(izz,2,i)+z(izz+1,2,i)
c
c     loop through cross-flow planes
c
 1001 continue
      t1 = -1.e0/24.e0
      do 1000 i=1,idim1
c
c      average point in volume
c
      l = i+1
      n = jdim*kdim1-1
cdir$ ivdep
      do 1002 izz=1,n
      t(izz,10) = x(izz,1,l)+x(izz+1,1,l)+x(izz,2,l)+x(izz+1,2,l)
     .           +t(izz,13)
      t(izz,11) = y(izz,1,l)+y(izz+1,1,l)+y(izz,2,l)+y(izz+1,2,l)
     .           +t(izz,14)
      t(izz,12) = z(izz,1,l)+z(izz+1,1,l)+z(izz,2,l)+z(izz+1,2,l)
     .           +t(izz,15)
c
c     accumulate volume as sum of pentahedrons
c
      vol(izz,1,i) = (2.e0*t(izz,13)-t(izz,10))*t(izz,7)
     .              +(2.e0*t(izz,14)-t(izz,11))*t(izz,8)
     .              +(2.e0*t(izz,15)-t(izz,12))*t(izz,9)
 1002 continue
c
c     directed area (j=constant face)
c
      n = jdim*kdim1
c
cdir$ ivdep
      do 1003 izz=1,n
      t(izz,7) = sj(izz,i,1)*sj(izz,i,4)
      t(izz,8) = sj(izz,i,2)*sj(izz,i,4)
      t(izz,9) = sj(izz,i,3)*sj(izz,i,4)
c
c     average point in j=constant face
c
      t(izz,4) = x(izz,1,i)+x(izz,1,i+1)+x(izz,2,i)+x(izz,2,i+1)
      t(izz,5) = y(izz,1,i)+y(izz,1,i+1)+y(izz,2,i)+y(izz,2,i+1)
      t(izz,6) = z(izz,1,i)+z(izz,1,i+1)+z(izz,2,i)+z(izz,2,i+1)
c
 1003 continue
      n = n-1
cdir$ ivdep
      do 1004 izz=1,n
      vol(izz,1,i) = vol(izz,1,i)+(2.e0*t(izz,4)  -t(izz,10))*t(izz,7)
     .                           +(2.e0*t(izz,5)  -t(izz,11))*t(izz,8)
     .                           +(2.e0*t(izz,6)  -t(izz,12))*t(izz,9)
     .                           -(2.e0*t(izz+1,4)-t(izz,10))*t(izz+1,7)
     .                           -(2.e0*t(izz+1,5)-t(izz,11))*t(izz+1,8)
     .                           -(2.e0*t(izz+1,6)-t(izz,12))*t(izz+1,9)
 1004 continue
c
c     directed area and magnitude (k=constant face)
c
      n = jdim*kdim-1
c
cdir$ ivdep
      do 1005 izz=1,n
      t(izz,7) = sk(izz,i,1)*sk(izz,i,4)
      t(izz,8) = sk(izz,i,2)*sk(izz,i,4)
      t(izz,9) = sk(izz,i,3)*sk(izz,i,4)
c
c     average point in k=constant face
c
      t(izz,4) = x(izz+1,1,i)+x(izz+1,1,i+1)+x(izz,1,i)+x(izz,1,i+1)
      t(izz,5) = y(izz+1,1,i)+y(izz+1,1,i+1)+y(izz,1,i)+y(izz,1,i+1)
      t(izz,6) = z(izz+1,1,i)+z(izz+1,1,i+1)+z(izz,1,i)+z(izz,1,i+1)
c
 1005 continue
      n = n-jdim
cdir$ ivdep
      do 1006 izz=1,n
      vol(izz,1,i) = vol(izz,1,i)
     .                     +(2.e0*t(izz,4)     -t(izz,10))*t(izz,7)
     .                     +(2.e0*t(izz,5)     -t(izz,11))*t(izz,8)
     .                     +(2.e0*t(izz,6)     -t(izz,12))*t(izz,9)
     .                     -(2.e0*t(izz+jdim,4)-t(izz,10))*t(izz+jdim,7)
     .                     -(2.e0*t(izz+jdim,5)-t(izz,11))*t(izz+jdim,8)
     .                     -(2.e0*t(izz+jdim,6)-t(izz,12))*t(izz+jdim,9)
 1006 continue
c
c     directed area and magnitude (i=constant face)
c
      l = i+1
      n = jdim*kdim1-1
c
cdir$ ivdep
      do 1007 izz=1,n
      t(izz,7) = si(izz,l,1)*si(izz,l,4)
      t(izz,8) = si(izz,l,2)*si(izz,l,4)
      t(izz,9) = si(izz,l,3)*si(izz,l,4)
c
c      average point in i=constant face
c
      t(izz,13) = t(izz,10)-t(izz,13)
      t(izz,14) = t(izz,11)-t(izz,14)
      t(izz,15) = t(izz,12)-t(izz,15)
c
      vol(izz,1,i) = vol(izz,1,i)-(2.e0*t(izz,13)-t(izz,10))*t(izz,7)
     .                           -(2.e0*t(izz,14)-t(izz,11))*t(izz,8)
     .                           -(2.e0*t(izz,15)-t(izz,12))*t(izz,9)
      vol(izz,1,i) = t1*vol(izz,1,i)
 1007 continue
      do 5000 kk=1,kdim1
      vol(jdim,kk,i) = vol(jdim1,kk,i)
 5000 continue
c
cdir$ ivdep
      do 1008 izz=1,jdim
      vol(izz,kdim,i) = vol(izz,kdim1,i)
 1008 continue
 1000 continue
c
      if (negvol .eq. 0) then
         do 64 i=1,idim1
         do 64 k=1,kdim
         do 64 j=1,jdim
         if (real(vol(j,k,i)).gt.0.0) go to 64
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),82) i,j,k,nbl
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
   64    continue
   82    format(39h stopping ... negative volume at i,j,k=,3i5,
     .           6h block,i5)
      else
         imin = 100000000
         imax = 0
         jmin = 100000000
         jmax = 0
         kmin = 100000000
         kmax = 0
         do 65 i=1,idim1
         do 65 k=1,kdim
         do 65 j=1,jdim
           if (real(vol(j,k,i)).lt.real(eps)) then 
             if(i.lt.imin) imin = i
             if(i.gt.imax) imax = i
             if(j.lt.jmin) jmin = j
             if(j.gt.jmax) jmax = j
             if(k.lt.kmin) kmin = k
             if(k.gt.kmax) kmax = k
             iflagv = 1
             vol(j,k,i) = ccabs(vol(j,k,i))
             if(iflagv1.eq.1) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),83) i,j,k,nbl
             end if
           end if
   65    continue
   83    format(39h WARNING  ... negative volume at i,j,k=,3i5,
     .             6h block,i5,14h not stopping!)
         imin = imin - 10
         jmin = jmin - 10
         kmin = kmin - 10
         imax = imax + 10
         jmax = jmax + 10
         kmax = kmax + 10
         if(imin.lt.1   ) imin=1
         if(jmin.lt.1   ) jmin=1
         if(kmin.lt.1   ) kmin=1
         if(imax.gt.idim) imax=idim
         if(jmax.gt.jdim) jmax=jdim
         if(kmax.gt.kdim) kmax=kdim
      end if
c
      return
      end
