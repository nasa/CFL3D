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
      subroutine dird(jdim,kdim,idim,x,y,z,sj,sk,si,snj0,snk0,sni0,
     .                snjm,snkm,snim,w,ivisc,nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c**********************************************************************
c      Purpose:  Evaluate directed distance from k=0 wall and i=0/j=0 
c      wall for use in evaluating the Baldwin-Lomax turbulence model.
c**********************************************************************
c
c      input arrays   : 
c       x,y,z         : grid point positions
c       sj,sk,si      : metrics
c                       (direction cosines,areas,speeds of cell faces)
c      output arrays  :
c       snj0,snjm     : directed distance from j=0,max wall 
c                        of cell-centers
c       snk0,snkm     : directed distance from k=0,max wall     
c                        of cell-centers
c       sni0,snim     : directed distance from i=0,max wall     
c                       of cell-centers
c      scratch arrays :
c       w(1-3)        : cell-centered positions
c       w(4)          : j=0/k=0/i=0 wall temporary
c
c**********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension x(jdim,kdim,idim),
     .          y(jdim,kdim,idim),
     .          z(jdim,kdim,idim)
      dimension sj(jdim,kdim,idim-1,5),
     .          sk(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5)
      dimension snj0(jdim-1,kdim-1,idim-1),
     .          snk0(jdim-1,kdim-1,idim-1),
     .          sni0(jdim-1,kdim-1,idim-1)
      dimension snjm(jdim-1,kdim-1,idim-1),
     .          snkm(jdim-1,kdim-1,idim-1),
     .          snim(jdim-1,kdim-1,idim-1)
c
      dimension w(jdim-1,kdim-1,idim-1,4),ivisc(3)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
c      cell-center locations
c 
      do 1000 i=1,idim1
      do 1000 k=1,kdim1
      do 1000 j=1,jdim1
      w(j,k,i,1) = (  x(j  ,k  ,i) + x(j  ,k  ,i+1) + 
     .                x(j  ,k+1,i) + x(j  ,k+1,i+1) + 
     .                x(j+1,k  ,i) + x(j+1,k  ,i+1) + 
     .                x(j+1,k+1,i) + x(j+1,k+1,i+1) ) 
      w(j,k,i,2) = (  y(j  ,k  ,i) + y(j  ,k  ,i+1) + 
     .                y(j  ,k+1,i) + y(j  ,k+1,i+1) + 
     .                y(j+1,k  ,i) + y(j+1,k  ,i+1) + 
     .                y(j+1,k+1,i) + y(j+1,k+1,i+1) ) 
      w(j,k,i,3) = (  z(j  ,k  ,i) + z(j  ,k  ,i+1) + 
     .                z(j  ,k+1,i) + z(j  ,k+1,i+1) + 
     .                z(j+1,k  ,i) + z(j+1,k  ,i+1) + 
     .                z(j+1,k+1,i) + z(j+1,k+1,i+1) ) 
 1000 continue
c
c*********************************************************
c      k=0 wall
c*********************************************************
c 
      if (ivisc(3).gt.1) then
c
      k = 1
      do 2000 i=1,idim1
      do 2000 j=1,jdim1
      w(j,k,i,4) = .25e0*((  x(j  ,k,i) + x(j  ,k,i+1) + 
     .                       x(j+1,k,i) + x(j+1,k,i+1) ) *sk(j,k,i,1)
     .                  + (  y(j  ,k,i) + y(j  ,k,i+1) + 
     .                       y(j+1,k,i) + y(j+1,k,i+1) ) *sk(j,k,i,2)
     .                  + (  z(j  ,k,i) + z(j  ,k,i+1) + 
     .                       z(j+1,k,i) + z(j+1,k,i+1) ) *sk(j,k,i,3) )
 2000 continue
c
      do 3000 i=1,idim1
      do 3000 k=1,kdim1
      do 3000 j=1,jdim1
      snk0(j,k,i) = ( w(j,k,i,1)*sk(j,1,i,1)
     .              + w(j,k,i,2)*sk(j,1,i,2)
     .              + w(j,k,i,3)*sk(j,1,i,3) )*.125e0 - w(j,1,i,4)
 3000 continue
c
      icount  = 0
      snk0mg  = 1.0e+06
      do 9900 i=1,idim1
      do 9900 k=1,kdim1
      snk0min = q8smin(jdim1,snk0(1,k,i))
c
      if (real(snk0min).le.0.) then
c     write(15,8001)
 8001 format(1x,43h *** probable miscue in dird *** check grid)
c
      if (icount .le. 100) then
         do 8902 j=1,jdim1
         if (real(snk0(j,k,i)) .gt. 0.e0) go to 8902
         icount = icount+1
c        if (icount.le.100) write(15,8903) j,k,i,real(snk0(j,k,i))
 8902    continue
      end if
 8903 format(1x,11h j,k,i,snk0,3i5,e12.5)
      end if
c
      snk0mg = ccmin(snk0min,snk0mg)
 9900 continue
c   ***ensure grid wrapping back in toward body doesn't mess up snk0
      if(icount .gt. 0) then
c       write(15,8111)
 8111   format(1x,39h *** The code is resetting some values.)
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),8112)
 8112   format(3x,49hnegative directed distances have been detected...)
         nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),8113) 
 8113   format(3x,50hthis is not usually a problem, but check output on)
         nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),8114)
 8114   format(3x,47hunit 15 to be sure (may need to check grid too))
        do 4967 i=1,idim1
        do 4967 j=1,jdim1
        do 4967 k=2,kdim1
          snk0(j,k,i)=ccmax(snk0(j,k,i),snk0(j,k-1,i))
 4967   continue
      end if
c     write(15,9901) real(snk0mg)
 9901 format(1x,10h snk0,min=,e12.5)
      end if
c*********************************************************
c      k=kdim wall
c*********************************************************
c
      if (ivisc(3).gt.1) then
c
      k  = kdim
      km = k - 1
      do 200 i=1,idim1
      do 200 j=1,jdim1
      w(j,km,i,4) = .25e0*((  x(j  ,k,i) + x(j  ,k,i+1) +
     .                        x(j+1,k,i) + x(j+1,k,i+1) ) *sk(j,k,i,1)
     .                   + (  y(j  ,k,i) + y(j  ,k,i+1) +
     .                        y(j+1,k,i) + y(j+1,k,i+1) ) *sk(j,k,i,2)
     .                   + (  z(j  ,k,i) + z(j  ,k,i+1) +
     .                        z(j+1,k,i) + z(j+1,k,i+1) ) *sk(j,k,i,3) )
  200 continue
c
      do 300 i=1,idim1
      do 300 k=1,kdim1
      do 300 j=1,jdim1
      snkm(j,k,i) =  w(j,kdim-1,i,4) - ( w(j,k,i,1)*sk(j,kdim,i,1)
     .              +                    w(j,k,i,2)*sk(j,kdim,i,2)
     .              +                    w(j,k,i,3)*sk(j,kdim,i,3) )
     .                                 * .125e0
  300 continue
c
c     write(15,11)
   11 format(2x,43hcomputed directed distance from k=kdim wall)
      do 600 i=1,idim1,5
      do 600 j=1,jdim1,5
c     write(15,1601) j,i,real(snkm(j,1,i)),real(snkm(j,kdim1,i))
 1601 format(1x,24h j,i,snkm(1),snkm(kdim1),2i5,e12.5,f12.3)
  600 continue
c
      icount  = 0
      snk0mg  = 1.0e+06
      do 990 i=1,idim1
      do 990 k=1,kdim1
      snk0min = q8smin(jdim1,snkm(1,k,i))
c
      if (real(snk0min).le.0.) then
c     write(15,801)
  801 format(1x,43h *** probable miscue in dird *** check grid)
c
      if (icount .le. 100) then
         do 892 j=1,jdim1
         if (real(snkm(j,k,i)) .gt. 0.e0) go to 892
         icount = icount+1
c        if (icount.le.100) write(15,893) j,k,i,real(snkm(j,k,i))
  892    continue
      end if
  893 format(1x,11h j,k,i,snkm,3i5,e12.5)
      end if
c
      snk0mg = ccmin(snk0min,snk0mg)
  990 continue
c   ***ensure grid wrapping back in toward body doesn't mess up snkm
      if(icount .gt. 0) then
c       write(15,811)
  811   format(1x,39h *** The code is resetting some values.)
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),812)
  812   format(3x,49hnegative directed distances have been detected...)
         nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),813)
  813   format(3x,50hthis is not usually a problem, but check output on)
         nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),814)
  814   format(3x,47hunit 15 to be sure (may need to check grid too))
        do 496 i=1,idim1
        do 496 j=1,jdim1
        do 496 k=1,kdim1-1
          snkm(j,k,i)=ccmax(snkm(j,k,i),snkm(j,k+1,i))
 496   continue
      end if
c     write(15,1991) real(snk0mg)
 1991 format(1x,10h snkm,min=,e12.5)
      end if
c
c*********************************************************
c      j=0 wall
c*********************************************************
c
      if (ivisc(2).gt.1) then
c
      j = 1
      do 2100 i=1,idim1
      do 2100 k=1,kdim1
      w(j,k,i,4) = .25e0*((  x(j  ,k,i) + x(j  ,k,i+1) + 
     .                       x(j,k+1,i) + x(j,k+1,i+1) ) *sj(j,k,i,1)
     .                  + (  y(j  ,k,i) + y(j  ,k,i+1) + 
     .                       y(j,k+1,i) + y(j,k+1,i+1) ) *sj(j,k,i,2)
     .                  + (  z(j  ,k,i) + z(j  ,k,i+1) + 
     .                       z(j,k+1,i) + z(j,k+1,i+1) ) *sj(j,k,i,3) )
2100  continue
c
      do 3100 i=1,idim1
      do 3100 k=1,kdim1
      do 3100 j=1,jdim1
      snj0(j,k,i) = ( w(j,k,i,1)*sj(1,k,i,1)
     .              + w(j,k,i,2)*sj(1,k,i,2)
     .              + w(j,k,i,3)*sj(1,k,i,3) )*.125e0 - w(1,k,i,4)
 3100 continue
c
      icount  = 0
      snj0mg  = 1.0e+06
      do 9700 i=1,idim1
      do 9700 k=1,kdim1
      snj0min = q8smin(jdim1,snj0(1,k,i))
c
      if (real(snj0min).le.0.) then
c     write(15,8001)
c
      if (icount .le. 100) then
         do 8702 j=1,jdim1
         if (real(snj0(j,k,i)) .gt. 0.e0) go to 8702
         icount = icount+1
c        if (icount.le.100) write(15,8703) j,k,i,real(snj0(j,k,i))
 8702    continue
      end if
 8703 format(1x,11h j,k,i,snj0,3i5,e12.5)
      end if
c
      snj0mg = ccmin(snj0min,snj0mg)
 9700 continue
c   ***ensure grid wrapping back in toward body doesn't mess up snj0
      if(icount .gt. 0) then
c       write(15,8111)
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),8112)
        do 4968 i=1,idim1
        do 4968 k=1,kdim1
        do 4968 j=2,jdim1
          snj0(j,k,i)=ccmax(snj0(j,k,i),snj0(j-1,k,i))
 4968   continue
      end if
c
c     write(15,9701) real(snj0mg)
 9701 format(1x,10h snj0,min=,e12.5)
c
      end if
c*********************************************************
c      j=jdim wall
c*********************************************************
c
      if (ivisc(2).gt.1) then
c
      j  = jdim
      jm = j - 1
      do 201 i=1,idim1
      do 201 k=1,kdim1
      w(jm,k,i,4) = .25e0*((  x(j  ,k,i) + x(j  ,k,i+1) +
     .                        x(j,k+1,i) + x(j,k+1,i+1) ) *sj(j,k,i,1)
     .                   + (  y(j  ,k,i) + y(j  ,k,i+1) +
     .                        y(j,k+1,i) + y(j,k+1,i+1) ) *sj(j,k,i,2)
     .                   + (  z(j  ,k,i) + z(j  ,k,i+1) +
     .                        z(j,k+1,i) + z(j,k+1,i+1) ) *sj(j,k,i,3) )
  201 continue
c
      do 301 i=1,idim1
      do 301 k=1,kdim1
      do 301 j=1,jdim1
      snjm(j,k,i) =  w(jdim-1,k,i,4) - ( w(j,k,i,1)*sj(jdim,k,i,1)
     .              +                    w(j,k,i,2)*sj(jdim,k,i,2)
     .              +                    w(j,k,i,3)*sj(jdim,k,i,3) )
     .                                 * .125e0
  301 continue
c
c     write(15,116)
  116 format(2x,43hcomputed directed distance from j=jdim wall)
      do 601 i=1,idim1,5
      do 601 k=1,kdim1,5
c     write(15,701) k,i,real(snjm(1,k,i)),real(snjm(jdim1,k,i))
  701 format(1x,24h k,i,snjm(1),snjm(jdim1),2i5,e12.5,f12.3)
  601 continue
c
      icount  = 0
      snj0mg  = 1.0e+06
      do 991 i=1,idim1
      do 991 j=1,jdim1
      snj0min = q8smin(kdim1,snjm(j,1,i))
c
      if (real(snj0min).le.0.) then
c     write(15,1811)
 1811 format(1x,43h *** probable miscue in dird *** check grid)
c
      if (icount .le. 100) then
         do 894 k=1,kdim1
         if (real(snjm(j,k,i)) .gt. 0.e0) go to 894
         icount = icount+1
c        if (icount.le.100) write(15,895) j,k,i,real(snjm(j,k,i))
  894    continue
      end if
  895 format(1x,11h j,k,i,snjm,3i5,e12.5)
      end if
c
      snj0mg = ccmin(snj0min,snj0mg)
  991 continue
c   ***ensure grid wrapping back in toward body doesn't mess up snjm
      if(icount .gt. 0) then
c       write(15,6811)
 6811   format(1x,39h *** The code is resetting some values.)
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),6812)
 6812   format(3x,49hnegative directed distances have been detected...)
         nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),6813)
 6813   format(3x,50hthis is not usually a problem, but check output on)
         nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),6814)
 6814   format(3x,47hunit 15 to be sure (may need to check grid too))
        do 498 i=1,idim1
        do 498 k=1,kdim1
        do 498 j=1,jdim1-1
          snjm(j,k,i)=ccmax(snjm(j,k,i),snjm(j+1,k,i))
 498   continue
      end if
c     write(15,6991) real(snj0mg)
 6991 format(1x,10h snjm,min=,e12.5)
      end if
c
c*********************************************************
c      i=0 wall
c*********************************************************
c
      if (ivisc(1).gt.1) then
c 
      i = 1
      do 4000 k=1,kdim1
      do 4000 j=1,jdim1
      w(j,k,i,4) = .25e0*((  x(j  ,k,i) + x(j  ,k+1,i) + 
     .                       x(j+1,k,i) + x(j+1,k+1,i) ) *si(j,k,i,1)
     .                  + (  y(j  ,k,i) + y(j  ,k+1,i) + 
     .                       y(j+1,k,i) + y(j+1,k+1,i) ) *si(j,k,i,2)
     .                  + (  z(j  ,k,i) + z(j  ,k+1,i) + 
     .                       z(j+1,k,i) + z(j+1,k+1,i) ) *si(j,k,i,3) )
 4000 continue
c
      do 5000 i=1,idim1
      do 5000 k=1,kdim1
      do 5000 j=1,jdim1
      sni0(j,k,i) = ( w(j,k,i,1)*si(j,k,1,1)
     .              + w(j,k,i,2)*si(j,k,1,2)
     .              + w(j,k,i,3)*si(j,k,1,3) )*.125e0 - w(j,k,1,4)
 5000 continue
c
      icount  = 0
      sni0mg  = 1.0e+06
      do 9000 i=1,idim1
      do 9000 k=1,kdim1
      sni0min = q8smin(jdim1,sni0(1,k,i))
c
      if (real(sni0min).le.0.) then
c     write(15,8001)
c
      if (icount .le. 100) then
         do 8002 j=1,jdim1
         if (real(sni0(j,k,i)) .gt. 0.e0 ) go to 8002
         icount = icount+1
c        if (icount.le.100) write(15,8003) j,k,i,real(sni0(j,k,i))
 8002    continue
      end if
 8003 format(1x,11h j,k,i,sni0,3i5,e12.5)
      end if
c
      sni0mg = ccmin(sni0min,sni0mg)
 9000 continue
c   ***ensure grid wrapping back in toward body doesn't mess up sni0
      if(icount .gt. 0) then
c       write(15,8111)
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),8112)
        do 4969 k=1,kdim1
        do 4969 j=1,jdim1
        do 4969 i=2,idim1
          sni0(j,k,i)=ccmax(sni0(j,k,i),sni0(j,k,i-1))
 4969   continue
      end if
c
c     write(15,9001) real(sni0mg)
 9001 format(1x,10h sni0,min=,e12.5)
      end if
c*********************************************************
c      i=idim wall
c*********************************************************
c
      if (ivisc(1).gt.1) then
c
      i  = idim
      im = i - 1
      do 207 k=1,kdim1
      do 207 j=1,jdim1
      w(j,k,im,4) = .25e0*((  x(j  ,k,i) + x(j  ,k+1,i) +
     .                        x(j+1,k,i) + x(j+1,k+1,i) ) *si(j,k,i,1)
     .                   + (  y(j  ,k,i) + y(j  ,k+1,i) +
     .                        y(j+1,k,i) + y(j+1,k+1,i) ) *si(j,k,i,2)
     .                   + (  z(j  ,k,i) + z(j  ,k+1,i) +
     .                        z(j+1,k,i) + z(j+1,k+1,i) ) *si(j,k,i,3) )
  207 continue
c
      do 307 i=1,idim1
      do 307 k=1,kdim1
      do 307 j=1,jdim1
      snim(j,k,i) =  w(j,k,idim-1,4) - ( w(j,k,i,1)*si(j,k,idim,1)
     .              +                    w(j,k,i,2)*si(j,k,idim,2)
     .              +                    w(j,k,i,3)*si(j,k,idim,3) )
     .                                 * .125e0
  307 continue
c
c     write(15,126)
  126 format(2x,43hcomputed directed distance from i=idim wall)
      do 607 k=1,kdim1,5
      do 607 j=1,jdim1,5
c     write(15,704) j,k,real(snim(j,k,1)),real(snjm(j,k,idim1))
  704 format(1x,24h j,k,snim(1),snim(idim1),2i5,e12.5,f12.3)
  607 continue
c
      icount  = 0
      sni0mg  = 1.0e+06
      do 997 i=1,idim1
      do 997 k=1,kdim1
      sni0min = q8smin(jdim1,snim(1,k,i))
c
      if (real(sni0min).le.0.) then
c     write(15,817)
  817 format(1x,43h *** probable miscue in dird *** check grid)
c
      if (icount .le. 100) then
         do 897 j=1,jdim1
         if (real(snim(j,k,i)) .gt. 0.e0) go to 897
         icount = icount+1
c        if (icount.le.100) write(15,899) j,k,i,real(snim(j,k,i))
  897    continue
      end if
  899 format(1x,11h j,k,i,snim,3i5,e12.5)
      end if
c
      sni0mg = ccmin(sni0min,sni0mg)
  997 continue
c   ***ensure grid wrapping back in toward body doesn't mess up snim
      if(icount .gt. 0) then
c       write(15,6711)
 6711   format(1x,39h *** The code is resetting some values.)
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),6712)
 6712   format(3x,49hnegative directed distances have been detected...,
     .  /3x,50hthis is not usually a problem, but check output on,
     .  /3x,47hunit 15 to be sure (may need to check grid too))
        do 798 k=1,kdim1
        do 798 j=1,jdim1
        do 798 i=1,idim1-1
          snim(j,k,i)=ccmax(snim(j,k,i),snim(j,k,i+1))
 798   continue
      end if
c     write(15,6997) real(sni0mg)
 6997 format(1x,10h snim,min=,e12.5)
      end if
      return
      end
