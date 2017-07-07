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
      subroutine metric(jdim,kdim,idim,x,y,z,sj,sk,si,t,t1,nbl,iflag,
     .                  icnt,nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                  nbckdim,ibcinfo,jbcinfo,kbcinfo,maxbl,maxseg,
     .                  nou,bou,nbuf,ibufdim,myid,mblk2nd)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Calculate the cell-interface directed areas.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension sj(jdim*kdim,idim-1,5),sk(jdim*kdim,idim-1,5),
     .          si(jdim*kdim,idim,5),mblk2nd(maxbl)
      dimension t(jdim*kdim,6),t1(jdim*kdim,idim,5)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
c
      common /sklton/ isklton
      common /twod/ i2d
      common /singular/ atol
      common /zero/iexp
c
c***********************************************************************
c  metrics:
c          si(1-3)...components of unit normal to i-face (directed area)
c          si( 4 )...area of  i-face
c          si( 5 ).. speed of i-face
c
c          sj(1-3)...components of unit normal to j-face (directed area)
c          sj( 4 )...area of  j-face
c          sj( 5 ).. speed of j-face
c
c          sk(1-3)...components of unit normal to k-face (directed area)
c          sk( 4 )...area of  k-face
c          sk( 5 ).. speed of k-face
c
c  notes:
c        1) the normal to a cell face is determined via the cross product
c           of two diagonal vectors connecting oposite corners of the cell
c           face. 
c        2) a unit normal is obtained by dividing by the magnitude 
c           of the inner product. 
c        3) The area of the cell face is given by  one-half the magnitude 
c           of the cross product of the two diagonal vectors
c        4) in a jdim*kdim*idim grid there are really only 
c                   idim*(jdim-1)*(kdim-1) i-faces
c                   jdim*(idim-1)*(kdim-1) j-faces
c                   kdim*(idim-1)*(jdim-1) k-faces
c           however the metric arrays are dimensioned larger than this -
c           for efficient vectorization, loops throughout the code run 
c           over these fictitious faces. thus, for safety, the metrics 
c           for these fictitious faces are set to the corresponding metrics 
c           for the neighboring face.
c
c        5) unsteady metric terms si(5), sj(5), sk(5) are set to zero
c           here, and must be reset elsewhere if the mesh is actually 
c           in motion
c
c        6) tolerance to automatically detect collapsed (singular) metrics
c           is set by atol (now set in subroutine readkey, and may be
c           overwritten via keyword-driven input
c
c***********************************************************************
c
      icntmax = 10
      atol2   = 10.**(-iexp+1)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
c***********************************************************************
c
c     metrics for i=constant surfaces
c
c***********************************************************************
c
c     *** interior faces ***
c
      n = jdim*kdim-jdim
c
      do 1040 i=2,idim1
cdir$ ivdep
      do 1050 izz=1,n
c
c     components of vectors connecting opposite corners of cell j,k
      t(izz,1) = x(izz+1,1,i)-x(izz,2,i)
      t(izz,2) = y(izz+1,1,i)-y(izz,2,i)
      t(izz,3) = z(izz+1,1,i)-z(izz,2,i)
      t(izz,4) = x(izz+1,2,i)-x(izz,1,i)
      t(izz,5) = y(izz+1,2,i)-y(izz,1,i)
      t(izz,6) = z(izz+1,2,i)-z(izz,1,i)
c
c     cross product of vectors
      si(izz,i,1) =  t(izz,2)*t(izz,6)-t(izz,3)*t(izz,5)
      si(izz,i,2) = -t(izz,1)*t(izz,6)+t(izz,3)*t(izz,4)
      si(izz,i,3) =  t(izz,1)*t(izz,5)-t(izz,2)*t(izz,4)
c
c     magnitude of cross product
      si(izz,i,4) = si(izz,i,1)*si(izz,i,1)+si(izz,i,2)*si(izz,i,2)+
     .              si(izz,i,3)*si(izz,i,3)
 1050 continue
c
c     due to the ordering in izz, the cross product above is incorrect
c     at the fictitious interfaces izz=k*jdim (k=1...kdim-1), and may
c     have a zero value. for safety in the loop below, first set these
c     fictitious interface values to temporary safe values.
c 
      do 1055 k=1,kdim1
      izz = k*jdim
      si(izz,i,1) = 0.
      si(izz,i,2) = 0.
      si(izz,i,3) = 0.
      si(izz,i,4) = 4.
 1055 continue
c
cdir$ ivdep
c     store metrics:
      do 1060 izz=1,n
      si(izz,i,4) = 1.e0/sqrt(si(izz,i,4))
      si(izz,i,1) = si(izz,i,1)*si(izz,i,4)
      si(izz,i,2) = si(izz,i,2)*si(izz,i,4)
      si(izz,i,3) = si(izz,i,3)*si(izz,i,4)
      si(izz,i,4) = 0.5e0/si(izz,i,4)
      si(izz,i,5) = 0.0
 1060 continue
c
 1040 continue
c
c     *** i=1/idim faces ***
c
      do 1000 m=1,2
c
      if(m.eq.1) then
        i=1
        nseg = nbci0(nbl)
      else
        i=idim
        nseg = nbcidim(nbl)
      end if
c
      do 1010 ns=1,nseg
      js = ibcinfo(nbl,ns,2,m)
      je = ibcinfo(nbl,ns,3,m) - 1
      ks = ibcinfo(nbl,ns,4,m)
      ke = ibcinfo(nbl,ns,5,m) - 1
      mtyp = ibcinfo(nbl,ns,1,m)
c
      asum = 0.
c
      do 1020 j=js,je
      do 1020 k=ks,ke
      izz = (k-1)*jdim+j
c
c     components of vectors connecting opposite corners of cell j,k
      t(izz,1) = x(j+1,k,i)   - x(j,k+1,i)
      t(izz,2) = y(j+1,k,i)   - y(j,k+1,i)
      t(izz,3) = z(j+1,k,i)   - z(j,k+1,i)
      t(izz,4) = x(j+1,k+1,i) - x(j,k,i)
      t(izz,5) = y(j+1,k+1,i) - y(j,k,i)
      t(izz,6) = z(j+1,k+1,i) - z(j,k,i)
c
c     cross product of vectors
      si(izz,i,1) =  t(izz,2)*t(izz,6)-t(izz,3)*t(izz,5)
      si(izz,i,2) = -t(izz,1)*t(izz,6)+t(izz,3)*t(izz,4)
      si(izz,i,3) =  t(izz,1)*t(izz,5)-t(izz,2)*t(izz,4)
c
c     magnitude of cross product
      si(izz,i,4) = si(izz,i,1)*si(izz,i,1)+si(izz,i,2)*si(izz,i,2)+
     .              si(izz,i,3)*si(izz,i,3)
c
      asum = asum + sqrt(si(izz,i,4))
      if (si(izz,i,4) .eq. 0.) si(izz,i,4) = atol2
c
 1020 continue
c
      asum = 0.5*asum
      if (real(asum) .lt. real(atol)) then
c
c       collapsed metrics
c
        if (isklton .gt. 0) then
           if (m.eq.1) then
              if (isklton .eq. 1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),101) js,je+1,ks,ke+1
              else
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),103) js,je+1,ks,ke+1
              end if
           else
              if (isklton .eq. 1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),102) js,je+1,ks,ke+1
              else
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),104) js,je+1,ks,ke+1
              end if
           end if
        end if
c  
        izzndx = i-(2*m)+3
        do 1030 j=js,je
        do 1030 k=ks,ke
        izz = (k-1)*jdim+j
c       set directions on collapsed face equal to those at closest interior face
        si(izz,i,1) = si(izz,izzndx,1)
        si(izz,i,2) = si(izz,izzndx,2)
        si(izz,i,3) = si(izz,izzndx,3)
        si(izz,i,4) = 0.e0
        si(izz,i,5) = si(izz,izzndx,5)
 1030   continue
c
      else
c
c       non-singular metrics
c
        do 1045 j=js,je
        do 1045 k=ks,ke
        izz = (k-1)*jdim+j
        si(izz,i,4) = 1.e0/sqrt(si(izz,i,4))
        si(izz,i,1) = si(izz,i,1)*si(izz,i,4)
        si(izz,i,2) = si(izz,i,2)*si(izz,i,4)
        si(izz,i,3) = si(izz,i,3)*si(izz,i,4)
        si(izz,i,4) = 0.5e0/si(izz,i,4)
        si(izz,i,5) = 0.0
 1045   continue
c
      end if
c
 1010 continue
c
 1000 continue
c
c     fill in extra values of si for safety
      do 1070 i=1,idim
c     set metrics at jdim
      do 1080 k=1,kdim-1
      izz  = (k-1)*jdim+jdim
      izz1 = izz - 1
      si(izz,i,1) = si(izz1,i,1)
      si(izz,i,2) = si(izz1,i,2)
      si(izz,i,3) = si(izz1,i,3)
      si(izz,i,4) = si(izz1,i,4)
      si(izz,i,5) = si(izz1,i,5)
 1080 continue  
c     set metrics at kdim
cdir$ ivdep
      do 1090 j=1,jdim
      izz  = jdim*(kdim-1) + j
      izz1 = izz - jdim
      si(izz,i,1) = si(izz1,i,1)
      si(izz,i,2) = si(izz1,i,2)
      si(izz,i,3) = si(izz1,i,3)
      si(izz,i,4) = si(izz1,i,4)
      si(izz,i,5) = si(izz1,i,5)
 1090 continue
 1070 continue
c
c   Don't do the following checks if iflag=-1:
      if(iflag .ne. -1) then
c
c   Check for non-planar i-planes when 2-D:
      if (i2d .eq. 1) then
      izz1=jdim*kdim
      do 1235 i=1,idim
      do 1235 m=1,3
      temp=si(1,i,m)
      do 1234 izz=2,izz1
        if (abs(real(si(izz,i,m)-temp)) .gt. 1.e-5) then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' Error. Grid not planar at i,nbl='',
     .     2i5)') i,nbl
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' The i-planes need to be planar'',
     .     '' for 2-D (i2d=1)'')')
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'(''   ... or else the problem may'',
     .     '' be crossing lines (negative volumes) in the grid'')')
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
 1234 continue
 1235 continue
      do j=1,jdim
      do k=1,kdim
        if (abs(x(j,k,2)-x(j,k,1)) .gt. 1.e-5 .or.
     .      abs(z(j,k,2)-z(j,k,1)) .gt. 1.e-5) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' The 2 i-planes need to be'',
     .       '' identical (in their planar values) for 2-D (i2d=1)'')')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''   ... or else the problem may'',
     .       '' be that IALPH is set incorrectly for your grid'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
      enddo
      enddo
      end if
c   Check for large changes in metric from one i to the next,
c   indicating problem with the grid
      do k=1,kdim-1
      do j=1,jdim-1
        izz=(k-1)*jdim+j
        do i=1,idim-1 
          delta1 = ccabs(si(izz,i+1,1)-si(izz,i,1))
          delta2 = ccabs(si(izz,i+1,2)-si(izz,i,2))
          delta3 = ccabs(si(izz,i+1,3)-si(izz,i,3))
          if(real(delta1).gt.1.99 .or. real(delta2).gt.1.99 .or.
     +       real(delta3).gt.1.99) then
            iflag=1
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''   FATAL si grid normal direction'',
     +       '' change near j,k,i,i+1='',4i5)')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''     ... suspect bad grid'',4i5)') 
     +       j,k,i,i+1
          else if(real(delta1).gt.1.5 .or. real(delta2).gt.1.5 .or.
     +       real(delta3).gt.1.5) then
             icnt = icnt + 1
             if (icnt.le.icntmax) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'(''   WARNING: Dramatic si grid'',
     +          '' norm direction change (>120deg)'')')
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'(''     near j,k,i,i+1='',4i5)') 
     +          j,k,i,i+1
             else if (icnt.eq.icntmax+1) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'(''   NOTE: no further WARNINGS'',
     +          '' will be output...better check this grid!'')')
             end if
          end if
        enddo
      enddo
      enddo
      end if
c
c***********************************************************************
c
c     metrics for j=constant surfaces
c
c***********************************************************************
c
c     *** interior faces ***
c
c     note: interior metrics for j=constant surfaces are initially 
c           calculated with a reordered izz index (compared to si and  
c           sk metrics) to permit vectorization and yet skip over j=1/jdim,
c           which may have collapsed metrics. after the metrics are 
c           evaluated, the proper order is reset.
c
      do 2035 i=1,idim
      do 2035 j=1,jdim
      do 2035 k=1,kdim
      izz = (j-1)*kdim+k
      t1(izz,i,1) = x(j,k,i)
      t1(izz,i,2) = y(j,k,i)
      t1(izz,i,3) = z(j,k,i)
 2035 continue
c
      ns = kdim+1
      ne = (jdim-1)*kdim
c
      do 2040 i=1,idim1
cdir$ ivdep
      do 2050 izz=ns,ne
c
c     components of vectors connecting opposite corners of cell i,k
      t(izz,1) = t1(izz+1,i,1)-t1(izz,i+1,1)
      t(izz,2) = t1(izz+1,i,2)-t1(izz,i+1,2)
      t(izz,3) = t1(izz+1,i,3)-t1(izz,i+1,3)
      t(izz,4) = t1(izz+1,i+1,1)-t1(izz,i,1)
      t(izz,5) = t1(izz+1,i+1,2)-t1(izz,i,2)
      t(izz,6) = t1(izz+1,i+1,3)-t1(izz,i,3)
c
c     cross product of vectors
      sj(izz,i,1) =  t(izz,2)*t(izz,6)-t(izz,3)*t(izz,5)
      sj(izz,i,2) = -t(izz,1)*t(izz,6)+t(izz,3)*t(izz,4)
      sj(izz,i,3) =  t(izz,1)*t(izz,5)-t(izz,2)*t(izz,4)
c
c     magnitude of cross product
      sj(izz,i,4) = sj(izz,i,1)*sj(izz,i,1)+sj(izz,i,2)*sj(izz,i,2)+
     .              sj(izz,i,3)*sj(izz,i,3)
 2050 continue
c
c     due to the ordering in izz, the cross product above is incorrect
c     at the fictitious interfaces izz=j*kdim (j=1...jdim-1), and may
c     have a zero value. for safety in the loop below, first set these
c     fictitious interface values to temporary safe values.
c
      do 2055 j=1,jdim1
      izz = kdim*j
      sj(izz,i,1) = 0.
      sj(izz,i,2) = 0.
      sj(izz,i,3) = 0.
      sj(izz,i,4) = 4.
 2055 continue
c
c     store temporary metrics (these are not in correct order): 
      do 2060 izz=ns,ne
      t1(izz,i,4) = 1.e0/sqrt(sj(izz,i,4))
      t1(izz,i,1) = sj(izz,i,1)*t1(izz,i,4)
      t1(izz,i,2) = sj(izz,i,2)*t1(izz,i,4)
      t1(izz,i,3) = sj(izz,i,3)*t1(izz,i,4)
      t1(izz,i,4) = 0.5e0/t1(izz,i,4)
      t1(izz,i,5) = 0.0
 2060 continue
c
c     store metrics in correct order
      do 2061 j=2,jdim1
      do 2062 k=1,kdim
      izz1 = (j-1)*kdim+k
      izz = (k-1)*jdim+j
      sj(izz,i,1) = t1(izz1,i,1)
      sj(izz,i,2) = t1(izz1,i,2)
      sj(izz,i,3) = t1(izz1,i,3)
      sj(izz,i,4) = t1(izz1,i,4)
      sj(izz,i,5) = t1(izz1,i,5)
 2062 continue
 2061 continue
c
 2040 continue
c
c     *** j=1/jdim faces ***
c
      do 2000 m=1,2
c
      if(m.eq.1) then
        j=1
        nseg = nbcj0(nbl)
      else
        j=jdim
        nseg = nbcjdim(nbl)
      end if
c
      do 2010 ns=1,nseg
      is = jbcinfo(nbl,ns,2,m)
      ie = jbcinfo(nbl,ns,3,m) - 1
      ks = jbcinfo(nbl,ns,4,m)
      ke = jbcinfo(nbl,ns,5,m) - 1
      mtyp = jbcinfo(nbl,ns,1,m)
c
      asum = 0.
c
      do 2020 i=is,ie
      do 2020 k=ks,ke
      izz = (k-1)*jdim+j
c
c     components of vectors connecting opposite corners of cell i,k
      t(izz,1) = x(j,k+1,i) - x(j,k,i+1)
      t(izz,2) = y(j,k+1,i) - y(j,k,i+1)
      t(izz,3) = z(j,k+1,i) - z(j,k,i+1)
      t(izz,4) = x(j,k+1,i+1) - x(j,k,i)
      t(izz,5) = y(j,k+1,i+1) - y(j,k,i)
      t(izz,6) = z(j,k+1,i+1) - z(j,k,i)
c
c     cross product of vectors
      sj(izz,i,1) =  t(izz,2)*t(izz,6)-t(izz,3)*t(izz,5)
      sj(izz,i,2) = -t(izz,1)*t(izz,6)+t(izz,3)*t(izz,4)
      sj(izz,i,3) =  t(izz,1)*t(izz,5)-t(izz,2)*t(izz,4)
c
c     magnitude of cross product
      sj(izz,i,4) = sj(izz,i,1)*sj(izz,i,1)+sj(izz,i,2)*sj(izz,i,2)+
     .              sj(izz,i,3)*sj(izz,i,3)
c
      asum = asum + sqrt(sj(izz,i,4))
      if (sj(izz,i,4) .eq. 0.) sj(izz,i,4) = atol2
c
 2020 continue
c
      asum = 0.5*asum 
      if (real(asum) .lt. real(atol)) then
c
c       collapsed metrics
c
        if (isklton .gt. 0) then
           if (m.eq.1) then
              if (isklton .eq. 1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),201) is,ie+1,ks,ke+1
              else
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),203) is,ie+1,ks,ke+1
              end if
           else
              if (isklton .eq. 1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),202) is,ie+1,ks,ke+1
              else
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),204) is,ie+1,ks,ke+1
              end if
           end if
        end if
c
        do 2030 i=is,ie
        do 2030 k=ks,ke
        izz = (k-1)*jdim+j
        izzndx = izz-(2*m)+3
c       set directions on collapsed face equal to those at closest interior face
        sj(izz,i,1) = sj(izzndx,i,1)
        sj(izz,i,2) = sj(izzndx,i,2)
        sj(izz,i,3) = sj(izzndx,i,3)
        sj(izz,i,4) = 0.e0
        sj(izz,i,5) = sj(izzndx,i,5)
 2030   continue
c
      else
c
c       non-singular metrics
c
        do 2045 i=is,ie
        do 2045 k=ks,ke
        izz = (k-1)*jdim+j
        sj(izz,i,4) = 1.e0/sqrt(sj(izz,i,4))
        sj(izz,i,1) = sj(izz,i,1)*sj(izz,i,4)
        sj(izz,i,2) = sj(izz,i,2)*sj(izz,i,4)
        sj(izz,i,3) = sj(izz,i,3)*sj(izz,i,4)
        sj(izz,i,4) = 0.5e0/sj(izz,i,4)
        sj(izz,i,5) = 0.0
 2045   continue
c
      end if
c
 2010 continue
c
 2000 continue
c
c     fill in extra values of sj for safety
      do 2070 i=1,idim-1
c     set metrics at kdim
cdir$ ivdep
      do 2090 j=1,jdim
      izz  = jdim*(kdim-1) + j
      izz1 = izz - jdim
      sj(izz,i,1) = sj(izz1,i,1)
      sj(izz,i,2) = sj(izz1,i,2)
      sj(izz,i,3) = sj(izz1,i,3)
      sj(izz,i,4) = sj(izz1,i,4)
      sj(izz,i,5) = sj(izz1,i,5)
 2090 continue
 2070 continue
c
c   Don't do the following checks if iflag=-1:
      if(iflag .ne. -1) then
c
c   Check for large changes in metric from one j to the next,
c   indicating problem with the grid
      do i=1,idim-1
      do k=1,kdim-1
        do j=1,jdim-1
          izz=(k-1)*jdim+j
          izzp1=(k-1)*jdim+j+1
          delta1 = ccabs(sj(izzp1,i,1)-sj(izz,i,1))
          delta2 = ccabs(sj(izzp1,i,2)-sj(izz,i,2))
          delta3 = ccabs(sj(izzp1,i,3)-sj(izz,i,3))
          if(real(delta1).gt.1.99 .or. real(delta2).gt.1.99 .or.
     +       real(delta3).gt.1.99) then
            iflag=1
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''   FATAL sj grid normal direction'',
     +       '' change near j,j+1,k,i='',4i5)')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''     ... suspect bad grid'',4i5)')
     +       j,j+1,k,i
          else if(real(delta1).gt.1.5 .or. real(delta2).gt.1.5 .or.
     +       real(delta3).gt.1.5) then
             icnt = icnt + 1
             if (icnt.le.icntmax) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'(''   WARNING: Dramatic sj grid'',
     +          '' norm direction change (>120deg)'')')
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'(''     near j,j+1,k,i='',4i5)') 
     +          j,j+1,k,i
             else if (icnt.eq.icntmax+1) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'(''   NOTE: no further WARNINGS'',
     +          '' will be output...better check this grid!'')')
             end if
          end if
        enddo
      enddo
      enddo
      end if
c
c***********************************************************************
c
c     metrics for k=constant surfaces
c
c***********************************************************************
c
c     *** interior faces ***
c
      ns = jdim+1
      ne = jdim*(kdim-1)
c
      do 3040 i=1,idim1
cdir$ ivdep
      do 3050 izz=ns,ne
c
c     components of vectors connecting opposite corners of cell i,j
      t(izz,1) = x(izz+1,1,i)-x(izz,1,i+1)
      t(izz,2) = y(izz+1,1,i)-y(izz,1,i+1)
      t(izz,3) = z(izz+1,1,i)-z(izz,1,i+1)
      t(izz,4) = x(izz,1,i)-x(izz+1,1,i+1)
      t(izz,5) = y(izz,1,i)-y(izz+1,1,i+1)
      t(izz,6) = z(izz,1,i)-z(izz+1,1,i+1)
c
c     cross product of vectors
      sk(izz,i,1) =  t(izz,2)*t(izz,6)-t(izz,3)*t(izz,5)
      sk(izz,i,2) = -t(izz,1)*t(izz,6)+t(izz,3)*t(izz,4)
      sk(izz,i,3) =  t(izz,1)*t(izz,5)-t(izz,2)*t(izz,4)
c
c     magnitude of cross product
      sk(izz,i,4) = sk(izz,i,1)*sk(izz,i,1)+sk(izz,i,2)*sk(izz,i,2)+
     .              sk(izz,i,3)*sk(izz,i,3)
 3050 continue
c
c     due to the ordering in izz, the cross product above is incorrect
c     at the fictitious interfaces izz=j*kdim (j=1...jdim-1), and may
c     have a zero value. for safety in the loop below, first set these
c     fictitious interface values to temporary safe values.
c
      do 3055 k=1,kdim
      izz = jdim*k
      sk(izz,i,1) = 0.e0
      sk(izz,i,2) = 0.e0
      sk(izz,i,3) = 0.e0
      sk(izz,i,4) = 4.e0
 3055 continue
c
c     store metrics:
      do 3060 izz=ns,ne
      sk(izz,i,4) = 1.e0/sqrt(sk(izz,i,4))
      sk(izz,i,1) = sk(izz,i,1)*sk(izz,i,4)
      sk(izz,i,2) = sk(izz,i,2)*sk(izz,i,4)
      sk(izz,i,3) = sk(izz,i,3)*sk(izz,i,4)
      sk(izz,i,4) = 0.5e0/sk(izz,i,4)
      sk(izz,i,5) = 0.0
 3060 continue
c
 3040 continue
c
c     *** k=1/kdim faces ***
c
      do 3000 m=1,2
c
      if(m.eq.1) then
        k=1
        nseg = nbck0(nbl)
      else
        k=kdim
        nseg = nbckdim(nbl)
      end if
c
      do 3010 ns=1,nseg
      is = kbcinfo(nbl,ns,2,m)
      ie = kbcinfo(nbl,ns,3,m) - 1
      js = kbcinfo(nbl,ns,4,m)
      je = kbcinfo(nbl,ns,5,m) - 1
      mtyp = kbcinfo(nbl,ns,1,m)
c
      asum = 0.
c
      do 3020 i=is,ie
      do 3020 j=js,je
      izz = (k-1)*jdim+j
c
c     components of vectors connecting opposite corners of cell i,j
      t(izz,1) = x(j+1,k,i) - x(j,k,i+1)
      t(izz,2) = y(j+1,k,i) - y(j,k,i+1)
      t(izz,3) = z(j+1,k,i) - z(j,k,i+1)
      t(izz,4) = x(j,k,i) - x(j+1,k,i+1)
      t(izz,5) = y(j,k,i) - y(j+1,k,i+1)
      t(izz,6) = z(j,k,i) - z(j+1,k,i+1)
c
c     cross product of vectors
      sk(izz,i,1) =  t(izz,2)*t(izz,6)-t(izz,3)*t(izz,5)
      sk(izz,i,2) = -t(izz,1)*t(izz,6)+t(izz,3)*t(izz,4)
      sk(izz,i,3) =  t(izz,1)*t(izz,5)-t(izz,2)*t(izz,4)
c
c     magnitude of cross product
      sk(izz,i,4) = sk(izz,i,1)*sk(izz,i,1)+sk(izz,i,2)*sk(izz,i,2)+
     .              sk(izz,i,3)*sk(izz,i,3)
c
      asum = asum + sqrt(sk(izz,i,4))
      if (sk(izz,i,4) .eq. 0.) sk(izz,i,4) = atol2
c
 3020 continue
c
      asum = 0.5*asum
      if (real(asum) .lt. real(atol)) then
c
c       collapsed metrics
c
        if (isklton .gt. 0) then
           if (m.eq.1) then
              if (isklton .eq. 1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),301) is,ie+1,js,je+1
              else
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),303) is,ie+1,js,je+1
              end if
           else
              if (isklton .eq. 1) then
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),302) is,ie+1,js,je+1
              else
                 nou(1) = min(nou(1)+1,ibufdim)
                 write(bou(nou(1),1),304) is,ie+1,js,je+1
              end if
           end if
        end if
c 
        do 3030 i=is,ie
        do 3030 j=js,je
        izz = (k-1)*jdim+j
        izzndx = izz-(2*jdim*m)+(3*jdim)
c       set directions on collapsed face equal to those at closest interior face
        sk(izz,i,1) = sk(izzndx,i,1)
        sk(izz,i,2) = sk(izzndx,i,2)
        sk(izz,i,3) = sk(izzndx,i,3)
        sk(izz,i,4) = 0.e0
        sk(izz,i,5) = sk(izzndx,i,5)
 3030   continue
c
      else
c
c       non-singular metrics
c
        do 3045 i=is,ie
        do 3045 j=js,je
        izz = (k-1)*jdim+j
        sk(izz,i,4) = 1.e0/sqrt(sk(izz,i,4))
        sk(izz,i,1) = sk(izz,i,1)*sk(izz,i,4)
        sk(izz,i,2) = sk(izz,i,2)*sk(izz,i,4)
        sk(izz,i,3) = sk(izz,i,3)*sk(izz,i,4)
        sk(izz,i,4) = 0.5e0/sk(izz,i,4)
        sk(izz,i,5) = 0.0
 3045   continue
c
      end if
c
 3010 continue
c
 3000 continue
c
c     fill in extra values of sk for safety
      do 3070 i=1,idim-1
c     set metrics at jdim
cdir$ ivdep
      do 3090 k=1,kdim
      izz  = jdim*(k-1) + jdim
      izz1 = izz - 1
      sk(izz,i,1) = sk(izz1,i,1)
      sk(izz,i,2) = sk(izz1,i,2)
      sk(izz,i,3) = sk(izz1,i,3)
      sk(izz,i,4) = sk(izz1,i,4)
      sk(izz,i,5) = sk(izz1,i,5)
 3090 continue
 3070 continue
c
c   Don't do the following checks if iflag=-1:
      if(iflag .ne. -1) then
c
c   Check for large changes in metric from one k to the next,
c   indicating problem with the grid
      do i=1,idim-1
      do j=1,jdim-1
        do k=1,kdim-1
          izz=(k-1)*jdim+j
          izzp1=k*jdim+j
          delta1 = ccabs(sk(izzp1,i,1)-sk(izz,i,1))
          delta2 = ccabs(sk(izzp1,i,2)-sk(izz,i,2))
          delta3 = ccabs(sk(izzp1,i,3)-sk(izz,i,3))
          if(real(delta1).gt.1.99 .or. real(delta2).gt.1.99 .or.
     +       real(delta3).gt.1.99) then
            iflag=1
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''   FATAL sk grid normal direction'',
     +       '' change near j,k,k+1,i='',4i5)')
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''     ... suspect bad grid'',4i5)')
     +       j,k,k+1,i
          else if(real(delta1).gt.1.5 .or. real(delta2).gt.1.5 .or.
     +       real(delta3).gt.1.5) then
             icnt = icnt + 1
             if (icnt.le.icntmax) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'(''   WARNING: Dramatic sk grid'',
     +          '' norm direction change (>120deg)'')')
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'(''     near j,k,k+1,i='',4i5)') 
     +          j,k,k+1,i
             else if (icnt.eq.icntmax+1) then
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'(''   NOTE: no further WARNINGS'',
     +          '' will be output...better check this grid!'')')
             end if
          end if
        enddo
      enddo
      enddo
      end if
c
  101 format(3x,'singular metrics: i=1   ', 
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
  102 format(3x,'singular metrics: i=idim',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
  103 format(6x,'singular metrics: i=1   ',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
  104 format(6x,'singular metrics: i=idim',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
  201 format(3x,'singular metrics: j=1   ',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
  202 format(3x,'singular metrics: j=jdim',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
  203 format(6x,'singular metrics: j=1   ',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
  204 format(6x,'singular metrics: j=jdim',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
  301 format(3x,'singular metrics: k=1   ',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
  302 format(3x,'singular metrics: k=kdim',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
  303 format(6x,'singular metrics: k=1   ',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
  304 format(6x,'singular metrics: k=kdim',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
      return
      end
