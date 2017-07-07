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
      program maggie
c
c     $Id$
c
c
c***********************************************************************
c
c     establishes zone-to-zone communication file for cfl3d multi-block
c     applications with the chimera oveset grid option
c
c     If you find a coding error or have suggestions for improved
c     implementation, contact Bob Biedron, 757-864-2156
c                                          r.t.biedron@larc.nasa.gov
c
c***********************************************************************
c
c  input files:
c 
c               grid.bin   input grids in CFL3D format 
c                          x(jdim,kdim,ldim),y,z (binary)
c
c               unit 5     input file which contains grid
c                          connection information (ascii)
c
c  output files:
c
c               ovrlp.bin  interpolation data (binary)
c                          required as an input file for CFL3D  
c
c               unit 6     primary information file (ascii)
c              
c               unit 7     secondary information file (ascii)
c                          contains search routine history from search
c                          for interpolation stencils. turned on or off
c                          in subroutines hintpt and ointpt as follows:
c                            itrace < 0, do not write search history 
c                              for current fringe or boundary pt
c                            itrace = 0, overwrite history from previous 
c                              point  with current (small file) 
c                            itrace = 1, retain the search history for ALL 
c                              points (may get huge file)
c                          itrace < 0 recomended unless error in search is
c                          suspected. itrace option slows down search.
c                          itrace is set near the top of subroutines HOLE
c                          and OUTER, and may be turned on in one or both
c                          routines
c                         
c               unit 8     terciary information file (ascii)  
c
c               unit 9     plot3d diagnostic file (binary)
c                          contains points segregated into field, fringe, 
c                          or boundary points, as well as orphan points 
c                          (ie points requiring extrapolation stencils) 
c                          and points having illegal stencils...file may
c                          be large. set iplt3d > 0 in input file to get 
c                          this file
c                         
c              temp_grd.n  grid data, where n = mesh number (binary)
c                          ...deleted after successful run
c 
c              temp_cen.n  cell-center grid data, where n = mesh number (binary)
c                          ...deleted after successful run
c  
c              temp_ibl.n  iblank files corresponding to each grid, 
c                          where n = mesh number (binary)
c                          ...deleted after successful run
c 
c              temp_int.n  interpolation data, where n = mesh number (binary) 
c                          ...deleted after successful run
c 
c              unit 30 is used for all temporary I/O
c
c             hole_grd.nh  inital hole boundary, where nh = hole number (ascii)
c                          useful for verifying hole input data.
c                          to get this file, must have "hole no." set 
c                          to "-hole no." in input for hole definition
c
c             hole_nrm.nh  outward pointing normals to the inital hole 
c                          boundary, where nh = hole number (ascii)
c                          useful for verifying hole input data. 
c                          to get this file, must have "hole no." set 
c                          to "-hole no." in input for hole definition
c
c             unit 35 is used for hole boundary/normal output
c
c             NOTE: hole_grd.nh and hole_nrm.nh are NOT created for hole nh
c                   if hole nh in mesh m is formed only from coordinate
c                   surfaces in mesh m itself, regardless of the sign of
c                   hole no.
c
      include 'mag1.h'
c
      character*80 grid,outpt,ovrlp,plt3d
c 
      common /pltpt/  xbnd(idim),ybnd(idim),zbnd(idim),xorph(idim),
     .                yorph(idim),zorph(idim),xill(idim),yill(idim),
     .                zill(idim),nbnd,norph,nill
      common /files/  grid,outpt,ovrlp,plt3d 
c
      common /tol/    epsc
c
      itime = 0
      call cputim(itime)
c
      nbnd  = 0
      norph = 0
      nill  = 0
c 
c     epsc: tolerance on interpolation coefficients 
c
      epsc=.25e-3
c
c***********************************************************************
c                    open files
c***********************************************************************
c
c
      read(5,*)
c
      read(5,'(a60)')grid
c     following line can be uncommented only on CRAY machines
c     call asnfile(grid, '-F f77 -N ieee', IER)
      open(unit=10,file=grid,form='unformatted',status='old')
c
      read(5,'(a60)')outpt
      open(unit=6,file=outpt,form='formatted',status='unknown')
c
      read(5,'(a60)')ovrlp
      open(unit=2,file=ovrlp,form='unformatted',status='unknown')
c
      read(5,'(a60)')plt3d
c     following line can be uncommented only on CRAY machines
c     call asnfile(plt3d, '-F f77 -N ieee', IER)
      open(unit=9,file=plt3d,form='unformatted',status='unknown')
c
c     output banner
c
      write(6,83)
      write(6,83)
      write(6,87)
      write(6,9900)
 9900 format(2(2h *),38h           MAGGIE - CFL3D OVERSET-GRID,
     .13h PREPROCESSOR,11x,2(2h *))
      write(6,87)
      write(6,9990)
 9990 format(2(2h *),43h   VERSION 6.7 :  Computational Fluids Lab,,
     .15h Mail Stop 128,,4x,2(2h *),
     ./2(2h *),18x,41hNASA Langley Research Center, Hampton, VA,
     .3x,2(2h *),/2(2h *),18x,33hRelease Date:  February  1, 2017.,
     .11x,2(2h *))
      write(6,87)
      write(6,83)
      write(6,83)
   83 format(35(2h *))
   87 format(2(2h *),62x,2(2h *))
c
      write(6,88)
   88 format(/19hinput/output files:)
c
      write(6,'(''  '',a60)')grid
      write(6,'(''  '',a60)')outpt
      write(6,'(''  '',a60)')ovrlp
      write(6,'(''  '',a60)')plt3d
c
c***********************************************************************
c                    set up and initialize
c***********************************************************************
c 
      call blckdat
      call setup 
c 
c***********************************************************************
c                   compute interpolation data
c***********************************************************************
c 
      call hole 
c
      call outer
c 
c***********************************************************************
c                   diagnostic checks 
c***********************************************************************
c 
      call diagnos 
c
c***********************************************************************
c                   output connection file for flow solver  
c***********************************************************************
c 
      call output 
c
      itime = -1
      call cputim(itime)
c
      stop
      end 
c
      subroutine blckdat 
c 
c***********************************************************************
c     Purpose: the variables and scalars are given their initial or 
c     default values 
c***********************************************************************
c 
      include 'mag1.h'
      parameter (md2=mdim*mdim, mid=mdim*idim, md3=md2*mhldim,
     .           md4=mdim*mhldim)
c 
      common /book1/  ipntr(mdim,mhldim*mdim), npntr(mdim,mhldim*mdim),
     .                mhbs(mdim,mdim), mobs(mdim,mdim), nsets(mdim)
      common /book2/  dxint(idim), dyint(idim), dzint(idim),
     .                ibpts(mdim), jbpt(idim), kbpt(idim), lbpt(idim),
     .                ji(idim),ki(idim), li(idim)
      common /book3/  iblank(mlen) 
      common /bound1/ ihbtyp(mdim,mdim),jh1(ipmax),jh2(ipmax),
     .                kh1(ipmax),kh2(ipmax),lh1(ipmax),lh2(ipmax),
     .                ip1(ipmax),ip2(ipmax),mh(ipmax),ihplt(ipmax) 
      common /bound2/ iobtyp(mdim),nobtyp(mdim),
     .                jo1(ipmax),jo2(ipmax),ko1(ipmax),ko2(ipmax),
     .                lo1(ipmax),lo2(ipmax)
      common /grdlim/ jsrs(mdim), jsre(mdim), ksrs(mdim), ksre(mdim),
     .                lsrs(mdim), lsre(mdim) 
      common /grid1/  x(mlen), y(mlen), z(mlen)
      common /grid2/  mjmax(mdim), mkmax(mdim), mlmax(mdim)
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
c
      do 10 m=1,mdim
10    nsets(m) = 0 
c
      do 20 n=1,mhldim
      do 20 m=1,mdim
      ipntr(m,n) = 0
20    npntr(m,n) = 0
c
      do 30 m=1,mdim
30    ibpts(m) = 0
c
      do 40 i=1,idim
      jbpt(i)  = 0
      kbpt(i)  = 0
      lbpt(i)  = 0
      ji(i)    = 0
      ki(i)    = 0
      li(i)    = 0
      dxint(i) = 0.0
      dyint(i) = 0.0
40    dzint(i) = 0.0
c
      do 50 i=1,mlen
50    iblank(i) = 1
c
      do 60 m=1,mdim
      nhole(m) = 0
60    noutr(m) = 0 
c
      do 70 m=1,mdim 
      do 70 n=1,mhldim
      moutr(m,n) = 0
70    mhole(m,n) = 0
c  
      do 80 m=1,mdim
      do 80 n=1,mhldim
      do 80 nn=1,mdim
      loutr(m,n,nn) = 0
80    lhole(m,n,nn) = 0
c
      do 90 m=1,mdim
      jsrs(m) = 0
      jsre(m) = 0
      ksrs(m) = 0
      ksre(m) = 0 
      lsrs(m) = 0
      lsre(m) = 0
90    continue
c  
      do 100 m=1,mlen
      x(m) = 0.0
      y(m) = 0.0
100   z(m) = 0.0
c
      do 110 m=1,mdim
      mjmax(m) = 0
      mkmax(m) = 0
110   mlmax(m) = 0
c
      do 120 n=1,mdim
      do 120 m=1,mdim
120   ihbtyp(m,n) = 0
c
      do 130 ip=1,ipmax
      jh1(ip) = 0
      jh2(ip) = 0
      kh1(ip) = 0
      kh2(ip) = 0
      lh1(ip) = 0
      lh2(ip) = 0 
      jo1(ip) = 0
      jo2(ip) = 0
      ko1(ip) = 0
      ko2(ip) = 0
      lo1(ip) = 0
130   lo2(ip) = 0 
c
      do 140 m=1,mdim
      iobtyp(m) = 0
140   nobtyp(m) = 0 
c
      return 
      end 
      subroutine bound(j,k,l,x1,y1,z1,jd,kd,ld,x,y,z,ife)
c***************************************************************** 
c     Purpose: create boundary points to be interpolated
c***************************************************************** 
c
      include 'mag1.h'
c
      dimension x(jd,kd,ld),y(jd,kd,ld),z(jd,kd,ld)
      dimension ife(6)
c
      eps    = 1.0e-04
c
c     to reduce the distance that the outer boundary cells
c     are extrapolated, set t2 = 0 and t1 equal to the desired
c     fraction of the adjacent interior cell. Standard values
c     are t1 = 1.4 and t2 = 1./1.4, and generally give a nice
c     continuation of the interior cells, but the resulting
c     exprapolated boundary cells may be too large for tight
c     sitiuations
c
      t1     = 1.4
      t2     = 1./1.4
c---rtb
c     t1     = 0.5
c     t2     = 0.
c---rtb
c
      if(j.le.jd .and. k.le.kd .and. l.le.ld .and.
     &   j.ge.1  .and. k.ge.1  .and. l.ge.1) then
                 x1     = x(j,k,l)
                 y1     = y(j,k,l)
                 z1     = z(j,k,l)
      else if(j.gt.jd) then
                 ife(2) = ife(2) + 1
        dx1    = x(jd,k,l) - x(jd-1,k,l)
        dy1    = y(jd,k,l) - y(jd-1,k,l)
        dz1    = z(jd,k,l) - z(jd-1,k,l)
        rx1    = min( t1 , (dx1 + eps)/(x(jd-1,k,l)-x(jd-2,k,l)+eps) ) 
        ry1    = min( t1 , (dy1 + eps)/(y(jd-1,k,l)-y(jd-2,k,l)+eps) ) 
        rz1    = min( t1 , (dz1 + eps)/(z(jd-1,k,l)-z(jd-2,k,l)+eps) ) 
        rx1    = max( t2 , rx1 )
        ry1    = max( t2 , ry1 )
        rz1    = max( t2 , rz1 )
                 x1     = x(jd,k,l) + rx1*(1.+rx1*float(j-jd-1))*dx1
                 y1     = y(jd,k,l) + ry1*(1.+ry1*float(j-jd-1))*dy1
                 z1     = z(jd,k,l) + rz1*(1.+rz1*float(j-jd-1))*dz1
      else if(j.lt.1) then
                 ife(1) = ife(1) + 1
        dx1    = x(1,k,l) - x(2,k,l)
        dy1    = y(1,k,l) - y(2,k,l)
        dz1    = z(1,k,l) - z(2,k,l)
        rx1    = min( t1 , (dx1 + eps)/(x(2,k,l)-x(3,k,l)+eps) ) 
        ry1    = min( t1 , (dy1 + eps)/(y(2,k,l)-y(3,k,l)+eps) ) 
        rz1    = min( t1 , (dz1 + eps)/(z(2,k,l)-z(3,k,l)+eps) ) 
        rx1    = max( t2 , rx1 )
        ry1    = max( t2 , ry1 )
        rz1    = max( t2 , rz1 )
                 x1     = x(1,k,l) + rx1*(1.+rx1*float(-j))*dx1
                 y1     = y(1,k,l) + ry1*(1.+ry1*float(-j))*dy1
                 z1     = z(1,k,l) + rz1*(1.+rz1*float(-j))*dz1
      else if(k.gt.kd) then
                 ife(4) = ife(4) + 1
        dx1    = x(j,kd,l) - x(j,kd-1,l)
        dy1    = y(j,kd,l) - y(j,kd-1,l)
        dz1    = z(j,kd,l) - z(j,kd-1,l)
        rx1    = min( t1 , (dx1 + eps)/(x(j,kd-1,l)-x(j,kd-2,l)+eps) ) 
        ry1    = min( t1 , (dy1 + eps)/(y(j,kd-1,l)-y(j,kd-2,l)+eps) ) 
        rz1    = min( t1 , (dz1 + eps)/(z(j,kd-1,l)-z(j,kd-2,l)+eps) ) 
        rx1    = max( t2 , rx1 )
        ry1    = max( t2 , ry1 )
        rz1    = max( t2 , rz1 )
                 x1     = x(j,kd,l) + rx1*(1.+rx1*float(k-kd-1))*dx1
                 y1     = y(j,kd,l) + ry1*(1.+ry1*float(k-kd-1))*dy1
                 z1     = z(j,kd,l) + rz1*(1.+rz1*float(k-kd-1))*dz1
      else if(k.lt.1) then
                 ife(3) = ife(3) + 1
        dx1    = x(j,1,l) - x(j,2,l)
        dy1    = y(j,1,l) - y(j,2,l)
        dz1    = z(j,1,l) - z(j,2,l)
        rx1    = min( t1 , (dx1 + eps)/(x(j,2,l)-x(j,3,l)+eps) ) 
        ry1    = min( t1 , (dy1 + eps)/(y(j,2,l)-y(j,3,l)+eps) ) 
        rz1    = min( t1 , (dz1 + eps)/(z(j,2,l)-z(j,3,l)+eps) ) 
        rx1    = max( t2 , rx1 )
        ry1    = max( t2 , ry1 )
        rz1    = max( t2 , rz1 )
                 x1     = x(j,1,l) + rx1*(1.+rx1*float(-k))*dx1
                 y1     = y(j,1,l) + ry1*(1.+ry1*float(-k))*dy1
                 z1     = z(j,1,l) + rz1*(1.+rz1*float(-k))*dz1
      else if(l.gt.ld) then
                 ife(6) = ife(6) + 1
        dx1    = x(j,k,ld) - x(j,k,ld-1)
        dy1    = y(j,k,ld) - y(j,k,ld-1)
        dz1    = z(j,k,ld) - z(j,k,ld-1)
        rx1    = min( t1 , (dx1 + eps)/(x(j,k,ld-1)-x(j,k,ld-2)+eps) ) 
        ry1    = min( t1 , (dy1 + eps)/(y(j,k,ld-1)-y(j,k,ld-2)+eps) ) 
        rz1    = min( t1 , (dz1 + eps)/(z(j,k,ld-1)-z(j,k,ld-2)+eps) ) 
        rx1    = max( t2 , rx1 )
        ry1    = max( t2 , ry1 )
        rz1    = max( t2 , rz1 )
                 x1     = x(j,k,ld) + rx1*(1.+rx1*float(l-ld-1))*dx1
                 y1     = y(j,k,ld) + ry1*(1.+ry1*float(l-ld-1))*dy1
                 z1     = z(j,k,ld) + rz1*(1.+rz1*float(l-ld-1))*dz1
      else if(l.lt.1) then
                 ife(5) = ife(5) + 1
        dx1    = x(j,k,1) - x(j,k,2)
        dy1    = y(j,k,1) - y(j,k,2)
        dz1    = z(j,k,1) - z(j,k,2)
        rx1    = min( t1 , (dx1 + eps)/(x(j,k,2)-x(j,k,3)+eps) ) 
        ry1    = min( t1 , (dy1 + eps)/(y(j,k,2)-y(j,k,3)+eps) ) 
        rz1    = min( t1 , (dz1 + eps)/(z(j,k,2)-z(j,k,3)+eps) ) 
        rx1    = max( t2 , rx1 )
        ry1    = max( t2 , ry1 )
        rz1    = max( t2 , rz1 )
                 x1     = x(j,k,1) + rx1*(1.+rx1*float(-l))*dx1
                 y1     = y(j,k,1) + ry1*(1.+ry1*float(-l))*dy1
                 z1     = z(j,k,1) + rz1*(1.+rz1*float(-l))*dz1

            end if
c
            return
            end
      subroutine branch(m1,jd,kd,ld,jimage,kimage,limage,
     .                 x1,y1,z1)
c***********************************************************************
c     Purpose: Check for branch cuts in mesh m1 and set up connection arrays 
c     jimage,kimage,limage to allow search routine to pass across branch cuts 
c     in the grid. 
c***********************************************************************
c
      dimension x1(jd,kd,ld),y1(jd,kd,ld),z1(jd,kd,ld)
      dimension jimage(jd-1,kd-1,ld-1),kimage(jd-1,kd-1,ld-1),
     .          limage(jd-1,kd-1,ld-1),ibrnch(6)
c
      character*8 face(6)
c
c     eps sets threshold for detection of branch cuts...points that are
c     within eps of one another are considered to lie on a branch cut
c     this parameter may need to be adjusted depending on the accuracy
c     used in grid generation
c
      eps = 1.e-7
c
      face(1) = 'j = 1'
      face(2) = 'j = jdim'
      face(3) = 'k = 1'
      face(4) = 'k = kdim'
      face(5) = 'l = 1'
      face(6) = 'l = ldim'
c      
c     notes: the indicies of the image points found in this subroutine point 
c     to the lower left hand corner of the image cell. the test for branch 
c     cuts is performed on the grid of node values, while the image indicies
c     pertain to cell centers. for a grid of dimensions jd,kd,ld, the cell-
c     center mesh is jd-1,kd-1,ld-1; however, the range of stencil points 
c     extends only to jd-2,kd-2,ld-2.  Thus image point arrays are set at
c     the j,k,l=1 and j,k,l=jd-2,kd-2,ld-2. 
c
      write(6,*)'    checking mesh ',m1,' for branch cuts on',
     .' boundaries'
c
c     ibrnch: idicates type of branch cut on each of the 6 faces of mesh m1
c
c     face numbering: 1...j=1
c                     2...j=jd
c                     3...k=1
c                     4...k=kd
c                     5...l=1
c                     6...l=ld
c
c     ibrnch = 0 no branch cut
c            = 10 C-type in j
c            = 20 C-type in k
c            = 30 C-type in l
c            = 40 O-type in j
c            = 50 O-type in k
c            = 60 O-type in l
c
c     default to no branch cuts
c
      do 1 iface=1,6
1     ibrnch(iface) = 0
      do 5 j=1,jd-1
      do 5 k=1,kd-1
      do 5 l=1,ld-1
      jimage(j,k,l) = j
      kimage(j,k,l) = k  
      limage(j,k,l) = l 
5     continue
c***
c      mskip = 2
c      if(m1.eq.mskip) then
c        write(6,*)'        skipping branch cut check on this mesh'
c        return
c      end if
c*** 
c      
c     check for C-grid type branch cut on k=1 and k=kd boundaries
c
      if(kd.gt.2) then 
c
c       branch cut in j-direction
c
        if(jd.gt.2) then        
c
        do 100 k=1,kd,kd-1
        iface = 3
        if(k.eq.kd) iface = 4
        if(ibrnch(iface).gt.0) go to 100
        kcheck = 0
        js0 = 100000
        jq0 = 0
        ls0 = 100000
        lq0 = 0
        do 200 l=1,ld-1
        do 200 js=1,jd-2
        do 300 jq=js+1,jd 
        delxyz = abs( x1(js,k,l)-x1(jq,k,l) ) 
     .         + abs( y1(js,k,l)-y1(jq,k,l) )
     .         + abs( z1(js,k,l)-z1(jq,k,l) )
        if (delxyz.lt.eps) then
c          grid points js and jq lie on opposite sides of a branch cut
           del2 = abs( x1(js+1,k,l)-x1(jq-1,k,l) ) 
     .          + abs( y1(js+1,k,l)-y1(jq-1,k,l) )
     .          + abs( z1(js+1,k,l)-z1(jq-1,k,l) )
           if(del2.lt.eps)then
c            cell centers js and jq-1 lie on opposite sides of a branch cut
             if (k.eq.1) then
                if(jimage(js,1,l).eq.js) then
                  kcheck = kcheck+1
                  jimage(js,1,l)      = jq - 2
                  jimage(jq-2,1,l)    = js
                end if
             end if
             if (k.eq.kd) then
                if(jimage(js,kd-1,l).eq.js) then
                kcheck = kcheck+1
                jimage(js,kd-2,l)   = jq - 2
                jimage(jq-2,kd-2,l) = js
                end if
             end if
             if(js.lt.js0)js0 = js
             if(jq.gt.jq0)jq0 = jq
             if(l.lt.ls0) ls0 = l
             if(l.gt.lq0) lq0 = l
           end if
        end if
  300   continue
  200   continue
c
c       amalgamate points to avoid gap in branch cut mapping at tip...the 
c       above criteria for C-type branch cut fails at the very tip where
c       the upper surface becomes the lower and vice versa.
c 
        if(k.eq.1) kk = 1
        if(k.eq.kd) kk = kd-2
        do 201 l=1,ld-1
        do 301 j=2,jd-2
        if(jimage(j,kk,l).eq.j) then
          if(jimage(j+1,kk,l).ne.j+1.and.jimage(j-1,kk,l).ne.j-1) then
            jimage(j,kk,l) = j-1
c           write(6,*)' fixing hole in branch cut array at j,k,i = ',
c    .      j,k,l,' jimage = ',jimage(j,kk,l)
          end if 
        end if
301     continue
201     continue
c
        if(kcheck.gt.0) then
          ibrnch(iface) = 10
          write(6,*)'      C-type branch cut detected on k = ',k,
     .    ' boundary'
          write(6,*)'        for j = ',js0,' to j = ',jq0,
     .    ' and for l = ',ls0,' to ',' l = ',lq0
        end if
  100   continue
c
        end if
c
c       branch cut in l-direction
c        
        if(ld.gt.2) then
c
        do 110 k=1,kd,kd-1
        iface = 3
        if(k.eq.kd) iface = 4
        if(ibrnch(iface).gt.0) go to 110
        kcheck = 0
        js0 = 100000
        jq0 = 0
        ls0 = 100000
        lq0 = 0
        do 210 j=1,jd-1
        do 210 ls=1,ld-2
        do 310 lq=ls+1,ld 
        delxyz = abs( x1(j,k,ls)-x1(j,k,lq) ) 
     .         + abs( y1(j,k,ls)-y1(j,k,lq) )
     .         + abs( z1(j,k,ls)-z1(j,k,lq) )
        if (delxyz.lt.eps) then
c          grid points ls and lq lie on opposite sides of a branch cut
           del2 = abs( x1(j,k,ls+1)-x1(j,k,lq-1) ) 
     .          + abs( y1(j,k,ls+1)-y1(j,k,lq-1) )
     .          + abs( z1(j,k,ls+1)-z1(j,k,lq-1) )
           if(del2.lt.eps)then
c            cell centers ls and lq-1 lie on opposite sides of a branch cut
             if (k.eq.1) then
                if(limage(j,1,ls).eq.ls) then
                kcheck = kcheck+1
                limage(j,1,ls)   = lq - 2
                limage(j,1,lq-2) = ls
                end if
             end if
             if (k.eq.kd) then
                if(limage(j,kd-1,ls).eq.ls) then
                kcheck = kcheck+1
                limage(j,kd-2,ls)   = lq - 2
                limage(j,kd-2,lq-2) = ls
                end if
             end if
             if(ls.lt.ls0)ls0 = ls
             if(lq.gt.lq0)lq0 = lq
             if(j.lt.js0) js0 = j
             if(j.gt.jq0) jq0 = j
           end if
        end if
  310   continue
  210   continue
c
c       amalgamate points to avoid gap in branch cut mapping at tip...the
c       above criteria for C-type branch cut fails at the very tip where
c       the upper surface becomes the lower and vice versa.
c
        if(k.eq.1) kk = 1
        if(k.eq.kd) kk = kd-2
        do 211 j=1,jd-1
        do 311 l=2,ld-2
        if(limage(j,kk,l).eq.l) then
          if(limage(j,kk,l+1).ne.l+1.and.limage(j,kk,l-1).ne.l-1) then
            limage(j,kk,l) = l-1
c           write(6,*)' fixing hole in branch cut array at j,k,i = ',
c    .      j,k,l,' limage = ',limage(j,kk,l)
          end if
        end if
311     continue
211     continue
c
        if(kcheck.gt.0) then
          ibrnch(iface) = 30
          write(6,*)'      C-type branch cut detected on k = ',k,
     .    ' boundary'
          write(6,*)'        for j = ',js0,' to j = ',jq0,
     .    ' and for l = ',ls0,' to ',' l = ',lq0
        end if
  110   continue
c
        end if
c
      end if
c      
c     check for C-grid type branch cut on j=1 and j=jd boundaries 
c
      if(jd.gt.2) then
c
c       branch cut in k-direction
c
        if(kd.gt.2) then
c        
        do 120 j=1,jd,jd-1
        iface = 1
        if(j.eq.jd) iface = 2
        ls0 = 100000
        lq0 = 0
        ks0 = 100000
        kq0 = 0
        jcheck = 0
        do 220 l=1,ld-1
        do 220 ks=1,kd-2
        do 320 kq=ks+1,kd 
        delxyz = abs( x1(j,ks,l)-x1(j,kq,l) ) 
     .         + abs( y1(j,ks,l)-y1(j,kq,l) )
     .         + abs( z1(j,ks,l)-z1(j,kq,l) )
        if (delxyz.lt.eps) then
c          grid points ks and kq lie on opposite sides of a branch cut
           del2 = abs( x1(j,ks+1,l)-x1(j,kq-1,l) ) 
     .          + abs( y1(j,ks+1,l)-y1(j,kq-1,l) )
     .          + abs( z1(j,ks+1,l)-z1(j,kq-1,l) )
           if(del2.lt.eps)then
c            cell centers ks and kq-1 lie on opposite sides of a branch cut
             if (j.eq.1) then
                if(kimage(1,ks,l).eq.ks) then
                jcheck = jcheck+1
                kimage(1,ks,l)   = kq - 2
                kimage(1,kq-2,l) = ks
                end if
             end if
             if (j.eq.jd) then
                if(kimage(jd-1,ks,l).eq.ks) then
                jcheck = jcheck+1
                kimage(jd-2,ks,l)   = kq - 2
                kimage(jd-2,kq-2,l) = ks
                end if
             end if
             if(ks.lt.ks0)ks0 = ks
             if(kq.gt.kq0)kq0 = kq
             if(l.lt.ls0) ls0 = l
             if(l.gt.lq0) lq0 = l
           end if
        end if
  320   continue
  220   continue
c
c       amalgamate points to avoid gap in branch cut mapping at tip...the
c       above criteria for C-type branch cut fails at the very tip where
c       the upper surface becomes the lower and vice versa.
c
        if(j.eq.1) jj = 1
        if(j.eq.jd) jj = jd-2
        do 221 l=1,ld-1
        do 321 k=2,kd-2
        if(kimage(jj,k,l).eq.k) then
          if(kimage(jj,k+1,l).ne.k+1.and.kimage(jj,k-1,l).ne.k-1) then
            kimage(jj,k,l) = k-1
c           write(6,*)' fixing hole in branch cut array at j,k,i = ',
c    .      j,k,l,' kimage = ',kimage(jj,k,l)
          end if
        end if
321     continue
221     continue
c
        if(jcheck.gt.0) then
          ibrnch(iface) = 20
          write(6,*)'      C-type branch cut detected on j = ',j,
     .    ' boundary'
          write(6,*)'        for k = ',ks0,' to k = ',kq0,
     .    ' and for l = ',ls0,' to ',' l = ',lq0
        end if
  120   continue
c
        end if
c
c       branch cut in l-direction
c
        if(ld.gt.2) then
c        
        do 130 j=1,jd,jd-1
        iface = 1
        if(j.eq.jd) iface = 2
        if(ibrnch(iface).gt.0) go to 130
        ls0 = 100000
        lq0 = 0
        ks0 = 100000
        kq0 = 0
        jcheck = 0
        do 230 k=1,kd-1
        do 230 ls=1,ld-2
        do 330 lq=ls+1,ld 
        delxyz = abs( x1(j,k,ls)-x1(j,k,lq) ) 
     .         + abs( y1(j,k,ls)-y1(j,k,lq) )
     .         + abs( z1(j,k,ls)-z1(j,k,lq) )
        if (delxyz.lt.eps) then
c          grid points ls and lq lie on opposite sides of a branch cut
           del2 = abs( x1(j,k,ls+1)-x1(j,k,lq-1) ) 
     .          + abs( y1(j,k,ls+1)-y1(j,k,lq-1) )
     .          + abs( z1(j,k,ls+1)-z1(j,k,lq-1) )
           if(del2.lt.eps)then
c            cell centers ls and lq-1 lie on opposite sides of a branch cut
             if (j.eq.1) then
                if(limage(1,k,ls).eq.ls) then
                jcheck = jcheck+1
                limage(1,k,ls)   = lq - 2
                limage(1,k,lq-2) = ls
                end if
             end if
             if (j.eq.jd) then
                if(limage(jd-1,k,ls).eq.ls) then
                jcheck = jcheck+1
                limage(jd-2,k,ls)   = lq - 2
                limage(jd-2,k,lq-2) = ls
                end if
             end if
             if(ls.lt.ls0)ls0 = ls
             if(lq.gt.lq0)lq0 = lq
             if(k.lt.ks0) ks0 = k
             if(k.gt.kq0) kq0 = k
           end if
        end if
  330   continue
  230   continue
c
c       amalgamate points to avoid gap in branch cut mapping  at tip...the
c       above criteria for C-type branch cut fails at the very tip where
c       the upper surface becomes the lower and vice versa.
c
        if(j.eq.1) jj = 1
        if(j.eq.jd) jj = jd-2
        do 231 k=1,kd-1
        do 331 l=2,ld-2
        if(limage(jj,k,l).eq.l) then
          if(limage(jj,k,l+1).ne.l+1.and.limage(jj,k,l-1).ne.l-1) then
            limage(jj,k,l) = l-1
c           write(6,*)' fixing hole in branch cut array at j,k,i = ',
c    .      j,k,l,' limage = ',limage(jj,k,l)
          end if
        end if
331     continue
231     continue
c
        if(jcheck.gt.0) then
          ibrnch(iface) = 30
          write(6,*)'      C-type branch cut detected on j = ',j,
     .    ' boundary'
          write(6,*)'        for k = ',ks0,' to k = ',kq0,
     .    ' and for l = ',ls0,' to ',' l = ',lq0
        end if
  130   continue
c
        end if
c
      end if
c      
c     check for C-grid type branch cut on l=1 and l=ld boundaries 
c
      if(ld.gt.2) then
c
c       branch cut in j-direction
c
        if(jd.gt.2) then
c        
        do 140 l=1,ld,ld-1
        iface = 5
        if(l.eq.ld) iface = 6
        js0 = 100000
        jq0 = 0
        ks0 = 100000
        kq0 = 0
        lcheck = 0
        do 240 k=1,kd-1
        do 240 js=1,jd-2
        do 340 jq=js+1,jd 
        delxyz = abs( x1(js,k,l)-x1(jq,k,l) ) 
     .         + abs( y1(js,k,l)-y1(jq,k,l) )
     .         + abs( z1(js,k,l)-z1(jq,k,l) )
        if (delxyz.lt.eps) then
c          grid points js and jq lie on opposite sides of a branch cut
           del2 = abs( x1(js+1,k,l)-x1(jq-1,k,l) ) 
     .          + abs( y1(js+1,k,l)-y1(jq-1,k,l) )
     .          + abs( z1(js+1,k,l)-z1(jq-1,k,l) )
           if(del2.lt.eps)then
c            cell centers js and jq-1 lie on opposite sides of a branch cut
             if (l.eq.1) then
                if(jimage(js,k,1).eq.js) then
                lcheck = lcheck+1
                jimage(js,k,1)   = jq - 2
                jimage(jq-2,k,1) = js
                end if
             end if
             if (l.eq.ld) then
                if(jimage(js,k,ld-1).eq.js) then
                lcheck = lcheck+1
                jimage(js,k,ld-2)   = jq - 2
                jimage(jq-2,k,ld-2) = js
                end if
             end if
             if(js.lt.js0)js0 = js
             if(jq.gt.jq0)jq0 = jq
             if(k.lt.ks0) ks0 = k
             if(k.gt.kq0) kq0 = k
           end if
        end if
  340   continue
  240   continue
c
c       amalgamate points to avoid gap in branch cut mapping at tip...the
c       above criteria for C-type branch cut fails at the very tip where
c       the upper surface becomes the lower and vice versa.
c
        if(l.eq.1) ll = 1
        if(l.eq.ld) ll = ld-2
        do 241 k=1,kd-1
        do 341 j=2,jd-2
        if(jimage(j,k,ll).eq.j) then
          if(jimage(j+1,k,ll).ne.j+1.and.jimage(j-1,k,ll).ne.j-1) then
            jimage(j,k,ll) = j-1
c           write(6,*)' fixing hole in branch cut array at j,k,i = ',
c    .      j,k,l,' jimage = ',jimage(j,k,ll)
          end if
        end if
341     continue
241     continue
c
        if(lcheck.gt.0) then
          ibrnch(iface) = 10
          write(6,*)'      C-type branch cut detected on l = ',l,
     .    ' boundary'
          write(6,*)'        for j = ',js0,' to j = ',jq0,
     .    ' and for k = ',ks0,' to ',' k = ',kq0
        end if
  140   continue
c
        end if
c
c       branch cut in k-direction
c
        if(kd.gt.2) then
c        
        do 150 l=1,ld,ld-1
        iface = 5
        if(l.eq.ld) iface = 6
        if(ibrnch(iface).gt.0) go to 150
        js0 = 100000
        jq0 = 0
        ks0 = 100000
        kq0 = 0
        lcheck = 0
        do 250 j=1,jd-1
        do 250 ks=1,kd-2
        do 350 kq=ks+1,kd 
        delxyz = abs( x1(j,ks,l)-x1(j,kq,l) ) 
     .         + abs( y1(j,ks,l)-y1(j,kq,l) )
     .         + abs( z1(j,ks,l)-z1(j,kq,l) )
        if (delxyz.lt.eps) then
c          grid points ks and kq lie on opposite sides of a branch cut
           del2 = abs( x1(j,ks+1,l)-x1(j,kq-1,l) ) 
     .          + abs( y1(j,ks+1,l)-y1(j,kq-1,l) )
     .          + abs( z1(j,ks+1,l)-z1(j,kq-1,l) )
           if(del2.lt.eps)then
c            cell centers ks and kq-1 lie on opposite sides of a branch cut
             if (l.eq.1) then
                if(kimage(j,ks,1).eq.ks) then
                lcheck = lcheck+1
                kimage(j,ks,1)   = kq - 2
                kimage(j,kq-2,1) = ks
                end if
             end if
             if (l.eq.ld) then
                if(kimage(j,ks,ld-1).eq.ks) then
                lcheck = lcheck+1
                kimage(j,ks,ld-2)   = kq - 2
                kimage(j,kq-2,ld-2) = ks
                end if
             end if
             if(ks.lt.ks0)ks0 = ks
             if(kq.gt.kq0)kq0 = kq
             if(j.lt.js0) js0 = j
             if(j.gt.jq0) jq0 = j
           end if
        end if
  350   continue
  250   continue
c
c       amalgamate points to avoid gap in branch cut mapping at tip...the
c       above criteria for C-type branch cut fails at the very tip where
c       the upper surface becomes the lower and vice versa.
c
        if(l.eq.1) ll = 1
        if(l.eq.ld) ll = ld-2
        do 251 j=1,jd-1
        do 351 k=2,kd-2
        if(kimage(j,k,ll).eq.k) then
          if(kimage(j,k+1,ll).ne.k+1.and.kimage(j,k-1,ll).ne.k-1) then
            kimage(j,k,ll) = k-1
c           write(6,*)' fixing hole in branch cut array at j,k,i = ',
c    .      j,k,l,' kimage = ',kimage(j,k,ll)
          end if
        end if
351     continue
251     continue
c
        if(lcheck.gt.0) then
          ibrnch(iface) = 20
          write(6,*)'      C-type branch cut detected on l = ',l,
     .    ' boundary'
          write(6,*)'        for j = ',js0,' to j = ',jq0,
     .    ' and for k = ',ks0,' to ',' k = ',kq0
        end if
  150   continue
c
        end if
c
      end if
c
c     iflag = 1 causes a message to be printed if some points are detected 
c     as being O-type branch points, but the cut does not extend over the 
c     entire face...this is usually a spurious message
c
      iflag = 0
c     
c     check for O-type branch cut between k=1 and k=kd boundaries
c
      if(kd.gt.2) then
c
      if(ibrnch(3).eq.0 .and. ibrnch(4).eq.0) then 
        jcount = 0
        do 400 l=1,ld
        do 400 j=1,jd
        delxyz = abs( x1(j,1,l)-x1(j,kd,l) ) 
     .         + abs( y1(j,1,l)-y1(j,kd,l) )
     .         + abs( z1(j,1,l)-z1(j,kd,l) )
        if (delxyz.lt.eps) jcount = jcount+1
  400   continue
c
        if(jcount.eq.jd*ld) then
c
          ibrnch(3) = 50 
          ibrnch(4) = 50
c 
          do 500 j=1,jd-1
          do 500 l=1,ld-1
          kimage(j,1,l)    = kd - 2
          kimage(j,kd-2,l) = 1
  500     continue
c
          write(6,*)'      O-type branch cut detected between k = 1',
     .    ' and k = ',kd,' boundaries'
          write(6,*)'        for l = 1 to l = ',ld,' and for j = 1',
     .    ' to j = ',jd
c
        else if(jcount.gt.0 .and. jcount.lt.jd*ld) then
          if(iflag.gt.0) then
          write(6,*)'      WARNING... ',jcount,' O-type branch cut',
     .    ' points detected between k = 1 and k = ',kd
          write(6,*)'      however, this is not the entire j and l',
     .    ' range as required'
          write(6,*)'      O-type image points not set...verify as',
     .    ' correct'  
          end if
        end if
      end if
c
      end if
c     
c     check for O-type branch cut between j=1 and j=jd boundaries
c
      if(jd.gt.2) then
c
      if(ibrnch(1).eq.0 .and. ibrnch(2).eq.0) then 
        kcount = 0
        do 410 l=1,ld
        do 410 k=1,kd
        delxyz = abs( x1(1,k,l)-x1(jd,k,l) ) 
     .         + abs( y1(1,k,l)-y1(jd,k,l) )
     .         + abs( z1(1,k,l)-z1(jd,k,l) )
        if (delxyz.lt.eps) kcount = kcount+1
  410   continue
c
        if(kcount.eq.kd*ld) then
c
          ibrnch(1) = 40 
          ibrnch(2) = 40 
c 
          do 510 k=1,kd-1
          do 510 l=1,ld-1
          jimage(1,k,l)    = jd - 2
          jimage(jd-2,k,l) = 1
  510     continue
c
          write(6,*)'      O-type branch cut detected between j = 1',
     .    ' and j = ',jd,' boundaries'
          write(6,*)'        for l = 1 to l = ',ld,' and for k = 1',
     .    ' to k = ',kd
c
        else if(kcount.gt.0 .and .kcount.lt.kd*ld) then
          if(iflag.gt.0) then
          write(6,*)'      WARNING... ',kcount,' O-type branch cut',
     .    ' points detected between j = 1 and j = ',jd
          write(6,*)'      however, this is not the entire k and l',
     .    ' range as required'
          write(6,*)'      O-type image points not set...verify as',
     .    ' correct'  
          end if
        end if      
      end if
c
      end if
c     
c     check for O-type branch cut between l=1 and l=ld boundaries
c
      if(ld.gt.2) then
c
      if(ibrnch(5).eq.0 .and. ibrnch(6).eq.0) then 
        lcount = 0 
        do 420 j=1,jd
        do 420 k=1,kd
        delxyz = abs( x1(j,k,1)-x1(j,k,ld) ) 
     .         + abs( y1(j,k,1)-y1(j,k,ld) )
     .         + abs( z1(j,k,1)-z1(j,k,ld) )
        if (delxyz.lt.eps) lcount = lcount+1
  420   continue
c
        if(lcount.eq.jd*kd) then
c
          ibrnch(5) = 60 
          ibrnch(6) = 60
c 
          do 520 j=1,jd-1
          do 520 k=1,kd-1
          limage(j,k,1)    = ld - 2
          limage(j,k,ld-2) = 1
  520     continue
c
          write(6,*)'      O-type branch cut detected between l = 1',
     .    ' and l = ',ld,' boundaries'
          write(6,*)'        for k = 1 to k = ',kd,' and for j = 1',
     .    ' to j = ',jd
c
        else if(lcount.gt.0 .and. lcount.lt.jd*kd) then
          if(iflag.gt.0) then
          write(6,*)'      WARNING... ',lcount,' O-type branch cut',
     .    ' points detected between l = 1 and l = ',ld
          write(6,*)'      however, this is not the entire j and k',
     .    ' range as required'
          write(6,*)'      O-type image points not set...verify as',
     .    ' correct'
          end if  
        end if      
      end if
c
      end if
c  
c     message for no branch cuts on block boundaries
c
      icount = 0
      do 600 iface=1,6
      if(ibrnch(iface).ne.0) icount = icount + 1
600   continue
      if(icount.eq.0) then
        write(6,*) '      no branch cuts detected in this mesh'
        return
      end if
      return
      end 
      subroutine bsurf(js,je,ks,ke,ls,le,jd,kd,ld,x,y,z,m1) 
c 
c*********************************************************************** 
c     Purpose: calculate normals to the coordinate surface defined by
c     js,je...ls,le in mesh m1. this coordinate surface and its normals 
c     are then loaded into 1-d boundary arrays.
c***********************************************************************
c 
c 
      include 'mag1.h'
c 
      common /surf/ xb(ibdim),yb(ibdim),zb(ibdim),ibmax
      common /norm/ vnx(ibdim),vny(ibdim),vnz(ibdim) 
      common /temp/ temp(ndim,ndim,6)
c
      dimension x(jd,kd,ld), y(jd,kd,ld), z(jd,kd,ld) 
c
      eps = 1.e-2
c
      if(ls.eq.le) then
      icase = 1
      sgn = float(ls)/float(iabs(ls))
      end if      
      if(ks.eq.ke) then
      icase = 2
      sgn = float(ks)/float(iabs(ks))
      end if      
      if(js.eq.je) then
      icase = 3
      sgn = float(js)/float(iabs(js))
      end if      
      ls = iabs(ls)
      le = iabs(le)
      ks = iabs(ks)
      ke = iabs(ke)
      js = iabs(js)
      je = iabs(je)
c 
c     l = constant plane
c
c
      if(icase.eq.1) then
c
      if((je-js).gt.ndim .or. (ke-ks).gt.ndim) then
        write(6,*)'  stopping...parameter ndim too small'
        stop
      end if
c
      l=ls
      jskip = je - js
      kskip = ke - ks
c
c     components of normals at cell centers stored in temp(i,j,1-3)
c     components of normals at grid point stored in temp(i,j,1-3)
c
      do 15 j = js,je-1
      do 15 k = ks,ke-1
c
      t1x = 0.5*(x(j+1,k+1,l) - x(j,k,l))
      t1y = 0.5*(y(j+1,k+1,l) - y(j,k,l))
      t1z = 0.5*(z(j+1,k+1,l) - z(j,k,l))
      t2x = 0.5*(x(j,k+1,l) - x(j+1,k,l))
      t2y = 0.5*(y(j,k+1,l) - y(j+1,k,l))
      t2z = 0.5*(z(j,k+1,l) - z(j+1,k,l))
c
      cnx = t1y*t2z - t2y*t1z
      cny = t2x*t1z - t1x*t2z
      cnz = t1x*t2y - t2x*t1y
c
      d=sqrt(cnx*cnx + cny*cny + cnz*cnz)
      if(d.gt.0) then
        temp(j,k,1) = sgn*cnx/d
        temp(j,k,2) = sgn*cny/d
        temp(j,k,3) = sgn*cnz/d
      else
        write(6,101) j,k,m1
101     format(3x,'singular normal at cell center j,k =',i3,i3,
     .  ' of mesh',i2)
        stop
      end if
15    continue
c     
c     average surrounding cell center normals to get grid point normals
c
c     interior points
c
      if(js+1.le.je-1 .and. ks+1.le.ke-1) then
        do 16 j = js+1,je-1
        do 16 k = ks+1,ke-1
        temp(j,k,4) = 0.25*(temp(j-1,k-1,1) + temp(j-1,k,1)
     .                    + temp(j,k,1) + temp(j,k-1,1))
        temp(j,k,5) = 0.25*(temp(j-1,k-1,2) + temp(j-1,k,2)
     .                    + temp(j,k,2) + temp(j,k-1,2))
        temp(j,k,6) = 0.25*(temp(j-1,k-1,3) + temp(j-1,k,3)
     .                    + temp(j,k,3) + temp(j,k-1,3))
16      continue
      end if
c
c     edge points
c
      if(js+1.le.je-1) then
        jskip = je - js
        kskip = ke - ks
        do 17 k = ks,ke,kskip
        k1 = k
        if(k.eq.ke) k1 = k-1
        do 17 j = js+1,je-1
        temp(j,k,4) = 0.5*(temp(j-1,k1,1) + temp(j,k1,1))
        temp(j,k,5) = 0.5*(temp(j-1,k1,2) + temp(j,k1,2))
        temp(j,k,6) = 0.5*(temp(j-1,k1,3) + temp(j,k1,3))
17      continue
      end if
c
      if(ks+1.le.ke-1) then 
        do 18 j = js,je,jskip
        j1 = j
        if(j.eq.je) j1 = j-1
        do 18 k = ks+1,ke-1
        temp(j,k,4) = 0.5*(temp(j1,k-1,1) + temp(j1,k,1))
        temp(j,k,5) = 0.5*(temp(j1,k-1,2) + temp(j1,k,2))
        temp(j,k,6) = 0.5*(temp(j1,k-1,3) + temp(j1,k,3))
18      continue
      end if 
c
c     corner points
c
      do 19 j = js,je,jskip
      jinc = 0
      if(j.eq.je)jinc = -1
      do 19 k = ks,ke,kskip
      kinc = 0
      if(k.eq.ke)kinc = -1
      temp(j,k,4) = temp(j+jinc,k+kinc,1)
      temp(j,k,5) = temp(j+jinc,k+kinc,2)
      temp(j,k,6) = temp(j+jinc,k+kinc,3)
19    continue
c
c     recess edge values to prevent double-valued corners
c
      do 170 j = js,je
      do 170 k = ks,ke
      temp(j,k,1) = x(j,k,l)
      temp(j,k,2) = y(j,k,l)
      temp(j,k,3) = z(j,k,l)
170   continue
      do 171 k = ks,ke,kskip
      kdir = 1
      if(k.eq.ke) kdir = -1
      do 171 j = js,je
      dx = temp(j,k+kdir,1) - temp(j,k,1)
      dy = temp(j,k+kdir,2) - temp(j,k,2)
      dz = temp(j,k+kdir,3) - temp(j,k,3)
      temp(j,k,1) = temp(j,k,1) + eps*dx
      temp(j,k,2) = temp(j,k,2) + eps*dy
      temp(j,k,3) = temp(j,k,3) + eps*dz
171   continue
      do 181 j = js,je,jskip
      jdir = 1
      if(j.eq.je) jdir = -1
      do 181 k = ks,ke
      dx = temp(j+jdir,k,1) - temp(j,k,1)
      dy = temp(j+jdir,k,2) - temp(j,k,2)
      dz = temp(j+jdir,k,3) - temp(j,k,3)
      temp(j,k,1) = temp(j,k,1) + eps*dx
      temp(j,k,2) = temp(j,k,2) + eps*dy
      temp(j,k,3) = temp(j,k,3) + eps*dz
181   continue
c
c     load boundary points and normals into 1-d arrays
c
      do 20 j = js,je
      do 20 k = ks,ke
      ibmax = ibmax + 1
      xb(ibmax)  = temp(j,k,1)
      yb(ibmax)  = temp(j,k,2)
      zb(ibmax)  = temp(j,k,3)
      vnx(ibmax) = temp(j,k,4)
      vny(ibmax) = temp(j,k,5)
      vnz(ibmax) = temp(j,k,6)
20    continue
c
c     k = constant plane
c
      else if(icase.eq.2) then
c
      if((je-js).gt.ndim .or. (le-ls).gt.ndim) then
        write(6,*)'  stopping...parameter ndim too small'
        stop
      end if
c
      k=ks
      lskip = le - ls
      jskip = je - js
c
c     components of normals at cell centers stored in temp(l,j,1-3)
c     components of normals at grid point stored in temp(l,j,1-3)
c
      do 25 l = ls,le-1
      do 25 j = js,je-1
c
      t1x = 0.5*(x(j+1,k,l+1) - x(j,k,l))
      t1y = 0.5*(y(j+1,k,l+1) - y(j,k,l))
      t1z = 0.5*(z(j+1,k,l+1) - z(j,k,l))
      t2x = 0.5*(x(j+1,k,l) - x(j,k,l+1))
      t2y = 0.5*(y(j+1,k,l) - y(j,k,l+1))
      t2z = 0.5*(z(j+1,k,l) - z(j,k,l+1))
c
      cnx = t1y*t2z - t2y*t1z
      cny = t2x*t1z - t1x*t2z
      cnz = t1x*t2y - t2x*t1y
c
      d=sqrt(cnx*cnx + cny*cny + cnz*cnz)

      if(d.gt.0) then
        temp(l,j,1) = sgn*cnx/d
        temp(l,j,2) = sgn*cny/d
        temp(l,j,3) = sgn*cnz/d

      else
        write(6,102) l,j,m1
102     format(3x,'singular normal at cell center l,j =',i3,i3,
     .  ' of mesh',i2)
        stop
      end if
25    continue
c     
c     average surrounding cell center normals to get grid point normals
c
c     interior points
c
      if(ls+1.le.le-1 .and. js+1.le.je-1) then
        do 26 l = ls+1,le-1
        do 26 j = js+1,je-1
        temp(l,j,4) = 0.25*(temp(l-1,j-1,1) + temp(l-1,j,1)
     .                    + temp(l,j,1) + temp(l,j-1,1))
        temp(l,j,5) = 0.25*(temp(l-1,j-1,2) + temp(l-1,j,2)
     .                    + temp(l,j,2) + temp(l,j-1,2))
        temp(l,j,6) = 0.25*(temp(l-1,j-1,3) + temp(l-1,j,3)
     .                    + temp(l,j,3) + temp(l,j-1,3))

26      continue
      end if
c
c     edge points
c
      if(ls+1.le.le-1) then
        lskip = le - ls
        jskip = je - js
        do 27 j = js,je,jskip
        j1 = j
        if(j.eq.je) j1 = j-1
        do 27 l = ls+1,le-1
        temp(l,j,4) = 0.5*(temp(l-1,j1,1) + temp(l,j1,1))
        temp(l,j,5) = 0.5*(temp(l-1,j1,2) + temp(l,j1,2))
        temp(l,j,6) = 0.5*(temp(l-1,j1,3) + temp(l,j1,3))
27      continue
      end if
      if(js+1.le.je-1) then
        do 28 l = ls,le,lskip
        l1 = l
        if(l.eq.le) l1 = l-1
        do 28 j = js+1,je-1
        temp(l,j,4) = 0.5*(temp(l1,j-1,1) + temp(l1,j,1))
        temp(l,j,5) = 0.5*(temp(l1,j-1,2) + temp(l1,j,2))
        temp(l,j,6) = 0.5*(temp(l1,j-1,3) + temp(l1,j,3))
28      continue
      end if
c
c     corner points
c
      do 29 l = ls,le,lskip
      linc = 0
      if(l.eq.le)linc = -1
      do 29 j = js,je,jskip
      jinc = 0
      if(j.eq.je)jinc = -1
      temp(l,j,4) = temp(l+linc,j+jinc,1)
      temp(l,j,5) = temp(l+linc,j+jinc,2)
      temp(l,j,6) = temp(l+linc,j+jinc,3)
29    continue
c
c     recess edge values to prevent double-valued corners
c
      do 270 l = ls,le
      do 270 j = js,je
      temp(l,j,1) = x(j,k,l)
      temp(l,j,2) = y(j,k,l)
      temp(l,j,3) = z(j,k,l)
270   continue
      do 271 j=js,je,jskip
      jdir = 1
      if(j.eq.je) jdir = -1
      do 271 l = ls,le
      dx = temp(l,j+jdir,1) - temp(l,j,1)
      dy = temp(l,j+jdir,2) - temp(l,j,2)
      dz = temp(l,j+jdir,3) - temp(l,j,3)
      temp(l,j,1) = temp(l,j,1) + eps*dx
      temp(l,j,2) = temp(l,j,2) + eps*dy
      temp(l,j,3) = temp(l,j,3) + eps*dz
271   continue
      do 281 l = ls,le,lskip
      ldir = 1
      if(l.eq.le) ldir = -1
      do 281 j = js,je
      dx = temp(l+ldir,j,1) - temp(l,j,1)
      dy = temp(l+ldir,j,2) - temp(l,j,2)
      dz = temp(l+ldir,j,3) - temp(l,j,3)
      temp(l,j,1) = temp(l,j,1) + eps*dx
      temp(l,j,2) = temp(l,j,2) + eps*dy
      temp(l,j,3) = temp(l,j,3) + eps*dz
281   continue
c
c     load boundary points and normals into 1-d arrays
c
      do 30 l = ls,le
      do 30 j = js,je
      ibmax = ibmax + 1
      xb(ibmax)  = temp(l,j,1)
      yb(ibmax)  = temp(l,j,2)
      zb(ibmax)  = temp(l,j,3)
      vnx(ibmax) = temp(l,j,4)
      vny(ibmax) = temp(l,j,5)
      vnz(ibmax) = temp(l,j,6)
30    continue
c
c     j = constant plane
c
      else if(icase.eq.3) then
c
      if((le-ls).gt.ndim .or. (ke-ks).gt.ndim) then
        write(6,*)'  stopping...parameter ndim too small'
        stop
      end if
c
      j=js
      kskip = ke - ks
      lskip = le - ls
c
c     components of normals at cell centers stored in temp(i,j,1-3)
c     components of normals at grid point stored in temp(i,j,1-3)
c
      do 35 k = ks,ke-1
      do 35 l = ls,le-1
c
      t1x = 0.5*(x(j,k+1,l+1) - x(j,k,l))
      t1y = 0.5*(y(j,k+1,l+1) - y(j,k,l))
      t1z = 0.5*(z(j,k+1,l+1) - z(j,k,l))
      t2x = 0.5*(x(j,k,l+1) - x(j,k+1,l))
      t2y = 0.5*(y(j,k,l+1) - y(j,k+1,l))
      t2z = 0.5*(z(j,k,l+1) - z(j,k+1,l))
c
      cnx = t1y*t2z - t2y*t1z
      cny = t2x*t1z - t1x*t2z
      cnz = t1x*t2y - t2x*t1y
c
      d=sqrt(cnx*cnx + cny*cny + cnz*cnz)
      if(d.gt.0) then
        temp(k,l,1) = sgn*cnx/d
        temp(k,l,2) = sgn*cny/d
        temp(k,l,3) = sgn*cnz/d
      else
        write(6,103) k,l,m1
103     format(3x,'singular normal at cell center k,l =',i3,i3,
     .  ' of mesh',i2)
        stop
      end if
35    continue
c     
c     average surrounding cell center normals to get grid point normals
c
c     interior points
c
      if(ks+1.le.ke-1 .and. ls+1.le.le-1) then
        do 36 k = ks+1,ke-1
        do 36 l = ls+1,le-1
        temp(k,l,4) = 0.25*(temp(k-1,l-1,1) + temp(k-1,l,1)
     .                    + temp(k,l,1) + temp(k,l-1,1))
        temp(k,l,5) = 0.25*(temp(k-1,l-1,2) + temp(k-1,l,2)
     .                    + temp(k,l,2) + temp(k,l-1,2))
        temp(k,l,6) = 0.25*(temp(k-1,l-1,3) + temp(k-1,l,3)
     .                    + temp(k,l,3) + temp(k,l-1,3))
36      continue
      end if
c
c     edge points
c
      if(ks+1.le.ke-1) then
        kskip = ke - ks
        lskip = le - ls
        do 37 l = ls,le,lskip
        l1 = l
        if(l.eq.le) l1 = l-1
        do 37 k = ks+1,ke-1
        temp(k,l,4) = 0.5*(temp(k-1,l1,1) + temp(k,l1,1))
        temp(k,l,5) = 0.5*(temp(k-1,l1,2) + temp(k,l1,2))
        temp(k,l,6) = 0.5*(temp(k-1,l1,3) + temp(k,l1,3))
37      continue
      end if
      if(ls+1.le.le-1) then
        do 38 k = ks,ke,kskip
        k1 = k
        if(k.eq.ke) k1 = k-1
        do 38 l = ls+1,le-1
        temp(k,l,4) = 0.5*(temp(k1,l-1,1) + temp(k1,l,1))
        temp(k,l,5) = 0.5*(temp(k1,l-1,2) + temp(k1,l,2))
        temp(k,l,6) = 0.5*(temp(k1,l-1,3) + temp(k1,l,3))
38      continue
      end if
c
c     corner points
c
      do 39 k = ks,ke,kskip
      kinc = 0
      if(k.eq.ke)kinc = -1
      do 39 l = ls,le,lskip
      linc = 0
      if(l.eq.le)linc = -1
      temp(k,l,4) = temp(k+kinc,l+linc,1)
      temp(k,l,5) = temp(k+kinc,l+linc,2)
      temp(k,l,6) = temp(k+kinc,l+linc,3)
39    continue
c
c     recess edge values to prevent double-valued corners
c
      do 370 k = ks,ke
      do 370 l = ls,le
      temp(k,l,1) = x(j,k,l)
      temp(k,l,2) = y(j,k,l)
      temp(k,l,3) = z(j,k,l)
370   continue
      do 371 l = ls,le,lskip
      ldir = 1
      if(l.eq.le) ldir = -1
      do 371 k = ks,ke
      dx = temp(k,l+ldir,1) - temp(k,l,1)
      dy = temp(k,l+ldir,2) - temp(k,l,2)
      dz = temp(k,l+ldir,3) - temp(k,l,3)
      temp(k,l,1) = temp(k,l,1) + eps*dx
      temp(k,l,2) = temp(k,l,2) + eps*dy
      temp(k,l,3) = temp(k,l,3) + eps*dz
371   continue
      do 381 k = ks,ke,kskip
      kdir = 1
      if(k.eq.ke) kdir = -1
      do 381 l = ls,le
      dx = temp(k+kdir,l,1) - temp(k,l,1)
      dy = temp(k+kdir,l,2) - temp(k,l,2)
      dz = temp(k+kdir,l,3) - temp(k,l,3)
      temp(k,l,1) = temp(k,l,1) + eps*dx
      temp(k,l,2) = temp(k,l,2) + eps*dy
      temp(k,l,3) = temp(k,l,3) + eps*dz
381   continue
c
c     load boundary  points and normals into 1-d arrays
c
      do 40 k = ks,ke
      do 40 l = ls,le
      ibmax = ibmax + 1
      xb(ibmax)  = temp(k,l,1)
      yb(ibmax)  = temp(k,l,2)
      zb(ibmax)  = temp(k,l,3)
      vnx(ibmax) = temp(k,l,4)
      vny(ibmax) = temp(k,l,5)
      vnz(ibmax) = temp(k,l,6)
40    continue
c
      end if
c 
c     make sure ibdim is large enough for this case
c
      if(ibmax.gt.ibdim) then
        write(6,*) '  stopping in bsurf...ibmax = ',ibmax,
     .  ' ibdim = ',ibdim
        write(6,*) '  increase parameter ibdim'
        stop
      end if
c
      return
      end
      subroutine cellcen(x,y,z,jd,kd,ld)
c
c***************************************************************** 
c     Purpose: find the cell centers of the grid 
c***************************************************************** 
c
      include 'mag1.h'
c 
      dimension x(jd,kd,ld),y(jd,kd,ld),z(jd,kd,ld)
      common/temp1/   x1(jdim,kdim,ldim),y1(jdim,kdim,ldim),
     .                z1(jdim,kdim,ldim) 
c 
      jmax = jd
      kmax = kd
      lmax = ld
c 
      do 10 l=1,lmax-1 
      do 20 k=1,kmax-1 
      do 30 j=1,jmax-1
      x1(j,k,l)=(x(j,k,l)+x(j+1,k,l)+x(j,k+1,l)+x(j+1,k+1,l)+
     .           x(j,k,l+1)+x(j,k+1,l+1)+x(j+1,k,l+1)+ 
     .           x(j+1,k+1,l+1))/8.0 
      y1(j,k,l)=(y(j,k,l)+y(j+1,k,l)+y(j,k+1,l)+y(j+1,k+1,l)+
     .           y(j,k,l+1)+y(j,k+1,l+1)+y(j+1,k,l+1)+ 
     .           y(j+1,k+1,l+1))/8.0 
      z1(j,k,l)=(z(j,k,l)+z(j+1,k,l)+z(j,k+1,l)+z(j+1,k+1,l)+
     .           z(j,k,l+1)+z(j,k+1,l+1)+z(j+1,k,l+1)+ 
     .           z(j+1,k+1,l+1))/8.0 
30    continue
20    continue 
10    continue 
c
      do 140 l=1,lmax-1
      do 150 k=1,kmax-1 
      do 160 j=1,jmax-1
      x(j,k,l) = x1(j,k,l) 
      y(j,k,l) = y1(j,k,l) 
      z(j,k,l) = z1(j,k,l) 
160   continue 
150   continue
140   continue 
c
      return 
      end
      subroutine chkstn( iblank ) 
c
c***********************************************************************
c     Purpose: the points in the interpolation/extrapolation stencil are 
c     checked to determine if illegal communication between 
c     mesh m and mesh m1 will result. there are four types of illegal
c     communication transfers:
c 
c       1. data transfer from hole or fringe point to fringe point 
c       2. data transfer from (interpolated/extrapolated) boundary point 
c          to fringe point 
c       3. data transfer from hole or fringe point to (interpolated/
c          extrapolated) boundary point 
c       4. from (interpolated/extrapolated) boundary point to (interpolated/
c          extrapolated) boundary point 
c
c     the flow solver must not be run if illegal transfers are detected and
c     cannot be corrected. provided that iholcor is set > 0 below, stencils 
c     that contain illegal points, but have at least one legal (i.e. field) 
c     point, are corrected using zeroth order interpolation.  this is  
c     accomplished by identifying the nearest legal (field) point in the 
c     stencil to the point to be interpolated.  This nearest legal point is 
c     then defined as the target point, and all weight of the interpolation 
c     is thrown into this target point, i.e., dxint=dyint=dzint=0
c
c     also: 
c       checks to make sure that stencils have been found
c       for all points which need them, and that all target points lie
c       in the (cell center) range 1 < ji/ki/li < j/k/lmax-1 (where
c       j/k/lmax are the dimensions of the cell center mesh)
c
c***********************************************************************
c 
      include 'mag1.h'
      parameter ( mdim2 = mdim*mdim ) 
c 
      logical iholcor 
c 
      common /book1/  ipntr(mdim,mhldim*mdim), npntr(mdim,mhldim*mdim),
     .                mhbs(mdim,mdim), mobs(mdim,mdim), nsets(mdim)
      common /book2/  dxint(idim), dyint(idim), dzint(idim),
     .                ibpts(mdim), jbpt(idim), kbpt(idim), lbpt(idim),
     .                ji(idim),ki(idim), li(idim)
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
      common /grid1 / x(mlen), y(mlen), z(mlen)
      common /image/  jimage(mlen), kimage(mlen), limage(mlen) 
      common /grid2/  mjmax(mdim), mkmax(mdim), mlmax(mdim)
      common /iwork/   ibck(iwrdim) 
      common /pltpt/  xbnd(idim),ybnd(idim),zbnd(idim),xorph(idim),
     .                yorph(idim),zorph(idim),xill(idim),yill(idim),
     .                zill(idim),nbnd,norph,nill
      common /surf/   xb(ibdim), yb(ibdim), zb(ibdim),ibmax
      common /chkst/ nwr 
c 
      dimension      iblank(*),list(mdim)
      dimension      mhhtrn(mdim,mdim), mhotrn(mdim,mdim),
     .               mohtrn(mdim,mdim), mootrn(mdim,mdim),
     .               ife(6) 
      dimension      jinew(idim),kinew(idim),linew(idim)
c 
c 
      data           mhhtrn / mdim2*0 /, mhotrn / mdim2*0 /,    
     .               mohtrn / mdim2*0 /, mootrn / mdim2*0 / 
c
c     perform correction to stencils containing hole,fringe or boundary points
c     (if possible) when iholcor is true
c
      iholcor = .true.
c 
      iflag1 = 0
      iflag2 = 0
c
      do ll=1,6
         ife(ll) = 0
      end do
c
      do 31 m = 1,nmesh  
c 
c***************************************************************************
c     check to determine if any stencils for mesh m fringe points contain 
c     illegal points
c***************************************************************************
c
      nh = nhole(m) 
c
      if (nh.gt.0) then
        jmaxm = mjmax(m)
        kmaxm = mkmax(m)
        lmaxm = mlmax(m)
        do 14 n = 1,nh
        nserch = mhole(m,n)
        do 141 nn = 1,nserch
        call getgrd( m,x,y,z,jimage,kimage,limage,jmaxm,kmaxm,lmaxm )
        ioh = 0
        ihh = 0
        m1  = lhole(m,n,nn) 
        call getint( m1,ji,ki,li,jbpt,kbpt,lbpt,
     .       dxint,dyint,dzint )
        iset = mhbs(m,m1) 
        is   = ipntr(m1,iset) 
        ie   = npntr(m1,iset) 
        ic   = ie - is + 1 
        if(ic.gt.0) then 
c
c         make sure stencils of fringe points have been set and that all target 
c         points lie in the range 1 < ji,ki,li < j/k/lmax-1 
c         (where j/k/lmax is the dimension of cell center grid). 
c
          jmax = mjmax(m1)
          kmax = mkmax(m1)
          lmax = mlmax(m1)
          jmaxx = jmax - 1
          kmaxx = kmax - 1
          lmaxx = lmax - 1
          if(jmax.lt.2) jmaxx = 1
          if(kmax.lt.2) kmaxx = 1
          if(lmax.lt.2) lmaxx = 1
          do 1600 i=is,ie
c         points with zeroth order interpolation can actually use the range
c         1 < j/k/lmax instead of the range 1 < j,k,lmax-1 for points with
c         first order interpolation/extrapolation
          if(dxint(i).eq.0. .and. dyint(i).eq.0. .and. 
     .    dzint(i).eq.0.) then
            jmaxx = jmax
            kmaxx = kmax
            lmaxx = lmax
          end if
          if(ji(i).lt.1 .or. ji(i).gt.jmaxx .or. 
     .       ki(i).lt.1 .or. ki(i).gt.kmaxx .or.
     .       li(i).lt.1 .or. li(i).gt.lmaxx) then
            write(6,*) ' stopping...illegal target cell for fringe',
     .      ' point ',jbpt(i),kbpt(i),lbpt(i),' of mesh ',m
            write(6,*) '   target cell is ',ji(i),ki(i),li(i),
     .      ' in mesh ',m1
            write(6,*) '   legal range in mesh ',m1,' is:',
     .      '  1<j<',jmaxx,'  1<k<',kmaxx,'  1<l<',lmaxx
            write(6,*) '   possible miscue in search routine'
            stop
          end if 
1600      continue
c
          do 131 i = 1,ic
          ii     = i + is-1 
          call bound(jbpt(ii),kbpt(ii),lbpt(ii),x0,y0,z0,jmaxm,kmaxm,
     .    lmaxm,x,y,z,ife)
          xb(i) = x0
          yb(i) = y0
          zb(i) = z0
131       continue
          call getibl( m1,iblank,jmax,kmax,lmax ) 
          call getgrd( m1,x,y,z,jimage,kimage,limage,jmax,kmax,lmax )
          call ibchk( ibck,iblank,ic,is,ji,ki,li,m1,jmax,kmax,lmax,
     .    xb,yb,zb,jinew,kinew,linew,x,y,z )
          do 10 i = 1,ic
          if( ibck(i).gt.0 ) then 
c
c           stencil contains interpolated boundary points
c
            ioh    = ioh + 1
            iflag1 = iflag1 + 1
            ii     = i + is-1 
            nill = nill + 1
            xill(nill) = xb(i)
            yill(nill) = yb(i)
            zill(nill) = zb(i)
c
            write(8,702)jbpt(ii),kbpt(ii),lbpt(ii),m,ji(ii),ki(ii),
     .      li(ii),m1
 702        format(' the fringe point jb,kb,lb ',3i4,' of mesh ',
     .      i3,/,'    has a stencil from point ji,ki,li ',3i4,' which',
     .      ' lies on an interp bndry of mesh ',i3)
c
c           if possible (i.e. target cell itself is not an interpolated 
c           boundary  point), use zeroth order transfer to correct this 
c           illegal stencil.
c
            if ( iholcor ) then
              if(ibck(i) .gt. 1) then
c               zeroth order transfer possible
                write(8,*) '      this illegal transfer was',
     .          ' circumvented using zeroth order data transfer' 
                write(8,*) '        old ji,ki,li: ',ji(ii),ki(ii),li(ii)
                write(8,*) '        new ji,ki,li: ',jinew(ii),kinew(ii),
     .          linew(ii)
                dxint(ii) = 0.
                dyint(ii) = 0.
                dzint(ii) = 0.
                ji(ii) = jinew(ii)
                ki(ii) = kinew(ii)
                li(ii) = linew(ii)
                iflag2    = iflag2 + 1
              end if
            end if
          else if( ibck(i).lt.0 ) then 
c
c           stencil contains hole and/or fringe points
c
            ihh    = ihh + 1
            iflag1 = iflag1 + 1
            ii     = i + is-1 
            nill = nill + 1
            xill(nill) = xb(i)
            yill(nill) = yb(i)
            zill(nill) = zb(i)
c
            write(8,802)jbpt(ii),kbpt(ii),lbpt(ii),m,ji(ii),ki(ii),
     .      li(ii),m1
 802        format(' the fringe point jb,kb,lb ',3i4,' of mesh ',
     .      i3,/,'    has a stencil point ji,ki,li ',3i4,' which lies',
     .      ' in a hole of mesh ',i3)
c
c           if possible (i.e. target cell itself is not a hole point), use
c           zeroth order transfer to correct this illegal stencil.
c
            if ( iholcor ) then
              if(ibck(i) .lt. -1) then
c               zeroth order transfer possible
                write(8,*) '      this illegal transfer was',
     .          ' circumvented using zeroth order data transfer' 
                write(8,*) '        old ji,ki,li: ',ji(ii),ki(ii),li(ii)
                write(8,*) '        new ji,ki,li: ',jinew(ii),kinew(ii),
     .          linew(ii)
                dxint(ii) = 0.
                dyint(ii) = 0.
                dzint(ii) = 0.
                ji(ii) = jinew(ii)
                ki(ii) = kinew(ii)
                li(ii) = linew(ii)
                iflag2    = iflag2 + 1
              end if
            end if
          end if
   10     continue
          call putint(nwr,m1,ji,ki,li,jbpt,kbpt,lbpt,dxint,dyint,dzint) 
        end if
        mhhtrn(m,m1) = ihh + mhhtrn(m,m1)
        mohtrn(m,m1) = ioh + mohtrn(m,m1)
  141   continue
  14    continue
      end if
c***************************************************************************
c     check to determine if any stencils for mesh m boundary points contain 
c     illegal points
c***************************************************************************
c
      no = noutr(m) 
c 
      if (no.gt.0) then
        jmaxm = mjmax(m)
        kmaxm = mkmax(m)
        lmaxm = mlmax(m)
        do 24 n = 1,no
        nnn = 1
        nserch=moutr(m,nnn)
        do 241 nn=1,nserch
        call getgrd( m,x,y,z,jimage,kimage,limage,jmaxm,kmaxm,lmaxm )
        ioo = 0
        iho = 0
        nnn = 1
        m1 = loutr(m,nnn,nn)
        call getint( m1,ji,ki,li,jbpt,kbpt,lbpt,
     .       dxint,dyint,dzint )
        iset = mobs(m,m1) 
        is   = ipntr(m1,iset) 
        ie   = npntr(m1,iset) 
        ic   = ie - is + 1 
        if(ic.gt.0) then 
c
c         make sure stencils of boundary points have been set and that all 
c         target points lie in the range 1 < ji,ki,li < j/k/lmax-1 
c         (where j/k/lmax is the dimension of cell center grid). 
c
          jmax = mjmax(m1)
          kmax = mkmax(m1)
          lmax = mlmax(m1)
          jmaxx = jmax - 1
          kmaxx = kmax - 1
          lmaxx = lmax - 1
          if(jmax.lt.2) jmaxx = 1
          if(kmax.lt.2) kmaxx = 1
          if(lmax.lt.2) lmaxx = 1
          do 2600 i=is,ie
c         points with zeroth order interpolation can actually use the range
c         1 < j/k/lmax instead of the range 1 < j,k,lmax-1 for points with
c         first order interpolation/extrapolation
          if(dxint(i).eq.0. .and. dyint(i).eq.0. .and. 
     .    dzint(i).eq.0.) then
            jmaxx = jmax
            kmaxx = kmax
            lmaxx = lmax
          end if
          if(ji(i).lt.1 .or. ji(i).gt.jmaxx .or. 
     .       ki(i).lt.1 .or. ki(i).gt.kmaxx .or.
     .       li(i).lt.1 .or. li(i).gt.lmaxx) then
            write(6,*) ' stopping...illegal target cell for boundary',
     .      ' point ',jbpt(i),kbpt(i),lbpt(i),' of mesh ',m
            write(6,*) '   target cell is ',ji(i),ki(i),li(i),
     .      ' in mesh ',m1
            write(6,*) '   legal range in mesh ',m1,' is:',
     .      '  1<j<',jmaxx,'  1<k<',kmaxx,'  1<l<',lmaxx
            write(6,*) '   possible miscue in search routine'
            stop
          end if 
2600      continue
c
          do 231 i = 1,ic
          ii     = i + is-1 
          call bound(jbpt(ii),kbpt(ii),lbpt(ii),x0,y0,z0,jmaxm,kmaxm,
     .    lmaxm,x,y,z,ife)
          xb(i) = x0
          yb(i) = y0
          zb(i) = z0
231       continue
          call getibl( m1,iblank,jmax,kmax,lmax ) 
          call getgrd( m1,x,y,z,jimage,kimage,limage,jmax,kmax,lmax )
          call ibchk( ibck,iblank,ic,is,ji,ki,li,m1,jmax,kmax,lmax,
     .    xb,yb,zb,jinew,kinew,linew,x,y,z )
          do 20 i = 1,ic
          if( ibck(i).gt.0 ) then 
c
c           stencil contains interpolated boundary points
c
            ioo    = ioo + 1
            iflag1 = iflag1 + 1
            ii     = i + is-1 
            nill = nill + 1
            xill(nill) = xb(i)
            yill(nill) = yb(i)
            zill(nill) = zb(i)
c
            write(8,502)jbpt(ii),kbpt(ii),lbpt(ii),m,ji(ii),ki(ii),
     .      li(ii),m1
 502        format(' the boundary point jb,kb,lb ',3i4,' of mesh ',
     .      i3,/,'    has a stencil from point ji,ki,li ',3i4,' which',
     .      ' lies on an interp bndry of mesh ',i3)
c
c           if possible (i.e. target cell itself is not an interpolated 
c           boundary  point), use zeroth order transfer to correct this 
c           illegal stencil.
c
            if ( iholcor ) then
              if(ibck(i) .gt. 1) then
c               zeroth order transfer possible
                write(8,*) '      this illegal transfer was',
     .          ' circumvented using zeroth order data transfer' 
                write(8,*) '        old ji,ki,li: ',ji(ii),ki(ii),li(ii)
                write(8,*) '        new ji,ki,li: ',jinew(ii),kinew(ii),
     .          linew(ii)
                dxint(ii) = 0.
                dyint(ii) = 0.
                dzint(ii) = 0.
                ji(ii) = jinew(ii)
                ki(ii) = kinew(ii)
                li(ii) = linew(ii)
                iflag2    = iflag2 + 1
              end if
            end if
          else if( ibck(i).lt.0 ) then 
c
c           stencil contains hole and/or fringe points
c
            iho    = iho + 1
            iflag1 = iflag1 + 1
            ii     = i + is-1 
            nill = nill + 1
            xill(nill) = xb(i)
            yill(nill) = yb(i)
            zill(nill) = zb(i)
c
            write(8,602)jbpt(ii),kbpt(ii),lbpt(ii),m,ji(ii),ki(ii),
     .      li(ii),m1
 602        format(' the boundary point jb,kb,lb ',3i4,' of mesh ',
     .      i3,/,'    has a stencil point ji,ki,li ',3i4,' which lies',
     .      ' in a hole of mesh ',i3)
c
c           if possible (i.e. target cell itself is not a hole point), use
c           zeroth order transfer to correct this illegal stencil.
c
            if ( iholcor ) then
              if(ibck(i) .lt. -1) then
c               zeroth order transfer possible
                write(8,*) '      this illegal transfer was',
     .          ' circumvented using zeroth order data transfer' 
                write(8,*) '        old ji,ki,li: ',ji(ii),ki(ii),li(ii)
                write(8,*) '        new ji,ki,li: ',jinew(ii),kinew(ii),
     .          linew(ii)
                dxint(ii) = 0.
                dyint(ii) = 0.
                dzint(ii) = 0.
                ji(ii) = jinew(ii)
                ki(ii) = kinew(ii)
                li(ii) = linew(ii)
                iflag2    = iflag2 + 1
              end if
            end if
          end if
   20     continue
          call putint(nwr,m1,ji,ki,li,jbpt,kbpt,lbpt,dxint,dyint,dzint) 
        end if
        mhotrn(m,m1) = iho + mhotrn(m,m1)
        mootrn(m,m1) = ioo + mootrn(m,m1)
  241   continue
  24    continue
      end if
   31 continue
c 
c***********************************************************************
c            output summary of stencil checks
c***********************************************************************
c
      write(6,*)
      write(6,*) '    checking for four types of illegal',
     .' transfer of information:'
      write(6,*) '    1. fringe points receive data from hole/fringe',
     .' points'
      write(6,*) '    2. fringe points receive data from boundary',
     .' points'
      write(6,*) '    3. boundary points receive data from hole/fringe',
     .' points'
      write(6,*) '    4. boundary pts receive data from boundary points'      
c
      do 401 m=1,nmesh
      write(6,*)
      write(6,*) '    checking stencils for mesh ',m
c
      write(6,*) '      checking stencils of fringe points'
      nh = nhole(m)
      do 512 m1=1,nmesh
      list(m1) = 0
512   continue
      if(nh .le. 0) then
        write(6,*)'        there are no fringe pts in mesh ',m,
     .  ' to check'
      else
c
c       generate list of all meshes which got searched to find the fringe 
c       points of all holes in the current mesh
c
        do 612 n = 1,nh 
        nserch = mhole(m,n)
        do 612 nn=1,nserch
        m1 = lhole(m,n,nn) 
        if(list(m1).eq.0) list(m1) = 1
612     continue
c
        do 501 m1=1,nmesh
        if(list(m1).gt.0) then
          if(nhole(m1).gt.0) then
            write(6,*)'        1. ',mhhtrn(m,m1),' fringe pts in',
     .      ' mesh ',m,' receive data from hole/fringe pts in mesh ',m1
          else
            write(6,*)'        1. not applicable, no hole/fringe pts',
     .      ' in mesh ',m1
          end if
        end if
501     continue
        do 701 m1=1,nmesh
        if(list(m1).gt.0) then
          if(noutr(m1).gt.0) then
            write(6,*)'        2. ',mohtrn(m,m1),' fringe pts in',
     .      ' mesh ',m,' receive data from boundary pts in mesh ',m1
          else
            write(6,*)'        2. not applicable, no boundary pts in',
     .     ' mesh ',m1
          end if
        end if
701     continue
      end if
c
      write(6,*) '      checking stencils of boundary points'
      no = noutr(m)
      do 712 m1=1,nmesh
      list(m1) = 0
712   continue
      if(no .le. 0) then
        write(6,*)'        there are no boundary pts in mesh ',m,
     .  ' to check'
      else
        do 812 n = 1,no 
        nnn = 1
        nserch=moutr(m,nnn)
        do 812 nn=1,nserch
        nnn = 1
        m1 = loutr(m,nnn,nn)
        if(list(m1).eq.0) list(m1) = 1
812     continue
        do 801 m1=1,nmesh
        if(list(m1).gt.0) then
          if(nhole(m1).gt.0) then
            write(6,*)'        3. ',mhotrn(m,m1),' bndry pts in mesh ',
     .      m,' receive data from hole/fringe pts in mesh ',m1
          else
            write(6,*)'        3. not applicable, no hole/fringe pts',
     .      ' in mesh ',m1
          end if
        end if
801     continue
        do 901 m1=1,nmesh
        if(list(m1).gt.0) then
          if(noutr(m1).gt.0) then
            write(6,*)'        4. ',mootrn(m,m1),' bndry pts in mesh ',
     .      m,' receive data from bndry pts in mesh ',m1
          else
            write(6,*)'        4. not applicable, no boundary pts in',
     .      ' mesh ',m1
          end if
        end if
901     continue
      end if
401   continue
c
c
      if(iflag1.gt.0) then
        write(6,*)
        write(6,*) '    ***CAUTION***CAUTION***CAUTION***CAUTION',
     .  '***CAUTION***CAUTION***CAUTION***'
        write(6,*) '    of the ',iflag1,' illegal stencils',
     .  ' detected above, ',iflag2,' were circumvented'
        write(6,*) '    by  using zeroth order data transfer.',
     .  ' illegal transfer of information'
        write(6,*) '    will no longer occur with these stencils,',
     .  ' but accuracy and convergence'
        write(6,*) '    may suffer. a large geometric mismatch', 
     .  ' (see below, and also file 8)'
        write(6,*) '    indicates loss of accuracy. if possible,',
     .  ' adjust maggie input parameters'
        write(6,*) '    and/or component grids to reduce',
     .  ' the number of stencils needing correction'
        write(6,*) '    ***CAUTION***CAUTION***CAUTION***CAUTION',
     .  '***CAUTION***CAUTION***CAUTION***'
      end if 
c
      if(iflag2.ne.iflag1) then
        write(6,*)
        write(6,*) '    ***WARNINIG***WARNING***WARNING***WARNING',
     .  '***WARNING***WARNING***WARNING***'
        write(6,*) '                           DO NOT RUN FLOW',
     .  ' SOLVER!!!'
        write(6,*) '    ',iflag1-iflag2,' illegal stencils',
     .  ' could not be circumvented...must either adjust' 
        write(6,*) '    maggie input parameters and/or component grids',
     .  ' to rectify this problem.'
        write(6,*) '                           DO NOT RUN FLOW',
     .  ' SOLVER!!!'
        write(6,*) '    ***WARNINIG***WARNING***WARNING***WARNING',
     .  '***WARNING***WARNING***WARNING***'
      end if
c 
      return
      end
      subroutine chkxyz
c
c***********************************************************************
c     Purpose: check for mismatch between 1) x,y,z of the cell centers for 
c     fringe and boundary points as determined from the input grids and 
c     2) x,y,z of the same points as found by trilinear interpolation using 
c     the transfer stencil information calculated by maggie. Substantial
c     mismatch indicates an error in the transfer stencil.
c 
c***********************************************************************
c 
      include 'mag1.h'
c 
      common /book1/  ipntr(mdim,mhldim*mdim), npntr(mdim,mhldim*mdim),
     .                mhbs(mdim,mdim), mobs(mdim,mdim), nsets(mdim)
      common /book2/  dxint(idim), dyint(idim), dzint(idim),
     .                ibpts(mdim), jbpt(idim), kbpt(idim), lbpt(idim),
     .                ji(idim),ki(idim), li(idim)
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
      common /grid1 / x(mlen), y(mlen), z(mlen)
      common /grid2/  mjmax(mdim), mkmax(mdim), mlmax(mdim)
      common /chkpt/  xcc0(idim),ycc0(idim),zcc0(idim),xccint(idim),
     .                yccint(idim),zccint(idim),
     .                xscale(idim),yscale(idim),zscale(idim)
      common /image/  jimage(mlen), kimage(mlen), limage(mlen) 
c
      dimension ife(6)
c
      do ll=1,6
         ife(ll) = 0
      end do
c
      write(6,*)
      write(6,*) '    checking for mismatch between cell',
     .' centers of fringe and boundary pts'
      write(6,*) '    as calculated from 1) input meshes and 2)',
     .' data transfer stencils' 
c
      sigdif = .1
c
      do 1000 m=1,nmesh
c
      itotm = 0
      itotn = 0
      errmax = 0.
      jerrmx = 1
      kerrmx = 1
      lerrmx = 1
      nsig   = 0
      nh     = nhole(m) 
      no     = noutr(m) 
      jmax   = mjmax(m)
      kmax   = mkmax(m)
      lmax   = mlmax(m)
c
c     index increments for 2d or 3d case
c  
      jinc = 1
      kinc = 1
      linc = 1
      if(jmax.lt.2) jinc = 0
      if(kmax.lt.2) kinc = 0
      if(lmax.lt.2) linc = 0
c 
c     check fringe points
c
      if(nh.gt.0) then
        do 14 n = 1,nh
        nserch = mhole(m,n)
        do 14 nn = 1,nserch
        m1  = lhole(m,n,nn)
        jmax1 = mjmax(m1)
        kmax1 = mkmax(m1)
        lmax1 = mlmax(m1)
        call getint( m1,ji,ki,li,jbpt,kbpt,lbpt,dxint,dyint,dzint )
        iset = mhbs(m,m1) 
        is   = ipntr(m1,iset) 
        ie   = npntr(m1,iset)  
        call getgrd( m,x,y,z,jimage,kimage,limage,jmax,kmax,lmax )
c
        do 141 i=is,ie 
c
c       cell center value of point jbpt,kbpt,lbpt of mesh m from cell center
c       grid
c
        call bound(jbpt(i),kbpt(i),lbpt(i),x0,y0,z0,jmax,kmax,lmax,
     .  x,y,z,ife)
c
c       set x,y,z length scales from cell with jbpt,kbpt,lbpt as lower 
c       corner point (for extrapolated points, i.e. those that don't exist
c       in the input grid, use a nearby neighbor that does exist in the
c       input grid to set the length scales)
c
        j00 = jbpt(i)
        k00 = kbpt(i)
        l00 = lbpt(i)
        if(j00.ge.jmax) j00 = jmax - 1
        if(k00.ge.kmax) k00 = kmax - 1
        if(l00.ge.lmax) l00 = lmax - 1
        if(j00.le.0) j00 = 1
        if(k00.le.0) k00 = 1
        if(l00.le.0) l00 = 1
c
        xscal = 0.
        yscal = 0.
        zscal = 0.
        do 10 j=j00,j00+jinc
        do 10 k=k00,k00+kinc
        do 10 l=l00,l00+linc
        ii = j +(k-1)*jmax +(l-1)*jmax*kmax 
        dx = abs(x(ii) - x0)
        dy = abs(y(ii) - y0)
        dz = abs(z(ii) - z0)
        if(dx.gt.xscal) xscal = dx
        if(dy.gt.yscal) yscal = dy
        if(dz.gt.zscal) zscal = dz
10      continue
c
        if(xscal.le.0.) xscal = 1.
        if(yscal.le.0.) yscal = 1.
        if(zscal.le.0.) zscal = 1.
c
        itotn = itotn + 1
        xcc0(itotn) = x0 
        ycc0(itotn) = y0 
        zcc0(itotn) = z0 
        xscale(itotn) = xscal
        yscale(itotn) = yscal
        zscale(itotn) = zscal
c
141     continue
c
        call getgrd( m1,x,y,z,jimage,kimage,limage,
     .  jmax1,kmax1,lmax1 ) 
c
        do 13 i = is,ie
c
c       calculate the cell centers using interpolation stencils
c
        call interp(jmax1,kmax1,lmax1,x,y,z,ji(i),ki(i),li(i),
     .              dxint(i),dyint(i),dzint(i),xintrp,yintrp,
     .              zintrp)
c
        itotm = itotm + 1
        x0    = xcc0(itotm)
        y0    = ycc0(itotm)
        z0    = zcc0(itotm)
        xscal = xscale(itotm)
        yscal = yscale(itotm)
        zscal = zscale(itotm)
        xccint(itotm) = xintrp
        yccint(itotm) = yintrp
        zccint(itotm) = zintrp
c
c       compare the two cell center coordinates.  a sigdif difference
c       (scaled with cell size xscal,yscal,zscal) is deemed significant
c
        xerr  = abs(xintrp - x0)/xscal
        yerr  = abs(yintrp - y0)/yscal
        zerr  = abs(zintrp - z0)/zscal
c
	error = xerr
        if(yerr .gt. error) error = yerr
        if(zerr .gt. error) error = zerr
        if(error .gt. errmax) then
          errmax = error
          jerrmx = jbpt(i)
          kerrmx = kbpt(i)
          lerrmx = lbpt(i)
        end if
c
        if(xerr.gt.sigdif .or. yerr.gt.sigdif .or. zerr.gt.sigdif) then
          nsig = nsig + 1
          write(8,*)'interpolated coordinates of fringe pt j,k,i = ',
     .    jbpt(i),kbpt(i),lbpt(i),' mesh = ',m
          write(8,*)'differ by ',100*error,' percent',
     .    ' from the actual coordinates:'
          write(8,*)'        x0,y0,z0 = ',x0,y0,z0
          write(8,*)'  xint,yint,zint = ',xintrp,yintrp,zintrp
          write(8,*)'  xsca,ysca,zsca = ',xscal,yscal,zscal
          write(8,*)'  target cell: j,k,l = ',ji(i),ki(i),li(i),
     .    ' mesh = ',m1
          write(8,*)'  xie,eta,zeta = ',dxint(i),dyint(i),dzint(i)
        end if
c
   13   continue
   14   continue
      end if
c 
c     check boundary points
c
      if(no.gt.0) then
        do 16 n = 1,no
        nnn = 1
        nserch=moutr(m,nnn)
        do 16 nn = 1,nserch
        nnn = 1
        m1 = loutr(m,nnn,nn)
        jmax1 = mjmax(m1)
        kmax1 = mkmax(m1)
        lmax1 = mlmax(m1)
        call getint( m1,ji,ki,li,jbpt,kbpt,lbpt,dxint,dyint,dzint )
        iset = mobs(m,m1) 
        is   = ipntr(m1,iset) 
        ie   = npntr(m1,iset)  
        call getgrd( m,x,y,z,jimage,kimage,limage,jmax,kmax,lmax )
c
        do 161 i=is,ie 
c
c       cell center value of point jbpt,kbpt,lbpt of mesh m from cell center
c       grid
c
        call bound(jbpt(i),kbpt(i),lbpt(i),x0,y0,z0,jmax,kmax,lmax,
     .  x,y,z,ife)
c
c       set x,y,z length scales from cell with jbpt,kbpt,lbpt as lower 
c       corner point (for extrapolated points, i.e. those that don't exist
c       in the input grid, use a nearby neighbor that does exist in the
c       input grid to set the length scales)
c
        j00 = jbpt(i)
        k00 = kbpt(i)
        l00 = lbpt(i)
        if(j00.ge.jmax) j00 = jmax - 1
        if(k00.ge.kmax) k00 = kmax - 1
        if(l00.ge.lmax) l00 = lmax - 1
        if(j00.le.0) j00 = 1
        if(k00.le.0) k00 = 1
        if(l00.le.0) l00 = 1
c
        xscal = 0.
        yscal = 0.
        zscal = 0.
        do 11 j=j00,j00+jinc
        do 11 k=k00,k00+kinc
        do 11 l=l00,l00+linc
        ii = j +(k-1)*jmax +(l-1)*jmax*kmax 
        dx = abs(x(ii) - x0)
        dy = abs(y(ii) - y0)
        dz = abs(z(ii) - z0)
        if(dx.gt.xscal) xscal = dx
        if(dy.gt.yscal) yscal = dy
        if(dz.gt.zscal) zscal = dz
11      continue
c
        if(xscal.le.0.) xscal = 1.
        if(yscal.le.0.) yscal = 1.
        if(zscal.le.0.) zscal = 1.
c
        itotn = itotn + 1
        xcc0(itotn) = x0 
        ycc0(itotn) = y0 
        zcc0(itotn) = z0 
        xscale(itotn) = xscal
        yscale(itotn) = yscal
        zscale(itotn) = zscal
c
161     continue
c
        call getgrd( m1,x,y,z,jimage,kimage,limage,
     .  jmax1,kmax1,lmax1 ) 
c
        do 15 i = is,ie
c
c       calculate the cell centers using interpolation stencils
c
        call interp(jmax1,kmax1,lmax1,x,y,z,ji(i),ki(i),li(i),
     .              dxint(i),dyint(i),dzint(i),xintrp,yintrp,
     .              zintrp)
c
        itotm = itotm + 1
        x0    = xcc0(itotm)
        y0    = ycc0(itotm)
        z0    = zcc0(itotm)
        xscal = xscale(itotm)
        yscal = yscale(itotm)
        zscal = zscale(itotm)
        xccint(itotm) = xintrp
        yccint(itotm) = yintrp
        zccint(itotm) = zintrp
c
c       compare the two cell center coordinates.  a sigdif difference
c       (scaled with cell size xscal,yscal,zscal) is deemed significant
c
        xerr  = abs(xintrp - x0)/xscal
        yerr  = abs(yintrp - y0)/yscal
        zerr  = abs(zintrp - z0)/zscal
c
	error = xerr
        if(yerr .gt. error) error = yerr
        if(zerr .gt. error) error = zerr
        if(error .gt. errmax) then
          errmax = error
          jerrmx = jbpt(i)
          kerrmx = kbpt(i)
          lerrmx = lbpt(i)
        end if
c
        if(xerr.gt.sigdif .or. yerr.gt.sigdif .or. zerr.gt.sigdif) then
          nsig = nsig + 1
          write(8,*)'interpolated coordinates of boundary pt j,k,i = ',
     .    jbpt(i),kbpt(i),lbpt(i),' mesh = ',m
          write(8,*)'differ by ',100*error,' percent',
     .    ' from the actual coordinates:'
          write(8,*)'        x0,y0,z0 = ',x0,y0,z0
          write(8,*)'  xint,yint,zint = ',xintrp,yintrp,zintrp
          write(8,*)'  xsca,ysca,zsca = ',xscal,yscal,zscal
          write(8,*)'  target cell: j,k,l = ',ji(i),ki(i),li(i),
     .    ' mesh = ',m1
          write(8,*)'  xie,eta,zeta = ',dxint(i),dyint(i),dzint(i)
        end if
c
   15   continue
   16   continue
      end if
c
      if(nh+no .gt.0) then
        write(6,*)
        write(6,*) '      ',nsig,' cell centers of mesh ',m,
     .  ' calculated using the transfer stencils'
        write(6,*) '      differ by more than ',int(sigdif*100),
     .  ' percent from the actual cell centers'  
        write(6,*) '      the maximum deviation ',100*errmax,
     .  ' percent occurs at ',jerrmx,kerrmx,lerrmx 
c
      end if
c
1000  continue
      return
      end
      subroutine cindex (ibc,ibpnts,iipnts,jb,kb,lb,mesh )
c 
c***********************************************************************
c     Purpose: the flow solver cfl3d expects interpolation data in the
c     form of list of jb,kb,lb for a grid at a time.thus, 
c     construct the cross index array ibc and pointers for transfer
c     of interpolated values to the appropriate boundary points. 
c     these points are located in different meshes from those in 
c     which they were interpolated: thus the need for a cross index. 
c***********************************************************************
c 
      include 'mag1.h'
c 
      common /book1/  ipntr(mdim,mhldim*mdim), npntr(mdim,mhldim*mdim),
     .                mhbs(mdim,mdim), mobs(mdim,mdim), nsets(mdim)
      common /book2/  dxint(idim), dyint(idim), dzint(idim),
     .                ibpts(mdim), jbpt(idim), kbpt(idim), lbpt(idim),
     .                ji(idim),ki(idim), li(idim)
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
c 
      dimension ipts(mdim),list(mdim)
      dimension ibc(*),jb(*),kb(*),lb(*)
c  
c     set up offset counter 
c 
      ipts(1) = 0
      do 11 m = 2,nmesh 
      ipts(m) = ipts(m-1) + ibpts(m-1) 
   11 continue
c 
c     set up cross index and pointers 
c 
      m  = mesh 
      ic = 0
c 
c     outer boundaries 
c 
      no = noutr(m) 
c 
      if( no.ne.0 ) then
        do 312 m1=1,nmesh
        list(m1) = 0
312     continue
        do 412 n = 1,no 
        nnn = 1
        nserch=moutr(m,nnn)
        do 412 nn=1,nserch
        nnn = 1
        m1 = loutr(m,nnn,nn)
        if(list(m1).eq.0) list(m1) = 1
412     continue
        do 23 m1 = 1,nmesh
        if(list(m1) .gt. 0) then
          iset       = mobs(m,m1) 
          is         = ipntr(m1,iset) 
          ie         = npntr(m1,iset) 
          call getint( m1,ji,ki,li,jbpt,kbpt,lbpt,
     .                 dxint,dyint,dzint )
          do 22 i = is,ie 
          ic         = ic + 1 
          ibc(ic)    = i + ipts(m1) 
          jb(ic)     = jbpt(i)
          kb(ic)     = kbpt(i)
          lb(ic)     = lbpt(i)
   22     continue
        end if
   23   continue
      end if
c 
c     hole boundary (fringe point)
c 
      nh = nhole(m) 
c 
      if( nh.ne.0 ) then
        do 512 m1=1,nmesh
        list(m1) = 0
512     continue
        do 612 n = 1,nh 
        nserch = mhole(m,n)
        do 612 nn=1,nserch
        m1 = lhole(m,n,nn) 
        if(list(m1).eq.0) list(m1) = 1
612     continue
        do 32 m1 = 1,nmesh
        if(list(m1).gt.0) then
          iset       = mhbs(m,m1) 
          is         = ipntr(m1,iset) 
          ie         = npntr(m1,iset) 
          call getint( m1,ji,ki,li,jbpt,kbpt,lbpt,
     .                 dxint,dyint,dzint )
          do 31 i = is,ie 
          ic      = ic + 1 
          ibc(ic) = i + ipts(m1) 
          jb(ic)  = jbpt(i)
          kb(ic)  = kbpt(i)
          lb(ic)  = lbpt(i)
   31     continue
        end if
   32   continue
      end if
c 
c     set pointers 
c 
      ibpnts = ic 
      iipnts = ibpts(m) 
c 
c     check to make sure that ic is within dimensions of arrays
c 
      if( ibpnts.gt.idim ) then 
        write(6,601) mesh,ibpnts,idim 
  601   format( /,/,' ',10x,'failure in cindex the number of boundary',
     .  1x,'points in mesh ',i3,1x,'is ibpnts = ',i5, 
     .  1x,'and the dimension is idim = ',i5  ) 
        stop 'cindex' 
      end if
c 
      return
      end 
      subroutine diagnos 
c 
c***********************************************************************
c     Purpose: perform diagnostic checks
c***********************************************************************
c 
      include 'mag1.h'
c 
      character*80 grid,outpt,ovrlp,plt3d
c
      common /book3/  iblank(mlen) 
      common /grid1/  x(mlen), y(mlen), z(mlen)
      common /grid2/  mjmax(mdim), mkmax(mdim), mlmax(mdim)
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
      common /image/  jimage(mlen), kimage(mlen), limage(mlen) 
      common /pltpt/  xbnd(idim),ybnd(idim),zbnd(idim),xorph(idim),
     .                yorph(idim),zorph(idim),xill(idim),yill(idim),
     .                zill(idim),nbnd,norph,nill
      common /diagno/ jp(mdim*3+3),kp(mdim*3+3),lp(mdim*3+3),iplt3d
      common /files/  grid,outpt,ovrlp,plt3d
c 
      write(6,*)
      write(6,*) '  *** beginning diagnostic checks of interpolation',
     .'/extrapolation stencils ***'
c 
c***********************************************************************
c    check for hole points or interpolated points as part of interpolation 
c    stencil
c***********************************************************************
c 
      call chkstn( iblank ) 
c
c***********************************************************************
c     check for mismatch between 1) x,y,z of the cell centers for fringe 
c     and boundary points as determined from the input grids and 2) x,y,z
c     of the same points as found by trilinear interpolation using the 
c     interpolation stencil information calculated by maggie. Substantial
c     mismatch indicates an inappropriate interpolation stencil.
c***********************************************************************
c
      call chkxyz
c
c***********************************************************************
c     output plot3d file to unit 9 to provide a visual check of point
c     type classification (i.e. field, hole fringe, boundary) as well 
c     as points for which extrapolation rather than interpolation was
c     used (i.e. orphan points) and points having illegal interpolation
c     stencils. 
c***********************************************************************
c
      if(iplt3d.gt.0) then
        write(6,*)
        write(6,'(''  plot3d diagnostic file (mg/blank) being'',
     .  '' written to '',a60)')plt3d
        write(6,*)'        the file has the following',
     .  ' structure:'    
        nmesh3=nmesh*3
        mm = 0
        do 3 m=1,nmesh3,3
        mm = mm+1
        write(6,*)
        write(6,101) m,mm
        write(6,102) m+1,mm
        write(6,103) m+2,mm
3       continue
        if(nbnd.gt.0)then  
          write(6,104) nmesh3+1
          write(6,1012)
        end if
        if(norph.gt.0) then
          if(nbnd.gt.0) then
            write(6,105) nmesh3+2
          else
            write(6,105) nmesh3+1
          end if
          write(6,1012)
        end if
        if(nill.gt.0) then
          if(nbnd.gt.0 .and. norph.gt.0) then
            write(6,106) nmesh3+3
          else if(nbnd.gt.0 .or. norph.gt.0) then
            write(6,106) nmesh3+2
          else
            write(6,106) nmesh3+1
          end if
          write(6,1013)
          write(6,1012)
        end if
101     format(8x,'plot3d grid ',i2,'.....mesh ',i2,' field points')
102     format(8x,'plot3d grid ',i2,'.....mesh ',i2,' hole points')
103     format(8x,'plot3d grid ',i2,'.....mesh ',i2,' fringe points')
104     format(/,8x,'plot3d grid ',i2,'.....interpolated',
     .  ' boundary points (all meshes)')
105     format(/,8x,'plot3d grid ',i2,'.....extrapolated points',
     .  ' (orphans) (all meshes)')
106     format(/,8x,'plot3d grid ',i2,'.....points with illegal',
     .  ' stencils (all meshes)')
1012    format(22x,'.....must plot as POINTS in plot3d')
1013    format(22x,'.....(does not reflect any illegal stencils',
     .  /,22x,'.....corrected by zeroth order data transfer)')
        do 5 m=1,nmesh
        jd=mjmax(m) 
        kd=mkmax(m) 
        ld=mlmax(m) 
        call getgrd( m,x,y,z,jimage,kimage,limage,jd,kd,ld ) 
        call getibl( m,iblank,jd,kd,ld ) 
        call pltpts( nmesh,m,jd,kd,ld,x,y,z,jimage,iblank )
5       continue
      end if
c
      return
      end
      subroutine dsmin(jmax,kmax,lmax,x,y,z,xp,yp,zp,jp,kp,lp,js,
     .                 je,ks,ke,ls,le,jskip,kskip,lskip,dmin)
c
c***********************************************************************
c      Purpose: Find closest point in grid to point (xp,yp,zp) to serve 
c      a staring point for the search routine.
c      note that in many cases, the cell associated with the minimum 
c      distance point is not the cell which actually contains xp,yp,zp;
c      a better esimate of the cell that surrounds xp,yp,zp is therefore
c      made. 
c***********************************************************************
c
      implicit real (a-h,o-z)
c
      include 'mag1.h'
c
      common/work/ rm(iwrdim) 
c
      dimension x(jmax,kmax,lmax),y(jmax,kmax,lmax),
     .          z(jmax,kmax,lmax)
c
      dmin = 1.0e+20
c
      jp = 0
      kp = 0	
      lp = 0
      num = (je - js)/jskip + 1
      if( num.gt.iwrdim ) then
        write( 6,601 ) num,iwrdim 
601     format( ' ', 10x,'failure in dsmin: num exceeds iwrdim',
     .  1x,'num = ',i5,1x,'iwrdim = ',i5  ) 
        stop 'dsmin' 
      end if 
c
      do 1234 k=ks,ke,kskip
      do 1234 l=ls,le,lskip
c
      jj = 0
      do 1235 j=js,je,jskip
      jj = jj + 1
      dx = xp-x(j,k,l)
      dy = yp-y(j,k,l)
      dz = zp-z(j,k,l)
      rm(jj) = dx*dx + dy*dy + dz*dz
 1235 continue

      jmin = ismin( num,rm,1 ) 
      d1   = rm(jmin)
      if (d1.lt.dmin) then
         jp   = (jmin-1)*jskip + js 
         kp   = k
         lp   = l
         dmin = d1
      end if
 1234 continue
c
      return
      end
      subroutine dsmin2(jmax,kmax,lmax,x,y,z,xp,yp,zp,jp,kp,lp,js,
     .                 je,ks,ke,ls,le,dmin)
c
c***********************************************************************
c      Purpose: Find closest point in grid faces to point (xp,yp,zp) to 
c      serve a staring point for the search routine when extrapolating 
c      for orphan points.
c***********************************************************************
c
      implicit real (a-h,o-z)
c
      include 'mag1.h'
c
      common/work/ rm(iwrdim) 
c
      dimension x(jmax,kmax,lmax),y(jmax,kmax,lmax),
     .          z(jmax,kmax,lmax)
c
      dmin = 1.0e+20
c
c     search only faces of mesh, using original grid density 
c
      jskip = je - js 
      kskip = ke - ks
      lskip = le - ls
      if(jskip.lt.1)jskip = 1
      if(kskip.lt.1)kskip = 1
      if(lskip.lt.1)lskip = 1
c
      num = je - js + 1
      if( num.gt.iwrdim ) then
        write( 6,601 ) num,iwrdim 
601     format( ' ', 10x,'failure in dsmin: num exceeds iwrdim',
     .  1x,'num = ',i5,1x,'iwrdim = ',i5  ) 
        stop 'dsmin2' 
      end if 
c
c     k-faces
c
      do 2234 k=ks,ke,kskip
      do 2234 l=ls,le
c
      jj = 0
      do 2235 j=js,je
      jj = jj + 1
      dx = xp-x(j,k,l)
      dy = yp-y(j,k,l)
      dz = zp-z(j,k,l)
      rm(jj) = dx*dx + dy*dy + dz*dz
 2235 continue

      jmin = ismin( num,rm,1 ) 
      d1   = rm(jmin)
      if (d1.lt.dmin) then
         jp   = jmin + js - 1 
         kp   = k
         lp   = l
         dmin = d1
      end if
 2234 continue
c
c     l-faces
c
      do 3234 l=ls,le,lskip
      do 3234 k=ks,ke
c
      jj = 0
      do 3235 j=js,je
      jj = jj + 1
      dx = xp-x(j,k,l)
      dy = yp-y(j,k,l)
      dz = zp-z(j,k,l)
      rm(jj) = dx*dx + dy*dy + dz*dz
 3235 continue
      jmin = ismin( num,rm,1 ) 
      d1   = rm(jmin)
      if (d1.lt.dmin) then
         jp   = jmin + js - 1 
         kp   = k
         lp   = l
         dmin = d1
      end if
 3234 continue
c
c     j-faces
c
      num = ke - ks + 1
      if( num.gt.iwrdim ) then
        write( 6,601 ) num,iwrdim 
        stop 'dsmin2' 
      end if 
      do 4234 j=js,je,jskip
      do 4234 l=ls,le
c
      kk = 0
      do 4235 k=ks,ke
      kk = kk + 1
      dx = xp-x(j,k,l)
      dy = yp-y(j,k,l)
      dz = zp-z(j,k,l)
      rm(kk) = dx*dx + dy*dy + dz*dz
 4235 continue
      kmin = ismin( num,rm,1 ) 
      d1   = rm(kmin)
      if (d1.lt.dmin) then
         kp   = kmin + ks - 1 
         jp   = j
         lp   = l
         dmin = d1
      end if
 4234 continue
c
      return
      end
      subroutine extrap(i,m,m1,jd1,kd1,ld1,xm1,ym1,zm1,xp,yp,zp,iok,
     .              jp,kp,lp,jimage,kimage,limage)     
c
c******************************************************************
c
c     Purpose: search over mesh m1 to find the cell nearest the 
c     mesh m point xp,yp,zp and compute the extrapolation data
c******************************************************************
c
      include 'mag1.h'
c
c     parameter mfroz should be set at least as large as itmax
c
      parameter (mfroz=999)
c
      common /intrp1/ xi(idim), yi(idim), zi(idim) 
      common /intrp2/ jb(idim), kb(idim), lb(idim), jn(idim), kn(idim),
     .                ln(idim), itotal 
      common /tol/    epsc
c 
      dimension xm1(jd1,kd1,ld1),ym1(jd1,kd1,ld1),zm1(jd1,kd1,ld1)
      dimension jimage(jd1,kd1,ld1),kimage(jd1,kd1,ld1),
     .          limage(jd1,kd1,ld1)
      dimension jfroz(mfroz),kfroz(mfroz),lfroz(mfroz)
      dimension xfroz(mfroz),yfroz(mfroz),zfroz(mfroz)
c
      idum1 = 0
      idum2 = 0
      idum3 = 0
      idum4 = 0
      dum1  = 0.
      dum2  = 0.
      dum3  = 0.
c
      limit  = 1
      itmax  = 200
      iok    = 0
      ijump  = 0
      xie    = .5
      eta    = .5
      zeta   = .5
c***
c      itrace = 1
c***
c
      call trace(11,i,idum2,idum3,idum4,xp,yp,zp)
c
      do 5555 intern=1,itmax
c
      jfroz(intern)  = jp
      kfroz(intern)  = kp
      lfroz(intern)  = lp
c
      call trace(3,intern,idum2,idum3,idum4,dum1,dum2,dum3)
c
c     find local xie,eta,zeta via Newton iteraton in current target 
c     cell jp,kp,lp
c
      call trace(4,jp,kp,lp,m1,dum1,dum2,dum3)
c
      call xe(jd1,kd1,ld1,xm1,ym1,zm1,jp,kp,lp,xp,yp,zp, 
     .              xie,eta,zeta,imiss)
c
      xfroz(intern) = xie
      yfroz(intern) = eta
      zfroz(intern) = zeta
c
c     current target cell correct if imiss = 0 (i.e. interpolation will
c     be done, rather than extrapolation)
c
      if(imiss.eq.0) then 
        call trace(5,idum1,idum2,idum3,idum4,xie,eta,zeta)
c
        go to 5556
      end if
c
c     current block m1 a poor guess if xie,eta or zeta become large. exit
c     to try another mesh m1
c
      huge = 1.e5
      if(abs(xie).ge.huge .or. abs(eta).ge.huge .or. 
     .abs(zeta).ge.huge) then
        iok = 0
        return
      end if
 
c
      call trace(5,idum1,idum2,idum3,idum4,xie,eta,zeta)
c
c     update current guess for target cell based on result of Newton
c     iteration, with max allowable change set by limit
c
      jsav = jp
      ksav = kp
      lsav = lp
c
      if (xie.ge.0) jinc = abs(xie)
      if (xie.lt.0) jinc = abs(xie-1)
      if (eta.ge.0) kinc = abs(eta)
      if (eta.lt.0) kinc = abs(eta-1)
      if (zeta.ge.0) linc = abs(zeta)
      if (zeta.lt.0) linc = abs(zeta-1)
c
      jinc = min0( jinc , limit )
      kinc = min0( kinc , limit )
      linc = min0( linc , limit )
c
      if (xie.gt.1.0) then
         jp = jp + jinc  
      else if (xie.lt.0.) then
         jp = jp - jinc  
      end if
      if (eta.gt.1.0) then
         kp = kp + kinc
      else if (eta.lt.0.) then
         kp = kp - kinc 
      end if
      if (zeta.gt.1.0) then
         lp = lp + linc
      else if (zeta.lt.0.) then
         lp = lp - linc 
      end if
c
      xieg  = float(jp)
      etag  = float(kp)
      zetag = float(lp)
c
c     keep within bounds of mesh m1
c
      jp = min0( jp , jd1-1 )
      kp = min0( kp , kd1-1 )
      lp = min0( lp , ld1-1 )
      jp = max0( 1 , jp )
      kp = max0( 1 , kp )
      lp = max0( 1 , lp )
c
c     account for any branch cuts
c
      jpc = jp
      kpc = kp
      lpc = lp
c
      if(xieg .lt. 1. .or. xieg .gt.jd1-1) then
        jpc = jimage(jp,kp,lp)
        kpc = kimage(jp,kp,lp)
        lpc = limage(jp,kp,lp)
      end if
      if(etag .lt. 1. .or. etag .gt.kd1-1) then
        jpc = jimage(jp,kp,lp)
        kpc = kimage(jp,kp,lp)
        lpc = limage(jp,kp,lp)
      end if
      if(zetag.lt. 1. .or. zetag.gt.ld1-1) then
        jpc = jimage(jp,kp,lp)
        kpc = kimage(jp,kp,lp)
        lpc = limage(jp,kp,lp)
      end if
c
      if(jpc.ne.jp .or. kpc.ne.kp .or. lpc.ne.lp) then
        ijump = ijump +1
c       allow only two jumps across branch cut
        if(ijump.le.2) then
          call trace(23,m1,jpc,kpc,lpc,dum1,dum2,dum3)
          call trace(24,idum1,jp,kp,lp,dum1,dum2,dum3)
          jp = jpc
          kp = kpc
          go to 5555
        end if
      end if
c
c     check for frozen convergence: search routine keeps returning
c     to the same point, without 0 < xie,eta,zeta < 1 at that point.
c     of those cells searched in the frozen cycle, use the point with
c     with the minimum extrapolation coefficient as the point to 
c     extrapolate from, unless that requires extrapolation in two
c     or three directions.
c     note: extrapolation measured from 0.0 for negative values and 
c     1.0 for positive values.
c
      ifroz = 0
      do 77 ii=1,intern
      int = intern-ii+1
      if (jp.eq.jfroz(int) .and. kp.eq.kfroz(int)
     .    .and. lfroz(int).eq.lp) ifroz = 1 
77    continue
      if(ifroz.gt.0) then
c       frozen convervence...search is cyclic
        xi0 = xfroz(1)
        yi0 = yfroz(1) 
        zi0 = zfroz(1)
        iext = 1
        nex0 = 0
c
        if(xi0-1.0 .gt. epsc) then
          nex0 = nex0 + 1
          ximod0 = abs(xi0-1.0)
        else if(xi0 .lt. -epsc) then
          nex0 = nex0 + 1
          ximod0 = abs(xi0)
        else
          ximod0 = 0.
        end if
        if(yi0-1.0 .gt. epsc) then
          nex0 = nex0 + 1
          yimod0 = abs(yi0-1.0)
        else if(yi0 .lt. -epsc) then
          nex0 = nex0 + 1
          yimod0 = abs(yi0)
        else
          yimod0 = 0.
        end if
        if(zi0-1.0 .gt. epsc) then
          nex0 = nex0 + 1
          zimod0 = abs(zi0-1.0)
        else if(zi0 .lt. -epsc) then
          nex0 = nex0 + 1
          zimod0 = abs(zi0)
        else
          zimod0 = 0.
        end if
        exmax0 = max(ximod0,yimod0,zimod0)
c
        do 78 ii=2,intern
        xi1 = xfroz(ii)
        yi1 = yfroz(ii)
        zi1 = zfroz(ii)
        nex1 = 0
        if(xi1-1.0 .gt. epsc) then
          nex1 = nex1 + 1
          ximod1 = abs(xi1-1.0)
        else if(xi1 .lt. -epsc) then
          nex1 = nex1 + 1
          ximod1 = abs(xi1)
        else
          ximod1 = 0.
        end if
        if(yi1-1.0 .gt. epsc) then
          nex1 = nex1 + 1
          yimod1 = abs(yi1-1.0)
        else if(yi1 .lt. -epsc) then
          nex1 = nex1 + 1
          yimod1 = abs(yi1)
        else
          yimod1 = 0.
        end if
        if(zi1-1.0 .gt. epsc) then
          nex1 = nex1 + 1
          zimod1 = abs(zi1-1.0)
        else if(zi1 .lt. -epsc) then
          nex1 = nex1 + 1
          zimod1 = abs(zi1)
        else
          zimod1 = 0.
        end if
        exmax1 = max(ximod1,yimod1,zimod1)
c
c       prefer the stencil with extrapolation in fewest directions
        if(nex1.lt.nex0) then
          nex0 = nex1
          ximod0 = ximod1
          yimod0 = yimod1
          zimod0 = zimod1
          exmax0 = exmax1
          iext = ii
c       prefer the stencil with smaller extrapolation coefficient, provided
c       it does not extrapolate in more directions than the previous one
        else if(exmax1.lt.exmax0 .and. nex1.eq.nex0) then
          nex0 = nex1
          ximod0 = ximod1
          yimod0 = yimod1
          zimod0 = zimod1
          exmax0 = exmax1
          iext = ii
        end if
c
78      continue
c
        jp   = jfroz(iext)
        kp   = kfroz(iext)
        lp   = lfroz(iext)
        xie  = xfroz(iext)
        eta  = yfroz(iext)
        zeta = zfroz(iext)
        call trace(17,idum1,idum2,idum3,idum4,dum1,dum2,dum3)
        call trace(4,jp,kp,lp,m1,dum1,dum2,dum3)
        call trace(5,idum1,idum2,idum3,idum4,xie,eta,zeta)
        go to 5556
      end if
c
5555  continue
5556  continue
c
      jn(i)  = jp
      kn(i)  = kp
      ln(i)  = lp
      xi(i)  = xie 
      yi(i)  = eta
      zi(i)  = zeta
      iok = 1
c
c***
c      itrace = -1
c***     
      return
      end     
      subroutine frnge( iblank,jd,kd,ld,m1 )
c 
c***********************************************************************
c     Purpose: construct the  fringe boundary surrounding a hole.
c     specify fringe points by setting 
c     iblank = -( m1 + mdim ). Check for two points on either side
c     of a field point to give enough data for second order accurate 
c     solution
c***********************************************************************
c 
      include 'mag1.h'
c 
      dimension iblank(jd,kd,ld)
c 
      jmax       = jd 
      kmax       = kd 
      lmax       = ld 
c 
c     loop through cell centers
c
      do 11 j = 1,jmax
      do 11 k = 1,kmax
      do 11 l = 1,lmax
      if( iblank(j,k,l).eq.1 ) then 
c
c       point is field point. determine if it is adjacent to a 
c       hole point and hence should be identified as a fringe point
c 
c       keep within grid boundaries
c 
        jp   = min( j+1,jmax )
        kp   = min( k+1,kmax )
        lp   = min( l+1,lmax )
        jm   = max( j-1,1 ) 
        km   = max( k-1,1 ) 
        lm   = max( l-1,1 ) 
        jp1  = min( j+2,jmax )
        kp1  = min( k+2,kmax )
        lp1  = min( l+2,lmax )
        jm1  = max( j-2,1 ) 
        km1  = max( k-2,1 ) 
        lm1  = max( l-2,1 )
        imin = min( 
     .         abs( iblank(jm1,k,l) ),  abs( iblank(jm,k,l) ),
     .         abs( iblank(jp1,k,l) ),  abs( iblank(jp,k,l) ),
     .         abs( iblank(j,km1,l) ),  abs( iblank(j,km,l) ),
     .         abs( iblank(j,kp1,l) ),  abs( iblank(j,kp,l) ),
     .         abs( iblank(j,k,lm1) ),  abs( iblank(j,k,lm) ),
     .         abs( iblank(j,k,lp1) ),  abs( iblank(j,k,lp) ) )
        if( imin.eq.0 ) then
c         fringe point
          iblank(j,k,l) = -( m1 +mdim ) 
        end if
      end if
11    continue
      return
      end 
      subroutine getgrd( m,x,y,z,jimage,kimage,limage,jd,kd,ld )
c
c***********************************************************************
c     Purpose: get a copy of cell-center grid for mesh m from
c     the file temp_cen.m, along with branch-cut arrays
c***********************************************************************
c
      character*20 titl
c
      dimension x(jd,kd,ld), y(jd,kd,ld), z(jd,kd,ld) 
      dimension jimage(jd,kd,ld),kimage(jd,kd,ld),limage(jd,kd,ld)
c
      iunit = 30
      if (m.gt.99) then
          len = 12
          write (titl,'("temp_cen.",i3)') m
      else if (m.gt.9) then
          len = 11
          write (titl,'("temp_cen.",i2)') m
      else
          len = 10
          write (titl,'("temp_cen.",i1)') m
      endif
      do i = len+1, 20
          titl(i:i) = ' '
      end do
      open(iunit,file=titl(1:len),form='unformatted',
     .status='unknown')
c
      read(iunit) x,y,z 
      read(iunit) jimage,kimage,limage
      rewind iunit
      close(iunit)
      return
      end 
      subroutine getgrd2( m,x,y,z,jd,kd,ld ) 
c
c***********************************************************************
c     Purpose: get a copy of grid points (node points) for mesh m from
c     file temp_grd.m
c***********************************************************************
c
      character*20 titl
c
      dimension x(jd,kd,ld), y(jd,kd,ld), z(jd,kd,ld) 
c
      iunit = 30 
      if (m.gt.99) then
          len = 12
          write (titl,'("temp_grd.",i3)') m
      else if (m.gt.9) then
          len = 11
          write (titl,'("temp_grd.",i2)') m
      else
          len = 10
          write (titl,'("temp_grd.",i1)') m
      endif
      do i = len+1, 20
          titl(i:i) = ' '
      end do
      open(iunit,file=titl(1:len),form='unformatted',
     .status='unknown')
c 
      read(iunit) x,y,z
      rewind iunit
      close(iunit)
c 
      return
      end 
      subroutine getibl( m,iblank,jd,kd,ld ) 
c
c***********************************************************************
c     Purpose: get a copy of the iblank array for mesh m from
c     file temp_ibl.m
c***********************************************************************
c
      character*20 titl
c
      dimension iblank(jd,kd,ld)
c
      iunit = 30
      if (m.gt.99) then
          len = 12
          write (titl,'("temp_ibl.",i3)') m
      else if (m.gt.9) then
          len = 11
          write (titl,'("temp_ibl.",i2)') m
      else
          len = 10
          write (titl,'("temp_ibl.",i1)') m
      endif
      do i = len+1, 20
          titl(i:i) = ' '
      end do
      open(iunit,file=titl(1:len),form='unformatted',
     .status='unknown')
      read(iunit) iblank
      rewind iunit
      close(iunit)
      return
      end 
      subroutine getint( m,ji,ki,li,jbpt,kbpt,lbpt, 
     &                   dxint,dyint,dzint )
c 
c***********************************************************************
c     Purpose: get a copy of the interpolation data for mesh m from
c     file temp_int.m
c***********************************************************************
c 
      include 'mag1.h'
c
      character*20 titl
c 
      dimension ji(idim),ki(idim),li(idim), 
     .          jbpt(idim),kbpt(idim),lbpt(idim), 
     .          dxint(idim),dyint(idim),dzint(idim)
c
      common /chkst/ nwr 
c
      iunit = 30
      if (m.gt.99) then
          len = 12
          write (titl,'("temp_int.",i3)') m
      else if (m.gt.9) then
          len = 11
          write (titl,'("temp_int.",i2)') m
      else
          len = 10
          write (titl,'("temp_int.",i1)') m
      endif
      do i = len+1, 20
          titl(i:i) = ' '
      end do
      open(iunit,file=titl(1:len),form='unformatted',
     .status='unknown')
      read(iunit) nwr
      if(nwr.gt.0) then
        read(iunit) (ji(i),ki(i),li(i),i=1,nwr),
     .              (jbpt(i),kbpt(i),lbpt(i),i=1,nwr),
     .              (dxint(i),dyint(i),dzint(i),i=1,nwr)
      end if
      rewind iunit
      close(iunit)
      return
      end 
      subroutine hlocat( iblank,x,y,z,jd,kd,ld)
c 
c***********************************************************************
c     Purpose: locate the points interior to the hole boundary in 
c     mesh m due to the grid m1. 
c
c     interior points are identified by forming the dot product, d = (rp,vn) 
c     where rp is the position vector between the nearest point on the 
c     boundary surface to a field point of m, and where vn is the outward 
c     normal vector to the surface. d < 0 indicates an interior point
c***********************************************************************
c 
c 
      include 'mag1.h'
c 
      common /norm/   vnx(ibdim), vny(ibdim),vnz(ibdim) 
      common /surf/   xb(ibdim), yb(ibdim), zb(ibdim),ibmax
      common /work/   rm(iwrdim) 
c 
      dimension iblank(jd,kd,ld),x(jd,kd,ld),y(jd,kd,ld),z(jd,kd,ld) 
c 
      data epslon /1.e-04/
c 
      ijmax = ibmax
      if( ijmax.gt.iwrdim ) then
        write( 6,601 ) ijmax,iwrdim 
601     format( ' ', 10x,'failure in hlocat: ijmax exceeds iwrdim',
     .  1x,'ijmax = ',i5,1x,'iwrdim = ',i5  ) 
        stop 'hlocat' 
      end if 
c 
c     for radius test, compute origin for boundary of mesh m1 and
c     max radius from origin to boundary
c 
      ii     = 0
      do 11 ib = 1,ibmax
      ii     = ii + 1 
      rm(ii) = xb(ib)
   11 continue
      i      = ismax( ijmax,rm,1 ) 
      xmax   = rm(i)
      i      = ismin( ijmax,rm,1 ) 
      xmin   = rm(i) 
      ii     = 0
      do 12 ib = 1,ibmax
      ii     = ii + 1 
      rm(ii) = yb(ib)
   12 continue
      i      = ismax( ijmax,rm,1 ) 
      ymax   = rm(i)
      i      = ismin( ijmax,rm,1 ) 
      ymin   = rm(i)
      ii     = 0
      do 13 ib = 1,ibmax
      ii     = ii + 1 
      rm(ii) = zb(ib)
   13 continue
      i      = ismax( ijmax,rm,1 ) 
      zmax   = rm(i)
      i      = ismin( ijmax,rm,1 ) 
      zmin   = rm(i)
      xorg   = ( xmax +xmin )/2.0 
      yorg   = ( ymax +ymin )/2.0 
      zorg   = ( zmax +zmin )/2.0 
      ii     = 0
      do 14 ib = 1,ibmax
      ii     = ii + 1 
      rm(ii) = (xb(ib)-xorg)**2 + (yb(ib)-yorg)**2 
     .       + (zb(ib)-zorg)**2
   14 continue
      i     = ismax( ijmax,rm,1 ) 
      rmax  = rm(i)
c 
c      write(6,*) ' ........in subroutine hlocat.................'
c      write(6,*)'  .............xmin,xmax ',xmin,xmax
c      write(6,*)'  .............ymin,ymax ',ymin,ymax
c      write(6,*)'  .............zmin,zmax ',zmin,zmax
c      write(6,*)'  .............rmax ',rmax
c 
c     locate any points within a hole
c
      do 34 j = 1,jd
      do 34 k = 1,kd
      do 34 l = 1,ld
c 
c     first check if point j,k,l of mesh m lies within the minimum and
c     maximum boundary points of mesh m1. If not then j,k,l does not lie
c     inside mesh m1
c
      if ((x(j,k,l) .lt. xmin).or.
     .(x(j,k,l) .gt. xmax) .or.
     .(y(j,k,l) .lt. ymin) .or.
     .(y(j,k,l) .gt. ymax) .or.
     .(z(j,k,l) .lt. zmin) .or.
     .(z(j,k,l) .gt. zmax)) go to 34
c
      dx = x(j,k,l)-xorg
      dy = y(j,k,l)-yorg
      dz = z(j,k,l)-zorg
      r =  dx*dx + dy*dy +dz*dz 
c
c     check if radius r is smaller than rmax 
c 
c 
      if( r.gt.rmax ) then
c       point not located in mesh m1; iblank does not change 
      else
c 
c       locate nearest point on hole boundary
c 
        ii = 0
        do 31 ib = 1,ibmax
        ii = ii + 1 
        dx = x(j,k,l) - xb(ib) 
        dy = y(j,k,l) - yb(ib) 
        dz = z(j,k,l) - zb(ib) 
        rm(ii) = dx**2 + dy**2 + dz**2
31      continue
        in = ismin( ijmax,rm,1 ) 
c
c       locate 2nd nearest point on hole boundary
c
        rm(in) = 1.e30
        in2 = ismin( ijmax,rm,1 ) 
c
c       locate 3rd nearest point on hole boundary
c
        rm(in2) = 1.e30
        in3 = ismin( ijmax,rm,1 ) 
c 
c       form dot products
c 
        rpx = x(j,k,l) - xb(in) 
        rpy = y(j,k,l) - yb(in) 
        rpz = z(j,k,l) - zb(in)         
        d   = rpx*vnx(in) + rpy*vny(in) + rpz*vnz(in) 
c
        rpx = x(j,k,l) - xb(in2) 
        rpy = y(j,k,l) - yb(in2) 
        rpz = z(j,k,l) - zb(in2) 
        d2  = rpx*vnx(in2) + rpy*vny(in2) + rpz*vnz(in2) 
c 
        rpx = x(j,k,l) - xb(in3) 
        rpy = y(j,k,l) - yb(in3) 
        rpz = z(j,k,l) - zb(in3) 
        d3  = rpx*vnx(in3) + rpy*vny(in3) + rpz*vnz(in3) 
c
c       interior point tests...use inner product test on nearest 3 
c       boundary points...if two of these pass the test, identify
c       point as interior...usually, just the nearest point is put 
c       through the inner product test - three points checked here
c       to try to reduce occurances of spurious hole point generation
c 
c        if(d .le. epslon) then
        if( (d .le. epslon .and. d2 .le. epslon) .or.
     .      (d .le. epslon .and. d3 .le. epslon)) then 
c         interior point to hole boundary; set iblank to zero
          iblank(j,k,l) = 0
c         write(8,*)' hole point=> j,k,l ',j,k,l,iblank(j,k,l)
c         write(8,*)' x,y,z', x(j,k,l),y(j,k,l),z(j,k,l)
        end if
c
      end if
34    continue
c 
      return
      end 
      subroutine hlocat2(m,jd,kd,ld,jmn,jmx,kmn,kmx,lmn,lmx,iblank)
c
c***********************************************************************
c     Purpose: locate the points interior to the hole boundary in
c     mesh m due to the grid m1.
c
c     simplified version of subroutine hlocat, wherein a hole in mesh m
c     is cut only by coordinate surfaces in mesh m (i.e. the same mesh).
c     thus, iblank values can be set simply by zeroing out the
c     appropriate index ranges, without having to resort to the more
c     complex algorithm in subroutine hlocat.
c***********************************************************************
c
      dimension iblank(jd,kd,ld)
c
      write(6,*)
      write(6,*)'    hole limits:'
      write(6,*)'    jmn,jmx ',jmn,jmx
      write(6,*)'    kmn,kmx ',kmn,kmx
      write(6,*)'    lmn,lmx ',lmn,lmx
c
      do 10 j=jmn,jmx
      do 10 k=kmn,kmx
      do 10 l=lmn,lmx
      iblank(j,k,l) = 0
10    continue
c
      return
      end
      subroutine hole
c 
c***********************************************************************
c     Purpose: construct hole in mesh m introduced by the presence of 
c     mesh m1 and compute the interpolation stencils for mesh m fringe 
c     points surrounding the hole. 
c 
c     iblank = 0 for points interior to hole
c     iblank < 0 for points surrounding hole (fringe points)
c     iblank = 1 for regular field points
c
c     if ihplt(nh) > 0  plot3d files for the initial hole and normals are
c     output. hole (grid file) is in hole_grd.nh, components of normal vector 
c     (q file) are in hole_nrm.nh. plot hole points with "points" option
c     in plot3d. plot normals using vector velocity. nh is the hole number.
c***********************************************************************
c 
c 
      include 'mag1.h'
c 
      common /book3 / iblank(mlen) 
      common /bound1/ ihbtyp(mdim,mdim),jh1(ipmax),jh2(ipmax),
     .                kh1(ipmax),kh2(ipmax),lh1(ipmax),lh2(ipmax),
     .                ip1(ipmax),ip2(ipmax),mh(ipmax),ihplt(ipmax) 
      common /grid1 / x(mlen), y(mlen), z(mlen)
      common /grid2 / mjmax(mdim), mkmax(mdim), mlmax(mdim)
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
      common /surf/   xb(ibdim), yb(ibdim), zb(ibdim), ibmax
      common /surf2/  xbo(idim),ybo(idim),zbo(idim)
      common /norm/   vnx(ibdim), vny(ibdim),vnz(ibdim) 
      common /intrp2/ jb(idim), kb(idim), lb(idim), jn(idim), kn(idim),
     .                ln(idim), itotal 
      common /intrp1/ xi(idim), yi(idim), zi(idim) 
      common /where/  nblkpt(idim)
      common /image/  jimage(mlen), kimage(mlen), limage(mlen) 
      common /tol/    epsc
      common /pltpt/  xbnd(idim),ybnd(idim),zbnd(idim),xorph(idim),
     .                yorph(idim),zorph(idim),xill(idim),yill(idim),
     .                zill(idim),nbnd,norph,nill
      common /trace1/ itrace
      common /snafu/  iholeh(mdim,mhldim),iholeo(mdim,mhldim),
     .                iorphh(mdim,mhldim),iorpho(mdim,mhldim)
c
      character*20 titlhlgrd, titlhlnrm
c 
      dimension jbtmp(idim),kbtmp(idim),lbtmp(idim),
     .          jntmp(idim),kntmp(idim),lntmp(idim),xitmp(idim),
     .          yitmp(idim),zitmp(idim)
      dimension list(mdim),iskip(idim)
c
      character*4 flag(3)
c  
c     itrace < 0, do not write search history for current fringe point
c     itrace = 0, overwrite history from previous point with current 
c     itrace = 1, retain the search history for ALL points (may get huge file)
c     trace output found in unit 7
c
      itrace = -1
c
      flag(1) = 'xie' 
      flag(2) = 'eta' 
      flag(3) = 'zeta' 
c

      nh = 0
      nhplt = 0
c  
c     construct hole boundaries and locate fringe points of mesh m
c  
      do 13 m = 1,nmesh 
      jd = mjmax(m) 
      kd = mkmax(m) 
      ld = mlmax(m) 
      ntot = nhole(m)
      itotal = 0
c  
      if( ntot.gt.0 ) then
c
        write(6,*)
        write(6,*)
        write(6,*) '  *** beginning hole construction for mesh ',
     .  m,' ***'
c
        call getibl( m,iblank,jd,kd,ld )
c 
c       loop over all holes in mesh m
c
        do 11 n = 1,ntot 
        nh  = nh+1
c
        write(6,*)
        write(6,*) '    creating hole ',n,' in mesh ',m,
     .  ' using ',ip2(nh)-ip1(nh)+1,' coordinate surfaces'
c
c       loop through coordinate surfaces to define initial 
c       hole boundary
c
        m1last = 0
        mseg   = 0
        ibmax  = 0      
c
c       test for holein mesh m defined only by coorinate surfaces 
c       of mesh m - then can simply blank out the appropriate j,k,l
c       range without the using hlocat (hlocat is less reliable since
c       ocassional, it generates spurious hole pts. due to corners)
c
        ihtest = 1
        do 117 ip = ip1(nh),ip2(nh)
        if (mh(ip) .ne. m) ihtest = 0
117     continue
c
        if (ihtest .ne. 0) then
c
           jmn = 10000000
           jmx = 0
           kmn = 10000000
           kmx = 0
           lmn = 10000000
           lmx = 0
           do 118 ip = ip1(nh),ip2(nh)
           m1 = mh(ip)
           js = jh1(ip)
           je = jh2(ip)
           ks = kh1(ip)
           ke = kh2(ip)
           ls = lh1(ip)
           le = lh2(ip)
           jmn = min(jmn,abs(js))
           jmn = min(jmn,abs(je))
           kmn = min(kmn,abs(ks))
           kmn = min(kmn,abs(ke))
           lmn = min(lmn,abs(ls))
           lmn = min(lmn,abs(le))
           jmx = max(jmx,abs(js))
           jmx = max(jmx,abs(je))
           kmx = max(kmx,abs(ks))
           kmx = max(kmx,abs(ke))
           lmx = max(lmx,abs(ls))
           lmx = max(lmx,abs(le))
c          cell center range (js,je, etc. based on grid points):
           if (jmx .gt. 1) jmx = jmx-1
           if (kmx .gt. 1) kmx = kmx-1
           if (lmx .gt. 1) lmx = lmx-1
118        continue 
c
           call hlocat2(m,jd,kd,ld,jmn,jmx,kmn,kmx,lmn,lmx,iblank)
c
        else
c
           do 14 ip = ip1(nh),ip2(nh)
           m1 = mh(ip)
           js = jh1(ip)
           je = jh2(ip)
           ks = kh1(ip)
           ke = kh2(ip)
           ls = lh1(ip)
           le = lh2(ip)
           if (m1.ne.m1last) then
              jd1 = mjmax(m1)+1
              kd1 = mkmax(m1)+1
              ld1 = mlmax(m1)+1
              call getgrd2( m1,x,y,z,jd1,kd1,ld1 )
              m1last = m1
           end if
c
           mseg = mseg + 1
           write(6,*) '      coordinate surface number ',mseg,
     .     '...from mesh ',m1
           write(6,*) '      ......jh1,jh2,kh1,kh2,lh1,lh2 = ',
     .     js,je,ks,ke,ls,le
c
c          determine normals to the current surface ip of mesh m1
c
           call bsurf(js,je,ks,ke,ls,le,jd1,kd1,ld1,x,y,z,m1)
c
14         continue
c
c          plot files for initial hole boundaries
c
           nhplt = nhplt+1
           if (ihplt(nhplt).gt.0) then
              if (nhplt.gt.99) then
                 len = 12
                 write (titlhlgrd,'("hole_grd.",i3)') nhplt
                 write (titlhlnrm,'("hole_nrm.",i3)') nhplt
              else if (nhplt.gt.9) then
                len = 11
                write (titlhlgrd,'("hole_grd.",i2)') nhplt
                write (titlhlnrm,'("hole_nrm.",i2)') nhplt
              else
                len = 10
                write (titlhlgrd,'("hole_grd.",i1)') nhplt
                write (titlhlnrm,'("hole_nrm.",i1)') nhplt
              endif
              do i = len+1, 20
                  titlhlgrd(i:i) = ' '
                  titlhlnrm(i:i) = ' '
              end do
              iunit = 35
c
              open(iunit,file=titlhlgrd(1:len),
     .        form='formatted',status='unknown')
              write(6,*)
              write(6,*)'    ',titlhlgrd(1:len),
     .        ' and ',titlhlnrm(1:len),' contain plot3d files'
              write(6,*)'    (ascii) for this initial hole boundary',
     .        '...plot using'
              write(6,*)'    POINTS for grid and VECTOR',
     .        ' VELOCITY for q (normal vectors)'
              write(iunit,*)ibmax,1,1
              write(iunit,*)(xb(ib),ib=1,ibmax),
     .                    (yb(ib),ib=1,ibmax),
     .                    (zb(ib),ib=1,ibmax)
c
              open(iunit,file=titlhlnrm(1:len),
     .        form='formatted',status='unknown')
              write(iunit,*)ibmax,1,1
c             (dummy) mach, alpha, re, time
              write(iunit,*)0.1,0.0,0.0,0.0
              write(iunit,*)(1.,ib=1,ibmax),
     .                    (vnx(ib),ib=1,ibmax),
     .                    (vny(ib),ib=1,ibmax),
     .                    (vnz(ib),ib=1,ibmax),
     .                    (1.,ib=1,ibmax)
           end if
c
c          determine actual hole in mesh m as being all points of mesh m 
c          that lie inside of initial hole boundary 
c
           call getgrd( m,x,y,z,jimage,kimage,limage,jd,kd,ld ) 
           call hlocat( iblank,x,y,z,jd,kd,ld) 
c
        end if
c 
c       identify fringe points in mesh m as points adjacent to hole n in mesh m
c       such points are identified with iblank = -(n+mdim)
c
        call frnge( iblank,jd,kd,ld,n )
c
11      continue
c
        call putibl( m,iblank,jd,kd,ld )
c
      end if
c  
13    continue
      do 1313 m = 1,nmesh 
      jd = mjmax(m) 
      kd = mkmax(m) 
      ld = mlmax(m) 
      ntot = nhole(m)
      itotal = 0
c  
      if( ntot.gt.0 ) then
c
        do 610 mm=1,nmesh
        list(mm) = 0
610     continue
c
        write(6,*)
        write(6,*)
        write(6,*) '  *** beginning determination of stencils for',
     .  ' fringe points of mesh ',m,' ***'
c 
c       loop over all holes in mesh m
c
        do 1111 n = 1,ntot 
        nserch = mhole(m,n) 
        iorph = iorphh(m,n)
c
        call getgrd( m,x,y,z,jimage,kimage,limage,jd,kd,ld ) 
        call getibl( m,iblank,jd,kd,ld )
c 
c       place fringe points of mesh m into  1-d lists jb,kb,lb of points to be
c       interpolated
c
        i1 = itotal + 1
        do 111 l = 1,ld
        do 111 k = 1,kd
        do 111 j = 1,jd
        i = j +(k-1)*jd +(l-1)*jd*kd
        if(iblank(i).eq.-(n+mdim)) then
          itotal      = itotal + 1
          if( itotal .gt.idim ) then
            write(6,602) itotal,idim 
602         format('0',130('*'),/,/,5x,'failure in hole: total' 
     .      ,1x,'number of fringe points exceeds available storage', 
     .      /,20x,'itotal = ',i6,2x,'idim = ',i6, 
     .      /,/,130('*') )
            stop 'hole'
          end if
          jb(itotal)  = j
          kb(itotal)  = k
          lb(itotal)  = l
          xbo(itotal) = x(i)
          ybo(itotal) = y(i)
          zbo(itotal) = z(i)
        end if
111     continue
c
        write(6,*)
        write(6,*) '    there are ',itotal-i1+1,' fringe points',
     .  ' due to hole ',n
c 
c       initialize block pointer array nblkpt
c       ...nblkpt(i) = 0 if interp. stencil for fringe point i not yet found
c       ...nblkpt(i) = -(m1+mdim)  if interp. stencil for fringe point i 
c                                  has been found in mesh m1
        do 112 i=i1,itotal
        nblkpt(i) = 0 
112     continue
c
c       search over all meshes in search list to find interpolation
c       stencils for fringe points of mesh m associated with hole n
c
        do 12 nn = 1,nserch
        m1  = lhole( m,n,nn ) 
        jd1 = mjmax(m1)
        kd1 = mkmax(m1)
        ld1 = mlmax(m1)
        call getgrd( m1,x,y,z,jimage,kimage,limage,jd1,kd1,ld1 ) 
        call intpt( -1,jimage,kimage,limage,jd,kd,ld, 
     .  m1,x,y,z,jd1,kd1,ld1,m,i1 ) 
c       reject stencils that contain hole/fringe pts.
          call getibl( m1,iblank,jd1,kd1,ld1 )
          do 1114 i=i1,itotal
          if(nblkpt(i) .eq. -(m1+mdim)) then
             j     = jn(i)
             k     = kn(i)
             l     = ln(i)
             jp1   = min( j+1,jd1 )
             kp1   = min( k+1,kd1 )
             lp1   = min( l+1,ld1 )
             ii1 = j   + (k-1)*jd1   + (l-1)*jd1*kd1
             ii2 = jp1 + (k-1)*jd1   + (l-1)*jd1*kd1
             ii3 = jp1 + (kp1-1)*jd1 + (l-1)*jd1*kd1
             ii4 = j   + (kp1-1)*jd1 + (l-1)*jd1*kd1
             ii5 = j   + (k-1)*jd1   + (lp1-1)*jd1*kd1
             ii6 = jp1 + (k-1)*jd1   + (lp1-1)*jd1*kd1
             ii7 = jp1 + (kp1-1)*jd1 + (lp1-1)*jd1*kd1
             ii8 = j   + (kp1-1)*jd1 + (lp1-1)*jd1*kd1
             if(iblank(ii1) .le. 0 .or.
     .          iblank(ii2) .le. 0 .or.
     .          iblank(ii3) .le. 0 .or.
     .          iblank(ii4) .le. 0 .or.
     .          iblank(ii5) .le. 0 .or.
     .          iblank(ii6) .le. 0 .or.
     .          iblank(ii7) .le. 0 .or.
     .          iblank(ii8) .le. 0) then
                  nblkpt(i) = 0
                  write(88,*)'in hole, rejecting pt. i = ',i,
     .            ' mesh ',m
             end if
          end if
 1114     continue
12      continue
c 
c       check to see if interpolation stencils were not found for any 
c       mesh m fringe points. for such "orphan" points, use either the nearest
c       point from one of the meshes in the search list, with xi=eta=zeta=0 
c       (zeroth order interpolation), or, use extrapolation from the nearest 
c       point.  extrapolation is distinguished from interpolation in that an
c       interpolation stencil has 0 < xie,eta,zeta < 1 (to within a tolerance
c       epsc) while for extrapolation, either xie, eta or zeta is < 0 or > 1
c       use of the nearest neighbor or extrapolation should occur only from 
c       boundaries of mesh m1; otherwise search routine has failed in some way
c       
        notok = 0
        do 113 i=i1,itotal
        if(nblkpt(i).eq.0) notok = notok + 1
113     continue
c
        if(notok.gt.0) then
          iflg = 0
          call orphan(m,n,nserch,i1,iorph,iflg)
        end if    
c
c       check to see if any interpolation stencils contain hole points 
c       if any do, check all other meshes in search list (if more than one)
c       to see if a valid stencil can be found elsewhere
c       
        if(iholeh(m,n).gt.0) then
          iflg = 0
          call reserch(m,n,nserch,i1,iskip,iorph,iflg)
        end if
1111    continue 
c
        do 612 n = 1,ntot 
        nserch = mhole(m,n)
        do 612 nn=1,nserch
        m1 = lhole(m,n,nn) 
        if(list(m1).eq.0) list(m1) = 1
612     continue
c 
c       summary of search routine results
c      
        write(6,*)
        write(6,*) '    summary of search routine results for mesh ',m,
     .  ' stencils:'
c
        nintot = 0
        nextot = 0
        nzerot = 0
c
        do 114 m1=1,nmesh
        if(list(m1).gt.0) then
          xtrap = 0.
          xtrap1 = 0.
          jd1    = mjmax(m1)
          kd1    = mkmax(m1)
          ld1    = mlmax(m1)
          nin    = 0
          nex    = 0
          norphb = 0
          norphc = 0 
          norphd = 0
          nzero  = 0
          do 115 i=1,itotal
          if(nblkpt(i).eq.-(m1+mdim)) then
            j    = jb(i)
            k    = kb(i)
            l    = lb(i)      
            xp   = xbo(i)
            yp   = ybo(i)
            zp   = zbo(i)
            jp   = jn(i)
            kp   = kn(i)
            lp   = ln(i)
            xie  = xi(i)
            eta  = yi(i)
            zeta = zi(i)
c
c           (orphan) points using zeroth order interpolation
c
            if(xie.eq.0. .and. eta.eq.0. .and. zeta.eq.0.)then
              nzero = nzero + 1
c             for stencils weighted totally at nearest neighbor,
c             keep all points in stencil within bounds of cell center grid,
c             grid,  1 .ge. j/k/l .le. jd1-1/kd-1/ld1-1  
c             for example, if nearest neighbor is at last cell center in
c             l-direction, say, l=ld1, then switch stencil from
c             jn,kn,ln with xie,eta,zeta=0 to jn,kn,ln-1 with xie=eta=0,
c             and zeta=1. The two stencils are equivalent
              if(jn(i).eq.jd1) then
                jn(i) = jn(i) -1
                xi(i) = 1.
                jp    = jn(i)
                xie   = xi(i)
              end if
              if(kn(i).eq.kd1) then
                kn(i) = kn(i) -1
                yi(i) = 1.
                kp    = kn(i)
                eta   = yi(i)
              end if
              if(ln(i).eq.ld1) then
c               check for 2d case
                if(ld1.gt.1) then
                  ln(i) = ln(i)-1
                end if
                zi(i) = 1.
                lp    = ln(i)
                zeta  = zi(i)
              end if
c
              if(jn(i).eq.1.or.jn(i).eq.jd1-1 .or.
     .        kn(i).eq.1.or.kn(i).eq.kd1-1 .or.
     .        ln(i).eq.1.or.ln(i).eq.ld1-1) norphb = norphb + 1
c
              write(8,*) 'for the mesh ',m,' fringe point',
     .        '  j,k,l = ',j,k,l
              write(8,*) '  will use nearest point in',
     .        ' mesh ',m1,': j,k,l = ',jp,kp,lp
              write(8,*) '    with xie,eta,zeta = ',xie,eta,zeta
c
c           points using interpolation
c
            else if(xie .ge.-epsc .and. xie .le.1.+epsc .and. 
     .      eta .ge.-epsc .and. eta .le.1.+epsc .and. 
     .      zeta.ge.-epsc .and. zeta.le.1.+epsc) then 
              nin = nin + 1
            end if
c
c           (orphan) points using extrapolation
c
            if(xie .lt.-epsc .or. xie .gt.1.+epsc .or. 
     .      eta .lt.-epsc .or. eta .gt.1.+epsc .or. 
     .      zeta.lt.-epsc .or. zeta.gt.1.+epsc) then 
              iflagg = 0
              norph        = norph + 1
              xorph(norph) = xp
              yorph(norph) = yp
              zorph(norph) = zp
              nex = nex + 1
              if(jn(i).eq.1.or.jn(i).eq.jd1-1 .or.
     .        kn(i).eq.1.or.kn(i).eq.kd1-1 .or.
     .        ln(i).eq.1.or.ln(i).eq.ld1-1) norphb = norphb + 1
c
              write(8,*) 'for the mesh ',m,' fringe point',
     .        '  j,k,l = ',j,k,l
              write(8,*) '  will use extrapolation from',
     .        ' mesh ',m1,' point j,k,l = ',jp,kp,lp
              write(8,*) '    with xie,eta,zeta = ',xie,eta,zeta
c
              if(xie.lt.-epsc .or. xie.gt.1.+epsc) then
                iflagg = iflagg+1
                if(xie .gt. xtrap) then
                  xtrap  = xie
                  ixtrap = i
                  mxtrap = m1
                  iflag  = 1
                end if
                if(xie .lt. xtrap1) then
                  xtrap1  = xie
                  ixtrap1 = i
                  mxtrap1 = m1
                  iflag1  = 1
                end if
              end if
              if(eta.lt.-epsc .or. eta.gt.1.+epsc) then
                iflagg = iflagg+1
                if(eta .gt. xtrap) then
                  xtrap  = eta
                  ixtrap = i
                  mxtrap = m1
                  iflag  = 2
                end if
                if(eta .lt. xtrap1) then
                  xtrap1  = eta
                  ixtrap1 = i
                  mxtrap1 = m1
                  iflag1  = 2
                end if
              end if
              if(zeta.lt.-epsc .or. zeta.gt.1.+epsc) then
                iflagg = iflagg+1
                if(zeta .gt. xtrap) then
                  xtrap  = zeta
                  ixtrap = i
                  mxtrap = m1
                  iflag = 3
                end if
                if(zeta .lt. xtrap1) then
                  xtrap1  = zeta
                  ixtrap1 = i
                  mxtrap1 = m1
                  iflag1  = 3
                end if
              end if
              if(xie.lt.-0.5-epsc .or. xie.gt.1.5+epsc) 
     .        norphc = norphc+1
              if(eta.lt.-0.5-epsc .or. eta.gt.1.5+epsc) 
     .        norphc = norphc+1
              if(zeta.lt.-0.5-epsc .or. zeta.gt.1.5+epsc) 
     .        norphc = norphc+1
              if(iflagg .gt. 1) norphd = norphd+1
            end if
          end if
115       continue
c
          write(6,*)
          write(6,*) '      ',nin,' fringe points of mesh ',m,
     .    ' are interpolated from mesh ',m1
          if(nex.gt.0) then  
          write(6,*) '      ',nex,' orphaned fringe points of mesh ',
     .    m,' are extrapolated from mesh ',m1
          write(6,*) '      ',norphc,' of these orphans have',
     .    ' extrapolation coefficients <-0.5 or >1.5'
          if(norphd.gt.0) then
            write(6,*)'      ',norphd,' of these orphans',
     .      ' are extrapolated in more than one direction'
          end if
          if(xtrap.gt.1.+epsc) then
          write(6,*)'        the maximum extrapolation coefficient is ',
     .    flag(iflag),' = ',xtrap
          write(6,*)'        to the fringe pt ',jb(ixtrap),kb(ixtrap),
     .    lb(ixtrap),' from the mesh ',mxtrap,' target pt ',
     .    jn(ixtrap),kn(ixtrap),ln(ixtrap)
          end if
          if(xtrap1.lt.-epsc) then
          write(6,*)'        the minimum extrapolation coefficient is ',
     .    flag(iflag1),' = ',xtrap1
          write(6,*)'        to the fringe pt ',jb(ixtrap1),kb(ixtrap1),
     .    lb(ixtrap1),' from the mesh ',mxtrap1,' target pt ',
     .    jn(ixtrap1),kn(ixtrap1),ln(ixtrap1)
          end if
          write(6,*)'        check file 8 for more details' 
          end if
          if(nzero.gt.0) then
          write(6,*) '      ',nzero,' orphaned fringe points of mesh ',
     .    m,' use the nearest point in mesh ',m1
          write(6,*)'        check file 8 for more details' 
          end if
          if(nex+nzero-norphb.ne.0)then
              write(6,*)'        WARNING: ',nex-norphb,' of these',
     .        ' orphans are extrapolated from an'
              write(6,*)'          interior point of mesh ',m1,
     .        '...possible miscue in search routine'
          end if
c
          nintot = nintot + nin
          nextot = nextot + nex
          nzerot = nzerot + nzero
c
        end if
c
114     continue
c
        if(nintot+nextot+nzerot .ne. itotal) then
          write(6,*)'  stopping...unable to compute ',
     .    itotal-nintot-nextot-nzerot,' stencils for mesh ',m,
     .    ' fringe points'
          write(6,*)'     possible miscue in search routine'
          stop
        end if
c
c       put data back in original order (i.e. all points 
c       interpolated/extrapolated from mesh m1 are sequential in list jb,kb,lb
c
        itotl1 = itotal
        do 120 i=1,itotl1
        jbtmp(i) = jb(i)
        kbtmp(i) = kb(i)
        lbtmp(i) = lb(i)
        jntmp(i) = jn(i)
        kntmp(i) = kn(i)
        lntmp(i) = ln(i)
        xitmp(i) = xi(i)
        yitmp(i) = yi(i)
        zitmp(i) = zi(i)
120     continue
c
        do 125 m1 = 1,nmesh 
        if(list(m1).gt.0) then
          ii = 0
          do 126 i=1,itotl1
          if(nblkpt(i) .eq.-(m1+mdim)) then
            ii = ii + 1
            jb(ii) = jbtmp(i)
            kb(ii) = kbtmp(i)
            lb(ii) = lbtmp(i)
            jn(ii) = jntmp(i)
            kn(ii) = kntmp(i)
            ln(ii) = lntmp(i)
            xi(ii) = xitmp(i)
            yi(ii) = yitmp(i)
            zi(ii) = zitmp(i)
          end if
126       continue
c   
c         set interpolation pointers and load lists for input to flow solver
c
          itotal = ii
          icase  = 1
          call setptr( m,m1,icase )
        end if
125     continue
c  
      end if
c
1313  continue
c
c     reset itrace to "off"
      itrace = -1
c
      return
      end 
c
      subroutine ibchk( ib,iblank,ic,is,ji,ki,li,m1,jd,kd,ld, 
     . xb,yb,zb,jinew,kinew,linew,x,y,z )
c 
c***********************************************************************
c     Purpose: set up array to check for fringe points, hole points or 
c     interpolated boundary points in the interpolation stencil made up 
c     of mesh m1 points. 
c
c     if all points in the interpolation stencil are field points (valid 
c     stencil), ib(i) is returned as zero.
c     
c     if the stencil contains hole and/or fringe points (but no interpolated)
c      
c

c     if all points in the stencil are invalid (i.e. each stencil point
c     is either a hole, fringe or interpolated boundary point), ib(i) 
c     is returned as zero
c
c     if at least one point in the stencil is a field point, then 
c     ib(i) is returned as negative; also returned are the indicies 
c     jinew,kinew,linew of the nearest field point in the stencil to 
c     the point to be interpolated (xb,yb,zb) 
c
c  
c***********************************************************************
c 
      dimension ib(*), iblank(jd,kd,ld),ji(*),ki(*),li(*) 
      dimension x(jd,kd,ld),y(jd,kd,ld),z(jd,kd,ld),jinew(*),kinew(*),
     .linew(*),xb(*),yb(*),zb(*) 
c 
c***********************************************************************
c 
      dmin0 = 1.e30
c
      ii = is - 1 
      do 11 i = 1,ic
      dmin  = dmin0
      ib(i) = 0
      ii    = ii + 1 
      j     = ji(ii) 
      k     = ki(ii) 
      l     = li(ii) 
      jp1   = min( j+1,jd )
      kp1   = min( k+1,kd )
      lp1   = min( l+1,ld )
      i1   = iblank(j,k,l)
      i2   = iblank(jp1,k,l)
      i3   = iblank(jp1,kp1,l)
      i4   = iblank(j,kp1,l)
      i5   = iblank(j,k,lp1)
      i6   = iblank(jp1,k,lp1)
      i7   = iblank(jp1,kp1,lp1)
      i8   = iblank(j,kp1,lp1)  
      ib1  = min( i1,i2,i3,i4,i5,i6,i7,i8)
      ib2  = max( i1,i2,i3,i4,i5,i6,i7,i8)
c
      if(ib1*ib2 .ne. 1) then
c
c       stencil contains at least one hole, fringe, or interpolated 
c       boundary point
c
        if(i1.eq.1) then
          d1 = (x(j,k,l)-xb(ii))**2+(y(j,k,l)-yb(ii))**2
     .    +(z(j,k,l)-zb(ii))**2
          if(d1.lt.dmin) then
            dmin = d1
            jmin = j
            kmin = k 
            lmin = l
          end if
        end if
        if(i2.eq.1) then
          d1 = (x(jp1,k,l)-xb(ii))**2+(y(jp1,k,l)-yb(ii))**2
     .    +(z(jp1,k,l)-zb(ii))**2
          if(d1.lt.dmin) then
            dmin = d1
            jmin = jp1
            kmin = k 
            lmin = l
          end if
        end if
        if(i3.eq.1) then
          d1 = (x(jp1,kp1,l)-xb(ii))**2+(y(jp1,kp1,l)-yb(ii))**2
     .    +(z(jp1,kp1,l)-zb(ii))**2
          if(d1.lt.dmin) then
            dmin = d1
            jmin = jp1
            kmin = kp1
            lmin = l
          end if
        end if
        if(i4.eq.1) then
          d1 = (x(j,kp1,l)-xb(ii))**2+(y(j,kp1,l)-yb(ii))**2
     .    +(z(j,kp1,l)-zb(ii))**2
          if(d1.lt.dmin) then
            dmin = d1
            jmin = j
            kmin = kp1 
            lmin = l
          end if
        end if
        if(i5.eq.1) then
          d1 = (x(j,k,lp1)-xb(ii))**2+(y(j,k,lp1)-yb(ii))**2
     .    +(z(j,k,lp1)-zb(ii))**2
          if(d1.lt.dmin) then
            dmin = d1
            jmin = j
            kmin = k 
            lmin = lp1
          end if
        end if
        if(i6.eq.1) then
          d1 = (x(jp1,k,lp1)-xb(ii))**2+(y(jp1,k,lp1)-yb(ii))**2
     .    +(z(jp1,k,lp1)-zb(ii))**2
          if(d1.lt.dmin) then
            dmin = d1
            jmin = jp1
            kmin = k 
            lmin = lp1
          end if
        end if
        if(i7.eq.1) then
          d1 = (x(jp1,kp1,lp1)-xb(ii))**2
     .    +(y(jp1,kp1,lp1)-yb(ii))**2+(z(jp1,kp1,lp1)-zb(ii))**2
          if(d1.lt.dmin) then
            dmin = d1
            jmin = jp1
            kmin = kp1
            lmin = lp1
          end if
        end if
        if(i8.eq.1) then
          d1 = (x(j,kp1,lp1)-xb(ii))**2+(y(j,kp1,lp1)-yb(ii))**2
     .    +(z(j,kp1,lp1)-zb(ii))**2
          if(d1.lt.dmin) then
            dmin = d1
            jmin = j
            kmin = kp1 
            lmin = lp1
          end if
        end if
c
        if(ib1.le.0) then
c
c         stencil contains hole/fringe points
c
          if(dmin .lt. dmin0) then
            ib(i) = -2
            jinew(ii) = jmin
            kinew(ii) = kmin
            linew(ii) = lmin
          else
            ib(i) = -1
          end if
        end if
c
        if(ib2.gt.1) then
c
c         stencil contains interpolated boundary points
c
          if(dmin .lt. dmin0) then
            ib(i) = 2
            jinew(ii) = jmin
            kinew(ii) = kmin
            linew(ii) = lmin
          else
            ib(i) = 1
          end if
        end if
c
      end if
c
  11  continue
c 
      return
      end 
      subroutine iblcon( iblank,jd,kd,ld )
c
c***********************************************************************
c     Purpose: remove the connection data from iblank by setting 
c     iblank = 0 if iblank .ne. 1.
c***********************************************************************
c 
      dimension iblank(jd,kd,ld)
c 
      do 11 l = 1,ld
      do 11 k = 1,kd
      do 11 j = 1,jd
       if( iblank(j,k,l).ne.1 )  iblank(j,k,l)=0
   11 continue
c 
      return
      end 
      subroutine iblplt(m,iblank,jd,kd,ld)
c********************************************************************** 
c     Purpose: print the iblank array for each mesh. 
c********************************************************************** 
c 
      dimension iblank(jd,kd,ld)
c 
      write(6,299) m 
c 
      do 10 l=1,ld
        icell =0
        do 30 k=1,kd
        do 30 j=1,jd
        if(iblank(j,k,l).ne.1) icell = icell+1
 30   continue
        if(icell.eq.0) then
              write(6,302) l
        else
              write(6,301) l 
        do 20 k=1,kd
              write(6,305)(abs(iblank(j,k,l)),j=1,jd) 
 20     continue
        end if
 10   continue
c 
 299   format(/,/,' iblank array for cell centers of mesh=',i3) 
 301   format(' .....l=constant= ',i3)
 302   format(' .....all values unity for l=constant= ',i3)
 305   format(' ',130i1)
c 
      return 
      end
      subroutine initia 
c
c***********************************************************************
c     Purpose: read connection and grid parameters
c***********************************************************************
c 
      include 'mag1.h'
c 
      common /book2/  dxint(idim), dyint(idim), dzint(idim),
     .                ibpts(mdim), jbpt(idim), kbpt(idim), lbpt(idim),
     .                ji(idim),ki(idim), li(idim)
      common /book3/  iblank(mlen) 
      common /grdlim/ jsrs(mdim),jsre(mdim),ksrs(mdim),ksre(mdim),
     .                lsrs(mdim),lsre(mdim) 
      common /grid2/  mjmax(mdim),mkmax(mdim),mlmax(mdim)
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
      common /bound1/ ihbtyp(mdim,mdim),jh1(ipmax),jh2(ipmax),
     .                kh1(ipmax),kh2(ipmax),lh1(ipmax),lh2(ipmax),
     .                ip1(ipmax),ip2(ipmax),mh(ipmax),ihplt(ipmax) 
      common /bound2/ iobtyp(mdim),nobtyp(mdim),jo1(ipmax),jo2(ipmax),
     .                ko1(ipmax),ko2(ipmax),lo1(ipmax),lo2(ipmax)
      common /iseq/   iseq(mdim),icount(mdim)
      common /diagno/ jp(mdim*3+3),kp(mdim*3+3),lp(mdim*3+3),iplt3d
      common /snafu/  iholeh(mdim,mhldim),iholeo(mdim,mhldim),
     .                iorphh(mdim,mhldim),iorpho(mdim,mhldim)
      common /igrdtyp/ ip3dgrd,ialph
c 
      dimension title(20)
c
c ***************************************************************
c
      read(5,10)  (title(i),i=1,20)
      write(6,11) (title(i),i=1,20)
   10 format(20a4)
   11 format(/,1x,20a4)
c
      read(5,10)
      write(6,21)
      read(5,*)  nmesh,iplt3d
      write(6,20) nmesh,iplt3d
   20 format(8i10)
   21 format(5x,5hnmesh,4x,6hiplt3d)
c
      ip3dgrd = 0
      if (nmesh.lt.0) then
         ip3dgrd = 1
         nmesh =  iabs(nmesh)
      end if
c
      if(nmesh.gt.mdim) then
        write(6,5) mdim,nmesh
5       format(5x,5hmdim=,i2,2x,6hnmesh=,i2,16hincrease mdim...,
     .  25hmust have mdim .ge. nmesh)
        stop
      end if
c
      read(5,10)
      write(6,101)
      nout    = 0 
      nob     = 0
      ihotest = 0
      do 100 n=1,nmesh
      read(5,91) mjmax(n),mkmax(n),mlmax(n),noutr(n),nhole(n),nobtyp(n)
      write(6,91)mjmax(n),mkmax(n),mlmax(n),noutr(n),nhole(n),nobtyp(n)
      mjmax(n) = mjmax(n)-1
      mkmax(n) = mkmax(n)-1
      mlmax(n) = mlmax(n)-1
      nout     = nout + min(noutr(n),1)
      nob      = nob  + nobtyp(n)
      if (nhole(n).gt.0) ihotest = ihotest + nhole(n)
      if(nhole(n).gt.mhldim .or. noutr(n).gt.mhldim) then
        write(6,92)
92      format(2x,'stopping...number of holes/outer boundaries',
     .  ' exceeds parameter mhldim')
        stop
      end if
  100 continue
   91 format(7i10,l10)
  101 format(6x,4hjmax,6x,4hkmax,6x,4hlmax,5x,5hnoutr,5x,5hnhole,4x,
     .6hnobtyp)
c
      if(nob.gt.ipmax) then
        write(6,4051)
 4051   format(2x,'stopping...number of planes defining outer',
     .  ' boundary points to be interpolated exceeds ipmax')
        stop
      end if
c
      read(5,10)
      write(6,151)
      do 150 n=1,nmesh
      read(5,*)  jsrs(n),jsre(n),ksrs(n),ksre(n),lsrs(n),lsre(n)
c
      if (jsrs(n).le.0)                          jsrs(n) = 1
      if (jsre(n).le.0 .or. jsre(n).gt.mjmax(n)) jsre(n) = mjmax(n)
      if (ksrs(n).le.0)                          ksrs(n) = 1
      if (ksre(n).le.0 .or. ksre(n).gt.mkmax(n)) ksre(n) = mkmax(n)
      if (lsrs(n).le.0)                          lsrs(n) = 1
      if (lsre(n).le.0 .or. lsre(n).gt.mlmax(n)) lsre(n) = mlmax(n)
c
      write(6,20) jsrs(n),jsre(n),ksrs(n),ksre(n),lsrs(n),lsre(n)
  150 continue
  151 format(6x,4hjsrs,6x,4hjsre,6x,4hksrs,6x,4hksre,6x,
     .       4hlsrs,6x,4hlsre)
c
      do 200 n=1,nmesh
      do 200 m=1,nmesh
      ihbtyp(m,n) = 0
  200 continue
c
      if (ihotest.gt.0) then
      nt    = 0
      ntest = 0
      jp2   = 0
      ihtot = 0
      do 600 i=1,ihotest
c
      do 350 n=1,nmesh
      icount(n) = 0
  350 continue
c
      read(5,10)
      write(6,399)
c
      read(5,10)
      write(6,401)
      read(5,*)  nn,nhle,(ihbtyp(nn,m),m=1,nmesh)
c
      if (nn.lt.ntest) then
         write(6,340)nn,ntest
         stop
      else
         ntest = nn
         ihplt(i) = 0
         if(nhle.lt.0) then
           ihplt(i) = 1
           nhle = iabs(nhle)
         end if
      end if
  399 format(1x,16hhole definition:)
  340 format(/,/1x,46hstopping ... incorrect order for IHBTYP lines;,
     .  /1x,57hlines must be in ascending order according to mesh number
     .  /1x,32hcurrent mesh number being read =,i3,
     .   1x,24h; previous mesh number =,i3)
      write(6,20) nn,nhle,(ihbtyp(nn,m),m=1,nmesh)
  401 format(6x,4hmesh,2x,8hhole no.,1x,9hihbtyp(1),1x,9hihbtyp(2),1x,
     .       9hihbtyp(3),1x,9hihbtyp(4),1x,9hihbtyp(5))
c
      ih = 0
      do 360 m=1,nmesh
      ih = ih+ihbtyp(nn,m)
  360 continue
c
      jp1    = jp2 + 1
      jp2    = jp2 + ih
      ip1(i) = jp1
      ip2(i) = jp2
c
c
        ihtot = ihtot + ih
	if(ihtot.gt.ipmax) then
	  write(6,6051)
 6051   format(2x,'stopping...number of planes defining holes',
     .  ' exceeds ipmax')
	  stop
	end if
c
      iseq(1) = nt
      do 400 m=2,nmesh
      iseq(m) = ihbtyp(nn,m-1)+iseq(m-1)
  400 continue
c
      read(5,10)
      write(6,605)
c
      do 500 n=1,ih
      read(5,*)  mm,j1,j2,k1,k2,l1,l2
c
      if (icount(mm).gt.ihbtyp(nn,mm)) then
         write(6,501)mm,nn
         stop
      end if
  501 format(/,/1x,43hstopping ... too many holes defined by mesh,i3,
     .       1x,10hfor mesh =,i3)
c
      if (ihbtyp(nn,mm).eq.0) then
         write(6,222)nn
         write(6,502)mm
         stop
      else
         icount(mm) = icount(mm)+1
      end if
  222 format(/,/1x,44hstopping ... inconsistent input for mesh m =,i3)
  502 format(/1x,50hattempted entry of JH1, JH2... values for mesh n =,
     .       i3,1x,19hbut IHBTYP(m,n) = 0)
c
      if (j1.eq.0) j1 = 1
      if (j2.eq.0) j2 = mjmax(mm)+1
      if (k1.eq.0) k1 = 1
      if (k2.eq.0) k2 = mkmax(mm)+1
      if (l1.eq.0) l1 = 1
      if (l2.eq.0) l2 = mlmax(mm)+1
c
c
c     check for inappropriate indices on components of hole surface 
c
      if (abs(j1) .gt. mjmax(mm)+1) then
         write(6,440) j1,mjmax(mm)+1
440      format(1x,17h  stopping: j1 = ,i4,8h jmax = ,i4)
         stop
      end if
      if (abs(j2) .gt. mjmax(mm)+1) then
         write(6,441) j2,mjmax(mm)+1
441      format(1x,17h  stopping: j2 = ,i4,8h jmax = ,i4)
         stop
      end if
      if (abs(k1) .gt. mkmax(mm)+1) then
         write(6,442) k1,mkmax(mm)+1
442      format(1x,17h  stopping: k1 = ,i4,8h kmax = ,i4)
         stop
      end if
      if (abs(k2) .gt. mkmax(mm)+1) then
         write(6,443) k2,mkmax(mm)+1
443      format(1x,17h  stopping: k2 = ,i4,8h kmax = ,i4)
         stop
      end if

      if (abs(l1) .gt. mlmax(mm)+1) then
         write(6,444) l1,mlmax(mm)+1
444      format(1x,17h  stopping: l1 = ,i4,8h lmax = ,i4)
         stop
      end if
      if (abs(l2) .gt. mlmax(mm)+1) then
         write(6,445) l2,mlmax(mm)+1
445      format(1x,17h  stopping: l2 = ,i4,8h lmax = ,i4)
         stop
      end if
c
      write(6,20) mm,j1,j2,k1,k2,l1,l2
      iseq(mm)      = iseq(mm)+1
c
      mh(iseq(mm))  = mm
      jh1(iseq(mm)) = j1
      jh2(iseq(mm)) = j2
      kh1(iseq(mm)) = k1
      kh2(iseq(mm)) = k2
      lh1(iseq(mm)) = l1
      lh2(iseq(mm)) = l2
c
  500 continue
c
      read(5,10)
      read(5,10)
      write(6,129)
  129 format(6x,41hsearch list for associated fringe points:)
      write(6,131)
  131 format(5x,5hnlist,4x,6hlst(1),4x,6hlst(2),
     .4x,6hlst(3),4x,6hlst(4),4x,10h...  nlist)
      read(5,*)mhole(nn,nhle),(lhole(nn,nhle,ll),ll=1,mhole(nn,nhle))
      write(6,20)mhole(nn,nhle),(lhole(nn,nhle,ll),ll=1,mhole(nn,nhle))
      do 1502 ll=1,mhole(nn,nhle)
      if(lhole(nn,nhle,ll).eq.nn) then
        write(6,123)nn,nn
  123 format(4x,16hstopping...mesh ,i3,27h not allowed in search list,
     .10h for mesh ,i3,14h fringe points)
        stop
      end if
 1502 continue
c
      read(5,10)
      read(5,10)
      write(6,329)
  329 format(6x,22hproblem point options:)
      write(6,331)
  331 format(5x,5hiorph,5x,5hihole)
      read(5,*)  iorphh(nn,nhle),iholeh(nn,nhle)
      write(6,20) iorphh(nn,nhle),iholeh(nn,nhle)
c
      nt = nt+ih
      do 580 n=1,nmesh
      if (icount(n).ne.ihbtyp(nn,n)) then
         write(6,222)nn
         write(6,581)n,icount(n),ihbtyp(nn,n)
         stop
      end if
  580 continue
  581 format(/1x,52hnumber of entries of JH1, JH2... values for mesh n =
     .       ,i3,1x,3his ,i3,1x,17hbut IHBTYP(m,n) =,i3)
  600 continue
c
      else
c
      read(5,10)
      write(6,399)
      read(5,10)
      write(6,401)
      read(5,10)
      write(6,605)
      read(5,10)
      write(6,129)
      read(5,10)
      write(6,131)
      read(5,10)
      write(6,329)
      read(5,10)
      write(6,331)
c
      end if
  605 format(5x,5himesh,7x,3hjh1,7x,3hjh2,7x,3hkh1,7x,3hkh2,
     .       7x,3hlh1,7x,3hlh2)
c
      if (nob.gt.0) then
         iseq(1) = 0
         do 650 m=2,nmesh
         iseq(m) = nobtyp(m-1)+iseq(m-1)
  650    continue
         do 700 nn=1,nmesh
         if(nobtyp(nn).gt.0)then
           read(5,10)
           read(5,10)
           write(6,699)
           write(6,701)
           do 7001 n=1,nobtyp(nn)
           read(5,*)  mm,j1,j2,k1,k2,l1,l2
           if (nobtyp(mm).eq.0) then
              write(6,222)mm
              write(6,702)
              stop
           end if
c
           if (j1.eq.0) j1 = 1
           if (j2.eq.0) j2 = mjmax(mm)
           if (k1.eq.0) k1 = 1
           if (k2.eq.0) k2 = mkmax(mm)
           if (l1.eq.0) l1 = 1
           if (l2.eq.0) l2 = mlmax(mm)
c 
           if (j1.lt.0) j2 = 0
           if (k1.lt.0) k2 = 0
           if (l1.lt.0) l2 = 0
c  
           write(6,20) mm,j1,j2,k1,k2,l1,l2
           iseq(mm)      = iseq(mm)+1
c
           jo1(iseq(mm)) = j1
           jo2(iseq(mm)) = j2
           ko1(iseq(mm)) = k1
           ko2(iseq(mm)) = k2
           lo1(iseq(mm)) = l1
           lo2(iseq(mm)) = l2
c
 7001      continue
c
           read(5,10)
           read(5,10)
           write(6,229)
  229      format(4x,38hsearch list for outer boundary points:)
           write(6,231)
  231      format(5x,5hnlist,4x,6hlst(1),4x,6hlst(2),
     .     4x,6hlst(3),4x,6hlst(4),4x,10h...  nlist)
c
           n=1
c
           read(5,*)moutr(nn,n),(loutr(nn,n,ll),ll=1,moutr(nn,n))
           write(6,20)moutr(nn,n),(loutr(nn,n,ll),ll=1,moutr(nn,n))
           do 2502 ll=1,moutr(nn,n)
           if(loutr(nn,n,ll).eq.nn) then
             write(6,223)nn,nn
  223        format(4x,16hstopping...mesh ,i3,
     .       27h not allowed in search list,
     .       10h for mesh ,i3,16h boundary points)
             stop
           end if
 2502      continue
c
           read(5,10)
           read(5,10)
           write(6,329)
           write(6,331)
           read(5,*)  iorpho(nn,n),iholeo(nn,n)
           write(6,20) iorpho(nn,n),iholeo(nn,n)
c
         end if
c
  700    continue
      else
c         read(5,10)
c         read(5,10)
c         write(6,699)
c         write(6,701)
      end if
  699 format(1x,49hdefinition of outer boundary points (cell center),
     .19hto be interpolated:)
  701 format(6x,4hmesh,7x,3hjo1,7x,3hjo2,7x,3hko1,7x,3hko2,
     .       7x,3hlo1,7x,3hlo2)
  702 format(/1x,37hattempted entry of JO1, JO2... values,
     .       1x,17Hbut NOBTYP(m) = 0)
c
      mb = 0
      do 750 m=1,nmesh 
      iobtyp(m) = 0
      if (nobtyp(m).ne.0) then 
         iobtyp(m) = mb+1
         do 725 nb=1,nobtyp(m)
         mb = mb+1
         if (jo1(mb).gt.mjmax(m)+2 .or. jo2(mb).gt.mjmax(m)+2 .or. 
     .       ko1(mb).gt.mkmax(m)+2 .or. ko2(mb).gt.mkmax(m)+2 .or. 
     .       lo1(mb).gt.mlmax(m)+2 .or. lo2(mb).gt.mlmax(m)+2 .or. 
     .       jo1(mb).lt.-1         .or. jo2(mb).lt.-1         .or. 
     .       ko1(mb).lt.-1         .or. ko2(mb).lt.-1         .or. 
     .       lo1(mb).lt.-1         .or. lo2(mb).lt.-1) then
             write(6,755) m
             stop
         end if
  725    continue
      end if
  750 continue
  755 format(/,/5x,'value of outer boundary index exceeds',
     .       1x,'constrained values on mesh ',i3)
c
c
c***********************************************************************
c                initialize iblank to 1 for each grid 
c***********************************************************************
c 
      do 1000 m=1,nmesh 
      jd = mjmax(m) 
      kd = mkmax(m) 
      ld = mlmax(m) 
      call putibl( m,iblank,jd,kd,ld )
 1000 continue
c 
c***********************************************************************
c                initialize the interpolation arrays for each grid
c***********************************************************************
c 
      nwr = 0
      do 1100 m=1,nmesh 
      call putint( nwr,m,ji,ki,li,jbpt,kbpt,lbpt,dxint,dyint,dzint )
 1100 continue
c 
c***********************************************************************
c                         consistency check
c***********************************************************************
c 
c     check total number of points less than mlen 
c 
      istop = 0
      do 1200 m=1,nmesh 
      ntotal = (mjmax(m)+1)*(mkmax(m)+1)*(mlmax(m)+1)
      if (ntotal.gt.mlen) then 
         write(6,1205) m,ntotal,mlen
         istop = 1 
      end if
 1200 continue
 1205 format('0',10x,'fatal error in mesh = ',i3,/, 
     .           10x,'the total number of points in', 
     .            1x,'mesh is ',i6,1x,'max dimension is ',i6  )
c 
c     check each dimension
c 
      if (nmesh.gt.mdim) then
         write(6,1305) nmesh,mdim 
         istop = 1
      else
         do 1300 m=1,nmesh 
         if (mjmax(m)+1.gt.jdim .or. mkmax(m)+1.gt.kdim .or.
     .       mlmax(m)+1.gt.ldim)  then
             write(6,1310) m,mjmax(m)+1,jdim,mkmax(m)+1,kdim, 
     .                    mlmax(m)+1,ldim
             istop = 1
         end if 
 1300    continue
      end if
 1305 format(1h0,130('*'),/,/,5x,'failure in initia,number of' 
     .       ,1x,'grids input,nmesh=',i3,3x,'is greater than dimension,'
     .       ,1x,'mdim =',i3,/,/,130('*') ) 
 1310 format('0',130('*'),/,/,5x,'input failure in initia',
     .        1x,'exceeds maximum dimension', /, 
     .       25x,'mesh =',i5,3x,'jmax,jdim',2i5,5x,'kmax,kdim',2i5, 
     .        5x,'lmax,ldim',2i5,/,/,130('*') ) 

c 
      if (istop.eq.1) stop 'initia'
      return
      end 
      subroutine interp(jmax1,kmax1,lmax1,x11,y11,z11,ji,ki,li,
     .                  xie,eta,zeta,xintrp,yintrp,zintrp)
c
c***********************************************************************
c    Purpose:  determine cell center coordinates via interpolation stencils
c***********************************************************************
c
      dimension x11(jmax1,kmax1,lmax1),y11(jmax1,kmax1,lmax1),
     .          z11(jmax1,kmax1,lmax1)
c
c     index increments for 2d or 3d case
c
      jinc1 = 1
      kinc1 = 1
      linc1 = 1
      if(jmax1.lt.2) jinc1 = 0
      if(kmax1.lt.2) kinc1 = 0
      if(lmax1.lt.2) linc1 = 0
c
c     cell center value of point in mesh m via trilinear 
c     interpolation in target cell ji,ki,li of mesh m1
c
      x1 = x11(ji,ki,li)
      y1 = y11(ji,ki,li)
      z1 = z11(ji,ki,li)
      x2 = x11(ji+jinc1,ki,li)
      y2 = y11(ji+jinc1,ki,li)
      z2 = z11(ji+jinc1,ki,li)
      x4 = x11(ji,ki+kinc1,li)
      y4 = y11(ji,ki+kinc1,li)
      z4 = z11(ji,ki+kinc1,li)
      x3 = x11(ji+jinc1,ki+kinc1,li)
      y3 = y11(ji+jinc1,ki+kinc1,li)
      z3 = z11(ji+jinc1,ki+kinc1,li)
      x5 = x11(ji,ki,li+linc1)
      y5 = y11(ji,ki,li+linc1)
      z5 = z11(ji,ki,li+linc1)
      x6 = x11(ji+jinc1,ki,li+linc1)
      y6 = y11(ji+jinc1,ki,li+linc1)
      z6 = z11(ji+jinc1,ki,li+linc1)
      x7 = x11(ji+jinc1,ki+kinc1,li+linc1)
      y7 = y11(ji+jinc1,ki+kinc1,li+linc1)
      z7 = z11(ji+jinc1,ki+kinc1,li+linc1)
      x8 = x11(ji,ki+kinc1,li+linc1)
      y8 = y11(ji,ki+kinc1,li+linc1)
      z8 = z11(ji,ki+kinc1,li+linc1)
c
      dx2 = x2 - x1
      dy2 = y2 - y1
      dz2 = z2 - z1
      dx3 = x3 - x1
      dy3 = y3 - y1
      dz3 = z3 - z1
      dx4 = x4 - x1
      dy4 = y4 - y1
      dz4 = z4 - z1
      dx5 = x5 - x1
      dy5 = y5 - y1
      dz5 = z5 - z1
      dx6 = x6 - x1
      dy6 = y6 - y1
      dz6 = z6 - z1
      dx7 = x7 - x1
      dy7 = y7 - y1
      dz7 = z7 - z1
      dx8 = x8 - x1
      dy8 = y8 - y1
      dz8 = z8 - z1
c
c     coefficients for tri-linear interpolation
c
      a2 = dx2
      a3 = dx4
      a4 = dx5
      a5 = dx3 - dx2 - dx4
      a6 = dx6 - dx2 - dx5
      a7 = dx8 - dx4 - dx5
      a8 = dx7 - dx3 - dx6 + dx2 - dx8 + dx4 + dx5
c
      b2 = dy2
      b3 = dy4
      b4 = dy5
      b5 = dy3 - dy2 - dy4
      b6 = dy6 - dy2 - dy5
      b7 = dy8 - dy4 - dy5
      b8 = dy7 - dy3 - dy6 + dy2 - dy8 + dy4 + dy5
c
      c2 = dz2
      c3 = dz4
      c4 = dz5
      c5 = dz3 - dz2 - dz4
      c6 = dz6 - dz2 - dz5
      c7 = dz8 - dz4 - dz5
      c8 = dz7 - dz3 - dz6 + dz2 - dz8 + dz4 + dz5
c
      xintrp = x1 + xie*( a2 + a5*eta + a6*zeta + a8*eta*zeta )
     .        + eta*( a3 + a7*zeta ) + a4*zeta 
      yintrp = y1 + xie*( b2 + b5*eta + b6*zeta + b8*eta*zeta )
     .        + eta*( b3 + b7*zeta ) + b4*zeta 
      zintrp = z1 + xie*( c2 + c5*eta + c6*zeta + c8*eta*zeta )
     .        + eta*( c3 + c7*zeta ) + c4*zeta 
c
      return
      end
      subroutine intpt( itr,jimage,kimage,limage,jd,kd,
     .                   ld,m1,xm1,ym1,zm1,jd1,kd1,ld1,m,i1)
c
c*********************************************************************** 
c     Purpose: find interpolation coefficients for interpolation to fringe 
c     or boundary points of mesh m from mesh m1
c***********************************************************************
c 
      include 'mag1.h'
c 
      common /intrp2/ jb(idim), kb(idim), lb(idim), jn(idim), kn(idim),
     .                ln(idim), itotal 
      common /surf2/  xbo(idim),ybo(idim),zbo(idim)
      common /where/  nblkpt(idim)
c 
      dimension xm1(jd1,kd1,ld1), ym1(jd1,kd1,ld1), zm1(jd1,kd1,ld1)
      dimension jimage(jd1,kd1,ld1),kimage(jd1,kd1,ld1),
     .          limage(jd1,kd1,ld1)
c
      idum1 = 0
      idum2 = 0
      idum3 = 0
      idum4 = 0
      dum1  = 0.
      dum2  = 0.
      dum3  = 0.
c
      if( itotal .gt. idim) then 
        write(6,601) m,idim,itotal
  601   format('0', ' dimension exceeded on interpolation arrays for' 
     .  ,1x,'mesh =',i5,2x,'idim =',i5,2x,'itotal = ',i5) 
        stop 'intpt' 
      end if
c
      call trace(itr,m1,idum2,idum3,idum4,dum1,dum2,dum3)
c
      icall = 0
c
      do 12 i = i1,itotal
c
      if(nblkpt(i).eq.0) then
c
        xp = xbo(i) 
        yp = ybo(i) 
        zp = zbo(i) 
c
c       search over mesh m1 to find the cell surrounding the mesh m 
c       point xp,yp,zp and compute the required interpolation data
c
        icall = icall + 1
        call search2(i,m,m1,jd1,kd1,ld1,xm1,ym1,zm1,jimage,kimage,
     .  limage,xp,yp,zp,iok,icall)
        if(iok.eq.1) nblkpt(i) = -(m1+mdim)
c
      end if
c      
   12 continue
c
      return
      end 
      subroutine obibl( m,n,iblank,jd,kd,ld )
c 
c***********************************************************************
c     Purpose: outer boundary points that are interpolated from another grid,
c     m1, are flagged by setting iblank = m1 +mdim. this will keep
c     mesh = 1 and iblank= 1 unique. these points will be 
c     set to zero later
c***********************************************************************
c 
      include 'mag1.h'
c 
      common /book1/  ipntr(mdim,mhldim*mdim), npntr(mdim,mhldim*mdim),
     &               mhbs(mdim,mdim), mobs(mdim,mdim), nsets(mdim)
      common /book2/  dxint(idim), dyint(idim), dzint(idim),
     .                ibpts(mdim), jbpt(idim), kbpt(idim), lbpt(idim),
     .                ji(idim),ki(idim), li(idim)
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
c 
      dimension iblank(jd,kd,ld)
c 
c 
        nnn = 1
        nserch=moutr(m,nnn)
      do 21 nn = 1,nserch
        nnn        = 1
        m1         = loutr(m,nnn,nn)
        ns         = mobs(m,m1) 
        ip         = ipntr(m1,ns) 
        np         = npntr(m1,ns) 
        call getint( m1,ji,ki,li,jbpt,kbpt,lbpt,
     &               dxint,dyint,dzint )
        do 11 ib = ip,np
          j             = jbpt(ib)
          k             = kbpt(ib)
          l             = lbpt(ib)
          if(j.ge.1 .and. j.le.jd .and.
     &       k.ge.1 .and. k.le.kd .and.
     &       l.ge.1 .and. l.le.ld) then
          iblank(j,k,l) = m1 +mdim
          end if
   11   continue
   21 continue
c 
      return
      end 
      subroutine obtot( m,x,y,z,iblank,jd,kd,ld,i1 ) 
c 
c***********************************************************************
c     Purpose: any points on the outer boundary of mesh m which need to 
c     interpolated from other meshes are loaded into the 
c     1-d lists jb,kb,lb.  outer boundary points that are also hole 
c     points of mesh m are then removed from these lists
c***********************************************************************
c 
      include 'mag1.h'
c 
      common /bound2/ iobtyp(mdim),nobtyp(mdim),
     .                jo1(ipmax),jo2(ipmax),ko1(ipmax),ko2(ipmax),
     .                lo1(ipmax),lo2(ipmax)
      common /intrp2/ jb(idim), kb(idim), lb(idim), jn(idim), kn(idim),
     .                ln(idim), itotal 
      common /surf2/  xbo(idim),ybo(idim),zbo(idim)
      common /pltpt/  xbnd(idim),ybnd(idim),zbnd(idim),xorph(idim),
     .                yorph(idim),zorph(idim),xill(idim),yill(idim),
     .                zill(idim),nbnd,norph,nill
c 
      dimension iblank(jd,kd,ld),ife(6)
      dimension x(jd,kd,ld),y(jd,kd,ld),z(jd,kd,ld)
c 
      do ll=1,6
         ife(ll) = 0
      end do
c
      nb       = nobtyp(m)
      write(6,*) '    creating outer boundary for mesh ',m,
     .' using ',nb,' coordinate surfaces '
c
      ito = 0
      do 1000 mb  = 1,nb
      ist = iobtyp(m) + mb - 1
      js  = jo1(ist)
      je  = jo2(ist)
      ks  = ko1(ist)
      ke  = ko2(ist)
      ls  = lo1(ist)
      le  = lo2(ist)
c
      write(6,*) '      coordinate surface number ',mb
      write(6,*) '      ......jo1,jo2,ko1,ko2,lo1,lo2 = ',
     .js,je,ks,ke,ls,le
c
      do 100 l=ls,le
      do 100 k=ks,ke
      do 100 j=js,je      
      ito      = ito + 1
      if( itotal+ito.gt.idim ) then
         write(6,602) ito+itotal,idim 
602      format('0',130('*'),/,/,5x,'failure in obtot: total' 
     .   ,1x,'number of boundary points exceeds available storage', 
     .   /,20x,'itotal = ',i6,2x,'idim = ',i6, 
     .   /,/,130('*') )
         stop 'obtot'
      end if
      jb(ito) = j
      kb(ito) = k
      lb(ito) = l 
      call bound(j,k,l,x1,y1,z1,jd,kd,ld,x,y,z,ife)
      xbo(ito) = x1
      ybo(ito) = y1
      zbo(ito) = z1
  100 continue
c
 1000 continue
c
c     flags to indicate extrapolated mesh m outer boundary points
c
      if(ife(1).gt.0) write(6,*) '      extrapolated position of',
     .' boundary at j=   1 boundary for ',ife(1),' points'
      if(ife(2).gt.0) write(6,*) '      extrapolated position of',
     .' boundary at j=jdim boundary for ',ife(2),' points'
      if(ife(3).gt.0) write(6,*) '      extrapolated position of',
     .' boundary at k=   1 boundary for ',ife(3),' points'
      if(ife(4).gt.0) write(6,*) '      extrapolated position of',
     .' boundary at k=kdim boundary for ',ife(4),' points'
      if(ife(5).gt.0) write(6,*) '      extrapolated position of',
     .' boundary at l=   1 boundary for ',ife(5),' points'
      if(ife(6).gt.0) write(6,*) '      extrapolated position of',
     .' boundary at l=ldim boundary for ',ife(6),' points'
c
      write(6,*)
      write(6,*) '    there are ',ito,
     .' boundary points in mesh ',m,' for which'
        write(6,*) '    stencils must be determined' 
c
c     exclude hole points from list of outer boundary points
c 
      ni = itotal
      i1 = itotal + 1
      do 1001 i = 1,ito 
      j  = max(1,jb(i))
      k  = max(1,kb(i))
      l  = max(1,lb(i))
      j  = min(jd,j)
      k  = min(kd,k)
      l  = min(ld,l)
      if( iblank(j,k,l).eq.1 )then 
        ni        = ni + 1 
        jb(ni)    = jb(i) 
        kb(ni)    = kb(i) 
        lb(ni)    = lb(i)
        xbo(ni)   = xbo(i)
        ybo(ni)   = ybo(i)
        zbo(ni)   = zbo(i)
        nbnd      = nbnd + 1
        if(nbnd .gt.idim) then
          write(6,*) '   total number of boundary points exceeds idim'
          stop 'obtot'
        end if
        xbnd(nbnd) = xbo(i)
        ybnd(nbnd) = ybo(i)
        zbnd(nbnd) = zbo(i)
c      else
c        write(8,*) ' mesh ',m,'  boundary point eliminated  i,j,k,l= ',
c     .  i,jb(i),kb(i),lb(i) 
      end if
 1001 continue
c 
      if(ni.ne.itotal+ito) then
        ito        = ni-itotal 
        write(6,*)'    ...after omitting hole points, there',
     .  ' are ',ito,' points for which'
        write(6,*)'    stencils must be determined'
      end if
c
      itotal = itotal +ito
c 
      return
      end 
      subroutine orphan(m,n,nserch,i1,iorph,iflg)
c 
c***********************************************************************
c     Purpose: deal with the inevitable orphan points: if iorph = 0, use
c     nearest neighboring point (zeroth order interpolation) from the 
c     meshes in the search list. if iorph = 1, extrapolate from the nearest 
c     neighbor. if iorph = 2, use nearest neighbor only if extrapolation
c     occurs outside the usually acceptable range of -0.5 to 1.5.
c***********************************************************************
c
      include 'mag1.h'
c 
      common /grid1 / x(mlen), y(mlen), z(mlen)
      common /grid2 / mjmax(mdim), mkmax(mdim), mlmax(mdim)
      common /grdlim/ jsrs(mdim), jsre(mdim), ksrs(mdim), ksre(mdim),
     .                lsrs(mdim), lsre(mdim) 
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
      common /surf2/  xbo(idim),ybo(idim),zbo(idim)
      common /intrp2/ jb(idim), kb(idim), lb(idim), jn(idim), kn(idim),
     .                ln(idim), itotal 
      common /intrp1/ xi(idim), yi(idim), zi(idim) 
      common /where/  nblkpt(idim)
      common /image/  jimage(mlen), kimage(mlen), limage(mlen) 
      common /tol/    epsc
c 
      dimension ds(idim),nbltmp(idim)
c
      do 220 i=i1,itotal
      nbltmp(i) = 0
220   ds(i)     = 1.e30
c
      do 22 nn = 1,nserch
      if (iflg.eq.0) then
        m1  = lhole(m,n,nn) 
      else
        nnn = 1
        m1 = loutr(m,nnn,nn)
      end if      
      jd1 = mjmax(m1)
      kd1 = mkmax(m1)
      ld1 = mlmax(m1)
      call getgrd( m1,x,y,z,jimage,kimage,limage,
     .jd1,kd1,ld1 ) 
      do 23 i=i1,itotal
      if(nblkpt(i).eq.0) then
        j  = jb(i) 
        k  = kb(i)
	l  = lb(i)
        xp = xbo(i)
        yp = ybo(i)
        zp = zbo(i)
c
c       minimum distance search on mesh m1 faces only  
c
        js = jsrs(m1) 
        je = jsre(m1) 
        ks = ksrs(m1) 
        ke = ksre(m1) 
        ls = lsrs(m1) 
        le = lsre(m1) 
        call dsmin2(jd1,kd1,ld1,x,y,z,
     .  xp,yp,zp,jp,kp,lp,js,je,ks,ke,ls,le,dmin1)
c
        if(dmin1.lt.ds(i)) then
          ds(i) = dmin1
c
          jpold = jp
          kpold = kp
          lpold = lp
c
          if(iorph.eq.0) then
c           use nearest neighbor (zeroth order interpolation)
            jn(i) = jp
            kn(i) = kp
            ln(i) = lp
            xi(i) = 0.
            yi(i) = 0.
            zi(i) = 0.
            iok   = 1
          else
c           extrapolate
c           keep within bounds 1 < j,k,l < j/k/ld1-1 for interp/extrap stencils
            jp = min0( jp , jd1-1 )
            kp = min0( kp , kd1-1 )
            lp = min0( lp , ld1-1 )
            jp = max0( 1 , jp )
            kp = max0( 1 , kp )
            lp = max0( 1 , lp )
c
c           save extrapolation coeffients from the last block to have been
c           checked, in case it turns out that the extrapolation coefficients
c           from the last block are smaller than the current one, despite the
c           smaller physical distance to a point in the current mesh
c
            isav = 0
            if(nbltmp(i) .ne. 0) then
              isav  = 1
              jnsav = jn(i)
              knsav = kn(i)
              lnsav = ln(i)
              xisav = xi(i)
              yisav = yi(i)
              zisav = zi(i)
              nbsav = nbltmp(i)
            end if
c
            call extrap(i,m,m1,jd1,kd1,ld1,x,y,z,xp,yp,zp,iok,
     .      jp,kp,lp,jimage,kimage,limage)
            if(iok.eq.0) then
c              linear extrapolation not successful...use zeroth order
               jn(i) = jpold
               kn(i) = kpold
               ln(i) = lpold
               xi(i) = 0.
               yi(i) = 0.
               zi(i) = 0.
             end if
c
c           for iorph=2, check extrapolation coefficients to
c           see if they are outside usually acceptable range 
c           of -0.5 to 1.5 . if outside acceptable range, use 
c           nearest neighbor instead
c
            if(iorph.eq.2) then
              if((xi(i).lt.-0.5-epsc .or. xi(i).gt.1.5+epsc).or.
     .           (yi(i).lt.-0.5-epsc .or. yi(i).gt.1.5+epsc).or.
     .           (zi(i).lt.-0.5-epsc .or. zi(i).gt.1.5+epsc)) then
                    if(iok.gt.0) then
                      jn(i) = jpold
                      kn(i) = kpold
                      ln(i) = lpold
                      xi(i) = 0.
                      yi(i) = 0.
                      zi(i) = 0.
                    end if
               end if
            end if
c 
            iok = 1
c
          end if
      	  if(iok.gt.0)nbltmp(i) = -(m1+mdim)
c
c         check current stencil against previous one for best (smallest) to use
c
          if(isav.gt.0) then
c
c           saved stencil
            if(xisav-1.0 .gt. epsc) then
              ximod0 = abs(xisav-1.0)
            else if(xisav .lt. -epsc) then
              ximod0 = abs(xisav)
            else
              ximod0 = 0.
            end if
            if(yisav-1.0 .gt. epsc) then
              yimod0 = abs(yisav-1.0)
            else if(yisav .lt. -epsc) then
              yimod0 = abs(yisav)
            else
              yimod0 = 0.
            end if
            if(zisav-1.0 .gt. epsc) then
              zimod0 = abs(zisav-1.0)
            else if(zisav .lt. -epsc) then
              zimod0 = abs(zisav)
            else
              zimod0 = 0.
            end if
            exmax0 = max(ximod0,yimod0,zimod0)
c
c           new stencil
            if(xi(i)-1.0 .gt. epsc) then
              ximod1 = abs(xi(i)-1.0)
            else if(xi(i) .lt. -epsc) then
              ximod1 = abs(xi(i))
            else
              ximod1 = 0.
            end if
            if(yi(i)-1.0 .gt. epsc) then
              yimod1 = abs(yi(i)-1.0)
            else if(yi(i) .lt. -epsc) then
              yimod1 = abs(yi(i))
            else
              yimod1 = 0.
            end if
            if(zi(i)-1.0 .gt. epsc) then
              zimod1 = abs(zi(i)-1.0)
            else if(zi(i) .lt. -epsc) then
              zimod1 = abs(zi(i))
            else
              zimod1 = 0.
            end if
            exmax1 = max(ximod1,yimod1,zimod1)
c
           if(exmax0 .lt. exmax1) then
c            old stencil was actually better
c***
            write(99,*) 'i =',i,' using stencil from mesh ',
     .       -(mdim+nbsav),' not current mesh ',m1
             write(99,*) '  old xi,yi,zi: ',xisav,yisav,zisav
             write(99,*) '  new xi,yi,zi: ',xi(i),yi(i),zi(i)
c***
             xi(i) = xisav
             yi(i) = yisav
             zi(i) = zisav
             jn(i) = jnsav
             kn(i) = knsav
             ln(i) = lnsav
             nbltmp(i) = nbsav 
           end if
          end if
        end if
      end if
23    continue
22    continue
c
      do 221 i=i1,itotal
      if(nblkpt(i).eq.0) nblkpt(i) = nbltmp(i)
221   continue         
c
      return
      end
c     
      subroutine outer
c
c***********************************************************************
c     Purpose: compute the interpolation data for the outer boundaries of 
c     mesh  m 
c***********************************************************************
c 
      include 'mag1.h'
c 
      common /book3/  iblank(mlen) 
      common /grid1/  x(mlen), y(mlen), z(mlen)
      common /grid2/  mjmax(mdim), mkmax(mdim), mlmax(mdim)
      common /surf2/  xbo(idim),ybo(idim),zbo(idim)
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
      common /intrp1/ xi(idim), yi(idim), zi(idim) 
      common /intrp2/ jb(idim), kb(idim), lb(idim), jn(idim), kn(idim),
     .                ln(idim), itotal 
      common /where/  nblkpt(idim)
      common /image/  jimage(mlen), kimage(mlen), limage(mlen) 
      common /tol/    epsc
      common /pltpt/  xbnd(idim),ybnd(idim),zbnd(idim),xorph(idim),
     .                yorph(idim),zorph(idim),xill(idim),yill(idim),
     .                zill(idim),nbnd,norph,nill
      common /trace1/ itrace
      common /snafu/  iholeh(mdim,mhldim),iholeo(mdim,mhldim),
     .                iorphh(mdim,mhldim),iorpho(mdim,mhldim)
c 
      dimension jbtmp(idim),kbtmp(idim),lbtmp(idim),
     .          jntmp(idim),kntmp(idim),lntmp(idim),xitmp(idim),
     .          yitmp(idim),zitmp(idim)
      dimension list(mdim),iskip(idim)
c  
      character*4 flag(3)
c
c     itrace < 0, do not write search history for current boundary point
c     itrace = 0, overwrite history from previous point with current 
c     itrace = 1, retain the search history for ALL points (may get huge file)
c     trace output found in unit 7
c
      itrace =  -1
c
      flag(1) = 'xie ' 
      flag(2) = 'eta ' 
      flag(3) = 'zeta'
c 
      do 31 m = 1,nmesh
c
      jd = mjmax(m) 
      kd = mkmax(m) 
      ld = mlmax(m)
      ntot = noutr(m) 
      itotal = 0
c
      if(ntot .gt. 0) then 
c
        do 610 mm=1,nmesh
        list(mm) = 0
610     continue
c
        write(6,*)
        write(6,*)
        write(6,*) '  *** beginning determination of stencils for',
     .  ' boundary points of mesh ',m,' ***'
        write(6,*)
c
c       loop over all outer boundaries of mesh m
c
        do 1111 n=1,ntot
        nnn = 1
        nserch=moutr(m,nnn)
        iorph = iorpho(m,nnn)
c
        call getgrd( m,x,y,z,jimage,kimage,limage,jd,kd,ld ) 
        call getibl( m,iblank,jd,kd,ld )
c
c       place outer boundary points of mesh m into 1-d lists jb,kb,lb of
c       points to be interpolated
c 
        call obtot( m,x,y,z,iblank,jd,kd,ld,i1 )
c 
c       initialize block pointer array nblkpt
c       ...nblkpt(i) = 0 if interp. stencil for boundary point i not yet found
c       ...nblkpt(i) = -(m1+mdim)  if interp. stencil for boundary point i 
c                       has been found in mesh m1
        do 112 i=i1,itotal
        nblkpt(i) = 0 
112     continue
c
c       search over all meshes connected to  mesh m to find interpolation
c       stencils for outer boundary points of mesh m
c
        do 21 nn = 1,nserch
        nnn = 1
        m1  = loutr(m,nnn,nn)
        jd1 = mjmax(m1) 
        kd1 = mkmax(m1) 
        ld1 = mlmax(m1) 
        call getgrd( m1,x,y,z,jimage,kimage,limage,jd1,kd1,ld1 ) 
        call intpt( 0,jimage,kimage,limage,jd,kd,ld, 
     .  m1,x,y,z,jd1,kd1,ld1,m,i1 ) 
c       reject stencils that contain hole/fringe pts.
          call getibl( m1,iblank,jd1,kd1,ld1 )
          do 1114 i=i1,itotal
          if(nblkpt(i) .eq. -(m1+mdim)) then
             j     = jn(i)
             k     = kn(i)
             l     = ln(i)
             jp1   = min( j+1,jd1 )
             kp1   = min( k+1,kd1 )
             lp1   = min( l+1,ld1 )
             ii1 = j   + (k-1)*jd1   + (l-1)*jd1*kd1
             ii2 = jp1 + (k-1)*jd1   + (l-1)*jd1*kd1
             ii3 = jp1 + (kp1-1)*jd1 + (l-1)*jd1*kd1
             ii4 = j   + (kp1-1)*jd1 + (l-1)*jd1*kd1
             ii5 = j   + (k-1)*jd1   + (lp1-1)*jd1*kd1
             ii6 = jp1 + (k-1)*jd1   + (lp1-1)*jd1*kd1
             ii7 = jp1 + (kp1-1)*jd1 + (lp1-1)*jd1*kd1
             ii8 = j   + (kp1-1)*jd1 + (lp1-1)*jd1*kd1
             if(iblank(ii1) .le. 0 .or.
     .          iblank(ii2) .le. 0 .or.
     .          iblank(ii3) .le. 0 .or.
     .          iblank(ii4) .le. 0 .or.
     .          iblank(ii5) .le. 0 .or.
     .          iblank(ii6) .le. 0 .or.
     .          iblank(ii7) .le. 0 .or.
     .          iblank(ii8) .le. 0) then
                  nblkpt(i) = 0
                  write(88,*)'in outer, rejecting pt. i = ',i,
     .            ' mesh ',m
             end if
          end if
 1114     continue
21      continue
c 
c       check to see if interpolation stencils were not found for any 
c       mesh m boundary points. for such "orphan" points, use either the nearest
c       point from one of the meshes in the search list, with xi=eta=zeta=0 
c       (zeroth order interpolation), or, use extrapolation from the nearest 
c       point.  extrapolation is distinguished from interpolation in that an
c       interpolation stencil has 0 < xie,eta,zeta < 1 (to within a tolerance
c       epsc) while for extrapolation, either xie, eta or zeta is < 0 or > 1
c       use of the nearest neighbor or extrapolation should occur only from 
c       boundaries of mesh m1; otherwise search routine has failed in some way
c       
        notok = 0
        do 113 i=i1,itotal
        if(nblkpt(i).eq.0) notok = notok + 1
113     continue
c
        if(notok.gt.0) then
          iflg = 1
          call orphan(m,n,nserch,i1,iorph,iflg)
        end if    
c
c       check to see if any interpolation stencils contain hole points 
c       if any do, check all other meshes in search list (if more than one)
c       to see if a valid stencil can be found elsewhere
c 
        nnn = 1      
        if(iholeo(m,nnn).gt.0) then
          iflg = 1
          call reserch(m,n,nserch,i1,iskip,iorph,iflg)
        end if
1111    continue 
c
c       generate list of all meshes which got searched to find the boundary 
c       points of all outer boundaries in the current mesh
c
        do 612 n = 1,ntot 
        nnn = 1
        nserch=moutr(m,nnn)
        do 612 nn=1,nserch
        nnn = 1
        m1 = loutr(m,nnn,nn)
        if(list(m1).eq.0) list(m1) = 1
612     continue
c  
c       summary of search routine results
c      
        write(6,*)
        write(6,*) '    summary of search routine results for mesh ',m,
     .  ' stencils:'
c
        nintot = 0
        nextot = 0
        nzerot = 0
c
        do 114 m1=1,nmesh
        if(list(m1).gt.0) then
          xtrap = 0.
          xtrap1 = 0.
          jd1    = mjmax(m1)
          kd1    = mkmax(m1)
          ld1    = mlmax(m1)
          nin    = 0
          nex    = 0
          norphb = 0
          norphc = 0
          norphd = 0
          nzero  = 0
          do 115 i=1,itotal
          if(nblkpt(i).eq.-(m1+mdim)) then
            j    = jb(i)
            k    = kb(i)
            l    = lb(i)
            xp   = xbo(i)
            yp   = ybo(i)
            zp   = zbo(i)
            jp   = jn(i)
            kp   = kn(i)
            lp   = ln(i)
            xie  = xi(i)
            eta  = yi(i)
            zeta = zi(i)
c
c           (orphan) points using zeroth order interpolation
c
            if(xie.eq.0. .and. eta.eq.0. .and. zeta.eq.0.)then
              nzero = nzero + 1
c             for stencils weighted totally at nearest neighbor,
c             keep all points in stencil within bounds of cell center grid,
c             grid,  1 .ge. j/k/l .le. j.k.ld1-1  
c             for example, if nearest neighbor is at last cell center in
c             l-direction, say, l=ld1, then switch stencil from
c             jn,kn,ln with xie,eta,zeta=0 to jn,kn,ln-1 with xie=eta=0,
c             and zeta=1. The two stencils are equivalent
c
              if(jn(i).eq.jd1) then
                jn(i) = jn(i) -1
                xi(i) = 1.
                jp    = jn(i)
                xie   = xi(i)
              end if
              if(kn(i).eq.kd1) then
                kn(i) = kn(i) -1
                yi(i) = 1.
                kp    = kn(i)
                eta   = yi(i)
              end if
              if(ln(i).eq.ld1) then
c               check for 2d case
                if(ld1.gt.1) then
                  ln(i) = ln(i)-1
                end if
                zi(i) = 1.
                lp    = ln(i)
                zeta  = zi(i)
              end if
c
              if(jn(i).eq.1.or.jn(i).eq.jd1-1 .or.
     .        kn(i).eq.1.or.kn(i).eq.kd1-1 .or.
     .        ln(i).eq.1.or.ln(i).eq.ld1-1) norphb = norphb + 1
c
              write(8,*) 'for the mesh ',m,' boundary point',
     .        '  j,k,l = ',j,k,l
              write(8,*) '  will use nearest point in',
     .        ' mesh ',m1,': j,k,l = ',jp,kp,lp
              write(8,*) '    with xie,eta,zeta = ',xie,eta,zeta
c
c           points using interpolation
c
            else if(xie .ge.-epsc .and. xie .le.1.+epsc .and. 
     .         eta .ge.-epsc .and. eta .le.1.+epsc .and. 
     .         zeta.ge.-epsc .and. zeta.le.1.+epsc) then 
               nin = nin + 1
            end if
c
c           (orphan) points using extrapolation
c
            if(xie .lt.-epsc .or. xie .gt.1.+epsc .or. 
     .      eta .lt.-epsc .or. eta .gt.1.+epsc .or. 
     .      zeta.lt.-epsc .or. zeta.gt.1.+epsc) then 
              iflagg = 0
              norph        = norph + 1
              xorph(norph) = xp
              yorph(norph) = yp
              zorph(norph) = zp
              nex = nex + 1
              if(jn(i).eq.1.or.jn(i).eq.jd1-1 .or.
     .        kn(i).eq.1.or.kn(i).eq.kd1-1 .or.
     .        ln(i).eq.1.or.ln(i).eq.ld1-1) norphb = norphb + 1
c
              write(8,*) 'for the mesh ',m,' boundary point',
     .        '  j,k,l = ',j,k,l
              write(8,*) '  will use extrapolation from',
     .        ' mesh ',m1,' point j,k,l = ',jp,kp,lp
              write(8,*) '    with xie,eta,zeta = ',xie,eta,zeta
c
              if(xie.lt.-epsc .or. xie.gt.1.+epsc) then
                iflagg = iflagg+1
                if(xie .gt. xtrap) then
                  xtrap  = xie
                  ixtrap = i
                  mxtrap = m1
                  iflag  = 1
                end if
                if(xie .lt. xtrap1) then
                  xtrap1  = xie
                  ixtrap1 = i
                  mxtrap1 = m1
                  iflag1  = 1
                end if
              end if
              if(eta.lt.-epsc .or. eta.gt.1.+epsc) then
                iflagg = iflagg+1
                if(eta .gt. xtrap) then
                  xtrap  = eta
                  ixtrap = i
                  mxtrap = m1
                  iflag  = 2
                end if
                if(eta .lt. xtrap1) then
                  xtrap1  = eta
                  ixtrap1 = i
                  mxtrap1 = m1
                  iflag1  = 2
                end if
              end if
              if(zeta.lt.-epsc .or. zeta.gt.1.+epsc) then
                iflagg = iflagg+1
                if(zeta .gt. xtrap) then
                  xtrap  = zeta
                  ixtrap = i
                  mxtrap = m1
                  iflag = 3
                end if
                if(zeta .lt. xtrap1) then
                  xtrap1  = zeta
                  ixtrap1 = i
                  mxtrap1 = m1
                  iflag1  = 3
                end if
              end if
              if(xie.lt.-0.5-epsc .or. xie.gt.1.5+epsc)
     .        norphc = norphc+1
              if(eta.lt.-0.5-epsc .or. eta.gt.1.5+epsc)
     .        norphc = norphc+1
              if(zeta.lt.-0.5-epsc .or. zeta.gt.1.5+epsc)
     .        norphc = norphc+1
              if(iflagg .gt. 1) norphd = norphd+1
            end if
          end if
115       continue
c
          write(6,*)
          write(6,*) '      ',nin,' boundary points of mesh ',m,
     .    ' are interpolated from mesh ',m1
          if(nex.gt.0) then
          write(6,*) '      ',nex,' orphaned boundary points of mesh ',
     .    m,' are extrapolated from mesh ',m1
          write(6,*) '      ',norphc,' of these orphans have',
     .    ' extrapolation coefficients <-0.5 or >1.5'
          if(norphd.gt.0) then
            write(6,*)'      ',norphd,' of these orphans',
     .      ' are extrapolated in more than one direction'
          end if
          if(xtrap.gt.1.+epsc) then
          write(6,*)'        the maximum extrapolation coefficient is ',
     .    flag(iflag),' = ',xtrap
          write(6,*)'        to the boundary pt ',jb(ixtrap),kb(ixtrap),
     .    lb(ixtrap),' from the mesh ',mxtrap,' target pt ',
     .    jn(ixtrap),kn(ixtrap),ln(ixtrap)
          end if
          if(xtrap1.lt.-epsc) then
          write(6,*)'        the minimum extrapolation coefficient is ',
     .    flag(iflag1),' = ',xtrap1
          write(6,*)'        to the boundary pt ',jb(ixtrap1),
     .    kb(ixtrap1),
     .    lb(ixtrap1),' from the mesh ',mxtrap1,' target pt ',
     .    jn(ixtrap1),kn(ixtrap1),ln(ixtrap1)
          end if
          write(6,*)'        check file 8 for more details' 
          end if
          if(nzero.gt.0) then
          write(6,*) '      ',nzero,' orphaned boundary points of ',
     .    'mesh ',m,' use the nearest point in mesh ',m1
          write(6,*)'        check file 8 for more details' 
          end if
          if(nex+nzero-norphb.ne.0)then
            write(6,*)'        WARNING: ',nex+nzero-norphb,' of these',
     .      ' orphans use extrap./zeroth order interp.'
            write(6,*)'          from an interior point of mesh ',m1,
     .      '...possible miscue in search routine'
          end if
c
          nintot = nintot + nin
          nextot = nextot + nex
          nzerot = nzerot + nzero
c
        end if 
c
114     continue
c
        if(nintot+nextot+nzerot .ne. itotal) then
          write(6,*)'  stopping...unable to compute ',
     .    itotal-nintot-nextot,' stencils for mesh ',m,
     .    ' boundary points'
          write(6,*)'     possible miscue in search routine'
          stop
        end if
c
c       put data back in original order (i.e. all points 
c       interpolated/extrapolated from mesh m1 are sequential in list jb,kb,lb
c
        itotl1 = itotal
        do 120 i=1,itotl1
        jbtmp(i) = jb(i)
        kbtmp(i) = kb(i)
        lbtmp(i) = lb(i)
        jntmp(i) = jn(i)
        kntmp(i) = kn(i)
        lntmp(i) = ln(i)
        xitmp(i) = xi(i)
        yitmp(i) = yi(i)
        zitmp(i) = zi(i)
120     continue
        do 125 m1=1,nmesh
        if(list(m1).gt.0) then
          ii = 0
          do 126 i=1,itotl1
          if(nblkpt(i) .eq. -(m1+mdim)) then
            ii = ii +1
            jb(ii) = jbtmp(i)
            kb(ii) = kbtmp(i)
            lb(ii) = lbtmp(i)
            jn(ii) = jntmp(i)
            kn(ii) = kntmp(i)
            ln(ii) = lntmp(i)
            xi(ii) = xitmp(i)
            yi(ii) = yitmp(i)
            zi(ii) = zitmp(i)
          end if
126       continue
c   
c         set interpolation pointers and load lists for input to flow solver
c
          itotal = ii
          icase  = 2
          call setptr( m,m1,icase )
        end if
125     continue
c
      end if 
c
31    continue
c 
c     update and store iblank array for mesh m (set iblank on outer 
c     boundary equal to m1 + mdim)
c 
      do 41 m = 1,nmesh  
      ntot = noutr(m)
      if( ntot.gt.0 ) then
        do 42 n=1,ntot
        jd = mjmax(m) 
        kd = mkmax(m) 
        ld = mlmax(m) 
        call getibl( m,iblank,jd,kd,ld )
        call obibl( m,n,iblank,jd,kd,ld )
        call putibl( m,iblank,jd,kd,ld )
   42   continue
      end if
   41 continue 
c
c     reset itrace to "off"
      itrace = -1
c
      return
      end 
c
      subroutine output 
c 
c***********************************************************************
c     Purpose: write summary of overlapped mesh calculation to unit 6, 
c     and  write interpolation data to unit 2
c***********************************************************************
c 
      include 'mag1.h'
c 
      character*20 titl
c
      common /book2/  dxint(idim), dyint(idim), dzint(idim),
     .                ibpts(mdim), jbpt(idim), kbpt(idim), lbpt(idim),
     .                ji(idim),ki(idim), li(idim)
      common /book3/  iblank(mlen) 
      common /grid2/  mjmax(mdim), mkmax(mdim), mlmax(mdim)
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
      common /intrp2/ jb(idim), kb(idim), lb(idim), jn(idim), kn(idim),
     .                ln(idim), itotal 
c 
      dimension ibc(idim), ibpnts(mdim), iipnts(mdim)
      dimension intpts(4),iord(idim)
c 
c***********************************************************************
c     print out summary of overlapped mesh calculations
c***********************************************************************
c 
c     count number of hole/fringe points in each mesh; get total for all meshes
c
      iholpt = 0
      ifrgpt = 0
      npnts  = 0
c 
      do 12 m = 1,nmesh 
c 
      jd = mjmax(m) 
      kd = mkmax(m) 
      ld = mlmax(m)  
      npnts = npnts + jd*kd*ld 
      call getibl( m,iblank,jd,kd,ld )
c
      do 11 l = 1,ld
      do 11 k = 1,kd
      do 11 j = 1,jd
      i       = j +(k-1)*jd +(l-1)*jd*kd
      if(iblank(i) .eq. 0) iholpt = iholpt + 1
      if(iblank(i) .lt. 0) ifrgpt = ifrgpt + 1
   11 continue
   12 continue
c
      write(6,*) 
      write(6,*) '  *** summary of maggie preprocessor ***'
      write(6,*)
      write(6,*) '    there are ', nmesh,' grids in the composite mesh'
      write(6,*) '    with a total of ',npnts,' points'
      write(6,*) '    of which: ',iholpt,' are hole points'
      write(6,*) '              ',ifrgpt,' are fringe points'
c
c***********************************************************************
c               set up CFL3D data structures 
c***********************************************************************
c  
c
c     set iblank for flow solver (remove connection information)
c 
      do 60 m = 1,nmesh
      jd=mjmax(m) 
      kd=mkmax(m) 
      ld=mlmax(m) 
      call getibl( m,iblank,jd,kd,ld )
c*** 
c     output of iblank arrays WITH connection info: call subroutine 
c     putibl2 to overwrite units 10 + mesh with this data
c     for use with checkmag2 program
      call putibl2( m,iblank,jd,kd,ld ) 
c*** 
      call iblcon( iblank,jd,kd,ld ) 
      call putibl( m,iblank,jd,kd,ld ) 
60    continue
c
      iitot = 0
      do 61 m = 1,nmesh 
      iitot = iitot + ibpts(m) 
   61 continue
c 
c     set up cross index array ibc,associated pointers
c     and load boundary point arrays
c 
      do 63 m = 1,nmesh 
c 
      write(6,*)
      write(6,*) '      writing connection data for CFL3D to ',
     .           ' unit 2 for mesh ',m
c
      if( m.eq.1 ) then 
        iisptr     = 1
        iieptr     = ibpts(m) 
      else
        iisptr     = iisptr + ibpts(m-1)
        iieptr     = iieptr + ibpts(m)
      end if
c
      call cindex( ibc,ibpnts(m),iipnts(m),jb,kb,lb,m ) 
c
      write(6,*) '        ibpnts,iipnts,iieptr,iisptr ',
     .                  ibpnts(m),iipnts(m),iieptr,iisptr
c
      call reorder(idim,ibpnts(m),intpts,jb,kb,lb,iord,
     .             mjmax(m),mkmax(m),mlmax(m)) 
c
      write(6,*) '        intpts = ',intpts
c       
      write(2) mjmax(m),mkmax(m),mlmax(m)
      write(2) ibpnts(m),intpts,iipnts(m),iieptr,iisptr
      call getint( m,ji,ki,li,jbpt,kbpt,lbpt, 
     .             dxint,dyint,dzint )
c 
      write(2) ( ji(i),ki(i),li(i),dxint(i),dyint(i), 
     .           dzint(i), i=1,iipnts(m) )
      write(2) ( jb(iord(i)),kb(iord(i)),lb(iord(i)),ibc(iord(i)),
     .           i=1,ibpnts(m) )
c 
      jd = mjmax(m) 
      kd = mkmax(m) 
      ld = mlmax(m) 
c
c     add iblank array to file 2, as required by CFL3D 
c
      call getibl(m,iblank,jd,kd,ld)
      call wiblnk(  iblank,jd,kd,ld)
c      call iblplt(m,iblank,jd,kd,ld)
c 
   63 continue
c 
c 
c***********************************************************************
c          compute minimum sizes of parameters needed for overlaped mesh
c          mesh option in CFL3D (minimum size of iitot) 
c***********************************************************************
c 
c 
      iidmax = 0
      ibdmax = 0
c 
      do 71 m = 1,nmesh 
      iidmax = iidmax + iipnts(m) 
      ibdmax = ibdmax + ibpnts(m)
   71 continue
      iitotl = max(iidmax,ibdmax)
c 
      write(6,*) 
      write(6,*) '    minimum dimension for parameter IITOT',
     .' in CFL3D (module RHS): ',iitotl 
c  
c 
c***********************************************************************
c          free up file space no longer needed
c***********************************************************************
c
c 
      do m=1,nmesh
         if (m.gt.99) then
             len = 12
             write (titl,'("temp_grd.",i3)') m
         else if (m.gt.9) then
             len = 11
             write (titl,'("temp_grd.",i2)') m
         else
             len = 10
             write (titl,'("temp_grd.",i1)') m
         endif
         do i = len+1, 20
             titl(i:i) = ' '
         end do
         iunit = 30
         open(iunit,file=titl(1:len),form='unformatted',
     .   status='unknown')
         close(iunit,status='delete')
      end do
c
      do m=1,nmesh
         if (m.gt.99) then
             len = 12
             write (titl,'("temp_cen.",i3)') m
         else if (m.gt.9) then
             len = 11
             write (titl,'("temp_cen.",i2)') m
         else
             len = 10
             write (titl,'("temp_cen.",i1)') m
         endif
         do i = len+1, 20
             titl(i:i) = ' '
         end do
         iunit = 30
         open(iunit,file=titl(1:len),form='unformatted',
     .   status='unknown')
         close(iunit,status='delete')
      end do
c
      do m=1,nmesh
         if (m.gt.99) then
             len = 12
             write (titl,'("temp_ibl.",i3)') m
         else if (m.gt.9) then
             len = 11
             write (titl,'("temp_ibl.",i2)') m
         else
             len = 10
             write (titl,'("temp_ibl.",i1)') m
         endif
         do i = len+1, 20
             titl(i:i) = ' '
         end do
         iunit = 30
         open(iunit,file=titl(1:len),form='unformatted',
     .   status='unknown')
         close(iunit,status='delete')
      end do
c
      do m=1,nmesh
         if (m.gt.99) then
             len = 12
             write (titl,'("temp_int.",i3)') m
         else if (m.gt.9) then
             len = 11
             write (titl,'("temp_int.",i2)') m
         else
             len = 10
             write (titl,'("temp_int.",i1)') m
         endif
         do i = len+1, 20
             titl(i:i) = ' '
         end do
         iunit = 30
         open(iunit,file=titl(1:len),form='unformatted',
     .   status='unknown')
         close(iunit,status='delete')
      end do 
c 
      return
      end 
      subroutine pltpts(nmesh,m,jmax,kmax,lmax,xx,yy,zz,ib,ibb)
c
c**********************************************************************
c     Purpose: output plot3d file for visual check of classification of
c     points, i.e., field, hole, fringe, interpolated boundary, orphan
c     (if applicable), or point with illegal stencil (if applicable) 
c**********************************************************************
c
      include 'mag1.h'
c
      common /pltpt/  xxbnd(idim),yybnd(idim),zzbnd(idim),xxorph(idim),
     .                yyorph(idim),zzorph(idim),xxill(idim),
     .                yyill(idim),zzill(idim),nbnd,norph,nill
      common /diagno/ jp(mdim*3+3),kp(mdim*3+3),lp(mdim*3+3),iplt3d
      common /grid2/  mjmax(mdim), mkmax(mdim), mlmax(mdim)
      common /igrdtyp/ ip3dgrd,ialph
c
      dimension xx(jmax,kmax,lmax),yy(jmax,kmax,lmax),
     .          zz(jmax,kmax,lmax)
      dimension ib(jmax,kmax,lmax),ibb(jmax,kmax,lmax)
c
      real*4 x(jdim,kdim,ldim),y(jdim,kdim,ldim),z(jdim,kdim,ldim),
     .       xbnd(idim),ybnd(idim),zbnd(idim),xorph(idim),
     .       yorph(idim),zorph(idim),xill(idim),
     .       yill(idim),zill(idim)
c
      do 1 j=1,jmax
      do 1 k=1,kmax
      do 1 l=1,lmax
      x(j,k,l) = xx(j,k,l)
      y(j,k,l) = yy(j,k,l)
      z(j,k,l) = zz(j,k,l)
1     continue
c
      iflag1 = 0
      iflag2 = 0
      iflag3 = 0
c
      if(m .eq. 1) then
        nmesh3     = 3*nmesh
        i = 0
        do 5 imesh=1,nmesh3,3
        i = i + 1
        jp(imesh) = mjmax(i)
        kp(imesh) = mkmax(i)
        lp(imesh) = mlmax(i)
5       continue
        do 2 imesh = 1,nmesh3,3
        jp(imesh+1) = jp(imesh)
        jp(imesh+2) = jp(imesh)
        kp(imesh+1) = kp(imesh)
        kp(imesh+2) = kp(imesh)
        lp(imesh+1) = lp(imesh)
        lp(imesh+2) = lp(imesh)
    2   continue
c
        if(nbnd.gt.0) then
          nmesh3 = nmesh3 + 1
          jp(nmesh3)=nbnd
          kp(nmesh3)=1
          lp(nmesh3)=1
        end if
        if(norph.gt.0) then
          nmesh3 = nmesh3 + 1
          jp(nmesh3)=norph
          kp(nmesh3)=1
          lp(nmesh3)=1
        end if
        if(nill.gt.0) then
          nmesh3 = nmesh3 + 1
          jp(nmesh3)=nill
          kp(nmesh3)=1
          lp(nmesh3)=1
        end if
c
        write(9) nmesh3
        write(9) (jp(i),kp(i),lp(i),i=1,nmesh3)
      end if
c
      ibl0=0
      ibl1=0
      ibln=0
c
c     output field points
c
      do 10 l=1,lmax
      do 10 k=1,kmax
      do 10 j=1,jmax
      iblk = ibb(j,k,l)
      if(iblk.eq.1) then
        ib(j,k,l) = 1
        ibl1 = ibl1 + 1
      else if(iblk.eq.0) then
        ib(j,k,l) = 0
        ibl0=ibl0+1
        else
        ib(j,k,l) = 0
        ibln=ibln+1
      end if
   10 continue
c
      if(ialph.eq.0) then
        write(9) (((x(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((y(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((z(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((ib(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
      else
        write(9) (((x(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((z(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((-y(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((ib(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
      end if
c
c     output hole points
c
      do 15 l=1,lmax
      do 15 k=1,kmax
      do 15 j=1,jmax
      iblk = ibb(j,k,l)
      if(iblk.eq.1) then
        ib(j,k,l) = 0
      else if(iblk.eq.0) then
        ib(j,k,l) = 1
      else
        ib(j,k,l) = 0
      end if
   15 continue
      if(ialph.eq.0) then
        write(9) (((x(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((y(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((z(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((ib(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
      else
        write(9) (((x(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((z(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((-y(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((ib(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
      end if
c
c     output fringe points
c
      do 20 l=1,lmax
      do 20 k=1,kmax
      do 20 j=1,jmax
      iblk = ibb(j,k,l)
      if(iblk.eq.0 .or. iblk.eq.1) then
        ib(j,k,l) = 0
      else
        ib(j,k,l) = 1
      end if
   20 continue
c
      if(ialph.eq.0) then
        write(9) (((x(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((y(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((z(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((ib(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
      else
        write(9) (((x(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((z(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((-y(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),
     .           (((ib(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
      end if
c
      if(m .eq. nmesh) then
c
c     output (interpolated) boundary points
c
        if(nbnd.gt.0) then
          do 1234 i=1,nbnd
          xbnd(i) = xxbnd(i)
          ybnd(i) = yybnd(i)
          zbnd(i) = zzbnd(i)
          ib(i,1,1) = 1
 1234     continue
          if(ialph.eq.0) then
            write(9) (xbnd(i),i=1,nbnd),
     .               (ybnd(i),i=1,nbnd),
     .               (zbnd(i),i=1,nbnd),
     .               (ib(i,1,1),i=1,nbnd)
          else
            write(9) (xbnd(i),i=1,nbnd),
     .               (zbnd(i),i=1,nbnd),
     .               (-ybnd(i),i=1,nbnd),
     .               (ib(i,1,1),i=1,nbnd)
          end if
        end if
c
c     output extrapolated (orphan) points
c
        if(norph.gt.0) then
          do 1235 i=1,norph
          xorph(i) = xxorph(i)
          yorph(i) = yyorph(i)
          zorph(i) = zzorph(i)
          ib(i,1,1) = 1
 1235     continue
          if(ialph.eq.0) then
            write(9) (xorph(i),i=1,norph),
     .               (yorph(i),i=1,norph),
     .               (zorph(i),i=1,norph),
     .               (ib(i,1,1),i=1,norph)
          else
            write(9) (xorph(i),i=1,norph),
     .               (zorph(i),i=1,norph),
     .               (-yorph(i),i=1,norph),
     .               (ib(i,1,1),i=1,norph)
          end if
        end if
c
c     output points with illegal interpolation stencils
c
        if(nill.gt.0) then
          do 1236 i=1,nill
          xill(i) = xxill(i)
          yill(i) = yyill(i)
          zill(i) = zzill(i)
          ib(i,1,1) = 1
 1236     continue
          if(ialph.eq.0) then
            write(9) (xill(i),i=1,nill),
     .               (yill(i),i=1,nill),
     .               (zill(i),i=1,nill),
     .               (ib(i,1,1),i=1,nill)
          else
            write(9) (xill(i),i=1,nill),
     .               (zill(i),i=1,nill),
     .               (-yill(i),i=1,nill),
     .               (ib(i,1,1),i=1,nill)
          end if
        end if
      end if
c
      return
      end
c
      subroutine putibl( m,iblank,jd,kd,ld ) 
c 
c***********************************************************************
c     Purpose: put a copy of iblank array for mesh m onto file temp_ibl.m
c***********************************************************************
c 
      character*20 titl
c
      dimension iblank(jd,kd,ld)
c
      iunit = 30
      if (m.gt.99) then
          len = 12
          write (titl,'("temp_ibl.",i3)') m
      else if (m.gt.9) then
          len = 11
          write (titl,'("temp_ibl.",i2)') m
      else
          len = 10
          write (titl,'("temp_ibl.",i1)') m
      endif
      do i = len+1, 20
          titl(i:i) = ' '
      end do
      open(iunit,file=titl(1:len),form='unformatted',
     .status='unknown')
      write(iunit) iblank 
      rewind iunit
c 
      return
      end 
      subroutine putibl2( m,iblank,jd,kd,ld ) 
c 
c***********************************************************************
c     Purpose: put a copy of iblank array for mesh m into file temp_grd.m
c***********************************************************************
c 
      character*20 titl
c
      dimension iblank(jd,kd,ld)
c
      iunit = 30
      if (m.gt.99) then
          len = 12
          write (titl,'("temp_grd.",i3)') m
      else if (m.gt.9) then
          len = 11
          write (titl,'("temp_grd.",i2)') m
      else
          len = 10
          write (titl,'("temp_grd.",i1)') m
      endif
      do i = len+1, 20
          titl(i:i) = ' '
      end do
      open(iunit,file=titl(1:len),form='unformatted',
     .status='unknown')
      write(iunit) (((iblank(j,k,l),j=1,jd),k=1,kd),l=1,ld) 
      rewind iunit
c 
      return
      end 
      subroutine putint( nwr,m,ji,ki,li,jbpt,kbpt,lbpt, 
     &                   dxint,dyint,dzint )
c 
c***********************************************************************
c     Purpose: put interpolation data for mesh m, into file temp_int.m
c***********************************************************************
c 
      include 'mag1.h'
c 
      character*20 titl
c
      dimension ji(idim),ki(idim),li(idim), 
     &          jbpt(idim),kbpt(idim),lbpt(idim), 
     &          dxint(idim),dyint(idim),dzint(idim) 
c
      iunit = 30
      if (m.gt.99) then
          len = 12
          write (titl,'("temp_int.",i3)') m
      else if (m.gt.9) then
          len = 11
          write (titl,'("temp_int.",i2)') m
      else
          len = 10
          write (titl,'("temp_int.",i1)') m
      endif
      do i = len+1, 20
          titl(i:i) = ' '
      end do
      open(iunit,file=titl(1:len),form='unformatted',
     .status='unknown')
      write(iunit) nwr
      if(nwr.gt.0) then
      write(iunit) (ji(i),ki(i),li(i),i=1,nwr),
     .       (jbpt(i),kbpt(i),lbpt(i),i=1,nwr),
     .    (dxint(i),dyint(i),dzint(i),i=1,nwr)
      end if
      rewind iunit
c 
c 
      return
      end 
c
      subroutine reorder(idim,ibpnts,intpts,jb,kb,lb,iord,jd,kd,ld)
c 
c***********************************************************************
c     Purpose: reorder and modify interpolated values for QI0/QJ0/QK0
c     as required by cfl3d
c***********************************************************************
c 
      dimension  jb(idim),  kb(idim),  lb(idim),   iord(idim)
c 
      dimension intpts(4)
c 
      do 1 ip    = 1,4
    1 intpts(ip) = 0
c
      iseq       = 0
      iint       = 0
      do 100 i = 1,ibpnts
      if(jb(i).le.jd .and. kb(i).le.kd .and. lb(i).le.ld .and.
     &   jb(i).ge.1  .and. kb(i).ge.1  .and. lb(i).ge.1) then
      iint       = iint + 1
      if(i.eq.iint) iseq = iint
      iord(iint) = i
      end if
  100 continue
      intpts(1)  = iint
      write(6,*) '        ...reordering.......interior pts = ',intpts(1)
      write(6,*) '        .....iseq,ibpnts = ',iseq,ibpnts
c
      ist        = iseq+1
      if(ist.gt.ibpnts .or. iint.ge.ibpnts) go to 1000
      do 200 i = ist,ibpnts
      if(jb(i).lt.1 .or. jb(i).gt.jd) then
      iint       = iint + 1
      if(i.eq.iint) iseq = iint
      iord(iint) = i
      end if
  200 continue
      intpts(2)   = iint - intpts(1)
      write(6,*) '        ...reordering.......QJ0 pts = ',intpts(2)
      write(6,*) '        ......iseq,ibpnts = ',iseq,ibpnts
c
      ist        = iseq+1
      if(ist.gt.ibpnts  .or. iint.ge.ibpnts) go to 1000
      do 300 i = ist,ibpnts
      if(kb(i).lt.1 .or. kb(i).gt.kd) then
      iint       = iint + 1
      if(i.eq.iint) iseq = iint
      iord(iint) = i
      end if
  300 continue
      intpts(3)   = iint - intpts(2) - intpts(1)
      write(6,*) '        ...reordering.......QK0 pts = ',intpts(3)
      write(6,*) '        ......iseq,ibpnts = ',iseq,ibpnts
c
      ist        = iseq+1
      if(ist.gt.ibpnts .or. iint.ge.ibpnts) go to 1000
      do 400 i = ist,ibpnts
      if(lb(i).lt.1 .or. lb(i).gt.ld) then
      iint       = iint + 1
      if(i.eq.iint) iseq = iint
      iord(iint) = i
      end if
  400 continue
      intpts(4)   = iint - intpts(3) - intpts(2) - intpts(1)
      write(6,*) '        ...reordering.......QI0 pts = ',intpts(4)
      write(6,*) '        ......iseq,ibpnts = ',iseq,ibpnts
c
 1000 continue
c
      if(ibpnts.ne.iint) then
      write(6,*) '    stopping in reorder......ibpnts.ne.iint',
     &   ibpnts,iint,intpts
      stop
      end if
c
      return
      end 
      subroutine reserch(m,n,nserch,i1,iskip,iorph,iflg)
c 
c***********************************************************************
c     Purpose: the initial pass through the search routines has determined
c     interpolation/extrapolation  stencils for all fringe or boundary points;
c     check to see if these stencils contain hole/fringe points (illegal)  
c***********************************************************************
c 
      include 'mag1.h'
      common /book3 / iblank(mlen) 
      common /grid1 / x(mlen), y(mlen), z(mlen)
      common /grid2 / mjmax(mdim), mkmax(mdim), mlmax(mdim)
      common /grdlim/ jsrs(mdim), jsre(mdim), ksrs(mdim), ksre(mdim),
     .                lsrs(mdim), lsre(mdim) 
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
      common /surf2/  xbo(idim),ybo(idim),zbo(idim)
      common /intrp2/ jb(idim), kb(idim), lb(idim), jn(idim), kn(idim),
     .                ln(idim), itotal 
      common /intrp1/ xi(idim), yi(idim), zi(idim) 
      common /where/  nblkpt(idim)
      common /image/  jimage(mlen), kimage(mlen), limage(mlen) 
      common /tol/    epsc
c 
      dimension ds(idim),nbltmp(idim)
      dimension iskip(idim)
c
          notok = 0
          do 1110 i=i1,itotal
          iskip(i)  = 9999
          nbltmp(i) = 0
1110      continue
          do 1113 nn = 1,nserch
      if (iflg.eq.0) then
        m1  = lhole(m,n,nn) 
      else
        nnn = 1
        m1 = loutr(m,nnn,nn)
      end if      
          jd1 = mjmax(m1)
          kd1 = mkmax(m1)
          ld1 = mlmax(m1)
          call getibl( m1,iblank,jd1,kd1,ld1 )
          do 1114 i=i1,itotal
          if(nblkpt(i) .ne. -(m1+mdim)) go to 1115
c
          j     = jn(i) 
          k     = kn(i) 
          l     = ln(i) 
          jp1   = min( j+1,jd1 )
          kp1   = min( k+1,kd1 )
          lp1   = min( l+1,ld1 )
          ii1 = j   + (k-1)*jd1   + (l-1)*jd1*kd1
          ii2 = jp1 + (k-1)*jd1   + (l-1)*jd1*kd1
          ii3 = jp1 + (kp1-1)*jd1 + (l-1)*jd1*kd1
          ii4 = j   + (kp1-1)*jd1 + (l-1)*jd1*kd1
          ii5 = j   + (k-1)*jd1   + (lp1-1)*jd1*kd1
          ii6 = jp1 + (k-1)*jd1   + (lp1-1)*jd1*kd1
          ii7 = jp1 + (kp1-1)*jd1 + (lp1-1)*jd1*kd1
          ii8 = j   + (kp1-1)*jd1 + (lp1-1)*jd1*kd1
c
          if(iblank(ii1) .le. 0 .or.
     .       iblank(ii2) .le. 0 .or.
     .       iblank(ii3) .le. 0 .or.
     .       iblank(ii4) .le. 0 .or.
     .       iblank(ii5) .le. 0 .or.
     .       iblank(ii6) .le. 0 .or.
     .       iblank(ii7) .le. 0 .or.
     .       iblank(ii8) .le. 0) then
c
c            stencil contains hole/fringe point from mesh m1: do not search 
c            mesh m1 in next attempt
c
             iskip(i) = m1
             notok = notok + 1
          end if
1115      continue
1114      continue
1113      continue
c
	  if(notok.gt.0 .and. nserch.gt.1) then
            write(6,*)
            write(6,*) '    initial pass through search routine ',
     .      'yielded ',notok,' illegal stencils'
            write(6,*) '    ...will attempt to',
     .      ' extrapolate/interpolate from another mesh'
            do 2221 i=i1,itotal
2221        ds(i) = 1.e30
c
c              
c           try again for those points which had illegal stencils
c
            nok = 0
            do 222 nn = 1,nserch
           if (iflg.eq.0) then
             m1  = lhole(m,n,nn) 
           else
             nnn = 1
             m1 = loutr(m,nnn,nn)
           end if      
            jd1 = mjmax(m1)
            kd1 = mkmax(m1)
            ld1 = mlmax(m1)
            call getibl( m1,iblank,jd1,kd1,ld1 )
            call getgrd( m1,x,y,z,jimage,kimage,limage,
     .      jd1,kd1,ld1 ) 
            do 233 i=i1,itotal
            if(iskip(i).eq.9999) go to 234
            if(iskip(i).eq.m1) go to 234
            j  = jb(i)
            k  = kb(i)
  	    l  = lb(i)
            xp = xbo(i)
            yp = ybo(i)
            zp = zbo(i)
            jnsav = jn(i)
            knsav = kn(i)
            lnsav = ln(i)
            xisav = xi(i)
            yisav = yi(i)
            zisav = zi(i)
c
c           minimum distance search 
c
            js = jsrs(m1) 
            je = jsre(m1) 
            ks = ksrs(m1) 
            ke = ksre(m1) 
            ls = lsrs(m1) 
            le = lsre(m1) 
            call dsmin2(jd1,kd1,ld1,x,y,z,
     .      xp,yp,zp,jp,kp,lp,js,je,ks,ke,ls,le,dmin1)
c***
c           may need to use entire mesh to find true nearest point
c           dsmin2 called above only searches boundaries
c           jskip = 1
c           kskip = 1
c           lskip = 1
c           call dsmin(jd1,kd1,ld1,x,y,z,xp,yp,zp,jp,kp,lp,
c     .     js,je,ks,ke,ls,le,jskip,kskip,lskip,dmin1)
c***
c     
            if(dmin1.lt.ds(i)) then
              ds(i) = dmin1
c
                jpold = jp
                kpold = kp
                lpold = lp
c
              if(iorph.eq.0) then
c               use nearest neighbor (zeroth order interpolation)
                jn(i) = jp
                kn(i) = kp
                ln(i) = lp
                xi(i) = 0.
                yi(i) = 0.
                zi(i) = 0.
                iok   = 1
              else
c               keep within bounds 1 < j,k,l < j/k/ld1-1 for interp/extrap stencils
                jp = min0( jp , jd1-1 )
                kp = min0( kp , kd1-1 )
                lp = min0( lp , ld1-1 )
                jp = max0( 1 , jp )
                kp = max0( 1 , kp )
                lp = max0( 1 , lp )
                call extrap(i,m,m1,jd1,kd1,ld1,x,y,z,xp,yp,zp,iok,
     .          jp,kp,lp,jimage,kimage,limage)
                if(iok.eq.0) then 
c                 linear extrapolation not successful...use zeroth order
                  jn(i) = jpold
                  kn(i) = kpold
                  ln(i) = lpold
                  xi(i) = 0.
                  yi(i) = 0.
                  zi(i) = 0.
                end if
c               for iorph=2, check extrapolation coefficients to
c               see if they are outside usually acceptable range
c               of -0.5 to 1.5 . if outside acceptable range, use
c               nearest neighbor instead
c
                if(iorph.eq.2) then
                  if((xi(i).lt.-0.5-epsc .or. xi(i).gt.1.5+epsc).or.
     .               (yi(i).lt.-0.5-epsc .or. yi(i).gt.1.5+epsc).or.
     .               (zi(i).lt.-0.5-epsc .or. zi(i).gt.1.5+epsc)) then
                     if(iok.gt.0) then
                       jn(i) = jpold
                       kn(i) = kpold
                       ln(i) = lpold
                       xi(i) = 0.
                       yi(i) = 0.
                       zi(i) = 0.
                     end if
                  end if
                end if
c
                iok = 1
              end if
	      if(iok.gt.0) then
                j     = jn(i) 
                k     = kn(i) 
                l     = ln(i) 
                jp1   = min( j+1,jd1 )
                kp1   = min( k+1,kd1 )
                lp1   = min( l+1,ld1 )
                ii1 = j   + (k-1)*jd1   + (l-1)*jd1*kd1
                ii2 = jp1 + (k-1)*jd1   + (l-1)*jd1*kd1
                ii3 = jp1 + (kp1-1)*jd1 + (l-1)*jd1*kd1
                ii4 = j   + (kp1-1)*jd1 + (l-1)*jd1*kd1
                ii5 = j   + (k-1)*jd1   + (lp1-1)*jd1*kd1
                ii6 = jp1 + (k-1)*jd1   + (lp1-1)*jd1*kd1
                ii7 = jp1 + (kp1-1)*jd1 + (lp1-1)*jd1*kd1
                ii8 = j   + (kp1-1)*jd1 + (lp1-1)*jd1*kd1
c
c               only accept new stencil if it does not contain hole/fringe pts.
c                                      
                if(iblank(ii1) .gt. 0 .and.
     .             iblank(ii2) .gt. 0 .and.
     .             iblank(ii3) .gt. 0 .and.
     .             iblank(ii4) .gt. 0 .and.
     .             iblank(ii5) .gt. 0 .and.
     .             iblank(ii6) .gt. 0 .and.
     .             iblank(ii7) .gt. 0 .and.
     .             iblank(ii8) .gt. 0) then
                   if(nbltmp(i).eq.0)nok = nok + 1
                   nblkpt(i) = -(m1+mdim)
                   nbltmp(i) = 1
                else
                   jn(i) = jnsav
                   kn(i) = knsav
                   ln(i) = lnsav
                   xi(i) = xisav
                   yi(i) = yisav
                   zi(i) = zisav
                end if
              end if
            end if
234         continue
233         continue
222         continue
            write(6,*) '       ',nok,' of these illegal stencils ',
     .      'were successfully replaced'
c
          end if
c
          return
          end
c
      subroutine rgrid( m,x,y,z,jd,kd,ld )
c 
c***********************************************************************
c     Purpose: read a grid from external unit = 10 (unit 10 is opened 
c     elsewhere in the code and contains all  grids in cfl3d format)
c     Also check the consistency of the input dimensions, and output
c     individual grids to file temp_grd.m, where m=grid number
c***********************************************************************
c 
      character*20 titl
c
      dimension x(jd,kd,ld), y(jd,kd,ld), z(jd,kd,ld) 
c
      read(10) jmax,kmax,lmax
c 
      if( jd.ne.jmax .or. kd.ne.kmax .or. ld.ne.lmax ) then
        write(6,606) m,jd,jmax,kd,kmax,ld,lmax 
  606   format('0',130('*'),/,/,5x,'grid input failed ' 
     .  ,'due to inconsistent input dimensions:',/,
     .  25x,'mesh number ',i5,/,
     .  25x,'jmax, tape jmax',2i5,/,
     .  25x,'kmax, tape kmax',2i5,/,
     .  25x,'lmax, tape lmax',2i5,/,/,130('*') ) 
        stop 'rgrid' 
      end if
c 
c     read grid from input file and write coordinates for each grid 
c     to a separate unit 
c 
      read(10) x,y,z
c
      iunit = 30
      if (m.gt.99) then
          len = 12
          write (titl,'("temp_grd.",i3)') m
      else if (m.gt.9) then
          len = 11
          write (titl,'("temp_grd.",i2)') m
      else
          len = 10
          write (titl,'("temp_grd.",i1)') m
      endif
      do i = len+1, 20
          titl(i:i) = ' '
      end do
      open(iunit,file=titl(1:len),form='unformatted',
     .status='unknown')
      write(iunit) x,y,z
      rewind iunit
c 
      return
      end  
c
      subroutine rp3d( m,x,y,z,jd,kd,ld,ialph,nmesh )
c
c***********************************************************************
c     Purpose: read a plot3d/tlns3d grid from external unit = 10 
c     (unit 10 is opened elsewhere in the code and in this case 
c     contains all grids in plot3d/tlns3d format);
c     Also check the consistency of the input dimensions, and output
c     individual grids to file temp_grd.m, where m=grid number
c***********************************************************************
c
      include 'mag1.h'
c
      character*20 titl
c
      dimension x(jd,kd,ld), y(jd,kd,ld), z(jd,kd,ld)
c
      common/dimen/ ltest(mdim),jtest(mdim),ktest(mdim)
c
c     ialph - flag to interpret input geometry in plot3d-type grid file
c             = 0 alpha measured in x-z plane (cfl3d standard)
c             > 0 alpha measured in x-y plane (tlns3d standard)
c
      ialph = 1
c
      if(m .eq. 1) then
        read(10) ngrid
        if (ngrid .ne. nmesh) then
           write(6,605) nmesh,ngrid
  605      format(1x,30h stopping, inconsistent input:,
     .     16h nmesh, ngrid = ,i3,i3)
           stop 'rp3d'
        end if
        read(10) (ltest(mm),jtest(mm),ktest(mm),mm=1,nmesh)
      end if
c
      if( jd.ne.jtest(m) .or. kd.ne.ktest(m) .or. 
     .  ld.ne.ltest(m) ) then
        write(6,606) m,jd,jtest(m),kd,ktest(m),ld,ltest(m)
  606   format('0',130('*'),/,/,5x,'grid input failed '
     .  ,'due to inconsistent input dimensions:',/,
     .  25x,'mesh number ',i5,/,
     .  25x,'jmax, tape jmax',2i5,/,
     .  25x,'kmax, tape kmax',2i5,/,
     .  25x,'lmax, tape lmax',2i5,/,/,130('*') )
        stop 'rp3d'
      end if
c
c     read grid from input file and write coordinates for each grid
c     to a separate unit
c
      if (ialph.eq.0) then
         read(10) (((x(j,k,l),l=1,ld),j=1,jd),k=1,kd),
     .            (((y(j,k,l),l=1,ld),j=1,jd),k=1,kd),
     .            (((z(j,k,l),l=1,ld),j=1,jd),k=1,kd)
      else
         read(10) (((x(j,k,l),l=1,ld),j=1,jd),k=1,kd),
     .            (((z(j,k,l),l=1,ld),j=1,jd),k=1,kd),
     .            (((y(j,k,l),l=1,ld),j=1,jd),k=1,kd)
         do 10 l=1,ld
         do 10 j=1,jd
         do 10 k=1,kd
         y(j,k,l) = -y(j,k,l)
   10    continue
      end if
c
      iunit = 30
      if (m.gt.99) then
          len = 12
          write (titl,'("temp_grd.",i3)') m
      else if (m.gt.9) then
          len = 11
          write (titl,'("temp_grd.",i2)') m
      else
          len = 10
          write (titl,'("temp_grd.",i1)') m
      endif
      do i = len+1, 20
          titl(i:i) = ' '
      end do
      open(iunit,file=titl(1:len),form='unformatted',
     .status='unknown')
      write(iunit) x,y,z
      rewind iunit
c
      return
      end  
c
      subroutine search2(i,m,m1,jd1,kd1,ld1,xm1,ym1,zm1,
     .                   jimage,kimage,limage,xp,yp,zp,iok,icall)     
c
c******************************************************************
c     Purpose: search over mesh m1 to find the cell surrounding the 
c     mesh m point xp,yp,zp and compute the required interpolation data
c******************************************************************
c
      include 'mag1.h'
c
      common /grdlim/ jsrs(mdim), jsre(mdim), ksrs(mdim), ksre(mdim),
     .                lsrs(mdim), lsre(mdim) 
      common /intrp1/ xi(idim), yi(idim), zi(idim) 
      common /intrp2/ jb(idim), kb(idim), lb(idim), jn(idim), kn(idim),
     .                ln(idim), itotal 
      common /where/  nblkpt(idim)
c 
      dimension xm1(jd1,kd1,ld1), ym1(jd1,kd1,ld1), zm1(jd1,kd1,ld1),
     .          jimage(jd1,kd1,ld1),kimage(jd1,kd1,ld1),
     .          limage(jd1,kd1,ld1)
c
      idum1 = 0
      idum2 = 0
      idum3 = 0
      idum4 = 0
      dum1  = 0.
      dum2  = 0.
      dum3  = 0.
c
      itmax  = 50
      limit  = 5
      icount = 1
c 
c     search limits for minimum distance search
c
      js = jsrs(m1) 
      je = jsre(m1) 
      ks = ksrs(m1) 
      ke = ksre(m1) 
      ls = lsrs(m1) 
      le = lsre(m1) 
c
17085 continue
c
      call trace(1,i,idum2,idum3,idum4,xp,yp,zp)
c
c
c     use last successful solution point in mesh m1 as the initial 
c     guess for the cell which surrounds  point i if one exists; otherwise,
c     use minimum distance search to get a starting point  
c
      if(icall.eq.1) then
        idsmin = 1
      else
        do 555 ii = 1,i-1
        iii = i - ii
        if(jn(iii).gt.0 .and. 
     .  abs(nblkpt(iii))-mdim .eq. m1) then
          xiet  = xi(iii)
          etat  = yi(iii)
          zetat = zi(iii)
          jp    = jn(iii)
          kp    = kn(iii)
          lp    = ln(iii)
          idsmin = 0
          go to 556
        end if
555       continue
          idsmin = 1
556       continue
      end if
c
      call topol(m1,jd1,kd1,ld1,xm1,ym1,zm1,jimage,kimage,limage,
     .     xp,yp,zp,xiet,etat,zetat,jp,kp,lp,itmax,
     .     limit,iok,js,je,ks,ke,ls,le,idsmin)
c
      if(iok.eq.1) then
c
c       search successful
c
        jn(i)  = jp
        kn(i)  = kp
        ln(i)  = lp
        xi(i)  = xiet
        yi(i)  = etat
        zi(i)  = zetat
      else
        call trace(20,idum1,idum2,idum3,idum4,dum1,dum2,dum3)
      end if
c
      return
      end     
      subroutine setptr( m,m1,icase ) 
c 
c***********************************************************************
c     Purpose: set interpolation pointers and load lists
c 
c     note that jn,kn,ln belong to mesh m1 and 
c     jb,kb,lb belong to mesh m
c     thus the data must be added to m1 data 
c 
c               icase  = 1 indicates hole boundary
c                      = 2 indicates outer boundary 
c***********************************************************************
c 
      include 'mag1.h'
c 
      common /book1/  ipntr(mdim,mhldim*mdim), npntr(mdim,mhldim*mdim),
     .                mhbs(mdim,mdim), mobs(mdim,mdim), nsets(mdim)
      common /book2/  dxint(idim), dyint(idim), dzint(idim),
     .                ibpts(mdim), jbpt(idim), kbpt(idim), lbpt(idim),
     .                ji(idim),ki(idim), li(idim)
      common /intrp1/ xi(idim), yi(idim), zi(idim) 
      common /intrp2/ jb(idim), kb(idim), lb(idim), jn(idim), kn(idim),
     .                ln(idim), itotal 
c
      write(6,*)
      write(6,*) '    updating pointers for mesh ',m
c 
c     update pointers
c 
      call getint( m1,ji,ki,li,jbpt,kbpt,lbpt,
     .             dxint,dyint,dzint )
      nsets(m1) = nsets(m1) + 1
      n = nsets(m1)
      if( n.gt.mhldim*mdim ) then
        write( 6,601) n,mhldim*mdim,m1 
601     format( '0',10x,'failure in setptr. the number of sets, n = ',
     .  i5,1x,'exceeds the max dimension = ',i5,1x,'for mesh ',i3 )
        stop 'setptr '
      end if
c 
      ipntr(m1,n) = ibpts(m1) + 1 
      npntr(m1,n) = ipntr(m1,n) + itotal - 1
      ibpts(m1)   = ibpts(m1) + itotal 
c
      write(6,*) '      nsets(m1),m1= ',n,m1
      write(6,*) '      ibpts,ipntr(m1,n),npntr(m1,n)= ',
     .             ibpts(m1),ipntr(m1,n),npntr(m1,n)
c 
      if( ibpts(m1).gt.idim ) then
        write(6,602) ibpts(m1),idim,m1
602     format( '0',10x,'failure in setup. the number of points,',
     .  ' ibpts =',i5,1x,'exceeds the max dimension = ',i5,
     .  1x,'for mesh ',i3 )
        stop 'setptr'
      end if
c 
      if( icase.eq.1 ) then 
c       hole boundaries
        mhbs(m,m1) = n
        write(6,*) '      m,m1,mhbs(m,m1)= ',m,m1,mhbs(m,m1)
      else if( icase.eq.2 ) then
c       outer boundaries 
        mobs(m,m1) = n
        write(6,*) '      m,m1,mobs(m,m1)= ',m,m1,mobs(m,m1)
      else
c       a foul up
        stop 'setptr' 
      end if
c 
c     update lists 
c 
      ib = ipntr(m1,n) - 1
      do 11 i = 1,itotal
      ib         = ib +  1 
      jbpt(ib)   = jb(i)
      kbpt(ib)   = kb(i)
      lbpt(ib)   = lb(i)
      ji(ib)     = jn(i)
      ki(ib)     = kn(i)
      li(ib)     = ln(i)
      dxint(ib)  = xi(i)
      dyint(ib)  = yi(i)
      dzint(ib)  = zi(i)
   11 continue
      call putint( ib,m1,ji,ki,li,jbpt,kbpt,lbpt,dxint,dyint,dzint ) 
c 
      return
      end 
      subroutine setup 
c 
c***********************************************************************
c     Purpose:  set up for overlaped grid interpolation
c***********************************************************************
c 
c 
      include 'mag1.h'
c 
      common /grid1/  x(mlen), y(mlen), z(mlen)
      common /grid2/  mjmax(mdim), mkmax(mdim), mlmax(mdim)
      common /conec1/ nmesh, nhole(mdim), noutr(mdim),
     .                mhole(mdim,mhldim),moutr(mdim,mhldim),
     .                lhole(mdim,mhldim,mdim),loutr(mdim,mhldim,mdim)
      common /image/  jimage(mlen), kimage(mlen), limage(mlen) 
      common /igrdtyp/ ip3dgrd,ialph
c 
c     read input data and initialize
c 
      call initia 
c
      do 11 mesh = 1,nmesh
c 
      jdg = mjmax(mesh)+1
      kdg = mkmax(mesh)+1
      ldg = mlmax(mesh)+1
c 
c     read in grid
c
      if(ip3dgrd.eq.0) then
c       (cfl3d format)
        call rgrid( mesh,x,y,z,jdg,kdg,ldg )
      else
c       (plot3d/tlns3d format)
        call rp3d( mesh,x,y,z,jdg,kdg,ldg,ialph,nmesh )
      end if
c 
      npnts = jdg*kdg*ldg 
c 
      if( npnts.gt.mlen ) then
        write(6,601) npnts, mlen
  601   format(1h0,130('*'),/,/,10x,'failure in setup due to number',
     .  1x,'of total grid points exceeding allocate',
     .  1x,'storage',/,30x,'npnts =',i3,5x,'mlen =',i3,/,/,130('*') ) 
        stop 'setup' 
      end if
c
c     check for branch cuts in grid
c
      if(mesh.eq.1) then
        write(6,*)
        write(6,*) '  *** beginning checks to identify branch',
     .  ' cuts ***'
      end if
      write(6,*)
      call branch(mesh,jdg,kdg,ldg,jimage,kimage,limage,
     .                 x,y,z)
c
c     transform grid to cell centered grid. 
c 
      call cellcen(x,y,z,jdg,kdg,ldg)
c 
c     write cell centered grid 
c 
      call wcelgrd(mesh,x,y,z,jimage,kimage,limage,jdg,kdg,ldg) 
c 
   11 continue
c
      rewind 1
c 
      npnts = 0
      do 21 m = 1,nmesh 
      npnts = (mjmax(m)+1)*(mkmax(m)+1)*(mlmax(m)+1) +npnts
   21 continue
c 
      return
      end 
      subroutine topol(m1,jdim,kdim,ldim,x,y,z,jimage,kimage,limage,
     .                 xp,yp,zp,xie,eta,zeta,jp,kp,lp,itmax,limit,iok,
     .                 js,je,ks,ke,ls,le,idsmin)
c
c***********************************************************************
c     Purpose: search mesh m1 for  the mesh m point with 
c     coordinates xp,yp,zp; determine corresponding xie,eta,zeta
c***********************************************************************
c
      implicit real (a-h,o-z)
c
c     parameter mfroz should be set at least as large as itmax
c
      parameter (mfroz=500)
c
      dimension x(jdim,kdim,ldim),y(jdim,kdim,ldim),
     .          z(jdim,kdim,ldim)
      dimension jimage(jdim,kdim,ldim),kimage(jdim,kdim,ldim),
     .          limage(jdim,kdim,ldim)
      dimension jfroz(mfroz),kfroz(mfroz),lfroz(mfroz)
c
      ifr = 1
      ihu = 1
c
      ifroze = 0
      ifroz  = 0
      ihuge  = 0
      iok    = 0
c
      idum1 = 0
      idum2 = 0
      idum3 = 0
      idum4 = 0
      dum1  = 0.
      dum2  = 0.
      dum3  = 0.
c
c     find nearest point (with indicies jp,kp,lp) to xp,yp,zp; start with
c     minimum distance search if idsmin > 0
c
c     search only on coarser version of  mesh 
      jskip = 4
      kskip = 4
      lskip = 4
      if(je-js.lt.4) jskip = 2
      if(ke-ks.lt.4) kskip = 2
      if(le-ls.lt.4) lskip = 2
c
      if(idsmin.gt.0) call dsmin(jdim,kdim,ldim,x,y,z,
     .xp,yp,zp,jp,kp,lp,js,je,ks,ke,ls,le,jskip,kskip,lskip,dmin)
c
c     keep within bounds 1 < j,k,l < j/k/ldim-1 for interp/extrap stencils
c
      jp = min0( jp , jdim-1 )
      kp = min0( kp , kdim-1 )
      lp = min0( lp , ldim-1 )
      jp = max0( 1 , jp )
      kp = max0( 1 , kp )
      lp = max0( 1 , lp )
c
      do 5555 intern=1,itmax
c
      jfroz(intern) = jp
      kfroz(intern) = kp
      lfroz(intern) = lp
c
      call trace(3,intern,idum2,idum3,idum4,dum1,dum2,dum3)
c
c     find local xie,eta,zeta via Newton iteraton in current target 
c     cell jp,kp,lp
c
      call trace(4,jp,kp,lp,m1,dum1,dum2,dum3)
c
      call xe(jdim,kdim,ldim,x,y,z,jp,kp,lp,xp,yp,zp, 
     .              xie,eta,zeta,imiss)
c
      call trace(5,idum1,idum2,idum3,idum4,xie,eta,zeta)
c
c     current target cell correct if imiss = 0
c
      if (imiss.eq.0) go to 5556
c
c     update current guess for target cell based on result of Newton
c     iteration, with max allowable change set by limit
c
      if (xie.ge.0) jinc = abs(xie)
      if (xie.lt.0) jinc = abs(xie-1)
      if (eta.ge.0) kinc = abs(eta)
      if (eta.lt.0) kinc = abs(eta-1)
      if (zeta.ge.0) linc = abs(zeta)
      if (zeta.lt.0) linc = abs(zeta-1)
c
      jinc = min0( jinc , limit )
      kinc = min0( kinc , limit )
      linc = min0( linc , limit )
c
      if (xie.gt.1.0) then
         jp = jp + jinc  
      else if (xie.lt.0.) then
         jp = jp - jinc  
      end if
      if (eta.gt.1.0) then
         kp = kp + kinc
      else if (eta.lt.0.) then
         kp = kp - kinc 
      end if
      if (zeta.gt.1.0) then
         lp = lp + linc
      else if (zeta.lt.0.) then
         lp = lp - linc 
      end if
c
      xieg  = float(jp) 
      etag  = float(kp) 
      zetag = float(lp) 
c
c     keep within bounds of mesh m1
c
      jp = min0( jp , jdim-1 )
      kp = min0( kp , kdim-1 )
      lp = min0( lp , ldim-1 )
      jp = max0( 1 , jp )
      kp = max0( 1 , kp )
      lp = max0( 1 , lp )
c
c     account for any branch cuts
c
      jpc = jp
      kpc = kp
      lpc = lp
c
      if(xieg .lt. 1. .or. xieg .gt.jdim-1) then
        jpc = jimage(jp,kp,lp)
        kpc = kimage(jp,kp,lp)
        lpc = limage(jp,kp,lp)
      end if
      if(etag .lt. 1. .or. etag .gt.kdim-1) then
        jpc = jimage(jp,kp,lp)
        kpc = kimage(jp,kp,lp)
        lpc = limage(jp,kp,lp)
      end if
      if(zetag.lt. 1. .or. zetag.gt.ldim-1) then
        jpc = jimage(jp,kp,lp)
        kpc = kimage(jp,kp,lp)
        lpc = limage(jp,kp,lp)
      end if
c
      if(jpc.ne.jp .or. kpc.ne.kp .or. lpc.ne.lp) then
        call trace(23,m1,jpc,kpc,lpc,dum1,dum2,dum3)
        call trace(24,idum1,jp,kp,lp,dum1,dum2,dum3)
      end if
c
      jp  = jpc
      kp  = kpc
      lp  = lpc

c
c     search routine off track if xie,eta or zeta are huge...try to 
c     get back on track with minumum distance search 
c
      huge=1.e3
      if(abs(xie).ge.huge .or. abs(eta).ge.huge .or. 
     .abs(zeta).ge.huge) then
        ihuge = ihuge + 1
        if (ihuge.gt.ihu) go to 1000
        if (ifroz.gt.0) go to 1000
        call trace(6,idum1,idum2,idum3,idum4,dum1,dum2,dum3)
        call dsmin(jdim,kdim,ldim,x,y,z,
     .       xp,yp,zp,jp,kp,lp,js,je,ks,ke,ls,le,jskip,kskip,
     .       lskip,dmin)
c
c       keep within bounds 1 < j,k,l < j/k/ldim-1 for interp/extrap stencils
c
        jp = min0( jp , jdim-1 )
        kp = min0( kp , kdim-1 )
        lp = min0( lp , ldim-1 )
        jp = max0( 1 , jp )
        kp = max0( 1 , kp )
        lp = max0( 1 , lp )
c
        go to 5555
      end if
c
c     check for frozen convergence: search routine keeps returning
c     to the same point, without 0 < xie,eta,zeta < 1 at that point. if frozen, 
c     attempt to break out of cycle by using minimum distance search
c
      ifroz = 0
      do 77 ii=1,intern
      int = intern-ii+1
      if (jp.eq.jfroz(int) .and. kp.eq.kfroz(int)
     .    .and. lfroz(int).eq.lp) then 
          ifroz  = 1
          jpfroz = jp
          kpfroz = kp
          lpfroz = lp
      end if
77    continue
      if (ifroz.eq.1) then
         ifroze = ifroze + 1
         if (ifroze.gt.ifr) go to 1000
         if (ihuge.gt.0) go to 1000
         call trace(7,jp,kp,lp,idum4,dum1,dum2,dum3)
         call dsmin(jdim,kdim,ldim,x,y,z,
     .        xp,yp,zp,jp,kp,lp,js,je,ks,ke,ls,le,jskip,kskip,
     .        lskip,dmin)
c
c        keep within bounds 1 < j,k,l < j/k/ldim-1 for interp/extrap stencils
c
         jp = min0( jp , jdim-1 )
         kp = min0( kp , kdim-1 )
         lp = min0( lp , ldim-1 )
         jp = max0( 1 , jp )
         kp = max0( 1 , kp )
         lp = max0( 1 , lp )
      end if
c
5555  continue
c 
c     search routine has been unsuccessful
c
1000  continue
      iok = 0
      return     
c
5556  continue
c
c     search routine has been successful
c
      iok = 1
      return
      end
      subroutine trace(icall,idum1,idum2,idum3,idum4,dum1,dum2,dum3)
c
c***********************************************************************
c     Purpose: Writes the search routine history for the current fringe 
c     or boundary point to unit 7.  
c***********************************************************************
c
      common /trace1/itrace
c
c     itrace < 0, do not write history for current fringe or boundary point
c     itrace = 0, overwrite history from previous point with current 
c     itrace = 1, retain the search history for ALL points 
c
      iunit = 7
c
      if(itrace.lt.0) return
c
      if(icall.eq.-1) then
      write(iunit,99)idum1
   99 format(' ',2x,'*** Search History For Fringe Points',
     .' ***',/,11x,'*** Searching Mesh ',i2,' ***')
      end if
c
      if(icall.eq.0) then
      write(iunit,*)
      write(iunit,*)
      write(iunit,*)
      write(iunit,100)idum1
  100 format(' ',2x,'*** Search History For Boundary Points',
     .' ***',/,12x,'*** Searching Mesh ',i2,' ***')
      end if
c
      if(icall.eq.1) then
        if(itrace.eq.0) then
          rewind(iunit)
          write(iunit,100)
        end if
      end if
c
      if(icall.eq.1) then
      write(iunit,101)idum1,dum1,dum2,dum3
  101 format(' ',2x,'iterations for interpolation to point  i =',i6,
     .       /,2x,'  with coordinates xp,yp,zp=',
     .       e11.4,',',e11.4,',',e11.4)
      end if
c
      if(icall.eq.11) then
      write(iunit,1101)idum1,dum1,dum2,dum3
 1101 format(' ',2x,'iterations for extrapolation to point  i =',i6,
     .       /,2x,'  with coordinates xp,yp,zp=',
     .       e11.4,',',e11.4,',',e11.4)
      end if
c
      if(icall.eq.3) then
      write(iunit,103) idum1
  103 format(' ',4x,'intern=',i2)
      end if
c
      if(icall.eq.4) then
      write(iunit,104) idum1,idum2,idum3,idum4
  104 format(' ',4x,'searching in cell j,k,i = ',i3,i3,i3,
     .       ' of block ',i3)
      end if
c
      if(icall.eq.5) then
      write(iunit,105)dum1,dum2,dum3
  105 format(' ',8x,'xie,eta,zeta=',e10.3,',',e10.3,
     .       ',',e10.3,' (local values)')
      end if
c     
      if(icall.eq.6) then
      write(iunit,*)'search off track...using dsmin to get back on',
     .              ' track'
      end if
c
      if(icall.eq.7) then
      write(iunit,*)'frozen convergence...using dsmin to attempt to',
     .              ' break cycle'
      end if
c
      if(icall.eq.17) then
      write(iunit,*)'frozen convergence...will extrapolate from best ',
     .'cell'
      end if
c
      if(icall.eq.20) then
      write(iunit,120)
  120 format(' ',4x,'giving up...will try in another mesh or',/,
     .       6x,'attempt extrapolation if that fails')
      end if
c 
      if(icall.eq.21) then
      write(iunit,*) '      miscue in xe2: newton iteration not',
     .' converged but 0 < xie,eta,zeta < 1'
      write(iunit,*) '      xie,eta,zeta = ',dum1,dum2,dum3
      end if
c
      if(icall.eq.22) then
      write(iunit,*) '      iteration = ',idum1,' error = ',dum1,
     .' error tolerance = ',dum2
      write(iunit,*) '      will set xie,eta,zeta to large value to',
     .'trigger min. distance search'
      end if
c
      if(icall.eq.23) then
        write(7,*)'   mesh ',idum1,' branch cut: jpc,kpc,lpc = ',
     .  idum2,idum3,idum4
      end if
c
      if(icall.eq.24) then
        write(7,*)'                        jp,kp,lp = ',
     .  idum2,idum3,idum4
      end if
c
      return
      end
      subroutine wcelgrd(m,x,y,z,jimage,kimage,limage,jd,kd,ld)
c 
c******************************************************************** 
c     Purpose: write the cell center coordinates for mesh m into
c     file temp_cen.m
c********************************************************************
c 
      character*20 titl
c
      dimension x(jd,kd,ld),y(jd,kd,ld),z(jd,kd,ld)
      dimension jimage(jd-1,kd-1,ld-1),kimage(jd-1,kd-1,ld-1),
     .limage(jd-1,kd-1,ld-1)
c
      iunit = 30
      if (m.gt.99) then
          len = 12
          write (titl,'("temp_cen.",i3)') m
      else if (m.gt.9) then
          len = 11
          write (titl,'("temp_cen.",i2)') m
      else
          len = 10
          write (titl,'("temp_cen.",i1)') m
      endif
      do i = len+1, 20
          titl(i:i) = ' '
      end do
      open(iunit,file=titl(1:len),form='unformatted',
     .status='unknown')
c 
      write(iunit)(((x(j,k,l),j=1,jd-1),k=1,kd-1),l=1,ld-1),
     &            (((y(j,k,l),j=1,jd-1),k=1,kd-1),l=1,ld-1),
     &            (((z(j,k,l),j=1,jd-1),k=1,kd-1),l=1,ld-1)
      write(iunit)(((jimage(j,k,l),j=1,jd-1),k=1,kd-1),l=1,ld-1),
     &            (((kimage(j,k,l),j=1,jd-1),k=1,kd-1),l=1,ld-1),
     &            (((limage(j,k,l),j=1,jd-1),k=1,kd-1),l=1,ld-1)
      rewind iunit 
c 
      return 
      end
      subroutine wiblnk(iblank,jd,kd,ld) 
c 
c******************************************************************** 
c     Purpose: write iblank array in form expected by cfl3d
c********************************************************************
c 
      dimension iblank(jd,kd,ld)
c 
      write(2)(((float(iblank(j,k,l)),j=1,jd),k=1,kd),l=1,ld)
c
      return
      end 
      subroutine xe(jdim,kdim,ldim,x,y,z,jcell,kcell,lcell,
     .              xc,yc,zc,xie,eta,zeta,imiss)
c
c***********************************************************************
c     Purpose: set up points for subroutine xe2
c***********************************************************************
c
      implicit real (a-h,o-z)
c
      dimension x(jdim,kdim,ldim),y(jdim,kdim,ldim),
     .          z(jdim,kdim,ldim)
c
c
c     index increments for 2d or 3d case
c
      jinc = 1
      kinc = 1
      linc = 1
      if(jdim.lt.2) jinc = 0
      if(kdim.lt.2) kinc = 0
      if(ldim.lt.2) linc = 0
c
      x1 = x(jcell,kcell,lcell)
      y1 = y(jcell,kcell,lcell)
      z1 = z(jcell,kcell,lcell)
      x2 = x(jcell+jinc,kcell,lcell)
      y2 = y(jcell+jinc,kcell,lcell)
      z2 = z(jcell+jinc,kcell,lcell)
      x4 = x(jcell,kcell+kinc,lcell)
      y4 = y(jcell,kcell+kinc,lcell)
      z4 = z(jcell,kcell+kinc,lcell)
      x3 = x(jcell+jinc,kcell+kinc,lcell)
      y3 = y(jcell+jinc,kcell+kinc,lcell)
      z3 = z(jcell+jinc,kcell+kinc,lcell)
      x5 = x(jcell,kcell,lcell+linc)
      y5 = y(jcell,kcell,lcell+linc)
      z5 = z(jcell,kcell,lcell+linc)
      x6 = x(jcell+jinc,kcell,lcell+linc)
      y6 = y(jcell+jinc,kcell,lcell+linc)
      z6 = z(jcell+jinc,kcell,lcell+linc)
      x7 = x(jcell+jinc,kcell+kinc,lcell+linc)
      y7 = y(jcell+jinc,kcell+kinc,lcell+linc)
      z7 = z(jcell+jinc,kcell+kinc,lcell+linc)
      x8 = x(jcell,kcell+kinc,lcell+linc)
      y8 = y(jcell,kcell+kinc,lcell+linc)
      z8 = z(jcell,kcell+kinc,lcell+linc)
c
      call xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,y1,y2,y3,y4,y5,y6,y7,y8,yc,
     .         z1,z2,z3,z4,z5,z6,z7,z8,zc,xie,eta,zeta,imiss,
     .         jinc,kinc,linc)
      return
      end
      subroutine xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,y1,y2,y3,y4,y5,y6,y7,y8,
     .               yc,z1,z2,z3,z4,z5,z6,z7,z8,zc,xie,eta,zeta,
     .               imiss,jinc,kinc,linc)
c
c***********************************************************************
c     Purpose: Set up coefficients for (locally) fitting a linear
c     variation in the xie, eta, and zeta directions .
c***********************************************************************
c
      common /tol/ epsc
c
      idum1 = 0
      idum2 = 0
      idum3 = 0
      idum4 = 0
      dum1  = 0.
      dum2  = 0.
      dum3  = 0.
c
c     initial guess
c
      xie  = .5
      eta  = .5 
      zeta = .5
c
      limit = 50
      iter  = 0
c
      dx2 = x2 - x1
      dy2 = y2 - y1
      dz2 = z2 - z1
      dx3 = x3 - x1
      dy3 = y3 - y1
      dz3 = z3 - z1
      dx4 = x4 - x1
      dy4 = y4 - y1
      dz4 = z4 - z1
      dx5 = x5 - x1
      dy5 = y5 - y1
      dz5 = z5 - z1
      dx6 = x6 - x1
      dy6 = y6 - y1
      dz6 = z6 - z1
      dx7 = x7 - x1
      dy7 = y7 - y1
      dz7 = z7 - z1
      dx8 = x8 - x1
      dy8 = y8 - y1
      dz8 = z8 - z1
c
c     coefficients for tri-linear representation
c
      a2 = dx2
      a3 = dx4
      a4 = dx5
      a5 = dx3 - dx2 - dx4
      a6 = dx6 - dx2 - dx5
      a7 = dx8 - dx4 - dx5
      a8 = dx7 - dx3 - dx6 + dx2 - dx8 + dx4 + dx5
c
      b2 = dy2
      b3 = dy4
      b4 = dy5
      b5 = dy3 - dy2 - dy4
      b6 = dy6 - dy2 - dy5
      b7 = dy8 - dy4 - dy5
      b8 = dy7 - dy3 - dy6 + dy2 - dy8 + dy4 + dy5
c
      c2 = dz2
      c3 = dz4
      c4 = dz5
      c5 = dz3 - dz2 - dz4
      c6 = dz6 - dz2 - dz5
      c7 = dz8 - dz4 - dz5
      c8 = dz7 - dz3 - dz6 + dz2 - dz8 + dz4 + dz5
c
c     2d cases
c
      if(jinc .eq. 0) then
        a2 = 1.
        a5 = 0.
        a6 = 0.
        a8 = 0.
        b2 = 1.
        b5 = 0.
        b6 = 0.
        b8 = 0.
        c2 = 1.
        c5 = 0.
        c6 = 0.
        c8 = 0.
      end if
      if(kinc .eq. 0) then
        a3 = 1.
        a5 = 0.
        a7 = 0.
        a8 = 0.
        b3 = 1.
        b5 = 0.
        b7 = 0.
        b8 = 0.
        c3 = 1.
        c5 = 0.
        c7 = 0.
        c8 = 0.
      end if
      if(linc .eq. 0) then
        a4 = 1.
        a6 = 0.
        a7 = 0.
        a8 = 0.
        b4 = 1.
        b6 = 0.
        b7 = 0.
        b8 = 0.
        c4 = 1.
        c6 = 0.
        c7 = 0.
        c8 = 0.
      end if
c
c     newton iteration to invert for xie,eta,zeta
c
c     convergence criterion for Newton iteration...require 4 orders of 
c     reduction to starting error=abs(f1)+abs(f2)+abs(f3) 
c
      f1 = x1 + xie*( a2 + a5*eta + a6*zeta + a8*eta*zeta )
     .        + eta*( a3 + a7*zeta ) + a4*zeta - xc
      f2 = y1 + xie*( b2 + b5*eta + b6*zeta + b8*eta*zeta )
     .        + eta*( b3 + b7*zeta ) + b4*zeta - yc
      f3 = z1 + xie*( c2 + c5*eta + c6*zeta + c8*eta*zeta )
     .        + eta*( c3 + c7*zeta ) + c4*zeta - zc      
c
      error0 = abs(f1)+abs(f2)+abs(f3)
      epsf   = 1.0e-4*error0
      epsf1  = 1.0e-6
      if(epsf .lt. epsf1) epsf = epsf1
c
    2 continue
c
      f1 = x1 + xie*( a2 + a5*eta + a6*zeta + a8*eta*zeta )
     .        + eta*( a3 + a7*zeta ) + a4*zeta - xc
      f2 = y1 + xie*( b2 + b5*eta + b6*zeta + b8*eta*zeta )
     .        + eta*( b3 + b7*zeta ) + b4*zeta - yc
      f3 = z1 + xie*( c2 + c5*eta + c6*zeta + c8*eta*zeta )
     .        + eta*( c3 + c7*zeta ) + c4*zeta - zc      
c
c d(f1)/d(xie)
      a2b = a2 + a5*eta + zeta*( a6 + a8*eta )
c
c d(f1)/d(eta)
      a3b = a3 + a5*xie + zeta*( a7 + a8*xie ) 
c
c d(f1)/d(zeta)
      a4b = a4 + a6*xie + eta*( a7 + a8*xie )
c
c d(f2)/d(xie)
      b2b = b2 + b5*eta + zeta*( b6 + b8*eta )
c
c d(f2)/d(eta)
      b3b = b3 + b5*xie + zeta*( b7 + b8*xie )
c
c d(f2)/d(zeta)
      b4b = b4 + b6*xie + eta*( b7 + b8*xie )
c
c d(f3)/d(xie)
      c2b = c2 + c5*eta + zeta*( c6 + c8*eta )
c
c d(f3)/d(eta)
      c3b = c3 + c5*xie + zeta*( c7 + c8*xie )
c
c d(f3)/d(zeta)
      c4b = c4 + c6*xie + eta*( c7 + c8*xie )
c
      iter  = iter + 1
c
      det   = 1./( a2b*(b3b*c4b-b4b*c3b)
     .           - a3b*(b2b*c4b-b4b*c2b)
     .           + a4b*(b2b*c3b-b3b*c2b) )
c
      xie   =  xie - det*( f1*(b3b*c4b-b4b*c3b)
     .                   - f2*(a3b*c4b-a4b*c3b) 
     .                   + f3*(a3b*b4b-a4b*b3b) )
c
      eta   =  eta - det*(-f1*(b2b*c4b-b4b*c2b)
     .                   + f2*(a2b*c4b-a4b*c2b)
     .                   - f3*(a2b*b4b-a4b*b2b) )
c
      zeta  = zeta - det*(f1*(b2b*c3b-b3b*c2b)
     .                  - f2*(a2b*c3b-a3b*c2b)
     .                  + f3*(a2b*b3b-a3b*b2b) )
c
c     exit newton iteration if xie,eta,zeta become large...the current target
c     cell is a poor guess. set xie,eta,zeta to a fixed large value so as to 
c     trigger a minimum distance search upon exiting this routine 

c
      huge=1.e10
      imiss = 0
      if(abs(xie) .gt. huge) then
        xie = huge
        imiss = 1
      end if
      if(abs(eta) .gt. huge) then
        eta = huge
        imiss = 1
      end if
      if(abs(zeta) .gt. huge) then
        zeta = huge
        imiss = 1
      end if
      if(imiss .gt. 0) return
c
      error = abs(f1) + abs(f2) +abs(f3)
      if (error.gt.epsf .and. iter.le.limit) go to 2
c
c     check to ensure point is inside cell
c
999   imiss = 0
      if (xie.lt.-epsc .or. xie.gt.1.+epsc .or.
     .    eta.lt.-epsc .or. eta.gt.1.+epsc .or.
     .    zeta.lt.-epsc .or. zeta.gt.1.+epsc) imiss = 1
c
c     check for non-converged newton iteration which nevertheless manages
c     to have  0 < xie,eta,zeta < 1. set xie,eta,zeta to a fixed large 
c     value so as to trigger a minimum distance search upon exiting this 
c     routine 
c
      if(imiss.eq.0 .and. error.gt.epsf) then
        call trace(21,idum1,idum2,idum3,idum4,xie,eta,zeta)
        call trace(22,iter,idum2,idum3,idum4,error,epsf,dum3)
        imiss = 1
        xie   = huge
        eta   = huge
        zeta  = huge
      end if
      return
      end
