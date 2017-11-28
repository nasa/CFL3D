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
      subroutine calyplus(jdim,kdim,idim,nbl,q,qi0,qj0,qk0,x,y,z,
     .                   vist3d,iover,bcj,bck,bci,sj,sk,si,yplusj,
     .                   yplusk,yplusi,blankj,blankk,blanki,
     .                   dnj,dnk,dni,vistj,vistk,visti,ifunc,iunit,
     .                   vj0,vk0,vi0,maxbl,maxseg,nbci0,nbcj0,
     .                   nbck0,nbcidim,nbcjdim,nbckdim,myid,myhost,
     .                   mycomm,mblk2nd,ibcinfo,jbcinfo,kbcinfo,vol)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Calculate yplus in turbulent flows at solid surfaces
c     in a block, and calculates statistics for the yplus distribution 
c     (average yplus, max yplus and location, standard deviation)
c
c     yplus values are calculated at surface grid point, first gridpoint
c     above the wall
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
      dimension intyp(10),relyp(6)
#endif
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension vist3d(jdim,kdim,idim)
      dimension sk(jdim,kdim,idim-1,5),sj(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
      dimension yplusj(kdim,idim,2),yplusk(jdim,idim,2),
     .          yplusi(jdim,kdim,2),blankj(kdim,idim,2),
     .          blankk(jdim,idim,2),blanki(jdim,kdim,2),
     .          dnj(kdim,idim,2),dnk(jdim,idim,2),dni(jdim,kdim,2),
     .          vistj(kdim,idim,2),vistk(jdim,idim,2),visti(jdim,kdim,2)
      dimension vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension mblk2nd(maxbl)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /reyue/ reue,tinf,ivisc(3)
      common /twod/ i2d
      common /wallfun/ iwf(3)
      common /ypinfo/ ypsumb,ypsumsqb,ypmaxb,ypminb,dnmaxb,dnminb,ypchk,
     .                nptsb,jypmaxb,kypmaxb,iypmaxb,jypminb,kypminb,
     .                iypminb,jdnmaxb,kdnmaxb,idnmaxb,jdnminb,
     .                kdnminb,idnminb,nypchkb
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      ypsumb   = 0.
      ypsumsqb = 0.
      ypmaxb   = 0.
      iypmaxb  = 0
      jypmaxb  = 0
      kypmaxb  = 0
      ypminb   = 1.e9
      iypminb  = 0
      jypminb  = 0
      kypminb  = 0
      dnmaxb   = 0.
      idnmaxb  = 0
      jdnmaxb  = 0
      kdnmaxb  = 0
      dnminb   = 1.e9
      idnminb  = 0
      jdnminb  = 0
      kdnminb  = 0
      nypchkb  = 0
      nptsb    = 0
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset     = maxbl
      itag_x      = 1
      itag_y      = itag_x      + ioffset
      itag_z      = itag_y      + ioffset
      itag_blank  = itag_z      + ioffset
      itag_vist   = itag_blank  + ioffset
      itag_yp     = itag_vist   + ioffset
      itag_dn     = itag_yp     + ioffset
      itag_int    = itag_dn     + ioffset
      itag_real   = itag_int    + ioffset
      itag_ypcomr = itag_real   + ioffset
      itag_ypcomi = itag_ypcomr + ioffset
      nn          = maxbl
c
#endif
      if (ivisc(3).gt.1) then
c
c        k=1 and/or k=kdim surfaces
c
         do 6000 m=1,2
c
         iflag   = 0
         ypsum   = 0.
         ypsumsq = 0.
         ypmax   = 0.
         iypmax  = 0
         jypmax  = 0
         ypmin   = 1.e9
         iypmin  = 0
         jypmin  = 0
         dnmax   = 0.
         idnmax  = 0
         jdnmax  = 0
         dnmin   = 1.e9
         idnmin  = 0
         jdnmin  = 0
         nypchk  = 0
         npts    = 0
c
         if (m.eq.1) then
            mm  = 2
            nbc = nbck0(nbl)
            kd  = 1
            kdx = 1
            kd1 = kd+1
         else 
            mm  = 4
            nbc = nbckdim(nbl)
            kd  = kdim1
            kdx = kdim
            kd1 = kd-1
         end if
c
         do 5990 i=1,idim
         do 5990 j=1,jdim
         blankk(j,i,m) = 0.
         yplusk(j,i,m) = 0.
         dnk(j,i,m)    = 0.
         vistk(j,i,m)  = 0.
 5990    continue
c
         do 6001 nseg = 1,nbc
c
         if (abs(kbcinfo(nbl,nseg,1,m)) .eq. 2004 .or.
     .       abs(kbcinfo(nbl,nseg,1,m)) .eq. 2014 .or.
     .       abs(kbcinfo(nbl,nseg,1,m)) .eq. 2024 .or.
     .       abs(kbcinfo(nbl,nseg,1,m)) .eq. 2034 .or.
     .       abs(kbcinfo(nbl,nseg,1,m)) .eq. 2016) then
c
            iflag = 1
c
c   do not include endpoints of segment, to avoid nonphysical behavior
            i1 = kbcinfo(nbl,nseg,2,m)+1
            i2 = kbcinfo(nbl,nseg,3,m)-1
            j1 = kbcinfo(nbl,nseg,4,m)+1
            j2 = kbcinfo(nbl,nseg,5,m)-1
            if (i1 .gt. i2) i1=i2
            if (j1 .gt. j2) j1=j2
c
#if defined DIST_MPI
            if (myid .ne. myhost) then
c
#endif
            do 6002 i=i1,i2
            id   = i
            id1  = id-1
            if (id1.lt.1)    id1 = 1
            if (id.gt.idim1) id = idim1
            do 6002 j=j1,j2
            jd   = j
            jd1  = jd-1
            if (jd1.lt.1)    jd1 = 1
            if (jd.gt.jdim1) jd = jdim1
c
            blankk(j,i,m) = 1.
c
c           wall value - average in J and I
c
            q1k1 = .25*(qk0(jd,id,1,mm)  + qk0(jd1,id,1,mm)
     .           +      qk0(jd,id1,1,mm) + qk0(jd1,id1,1,mm))
            q2k1 = .25*(qk0(jd,id,2,mm)  + qk0(jd1,id,2,mm)
     .           +      qk0(jd,id1,2,mm) + qk0(jd1,id1,2,mm))
            q3k1 = .25*(qk0(jd,id,3,mm)  + qk0(jd1,id,3,mm)
     .           +      qk0(jd,id1,3,mm) + qk0(jd1,id1,3,mm))
            q4k1 = .25*(qk0(jd,id,4,mm)  + qk0(jd1,id,4,mm)
     .           +      qk0(jd,id1,4,mm) + qk0(jd1,id1,4,mm))
            q5k1 = .25*(qk0(jd,id,5,mm)  + qk0(jd1,id,5,mm)
     .           +      qk0(jd,id1,5,mm) + qk0(jd1,id1,5,mm))*gamma
            t1k1 = q5k1/q1k1
c
c           first cell center location - average in J and I
c
            q1k2 = .25*(q(jd,kd,id,1)   + q(jd1,kd,id,1)
     .           +      q(jd,kd,id1,1)  + q(jd1,kd,id1,1))
            q2k2 = .25*(q(jd,kd,id,2)   + q(jd1,kd,id,2)
     .           +      q(jd,kd,id1,2)  + q(jd1,kd,id1,2))
            q3k2 = .25*(q(jd,kd,id,3)   + q(jd1,kd,id,3)
     .           +      q(jd,kd,id1,3)  + q(jd1,kd,id1,3))
            q4k2 = .25*(q(jd,kd,id,4)   + q(jd1,kd,id,4)
     .           +      q(jd,kd,id1,4)  + q(jd1,kd,id1,4))
            q5k2 = .25*(q(jd,kd,id,5)   + q(jd1,kd,id,5)
     .           +      q(jd,kd,id1,5)  + q(jd1,kd,id1,5))*gamma
            vist = .25*(vist3d(jd,kd,id)  + vist3d(jd1,kd,id)
     .           +      vist3d(jd,kd,id1) + vist3d(jd1,kd,id1))
            t1k2 = q5k2/q1k2
c
            if (sk(jd ,kdx,id ,4) .eq. 0. .or.
     +          sk(jd1,kdx,id ,4) .eq. 0. .or.
     +          sk(jd ,kdx,id1,4) .eq. 0. .or.
     +          sk(jd1,kdx,id1,4) .eq. 0.) then
            dx    = x(j,kd+1,i) - x(j,kd,i)
            dy    = y(j,kd+1,i) - y(j,kd,i)
            dz    = z(j,kd+1,i) - z(j,kd,i)
            dn    = sqrt(dx**2 + dy**2 + dz**2)
            else
            dn = (vol(jd ,kd,id )/sk(jd ,kdx,id ,4)+
     +            vol(jd1,kd,id )/sk(jd1,kdx,id ,4)+
     +            vol(jd ,kd,id1)/sk(jd ,kdx,id1,4)+
     +            vol(jd1,kd,id1)/sk(jd1,kdx,id1,4))/4.
            end if
            if (real(dn).gt.real(dnmax)) then
               dnmax  = dn
               idnmax = i
               jdnmax = j
            end if
            if (real(dn).lt.real(dnmin)) then
               dnmin  = dn
               idnmin = i
               jdnmin = j
            end if
            dnk(j,i,m)   = dn
            vistk(j,i,m) = vist
            if (real(dn).eq.0.) then
               dn = 1.
            end if
            qa1   = sqrt(q2k1**2 + q3k1**2 + q4k1**2)
            qa2   = sqrt(q2k2**2 + q3k2**2 + q4k2**2)
c        Get turb viscosity at wall (0 unless wall fn used)
            if (iwf(3) .eq. 1) then
              avgmut=.25*(vk0(jd,id,1,mm)  +vk0(jd1,id,1,mm)
     .                  + vk0(jd,id1,1,mm) +vk0(jd1,id1,1,mm))
            else
             avgmut=0.
            end if
            emuka = (t1k1**1.5)*((1.0+cbar/tinf)/(t1k1+cbar/tinf))
            pres1 = q5k1
            cfy   = 2.*(qa2-qa1)/dn
            cfy   = 2.*(emuka+avgmut)*cfy/(reue*xmach)
c
            yplusk(j,i,m) = 0.
            if (abs(real(cfy)).gt.0.) yplusk(j,i,m) = 
     .         dn*reue*q1k1*sqrt(abs(real(cfy)*0.5/real(q1k1)))/emuka
            ypsum   = ypsum  + yplusk(j,i,m)
            if (real(yplusk(j,i,m)).gt.real(ypmax)) then
               ypmax  = yplusk(j,i,m)
               iypmax = i
               jypmax = j
            end if
            if (real(yplusk(j,i,m)).gt.0. .and.
     .         real(yplusk(j,i,m)).lt.real(ypmin)) then
               ypmin  = yplusk(j,i,m)
               iypmin = i
               jypmin = j
            end if
            if (real(yplusk(j,i,m)).gt.real(ypchk)) then
               nypchk = nypchk+1
            end if
c
 6002       continue
c
            npts  = npts + (i2-i1+1)*(j2-j1+1)
#if defined DIST_MPI
c
            end if
#endif
c
c           2d output on segment by segment basis; 3d output on a
c           block face basis is done below
            if (ifunc.gt.0. .and. i2d.eq.1) then
#if defined DIST_MPI
c
               nxyz = idim*jdim*kdim
               nypk = jdim*idim*2
               nypj = kdim*idim*2
               nypi = jdim*kdim*2           
               if (myid.eq.mblk2nd(nbl)) then
                  mytag = itag_x + nbl
                  call MPI_Send(x,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_z + nbl
                  call MPI_Send(z,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_vist + nbl
                  call MPI_Send(vistk,nypk,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_yp + nbl
                  call MPI_Send(yplusk,nypk,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_dn + nbl
                  call MPI_Send(dnk,nypk,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
               end if
               if (myid.eq.myhost) then
                  nd_srce = mblk2nd(nbl)
                  mytag = itag_x + nbl
                  call MPI_Recv(x,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_z + nbl
                  call MPI_Recv(z,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_vist + nbl
                  call MPI_Recv(vistk,nypk,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_yp + nbl
                  call MPI_Recv(yplusk,nypk,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_dn + nbl
                  call MPI_Recv(dnk,nypk,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
               end if
c
               if (myid.eq.myhost) then
#endif
               write(iunit,*) 'zone t = "  ", i = ',j2-j1+1,
     .                        ', j = ',1,', k = ',1
               i      = 1
               do 6004 j=j1,j2
               write(iunit,402) real(x(j,kd,i)),real(z(j,kd,i)),
     .                          real(yplusk(j,i,m)),
     .                          real(dnk(j,i,m)),real(vistk(j,i,m))
 6004          continue
#if defined DIST_MPI
               end if
#endif
            end if
c
         end if
c
 6001    continue
c
         if (iflag.gt.0) then
c
#if defined DIST_MPI
            if (myid .ne. myhost) then
c
#endif
            ypavg = ypsum/npts
            do 6003 i=1,idim
            do 6003 j=1,jdim
            diff    = yplusk(j,i,m) - ypavg
            ypsumsq = ypsumsq + diff*diff
 6003       continue
            ypsdev  = sqrt(ypsumsq/(npts-1))
c
#if defined DIST_MPI
            end if
c
            if (myid.eq.mblk2nd(nbl)) then
               intyp(1) = jypmax
               intyp(2) = iypmax
               intyp(3) = jypmin
               intyp(4) = iypmin
               intyp(5) = jdnmax
               intyp(6) = idnmax
               intyp(7) = jdnmin
               intyp(8) = idnmin
               intyp(9) = nypchk
               intyp(10) = npts
               relyp(1) = ypmax
               relyp(2) = ypmin
               relyp(3) = dnmax
               relyp(4) = dnmin
               relyp(5) = ypavg
               relyp(6) = ypsdev
               mytag = itag_int + nbl
               call MPI_Send(intyp,10,MPI_INTEGER,myhost,
     .                       mytag,mycomm,ierr)
               mytag = itag_real + nbl
               call MPI_Send(relyp,6,MY_MPI_REAL,myhost,
     .                       mytag,mycomm,ierr)
            end if
            if (myid.eq.myhost) then
               nd_srce = mblk2nd(nbl)
               mytag = itag_int + nbl
               call MPI_Recv(intyp,10,MPI_INTEGER,nd_srce,
     .                       mytag,mycomm,istat,ierr)
               mytag = itag_real + nbl
               call MPI_Recv(relyp,6,MY_MPI_REAL,nd_srce,
     .                       mytag,mycomm,istat,ierr)
               jypmax = intyp(1)
               iypmax = intyp(2)
               jypmin = intyp(3)
               iypmin = intyp(4)
               jdnmax = intyp(5)
               idnmax = intyp(6)
               jdnmin = intyp(7)
               idnmin = intyp(8)
               nypchk = intyp(9)
               npts   = intyp(10)
               ypmax  = relyp(1)
               ypmin  = relyp(2)
               dnmax  = relyp(3)
               dnmin  = relyp(4)
               ypavg  = relyp(5)
               ypsdev = relyp(6)
            end if
c
            if (myid.eq.myhost) then
#endif
            if (m.eq.1) then
               write(11,101)
            else
               write(11,102)
            end if
            write(11,103)
            write(11,106) real(ypmax),jypmax,iypmax,
     .                    real(ypmin),jypmin,iypmin
            write(11,104)
            write(11,106) real(dnmax),jdnmax,idnmax,
     .                    real(dnmin),jdnmin,idnmin
            write(11,105) int(real(ypchk))
            write(11,107) real(ypavg),real(ypsdev),nypchk,npts
#if defined DIST_MPI
            end if
#endif
c
            if (ifunc.gt.0. .and. i2d.eq.0) then
#if defined DIST_MPI
c
               nxyz = idim*jdim*kdim
               nypk = jdim*idim*2
               nypj = kdim*idim*2
               nypi = jdim*kdim*2           
               if (myid.eq.mblk2nd(nbl)) then
                  mytag = itag_x + nbl
                  call MPI_Send(x,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_y + nbl
                  call MPI_Send(y,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_z + nbl
                  call MPI_Send(z,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_blank + nbl
                  call MPI_Send(blankk,nypk,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_vist + nbl
                  call MPI_Send(vistk,nypk,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_yp + nbl
                  call MPI_Send(yplusk,nypk,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_dn + nbl
                  call MPI_Send(dnk,nypk,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
               end if
               if (myid.eq.myhost) then
                  nd_srce = mblk2nd(nbl)
                  mytag = itag_x + nbl
                  call MPI_Recv(x,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_y + nbl
                  call MPI_Recv(y,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_z + nbl
                  call MPI_Recv(z,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_blank + nbl
                  call MPI_Recv(blankk,nypk,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_vist + nbl
                  call MPI_Recv(vistk,nypk,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_yp + nbl
                  call MPI_Recv(yplusk,nypk,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_dn + nbl
                  call MPI_Recv(dnk,nypk,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
               end if
c
               if (myid.eq.myhost) then
#endif
               kk = 1
               if (m.eq.2) kk = kdim
               write(iunit,*) ((real(x(j,kk,i)),i=1,idim),j=1,jdim),
     .                        ((real(y(j,kk,i)),i=1,idim),j=1,jdim),
     .                        ((real(z(j,kk,i)),i=1,idim),j=1,jdim),
     .               ((int(real(blankk(j,i,m))),i=1,idim),j=1,jdim)
               write(iunit+1,*) ((real(yplusk(j,i,m)),i=1,idim),
     .                            j=1,jdim)
               write(iunit+2,*) ((real(dnk(j,i,m)),i=1,idim),
     .                            j=1,jdim)
               write(iunit+3,*) ((real(vistk(j,i,m)),i=1,idim),
     .                            j=1,jdim)
#if defined DIST_MPI
               end if
#endif
            end if
#if defined DIST_MPI
c
            if (myid .ne. myhost) then
#endif
c
c           block totals
c
            ypsumb   = ypsumb   + ypsum
            ypsumsqb = ypsumsqb + ypsumsq
            nptsb    = nptsb    + npts
            nypchkb  = nypchkb  + nypchk
            if (real(ypmax) .gt. real(ypmaxb)) then
               ypmaxb   = ypmax
               iypmaxb  = iypmax
               jypmaxb  = jypmax
               if (m.eq.1) then
                  kypmaxb = 1
               else
                  kypmaxb = kdim
               end if
            end if
            if (real(ypmin) .lt. real(ypminb)) then
               ypminb   = ypmin
               iypminb  = iypmin
               jypminb  = jypmin
               if (m.eq.1) then
                  kypminb = 1
               else
                  kypminb = kdim
               end if
            end if
            if (real(dnmax) .gt. real(dnmaxb)) then
               dnmaxb   = dnmax
               idnmaxb  = idnmax
               jdnmaxb  = jdnmax
               if (m.eq.1) then
                  kdnmaxb = 1
               else
                  kdnmaxb = kdim
               end if
            end if
            if (real(dnmin) .lt. real(dnminb)) then
               dnminb   = dnmin
               idnminb  = idnmin
               jdnminb  = jdnmin
               if (m.eq.1) then
                  kdnminb = 1
               else
                  kdnminb = kdim
               end if
            end if
#if defined DIST_MPI
c
            end if
c
#endif
         end if
c
 6000    continue
      end if
c
      if (ivisc(2).gt.1) then
c
c        j=1 and/or j=jdim surfaces
c
         do 7000 m=1,2
c
         iflag   = 0
         ypsum   = 0.
         ypsumsq = 0.
         ypmax   = 0.
         iypmax  = 0
         kypmax  = 0
         ypmin   = 1.e9
         iypmin  = 0
         kypmaxin= 0
         dnmax   = 0.
         idnmax  = 0
         kdnmax  = 0
         dnmin   = 1.e9
         idnmin  = 0
         kdnmin  = 0
         nypchk  = 0
         npts    = 0
c
         if (m.eq.1) then
            mm  = 2
            nbc = nbcj0(nbl)
            jd  = 1
            jdx = 1
            jd1 = jd+1
         else
            mm  = 4
            nbc = nbcjdim(nbl)
            jd  = jdim1
            jdx = jdim
            jd1 = jd-1
         end if
c
         do 6990 i=1,idim
         do 6990 k=1,kdim
         blankj(k,i,m) = 0.
         yplusj(k,i,m) = 0.
         dnj(k,i,m)    = 0.
         vistj(k,i,m)  = 0.
 6990    continue
c
         do 7001 nseg = 1,nbc
c
         if (abs(jbcinfo(nbl,nseg,1,m)) .eq. 2004 .or.
     .       abs(jbcinfo(nbl,nseg,1,m)) .eq. 2014 .or.
     .       abs(jbcinfo(nbl,nseg,1,m)) .eq. 2024 .or.
     .       abs(jbcinfo(nbl,nseg,1,m)) .eq. 2034 .or.
     .       abs(jbcinfo(nbl,nseg,1,m)) .eq. 2016) then
c
            iflag = 1
c
c   do not include endpoints of segment, to avoid nonphysical behavior
            i1 = jbcinfo(nbl,nseg,2,m)+1
            i2 = jbcinfo(nbl,nseg,3,m)-1
            k1 = jbcinfo(nbl,nseg,4,m)+1
            k2 = jbcinfo(nbl,nseg,5,m)-1
            if (i1 .gt. i2) i1=i2
            if (k1 .gt. k2) k1=k2
c
#if defined DIST_MPI
            if (myid .ne. myhost) then
c
#endif
            do 7002 i=i1,i2
            id   = i
            id1  = id-1
            if (id1.lt.1)    id1 = 1
            if (id.gt.idim1) id = idim1
            do 7002 k=k1,k2
            kd   = k
            kd1  = kd-1
            if (kd1.lt.1)    kd1 = 1
            if (kd.gt.kdim1) kd = kdim1
c
            blankj(k,i,m) = 1.
c
c           wall value - average in K and I
c
            q1j1 = .25*(qj0(kd,id,1,mm)  + qj0(kd1,id,1,mm)
     .           +      qj0(kd,id1,1,mm) + qj0(kd1,id1,1,mm))
            q2j1 = .25*(qj0(kd,id,2,mm)  + qj0(kd1,id,2,mm)
     .           +      qj0(kd,id1,2,mm) + qj0(kd1,id1,2,mm))
            q3j1 = .25*(qj0(kd,id,3,mm)  + qj0(kd1,id,3,mm)
     .           +      qj0(kd,id1,3,mm) + qj0(kd1,id1,3,mm))
            q4j1 = .25*(qj0(kd,id,4,mm)  + qj0(kd1,id,4,mm)
     .           +      qj0(kd,id1,4,mm) + qj0(kd1,id1,4,mm))
            q5j1 = .25*(qj0(kd,id,5,mm)  + qj0(kd1,id,5,mm)
     .           +      qj0(kd,id1,5,mm) + qj0(kd1,id1,5,mm))*gamma
            t1j1 = q5j1/q1j1
c
c           first cell center location - average in K and I
c
            q1j2 = .25*(q(jd,kd,id,1)   + q(jd,kd1,id,1)
     .           +      q(jd,kd,id1,1)  + q(jd,kd1,id1,1))
            q2j2 = .25*(q(jd,kd,id,2)   + q(jd,kd1,id,2)
     .           +      q(jd,kd,id1,2)  + q(jd,kd1,id1,2))
            q3j2 = .25*(q(jd,kd,id,3)   + q(jd,kd1,id,3)
     .           +      q(jd,kd,id1,3)  + q(jd,kd1,id1,3))
            q4j2 = .25*(q(jd,kd,id,4)   + q(jd,kd1,id,4)
     .           +      q(jd,kd,id1,4)  + q(jd,kd1,id1,4))
            q5j2 = .25*(q(jd,kd,id,5)   + q(jd,kd1,id,5)
     .           +      q(jd,kd,id1,5)  + q(jd,kd1,id1,5))*gamma
            vist = .25*(vist3d(jd,kd,id)  + vist3d(jd,kd1,id)
     .           +      vist3d(jd,kd,id1) + vist3d(jd,kd1,id1))
            t1j2 = q5j2/q1j2
c
            if (sj(jdx,kd ,id ,4) .eq. 0. .or.
     +          sj(jdx,kd1,id ,4) .eq. 0. .or.
     +          sj(jdx,kd ,id1,4) .eq. 0. .or.
     +          sj(jdx,kd1,id1,4) .eq. 0.) then
            dx    = x(jd+1,k,i) - x(jd,k,i)
            dy    = y(jd+1,k,i) - y(jd,k,i)
            dz    = z(jd+1,k,i) - z(jd,k,i)
            dn    = sqrt(dx**2 + dy**2 + dz**2)
            else
            dn = (vol(jd,kd ,id )/sj(jdx,kd ,id ,4)+
     +            vol(jd,kd1,id )/sj(jdx,kd1,id ,4)+
     +            vol(jd,kd ,id1)/sj(jdx,kd ,id1,4)+
     +            vol(jd,kd1,id1)/sj(jdx,kd1,id1,4))/4.
            end if
            if (real(dn).gt.real(dnmax)) then
               dnmax  = dn
               idnmax = i
               kdnmax = k
            end if
            if (real(dn).lt.real(dnmin)) then
               dnmin  = dn
               idnmin = i
               kdnmin = k
            end if
            dnj(k,i,m)   = dn
            vistj(k,i,m) = vist
            if (real(dn).eq.0.) then
               dn = 1.
            end if
            qa1   = sqrt(q2j1**2 + q3j1**2 + q4j1**2)
            qa2   = sqrt(q2j2**2 + q3j2**2 + q4j2**2)
c        Get turb viscosity at wall (0 unless wall fn used)
            if (iwf(2) .eq. 1) then
              avgmut=.25*(vj0(kd,id,1,mm)  +vj0(kd1,id,1,mm)
     .                  + vj0(kd,id1,1,mm) +vj0(kd1,id1,1,mm))
            else
             avgmut=0.
            end if
            emuka = (t1j1**1.5)*((1.0+cbar/tinf)/(t1j1+cbar/tinf))
            pres1 = q5j1
            cfy   = 2.*(qa2-qa1)/dn
            cfy   = 2.*(emuka+avgmut)*cfy/(reue*xmach)
c
            yplusj(k,i,m) = 0.
            if (abs(real(cfy)).gt.0.) yplusj(k,i,m) = 
     .         dn*reue*q1j1*sqrt(abs(real(cfy)*0.5/real(q1j1)))/emuka
            ypsum   = ypsum   + yplusj(k,i,m)
            if (real(yplusj(k,i,m)).gt.real(ypmax)) then
               ypmax  = yplusj(k,i,m)
               iypmax = i
               kypmax = k
            end if
            if (real(yplusj(k,i,m)).gt.0. .and.
     .         real(yplusj(k,i,m)).lt.real(ypmin)) then
               ypmin  = yplusj(k,i,m)
               iypmin = i
               kypmin = k
            end if
            if (real(yplusj(k,i,m)).gt.real(ypchk)) then
               nypchk = nypchk+1
            end if
c
 7002       continue
c
            npts  = npts + (i2-i1+1)*(k2-k1+1)
#if defined DIST_MPI
c
            end if
#endif
c
c           2d output on segment by segment basis; 3d output on a
c           block face basis is done below
            if (ifunc.gt.0. .and. i2d.eq.1) then
#if defined DIST_MPI
c
               nxyz = idim*jdim*kdim
               nypk = jdim*idim*2
               nypj = kdim*idim*2
               nypi = jdim*kdim*2           
               if (myid.eq.mblk2nd(nbl)) then
                  mytag = itag_x + nbl
                  call MPI_Send(x,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_z + nbl
                  call MPI_Send(z,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_vist + nbl
                  call MPI_Send(vistj,nypj,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_yp + nbl
                  call MPI_Send(yplusj,nypj,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_dn + nbl
                  call MPI_Send(dnj,nypj,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
               end if
               if (myid.eq.myhost) then
                  nd_srce = mblk2nd(nbl)
                  mytag = itag_x + nbl
                  call MPI_Recv(x,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_z + nbl
                  call MPI_Recv(z,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_vist + nbl
                  call MPI_Recv(vistj,nypj,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_yp + nbl
                  call MPI_Recv(yplusj,nypj,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_dn + nbl
                  call MPI_Recv(dnj,nypj,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
               end if
c
               if (myid.eq.myhost) then
#endif
               write(iunit,*) 'zone t = "  ", i = ',k2-k1+1,
     .                        ', j = ',1,', k = ',1
               i      = 1
               write(iunit,*) 'npts = ',k2-k1+1
               write(iunit,401)
               do 7004 k=k1,k2
               write(iunit,402) real(x(jd,k,i)),real(z(jd,k,i)),
     .                          real(yplusj(k,i,m)),
     .                          real(dnj(k,i,m)),real(vistj(k,i,m))
 7004          continue
#if defined DIST_MPI
               end if
#endif
            end if
c
         end if
c
 7001    continue
c
         if (iflag.gt.0) then
c
#if defined DIST_MPI
            if (myid .ne. myhost) then
c
#endif
            ypavg = ypsum/npts
            do 7003 i=i1,i2
            do 7003 k=k1,k2
            diff    = yplusj(k,i,m) - ypavg
            ypsumsq = ypsumsq + diff*diff
 7003       continue
            ypsdev  = sqrt(ypsumsq/(npts-1))
c
#if defined DIST_MPI
            end if
c
            if (myid.eq.mblk2nd(nbl)) then
               intyp(1) = kypmax
               intyp(2) = iypmax
               intyp(3) = kypmin
               intyp(4) = iypmin
               intyp(5) = kdnmax
               intyp(6) = idnmax
               intyp(7) = kdnmin
               intyp(8) = idnmin
               intyp(9) = nypchk
               intyp(10) = npts
               relyp(1) = ypmax
               relyp(2) = ypmin
               relyp(3) = dnmax
               relyp(4) = dnmin
               relyp(5) = ypavg
               relyp(6) = ypsdev
               mytag = itag_int + nn + nbl
               call MPI_Send(intyp,10,MPI_INTEGER,myhost,
     .                       mytag,mycomm,ierr)
               mytag = itag_real + nn + nbl
               call MPI_Send(relyp,6,MY_MPI_REAL,myhost,
     .                       mytag,mycomm,ierr)
            end if
            if (myid.eq.myhost) then
               nd_srce = mblk2nd(nbl)
               mytag = itag_int + nn + nbl
               call MPI_Recv(intyp,10,MPI_INTEGER,nd_srce,
     .                       mytag,mycomm,istat,ierr)
               mytag = itag_real + nn + nbl
               call MPI_Recv(relyp,6,MY_MPI_REAL,nd_srce,
     .                       mytag,mycomm,istat,ierr)
               kypmax = intyp(1)
               iypmax = intyp(2)
               kypmin = intyp(3)
               iypmin = intyp(4)
               kdnmax = intyp(5)
               idnmax = intyp(6)
               kdnmin = intyp(7)
               idnmin = intyp(8)
               nypchk = intyp(9)
               npts   = intyp(10)
               ypmax  = relyp(1)
               ypmin  = relyp(2)
               dnmax  = relyp(3)
               dnmin  = relyp(4)
               ypavg  = relyp(5)
               ypsdev = relyp(6)
            end if
c
            if (myid.eq.myhost) then
#endif
            if (m.eq.1) then
               write(11,201)
            else
               write(11,202)
            end if
            write(11,203)
            write(11,106) real(ypmax),kypmax,iypmax,
     .                    real(ypmin),kypmin,iypmin
            write(11,204)
            write(11,106) real(dnmax),kdnmax,idnmax,
     .                    real(dnmin),kdnmin,idnmin
            write(11,105) int(real(ypchk))
            write(11,107) real(ypavg),real(ypsdev),nypchk,npts
#if defined DIST_MPI
            end if
#endif
c
            if (ifunc.gt.0 .and. i2d.eq.0) then
#if defined DIST_MPI
c
               nxyz = idim*jdim*kdim
               nypk = jdim*idim*2
               nypj = kdim*idim*2
               nypi = jdim*kdim*2           
               if (myid.eq.mblk2nd(nbl)) then
                  mytag = itag_x + nbl
                  call MPI_Send(x,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_y + nbl
                  call MPI_Send(y,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_z + nbl
                  call MPI_Send(z,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_blank + nbl
                  call MPI_Send(blankj,nypj,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_vist + nbl
                  call MPI_Send(vistj,nypj,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_yp + nbl
                  call MPI_Send(yplusj,nypj,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_dn + nbl
                  call MPI_Send(dnj,nypj,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
               end if
               if (myid.eq.myhost) then
                  nd_srce = mblk2nd(nbl)
                  mytag = itag_x + nbl
                  call MPI_Recv(x,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_y + nbl
                  call MPI_Recv(y,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_z + nbl
                  call MPI_Recv(z,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_blank + nbl
                  call MPI_Recv(blankj,nypj,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_vist + nbl
                  call MPI_Recv(vistj,nypj,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_yp + nbl
                  call MPI_Recv(yplusj,nypj,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_dn + nbl
                  call MPI_Recv(dnj,nypj,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
               end if
c
               if (myid.eq.myhost) then
#endif
               jj = 1
               if (m.eq.2) jj = jdim
               write(iunit,*) ((real(x(jj,k,i)),i=1,idim),k=1,kdim),
     .                        ((real(y(jj,k,i)),i=1,idim),k=1,kdim),
     .                        ((real(z(jj,k,i)),i=1,idim),k=1,kdim),
     .               ((int(real(blankj(k,i,m))),i=1,idim),k=1,kdim)
               write(iunit+1,*) ((real(yplusj(k,i,m)),i=1,idim),
     .                            k=1,kdim)
               write(iunit+2,*) ((real(dnj(k,i,m)),i=1,idim),
     .                            k=1,kdim)
               write(iunit+3,*) ((real(vistj(k,i,m)),i=1,idim),
     .                            k=1,kdim)
#if defined DIST_MPI
               end if
#endif
            end if
#if defined DIST_MPI
c
            if (myid .ne. myhost) then
#endif
c
c           block totals
c
            ypsumb   = ypsumb   + ypsum
            ypsumsqb = ypsumsqb + ypsumsq
            nptsb    = nptsb    + npts
            nypchkb  = nypchkb  + nypchk
            if (real(ypmax) .gt. real(ypmaxb)) then
               ypmaxb   = ypmax
               iypmaxb  = iypmax
               kypmaxb  = kypmax
               if (m.eq.1) then
                  jypmaxb = 1
               else
                  jypmaxb = jdim
               end if
            end if
            if (real(ypmin) .lt. real(ypminb)) then
               ypminb   = ypmin
               iypminb  = iypmin
               kypminb  = kypmin
               if (m.eq.1) then
                  jypminb = 1
               else
                  jypminb = jdim
               end if
            end if
            if (real(dnmax) .gt. real(dnmaxb)) then
               dnmaxb   = dnmax
               idnmaxb  = idnmax
               kdnmaxb  = kdnmax
               if (m.eq.1) then
                  jdnmaxb = 1
               else
                  jdnmaxb = jdim
               end if
            end if
            if (real(dnmin) .lt. real(dnminb)) then
               dnminb   = dnmin
               idnminb  = idnmin
               kdnminb  = kdnmin
               if (m.eq.1) then
                  jdnminb = 1
               else
                  jdnminb = jdim
               end if
            end if
#if defined DIST_MPI
c
            end if
c
#endif
         end if
c
 7000    continue
c
      end if
c
      if (idim.eq.2) go to 999
c
      if (ivisc(1).gt.1) then
c
c        i=1 and/or i=idim surfaces
c
         do 8000 m=1,2
c
         iflag   = 0
         ypsum   = 0.
         ypsumsq = 0.
         ypmax   = 0.
         jypmax  = 0
         kypmax  = 0
         ypmin   = 1.e9
         jypmin  = 0
         kypmin  = 0
         dnmax   = 0.
         jdnmax  = 0
         kdnmax  = 0
         dnmin   = 1.e9
         jdnmin  = 0
         kdnmin  = 0
         nypchk  = 0
         npts    = 0
c
         if (m.eq.1) then
            mm  = 2
            nbc = nbci0(nbl)
            id  = 1
            idx = 1
            id1 = id+1
         else
            mm  = 4
            nbc = nbcidim(nbl)
            id  = idim1
            idx = idim
            id1 = id-1
         end if
c
         do 7990 k=1,kdim
         do 7990 j=1,jdim
         blanki(j,k,m) = 0.
         yplusi(j,k,m) = 0.
         dni(j,k,m)    = 0.
         visti(j,k,m)  = 0.
 7990    continue
c
         do 8001 nseg = 1,nbc
c
         if (abs(ibcinfo(nbl,nseg,1,m)) .eq. 2004 .or.
     .       abs(ibcinfo(nbl,nseg,1,m)) .eq. 2014 .or.
     .       abs(ibcinfo(nbl,nseg,1,m)) .eq. 2024 .or.
     .       abs(ibcinfo(nbl,nseg,1,m)) .eq. 2034 .or.
     .       abs(ibcinfo(nbl,nseg,1,m)) .eq. 2016) then
c
            iflag = 1
c
c   do not include endpoints of segment, to avoid nonphysical behavior
            j1 = ibcinfo(nbl,nseg,2,m)+1
            j2 = ibcinfo(nbl,nseg,3,m)-1
            k1 = ibcinfo(nbl,nseg,4,m)+1
            k2 = ibcinfo(nbl,nseg,5,m)-1
            if (j1 .gt. j2) j1=j2
            if (k1 .gt. k2) k1=k2
c
#if defined DIST_MPI
            if (myid .ne. myhost) then
c
#endif
            do 8002 j=j1,j2
            jd   = j
            jd1  = jd-1
            if (jd1.lt.1)    jd1 = 1
            if (jd.gt.jdim1) jd = jdim1
            do 8002 k=k1,k2
            kd   = k
            kd1  = kd-1
            if (kd1.lt.1)    kd1 = 1
            if (kd.gt.kdim1) kd = kdim1
c
            blanki(j,k,m) = 1.
c
c           wall value - average in J and k
c
            q1i1 = .25*(qi0(jd,kd,1,mm)  + qi0(jd1,kd,1,mm)
     .           +      qi0(jd,kd1,1,mm) + qi0(jd1,kd1,1,mm))
            q2i1 = .25*(qi0(jd,kd,2,mm)  + qi0(jd1,kd,2,mm)
     .           +      qi0(jd,kd1,2,mm) + qi0(jd1,kd1,2,mm))
            q3i1 = .25*(qi0(jd,kd,3,mm)  + qi0(jd1,kd,3,mm)
     .           +      qi0(jd,kd1,3,mm) + qi0(jd1,kd1,3,mm))
            q4i1 = .25*(qi0(jd,kd,4,mm)  + qi0(jd1,kd,4,mm)
     .           +      qi0(jd,kd1,4,mm) + qi0(jd1,kd1,4,mm))
            q5i1 = .25*(qi0(jd,kd,5,mm)  + qi0(jd1,kd,5,mm)
     .           +      qi0(jd,kd1,5,mm) + qi0(jd1,kd1,5,mm))*gamma
            t1i1 = q5i1/q1i1
c
c           first cell center location - average in J and K
c
            q1i2 = .25*(q(jd,kd,id,1)   + q(jd,kd1,id,1)
     .           +      q(jd1,kd,id,1)  + q(jd1,kd1,id,1))
            q2i2 = .25*(q(jd,kd,id,2)   + q(jd,kd1,id,2)
     .           +      q(jd1,kd,id,2)  + q(jd1,kd1,id,2))
            q3i2 = .25*(q(jd,kd,id,3)   + q(jd,kd1,id,3)
     .           +      q(jd1,kd,id,3)  + q(jd1,kd1,id,3))
            q4i2 = .25*(q(jd,kd,id,4)   + q(jd,kd1,id,4)
     .           +      q(jd1,kd,id,4)  + q(jd1,kd1,id,4))
            q5i2 = .25*(q(jd,kd,id,5)   + q(jd,kd1,id,5)
     .           +      q(jd1,kd,id,5)  + q(jd1,kd1,id,5))*gamma
            vist = .25*(vist3d(jd,kd,id)  + vist3d(jd,kd1,id)
     .           +      vist3d(jd1,kd,id) + vist3d(jd1,kd1,id))
            t1i2 = q5i2/q1i2
c
            if (si(jd ,kd ,idx,4) .eq. 0. .or.
     +          si(jd ,kd1,idx,4) .eq. 0. .or.
     +          si(jd1,kd ,idx,4) .eq. 0. .or.
     +          si(jd1,kd1,idx,4) .eq. 0.) then
            dx    = x(j,k,id+1) - x(j,k,id)
            dy    = y(j,k,id+1) - y(j,k,id)
            dz    = z(j,k,id+1) - z(j,k,id)
            dn    = sqrt(dx**2 + dy**2 + dz**2)
            else
            dn = (vol(jd ,kd ,id)/si(jd ,kd ,idx,4)+
     +            vol(jd ,kd1,id)/si(jd ,kd1,idx,4)+
     +            vol(jd1,kd ,id)/si(jd1,kd ,idx,4)+
     +            vol(jd1,kd1,id)/si(jd1,kd1,idx,4))/4.
            end if
            if (real(dn).gt.real(dnmax)) then
               dnmax  = dn
               jdnmax = j
               kdnmax = k
            end if
            if (real(dn).lt.real(dnmin)) then
               dnmin  = dn
               jdnmin = j
               kdnmin = k
            end if
            dni(j,k,m)   = dn
            visti(j,k,m) = vist
            if (real(dn).eq.0.) then
               dn = 1.
            end if
            qa1   = sqrt(q2i1**2 + q3i1**2 + q4i1**2)
            qa2   = sqrt(q2i2**2 + q3i2**2 + q4i2**2)
c        Get turb viscosity at wall (0 unless wall fn used)
            if (iwf(1) .eq. 1) then
              avgmut=.25*(vi0(jd,kd,1,mm)  +vi0(jd1,kd,1,mm)
     .                  + vi0(jd,kd1,1,mm) +vi0(jd1,kd1,1,mm))
            else
             avgmut=0.
            end if
            emuka = (t1i1**1.5)*((1.0+cbar/tinf)/(t1i1+cbar/tinf))
            pres1 = q5i1
            cfy   = 2.*(qa2-qa1)/dn
            cfy   = 2.*(emuka+avgmut)*cfy/(reue*xmach)
c
            yplusi(j,k,m) = 0.
            if (abs(real(cfy)).gt.0.) yplusi(j,k,m) = 
     .         dn*reue*q1i1*sqrt(abs(real(cfy)*0.5/real(q1i1)))/emuka
            ypsum   = ypsum   + yplusi(j,k,m)
            if (real(yplusi(j,k,m)).gt.real(ypmax)) then
               ypmax  = yplusi(j,k,m)
               jypmax = j
               kypmax = k
            end if
            if (real(yplusi(j,k,m)).gt.0. .and.
     .         real(yplusi(j,k,m)).lt.real(ypmin)) then
               ypmin  = yplusi(j,k,m)
               jypmin = j
               kypmin = k
            end if
            if (real(yplusi(j,k,m)).gt.real(ypchk)) then
               nypchk = nypchk+1
            end if
c
 8002       continue
c
            npts  = npts + (k2-k1+1)*(j2-j1+1)
#if defined DIST_MPI
c
            end if
#endif
c
         end if
c
 8001    continue
c
         if (iflag.gt.0) then
c
#if defined DIST_MPI
            if (myid .ne. myhost) then
c
#endif
            ypavg = ypsum/npts
            do 8003 j=j1,j2
            do 8003 k=k1,k2
            diff    = yplusi(j,k,m) - ypavg
            ypsumsq = ypsumsq + diff*diff
 8003       continue
            ypsdev  = sqrt(ypsumsq/(npts-1))
c
#if defined DIST_MPI
            end if
c
            if (myid.eq.mblk2nd(nbl)) then
c
               intyp(1) = jypmax
               intyp(2) = kypmax
               intyp(3) = jypmin
               intyp(4) = kypmin
               intyp(5) = jdnmax
               intyp(6) = kdnmax
               intyp(7) = jdnmin
               intyp(8) = kdnmin
               intyp(9) = nypchk
               intyp(10) = npts
               relyp(1) = ypmax
               relyp(2) = ypmin
               relyp(3) = dnmax
               relyp(4) = dnmin
               relyp(5) = ypavg
               relyp(6) = ypsdev
               mytag = itag_int + 2*nn + nbl
               call MPI_Send(intyp,10,MPI_INTEGER,myhost,
     .                       mytag,mycomm,ierr)
               mytag = itag_real + 2*nn + nbl
               call MPI_Send(relyp,6,MY_MPI_REAL,myhost,
     .                       mytag,mycomm,ierr)
            end if
            if (myid.eq.myhost) then
               nd_srce = mblk2nd(nbl)
               mytag = itag_int + 2*nn + nbl
               call MPI_Recv(intyp,10,MPI_INTEGER,nd_srce,
     .                       mytag,mycomm,istat,ierr)
               mytag = itag_real + 2*nn + nbl
               call MPI_Recv(relyp,6,MY_MPI_REAL,nd_srce,
     .                       mytag,mycomm,istat,ierr)
               jypmax = intyp(1)
               kypmax = intyp(2)
               jypmin = intyp(3)
               kypmin = intyp(4)
               jdnmax = intyp(5)
               kdnmax = intyp(6)
               jdnmin = intyp(7)
               kdnmin = intyp(8)
               nypchk = intyp(9)
               npts   = intyp(10)
               ypmax  = relyp(1)
               ypmin  = relyp(2)
               dnmax  = relyp(3)
               dnmin  = relyp(4)
               ypavg  = relyp(5)
               ypsdev = relyp(6)
            end if
c
            if (myid.eq.myhost) then
c
#endif
            if (m.eq.1) then
               write(11,301)
            else
               write(11,302)
            end if
            write(11,303)
            write(11,106) real(ypmax),jypmax,kypmax,
     .                    real(ypmin),jypmin,kypmin
            write(11,304)
            write(11,106) real(dnmax),jdnmax,kdnmax,
     .                    real(dnmin),jdnmin,kdnmin
            write(11,105) int(ypchk)
            write(11,107) real(ypavg),real(ypsdev),nypchk,npts
#if defined DIST_MPI
c
            end if
#endif
c
            if (ifunc.gt.0) then
#if defined DIST_MPI
c
               nxyz = idim*jdim*kdim
               nypk = jdim*idim*2
               nypj = kdim*idim*2
               nypi = jdim*kdim*2           
               if (myid.eq.mblk2nd(nbl)) then
                  mytag = itag_x + nbl
                  call MPI_Send(x,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_y + nbl
                  call MPI_Send(y,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_z + nbl
                  call MPI_Send(z,nxyz,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_blank + nbl
                  call MPI_Send(blanki,nypi,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_vist + nbl
                  call MPI_Send(visti,nypi,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_yp + nbl
                  call MPI_Send(yplusi,nypi,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
                  mytag = itag_dn + nbl
                  call MPI_Send(dni,nypi,MY_MPI_REAL,myhost,
     .                          mytag,mycomm,ierr)
               end if
               if (myid.eq.myhost) then
                  nd_srce = mblk2nd(nbl)
                  mytag = itag_x + nbl
                  call MPI_Recv(x,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_y + nbl
                  call MPI_Recv(y,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_z + nbl
                  call MPI_Recv(z,nxyz,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_blank + nbl
                  call MPI_Recv(blanki,nypi,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_vist + nbl
                  call MPI_Recv(visti,nypi,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_yp + nbl
                  call MPI_Recv(yplusi,nypi,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  mytag = itag_dn + nbl
                  call MPI_Recv(dni,nypi,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
               end if
c
               if (myid.eq.myhost) then
#endif
               ii = 1
               if (m.eq.2) ii = idim
               write(iunit,*) ((real(x(j,k,ii)),j=1,jdim),k=1,kdim),
     .                        ((real(y(j,k,ii)),j=1,jdim),k=1,kdim),
     .                        ((real(z(j,k,ii)),j=1,jdim),k=1,kdim),
     .               ((int(real(blanki(j,k,m))),j=1,jdim),k=1,kdim)
               write(iunit+1,*) ((real(yplusi(j,k,m)),j=1,jdim),
     .                            k=1,kdim)
               write(iunit+2,*) ((real(dni(j,k,m)),j=1,jdim),
     .                            k=1,kdim)
               write(iunit+3,*) ((real(visti(j,k,m)),j=1,jdim),
     .                            k=1,kdim)
#if defined DIST_MPI
               end if
#endif
            end if
#if defined DIST_MPI
c
            if (myid .ne. myhost) then
#endif
c
c           block totals
c
            ypsumb   = ypsumb   + ypsum
            ypsumsqb = ypsumsqb + ypsumsq
            nptsb    = nptsb    + npts
            nypchkb  = nypchkb  + nypchk
            if (real(ypmax) .gt. real(ypmaxb)) then
               ypmaxb   = ypmax
               kypmaxb  = kypmax
               jypmaxb  = jypmax
               if (m.eq.1) then
                  iypmaxb = 1
               else
                  iypmaxb = idim
               end if
            end if
            if (real(ypmin) .lt. real(ypminb)) then
               ypminb   = ypmin
               kypminb  = kypmin
               jypminb  = jypmin
               if (m.eq.1) then
                  iypminb = 1
               else
                  iypminb = idim
               end if
            end if
            if (real(dnmax) .gt. real(dnmaxb)) then
               dnmaxb   = dnmax
               kdnmaxb  = kdnmax
               jdnmaxb  = jdnmax
               if (m.eq.1) then
                  idnmaxb = 1
               else
                  idnmaxb = idim
               end if
            end if
            if (real(dnmin) .lt. real(dnminb)) then
               dnminb   = dnmin
               kdnminb  = kdnmin
               jdnminb  = jdnmin
               if (m.eq.1) then
                  idnminb = 1
               else
                  idnminb = idim
               end if
            end if
#if defined DIST_MPI
c
            end if
c
#endif
         end if
c
 8000    continue
c
      end if 
c
  101 format(/2x,'K=1    SURFACE:')
  102 format(/2x,'K=KDIM SURFACE:')
  103 format(5x,6hY+ MAX,3x,4hJLOC,3x,4hILOC,
     .       7x,6hY+ MIN,3x,4hJLOC,3x,4hILOC)
  104 format(5x,6hDN MAX,3x,4hJLOC,3x,4hILOC,
     .       7x,6hDN MIN,3x,4hJLOC,3x,4hILOC)
  105 format(5x,6hY+ AVG,4x,10hY+ STD DEV,6x,5hNY+ >,i2,
     .       3x,4hNPTS)
  106 format(1x,e10.3,4x,i3,4x,i3,3x,e10.3,4x,i3,4x,i3)
  107 format(1x,e10.3,4x,e10.3,7x,i6,1x,i6)
  201 format(/2x,'J=1    SURFACE:')
  202 format(/2x,'J=JDIM SURFACE:')
  203 format(5x,6hY+ MAX,3x,4hKLOC,3x,4hILOC,
     .       7x,6hY+ MIN,3x,4hKLOC,3x,4hILOC)
  204 format(5x,6hDN MAX,3x,4hKLOC,3x,4hILOC,
     .       7x,6hDN MIN,3x,4hKLOC,3x,4hILOC)
  301 format(/2x,'I=1    SURFACE:')
  302 format(/2x,'I=IDIM SURFACE:')
  303 format(5x,6hY+ MAX,3x,4hJLOC,3x,4hKLOC,
     .       7x,6hY+ MIN,3x,4hJLOC,3x,4hKLOC)
  304 format(5x,6hDN MAX,3x,4hJLOC,3x,4hKLOC,
     .       7x,6hDN MIN,3x,4hJLOC,3x,4hKLOC)
  401       format(7x,1hX,12x,1hY,12x,2hY+,11x,2hDN,6x,11hTURB. VISC.)
  402       format(5(2x,e11.5))
c
  999 continue
#if defined DIST_MPI
c
c     for parallel applications, pass the common block, YPINFO,
c     in real and integer parts for accumulating the global yplus
c     statistics in the calling routine, YPLUSOUT
c
      if (myid.eq.mblk2nd(nbl)) then
         mytag = itag_ypcomr + nbl
         call MPI_Send(ypsumb,7,MY_MPI_REAL,myhost,
     .                 mytag,mycomm,ierr)
         mytag = itag_ypcomi + nbl
         call MPI_Send(nptsb,14,MPI_INTEGER,myhost,
     .                 mytag,mycomm,ierr)
      end if
      if (myid.eq.myhost) then
         nd_srce = mblk2nd(nbl)
         mytag = itag_ypcomr + nbl
         call MPI_Recv(ypsumb,7,MY_MPI_REAL,nd_srce,
     .                 mytag,mycomm,istat,ierr)
         mytag = itag_ypcomi + nbl
         call MPI_Recv(nptsb,14,MPI_INTEGER,nd_srce,
     .                 mytag,mycomm,istat,ierr)
      end if
c
#endif
      return
      end
