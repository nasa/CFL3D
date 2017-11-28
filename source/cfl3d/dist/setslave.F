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
      subroutine setslave(lw,lw2,w,mgwk,wk,nwork,maxbl,mxbli,maxgr,
     .                    maxseg,nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                    nbckdim,ibcinfo,jbcinfo,kbcinfo,nblock,
     .                    idefrm,iadvance,nou,bou,nbuf,ibufdim,myid,
     .                    myhost,mycomm,mblk2nd,icsi,icsf,jcsi,
     .                    jcsf,kcsi,kcsf,islavept,nslave,nsegdfrm,
     .                    idfrmseg,maxsegdg,iwk,iwork,nmaster,ngrid,
     .                    jskip,kskip,iskip,nblg,levelg,lfgm,nblk,
     .                    limblk,isva,nblelst,nnodes,iskmax,jskmax,
     .                    kskmax,nbli) 
c
c     $Id$
c
c***********************************************************************
c     Purpose: Set up slave point data array for deforming meshes
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
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension w(mgwk),lw(65,maxbl),lw2(43,maxbl),wk(nwork),iwk(iwork)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),
     .          nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          ibcinfo(maxbl,maxseg,7,2),jbcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2)
      dimension idefrm(maxbl),nblg(maxgr)
      dimension iadvance(maxbl),mblk2nd(maxbl)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg)
      dimension iskip(maxbl,500),jskip(maxbl,500),kskip(maxbl,500)
      dimension islavept(nslave,nmaster,5),levelg(maxbl)
      dimension nblk(2,mxbli),limblk(2,6,mxbli),
     .          isva(2,2,mxbli)
      dimension nsgst1(2),nsgst(0:maxbl,2)
      dimension nskp1(3),nskp(3),in(3),im(3),ig(3),id(3),
     .          iss(3,2),ise(3,2),ijktot(2,3)
      dimension n11i(20,nslave),iimax(nslave),iimx1(nslave)
      dimension nblelst(maxbl,2)
      dimension iskmax(maxbl)
      dimension jskmax(maxbl)
      dimension kskmax(maxbl)
c
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /sklton/ isklton
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /elastic_ss/ idef_ss
      common /twod/ i2d
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt 
c
c
c
      islavept = 0
      if (iunst.gt.1 .or. idef_ss.gt.0) then
       do iseq=1,mseq
         do n = 1,nslave
           iimax(n)    = 1
         enddo
c
      lglobal = lfgm-(mseq-iseq)
c
c
c
         do nbl = 0,nblock
          do nn = 1,2
            nsgst(nbl,nn) = 0
          enddo
         enddo
         ivert     = 0
         nbl2      = 1
         do nbl=1,nblock
         if ((levelg(nbl).ge.lglobal .and.
     .      levelg(nbl).le.levelt(iseq))) then
c           if (myid.eq.mblk2nd(nbl)) then
               call lead(nbl,lw,lw2,maxbl)
               iskp = iskip(nbl,1)
               jskp = jskip(nbl,1)
               kskp = kskip(nbl,1)
               if(abs(isktyp).eq.1) then
                 iskmax(nbl) = (idim-1)/iskp + 1
                 jskmax(nbl) = (jdim-1)/jskp + 1
                 kskmax(nbl) = (kdim-1)/kskp + 1
                 nsgst(nbl2,1) = nsgst(nbl2-1,1) + 
     .              (jdim+jskp-1)*(kdim+kskp-1)*(idim+iskp-1)
     .                 /jskp/kskp/iskp
                 nsgst(nbl2,2) = nbl
                 nbl2  = nbl2 + 1
                 do i=1,idim,iskp 
                   do k=1,kdim,kskp
                     do j=1,jdim,jskp
                       ll = j-1+jdim*(k-1)+jdim*kdim*(i-1) 
                       ivert = ivert + 1
                       if (ivert.gt.nslave) then
                         nou(1) = min(nou(1)+1,ibufdim)
                         write(bou(nou(1),1),
     .                       '(''stopping...parameter nslave '',
     .                              ''too small'')')
                         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                       end if
                       islavept(ivert,1,iseq) = ll  
                       islavept(ivert,8,iseq) = 1
                       islavept(ivert,9,iseq) = nbl 
                       jm1 =-1
                       km1 =-(jdim+jskp-1)/jskp
                       im1 =-(jdim+jskp-1)*(kdim+kskp-1)/jskp/kskp
                       jp1 = 1
                       kp1 = (jdim+jskp-1)/jskp
                       ip1 = (jdim+jskp-1)*(kdim+kskp-1)/jskp/kskp
                       if(j.eq.   1) jm1 = 0 
                       if(j.eq.jdim) jp1 = 0 
                       if(k.eq.   1) km1 = 0 
                       if(k.eq.kdim) kp1 = 0 
                       if(i.eq.   1) im1 = 0 
                       if(i.eq.idim) ip1 = 0 
                       if(i2d.ne.0) goto 1000
                         do nseg = 1,nbci0(nbl)
                          ista = 1
                          iend = 1
                          jsta = ibcinfo(nbl,nseg,2,1)
                          jend = ibcinfo(nbl,nseg,3,1)
                          ksta = ibcinfo(nbl,nseg,4,1)
                          kend = ibcinfo(nbl,nseg,5,1)
                          if (ibcinfo(nbl,nseg,1,1).ne.2005) then     
                           if(ibcinfo(nbl,nseg,1,1).eq.1005.or. 
     .                        ibcinfo(nbl,nseg,1,1).eq.1006.or. 
     .                        ibcinfo(nbl,nseg,1,1).eq.2004.or.   
     .                        ibcinfo(nbl,nseg,1,1).eq.2014.or.   
     .                        ibcinfo(nbl,nseg,1,1).eq.2024.or.   
     .                        ibcinfo(nbl,nseg,1,1).eq.2034.or.   
     .                        ibcinfo(nbl,nseg,1,1).eq.2016.or.   
     .                        ibcinfo(nbl,nseg,1,1).eq.1002.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.1000.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.1003.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.1008.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2003.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2006.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2007.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2008.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2102.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2103) then
                            if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                        and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                         then 
                              islavept(ivert,8,iseq)=0
                            end if
                           end if
                          end if
                         enddo
                         do nseg = 1,nbcidim(nbl)
                          ista = idim
                          iend = idim
                          jsta = ibcinfo(nbl,nseg,2,2)
                          jend = ibcinfo(nbl,nseg,3,2)
                          ksta = ibcinfo(nbl,nseg,4,2)
                          kend = ibcinfo(nbl,nseg,5,2)
                          if (ibcinfo(nbl,nseg,1,2).ne.2005) then
                           if(ibcinfo(nbl,nseg,1,2).eq.1005.or. 
     .                        ibcinfo(nbl,nseg,1,2).eq.1006.or. 
     .                        ibcinfo(nbl,nseg,1,2).eq.2004.or.   
     .                        ibcinfo(nbl,nseg,1,2).eq.2014.or.   
     .                        ibcinfo(nbl,nseg,1,2).eq.2024.or.   
     .                        ibcinfo(nbl,nseg,1,2).eq.2034.or.   
     .                        ibcinfo(nbl,nseg,1,2).eq.2016.or.   
     .                        ibcinfo(nbl,nseg,1,2).eq.1002.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.1000.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.1003.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.1008.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2003.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2006.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2007.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2008.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2102.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2103) then
                            if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                        and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                         then 
                              islavept(ivert,8,iseq)=0
                            end if
                           end if
                          end if
                         enddo
1000                   continue
                       do nseg = 1,nbcj0(nbl)
                        ista = jbcinfo(nbl,nseg,2,1)
                        iend = jbcinfo(nbl,nseg,3,1)
                        jsta = 1
                        jend = 1
                        ksta = jbcinfo(nbl,nseg,4,1)
                        kend = jbcinfo(nbl,nseg,5,1)
                        if (jbcinfo(nbl,nseg,1,1).ne.2005) then
                         if(jbcinfo(nbl,nseg,1,1).eq.1005.or. 
     .                      jbcinfo(nbl,nseg,1,1).eq.1006.or. 
     .                      jbcinfo(nbl,nseg,1,1).eq.2004.or.   
     .                      jbcinfo(nbl,nseg,1,1).eq.2014.or.   
     .                      jbcinfo(nbl,nseg,1,1).eq.2024.or.   
     .                      jbcinfo(nbl,nseg,1,1).eq.2034.or.   
     .                      jbcinfo(nbl,nseg,1,1).eq.2016.or.   
     .                      jbcinfo(nbl,nseg,1,1).eq.1002.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.1000.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.1003.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.1008.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2003.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2006.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2007.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2008.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2102.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2103.or.
     .                     (jbcinfo(nbl,nseg,1,1).eq.1013.and.   
     .                      idim.eq.2)) then
                          if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                       and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                       then 
                            islavept(ivert,8,iseq)=0
                          end if
                         end if
                        end if
                       enddo
                       do nseg = 1,nbcjdim(nbl)
                        ista = jbcinfo(nbl,nseg,2,2)
                        iend = jbcinfo(nbl,nseg,3,2)
                        jsta = jdim
                        jend = jdim
                        ksta = jbcinfo(nbl,nseg,4,2)
                        kend = jbcinfo(nbl,nseg,5,2)
                        if (jbcinfo(nbl,nseg,1,2).ne.2005) then
                         if(jbcinfo(nbl,nseg,1,2).eq.1005.or. 
     .                      jbcinfo(nbl,nseg,1,2).eq.1006.or. 
     .                      jbcinfo(nbl,nseg,1,2).eq.2004.or.   
     .                      jbcinfo(nbl,nseg,1,2).eq.2014.or.   
     .                      jbcinfo(nbl,nseg,1,2).eq.2024.or.   
     .                      jbcinfo(nbl,nseg,1,2).eq.2034.or.   
     .                      jbcinfo(nbl,nseg,1,2).eq.2016.or.   
     .                      jbcinfo(nbl,nseg,1,2).eq.1002.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.1000.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.1003.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.1008.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2003.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2006.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2007.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2008.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2102.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2103.or.
     .                     (jbcinfo(nbl,nseg,1,2).eq.1013.and.   
     .                      idim.eq.2)) then
                          if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                       and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                       then 
                            islavept(ivert,8,iseq)=0
                          end if
                         end if
                        end if
                       enddo
                       do nseg = 1,nbck0(nbl)
                        ista = kbcinfo(nbl,nseg,2,1)
                        iend = kbcinfo(nbl,nseg,3,1)
                        jsta = kbcinfo(nbl,nseg,4,1)
                        jend = kbcinfo(nbl,nseg,5,1)
                        ksta = 1
                        kend = 1
                        if (kbcinfo(nbl,nseg,1,1).ne.2005) then
                         if(kbcinfo(nbl,nseg,1,1).eq.1005.or. 
     .                      kbcinfo(nbl,nseg,1,1).eq.1006.or. 
     .                      kbcinfo(nbl,nseg,1,1).eq.2004.or.   
     .                      kbcinfo(nbl,nseg,1,1).eq.2014.or.   
     .                      kbcinfo(nbl,nseg,1,1).eq.2024.or.   
     .                      kbcinfo(nbl,nseg,1,1).eq.2034.or.   
     .                      kbcinfo(nbl,nseg,1,1).eq.2016.or.   
     .                      kbcinfo(nbl,nseg,1,1).eq.1002.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.1000.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.1003.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.1008.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2003.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2006.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2007.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2008.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2102.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2103.or.
     .                     (kbcinfo(nbl,nseg,1,1).eq.1013.and.   
     .                      idim.eq.2)) then
                          if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                       and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                       then 
                            islavept(ivert,8,iseq)=0
                          end if
                         end if
                        end if
                       enddo
                       do nseg = 1,nbckdim(nbl)
                        ista = kbcinfo(nbl,nseg,2,2)
                        iend = kbcinfo(nbl,nseg,3,2)
                        jsta = kbcinfo(nbl,nseg,4,2)
                        jend = kbcinfo(nbl,nseg,5,2)
                        ksta = kdim
                        kend = kdim
                        if (kbcinfo(nbl,nseg,1,2).ne.2005) then
                         if(kbcinfo(nbl,nseg,1,2).eq.1005.or. 
     .                      kbcinfo(nbl,nseg,1,2).eq.1006.or. 
     .                      kbcinfo(nbl,nseg,1,2).eq.2004.or.   
     .                      kbcinfo(nbl,nseg,1,2).eq.2014.or.   
     .                      kbcinfo(nbl,nseg,1,2).eq.2024.or.   
     .                      kbcinfo(nbl,nseg,1,2).eq.2034.or.   
     .                      kbcinfo(nbl,nseg,1,2).eq.2016.or.   
     .                      kbcinfo(nbl,nseg,1,2).eq.1002.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.1000.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.1003.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.1008.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2003.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2006.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2007.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2008.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2102.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2103.or.
     .                     (kbcinfo(nbl,nseg,1,2).eq.1013.and.   
     .                      idim.eq.2)) then
                          if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                       and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                       then 
                            islavept(ivert,8,iseq)=0
                          end if
                         end if
                        end if
                       enddo
                       islavept(ivert,2,iseq) = ivert+jm1  
                       islavept(ivert,3,iseq) = ivert+jp1  
                       islavept(ivert,4,iseq) = ivert+km1  
                       islavept(ivert,5,iseq) = ivert+kp1  
                       islavept(ivert,6,iseq) = ivert+im1  
                       islavept(ivert,7,iseq) = ivert+ip1  
                     end do
                   end do
                 end do
               else
                 do n = 2,500
                   if(iskip(nbl,n).eq.0) then
                      iskmax(nbl) = n - 1
                      goto 10
                   end if
                 enddo
10               continue
                 do n = 2,500
                   if(jskip(nbl,n).eq.0) then
                      jskmax(nbl) = n - 1
                      goto 20
                   end if
                 enddo
20               continue
                 do n = 2,500
                   if(kskip(nbl,n).eq.0) then
                      kskmax(nbl) = n - 1
                      goto 30
                   end if
                 enddo
30               continue
                 itot = iskmax(nbl) 
                 jtot = jskmax(nbl) 
                 ktot = kskmax(nbl) 
                 nsgst(nbl2,1) = nsgst(nbl2-1,1) + itot*jtot*ktot 
                 nsgst(nbl2,2) = nbl
                 nbl2  = nbl2 + 1
                 do i1=1,itot
                   i = iskip(nbl,i1)
                   do k1=1,ktot
                     k = kskip(nbl,k1)
                     do j1=1,jtot
                       j = jskip(nbl,j1)
                       ll = j-1+jdim*(k-1)+jdim*kdim*(i-1) 
                       iwrap = 0
                       jwrap = 0
                       kwrap = 0
                       ivert = ivert + 1
                       if (ivert.gt.nslave) then
                         nou(1) = min(nou(1)+1,ibufdim)
                         write(bou(nou(1),1),
     .                       '(''stopping...parameter nslave '',
     .                              ''too small'')')
                         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                       end if
                       islavept(ivert,1,iseq) = ll  
                       islavept(ivert,8,iseq) = 1
                       islavept(ivert,9,iseq) = nbl 
                       jm1 =-1
                       km1 =-jtot
                       im1 =-jtot*ktot
                       jp1 = 1
                       kp1 = jtot
                       ip1 = jtot*ktot
                       if(j.eq.   1) jm1 = 0 
                       if(j.eq.jdim) jp1 = 0 
                       if(k.eq.   1) km1 = 0 
                       if(k.eq.kdim) kp1 = 0 
                       if(i.eq.   1) im1 = 0 
                       if(i.eq.idim) ip1 = 0 
                       if(i2d.ne.0) goto 2000
                         do nseg = 1,nbci0(nbl)
                          ista = 1
                          iend = 1
                          jsta = ibcinfo(nbl,nseg,2,1)
                          jend = ibcinfo(nbl,nseg,3,1)
                          ksta = ibcinfo(nbl,nseg,4,1)
                          kend = ibcinfo(nbl,nseg,5,1)
                          if (ibcinfo(nbl,nseg,1,1).ne.2005) then     
                           if(ibcinfo(nbl,nseg,1,1).eq.1005.or. 
     .                        ibcinfo(nbl,nseg,1,1).eq.1006.or. 
     .                        ibcinfo(nbl,nseg,1,1).eq.2004.or.   
     .                        ibcinfo(nbl,nseg,1,1).eq.2014.or.   
     .                        ibcinfo(nbl,nseg,1,1).eq.2024.or.   
     .                        ibcinfo(nbl,nseg,1,1).eq.2034.or.   
     .                        ibcinfo(nbl,nseg,1,1).eq.2016.or.   
     .                        ibcinfo(nbl,nseg,1,1).eq.1002.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.1000.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.1003.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.1008.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2003.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2006.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2007.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2008.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2102.or.
     .                        ibcinfo(nbl,nseg,1,1).eq.2103) then
                            if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                        and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                         then 
                              islavept(ivert,8,iseq)=0
                            end if
                           end if
                          end if
                         enddo
                         do nseg = 1,nbcidim(nbl)
                          ista = idim
                          iend = idim
                          jsta = ibcinfo(nbl,nseg,2,2)
                          jend = ibcinfo(nbl,nseg,3,2)
                          ksta = ibcinfo(nbl,nseg,4,2)
                          kend = ibcinfo(nbl,nseg,5,2)
                          if (ibcinfo(nbl,nseg,1,2).ne.2005) then
                           if(ibcinfo(nbl,nseg,1,2).eq.1005.or. 
     .                        ibcinfo(nbl,nseg,1,2).eq.1006.or. 
     .                        ibcinfo(nbl,nseg,1,2).eq.2004.or.   
     .                        ibcinfo(nbl,nseg,1,2).eq.2014.or.   
     .                        ibcinfo(nbl,nseg,1,2).eq.2024.or.   
     .                        ibcinfo(nbl,nseg,1,2).eq.2034.or.   
     .                        ibcinfo(nbl,nseg,1,2).eq.2016.or.   
     .                        ibcinfo(nbl,nseg,1,2).eq.1002.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.1000.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.1003.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.1008.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2003.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2006.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2007.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2008.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2102.or.
     .                        ibcinfo(nbl,nseg,1,2).eq.2103) then
                            if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                        and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                         then 
                              islavept(ivert,8,iseq)=0
                            end if
                           end if
                          end if
                         enddo
2000                   continue
                       do nseg = 1,nbcj0(nbl)
                        ista = jbcinfo(nbl,nseg,2,1)
                        iend = jbcinfo(nbl,nseg,3,1)
                        jsta = 1
                        jend = 1
                        ksta = jbcinfo(nbl,nseg,4,1)
                        kend = jbcinfo(nbl,nseg,5,1)
                        if (jbcinfo(nbl,nseg,1,1).ne.2005) then
                         if(jbcinfo(nbl,nseg,1,1).eq.1005.or. 
     .                      jbcinfo(nbl,nseg,1,1).eq.1006.or. 
     .                      jbcinfo(nbl,nseg,1,1).eq.2004.or.   
     .                      jbcinfo(nbl,nseg,1,1).eq.2014.or.   
     .                      jbcinfo(nbl,nseg,1,1).eq.2024.or.   
     .                      jbcinfo(nbl,nseg,1,1).eq.2034.or.   
     .                      jbcinfo(nbl,nseg,1,1).eq.2016.or.   
     .                      jbcinfo(nbl,nseg,1,1).eq.1002.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.1000.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.1003.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.1008.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2003.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2006.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2007.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2008.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2102.or.
     .                      jbcinfo(nbl,nseg,1,1).eq.2103.or.
     .                     (jbcinfo(nbl,nseg,1,1).eq.1013.and.   
     .                      idim.eq.2)) then
                          if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                       and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                       then 
                            islavept(ivert,8,iseq)=0
                          end if
                         end if
                        end if
                       enddo
                       do nseg = 1,nbcjdim(nbl)
                        ista = jbcinfo(nbl,nseg,2,2)
                        iend = jbcinfo(nbl,nseg,3,2)
                        jsta = jdim
                        jend = jdim
                        ksta = jbcinfo(nbl,nseg,4,2)
                        kend = jbcinfo(nbl,nseg,5,2)
                        if (jbcinfo(nbl,nseg,1,2).ne.2005) then
                         if(jbcinfo(nbl,nseg,1,2).eq.1005.or. 
     .                      jbcinfo(nbl,nseg,1,2).eq.1006.or. 
     .                      jbcinfo(nbl,nseg,1,2).eq.2004.or.   
     .                      jbcinfo(nbl,nseg,1,2).eq.2014.or.   
     .                      jbcinfo(nbl,nseg,1,2).eq.2024.or.   
     .                      jbcinfo(nbl,nseg,1,2).eq.2034.or.   
     .                      jbcinfo(nbl,nseg,1,2).eq.2016.or.   
     .                      jbcinfo(nbl,nseg,1,2).eq.1002.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.1000.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.1003.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.1008.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2003.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2006.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2007.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2008.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2102.or.
     .                      jbcinfo(nbl,nseg,1,2).eq.2103.or.
     .                     (jbcinfo(nbl,nseg,1,2).eq.1013.and.   
     .                      idim.eq.2)) then
                          if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                       and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                       then 
                            islavept(ivert,8,iseq)=0
                          end if
                         end if
                        end if
                       enddo
                       do nseg = 1,nbck0(nbl)
                        ista = kbcinfo(nbl,nseg,2,1)
                        iend = kbcinfo(nbl,nseg,3,1)
                        jsta = kbcinfo(nbl,nseg,4,1)
                        jend = kbcinfo(nbl,nseg,5,1)
                        ksta = 1
                        kend = 1
                        if (kbcinfo(nbl,nseg,1,1).ne.2005) then
                         if(kbcinfo(nbl,nseg,1,1).eq.1005.or. 
     .                      kbcinfo(nbl,nseg,1,1).eq.1006.or. 
     .                      kbcinfo(nbl,nseg,1,1).eq.2004.or.   
     .                      kbcinfo(nbl,nseg,1,1).eq.2014.or.   
     .                      kbcinfo(nbl,nseg,1,1).eq.2024.or.   
     .                      kbcinfo(nbl,nseg,1,1).eq.2034.or.   
     .                      kbcinfo(nbl,nseg,1,1).eq.2016.or.   
     .                      kbcinfo(nbl,nseg,1,1).eq.1002.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.1000.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.1003.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.1008.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2003.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2006.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2007.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2008.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2102.or.
     .                      kbcinfo(nbl,nseg,1,1).eq.2103.or.
     .                     (kbcinfo(nbl,nseg,1,1).eq.1013.and.   
     .                      idim.eq.2)) then
                          if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                       and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                       then 
                            islavept(ivert,8,iseq)=0
                          end if
                         end if
                        end if
                       enddo
                       do nseg = 1,nbckdim(nbl)
                        ista = kbcinfo(nbl,nseg,2,2)
                        iend = kbcinfo(nbl,nseg,3,2)
                        jsta = kbcinfo(nbl,nseg,4,2)
                        jend = kbcinfo(nbl,nseg,5,2)
                        ksta = kdim
                        kend = kdim
                        if (kbcinfo(nbl,nseg,1,2).ne.2005) then
                         if(kbcinfo(nbl,nseg,1,2).eq.1005.or. 
     .                      kbcinfo(nbl,nseg,1,2).eq.1006.or. 
     .                      kbcinfo(nbl,nseg,1,2).eq.2004.or.   
     .                      kbcinfo(nbl,nseg,1,2).eq.2014.or.   
     .                      kbcinfo(nbl,nseg,1,2).eq.2024.or.   
     .                      kbcinfo(nbl,nseg,1,2).eq.2034.or.   
     .                      kbcinfo(nbl,nseg,1,2).eq.2016.or.   
     .                      kbcinfo(nbl,nseg,1,2).eq.1002.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.1000.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.1003.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.1008.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2003.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2006.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2007.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2008.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2102.or.
     .                      kbcinfo(nbl,nseg,1,2).eq.2103.or.
     .                     (kbcinfo(nbl,nseg,1,2).eq.1013.and.   
     .                      idim.eq.2)) then
                          if(i.ge.ista.and.i.le.iend.and.j.ge.jsta.
     .                       and.j.le.jend.and.k.ge.ksta.and.k.le.kend)
     .                       then 
                            islavept(ivert,8,iseq)=0
                          end if
                         end if
                        end if
                       enddo
                       islavept(ivert,2,iseq) = ivert+jm1  
                       islavept(ivert,3,iseq) = ivert+jp1  
                       islavept(ivert,4,iseq) = ivert+km1  
                       islavept(ivert,5,iseq) = ivert+kp1  
                       islavept(ivert,6,iseq) = ivert+im1  
                       islavept(ivert,7,iseq) = ivert+ip1  
                     end do
                   end do
                 end do
               end if
c           end if
          end if
         end do 
         do nbl = nblock,1,-1
          nsgst(nbl,1) = nsgst(nbl-1,1)
         enddo
         ivert = 0
         do m = 1,nbli
          nbl1 = nblk(1,m)
          if ((levelg(nbl1).ge.lglobal .and.
     .         levelg(nbl1).le.levelt(iseq))) then
            do nbl = 1,nblock
              if(nsgst(nbl,2).eq.nblk(1,m)) then
                 nsgst12 = nsgst(nbl,2)
                 nsgst11 = nsgst(nbl,1)
              end if
              if(nsgst(nbl,2).eq.nblk(2,m)) then
                 nsgst22 = nsgst(nbl,2)
                 nsgst21 = nsgst(nbl,1)
              end if
            enddo
            call lead(nsgst12,lw,lw2,maxbl)
            idim1 = idim 
            jdim1 = jdim 
            kdim1 = kdim 
            call lead(nsgst22,lw,lw2,maxbl)
            idim2 = idim 
            jdim2 = jdim 
            kdim2 = kdim 
c
c  Use limblk(2...) data as boundary data for this block
c                
            do nn = 1,3
              iss(nn,1)= 0
              ise(nn,1)= 0
              iss(nn,2)= 0
              ise(nn,2)= 0
            enddo
            isva21 = isva(2,1,m)
            isva22 = isva(2,2,m)
            isva11 = isva(1,1,m)
            isva12 = isva(1,2,m)
            if(isva(2,1,m).eq.1.and.isva(2,2,m).eq.2) then
             isva23 = 3
            else if(isva(2,2,m).eq.1.and.isva(2,1,m).eq.2) then
             isva23 = 3
            else if(isva(2,1,m).eq.2.and.isva(2,2,m).eq.3) then
             isva23 = 1
            else if(isva(2,2,m).eq.2.and.isva(2,1,m).eq.3) then
             isva23 = 1
            else if(isva(2,1,m).eq.1.and.isva(2,2,m).eq.3) then
             isva23 = 2
            else if(isva(2,2,m).eq.1.and.isva(2,1,m).eq.3) then
             isva23 = 2
            end if
            if(isva(1,1,m).eq.1.and.isva(1,2,m).eq.2) then
             isva13 = 3
            else if(isva(1,2,m).eq.1.and.isva(1,1,m).eq.2) then
             isva13 = 3
            else if(isva(1,1,m).eq.2.and.isva(1,2,m).eq.3) then
             isva13 = 1
            else if(isva(1,2,m).eq.2.and.isva(1,1,m).eq.3) then
             isva13 = 1
            else if(isva(1,1,m).eq.1.and.isva(1,2,m).eq.3) then
             isva13 = 2
            else if(isva(1,2,m).eq.1.and.isva(1,1,m).eq.3) then
             isva13 = 2
            end if
            lmb2210= limblk(2,isva21  ,m) 
            lmb2213= limblk(2,isva21+3,m) 
            lmb2220= limblk(2,isva22  ,m) 
            lmb2223= limblk(2,isva22+3,m) 
            lmb2230= limblk(2,isva23  ,m) 
            lmb2233= limblk(2,isva23+3,m) 
            lmb1110= limblk(1,isva11  ,m) 
            lmb1113= limblk(1,isva11+3,m) 
            lmb1120= limblk(1,isva12  ,m) 
            lmb1123= limblk(1,isva12+3,m) 
            lmb1130= limblk(1,isva13  ,m) 
            lmb1133= limblk(1,isva13+3,m) 
            iss(isva23,2) = lmb2230    
            ise(isva23,2) = lmb2233 
            if(lmb2210.le.lmb2213) then
              if(lmb2220.le.lmb2223) then 
               iss(isva21,2) = lmb2210    
               ise(isva21,2) = lmb2213+1 
               iss(isva22,2) = lmb2220    
               ise(isva22,2) = lmb2223+1 
              else
               iss(isva21,2) = lmb2210    
               ise(isva21,2) = lmb2213+1 
               iss(isva22,2) = lmb2220+1    
               ise(isva22,2) = lmb2223 
              end if
            else
              if(lmb2220.le.lmb2223) then 
               iss(isva21,2) = lmb2210+1    
               ise(isva21,2) = lmb2213 
               iss(isva22,2) = lmb2220    
               ise(isva22,2) = lmb2223+1 
              else
               iss(isva21,2) = lmb2210+1    
               ise(isva21,2) = lmb2213 
               iss(isva22,2) = lmb2220+1    
               ise(isva22,2) = lmb2223 
              end if
            end if
c
c  Use limblk(1...) data as boundary data for this block
c                
            iss(isva13,1) = lmb1130    
            ise(isva13,1) = lmb1133 
            if(lmb1110.le.lmb1113) then 
              if(lmb1120.le.lmb1123) then 
               iss(isva11,1) = lmb1110    
               ise(isva11,1) = lmb1113+1 
               iss(isva12,1) = lmb1120    
               ise(isva12,1) = lmb1123+1 
              else
               iss(isva11,1) = lmb1110    
               ise(isva11,1) = lmb1113+1 
               iss(isva12,1) = lmb1120+1    
               ise(isva12,1) = lmb1123 
              end if
            else
              if(lmb1120.le.lmb1123) then 
               iss(isva11,1) = lmb1110+1    
               ise(isva11,1) = lmb1113 
               iss(isva12,1) = lmb1120    
               ise(isva12,1) = lmb1123+1 
              else
               iss(isva11,1) = lmb1110+1    
               ise(isva11,1) = lmb1113 
               iss(isva12,1) = lmb1120+1    
               ise(isva12,1) = lmb1123 
              end if
            end if
            if(abs(isktyp).eq.1) then 
              do ii = 1,3
                do jj = 1,2
                  if(iss(ii,jj).eq.0.or.ise(ii,jj).eq.0) then
                         nou(1) = min(nou(1)+1,ibufdim)
                         write(bou(nou(1),1),
     .          '(''stopping...deforming mesh control points do not '',
     .                              ''match blocking end points'')')
                         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                enddo
              enddo
              iskp0  = iskip(nsgst12,1)
              jskp0  = jskip(nsgst12,1) 
              kskp0  = kskip(nsgst12,1) 
              nskp(1)= iskip(nsgst12,1)
              nskp(2)= jskip(nsgst12,1) 
              nskp(3)= kskip(nsgst12,1) 
              iskp2  = iskip(nsgst22,1)
              jskp2  = jskip(nsgst22,1) 
              kskp2  = kskip(nsgst22,1) 
              nskp1(1)  = iskip(nsgst22,1)
              nskp1(2)  = jskip(nsgst22,1) 
              nskp1(3)  = kskip(nsgst22,1) 
              do nn = 1,3
               if(iss(nn,1).eq.ise(nn,1)) nskp(nn) = 0
               if(iss(nn,1).gt.ise(nn,1)) nskp(nn) = -nskp(nn) 
              enddo
              do nn = 1,3
               if(iss(nn,2).gt.ise(nn,2)) nskp1(nn) = -nskp1(nn) 
              enddo
              if     (iss(2,1).eq.ise(2,1).and.iss(2,1).eq.1    ) then
                inc1= 1
                jj1 = 2
              else if(iss(2,1).eq.ise(2,1).and.iss(2,1).eq.jdim1) then
                inc1=-1
                jj1 = 3
              else if(iss(3,1).eq.ise(3,1).and.iss(3,1).eq.1    ) then
                inc1= (jdim1+jskp0-1)/jskp0 
                jj1 = 4
              else if(iss(3,1).eq.ise(3,1).and.iss(3,1).eq.kdim1) then
                inc1=-(jdim1+jskp0-1)/jskp0 
                jj1 = 5
              else if(iss(1,1).eq.ise(1,1).and.iss(1,1).eq.1    ) then
                inc1= (jdim1+jskp0-1)*(kdim1+kskp0-1)/jskp0/kskp0 
                jj1 = 6
              else
                inc1=-(jdim1+jskp0-1)*(kdim1+kskp0-1)/jskp0/kskp0 
                jj1 = 7
              end if
              if     (iss(2,2).eq.ise(2,2).and.iss(2,2).eq.1    ) then
                inc2= 1
                jj2 = 2
              else if(iss(2,2).eq.ise(2,2).and.iss(2,2).eq.jdim2) then
                inc2=-1
                jj2 = 3
              else if(iss(3,2).eq.ise(3,2).and.iss(3,2).eq.1    ) then
                inc2= (jdim2+jskp2-1)/jskp2 
                jj2 = 4
              else if(iss(3,2).eq.ise(3,2).and.iss(3,2).eq.kdim2) then
                inc2=-(jdim2+jskp2-1)/jskp2 
                jj2 = 5
              else if(iss(1,2).eq.ise(1,2).and.iss(1,2).eq.1    ) then
                inc2= (jdim2+jskp2-1)*(kdim2+kskp2-1)/jskp2/kskp2 
                jj2 = 6
              else
                inc2=-(jdim2+jskp2-1)*(kdim2+kskp2-1)/jskp2/kskp2 
                jj2 = 7
              end if
              do nn = 1,2
               in(nn) = isva(1,nn,m)
               im(nn) = isva(2,nn,m)
              enddo
              in(3) = isva13
              im(3) = isva23
              id(in(1)) = iss(in(1),1)
              do i2= iss(im(1),2),ise(im(1),2),nskp1(im(1))
                ig(im(1)) = i2
                id(in(3)) = iss(in(3),1)
                do k2= iss(im(3),2),ise(im(3),2),nskp1(im(3)) 
                  ig(im(3)) = k2
                  id(in(2)) = iss(in(2),1)
                  do j2= iss(im(2),2),ise(im(2),2),nskp1(im(2)) 
                    ig(im(2)) = j2
                    ll1=(id(2)+jskp0-1)/jskp0 + 
     .                  (jdim1+jskp0-1)*(id(3)      -1)/jskp0/kskp0+  
     .                  (jdim1+jskp0-1)*(kdim1+kskp0-1)*
     .                  (id(1)      -1)/jskp0/kskp0/iskp0 
                    ll2=(ig(2)+jskp2-1)/jskp2 + 
     .                  (jdim2+jskp2-1)*(ig(3)      -1)/jskp2/kskp2+  
     .                  (jdim2+jskp2-1)*(kdim2+kskp2-1)*
     .                  (ig(1)      -1)/jskp2/kskp2/iskp2 
                    islavept(nsgst11+ll1,jj1,iseq)=nsgst21+ll2+inc2 
                    if(nsgst11+ll1.ne.nsgst21+ll2) then
                      iimax(nsgst11+ll1) = iimax(nsgst11+ll1) + 1
                      ii3 = 10 + iimax(nsgst11+ll1) 
                      islavept(nsgst11+ll1,ii3,iseq) = nsgst21+ll2 
                      iimax(nsgst21+ll2) = iimax(nsgst21+ll2) + 1
                      ii3 = 10 + iimax(nsgst21+ll2) 
                      islavept(nsgst21+ll2,ii3,iseq) = nsgst11+ll1 
                    end if
                    islavept(nsgst21+ll2,jj2,iseq)=nsgst11+ll1+inc1  
                    id(in(2)) = id(in(2)) + nskp(in(2))
                  enddo
                  id(in(3)) = id(in(3)) + nskp(in(3))
                enddo
                id(in(1)) = id(in(1)) + nskp(in(1))
              enddo    
            else
              do ii = 1,3
                do jj = 1,2
                  if(iss(ii,jj).eq.0.or.ise(ii,jj).eq.0) then
                         nou(1) = min(nou(1)+1,ibufdim)
                         write(bou(nou(1),1),
     .          '(''stopping...deforming mesh control points do not '',
     .                              ''match blocking end points'')')
                         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                enddo
              enddo
              do nn = 1,3
               nskp(nn) = 1
               if(iss(nn,1).eq.ise(nn,1)) nskp(nn) = 0
               if(iss(nn,1).gt.ise(nn,1)) nskp(nn) = -1
              enddo
              do nn = 1,3
               nskp1(nn) = 1
               if(iss(nn,2).gt.ise(nn,2)) nskp1(nn) = -1
              enddo
              ijktot(1,1)= iskmax(nsgst12)
              ijktot(1,2)= jskmax(nsgst12)
              ijktot(1,3)= kskmax(nsgst12)
              ijktot(2,1)= iskmax(nsgst22)
              ijktot(2,2)= jskmax(nsgst22)
              ijktot(2,3)= kskmax(nsgst22)
              if     (iss(2,1).eq.ise(2,1).and.iss(2,1).eq.1    ) then
                inc1= 1
                jj1 = 2
              else if(iss(2,1).eq.ise(2,1).and.iss(2,1).eq.jdim1) then
                inc1=-1
                jj1 = 3
              else if(iss(3,1).eq.ise(3,1).and.iss(3,1).eq.1    ) then
                inc1= ijktot(1,2) 
                jj1 = 4
              else if(iss(3,1).eq.ise(3,1).and.iss(3,1).eq.kdim1) then
                inc1=-ijktot(1,2) 
                jj1 = 5
              else if(iss(1,1).eq.ise(1,1).and.iss(1,1).eq.1    ) then
                inc1= ijktot(1,2)*ijktot(1,3)
                jj1 = 6
              else
                inc1=-ijktot(1,2)*ijktot(1,3)
                jj1 = 7
              end if
              if     (iss(2,2).eq.ise(2,2).and.iss(2,2).eq.1    ) then
                inc2= 1
                jj2 = 2
              else if(iss(2,2).eq.ise(2,2).and.iss(2,2).eq.jdim2) then
                inc2=-1
                jj2 = 3
              else if(iss(3,2).eq.ise(3,2).and.iss(3,2).eq.1    ) then
                inc2= ijktot(2,2)
                jj2 = 4
              else if(iss(3,2).eq.ise(3,2).and.iss(3,2).eq.kdim2) then
                inc2=-ijktot(2,2)
                jj2 = 5
              else if(iss(1,2).eq.ise(1,2).and.iss(1,2).eq.1    ) then
                inc2= ijktot(2,2)*ijktot(2,3)
                jj2 = 6
              else
                inc2=-ijktot(2,2)*ijktot(2,3)
                jj2 = 7
              end if
              do nn = 1,2
               in(nn) = isva(1,nn,m)
               im(nn) = isva(2,nn,m)
              enddo
              in(3) = isva13
              im(3) = isva23
              do ii = 1,ijktot(1,1)  
                 if(iss(1,1).eq.iskip(nsgst12,ii)) then
                   issim11= ii
                 end if
                 if(ise(1,1).eq.iskip(nsgst12,ii)) then
                   iseim11= ii
                 end if
              enddo
              do ii = 1,ijktot(2,1)  
                 if(iss(1,2).eq.iskip(nsgst22,ii)) then
                   issim12= ii
                 end if
                 if(ise(1,2).eq.iskip(nsgst22,ii)) then
                   iseim12= ii
                 end if
              enddo
              do ii = 1,ijktot(1,2)  
                 if(iss(2,1).eq.jskip(nsgst12,ii)) then
                   issim21= ii
                 end if
                 if(ise(2,1).eq.jskip(nsgst12,ii)) then
                   iseim21= ii
                 end if
              enddo
              do ii = 1,ijktot(2,2)  
                 if(iss(2,2).eq.jskip(nsgst22,ii)) then
                   issim22= ii
                 end if
                 if(ise(2,2).eq.jskip(nsgst22,ii)) then
                   iseim22= ii
                 end if
              enddo
              do ii = 1,ijktot(1,3)  
                 if(iss(3,1).eq.kskip(nsgst12,ii)) then
                   issim31= ii
                 end if
                 if(ise(3,1).eq.kskip(nsgst12,ii)) then
                   iseim31= ii
                 end if
              enddo
              do ii = 1,ijktot(2,3)  
                 if(iss(3,2).eq.kskip(nsgst22,ii)) then
                   issim32= ii
                 end if
                 if(ise(3,2).eq.kskip(nsgst22,ii)) then
                   iseim32= ii
                 end if
              enddo
              iss(1,1) = issim11 
              ise(1,1) = iseim11 
              iss(2,1) = issim21 
              ise(2,1) = iseim21
              iss(3,1) = issim31 
              ise(3,1) = iseim31 
              iss(1,2) = issim12 
              ise(1,2) = iseim12 
              iss(2,2) = issim22 
              ise(2,2) = iseim22
              iss(3,2) = issim32 
              ise(3,2) = iseim32 
              do ii = 1,3
                do jj = 1,2
                  if(iss(ii,jj).eq.0.or.ise(ii,jj).eq.0) then
                         nou(1) = min(nou(1)+1,ibufdim)
                         write(bou(nou(1),1),
     .          '(''stopping...deforming mesh control points do not '',
     .                              ''match blocking end points'')')
                         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                enddo
              enddo
              id(in(1)) = iss(in(1),1)
              do i2= iss(im(1),2),ise(im(1),2),nskp1(im(1))
                ig(im(1)) = i2
                id(in(3)) = iss(in(3),1)
                do k2= iss(im(3),2),ise(im(3),2),nskp1(im(3)) 
                  ig(im(3)) = k2
                  id(in(2)) = iss(in(2),1)
                  do j2= iss(im(2),2),ise(im(2),2),nskp1(im(2)) 
                    ig(im(2)) = j2
                    ll1=id(2)+(id(3)-1)*ijktot(1,2)+(id(1)-1)*
     .                                      ijktot(1,2)*ijktot(1,3) 
                    ll2=ig(2)+(ig(3)-1)*ijktot(2,2)+(ig(1)-1)*
     .                                      ijktot(2,2)*ijktot(2,3) 
                    islavept(nsgst11+ll1,jj1,iseq)=nsgst21+ll2+inc2 
                    if(nsgst11+ll1.ne.nsgst21+ll2) then
                      iimax(nsgst11+ll1) = iimax(nsgst11+ll1) + 1
                      ii3 = 10 + iimax(nsgst11+ll1) 
                      islavept(nsgst11+ll1,ii3,iseq) = nsgst21+ll2 
                      iimax(nsgst21+ll2) = iimax(nsgst21+ll2) + 1
                      ii3 = 10 + iimax(nsgst21+ll2) 
                      islavept(nsgst21+ll2,ii3,iseq) = nsgst11+ll1 
                    end if
                    islavept(nsgst21+ll2,jj2,iseq)=nsgst11+ll1+inc1  
                    id(in(2)) = id(in(2)) + nskp(in(2))
                  enddo
                  id(in(3)) = id(in(3)) + nskp(in(3))
                enddo
                id(in(1)) = id(in(1)) + nskp(in(1))
              enddo    
            end if
          end if
         enddo
c
c   Determine the number of coincident points at block interfaces, starting
c   with the first two coincident points.
c
         do n = 1,nslave
           if(iimax(n).eq.1) iimax(n)= 0
           islavept(n,11,iseq) = iimax(n)
           if(islavept(n,11,iseq).gt.0) then
             n11i(1,n)  = n
             do ii3 = 2,iimax(n) 
               n11i(ii3,n) = islavept(n,12+ii3-2,iseq)
             enddo
           end if
         enddo
         do n = 1,nslave
           if(islavept(n,11,iseq).gt.0) then
              iimx1(n) = iimax(n)
              do n2 = 1,nslave
c               if(n2.ne.n.and.iimax(n2).gt.1) then
                if(n2.ne.n) then
                 do ii3 = 1,iimax(n) 
                   do ii4 = 1,iimax(n2)-1 
                    if(n11i(ii3,n).eq.n11i(ii4,n2)) then
                      do ii5 = 2,iimax(n) 
                        if(n2.eq.n11i(ii5,n)) goto 2550
                      enddo
                      iimx1(n) = iimx1(n) + 1
                      n11i(iimx1(n),n) = n2
                      goto 2600
                    end if
2550               continue
                    if(n11i(ii3,n).eq.n2) then
                      do ii5 = 2,iimax(n) 
                        if(n11i(ii4,n2).eq.n11i(ii5,n)) goto 2600
                      enddo
                      iimx1(n) = iimx1(n)+ 1
                      n11i(iimx1(n),n) = n11i(ii4,n2)
                      goto 2600
                    end if
                   enddo
                 enddo
2600             continue
                 islavept(n,11,iseq) = iimx1(n) 
                 if(iimx1(n).gt.iimax(n).and.iimax(n).gt.0) then
                   ii4 = iimax(n)+1
                   do ii = 10+iimax(n)+1,10+iimx1(n) 
                     islavept(n,ii,iseq) = n11i(ii4,n)
                     ii4 = ii4 + 1
                   enddo
                 end if 
               end if
              enddo
2650          continue
           end if
         enddo
c
c   Establish the number of points in the ija and sa arrays for each
c   node point.
c
      do n = 1,nslave
       ipt1 = 0
       ipt2 = 0
       ii4  = 1
       if(islavept(n,8,iseq).eq.0) ii4 = 0
       if(ii4.ne.0) then
         iimax(n) = islavept(n,11,iseq)
         do ii2 = 2,iimax(n) 
           ni3 = islavept(n,12+ii2-2,iseq)
           if(islavept(ni3,8,iseq).eq.0) ii4 = 0
         enddo
       end if 
       if(ii4.ne.0) then
         ipt2 = 2
         if(islavept(n,2,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,3,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,4,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,5,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,6,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,7,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,2,iseq).ne.n.and.islavept(n,4,iseq).ne.n) ipt1= 
     .      ipt1+ 1
         if(islavept(n,2,iseq).ne.n.and.islavept(n,6,iseq).ne.n) ipt1=
     .      ipt1+ 1
         if(islavept(n,4,iseq).ne.n.and.islavept(n,6,iseq).ne.n) ipt1=
     .      ipt1+ 1
         if(islavept(n,2,iseq).ne.n.and.islavept(n,4,iseq).ne.n.and.
     .      islavept(n,6,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,2,iseq).ne.n.and.islavept(n,7,iseq).ne.n) ipt1=
     .      ipt1+ 1
         if(islavept(n,4,iseq).ne.n.and.islavept(n,7,iseq).ne.n) ipt1=
     .      ipt1+ 1
         if(islavept(n,2,iseq).ne.n.and.islavept(n,4,iseq).ne.n.and.
     .      islavept(n,7,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,3,iseq).ne.n.and.islavept(n,4,iseq).ne.n) ipt1=
     .      ipt1+ 1
         if(islavept(n,3,iseq).ne.n.and.islavept(n,6,iseq).ne.n) ipt1=
     .      ipt1+ 1
         if(islavept(n,3,iseq).ne.n.and.islavept(n,4,iseq).ne.n.and.
     .      islavept(n,6,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,3,iseq).ne.n.and.islavept(n,7,iseq).ne.n) ipt1=
     .      ipt1+ 1
         if(islavept(n,3,iseq).ne.n.and.islavept(n,4,iseq).ne.n.and.
     .      islavept(n,7,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,2,iseq).ne.n.and.islavept(n,5,iseq).ne.n) ipt1=
     .      ipt1+ 1
         if(islavept(n,5,iseq).ne.n.and.islavept(n,6,iseq).ne.n) ipt1=
     .      ipt1+ 1
         if(islavept(n,2,iseq).ne.n.and.islavept(n,5,iseq).ne.n.and.
     .      islavept(n,6,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,5,iseq).ne.n.and.islavept(n,7,iseq).ne.n) ipt1=
     .      ipt1+ 1
         if(islavept(n,2,iseq).ne.n.and.islavept(n,5,iseq).ne.n.and.
     .      islavept(n,7,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,3,iseq).ne.n.and.islavept(n,5,iseq).ne.n) ipt1=
     .      ipt1+ 1
         if(islavept(n,3,iseq).ne.n.and.islavept(n,5,iseq).ne.n.and.
     .      islavept(n,6,iseq).ne.n) ipt1= ipt1+ 1
         if(islavept(n,3,iseq).ne.n.and.islavept(n,5,iseq).ne.n.and.
     .      islavept(n,7,iseq).ne.n) ipt1= ipt1+ 1
       end if
         islavept(n,10,iseq) =  3*ipt1 + ipt2
      enddo
c
        nblelst      = 0
        ii           = 1
        nblelst(ii,1) = islavept(1,9,iseq)
        nblelst(ii,2) = 1
        do n = 2,nslave
          if(islavept(n,9,iseq).ne.islavept(n-1,9,iseq)) then
           ii            = ii + 1
           nblelst(ii,1) = islavept(n,9,iseq)
           nblelst(ii,2) = n
          end if
        enddo
        nblelst(ii+1,2) = nslave+1 
       end do
      end if
      if(myid.eq.myhost) then
        if(isktyp.gt.0.and.meshdef.eq.1) then
          write(1000+myid,31029)
31029     format(' Iteration   Total Resid   N@Maxres  Max Resid')
        else if(isktyp.gt.0.and.meshdef.eq.0) then
          write(1000+myid,31028)
31028     format(' Iterations  Total Resid')
        end if
      end if
c
      return
      end
