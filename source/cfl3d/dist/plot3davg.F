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
      subroutine plot3davg(jdim,kdim,idim,                               
     .                  i1,i2,i3,j1,j2,j3,k1,k2,k3,q,                    
     .                  qi0,qj0,qk0,x,y,z,xw,blank2,blank,xg,iflag,
     .                  vist3d,iover,nblk,nmap,bcj,bck,bci,
     .                  vj0,vk0,vi0,ifunc,iplot,jdw,kdw,idw,
     .                  nplots,jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                  jeg,keg,ieg,ninter,iindex,intmax,nsub1,
     .                  maxxe,nblkk,nbli,limblk,isva,nblon,mxbli,
     .                  thetay,maxbl,maxgr,myid,myhost,mycomm,
     .                  mblk2nd,inpl3d,nblock,nblkpt,xv,
     .                  sj,sk,si,vol,nset,qavg,q2avg,nt,movabs)          
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Write the output file at the grid points in PLOT3D
c     format and print solution data.
c
c     outputs grid/solution in single precision for use with FAST/PLOT3D
c***********************************************************************
c
#if defined ADP_OFF
#   ifdef CMPLX
#     ifdef DBLE_PRECSN
      implicit complex*8(a-h,o-z)
#     else
      implicit complex(a-h,o-z)
#     endif
#   else
#     ifdef DBLE_PRECSN
      implicit real*8 (a-h,o-z)
#     endif
#   endif
#else
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
#endif
c
#if defined DIST_MPI
#     include "mpif.h"
#   ifdef P3D_SINGLE
#     define MY_MPI_REAL MPI_REAL
#   else
#     define MY_MPI_REAL MPI_DOUBLE_PRECISION
#   endif
#endif
#if defined P3D_SINGLE
      real*4    xw(jdw,kdw,idw,nset),xg(jdw,kdw,idw,4)
      real*4    xmachw,alphww,reuew,timew
#else
      real xw(jdw,kdw,idw,nset),xg(jdw,kdw,idw,4)
#endif

      dimension xv(jdw,kdw,idw,5)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension vist3d(jdim,kdim,idim)
      dimension vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4)
      dimension blank(jdim,kdim,idim),blank2(jdim,kdim,idim,2)
      dimension nmap(maxbl),jdimg(maxbl),kdimg(maxbl),idimg(maxbl),
     .          nblcg(maxbl),jsg(maxbl),ksg(maxbl),isg(maxbl),
     .          jeg(maxbl),keg(maxbl),ieg(maxbl),
     .          iindex(intmax,6*nsub1+9),nblkpt(maxxe),nblkk(2,mxbli),
     .          limblk(2,6,mxbli),isva(2,2,mxbli),nblon(mxbli), 
     .          inpl3d(nplots,11),thetay(maxbl),mblk2nd(maxbl)
      dimension sk(jdim,kdim,idim-1,5),sj(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
      dimension qavg(jdim,kdim,idim,5)                                   
      dimension q2avg(jdim,kdim,idim,5)                                  
c
      common /bin/ ibin,iblnk,iblnkfr,ip3dgrad
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /reyue/ reue,tinf,ivisc(3)
      common /twod/ i2d
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /wallfun/ iwf(3)
      common /igrdtyp/ ip3dgrd,ialph
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
      common /conversion/ radtodeg
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset = maxbl
      itag_grid = 1
      itag_q    = itag_grid + ioffset
#endif
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
c     initialize xw, xv and xg arrays
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
      do j=1,jw
         do k=1,kw
            do i=1,iw
               do l= 1,5
               xw(j,k,i,l) = 0.
               xv(j,k,i,l) = 0.
               end do
               do l= 1,4
               xg(j,k,i,l) = 0.
               end do
            end do
         end do
      end do
c
c     assign single precision scalars
c
      alphaw = radtodeg*(alpha+thetay(nblk))
      xmachw = xmach
      alphww = alphaw
      reuew  = reue
      timew  = time
c
c     determine q values at grid points and load into xv array
c
      jdw = (j2-j1)/j3 + 1
      kdw = (k2-k1)/k3 + 1
      idw = (i2-i1)/i3 + 1
      ldw = 5
c
      call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,q,
     .            qi0,qj0,qk0,jdw,kdw,idw,xv,ldw) 
c
c   get iteration-averaged Q values                                
c   note: the qavg values are kept as primitive variables          


c     ! maxnum contains the max vals of each of the variables
      maxnum = 0
      minnum = 99999999

c     ! maxpos contains the indices of the maximum vals (4 is block #)

      maxpos = 0
      minpos = 0

      do i=1,idim                                                  
        do j=1,jdim                                                
          do k=1,kdim                                              
            qavg(j,k,i,1)=(qavg(j,k,i,1)*(xnumavg-1.)+             
     +        xv(j,k,i,1))/xnumavg                                 
            qavg(j,k,i,2)=(qavg(j,k,i,2)*(xnumavg-1.)+             
     +        xv(j,k,i,2))/xnumavg                                 
            qavg(j,k,i,3)=(qavg(j,k,i,3)*(xnumavg-1.)+             
     +        xv(j,k,i,3))/xnumavg                                 
            qavg(j,k,i,4)=(qavg(j,k,i,4)*(xnumavg-1.)+             
     +        xv(j,k,i,4))/xnumavg                                 
            qavg(j,k,i,5)=(qavg(j,k,i,5)*(xnumavg-1.)+             
     +        xv(j,k,i,5))/xnumavg                                 

            q2avg(j,k,i,1)=(q2avg(j,k,i,1)*(xnumavg2-1.)+          
     +        xv(j,k,i,1)**2)/xnumavg2                             
            q2avg(j,k,i,2)=(q2avg(j,k,i,2)*(xnumavg2-1.)+          
     +        xv(j,k,i,2)**2)/xnumavg2                             
            q2avg(j,k,i,3)=(q2avg(j,k,i,3)*(xnumavg2-1.)+          
     +        xv(j,k,i,3)**2)/xnumavg2                             
            q2avg(j,k,i,4)=(q2avg(j,k,i,4)*(xnumavg2-1.)+          
     +        xv(j,k,i,4)**2)/xnumavg2                             
            q2avg(j,k,i,5)=(q2avg(j,k,i,5)*(xnumavg2-1.)+          
     +        xv(j,k,i,5)**2)/xnumavg2                             

          enddo                                                    
        enddo                                                      
      enddo                                                        

      return
      end
