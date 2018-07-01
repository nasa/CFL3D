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
c=============================================================================80
c
c  All the following is for Method of Manufactured Solution (MMS)
c  currently coded in CFL3D according to the Lisbon Workshop series
c  (see ,e.g., AIAA 2007-4089 or Int J. Numer. Meth. fluids 54:119-154, 2007)
c  Currently only coded for use with Splart-Allmaras turbulence model
c
c  Currently coded to run "MS1"
c  To run, the grid should be a 2-D square grid 0.5<=x<=1, 0<=z<=0.5
c  (stretched in the z-direction, with lower wall clustering, is best)
c
c  The following is a typical cfl3d.inp file (with "c " in 1st columns):
c
c I/O FILES
c mms_stretch1_289_1153.p3d
c plot3dg.bin
c plot3dq.bin
c cfl3d.out
c cfl3d.res
c cfl3d.turres
c cfl3d.blomax
c cfl3d.out15
c cfl3d.prout
c cfl3d.out20
c ovrlp.bin
c patch.bin
c restart.bin
c >
c ifullns 1
c iturbord 2
c iexact_trunc 0
c iexact_disc 1
c iexact_ring 0
c <
c     turbulent flat plate (plate from j=17-65, prior to 17 is symmetry)
c      XMACH     ALPHA      BETA  REUE,MIL   TINF,DR     IALPH    IHSTRY
c     0.2000    00.000       0.0    01.000     540.0         0         0
c       SREF      CREF      BREF       XMC       YMC       ZMC
c    1.00000   1.00000    1.0000   0.00000      0.00      0.00
c         DT     IREST   IFLAGTS      FMAX     IUNST    CFLTAU
c     -5.000         1       000   05.0000         0       10.
c      NGRID   NPLOT3D    NPRINT    NWREST      ICHK       I2D    NTSTEP       ITA
c         -1         1         0      5000         0         1         1         1
c        NCG       IEM  IADVANCE    IFORCE  IVISC(I)  IVISC(J)  IVISC(K)
c          3         0         0       001         0        05        05
c     IDIM    JDIM    KDIM
c       02     289    1153
c     ILAMLO    ILAMHI    JLAMLO    JLAMHI    KLAMLO    KLAMHI
c          0         0         0         0         0         0
c      INEWG    IGRIDC        IS        JS        KS        IE        JE        KE
c          0         0         0         0         0         0         0         0
c   IDIAG(I)  IDIAG(J)  IDIAG(K)  IFLIM(I)  IFLIM(J)  IFLIM(K)
c          1         1         1         0         0         0
c    IFDS(I)   IFDS(J)   IFDS(K)  RKAP0(I)  RKAP0(J)  RKAP0(K)
c          1         1         1   0.33333   0.33333   0.33333
c       GRID     NBCI0   NBCIDIM     NBCJ0   NBCJDIM     NBCK0   NBCKDIM    IOVRLP
c          1         1         1         1         1         1         1         0
c I0:   GRID   SEGMENT    BCTYPE      JSTA      JEND      KSTA      KEND     NDATA
c          1         1      1001         0         0         0         0         0
c IDIM: GRID   SEGMENT    BCTYPE      JSTA      JEND      KSTA      KEND     NDATA
c          1         1      1002         0         0         0         0         0
c J0:   GRID   SEGMENT    BCTYPE      ISTA      IEND      KSTA      KEND     NDATA
c          1         1      9999         0         0         0         0         0
c JDIM: GRID   SEGMENT    BCTYPE      ISTA      IEND      KSTA      KEND     NDATA
c          1         1      9999         0         0         0         0         0
c K0:   GRID   SEGMENT    BCTYPE      ISTA      IEND      JSTA      JEND     NDATA
c          1         1      2004         0         0         0         0         2
c               TWTYPE        CQ
c                  -1.        0.
c KDIM: GRID   SEGMENT    BCTYPE      ISTA      IEND      JSTA      JEND     NDATA
c          1         1      9999         0         0         0         0         0
c       MSEQ    MGFLAG    ICONSF       MTT      NGAM
c          2         1         0         0        02
c       ISSC EPSSSC(1) EPSSSC(2) EPSSSC(3)      ISSR EPSSSR(1) EPSSSR(2) EPSSSR(3)
c          0       0.3       0.3       0.3         0       0.3       0.3      0.3
c       NCYC    MGLEVG     NEMGL     NITFO
c      00001        03        00       000
c      40000        03        00       000
c       MIT1      MIT2      MIT3      MIT4      MIT5      MIT6      MIT7     MIT8
c         01        01        01        01        01         1         1        1
c         01        01        01        01        01         1         1        1
c    1-1 BLOCKING DATA:
c       NBLI
c          0
c  NUMBER   GRID     :    ISTA   JSTA   KSTA   IEND   JEND   KEND  ISVA1  ISVA2
c  NUMBER   GRID     :    ISTA   JSTA   KSTA   IEND   JEND   KEND  ISVA1  ISVA2
c   PATCH SURFACE DATA:
c     NINTER
c          0
c   PLOT3D OUTPUT:
c    GRID IPTYPE ISTART   IEND   IINC JSTART   JEND   JINC KSTART   KEND   KINC
c       1      0      0      0      0      0      0      0      0      0      0
c  IMOVIE
c       0
c   PRINT OUT:
c    GRID IPTYPE ISTART   IEND   IINC JSTART   JEND   JINC KSTART   KEND   KINC
c   CONTROL SURFACE:
c   NCS
c     0
c    GRID ISTART   IEND   JSTART   JEND   KSTART   KEND  IWALL  INORM
c
c  To run a test on discretization error, use the above input file,
c  and converge the solution.
c  Results are given in fort.201, fort.301, and fort.302
c  To run a test on truncation error only, change:
c     iexact_trunc from 0 to 1
c     iexact_disc  from 1 to 0
c     NCG=0, MSEQ=1, MGFLAG=0, NCYC=1, MGLEVG=1
c     start from scratch (IREST=0)
c
c=============================================================================80
      subroutine exact_norm(jdim,kdim,idim,x,y,z,q,tursav,vol,res,smin,
     +   iexact_trunc,iexact_disc)
c==================================== EXACT_NORM =============================80
c
c Routine to get and print out norms
c
c=============================================================================80
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /reyue/ reue,tinf,ivisc(3)
c
      dimension x(jdim,kdim,idim), y(jdim,kdim,idim), z(jdim,kdim,idim),
     +          q(jdim,kdim,idim,5), tursav(jdim,kdim,idim)
      dimension vol(jdim,kdim,idim-1), res(jdim,kdim,idim-1,5)
      dimension smin(jdim-1,kdim-1,idim-1)
c   ***
c     dimension xtemp(jdim-1,kdim-1),ztemp(jdim-1,kdim-1)
c   ***
c
      dimension f(6)
c
c  delta Q's
      u_max=0.
      w_max=0.
      p_max=0.
      t_max=0.
      i_u_max=0
      j_u_max=0
      k_u_max=0
      i_w_max=0
      j_w_max=0
      k_w_max=0
      i_p_max=0
      j_p_max=0
      k_p_max=0
      i_t_max=0
      j_t_max=0
      k_t_max=0
      sum_u=0.
      sum_w=0.
      sum_p=0.
      sum_t=0.
c
      i_convert_q = 0
      i_forcing = 0
      i_gradient = 0
      nvar=6
      do i=1,idim-1
        do j=1,jdim-1
          do k=1,kdim-1
            do n=1,nvar
              f(n)=0.
            enddo
            distf=smin(j,k,i)
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            delta_u=ccabs(q(j,k,i,2)-f(2))/xmach
            delta_w=ccabs(q(j,k,i,4)-f(4))/xmach
            delta_p=ccabs(q(j,k,i,5)-f(5))*gamma/(xmach*xmach)
            delta_t=ccabs(tursav(j,k,i)-f(6))/reue
            if (real(delta_u) .gt. real(u_max)) then
              u_max=delta_u
              i_u_max=i
              j_u_max=j
              k_u_max=k
            end if
            if (real(delta_w) .gt. real(w_max)) then
              w_max=delta_w
              i_w_max=i
              j_w_max=j
              k_w_max=k
            end if
            if (real(delta_p) .gt. real(p_max)) then
              p_max=delta_p
              i_p_max=i
              j_p_max=j
              k_p_max=k
            end if
            if (real(delta_t) .gt. real(t_max)) then
              t_max=delta_t
              i_t_max=i
              j_t_max=j
              k_t_max=k
            end if
            sum_u=sum_u+delta_u
            sum_w=sum_w+delta_w
            sum_p=sum_p+delta_p
            sum_t=sum_t+delta_t
          enddo
        enddo
      enddo
      sum_u=sum_u/((jdim-1)*(kdim-1)*(idim-1))
      sum_w=sum_w/((jdim-1)*(kdim-1)*(idim-1))
      sum_p=sum_p/((jdim-1)*(kdim-1)*(idim-1))
      sum_t=sum_t/((jdim-1)*(kdim-1)*(idim-1))
c
      write(201,'(''# discretization error (compared to exact soln)'')')
      write(201,'(''variables="N","1/N","sqrt(1/N)",'',
     +            ''"u_max","i_u_max","j_u_max","k_u_max","sum_u",'')')
      write(201,'(''"w_max","i_w_max","j_w_max","k_w_max","sum_w",'',
     +            ''"p_max","i_p_max","j_p_max","k_p_max","sum_p",'')')
      write(201,'(''"t_max","i_t_max","j_t_max","k_t_max","sum_t"'')')
      inumtot=(jdim-1)*(kdim-1)
      hhhh=1./(float((jdim-1)*(kdim-1)))
      sqhhhh=sqrt(hhhh)
      write(201,'(i9,3e18.8,3i8,e18.8,e18.8,3i8,e18.8,e18.8,3i8,e18.8,
     +  e18.8,3i5,e18.8)') inumtot,hhhh,sqhhhh,
     +                     u_max,i_u_max,j_u_max,k_u_max,sum_u,
     +                     w_max,i_w_max,j_w_max,k_w_max,sum_w,
     +                     p_max,i_p_max,j_p_max,k_p_max,sum_p,
     +                     t_max,i_t_max,j_t_max,k_t_max,sum_t
c
c   residuals
      res_1_max=0.
      res_2_max=0.
      res_3_max=0.
      res_4_max=0.
      res_5_max=0.
      sum_1_res=0.
      sum_2_res=0.
      sum_3_res=0.
      sum_4_res=0.
      sum_5_res=0.
      i_res_3_max=0
      j_res_3_max=0
      k_res_3_max=0
      i=1
      do j=3,jdim-3
        do k=3,kdim-3
          if (real(ccabs(res(j,k,i,1)/vol(j,k,i))) .gt. 
     +       real(res_1_max)) then
            res_1_max = ccabs(res(j,k,i,1)/vol(j,k,i))
            i_res_1_max=i
            j_res_1_max=j
            k_res_1_max=k
            res_1_max_alone=res(j,k,i,1)
            vol_1_max_alone=vol(j,k,i)
          end if
          if (real(ccabs(res(j,k,i,2)/vol(j,k,i))) .gt. 
     +       real(res_2_max)) then
            res_2_max = ccabs(res(j,k,i,2)/vol(j,k,i))
            i_res_2_max=i
            j_res_2_max=j
            k_res_2_max=k
            res_2_max_alone=res(j,k,i,2)
            vol_2_max_alone=vol(j,k,i)
          end if
          if (real(ccabs(res(j,k,i,3)/vol(j,k,i))) .gt. 
     +       real(res_3_max)) then
            res_3_max = ccabs(res(j,k,i,3)/vol(j,k,i))
            i_res_3_max=i
            j_res_3_max=j
            k_res_3_max=k
            res_3_max_alone=res(j,k,i,3)
            vol_3_max_alone=vol(j,k,i)
          end if
          if (real(ccabs(res(j,k,i,4)/vol(j,k,i))) .gt. 
     +       real(res_4_max)) then
            res_4_max = ccabs(res(j,k,i,4)/vol(j,k,i))
            i_res_4_max=i
            j_res_4_max=j
            k_res_4_max=k
            res_4_max_alone=res(j,k,i,4)
            vol_4_max_alone=vol(j,k,i)
          end if
          if (real(ccabs(res(j,k,i,5)/vol(j,k,i))) .gt. 
     +       real(res_5_max)) then
            res_5_max = ccabs(res(j,k,i,5)/vol(j,k,i))
            i_res_5_max=i
            j_res_5_max=j
            k_res_5_max=k
            res_5_max_alone=res(j,k,i,5)
            vol_5_max_alone=vol(j,k,i)
          end if
          sum_1_res=sum_1_res+ccabs(res(j,k,i,1)/vol(j,k,i))
          sum_2_res=sum_2_res+ccabs(res(j,k,i,2)/vol(j,k,i))
          sum_3_res=sum_3_res+ccabs(res(j,k,i,3)/vol(j,k,i))
          sum_4_res=sum_4_res+ccabs(res(j,k,i,4)/vol(j,k,i))
          sum_5_res=sum_5_res+ccabs(res(j,k,i,5)/vol(j,k,i))
        enddo
      enddo
      sum_1_res=sum_1_res/((jdim-5)*(kdim-5)*(idim-1))
      sum_2_res=sum_2_res/((jdim-5)*(kdim-5)*(idim-1))
      sum_3_res=sum_3_res/((jdim-5)*(kdim-5)*(idim-1))
      sum_4_res=sum_4_res/((jdim-5)*(kdim-5)*(idim-1))
      sum_5_res=sum_5_res/((jdim-5)*(kdim-5)*(idim-1))
c
      write(301,'(''# mean flow residuals - truncation error'',
     +  '' if iexact_trunc=1'')')
      write(301,'(''variables="N","1/N","sqrt(1/N)",'')')
      write(301,'(''"res1_max","i_res1_max","j_res1_max",'',
     +  ''"k_res1_max","res1_L1norm","res_1_max_alone",'',
     +  ''"vol_1_max_alone",'')')
      write(301,'(''"res2_max","i_res2_max","j_res2_max",'',
     +  ''"k_res2_max","res2_L1norm","res_2_max_alone",'',
     +  ''"vol_2_max_alone",'')')
      write(301,'(''"res3_max","i_res3_max","j_res3_max",'',
     +  ''"k_res3_max","res3_L1norm","res_3_max_alone",'',
     +  ''"vol_3_max_alone",'')')
      write(301,'(''"res4_max","i_res4_max","j_res4_max",'',
     +  ''"k_res4_max","res4_L1norm","res_4_max_alone",'',
     +  ''"vol_4_max_alone",'')')
      write(301,'(''"res5_max","i_res5_max","j_res5_max",'',
     +  ''"k_res5_max","res5_L1norm","res_5_max_alone",'',
     +  ''"vol_5_max_alone"'')')
      inumtot=(jdim-1)*(kdim-1)
      hhhh=1./(float((jdim-1)*(kdim-1)))
      sqhhhh=sqrt(hhhh)
      write(301,'(i9,3e18.8,3i5,3e18.8,e18.8,3i5,3e18.8,e18.8,3i5,3e18.8
     +           ,e18.8,3i5,3e18.8,e18.8,3i5,3e18.8)')
     +         inumtot,hhhh,sqhhhh,
     +         res_1_max,i_res_1_max,j_res_1_max,k_res_1_max,sum_1_res,
     +         res_1_max_alone,vol_1_max_alone,
     +         res_2_max,i_res_2_max,j_res_2_max,k_res_2_max,sum_2_res,
     +         res_2_max_alone,vol_2_max_alone,
     +         res_3_max,i_res_3_max,j_res_3_max,k_res_3_max,sum_3_res,
     +         res_3_max_alone,vol_3_max_alone,
     +         res_4_max,i_res_4_max,j_res_4_max,k_res_4_max,sum_4_res,
     +         res_4_max_alone,vol_4_max_alone,
     +         res_5_max,i_res_5_max,j_res_5_max,k_res_5_max,sum_5_res,
     +         res_5_max_alone,vol_5_max_alone
c   ***
c     i=1
c     do j=1,jdim-1
c       do k=1,kdim-1
c         xtemp(j,k)=0.25*(x(j,k,i)+x(j+1,k,i)+x(j,k+1,i)+x(j+1,k+1,i))
c         ztemp(j,k)=0.25*(z(j,k,i)+z(j+1,k,i)+z(j,k+1,i)+z(j+1,k+1,i))
c       enddo
c     enddo
c     ngd=1
c     write(801,'(i5)') ngd
c     write(801,'(2i5)') jdim-1,kdim-1
c     write(801,'(5e15.5)') ((xtemp(j,k),j=1,jdim-1),k=1,kdim-1),
c    +                      ((ztemp(j,k),j=1,jdim-1),k=1,kdim-1)
c     ifive=5
c     write(802,'(i5)') ngd
c     write(802,'(3i5)') jdim-1,kdim-1,ifive
c     write(802,'(5e15.5)') (((res(j,k,i,n),j=1,jdim-1),k=1,kdim-1),
c    +                      n=1,ifive)
c   ***
      return
      end
c
      subroutine exact_turb_res(jdim,kdim,idim,vol,res,iexact_trunc,
     +  iexact_disc)
c==================================== EXACT_TURB_RES =========================80
c
c Routine to get and print out turb residuals
c
c=============================================================================80
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension vol(jdim,kdim,idim-1), res(jdim,kdim,idim)
c
c   residuals
      res_max=0.
      sum_res=0.
      i=1
      do j=3,jdim-3
        do k=3,kdim-3
          if (real(ccabs(res(j,k,i))) .gt. real(res_max)) then
            res_max = ccabs(res(j,k,i))
            i_res_max=i
            j_res_max=j
            k_res_max=k
            res_max_alone=res(j,k,i)
            vol_max_alone=vol(j,k,i)
          end if
          sum_res=sum_res+ccabs(res(j,k,i))
        enddo
      enddo
      sum_res=sum_res/((jdim-5)*(kdim-5)*(idim-1))
c
      write(302,'(''# turbulence residuals - truncation error'',
     +  '' if iexact_trunc=1'')')
      write(302,'(''variables="N","1/N","sqrt(1/N)",'')')
      write(302,'(''"res_max","i_res_max","j_res_max",'',
     +  ''"k_res_max","res_L1norm","res_max_alone",'',
     +  ''"vol_max_alone"'')')
      inumtot=(jdim-1)*(kdim-1)
      hhhh=1./(float((jdim-1)*(kdim-1)))
      sqhhhh=sqrt(hhhh)
      write(302,'(i9,3e18.8,3i5,3e18.8)')
     +         inumtot,hhhh,sqhhhh,
     +         res_max,i_res_max,j_res_max,k_res_max,sum_res,
     +         res_max_alone,vol_max_alone
      return
      end
c
      subroutine exact_flow_force(jdim,kdim,idim,x,y,z,vol,res,smin,
     +  iexact_trunc,iexact_disc)
c==================================== EXACT_FLOW_FORCE =======================80
c
c Driver routine to compute the exact forcing terms as a function of (x,y,z).
c
c Note: When distf is non-zero, the turbulent forcing terms are computed,
c            otherwise skipped.
c
c=============================================================================80
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension x(jdim,kdim,idim), y(jdim,kdim,idim), z(jdim,kdim,idim),
     +          vol(jdim,kdim,idim-1), res(jdim,kdim,idim-1,5)
      dimension smin(jdim-1,kdim-1,idim-1)
      dimension f(6)
c
      i_convert_q = 0
      i_forcing = 1
      i_gradient = 0
      nvar=6
      do i=1,idim-1
        do j=1,jdim-1
          do k=1,kdim-1
            do n=1,nvar
              f(n)=0.
            enddo
            distf=smin(j,k,i)
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              res(j,k,i,n)=res(j,k,i,n)-vol(j,k,i)*f(n)
            enddo
          enddo
        enddo
      enddo
      return
      end
c
      subroutine exact_turb_force(jdim,kdim,idim,x,y,z,vol,res,smin,
     +  iexact_trunc,iexact_disc)
c==================================== EXACT_TURB_FORCE =======================80
c
c Driver routine to compute the exact turb forcing terms as a function of (x,y,z).
c
c Note: When distf is non-zero, the turbulent forcing terms are computed,
c            otherwise skipped.
c
c=============================================================================80
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension x(jdim,kdim,idim), y(jdim,kdim,idim), z(jdim,kdim,idim),
     +          vol(jdim,kdim,idim-1), res(jdim,kdim,idim), 
     +          smin(jdim-1,kdim-1,idim-1)
      dimension f(6)
c
      i_convert_q = 0
      i_forcing = 1
      i_gradient = 0
      nvar=6
      do i=1,idim-1
        do j=1,jdim-1
          do k=1,kdim-1
            do n=1,nvar
              f(n)=0.
            enddo
            distf=smin(j,k,i)
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            res(j,k,i)=res(j,k,i)+f(6)
          enddo
        enddo
      enddo
      return
      end
c
      subroutine exact_flow_q(jdim,kdim,idim,x,y,z,q,
     +  iexact_trunc,iexact_disc)
c==================================== EXACT_FLOW_Q ===========================80
c
c Driver routine to compute the exact q terms as a function of (x,y,z).
c
c Note: When distf is non-zero, the turbulent forcing terms are computed,
c            otherwise skipped.
c
c=============================================================================80
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension x(jdim,kdim,idim), y(jdim,kdim,idim), z(jdim,kdim,idim),
     +          q(jdim,kdim,idim,5)
      dimension f(6)
c
      i_convert_q = 0
      i_forcing = 0
      i_gradient = 0
      nvar=6
      distf=0.
      do i=1,idim-1
        do j=1,jdim-1
          do k=1,kdim-1
            do n=1,nvar
              f(n)=0.
            enddo
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              q(j,k,i,n)=f(n)
            enddo
          enddo
        enddo
      enddo
      return
      end
c
      subroutine zero_resid_ring(jdim,kdim,idim,res,jj,kk,ii,numeq,
     +                           numrows,iexact_trunc,iexact_disc)
c==================================== ZERO_RESID_RING ========================80
c
c Routine to zero residual in outer numrow rings - valid for 2-D only
c
c=============================================================================80
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension res(jj,kk,ii,numeq)
c
      do n=1,numeq
        do i=1,idim-1
          do j=1,numrows
            do k=1,kdim-1
              res(j,k,i,n)=0.0
            enddo
          enddo
          do j=jdim-numrows,jdim-1
            do k=1,kdim-1
              res(j,k,i,n)=0.0
            enddo
          enddo
          do j=1,jdim-1
            do k=1,numrows
              res(j,k,i,n)=0.0
            enddo
          enddo
          do j=1,jdim-1
            do k=kdim-numrows,kdim-1
              res(j,k,i,n)=0.0
            enddo
          enddo
        enddo
      enddo
      return
      end
c
      subroutine exact_flow_q_ring(jdim,kdim,idim,x,y,z,q,
     +  iexact_trunc,iexact_disc)
c==================================== EXACT_FLOW_Q_RING ======================80
c
c Driver routine to compute the exact q terms as a function of (x,y,z).
c done only in the outer ring of 2 cells on all sides
c
c Note: When distf is non-zero, the turbulent forcing terms are computed,
c            otherwise skipped.
c
c=============================================================================80
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension x(jdim,kdim,idim), y(jdim,kdim,idim), z(jdim,kdim,idim),
     +          q(jdim,kdim,idim,5)
      dimension f(6)
c
      i_convert_q = 0
      i_forcing = 0
      i_gradient = 0
      nvar=6
      distf=0.0
      do i=1,idim-1
        do j=1,2
          do k=1,kdim-1
            do n=1,nvar
              f(n)=0.
            enddo
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              q(j,k,i,n)=f(n)
            enddo
          enddo
        enddo
        do j=jdim-2,jdim-1
          do k=1,kdim-1
            do n=1,nvar
              f(n)=0.
            enddo
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              q(j,k,i,n)=f(n)
            enddo
          enddo
        enddo
        do j=1,jdim-1
          do k=1,2
            do n=1,nvar
              f(n)=0.
            enddo
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              q(j,k,i,n)=f(n)
            enddo
          enddo
        enddo
        do j=1,jdim-1
          do k=kdim-2,kdim-1
            do n=1,nvar
              f(n)=0.
            enddo
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              q(j,k,i,n)=f(n)
            enddo
          enddo
        enddo
      enddo
      return
      end
c
      subroutine exact_turb_q(jdim,kdim,idim,x,y,z,q,smin,vist3d,
     +  iexact_trunc,iexact_disc)
c==================================== EXACT_TURB_Q ===========================80
c
c Driver routine to compute the exact turb q terms as a function of (x,y,z).
c
c Note: When distf is non-zero, the turbulent forcing terms are computed,
c            otherwise skipped.
c
c=============================================================================80
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension x(jdim,kdim,idim), y(jdim,kdim,idim), z(jdim,kdim,idim),
     +          q(jdim,kdim,idim), smin(jdim-1,kdim-1,idim-1)
      dimension vist3d(jdim,kdim,idim)
      dimension f(6)
c
      i_convert_q = 0
      i_forcing = 0
      i_gradient = 0
      nvar=6
      do i=1,idim-1
        do j=1,jdim-1
          do k=1,kdim-1
            do n=1,nvar
              f(n)=0.
            enddo
            distf = smin(j,k,i)
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            q(j,k,i)=f(6)
            vist3d(j,k,i)=xmut
          enddo
        enddo
      enddo
      return
      end
c
      subroutine exact_turb_q_ring(jdim,kdim,idim,x,y,z,q,smin,vist3d,
     +  iexact_trunc,iexact_disc)
c==================================== EXACT_TURB_Q_RING ======================80
c
c Driver routine to compute the exact turb q terms as a function of (x,y,z).
c done only in the outer ring of 2 cells on all sides
c
c Note: When distf is non-zero, the turbulent forcing terms are computed,
c            otherwise skipped.
c
c=============================================================================80
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension x(jdim,kdim,idim), y(jdim,kdim,idim), z(jdim,kdim,idim),
     +          q(jdim,kdim,idim), smin(jdim-1,kdim-1,idim-1)
      dimension vist3d(jdim,kdim,idim)
      dimension f(6)
c
      i_convert_q = 0
      i_forcing = 0
      i_gradient = 0
      nvar=6
      do i=1,idim-1
        do j=1,2
          do k=1,kdim-1
            do n=1,nvar
              f(n)=0.
            enddo
            distf = smin(j,k,i)
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            q(j,k,i)=f(6)
            vist3d(j,k,i)=xmut
          enddo
        enddo
        do j=jdim-2,jdim-1
          do k=1,kdim-1
            do n=1,nvar
              f(n)=0.
            enddo
            distf = smin(j,k,i)
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            q(j,k,i)=f(6)
            vist3d(j,k,i)=xmut
          enddo
        enddo
        do j=1,jdim-1
          do k=1,2
            do n=1,nvar
              f(n)=0.
            enddo
            distf = smin(j,k,i)
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            q(j,k,i)=f(6)
            vist3d(j,k,i)=xmut
          enddo
        enddo
        do j=1,jdim-1
          do k=kdim-2,kdim-1
            do n=1,nvar
              f(n)=0.
            enddo
            distf = smin(j,k,i)
            xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            q(j,k,i)=f(6)
            vist3d(j,k,i)=xmut
          enddo
        enddo
      enddo
      return
      end
c
      subroutine exact_flow_bc(jdim,kdim,idim,x,y,z,qi0,qj0,qk0,nface,
     +                           bcj,bck,bci,smin,tj0,tk0,ti0,
     +                           vj0,vk0,vi0,nummem,
     +                           iexact_trunc,iexact_disc)
c==================================== EXACT_FLOW_BC ==========================80
c
c Driver routine to compute the exact BC terms as a function of (x,y,z).
c
c Note: When distf is non-zero, the turbulent forcing terms are computed,
c            otherwise skipped.
c
c=============================================================================80
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      common /twod/ i2d
      common /maxiv/ ivmx
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
c
      dimension x(jdim,kdim,idim), y(jdim,kdim,idim), z(jdim,kdim,idim),
     +          qi0(jdim,kdim,5,4),qj0(kdim,idim-1,5,4),
     +          qk0(jdim,idim-1,5,4)
      dimension tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),
     .          ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),
     .          vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension smin(jdim-1,kdim-1,idim-1)
      dimension f(6)
c
      i_convert_q = 0
      i_forcing = 0
      i_gradient = 0
      nvar=6
c  QI0 lo
      if (i2d .ne. 1) then
      if (nface.eq.1) then
      do j=1,jdim-1
        do k=1,kdim-1
          bci(j,k,1)   = 0.0
          do n=1,nvar
            f(n)=0.
          enddo
            i=1
            distf1=smin(j,k,i)
            xx1=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy1=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz1=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            i=2
            distf2=smin(j,k,i)
            xx2=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy2=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz2=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            xx=2.*xx1-xx2
            yy=2.*yy1-yy2
            zz=2.*zz1-zz2
            distf=2.*distf1-distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qi0(j,k,n,1)=f(n)
c             qi0(j,k,n,2)=qi0(j,k,n,1)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                ti0(j,k,n,1)=f(6)
c               ti0(j,k,n,2)=ti0(j,k,n,1)
              enddo
            end if
            if (ivmx .ge. 2) then
              vi0(j,k,1,1)=xmut
c             vi0(j,k,1,2)=vi0(j,k,1,1)
            end if
            xx=3.*xx1-2.*xx2
            yy=3.*yy1-2.*yy2
            zz=3.*zz1-2.*zz2
            distf=3.*distf1-2.*distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qi0(j,k,n,2)=f(n)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                ti0(j,k,n,2)=f(6)
              enddo
            end if
            if (ivmx .ge. 2) then
              vi0(j,k,1,2)=xmut
            end if
        enddo
      enddo
      end if
c  QI0 hi
      if (nface.eq.2) then
      do j=1,jdim-1
        do k=1,kdim-1
          bci(j,k,2)   = 0.0
          do n=1,nvar
            f(n)=0.
          enddo
            i=idim-1
            distf1=smin(j,k,i)
            xx1=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy1=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz1=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            i=idim-2
            distf2=smin(j,k,i)
            xx2=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy2=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz2=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            xx=2.*xx1-xx2
            yy=2.*yy1-yy2
            zz=2.*zz1-zz2
            distf=2.*distf1-distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qi0(j,k,n,3)=f(n)
c             qi0(j,k,n,4)=qi0(j,k,n,3)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                ti0(j,k,n,3)=f(6)
c               ti0(j,k,n,4)=ti0(j,k,n,3)
              enddo
            end if
            if (ivmx .ge. 2) then
              vi0(j,k,1,3)=xmut
c             vi0(j,k,1,4)=vi0(j,k,1,3)
            end if
            xx=3.*xx1-2.*xx2
            yy=3.*yy1-2.*yy2
            zz=3.*zz1-2.*zz2
            distf=3.*distf1-2.*distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qi0(j,k,n,4)=f(n)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                ti0(j,k,n,4)=f(6)
              enddo
            end if
            if (ivmx .ge. 2) then
              vi0(j,k,1,4)=xmut
            end if
        enddo
      enddo
      end if
      end if
c  QJ0 lo
      if (nface .eq. 3) then
      do i=1,idim-1
        do k=1,kdim-1
          bcj(k,i,1)   = 0.0
          do n=1,nvar
            f(n)=0.
          enddo
            j=1
            distf1=smin(j,k,i)
            xx1=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy1=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz1=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            j=2
            distf2=smin(j,k,i)
            xx2=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy2=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz2=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            xx=2.*xx1-xx2
            yy=2.*yy1-yy2
            zz=2.*zz1-zz2
            distf=2.*distf1-distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qj0(k,i,n,1)=f(n)
c             qj0(k,i,n,2)=qj0(k,i,n,1)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                tj0(k,i,n,1)=f(6)
c               tj0(k,i,n,2)=tj0(k,i,n,1)
              enddo
            end if
            if (ivmx .ge. 2) then
              vj0(k,i,1,1)=xmut
c             vj0(k,i,1,2)=vj0(k,i,1,1)
            end if
            xx=3.*xx1-2.*xx2
            yy=3.*yy1-2.*yy2
            zz=3.*zz1-2.*zz2
            distf=3.*distf1-2.*distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qj0(k,i,n,2)=f(n)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                tj0(k,i,n,2)=f(6)
              enddo
            end if
            if (ivmx .ge. 2) then
              vj0(k,i,1,2)=xmut
            end if
        enddo
      enddo
      end if
c  QJO hi
      if (nface .eq. 4) then
      do i=1,idim-1
        do k=1,kdim-1
          bcj(k,i,2)   = 0.0
          do n=1,nvar
            f(n)=0.
          enddo
            j=jdim-1
            distf1=smin(j,k,i)
            xx1=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy1=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz1=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            j=jdim-2
            distf2=smin(j,k,i)
            xx2=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy2=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz2=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            xx=2.*xx1-xx2
            yy=2.*yy1-yy2
            zz=2.*zz1-zz2
            distf=2.*distf1-distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qj0(k,i,n,3)=f(n)
c             qj0(k,i,n,4)=qj0(k,i,n,3)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                tj0(k,i,n,3)=f(6)
c               tj0(k,i,n,4)=tj0(k,i,n,3)
              enddo
            end if
            if (ivmx .ge. 2) then
              vj0(k,i,1,3)=xmut
c             vj0(k,i,1,4)=vj0(k,i,1,3)
            end if
            xx=3.*xx1-2.*xx2
            yy=3.*yy1-2.*yy2
            zz=3.*zz1-2.*zz2
            distf=3.*distf1-2.*distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qj0(k,i,n,4)=f(n)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                tj0(k,i,n,4)=f(6)
              enddo
            end if
            if (ivmx .ge. 2) then
              vj0(k,i,1,4)=xmut
            end if
        enddo
      enddo
      end if
c  QK0 lo
      if (nface.eq.5) then
      do i=1,idim-1
        do j=1,jdim-1
          bck(j,i,1)   = 0.0
          do n=1,nvar
            f(n)=0.
          enddo
            k=1
            distf1=smin(j,k,i)
            xx1=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy1=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz1=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            k=2
            distf2=smin(j,k,i)
            xx2=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy2=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz2=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            xx=2.*xx1-xx2
            yy=2.*yy1-yy2
            zz=2.*zz1-zz2
            distf=2.*distf1-distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qk0(j,i,n,1)=f(n)
c             qk0(j,i,n,2)=qk0(j,i,n,1)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                tk0(j,i,n,1)=f(6)
c               tk0(j,i,n,2)=tk0(j,i,n,1)
              enddo
            end if
            if (ivmx .ge. 2) then
              vk0(j,i,1,1)=xmut
c             vk0(j,i,1,2)=vk0(j,i,1,1)
            end if
            xx=3.*xx1-2.*xx2
            yy=3.*yy1-2.*yy2
            zz=3.*zz1-2.*zz2
            distf=3.*distf1-2.*distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qk0(j,i,n,2)=f(n)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                tk0(j,i,n,2)=f(6)
              enddo
            end if
            if (ivmx .ge. 2) then
              vk0(j,i,1,2)=xmut
            end if
        enddo
      enddo
      end if
c  QKO hi
      if (nface.eq.6) then
      do i=1,idim-1
        do j=1,jdim-1
          bck(j,i,2)   = 0.0
          do n=1,nvar
            f(n)=0.
          enddo
            k=kdim-1
            distf1=smin(j,k,i)
            xx1=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy1=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz1=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            k=kdim-2
            distf2=smin(j,k,i)
            xx2=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
     +                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
     +                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
     +                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1))
            yy2=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
     +                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
     +                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
     +                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1))
            zz2=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
     +                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
     +                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
     +                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1))
            xx=2.*xx1-xx2
            yy=2.*yy1-yy2
            zz=2.*zz1-zz2
            distf=2.*distf1-distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qk0(j,i,n,3)=f(n)
c             qk0(j,i,n,4)=qk0(j,i,n,3)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                tk0(j,i,n,3)=f(6)
c               tk0(j,i,n,4)=tk0(j,i,n,3)
              enddo
            end if
            if (ivmx .ge. 2) then
              vk0(j,i,1,3)=xmut
c             vk0(j,i,1,4)=vk0(j,i,1,3)
            end if
            xx=3.*xx1-2.*xx2
            yy=3.*yy1-2.*yy2
            zz=3.*zz1-2.*zz2
            distf=3.*distf1-2.*distf2
#   ifdef CMPLX
c     analytic_compressible not compiled if complex used
      continue
#   else
            call analytic_compressible(xx,yy,zz,nvar,f,
     +        i_convert_q,i_forcing,i_gradient,distf,xmut,
     +        iexact_trunc,iexact_disc)
#   endif
            do n=1,5
              qk0(j,i,n,4)=f(n)
            enddo
            if (level .ge. lglobal .and. ivmx .ge. 4) then
              do n=1,nummem
                tk0(j,i,n,4)=f(6)
              enddo
            end if
            if (ivmx .ge. 2) then
              vk0(j,i,1,4)=xmut
            end if
        enddo
      enddo
      end if
      return
      end
c
#   ifdef CMPLX
c     do not compile all the exact routines if complex
#   else
      subroutine analytic_compressible(x, y, z, neq, q, i_convert_q,
     +  i_forcing, i_gradient, distf, xmut, iexact_trunc, iexact_disc)
c================================== ANALYTIC_COMPRESSIBLE  ===================80
c
c compressible perfect gas flow...analytic functions.
c
c=============================================================================80
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /reyue/ reue,tinf,ivisc(3)
      common /maxiv/ ivmx
      common /zero/ iexp
c
      dimension q(neq), f(neq), u_polyf(13), w_polyf(13), p_polyf(13)
c
c   Note: (10.**(-iexp) is machine zero)
      xminn=10.**(-iexp+1)
c
c Initialize q
      do n=1,neq
        q(n)=0.
      enddo 
c
c Spalart Constants
      sig    = 0.66667
      cv1    = 7.1
      ct3    = 1.2
      ct4    = 0.5
      vkar   = 0.41
      cb1    = 0.1355
      cb2    = 0.622
      cw1    = cb1/vkar**2+(1.+cb2)/sig
      cw2    = 0.3
      cw3    = 2.0
c
c Other variables
      xmre   = xmach / reue
      xmre_s = xmre / sig
      cgp    = xmre /(gm1*pr)
      cgpt   = xmre /(gm1*prt)
c
c Initialize
      call exact_zero(u, ux, uy, uz, uxx, uyy, uzz, uxy, uxz, uyz )
      call exact_zero(v, vx, vy, vz, vxx, vyy, vzz, vxy, vxz, vyz )
      call exact_zero(w, wx, wy, wz, wxx, wyy, wzz, wxy, wxz, wyz )
      call exact_zero(p, px, py, pz, pxx, pyy, pzz, pxy, pxz, pyz )
      call exact_zero(t, tx, ty, tz, txx, tyy, tzz, txy, txz, tyz )
c
c Choose manufactured solution
c 0=polynomial, 1=MS1(Lisbon), 2=MS2(Lisbon), 4=MS4(Lisbon)
      ims_sol=max(iexact_trunc, iexact_disc)
      if (ims_sol .eq. 999) then
        u_scale=xmach
        w_scale=xmach
        p_scale=1.0/gamma
        u_polyf(1)=1.
        u_polyf(2)=1.
        u_polyf(3)=0.
        u_polyf(4)=0.
        u_polyf(5)=0.
        u_polyf(6)=0.
        u_polyf(7)=0.
        u_polyf(8)=0.
        u_polyf(9)=0.
        u_polyf(10)=0.
        u_polyf(11)=0.
        u_polyf(12)=0.
        u_polyf(13)=0.
        w_polyf(1)=1.
        w_polyf(2)=1.
        w_polyf(3)=0.
        w_polyf(4)=0.
        w_polyf(5)=0.
        w_polyf(6)=0.
        w_polyf(7)=0.
        w_polyf(8)=0.
        w_polyf(9)=0.
        w_polyf(10)=0.
        w_polyf(11)=0.
        w_polyf(12)=0.
        w_polyf(13)=0.
        p_polyf(1)=1.
        p_polyf(2)=1.
        p_polyf(3)=0.
        p_polyf(4)=0.
        p_polyf(5)=0.
        p_polyf(6)=0.
        p_polyf(7)=0.
        p_polyf(8)=0.
        p_polyf(9)=0.
        p_polyf(10)=0.
        p_polyf(11)=0.
        p_polyf(12)=0.
        p_polyf(13)=0.
        call exact_polyfg(u_scale, u_polyf, x, y, z, u, ux, uy, uz, 
     +              uxx, uyy, uzz, uxy, uxz, uyz )
        call exact_polyfg(w_scale, w_polyf, x, y, z, w, wx, wy, wz, 
     +              wxx, wyy, wzz, wxy, wxz, wyz )
        call exact_polyfg(p_scale, p_polyf, x, y, z, p, px, py, pz, 
     +              pxx, pyy, pzz, pxy, pxz, pyz )
      else if (ims_sol .eq. 1 .or. ims_sol .eq. 2 .or.
     +         ims_sol .eq. 4) then
        xlisbon_p_scale = 1.0/gamma
        call exact_lisbon_ms2_u(xmach, x, z, u, ux, uz, uxx, uzz, uxz )
        call exact_lisbon_ms2_v(xmach, x, z, w, wx, wz, wxx, wzz, wxz )
        call exact_lisbon_ms2_p(xlisbon_p_scale, x, z, p, px, pz )
        p = p + xlisbon_p_scale
      else
        write(6,'('' MFG soln not chosen'')')
      end if
c
c Enforce constant adiabatic wall temperature (2D).
      call total_temperature_2d( t, tx, tz, txx, tzz, txz, 
     +                           u, ux, uz, uxx, uzz, uxz, 
     +                           w, wx, wz, wxx, wzz, wxz )
c
c Density
      rho  = gamma*( p/t )
      rhox = gamma*( px/t - p*tx/t**2 )
      rhoy = gamma*( py/t - p*ty/t**2 )
      rhoz = gamma*( pz/t - p*tz/t**2 )
c
c Molecular viscosity
      xmu    =  viscosity_law( t )
      dxmudt =  dviscosity_law( t )
      xmux  = dxmudt*tx
      xmuy  = dxmudt*ty
      xmuz  = dxmudt*tz
c
c Kinematic viscosity
      xnu   =  xmu/rho
      xnux  = (xmux - xnu*rhox)/rho
      xnuy  = (xmuy - xnu*rhoy)/rho
      xnuz  = (xmuz - xnu*rhoz)/rho
c
c Continuity equation
      d    = ux + vy + wz
      f(1) = rho*d + u*rhox + v*rhoy + w*rhoz
c
c Inviscid momentum equation
      f(2) = rho*( u*ux + v*uy + w*uz ) + px + u*f(1)
      f(3) = rho*( u*vx + v*vy + w*vz ) + py + v*f(1)
      f(4) = rho*( u*wx + v*wy + w*wz ) + pz + w*f(1)
c
c Inviscid energy equation
      h  =  t/gm1 + 0.5*(u**2 + v**2 + w**2)
      hx = tx/gm1 + u*ux + v*vx + w*wx
      hy = ty/gm1 + u*uy + v*vy + w*wy
      hz = tz/gm1 + u*uz + v*vz + w*wz
      rhoh  = rho*h
      rhohx = rhox*h + rho*hx
      rhohy = rhoy*h + rho*hy
      rhohz = rhoz*h + rho*hz
      f(5) = u*rhohx + v*rhohy + w*rhohz + rhoh*d
c
c   s2 = twice mean-strain-rate tensor
      s2_xx = ux + ux
      s2_xy = vx + uy
      s2_xz = wx + uz
      s2_yy = vy + vy
      s2_yz = wy + vz
      s2_zz = wz + wz
c
c   Laminar momentum equations
      xla_u = uxx + uyy + uzz
      xla_v = vxx + vyy + vzz
      xla_w = wxx + wyy + wzz
      d3   = 1./3.*d
      d3x  = 1./3.*( uxx + vxy + wxz )
      d3y  = 1./3.*( uxy + vyy + wyz )
      d3z  = 1./3.*( uxz + vyz + wzz )
c
c   Momentum equations - first term
      vf2 = - xmre*xmu*( xla_u + d3x )
      vf3 = - xmre*xmu*( xla_v + d3y )
      vf4 = - xmre*xmu*( xla_w + d3z )
c
c   Momentum equations - second term
      vf2 = vf2 - xmre*( xmux*s2_xx + xmuy*s2_xy + xmuz*s2_xz )
      vf3 = vf3 - xmre*( xmux*s2_xy + xmuy*s2_yy + xmuz*s2_yz )
      vf4 = vf4 - xmre*( xmux*s2_xz + xmuy*s2_yz + xmuz*s2_zz )
      vf2 = vf2 + xmre*2.0*d3*xmux
      vf3 = vf3 + xmre*2.0*d3*xmuy
      vf4 = vf4 + xmre*2.0*d3*xmuz
      f(2) = f(2) + vf2
      f(3) = f(3) + vf3
      f(4) = f(4) + vf4
c
c   Laminar energy equation
      xla_t = txx + tyy + tzz
      vf5  = u*vf2 + v*vf3 + w*vf4                                 
     +   - cgp*( xmu*xla_t + xmux*tx + xmuy*ty + xmuz*tz )         
     +   - xmre*xmu*( ux*s2_xx + uy*s2_xy + uz*s2_xz              
     +              + vx*s2_xy + vy*s2_yy + vz*s2_yz              
     +              + wx*s2_xz + wy*s2_yz + wz*s2_zz - 2.0*d3*d )
      f(5) = f(5) + vf5
c
c   Spalart turbulence model equation
      if (ivmx .eq. 5) then
c       Turbulence variable
        call exact_zero(vt, vtx, vty, vtz, vtxx, vtyy, vtzz, vtxy, 
     +                  vtxz, vtyz )
        term = 1.0
        if (ims_sol .eq. 1) then
          call exact_lisbon_ms1_vt(term, x, z, vt, vtx, vtz, vtxx, vtzz)
        else if (ims_sol .eq. 2) then
          call exact_lisbon_ms2_vt(term, x, z, vt, vtx, vtz, vtxx, vtzz)
        else if (ims_sol .eq. 4) then
          call exact_lisbon_ms4_vt(term, x, z, vt, vtx, vtz, vtxx, vtzz)
        else
          write(6,'('' Error... must have ims_sol=1,2, or 4 with SA'')')
        end if
        chi  =   vt/xnu
        chix = ( vtx - chi*xnux )/xnu
        chiy = ( vty - chi*xnuy )/xnu
        chiz = ( vtz - chi*xnuz )/xnu
        fv1  = chi**3 / (chi**3+cv1**3)
        fv1x = 3.0*chix*( 1.0 - fv1 )*chi**2/(chi**3+cv1**3)
        fv1y = 3.0*chiy*( 1.0 - fv1 )*chi**2/(chi**3+cv1**3)
        fv1z = 3.0*chiz*( 1.0 - fv1 )*chi**2/(chi**3+cv1**3)
        xmut  =  rho*fv1*vt
        xmutx = rhox*fv1*vt + rho*( fv1x*vt + fv1*vtx )
        xmuty = rhoy*fv1*vt + rho*( fv1y*vt + fv1*vty )
        xmutz = rhoz*fv1*vt + rho*( fv1z*vt + fv1*vtz )
c       The shear stress terms are taken as:
c
c       t_ij = (xmu + xmut)*[ 2S_ij - (2d/3)*delta_ij ]
c
c       Momentum equations - first term
        tf2 = - xmre*xmut*( xla_u + d3x )
        tf3 = - xmre*xmut*( xla_v + d3y )
        tf4 = - xmre*xmut*( xla_w + d3z )
c
c       Momentum equations - second term
        tf2 = tf2 - xmre*( xmutx*s2_xx + xmuty*s2_xy + xmutz*s2_xz )
        tf3 = tf3 - xmre*( xmutx*s2_xy + xmuty*s2_yy + xmutz*s2_yz )
        tf4 = tf4 - xmre*( xmutx*s2_xz + xmuty*s2_yz + xmutz*s2_zz )
        tf2 = tf2 + xmre*2.0*d3*xmutx
        tf3 = tf3 + xmre*2.0*d3*xmuty
        tf4 = tf4 + xmre*2.0*d3*xmutz
        f(2) = f(2) + tf2
        f(3) = f(3) + tf3
        f(4) = f(4) + tf4
c
c       Energy equation
        f(5) = f(5) + u*tf2 + v*tf3 + w*tf4                           
     +   - cgpt*( xmut*xla_t + xmutx*tx + xmuty*ty + xmutz*tz )     
     +   - xmre*xmut*( ux*s2_xx + uy*s2_xy + uz*s2_xz              
     +               + vx*s2_xy + vy*s2_yy + vz*s2_yz              
     +               + wx*s2_xz + wy*s2_yz + wz*s2_zz - 2.0*d3*d )
c
        if (distf .ne. 0.) then
          t_conv=1.0
          t_diff1=1.0
          t_diff2=1.0
          t_diff3=1.0
          t_prod=1.0
          t_dest=1.0
          distance = distf
          fv2 = 1.0 - chi/(1.0+chi*fv1)
          ft2 = ct3 * exp(-ct4*chi*chi)
          vkar2 = vkar * vkar
          bot = (vkar*distance)**2
          s   = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
          sw  = s + xmre*vt*fv2/bot
c   the following is slightly different between FUN3D and CFL3D:
c         sw  = max(sw,0.00001)
          sw  = max(sw,xminn)
c   ------------------------------------------------------------
          rr  = xmre*vt/(bot*sw)
c   the following is slightly different between FUN3D and CFL3D:
c         if(rr > 20.0) rr = 20.0
          if(rr > 10.0) rr = 10.0
c   ------------------------------------------------------------
          gg  = rr + cw2*(rr**6-rr)
c   the following is done in CFL3D, but not in FUN3D:
          gg = max(gg,xminn)
c   -------------------------------------------------
          fw  = gg * ((1.0+cw3**6)/(gg**6+cw3**6))**(1./6.)
          xla_vt = vtxx + vtyy + vtzz
c
c         Terms defined in accordance with CFL3D User's Manual
c         and prod/dest terms are grouped as in the S-A paper.
c         ...e.g., the prod term uses sw (s in CFL3D grouping).
c
c         Note: this form differs from Allmaras AIAA 1999-3336
c         by terms of order xmre_s*cb2*grad(xnu)*grad(vt)
          f_conv = u*vtx + v*vty + w*vtz
c   the following is slightly different between FUN3D and CFL3D:
c         f_prod = cb1*(1.0-ft2)*sw*vt
c         f_dest = xmre*( cb1*ft2/vkar2 - cw1*fw )*(vt/distance)**2
          f_prod = cb1*(1.0-ft2)*s*vt
          f_dest = xmre*( cb1*((1.-ft2)*fv2+ft2)/vkar2 - 
     +             cw1*fw )*(vt/distance)**2
c   ------------------------------------------------------------
c
c         ...Note: a model of pure diffusion (Laplacian operator)
c         ...is with t_diff1=1, t_diff2=t_diff3=0, xnu=constant
c         ...(Laplace equation is recovered with xmach=re=sig=1).
          cb21 = t_diff2 + t_diff3*cb2
          f_diff = xmre_s*( ( t_diff1*xnu + t_diff2*vt )*xla_vt  
     +                    + ( t_diff1*xnux + cb21*vtx )*vtx     
     +                    + ( t_diff1*xnuy + cb21*vty )*vty     
     +                    + ( t_diff1*xnuz + cb21*vtz )*vtz     )
          q(6) = vt
          f(6) = t_conv*f_conv                            
     +         - ( t_prod*f_prod + t_dest*f_dest )        
     +         - f_diff
        else
          q(6) = vt
          f(6) = 0.0
        end if
      end if
c
c Return requested result in q
      q(1) = rho
      q(2) = u
      q(3) = v
      q(4) = w
      q(5) = p
      if ( ivmx .eq. 5 ) q(6) = vt
      if (i_convert_q .eq. 1) then
        q(1) = rho
        q(2) = rho*u
        q(3) = rho*v
        q(4) = rho*w
        q(5) = p/gm1 + 0.5*rho*(u**2 + v**2 + w**2)
      else if (i_forcing .eq. 1) then
        q(1) = f(1)
        q(2) = f(2)
        q(3) = f(3)
        q(4) = f(4)
        q(5) = f(5)
        if ( ivmx .eq. 5 ) q(6) = f(6)
      else if (i_gradient .eq. 1) then
        q(1) = rhox
        q(2) = ux
        q(3) = vx
        q(4) = wx
        q(5) = px
        if ( ivmx .eq. 5 ) q(6) = vtx
      else if (i_gradient .eq. 2) then
        q(1) = rhoy
        q(2) = uy
        q(3) = vy
        q(4) = wy
        q(5) = py
        if ( ivmx .eq. 5 ) q(6) = vty
      else if(i_gradient .eq. 3) then
        q(1) = rhoz
        q(2) = uz
        q(3) = vz
        q(4) = wz
        q(5) = pz
        if ( ivmx .eq. 5 ) q(6) = vtz
      end if
c
      return
      end
c
      subroutine exact_zero(q, qx, qy, qz, qxx, qyy, qzz, qxy, qxz, qyz)
c=============================== EXACT_ZERO ==================================80
c
c Initialize manufactured solution and gradients to zero.
c
c=============================================================================80
      q = 0.0
      qx =  0.0
      qy =  0.0
      qz =  0.0
      qxx = 0.0
      qyy =  0.0
      qzz =  0.0
      qxy =  0.0
      qxz =  0.0
      qyz =  0.0
      return
      end
c
      subroutine total_temperature_2d(t, tx, ty, txx, tyy, txy, 
     +                                u, ux, uy, uxx, uyy, uxy,
     +                                v, vx, vy, vxx, vyy, vxy  )
c=============================== TOTAL_TEMPERATURE_2D ========================80
c
c Find temperature using constant total temperature.
c
c=============================================================================80

      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar

c   correction for consistency with CFL3D's TWTYPE=-1 for BC2004:
c     d = gm1*sqrt(pr)
      d = gm1
c   -------------------------------------------------------------

      t0 = 1.0 + 0.5*d*xmach*xmach

      t = t0 - 0.5*d*( u**2 + v**2 )

      tx = -d*( u*ux + v*vx )
      ty = -d*( u*uy + v*vy )

      txx = -d*( u*uxx + ux*ux + v*vxx + vx*vx )
      tyy = -d*( u*uyy + uy*uy + v*vyy + vy*vy )
      txy = -d*( u*uxy + uy*ux + v*vxy + vy*vx )
      return
      end
c
      function viscosity_law( t )
      common /reyue/ reue,tinf,ivisc(3)
      suth_const=198.6
      cstar=suth_const/tinf
      viscosity_law = (1.0+cstar)/(t+cstar)*t**1.5
      return
      end
c
      function dviscosity_law( t )
      common /reyue/ reue,tinf,ivisc(3)
      suth_const=198.6
      cstar=suth_const/tinf
      tstar  = t+cstar
      tstari = 1.0/tstar
      viscosity_law = (1.0+cstar)*tstari*t**1.5
      dviscosity_law = viscosity_law*( 0.5 + 1.5*cstar/t )*tstari
      return
      end
c
      subroutine exact_lisbon_ms2_forcing(beta, x, y, f )
c=============================== EXACT_LISBON_MS2_FORCING ====================80
c
c Add exact solution and gradients from Lisbon backward facing step ms2.
c
c=============================================================================80
      dimension f(5)
c   !...account for nondimensionalization in CFL3D
      emu=1.e-6
      f(1) = beta*( dudxms(x,y) + dvdyms(x,y) )
      f(2) = smxsams(x,y)
      f(3) = 0.0
      f(4) = smysams(x,y)
      f(5) = ssams(x,y)/emu
      return
      end
c
      subroutine exact_lisbon_ms2_u(scale, x, y, q, qx, qy, qxx, 
     +                              qyy, qxy)
c=============================== EXACT_LISBON_MS2_U ==========================80
c
c Add exact solution and gradients from Lisbon backward facing step ms2.
c
c=============================================================================80
      q = q + scale*ums(x,y)
      qx =  qx + scale*dudxms(x,y)
      qy =  qy + scale*dudyms(x,y)
      qxx = qxx +  scale*dudx2ms(x,y)
      qyy = qyy +  scale*dudy2ms(x,y)
      qxy = qxy +  scale*dudxyms(x,y)
      return
      end
c
      subroutine exact_lisbon_ms2_v(scale, x, y, q, qx, qy, qxx, 
     +                              qyy, qxy )
c=============================== EXACT_LISBON_MS2_V ==========================80
c
c Add exact solution and gradients from Lisbon backward facing step ms2.
c
c=============================================================================80
      q = q + scale*vms(x,y)
      qx =  qx + scale*dvdxms(x,y)
      qy =  qy + scale*dvdyms(x,y)
      qxx = qxx +  scale*dvdx2ms(x,y)
      qyy = qyy +  scale*dvdy2ms(x,y)
      qxy = qxy +  scale*dvdxyms(x,y)
      return
      end 
c
      subroutine exact_lisbon_ms2_p(scale, x, y, q, qx, qy )
c=============================== EXACT_LISBON_MS2_P ==========================80
c
c Add exact solution and gradients from Lisbon backward facing step ms2.
c
c=============================================================================80
      q = q + scale*pms(x,y)
      qx =  qx + scale*dpdxms(x,y)
      qy =  qy + scale*dpdyms(x,y)
      return
      end
c
      subroutine exact_lisbon_ms1_vt(scale, x, y, q, qx, qy, qxx, qyy )
c=============================== EXACT_LISBON_MS1_VT =========================80
c
c Add exact solution and gradients from Lisbon backward facing step ms1.
c
c=============================================================================80
      emu=1.e-6
      q = q + scale*em1ms(x,y)/emu
      qx =  qx + scale*dem1dxms(x,y)/emu
      qy =  qy + scale*dem1dyms(x,y)/emu
      qxx = qxx +  scale*dem1dx2ms(x,y)/emu
      qyy = qyy +  scale*dem1dy2ms(x,y)/emu
      return
      end
c
      subroutine exact_lisbon_ms2_vt(scale, x, y, q, qx, qy, qxx, qyy )
c=============================== EXACT_LISBON_MS2_VT =========================80
c
c Add exact solution and gradients from Lisbon backward facing step ms2.
c
c=============================================================================80
      emu=1.e-6
      q = q + scale*em2ms(x,y)/emu
      qx =  qx + scale*dem2dxms(x,y)/emu
      qy =  qy + scale*dem2dyms(x,y)/emu
      qxx = qxx +  scale*dem2dx2ms(x,y)/emu
      qyy = qyy +  scale*dem2dy2ms(x,y)/emu
      return
      end
c
      subroutine exact_lisbon_ms4_vt(scale, x, y, q, qx, qy, qxx, qyy )
c=============================== EXACT_LISBON_MS2_VT =========================80
c
c Add exact solution and gradients from Lisbon backward facing step ms2.
c
c=============================================================================80
      emu=1.e-6
      q = q + scale*em4ms(x,y)/emu
      qx =  qx + scale*dem4dxms(x,y)/emu
      qy =  qy + scale*dem4dyms(x,y)/emu
      qxx = qxx +  scale*dem4dx2ms(x,y)/emu
      qyy = qyy +  scale*dem4dy2ms(x,y)/emu
      return
      end
c
      function ums(x,y)
      sigma=4.0
      u1=1.0
      eta=sigma*y/x
      ums=user_derf(eta)
      return
      end
c
      function dudxms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      dudxms=-u1*sigma*y*2.0/(x*x*sqrt(pi))*exp(-eta*eta)
      return
      end
c
      function dudyms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      dudyms=u1*sigma*2.0/(x*sqrt(pi))*exp(-eta*eta)
      return
      end
c
      function dudx2ms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      t1=4.0*u1/sqrt(pi)
      x2=x*x
      dudx2ms=t1*(eta/x2)*exp(-eta*eta)*(1.0-eta*eta)
      return
      end
c
      function dudy2ms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      t1=4.0*u1/sqrt(pi)
      x2=sigma/x
      dudy2ms=-t1*x2*x2*eta*exp(-eta*eta)
      return
      end
c
      function dudxyms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      t1=2.0*u1/sqrt(pi)*sigma/x/x
      dudxyms=t1*exp(-eta*eta)*(2.0*eta*eta-1.0)
      return
      end
c
      function dudx2yms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      eta2=eta*eta
      t1=4.0*u1/sqrt(pi)*sigma/x/x/x
      dudx2yms=t1*exp(-eta2)*(1.0+2.0*eta2*(eta2-2.5))
      return
      end
c
      function dudy3ms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      x2=sigma/x
      t1=4.0*u1/sqrt(pi)*x2*x2*x2
      dudy3ms=t1*exp(-eta*eta)*(2.0*eta*eta-1.0)
      return
      end 
c
      function vms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      t1=u1/(sigma*sqrt(pi))
      vms=t1*(1.0-exp(-eta*eta))
      return
      end
c
      function dvdxms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      dvdxms=-u1*sigma*y*y*2.0/(x*x*x*sqrt(pi))*exp(-eta*eta)
      return
      end
c
      function dvdyms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      dvdyms=u1*sigma*y*2.0/(x*x*sqrt(pi))*exp(-eta*eta)
      return
      end
c
      function dvdx2ms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      t1=2.0*u1/sqrt(pi)
      x2=sigma*y*y/(x*x*x*x)
      dvdx2ms=t1*x2*exp(-eta*eta)*(3.0-2.0*eta*eta)
      return
      end
c
      function dvdy2ms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      t1=2.0*u1/sqrt(pi)
      x2=sigma/(x*x)
      dvdy2ms=t1*x2*exp(-eta*eta)*(1.0-2.0*eta*eta)
      return
      end
c
      function dvdxyms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      t1=4.0*u1/sqrt(pi)*sigma*y/x/x/x
      dvdxyms=t1*exp(-eta*eta)*(eta*eta-1.0)
      return
      end
c
      function dvdxy2ms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      eta2=eta*eta
      t1=4.0*u1/sqrt(pi)*sigma/x/x/x
      dvdxy2ms=t1*exp(-eta2)*(2.0*eta2*(2.5-eta2)-1.0)
      return
      end
c
      function dvdx3ms(x,y)
      sigma=4.0
      pi=acos(-1.)
      u1=1.0
      eta=sigma*y/x
      eta2=eta*eta
      t1=4.0*u1/sqrt(pi)
      x2=eta*y/(x*x*x*x)
      dvdx3ms=t1*x2*exp(-eta2)*(2.0*eta2*(4.5-eta2)-6.0)
      return
      end
c
      function pms(x,y)
      pms=0.5*log(2.0*x-x*x+0.25)*log(4.0*y*y*y-3.0*y*y+1.25)
      return
      end
c
      function dpdxms(x,y)
      dpdxms=0.5*(2.0-2.0*x)/(2.0*x-x*x+0.25)
     +       *log(4.0*y*y*y-3.0*y*y+1.25)
      return
      end
c
      function dpdyms(x,y)
      dpdyms=0.5*(12.0*y*y-6.0*y)/(4.0*y*y*y-3.0*y*y+1.25)
     +       *log(2.0*x-x*x+0.25)
      return
      end
c
      function strainms(x,y)
      u1x=dudxms(x,y)
      u1y=dudyms(x,y)
      u2x=dvdxms(x,y)
      u2y=dvdyms(x,y)
      strainms=2.0*(u1x*u1x+u2y*u2y)+(u1y+u2x)*(u1y+u2x)
      strainms=sqrt(strainms)
      return
      end
c
      function vortms(x,y)
      u1y=dudyms(x,y)
      u2x=dvdxms(x,y)
      vortms=(u1y-u2x)*(u1y-u2x)
      vortms=sqrt(vortms)
      return
      end
c
      function dvodxms(x,y)
      vort=sign(1.0,dudyms(x,y)-dvdxms(x,y))
      dvodxms=vort*(dudxyms(x,y)-dvdx2ms(x,y))
      return
      end
c
      function dvodyms(x,y)
      vort=sign(1.0,dudyms(x,y)-dvdxms(x,y))
      dvodyms=vort*(dudy2ms(x,y)-dvdxyms(x,y))
      return
      end
c
      function dvodx2ms(x,y)
      vort=sign(1.0,dudyms(x,y)-dvdxms(x,y))
      dvodx2ms=vort*(dudx2yms(x,y)-dvdx3ms(x,y))
      return
      end
c
      function dvody2ms(x,y)
      vort=sign(1.0,dudyms(x,y)-dvdxms(x,y))
      dvody2ms=vort*(dudy3ms(x,y)-dvdxy2ms(x,y))
      return
      end
c
      function smxsams(x,y)
c  one-equation turbulence models
c  spalart & allmaras
      emu=1.e-6
c  x momentum
      u1=ums(x,y)
      u2=vms(x,y)
      em=eddysams(x,y)
c  convection
      cv1=u1*dudxms(x,y)
      cv2=u2*dudyms(x,y)
      cvx=cv1+cv2
c  difusion
      df1=-(emu+em)*(dudx2ms(x,y)+dudy2ms(x,y))
      df2=-2.0*desadxms(x,y)*dudxms(x,y)
      df3=-desadyms(x,y)*(dudyms(x,y)+dvdxms(x,y))
      dfx=df1+df2+df3
c  pressure
      dpx=dpdxms(x,y)
c
      smxsams=cvx+dfx+dpx
      return
      end
c
      function smysams(x,y)
      emu=1.e-6
c  y momentum
      u1=ums(x,y)
      u2=vms(x,y)
      em=eddysams(x,y)
c  convection
      cv1=u1*dvdxms(x,y)
      cv2=u2*dvdyms(x,y)
      cvy=cv1+cv2
c  difusion
      df1=-(emu+em)*(dvdx2ms(x,y)+dvdy2ms(x,y))
      df2=-2.0*desadyms(x,y)*dvdyms(x,y)
      df3=-desadxms(x,y)*(dudyms(x,y)+dvdxms(x,y))
      dfy=df1+df2+df3
c  pressure
      dpy=dpdyms(x,y)
c
      smysams=cvy+dfy+dpy
      return
      end
c
      function smxmtms(x,y)
c  Menter model
      emu=1.e-6
c  x momentum
      u1=ums(x,y)
      u2=vms(x,y)
      em=eddymtms(x,y)
c  convection
      cv1=u1*dudxms(x,y)
      cv2=u2*dudyms(x,y)
      cvx=cv1+cv2
c  difusion
      df1=-(emu+em)*(dudx2ms(x,y)+dudy2ms(x,y))
      df2=-2.0*demtdxms(x,y)*dudxms(x,y)
      df3=-demtdyms(x,y)*(dudyms(x,y)+dvdxms(x,y))
      dfx=df1+df2+df3
c  pressure
      dpx=dpdxms(x,y)
c
      smxmtms=cvx+dfx+dpx
      return
      end
c
      function smymtms(x,y)
      emu=1.e-6
c  y momentum
      u1=ums(x,y)
      u2=vms(x,y)
      em=eddymtms(x,y)
c  convection
      cv1=u1*dvdxms(x,y)
      cv2=u2*dvdyms(x,y)
      cvy=cv1+cv2
c  difusion
      df1=-(emu+em)*(dvdx2ms(x,y)+dvdy2ms(x,y))
      df2=-2.0*demtdyms(x,y)*dvdyms(x,y)
      df3=-demtdxms(x,y)*(dudyms(x,y)+dvdxms(x,y))
      dfy=df1+df2+df3
c  pressure
      dpy=dpdyms(x,y)
c
      smymtms=cvy+dfy+dpy
      return
      end
c
c  turbulence models
c
c  dependent variable of the one-equation models
c
      function em1ms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      eta=sigeta*y/x
      t0=sqrt(2.0)*exp(0.5)*emmax
      em1ms=t0*eta*exp(-eta*eta)
      return
      end
c
      function em2ms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      eta=sigeta*y/x
      t1=exp(1.0)*emmax
      em2ms=t1*eta*eta*exp(-eta*eta)
      return
      end
c
      function em4ms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      eta=sigeta*y/x
      t1=0.25*exp(2.)*emmax
      em4ms=t1*eta*eta*eta*eta*exp(-eta*eta)
      return
      end
c
      function dem1dxms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      t0=sqrt(2.0)*exp(0.5)*emmax
      eta=sigeta*y/x
      eta_x = -eta/x
      dem1dxms=t0*eta_x*exp(-eta*eta)*( 1.0 - 2.0*eta*eta )
      return
      end
c
      function dem2dxms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      eta=sigeta*y/x
      t1=2.0*exp(1.0)*emmax/x
      dem2dxms=t1*eta*eta*exp(-eta*eta)*(eta*eta-1.0)
      return
      end
c
      function dem4dxms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      eta=sigeta*y/x
      t1=0.5*exp(2.)*emmax/x
      dem4dxms=t1*eta*eta*eta*eta*exp(-eta*eta)*(eta*eta-2.0)
      return
      end
c
      function dem1dyms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      t0=sqrt(2.0)*exp(0.5)*emmax
      eta   = sigeta*y/x
      eta_y = sigeta/x
      dem1dyms=t0*eta_y*exp(-eta*eta)*( 1.0 - 2.0*eta*eta )
      return
      end
c
      function dem2dyms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      eta=sigeta*y/x
      t1=2.0*exp(1.0)*emmax*sigeta/x
      dem2dyms=t1*eta*exp(-eta*eta)*(1.0-eta*eta)
      return
      end
c
      function dem4dyms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      eta=sigeta*y/x
      t1=0.5*exp(2.)*emmax*sigeta/x
      dem4dyms=t1*eta*eta*eta*exp(-eta*eta)*(2.0-eta*eta)
      return
      end
c
      function dem1dx2ms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      t0=sqrt(2.0)*exp(0.5)*emmax
      eta    = sigeta*y/x
      eta_x  = -eta/x
      eta_xx = -2.0*eta_x/x
      dem1dx2ms=t0*exp(-eta*eta)*(        eta_xx*(  1.0 - 2.0*eta*eta )
     +                       + eta*eta_x*eta_x*( -6.0 + 4.0*eta*eta ) )
      return
      end
c
      function dem2dx2ms(x,y)
      sigeta=1.e1
      eta=sigeta*y/x
      et2=eta*eta
      t1=2.0*em2ms(x,y)/x/x
      dem2dx2ms=t1*(et2*(2.0*et2-7.0)+3.0)
      return
      end
c
      function dem4dx2ms(x,y)
      sigeta=1.e1
      eta=sigeta*y/x
      et2=eta*eta
      t1=2.0*em4ms(x,y)/x/x
      dem4dx2ms=t1*(et2*(2.0*et2-11.0)+10.0)
      return
      end
c
      function dem1dy2ms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      t0=sqrt(2.0)*exp(0.5)*emmax
      eta    = sigeta*y/x
      eta_y  = sigeta/x
      eta_yy = 0.0
      dem1dy2ms=t0*exp(-eta*eta)*(        eta_yy*(  1.0 - 2.0*eta*eta )
     +                       + eta*eta_y*eta_y*( -6.0 + 4.0*eta*eta ) )
      return
      end
c
      function dem2dy2ms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      eta=sigeta*y/x
      et2=eta*eta
      t1=2.0*exp(1.0)*emmax*(sigeta/x)*(sigeta/x)*exp(-et2)
      dem2dy2ms=t1*(et2*(2.0*et2-5.0)+1.0)
      return
      end
c
      function dem4dy2ms(x,y)
      sigeta=1.e1
      emmax=1.e-3
      eta=sigeta*y/x
      et2=eta*eta
      t1=0.5*exp(2.)*emmax*et2*(sigeta/x)*(sigeta/x)*exp(-et2)
      dem4dy2ms=t1*(et2*(2.0*et2-9.0)+6.0)
      return
      end
c
      function ssams(x,y)
      emu=1.e-6
c  spalart & allmaras one-equation model
      yloc=max(1.e-40,y)
      yloc=yloc*yloc
      u1  =ums(x,y)
      u2  =vms(x,y)
      emsa=em2ms(x,y)
c  first and second derivatives of mu_tilde
      desadx  =dem2dxms(x,y)
      desady  =dem2dyms(x,y)
      desadx2 =dem2dx2ms(x,y)
      desady2 =dem2dy2ms(x,y)
c  convection
      cv1=u1*desadx+u2*desady
c  diffusion
      sig=2./3.
      df1=(emu+emsa)*(desadx2+desady2)
      df1=df1+1.622*(desadx*desadx+desady*desady)
      df1=-df1/sig
!  production
      emt=emsa/emu
      fv1=emsa**3/(emsa**3+(7.1*emu)**3)
      fv2=1.0-emt/(emt*fv1+1.0)
      dyp=emsa/0.41/0.41/yloc
      str=vortms(x,y)+dyp*fv2
      cp1=-0.1355*str*emsa
!  dissipation
      cw1=0.1355/0.41/0.41+1.622/sig
      r=dyp/str
      g=r*(1.0+0.3*(r**5-1.0))
      fw=65.0/(g**6+64.0)
      fw=g*(fw**(1.0/6.0))
      cd1=cw1*fw*emsa*emsa/yloc
      ssams=cv1+df1+cp1+cd1
      return
      end
c
      function eddysams(x,y)
      emsa=em2ms(x,y)
      eddysams=fv1sams(emsa)*emsa
      return
      end
c
      function desadxms(x,y)
      emsa=em2ms(x,y)
      desadxms=(fv1sams(emsa)+dfv1sams(emsa)*emsa)*dem2dxms(x,y)
      return
      end
c
      function desadyms(x,y)
      emsa=em2ms(x,y)
      desadyms=(fv1sams(emsa)+dfv1sams(emsa)*emsa)*dem2dyms(x,y)
      return
      end
c
      function fv1sams(em)
      emu=1.e-6
      cv1=(7.1*emu)**3
      fv1sams=em*em*em/(em*em*em+cv1)
      return
      end
c
      function dfv1sams(em)
      emu=1.e-6
      cv1=(7.1*emu)**3
      dfv1sams=3.0*em*em*cv1/(em*em*em+cv1)/(em*em*em+cv1)
      return
      end
c
      function smtms(x,y)
c menter one-equation turbulence model
      emu=1.e-6
      u1  =ums(x,y)
      u2  =vms(x,y)
      emmt=em2ms(x,y)
      emtu=eddymtms(x,y)
c  first and second derivatives of mu_tilde
      demtdx  =dem2dxms(x,y)
      demtdy  =dem2dyms(x,y)
      demtdx2 =dem2dx2ms(x,y)
      demtdy2 =dem2dy2ms(x,y)
c  convection
      cv1=u1*demtdx+u2*demtdy
c  diffusion
      sig=1.0
      df1=-(emu+emmt/sig)*(demtdx2+demtdy2)
      df1=df1-(demtdx*demtdx+demtdy*demtdy)/sig
c  production
      fd1=(emtu+emu)/(emmt+emu)
      str=strainms(x,y)
      cp1=-0.144*fd1*str*emmt
c  dissipation
      u1dx  =dudxms(x,y)
      u1dy  =dudyms(x,y)
      u1dx2 =dudx2ms(x,y)
      u1dy2 =dudy2ms(x,y)
      u1dxy2=dudxyms(x,y)
      u2dx  =dvdxms(x,y)
      u2dy  =dvdyms(x,y)
      u2dx2 =dvdx2ms(x,y)
      u2dy2 =dvdy2ms(x,y)
      u2dxy2=dvdxyms(x,y)
c
      dsdx=4.0*(u1dx*u1dx2+u2dy*u2dxy2)+2.0*(u1dy+u2dx)*(u1dxy2+u2dx2)
      dsdy=4.0*(u1dx*u1dxy2+u2dy*u2dy2)+2.0*(u1dy+u2dx)*(u1dy2+u2dxy2)
      dsdt=(dsdx*dsdx+dsdy*dsdy)/(4.0*str*str)
      ekep=emmt*emmt*dsdt/str/str
      epbb=demtdx*demtdx+demtdy*demtdy
      ep1e=7.0*epbb*tanh(ekep/(7.0*epbb))
      cd1 =1.862*ep1e
c
      smtms=cv1+df1+cp1+cd1
      return
      end
c
      function eddymtms(x,y)
      emmt=em2ms(x,y)
      eddymtms=d2mtms(emmt)*emmt
      return
      end
c
      function demtdxms(x,y)
      emmt=em2ms(x,y)
      demtdxms=(d2mtms(emmt)+dd2mtms(emmt)*emmt)*dem2dxms(x,y)
      return
      end
c
      function demtdyms(x,y)
      emmt=em2ms(x,y)
      demtdyms=(d2mtms(emmt)+dd2mtms(emmt)*emmt)*dem2dyms(x,y)
      return
      end
c
      function d2mtms(em)
      emu=1.e-6
      cv1 =13.0*0.41*emu
      cdum=(em/cv1)**2
      d2mtms=1.0-exp(-cdum)
      return
      end
c
      function dd2mtms(em)
      emu=1.e-6
      cv1    =13.0*0.41*emu
      cdum   =em/cv1
      dd2mtms=2.0*cdum*exp(-cdum*cdum)/cv1
      return
      end
c
      function user_derfc(x)
      parameter(pv= 1.26974899965115684e+01, 
     +          ph= 6.10399733098688199e+00, 
     +          p0= 2.96316885199227378e-01, 
     +          p1= 1.81581125134637070e-01, 
     +          p2= 6.81866451424939493e-02, 
     +          p3= 1.56907543161966709e-02, 
     +          p4= 2.21290116681517573e-03, 
     +          p5= 1.91395813098742864e-04, 
     +          p6= 9.71013284010551623e-06, 
     +          p7= 1.66642447174307753e-07)
      parameter(q0= 6.12158644495538758e-02, 
     +          q1= 5.50942780056002085e-01, 
     +          q2= 1.53039662058770397e+00, 
     +          q3= 2.99957952311300634e+00, 
     +          q4= 4.95867777128246701e+00, 
     +          q5= 7.41471251099335407e+00, 
     +          q6= 1.04765104356545238e+01, 
     +          q7= 1.48455557345597957e+01)
      y=x*x
      y=exp(-y)*x*(p7/(y+q7)+p6/(y+q6) 
     +  +p5/(y+q5)+p4/(y+q4)+p3/(y+q3) 
     +  +p2/(y+q2)+p1/(y+q1)+p0/(y+q0))
      if(x < ph) y=y+2/(exp(pv*x)+1)
      user_derfc=y
      return
      end
c
      function user_derf(x)
      parameter(p0= 1.12837916709551257e+00, 
     +          p1=-3.76126389031833602e-01, 
     +          p2= 1.12837916706621301e-01, 
     +          p3=-2.68661698447642378e-02, 
     +          p4= 5.22387877685618101e-03, 
     +          p5=-8.49202435186918470e-04)
      y=abs(x)
      if(y .gt. 0.125) then
        user_derf=sign(1-user_derfc(y),x)
      else
        y=x*x
        user_derf=(((((p5*y+p4)*y+p3)*y+p2)*y+p1)*y+p0)*x
      end if
      return
      end
c
      subroutine exact_polyfg(scale, polyfc, x, y, z, q, qx, qy, qz,
     +                    qxx, qyy, qzz, qxy, qxz, qyz )
c=============================== EXACT_POLYFG ================================80
c
c Add exact polynomial solution and gradients.
c
c=============================================================================80
      dimension polyfc(13)
c
      c     = polyfc( 1)
      cx    = polyfc( 2)
      cy    = polyfc( 3)
      cz    = polyfc( 4)
      cxx   = polyfc( 5)
      cyy   = polyfc( 6)
      czz   = polyfc( 7)
      cxxx  = polyfc( 8)
      cyyy  = polyfc( 9)
      czzz  = polyfc(10)
      cxxxx = polyfc(11)
      cyyyy = polyfc(12)
      czzzz = polyfc(13)

      q = q + scale*( c                                     
     +      + cx*x + 0.5*cxx*(x**2) + 1./6.*cxxx*(x**3) 
     +      + cy*y + 0.5*cyy*(y**2) + 1./6.*cyyy*(y**3) 
     +      + cz*z + 0.5*czz*(z**2) + 1./6.*czzz*(z**3) 
     +      + 1./12.*cxxxx*(x**4)                           
     +      + 1./12.*cyyyy*(y**4)                           
     +      + 1./12.*czzzz*(z**4) )

      qx =  qx + scale*( cx + cxx*x + 0.5*cxxx*(x**2) + 
     +   1./3.*cxxxx*(x**3) )
      qy =  qy + scale*( cy + cyy*y + 0.5*cyyy*(y**2) + 
     +   1./3.*cyyyy*(y**3) )
      qz =  qz + scale*( cz + czz*z + 0.5*czzz*(z**2) + 
     +   1./3.*czzzz*(z**3) )

      qxx = qxx + scale*( cxx + cxxx*x + cxxxx*(x**2) )
      qyy = qyy + scale*( cyy + cyyy*y + cyyyy*(y**2) )
      qzz = qzz + scale*( czz + czzz*z + czzzz*(z**2) )

      qxy = qxy + scale*( 0.0 )
      qxz = qxz + scale*( 0.0 )
      qyz = qyz + scale*( 0.0 )
      return
      end
#   endif
