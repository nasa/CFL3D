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
      subroutine initvist(nbl,jdim,kdim,idim,vist3d,tursav,smin,cmuv,
     .                    nummem,x,y,z)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Initialize the turbulent initial conditions on a mesh
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension tursav(jdim,kdim,idim,nummem),vist3d(jdim,kdim,idim),
     .          smin(jdim-1,kdim-1,idim-1),cmuv(jdim-1,kdim-1,idim-1)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /reyue/ reue,tinf,ivisc(3)
      common /mms/ iexact_trunc,iexact_disc,iexact_ring
      common /initfac/ scal_ic
      common /reystressmodel/ issglrrw2012,i_sas_rsm,i_yapterm
c
c   scaling factor for distance-based I.C.s:
      if (scal_ic .eq. 0.) then
c       Use freestream ICs everywhere
        iscal=0
      else
c       Scale turb-like IC in BL; lower scal_ic yields thinner layer
        refac=reue/scal_ic
        iscal=1
      end if
c
c     write(15,904)nbl
  904 format(39h turbulent initial conditions for block,i3)
c
      if (ivisc(1).eq.8  .or. ivisc(2).eq.8  .or. ivisc(3).eq.8  .or.
     +    ivisc(1).eq.9  .or. ivisc(2).eq.9  .or. ivisc(3).eq.9  .or.
     +    ivisc(1).eq.13 .or. ivisc(2).eq.13 .or. ivisc(3).eq.13 .or.
     +    ivisc(1).eq.14 .or. ivisc(2).eq.14 .or. ivisc(3).eq.14) then
      do 4020 i=1,idim-1
        do 4020 k=1,kdim-1
          do 4020 j=1,jdim-1
            cmuv(j,k,i)=-.09
 4020 continue
      end if
c   For field eqn turbulence models:
      if (ivisc(1).ge.4 .or. ivisc(2).ge.4 .or. ivisc(3).ge.4) then
        if (ivisc(1).eq.4 .or. ivisc(2).eq.4 .or. ivisc(3).eq.4) then
          do 2000 i=1,idim-1
            do 2000 k=1,kdim-1
              do 2000 j=1,jdim-1
                tursav(j,k,i,1)=tur10(1)
                vist3d(j,k,i)=0.09*tur10(1)
 2000     continue
        else if (ivisc(1).eq.5 .or. ivisc(2).eq.5 .or.
     .           ivisc(3).eq.5) then
          do 2001 i=1,idim-1
            do 2001 k=1,kdim-1
              do 2001 j=1,jdim-1
                tursav(j,k,i,1)=tur10(1)
                vist3d(j,k,i)=tur10(1)*(tur10(1)**3/
     .                        (tur10(1)**3+357.911))
 2001     continue
        else if (ivisc(1).eq.16 .or. ivisc(2).eq.16 .or.
     .           ivisc(3).eq.16) then
          do 2011 i=1,idim-1
            do 2011 k=1,kdim-1
              do 2011 j=1,jdim-1
                tursav(j,k,i,1)=tur10(1)
                tursav(j,k,i,2)=tur10(2)
                vist3d(j,k,i)=0.54772*tur10(1)/sqrt(tur10(2))
 2011     continue
        else if (ivisc(1).eq.11 .or. ivisc(2).eq.11 .or.
     .           ivisc(3).eq.11 .or.
     .           ivisc(1).eq.10 .or. ivisc(2).eq.10 .or.
     .           ivisc(3).eq.10 .or.
     .           ivisc(1).eq. 9 .or. ivisc(2).eq. 9 .or.
     .           ivisc(3).eq. 9 .or. ivisc(1).eq.13 .or.
     .           ivisc(2).eq.13 .or. ivisc(3).eq.13 .or.
     .           ivisc(1).eq.15 .or. ivisc(2).eq.15 .or.
     .           ivisc(3).eq.15) then
          do 2002 i=1,idim-1
            do 2002 k=1,kdim-1
              do 2002 j=1,jdim-1
                zk1=tur10(2)
                exponent=-471.*(ccabs(smin(j,k,i))*refac)+.47
                if(real(exponent) .lt. -40.) then
                  zk2=0.
                else
                  zk2=10.**(exponent)
                end if
                exponent=-37.5*(ccabs(smin(j,k,i))*refac)-3.7
                if(real(exponent) .lt. -40.) then
                  zk3=0.
                else
                  zk3=10.**(exponent)
                end if
                zk4=6.7*(ccabs(smin(j,k,i))*refac)
                tursav(j,k,i,2)=ccmin(zk2,zk3)
                tursav(j,k,i,2)=ccmax(tursav(j,k,i,2),zk1)
                tursav(j,k,i,2)=ccmin(tursav(j,k,i,2),zk4)
                tursav(j,k,i,2)=tursav(j,k,i,2)*iscal +
     +             tur10(2)*(1-iscal)
c
                ep1=tur10(1)
                exponent=-555.*(ccabs(smin(j,k,i))*refac)-6.
                if(real(exponent) .lt. -20.) then
                  ep2=0.
                else
                  ep2=10.**(exponent)
                end if
                exponent=-280.*(ccabs(smin(j,k,i))*refac)-9.2
                if(real(exponent) .lt. -20.) then
                  ep3=0.
                else
                  ep3=10.**(exponent)
                end if
                exponent=13333.*(ccabs(smin(j,k,i))*refac)-9.8
                if(real(exponent) .gt. 20.) then
                  ep4=1.e20
                else
                  ep4=10.**(exponent)
                end if
                tursav(j,k,i,1)=ccmin(ep2,ep3)
                tursav(j,k,i,1)=ccmax(tursav(j,k,i,1),ep1)
                tursav(j,k,i,1)=ccmin(tursav(j,k,i,1),ep4)
                tursav(j,k,i,1)=tursav(j,k,i,1)*iscal +
     +             tur10(1)*(1-iscal)
                vist3d(j,k,i)=.09*rho0*tursav(j,k,i,2)**2/
     .                        tursav(j,k,i,1)

 2002     continue
        else if (ivisc(1).eq.25 .or. ivisc(2).eq.25 .or.
     .           ivisc(3).eq.25) then
c         do nothing
          do i=1,idim-1
            do k=1,kdim-1
              do j=1,jdim-1
                tursav(j,k,i,1)=0.0
                vist3d(j,k,i)=0.0
              enddo
            enddo
          enddo
          continue
        elseif(ivisc(1).eq.72.or.ivisc(2).eq.72.or.ivisc(3).eq.72) then
          if (issglrrw2012 .eq. 6) then
          do i=1,idim-1
            do k=1,kdim-1
              do j=1,jdim-1
                tursav(j,k,i,1)=tur10(1)
                tursav(j,k,i,2)=tur10(2)
                tursav(j,k,i,3)=tur10(3)
                tursav(j,k,i,4)=tur10(4)
                tursav(j,k,i,5)=tur10(5)
                tursav(j,k,i,6)=tur10(6)
                tursav(j,k,i,7)=tur10(7)
                zkinf=-(tur10(1)+tur10(2)+tur10(3))/2.
                vist3d(j,k,i)=rho0*zkinf*(tur10(7)**2)
              enddo
            enddo
          enddo
          else
          do i=1,idim-1
            do k=1,kdim-1
              do j=1,jdim-1
                tursav(j,k,i,1)=tur10(1)
                tursav(j,k,i,2)=tur10(2)
                tursav(j,k,i,3)=tur10(3)
                tursav(j,k,i,4)=tur10(4)
                tursav(j,k,i,5)=tur10(5)
                tursav(j,k,i,6)=tur10(6)
                tursav(j,k,i,7)=tur10(7)
                zkinf=-(tur10(1)+tur10(2)+tur10(3))/2.
                vist3d(j,k,i)=rho0*zkinf/tur10(7)
              enddo
            enddo
          enddo
          end if
        else
          const1=45.8
          const2=1.68
          smax=const2/(2.*const1)
          tmax=-const1*smax**2 + const2*smax
          cmu=1.
          if (ivisc(1).eq. 8 .or. ivisc(2).eq. 8 .or.
     .        ivisc(3).eq. 8 .or.
     .        ivisc(1).eq.12 .or. ivisc(2).eq.12 .or.
     .        ivisc(3).eq.12 .or.
     .        ivisc(1).eq.14 .or. ivisc(2).eq.14 .or.
     .        ivisc(3).eq.14) cmu=.09
          v3dset=cmu*rho0*tur10(2)/tur10(1)
          do 2009 i=1,idim-1
            do 2009 k=1,kdim-1
              do 2009 j=1,jdim-1
                zk1=tur10(2)
                zk2=-const1*(ccabs(smin(j,k,i))*refac)**2 +
     +               const2*ccabs(smin(j,k,i))*refac
                tursav(j,k,i,2)=ccmax(zk1,zk2)
                tursav(j,k,i,2)=tursav(j,k,i,2)*iscal +
     +             tur10(2)*(1-iscal)
                v3d=tursav(j,k,i,2)*100./tmax
                v3d=ccmax(v3d,v3dset)
                om1=-12444.*(ccabs(smin(j,k,i))*refac) + .54
                om2=cmu*tursav(j,k,i,2)/v3d
                tursav(j,k,i,1)=ccmax(om1,om2)
                tursav(j,k,i,1)=tursav(j,k,i,1)*iscal +
     +             tur10(1)*(1-iscal)
                vist3d(j,k,i)=cmu*rho0*tursav(j,k,i,2)/
     .                        tursav(j,k,i,1)
 2009     continue
        end if
c
c       Special for 3-eqn transition model ivisc=30
        if (ivisc(1).eq.30.or. ivisc(2).eq.30.or.
     .           ivisc(3).eq.30) then
          do i=1,idim-1
            do k=1,kdim-1
              do j=1,jdim-1
                tursav(j,k,i,3)=tur10(3)
              enddo
            enddo
          enddo
        end if
c
c       Special for 4-eqn transition model ivisc=40
        if (ivisc(1).eq.40.or. ivisc(2).eq.40.or.
     .           ivisc(3).eq.40) then
          do i=1,idim-1
            do k=1,kdim-1
              do j=1,jdim-1
                tursav(j,k,i,3)=tur10(3)
                tursav(j,k,i,4)=tur10(4)
              enddo
            enddo
          enddo
        end if
c
      else
c
      do 4021 i=1,idim-1
        do 4021 k=1,kdim-1
          do 4021 j=1,jdim-1
            vist3d(j,k,i)=0.
 4021 continue
c
      end if
c
c   Overwrite with exact soln if doing MMS
      if (iexact_trunc .ne. 0 .or. iexact_disc .ne. 0) then
        call exact_turb_q(jdim,kdim,idim,x,y,z,tursav,smin,vist3d,
     +      iexact_trunc,iexact_disc)
      end if
c
      return
      end
