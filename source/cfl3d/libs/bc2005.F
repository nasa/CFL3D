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
      subroutine bc2005(jdim,kdim,idim,q,qj0,qk0,qi0,
     .                  ista,iend,jsta,jend,ksta,kend,nface,
     .                  tursav,tj0,tk0,ti0,vist3d,vj0,vk0,vi0,
     .                  mdim,ndim,bcdata,filname,qp,vp,tp,
     .                  jdimp,kdimp,idimp,qrotj,qrotk,qroti,nbl,
     .                  nblp,nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl,
     .                  nummem)
#   ifdef CMPLX
#   else
      use module_stm_2005, only: stm2k5_get_rotmat, stm2k5_bc
#   endif
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set periodic boundary conditions, given angular rotation 
c               angle to the periodic face and its block number
c
c     blnb   = block number of periodic face
c     dthtx  = angle about line parallel with x-axis
c     dthty  = angle about line parallel with y-axis
c     dthtz  = angle about line parallel with z-axis
c       (NOTE:  only one of the 3 angles can be used at a time; i.e.,
c               2 of them must be identically zero)
c       The angles should be measured using the right-hand rule; i.e.,
c       point your right thumb in the direction of the +N axis (for rotation
c       about N-axis, where N = x, y, or z) - 
c       the direction of finger curl is the direction
c       of positive angle.  If you are setting the angle for a particular
c       face (e.g., the i0 face), set the angle = the angle you 
c       have to move the PERIODIC FACE through to get to this face.
c
c     NOTE:  Currently, it is assumed that the current block and the
c            block it is periodic with are 1-to-1 at the corresponding
c            faces after rotation through the specified angle dthtN.
c            Also, the 2 blocks are assumed to be aligned similarly.
c            In other  words, i,j, and k must each run in the same
c            directions, respectively.  Therefore, if the periodic BC 
c            is being applied (on the current block) on the KMAX face, it 
c            is implicitly assumed that the corresponding surface it is 
c            periodic with is KMIN, with i and j running in the same 
c            directions as on the current block.  Also, therefore, the 2 
c            remaining dimensions on the periodic face (2 of idim, jdim, kdim)
c            of the current block must be identical to the same 2 dimensions
c            of the periodic block.
c            
c            This periodic BC also works for a 1-cell-in-the-periodic-dimension
c            grid that is periodic with itself.  Note, however, that if a
c            particular block is periodic with a DIFFERENT block, then
c            neither block should be only 1-cell wide.
c***********************************************************************
c     Description of variables:
c       jdim,kdim,idim    = dimensions of current block
c       q                 = primitive variables on current block
c       qj0,qk0,qi0       = BC values assigned for this block
c       ista,iend,etc.    = indices over which BC is applied
c       nface             = face number (1 = i=1  2 = i=idim
c                                        3 = j=1  4 = j=jdim
c                                        5 = k=1  6 = k=kdim)
c       tursav            = turb quantities on current block
c       tj0,tk0,ti0       = turb BC values assigned for this block
c       vist3d            = eddy viscosity on current block
c       vj0,vk0,vi0       = eddy viscosity BC values assigned for this block
c       mdim,ndim         = dimensions of bcdata
c       bcdata            = auxiliary data that goes with this BC
c       filename          = filename to read bcdata, if array values
c       qp                = primitive variables on periodic block
c       vp                = eddy viscosity on periodic block
c       tp                = turb quantities on periodic block
c       jdimp,kdimp,idimp = dimensions of periodic block
c       qrotj,qrotk,qroti = work arrays of q values on periodic block
c       nbl               = block number of current block
c       nblp              = block number of periodic block
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension qp(jdimp,kdimp,idimp,5),qrotj(2,kdimp,idimp,5),
     .          qrotk(jdimp,2,idimp,5),qroti(jdimp,kdimp,2,5)
      dimension vp(jdimp,kdimp,idimp),tp(jdimp,kdimp,idimp,nummem)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcdata(mdim,ndim,2,12)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
      dimension mblk2nd(maxbl)
c
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /conversion/ radtodeg
      common /igrdtyp/ ip3dgrd,ialph
      real :: rn(3,3), rnt(3,3)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c     this bc makes use of only one plane of data    
c
      ip    = 1
      dthtx = bcdata(1,1,ip,2)/radtodeg
      if (ialph.eq.0) then
         dthty = bcdata(1,1,ip,3)/radtodeg
         dthtz = bcdata(1,1,ip,4)/radtodeg
      else
         dthty = -bcdata(1,1,ip,4)/radtodeg
         dthtz =  bcdata(1,1,ip,3)/radtodeg
      end if
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=2005 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary             periodic boundary                   bctype 2005
c******************************************************************************
c
      if (nface.eq.3) then
c
c   Load qp values into qrotX array:
      do 2389 l=1,5
      do 2389 i=ista,iend1
      do 2389 k=ksta,kend1
        qrotj(1,k,i,l)=qp(jdimp-1,k,i,l)
        jmx=max(jdimp-2,1)
        qrotj(2,k,i,l)=qp(jmx,k,i,l)
 2389 continue
c   Rotate qrotX values
      if (jdimp .eq. 2) then
        call rotateq(2,kdimp,idimp,qrotj,qrotj,ista,iend1,1,1,
     .               ksta,kend1,dthtx,dthty,dthtz)
        call rotateq(2,kdimp,idimp,qrotj,qrotj,ista,iend1,2,2,
     .               ksta,kend1,2.*dthtx,2.*dthty,2.*dthtz)
      else
        call rotateq(2,kdimp,idimp,qrotj,qrotj,ista,iend1,1,2,
     .               ksta,kend1,dthtx,dthty,dthtz)
      end if
c
c   Apply periodic BCs
      do 100 l=1,5
      do 100 i=ista,iend1
      do 100 k=ksta,kend1
        qj0(k,i,l,1) = qrotj(1,k,i,l)
        qj0(k,i,l,2) = qrotj(2,k,i,l)
 100  continue
c
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
          vj0(k,i,1,1) = vp(jdimp-1,k,i)
          vj0(k,i,1,2) = 0.
  191   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.70.or. ivisc(2).ge.70.or. ivisc(1).ge.70) then
         
#   ifdef CMPLX
#   else
         call stm2k5_get_rotmat(dthtx, dthty, dthtz, rn,rnt)

         ! rotate the stress tensors
         do i=ista,iend1
            do k=ksta,kend1
               
               call stm2k5_bc(tp(jdimp-1,k,i,1:7), rn, rnt,
     $              tj0(k,i,1:7,1))
               if(jdimp==2) then
                  CALL stm2k5_bc(tj0(k,i,1:7,1), rn, rnt,
     $                 tj0(k,i,1:7,2))
               else
                  call stm2k5_bc(tp(jdimp-2,k,i,1:7), rn, rnt, 
     $                 tj0(k,i,1:7,2))

               endif
            enddo
         enddo
#   endif
      elseif (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        jdm2=max(jdimp-2,1)
        do l=1,nummem
        do 101 i=ista,iend1
        do 101 k=ksta,kend1
             tj0(k,i,l,1) = tp(jdimp-1,k,i,l)
             tj0(k,i,l,2) = tp(jdm2,k,i,l)
  101   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      j=jdim boundary          periodic boundary                   bctype 2005
c******************************************************************************
c
      if (nface.eq.4) then
c
c   Load qp values into qrotX array:
      do 3389 l=1,5
      do 3389 i=ista,iend1
      do 3389 k=ksta,kend1
        qrotj(1,k,i,l)=qp(1,k,i,l)
        jmx=min(jdimp-1,2)
        qrotj(2,k,i,l)=qp(jmx,k,i,l)
 3389 continue
c   Rotate qrotX values
      if (jdimp .eq. 2) then
        call rotateq(2,kdimp,idimp,qrotj,qrotj,ista,iend1,1,1,
     .               ksta,kend1,dthtx,dthty,dthtz)
        call rotateq(2,kdimp,idimp,qrotj,qrotj,ista,iend1,2,2,
     .               ksta,kend1,2.*dthtx,2.*dthty,2.*dthtz)
      else
        call rotateq(2,kdimp,idimp,qrotj,qrotj,ista,iend1,1,2,
     .               ksta,kend1,dthtx,dthty,dthtz)
      end if
c
c   Apply periodic BCs
      do 200 l=1,5
      do 200 i=ista,iend1
      do 200 k=ksta,kend1
        qj0(k,i,l,3) = qrotj(1,k,i,l)
        qj0(k,i,l,4) = qrotj(2,k,i,l)
 200  continue
c
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
          vj0(k,i,1,3) = vp(1,k,i)
          vj0(k,i,1,4) = 0.0
  291   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.70.or. ivisc(2).ge.70.or. ivisc(1).ge.70) then

#   ifdef CMPLX
#   else
         ! obtain the rotation matrix and its transpose
         call stm2k5_get_rotmat(dthtx, dthty, dthtz, rn,rnt)

         ! rotate the stress tensors
         do i=ista,iend1
            do k=ksta,kend1

               call stm2k5_bc(tp(1,k,i,1:7), rn, rnt, 
     $              tj0(k,i,1:7,3))
               if(jdimp==2) then
                 call stm2k5_bc(tj0(k,i,1:7,3), rn, rnt,
     $                tj0(k,i,1:7,4))
               else
                 call stm2k5_bc(tp(2,k,i,1:7), rn, rnt,
     $                tj0(k,i,1:7,4))
               endif
            enddo
         enddo
#   endif
      elseif (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        jdp2=min(2,jdimp-1)
        do l=1,nummem
        do 201 i=ista,iend1
        do 201 k=ksta,kend1
             tj0(k,i,l,3) = tp(1,k,i,l)
             tj0(k,i,l,4) = tp(jdp2,k,i,l)
  201   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=1 boundary             periodic boundary                   bctype 2005
c******************************************************************************
c
      if (nface.eq.5) then
c
c   Load qp values into qrotX array:
      do 4389 l=1,5
      do 4389 i=ista,iend1
      do 4389 j=jsta,jend1
        qrotk(j,1,i,l)=qp(j,kdimp-1,i,l)
        kmx=max(kdimp-2,1)
        qrotk(j,2,i,l)=qp(j,kmx,i,l)
 4389 continue
c   Rotate qrotX values
      if (kdimp .eq. 2) then
        call rotateq(jdimp,2,idimp,qrotk,qrotk,ista,iend1,jsta,jend1,
     .               1,1,dthtx,dthty,dthtz)
        call rotateq(jdimp,2,idimp,qrotk,qrotk,ista,iend1,jsta,jend1,
     .               2,2,2.*dthtx,2.*dthty,2.*dthtz)
      else
        call rotateq(jdimp,2,idimp,qrotk,qrotk,ista,iend1,jsta,jend1,
     .               1,2,dthtx,dthty,dthtz)
      end if
c
c   Apply periodic BCs
      do 300 l=1,5
      do 300 i=ista,iend1
      do 300 j=jsta,jend1
        qk0(j,i,l,1) = qrotk(j,1,i,l)
        qk0(j,i,l,2) = qrotk(j,2,i,l)
 300  continue
c
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
          vk0(j,i,1,1) = vp(j,kdimp-1,i)
          vk0(j,i,1,2) = 0.0
  391   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.70.or. ivisc(2).ge.70.or. ivisc(1).ge.70) then

#   ifdef CMPLX
#   else
         ! obtain the rotation matrix and its transpose
         call stm2k5_get_rotmat(dthtx, dthty, dthtz, rn,rnt)

         ! rotate the stress tensors
         do i=ista,iend1
            do j=jsta,jend1
               call stm2k5_bc(tp(j,kdimp-1,i,1:7), rn, rnt, 
     $              tk0(j,i,1:7,1))
               if(kdimp==2) then
                  call stm2k5_bc(tk0(j,i,1:7,1), rn, rnt, 
     $                 tk0(j,i,1:7,2))
            
               else
                  call stm2k5_bc(tp(j,kdimp-2,i,1:7), rn, rnt, 
     $                 tk0(j,i,1:7,2))
               endif

            enddo
         enddo
#   endif
      elseif (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        kdm2=max(kdimp-2,1)
        do l=1,nummem
        do 301 i=ista,iend1
        do 301 j=jsta,jend1
             tk0(j,i,l,1) = tp(j,kdimp-1,i,l)
             tk0(j,i,l,2) = tp(j,kdm2,i,l)
  301   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=kdim boundary          periodic boundary                   bctype 2005
c******************************************************************************
c
      if (nface.eq.6) then
c
c   Load qp values into qrotX array:
      do 5389 l=1,5
      do 5389 i=ista,iend1
      do 5389 j=jsta,jend1
        qrotk(j,1,i,l)=qp(j,1,i,l)
        kmx=min(kdimp-1,2)
        qrotk(j,2,i,l)=qp(j,kmx,i,l)
 5389 continue
c   Rotate qrotX values
      if (kdimp .eq. 2) then
        call rotateq(jdimp,2,idimp,qrotk,qrotk,ista,iend1,jsta,jend1,
     .               1,1,dthtx,dthty,dthtz)
        call rotateq(jdimp,2,idimp,qrotk,qrotk,ista,iend1,jsta,jend1,
     .               2,2,2.*dthtx,2.*dthty,2.*dthtz)
      else
        call rotateq(jdimp,2,idimp,qrotk,qrotk,ista,iend1,jsta,jend1,
     .               1,2,dthtx,dthty,dthtz)
      end if
c
c   Apply periodic BCs
      do 400 l=1,5
      do 400 i=ista,iend1
      do 400 j=jsta,jend1
        qk0(j,i,l,3) = qrotk(j,1,i,l)
        qk0(j,i,l,4) = qrotk(j,2,i,l)
 400  continue
c
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
          vk0(j,i,1,3) = vp(j,1,i)
          vk0(j,i,1,4) = 0.0
  491   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.70.or. ivisc(2).ge.70.or. ivisc(1).ge.70) then

#   ifdef CMPLX
#   else
         ! obtain the rotation matrix and its transpose
         call stm2k5_get_rotmat(dthtx, dthty, dthtz, rn,rnt)

         ! rotate the stress tensors
         do i=ista,iend1
            do j=jsta,jend1

               call stm2k5_bc(tp(j,1,i,1:7), rn, rnt, 
     $              tk0(j,i,1:7,3))
               if(kdimp==2) then
                  call stm2k5_bc(tk0(j,i,1:7,3), rn, rnt, 
     $                 tk0(j,i,1:7,4))
               else
                  call stm2k5_bc(tp(j,2,i,1:7), rn, rnt, 
     $                 tk0(j,i,1:7,4))
               endif

            enddo
         enddo
#   endif
      elseif (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        kdp2=min(2,kdimp-1)
        do l=1,nummem
        do 401 i=ista,iend1
        do 401 j=jsta,jend1
             tk0(j,i,l,3) = tp(j,1,i,l)
             tk0(j,i,l,4) = tp(j,kdp2,i,l)
  401   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=1 boundary             periodic boundary                   bctype 2005
c******************************************************************************
c
      if (nface.eq.1) then
c
c   Load qp values into qrotX array:
      do 6389 l=1,5
      do 6389 k=ksta,kend1
      do 6389 j=jsta,jend1
        qroti(j,k,1,l)=qp(j,k,idimp-1,l)
        imx=max(idimp-2,1)
        qroti(j,k,2,l)=qp(j,k,imx,l)
 6389 continue
c   Rotate qrotX values
      if (idimp .eq. 2) then
        call rotateq(jdimp,kdimp,2,qroti,qroti,1,1,jsta,jend1,
     .               ksta,kend1,dthtx,dthty,dthtz)
        call rotateq(jdimp,kdimp,2,qroti,qroti,2,2,jsta,jend1,
     .               ksta,kend1,2.*dthtx,2.*dthty,2.*dthtz)
      else
        call rotateq(jdimp,kdimp,2,qroti,qroti,1,2,jsta,jend1,
     .               ksta,kend1,dthtx,dthty,dthtz)
      end if
c
c   Apply periodic BCs
      do 500 l=1,5
      do 500 k=ksta,kend1
      do 500 j=jsta,jend1
        qi0(j,k,l,1) = qroti(j,k,1,l)
        qi0(j,k,l,2) = qroti(j,k,2,l)
 500  continue
c
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
          vi0(j,k,1,1) = vp(j,k,idimp-1)
          vi0(j,k,1,2) = 0.0
  591   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.70.or. ivisc(2).ge.70.or. ivisc(1).ge.70) then

#   ifdef CMPLX
#   else
         ! obtain the rotation matrix and its transpose
         call stm2k5_get_rotmat(dthtx, dthty, dthtz, rn,rnt)

         ! rotate the stress tensors
         do k=ksta,kend1
            do j=jsta,jend1

               call stm2k5_bc(tp(j,k,idimp-1,1:7), rn, rnt,
     $              ti0(j,k,1:7,1))
               if(idimp==2) then
                  call stm2k5_bc(ti0(j,k,1:7,1), rn, rnt,
     $                 ti0(j,k,1:7,2))
               else      
                  call stm2k5_bc(tp(j,k,idimp-2,1:7), rn, rnt,
     $                 ti0(j,k,1:7,2))
               endif

            enddo
         enddo
#   endif
      elseif (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        idm2=max(idimp-2,1)
        do l=1,nummem
        do 501 k=ksta,kend1
        do 501 j=jsta,jend1
             ti0(j,k,l,1) = tp(j,k,idimp-1,l)
             ti0(j,k,l,2) = tp(j,k,idm2,l)
  501   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=idim boundary          periodic boundary                   bctype 2005
c******************************************************************************
c
      if (nface.eq.2) then
c
c   Load qp values into qrotX array:
      do 7389 l=1,5
      do 7389 k=ksta,kend1
      do 7389 j=jsta,jend1
        qroti(j,k,1,l)=qp(j,k,1,l)
        imx=min(idimp-1,2)
        qroti(j,k,2,l)=qp(j,k,imx,l)
 7389 continue
c   Rotate qrotX values
      if (idimp .eq. 2) then
        call rotateq(jdimp,kdimp,2,qroti,qroti,1,1,jsta,jend1,
     .               ksta,kend1,dthtx,dthty,dthtz)
        call rotateq(jdimp,kdimp,2,qroti,qroti,2,2,jsta,jend1,
     .               ksta,kend1,2.*dthtx,2.*dthty,2.*dthtz)
      else
        call rotateq(jdimp,kdimp,2,qroti,qroti,1,2,jsta,jend1,
     .               ksta,kend1,dthtx,dthty,dthtz)
      end if
c
c   Apply periodic BCs
      do 600 l=1,5
      do 600 k=ksta,kend1
      do 600 j=jsta,jend1
        qi0(j,k,l,3) = qroti(j,k,1,l)
        qi0(j,k,l,4) = qroti(j,k,2,l)
 600  continue
c
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
          vi0(j,k,1,3) = vp(j,k,1)
          vi0(j,k,1,4) = 0.0
  691   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.70.or. ivisc(2).ge.70.or. ivisc(1).ge.70) then

#   ifdef CMPLX
#   else
         ! obtain the rotation matrix and its transpose
         call stm2k5_get_rotmat(dthtx, dthty, dthtz, rn,rnt)

         ! rotate the stress tensors
         do k=ksta,kend1
            do j=jsta,jend1

               call stm2k5_bc(tp(j,k,1,1:7), rn, rnt,
     $              ti0(j,k,1:7,3))
               if(idimp==2) then
                  call stm2k5_bc(ti0(j,k,1:7,3), rn, rnt,
     $                 ti0(j,k,1:7,4))
               else 
                  call stm2k5_bc(tp(j,k,2,1:7), rn, rnt,
     $                 ti0(j,k,1:7,4))
               endif

            enddo
         enddo
#   endif
      elseif (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        idp2=min(2,idimp-1)
        do l=1,nummem
        do 601 k=ksta,kend1
        do 601 j=jsta,jend1
             ti0(j,k,l,3) = tp(j,k,1,l)
             ti0(j,k,l,4) = tp(j,k,idp2,l)
  601   continue
        enddo
      end if
      end if
c
      end if
c
      return
      end
