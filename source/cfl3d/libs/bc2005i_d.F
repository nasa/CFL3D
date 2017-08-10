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
      subroutine bc2005i_d(jdim,kdim,idim,datai0,ista,iend,jsta,jend,
     .                  ksta,kend,nface,mdim,ndim,bcdata,filname,
     .                  jdimp,kdimp,idimp,datapi,nbl,nblp,ldim,
     .                  nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl,iflag)
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
c     iflag  = 0 for Qs, 1 for vist3d, 2 for turbulence quantities
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
c       data              = data on current block; corresponds to q
c                           if ldim=5, vist3d if ldim=1 and tursav
c                           if ldim=nummem
c       datai0            = BC values of data assigned for this block
c       ista,iend,etc.    = indices over which BC is applied
c       nface             = face number (1 = i=1  2 = i=idim)
c       mdim,ndim         = dimensions of bcdata
c       bcdata            = auxiliary data that goes with this BC
c       filename          = filename to read bcdata, if array values
c       jdimp,kdimp,idimp = dimensions of periodic block
c       datapi            = work array of data values on periodic block
c       nbl               = block number of current block
c       nblp              = block number of periodic block
c       ldim              = last dimension of data array, determines
c                           the interpretation of the info stored in
c                           the data array (see above)
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
      dimension datai0(jdim,kdim,ldim,4)
      dimension datapi(jdimp,kdimp,2,ldim)
      dimension bcdata(mdim,ndim,2,12),mblk2nd(maxbl)
c
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /conversion/ radtodeg
      common /igrdtyp/ ip3dgrd,ialph
      real :: rn(3,3), rnt(3,3) ! rotation matrix and its transpose
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
c      i=1 boundary             periodic boundary                   bctype 2005
c******************************************************************************
c
      if (nface.eq.1) then
c
      if (iflag.eq.0) then
c
c        rotate datapk values before setting datak0 values if data
c        corresponds to q
c
         if (idimp .eq. 2) then
            call rotateq(jdimp,kdimp,2,datapi,datapi,1,1,jsta,jend1,
     .                   ksta,kend1,dthtx,dthty,dthtz)
            call rotateq(jdimp,kdimp,2,datapi,datapi,2,2,jsta,jend1,
     .                   ksta,kend1,2.*dthtx,2.*dthty,2.*dthtz)
         else
            call rotateq(jdimp,kdimp,2,datapi,datapi,1,2,jsta,jend1,
     .                   ksta,kend1,dthtx,dthty,dthtz)
         end if
c
         do 500 k=ksta,kend1
         do 500 j=jsta,jend1
         do 500 l=1,5
           datai0(j,k,l,1) = datapi(j,k,1,l)
           datai0(j,k,l,2) = datapi(j,k,2,l)
 500     continue
c
      end if
c
      if (iflag.eq.1) then
c
         if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
           do 591 k=ksta,kend1
           do 591 j=jsta,jend1
             datai0(j,k,1,1) = datapi(j,k,1,1)
             datai0(j,k,1,2) = 0.0
  591      continue
         end if
c
      end if
c
      if (iflag.eq.2) then
c
c        only need to do advanced model turbulence B.C.s on finest grid
c
         if (level .ge. lglobal) then
            if(ivisc(3).ge.70 .or. ivisc(2).ge.70.or.ivisc(1).ge.70)then

#   ifdef CMPLX
#   else
               call stm2k5_get_rotmat(dthtx, dthty, dthtz,rn,rnt)
               do k=ksta,kend1
                  do  j=jsta,jend1
                 
                     call stm2k5_bc(datapi(j,k,1,1:7),rn,rnt, 
     $                              datai0(j,k,1:7,1))
                     if(idimp==2) then
                     call stm2k5_bc(datai0(j,k,1:7,1),rn,rnt, 
     $                              datai0(j,k,1:7,2))
                     else
                     call stm2k5_bc(datapi(j,k,2,1:7),rn,rnt, 
     $                              datai0(j,k,1:7,2))
                     endif 
                    
                  enddo
               enddo
#   endif
            elseif(ivisc(3).ge.4.or.ivisc(2).ge.4.or.ivisc(1).ge.4)then

               do 501 k=ksta,kend1
               do 501 j=jsta,jend1
               do 501 l=1,ldim
                  datai0(j,k,l,1) = datapi(j,k,1,l)
                  datai0(j,k,l,2) = datapi(j,k,2,l)
  501          continue
            end if
         end if
c
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
      if (iflag.eq.0) then
c
c        rotate datapk values before setting datak0 values if data
c        corresponds to q
c
         if (idimp .eq. 2) then
            call rotateq(jdimp,kdimp,2,datapi,datapi,1,1,jsta,jend1,
     .                   ksta,kend1,dthtx,dthty,dthtz)
            call rotateq(jdimp,kdimp,2,datapi,datapi,2,2,jsta,jend1,
     .                   ksta,kend1,2.*dthtx,2.*dthty,2.*dthtz)
         else
            call rotateq(jdimp,kdimp,2,datapi,datapi,1,2,jsta,jend1,
     .                   ksta,kend1,dthtx,dthty,dthtz)
          end if
c
          do 600 l=1,5
          do 600 k=ksta,kend1
          do 600 j=jsta,jend1
            datai0(j,k,l,3) = datapi(j,k,1,l)
            datai0(j,k,l,4) = datapi(j,k,2,l)
 600      continue
c
      end if
c
      if (iflag.eq.1) then
c
         if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
            do 691 k=ksta,kend1
            do 691 j=jsta,jend1
              datai0(j,k,1,3) = datapi(j,k,1,1)
              datai0(j,k,1,4) = 0.0
  691       continue
         end if
c
      end if
c
      if (iflag.eq.2) then
c
c        only need to do advanced model turbulence B.C.s on finest grid
c
         if (level .ge. lglobal) then
            if(ivisc(3).ge.70 .or. ivisc(2).ge.70.or.ivisc(1).ge.70)then

#   ifdef CMPLX
#   else
               call stm2k5_get_rotmat(dthtx, dthty, dthtz,rn,rnt)
               do k=ksta,kend1
               do j=jsta,jend1
                  call stm2k5_bc(datapi(j,k,1,1:7),rn, rnt,
     $                           datai0(j,k,1:7,3))
                  if(idimp==2) then
                  call stm2k5_bc(datai0(j,k,1:7,3),rn, rnt,
     $                           datai0(j,k,1:7,4))
                  else
                  call stm2k5_bc(datapi(j,k,2,1:7),rn, rnt,
     $                           datai0(j,k,1:7,4))
                  endif
               enddo
               enddo
#   endif
            elseif(ivisc(3).ge.4.or.ivisc(2).ge.4.or.ivisc(1).ge.4)then
               do 601 k=ksta,kend1
               do 601 j=jsta,jend1
               do 601 l=1,ldim
                 datai0(j,k,l,3) = datapi(j,k,1,l)
                 datai0(j,k,l,4) = datapi(j,k,2,l)
  601          continue
            end if
         end if
c
      end if
c
      end if
c
      return
      end
