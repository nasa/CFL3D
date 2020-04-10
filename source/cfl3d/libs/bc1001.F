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
      subroutine bc1001(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,iuns,nou,bou,nbuf,
     .                  ibufdim,nummem,x,z)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set symmetry plane boundary conditions
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
      dimension x(jdim,kdim,idim),z(jdim,kdim,idim)
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /specialtop_kmax1001/ i_specialtop_kmax1001,
     .        a_specialtop_kmax1001,xc_specialtop_kmax1001,
     .        sig_specialtop_kmax1001,vtp_specialtop_kmax1001,
     .        wc_specialtop_kmax1001,fac_specialtop_kmax1001,
     .        cc_specialtop_kmax1001,xerf_specialtop_kmax1001,
     .        sigerf_specialtop_kmax1001

      real :: rn(3,3),rnt(3,3) ! rotational matrices for ivisc>=70 (stress models)
      real :: tloc(6), tghost(6)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=1001 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary            symmetry plane                       bctype 1001
c******************************************************************************
c
      if (nface.eq.3) then
c
      do 38 i=ista,iend1
      do 38 k=ksta,kend1
c
      vcont1 =  q(1,k,i,2)*sj(1,k,i,1) +
     .          q(1,k,i,3)*sj(1,k,i,2) +
     .          q(1,k,i,4)*sj(1,k,i,3) + sj(1,k,i,5)
      vcont2 =  q(2,k,i,2)*sj(1,k,i,1) +
     .          q(2,k,i,3)*sj(1,k,i,2) +
     .          q(2,k,i,4)*sj(1,k,i,3) + sj(1,k,i,5)
c
      qj0(k,i,1,1) = q(1,k,i,1)
      qj0(k,i,2,1) = q(1,k,i,2) - 2.*vcont1*sj(1,k,i,1)
      qj0(k,i,3,1) = q(1,k,i,3) - 2.*vcont1*sj(1,k,i,2)
      qj0(k,i,4,1) = q(1,k,i,4) - 2.*vcont1*sj(1,k,i,3)
      qj0(k,i,5,1) = q(1,k,i,5)
c
      qj0(k,i,1,2) = q(2,k,i,1)
      qj0(k,i,2,2) = q(2,k,i,2) - 2.*vcont2*sj(1,k,i,1)
      qj0(k,i,3,2) = q(2,k,i,3) - 2.*vcont2*sj(1,k,i,2)
      qj0(k,i,4,2) = q(2,k,i,4) - 2.*vcont2*sj(1,k,i,3)
      qj0(k,i,5,2) = q(2,k,i,5)
c
      bcj(k,i,1) = 0.0
c
   38 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,2) = vist3d(2,k,i)
  191   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if(ivisc(3).ge.70.or.ivisc(2).ge.70.or.ivisc(1).ge.70) then

#   ifdef CMPLX
#   else
        do i=ista,iend1
           do k=ksta,kend1
! rn is the rotational matrix
              call get_n1n2n3(sj(1,k,i,1), sj(1,k,i,2),sj(1,k,i,3),rn)
              rnt(:,1) = rn(1,:)
              rnt(:,2) = rn(2,:)
              rnt(:,3) = rn(3,:)
              call tensor_rotate(tursav(1,k,i,1:6),rn, rnt,tloc)
              call tensor_sym(tloc, tghost)
              call tensor_rotate(tghost,rnt,rn,tloc)
              
              tj0(k,i,1:6,1) = tloc(1:6)
              tj0(k,i,7,1)   = tursav(1,k,i,7) 
              
              call tensor_rotate(tursav(2,k,i,1:6),rn, rnt,tloc)
              call tensor_sym(tloc, tghost)
              call tensor_rotate(tghost,rnt,rn,tloc)
              
              tj0(k,i,1:6,2) = tloc(1:6)
              tj0(k,i,7,2)   = tursav(2,k,i,7) 
           enddo
        enddo   
#   endif
      elseif(ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 101 i=ista,iend1
        do 101 k=ksta,kend1
          tj0(k,i,l,1) = tursav(1,k,i,l)
          tj0(k,i,l,2) = tursav(2,k,i,l)
  101   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      j=jdim boundary         symmetry plane                       bctype 1001
c******************************************************************************
c
      if (nface.eq.4) then
c
      do 39 i=ista,iend1
      do 39 k=ksta,kend1
c
      vcont1 =  q(jdim-1,k,i,2)*sj(jdim,k,i,1) +
     .          q(jdim-1,k,i,3)*sj(jdim,k,i,2) +
     .          q(jdim-1,k,i,4)*sj(jdim,k,i,3) + sj(jdim,k,i,5)
      vcont2 =  q(jdim-2,k,i,2)*sj(jdim,k,i,1) +
     .          q(jdim-2,k,i,3)*sj(jdim,k,i,2) +
     .          q(jdim-2,k,i,4)*sj(jdim,k,i,3) + sj(jdim,k,i,5)
c
      qj0(k,i,1,3) = q(jdim-1,k,i,1)
      qj0(k,i,2,3) = q(jdim-1,k,i,2) - 2.*vcont1*sj(jdim,k,i,1)
      qj0(k,i,3,3) = q(jdim-1,k,i,3) - 2.*vcont1*sj(jdim,k,i,2)
      qj0(k,i,4,3) = q(jdim-1,k,i,4) - 2.*vcont1*sj(jdim,k,i,3)
      qj0(k,i,5,3) = q(jdim-1,k,i,5)
c
      qj0(k,i,1,4) = q(jdim-2,k,i,1)
      qj0(k,i,2,4) = q(jdim-2,k,i,2) - 2.*vcont2*sj(jdim,k,i,1)
      qj0(k,i,3,4) = q(jdim-2,k,i,3) - 2.*vcont2*sj(jdim,k,i,2)
      qj0(k,i,4,4) = q(jdim-2,k,i,4) - 2.*vcont2*sj(jdim,k,i,3)
      qj0(k,i,5,4) = q(jdim-2,k,i,5)
c
      bcj(k,i,2) = 0.0
c
   39 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
          vj0(k,i,1,3) = vist3d(jdim-1,k,i)
          vj0(k,i,1,4) = vist3d(jdim-2,k,i)
  291   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if(ivisc(3)>=70.or.ivisc(2)>=70.or.ivisc(1)>=70) then

#   ifdef CMPLX
#   else
         DO i=ista,iend1
            DO k=ksta,kend1
               ! get local orthogonal coordinate system of sym-plane
               ! i.e., the rotational matrix, rn
               CALL get_n1n2n3(sj(jdim,k,i,1), 
     $                         sj(jdim,k,i,2),
     $                         sj(jdim,k,i,3),rn)
               ! rnt: the transpose of rn
               rnt(:,1) = rn(1,:)
               rnt(:,2) = rn(2,:)
               rnt(:,3) = rn(3,:)
               
               ! rotate the stress tensor to the local coordinate system
               ! tloc = RNT * T  * RN
               CALL tensor_rotate(tursav(JDIM-1,k,i,1:6),rn, rnt,tloc)
               ! apply symmetry condition within the local coordinate system
               CALL tensor_sym(tloc, tghost)
               ! rotate the tensor back the current coordinate system
               CALL tensor_rotate(tghost,rnt,rn,tloc)
               
               tj0(k,i,1:6,3) = tloc(1:6)
               tj0(k,i,7,3)   = tursav(JDIM-1,k,i,7) 
               
               CALL tensor_rotate(tursav(JDIM-2,k,i,1:6),rn, rnt,tloc)
               CALL tensor_sym(tloc, tghost)
               CALL tensor_rotate(tghost,rnt,rn,tloc)
               
               tj0(k,i,1:6,4) = tloc(1:6)
               tj0(k,i,7,4)   = tursav(JDIM-2,k,i,7) 
            ENDDO
         ENDDO
#   endif
      else if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4)then
        do l=1,nummem
        do 201 i=ista,iend1
        do 201 k=ksta,kend1
          tj0(k,i,l,3) = tursav(jdim-1,k,i,l)
          tj0(k,i,l,4) = tursav(jdim-2,k,i,l)
  201   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=1 boundary            symmetry plane                       bctype 1001
c******************************************************************************
c
      if (nface.eq.5) then
c
      do 48 i=ista,iend1
      do 48 j=jsta,jend1
c
      wcont1 =  q(j,1,i,2)*sk(j,1,i,1) +
     .          q(j,1,i,3)*sk(j,1,i,2) +
     .          q(j,1,i,4)*sk(j,1,i,3) + sk(j,1,i,5)
      wcont2 =  q(j,2,i,2)*sk(j,1,i,1) +
     .          q(j,2,i,3)*sk(j,1,i,2) +
     .          q(j,2,i,4)*sk(j,1,i,3) + sk(j,1,i,5)
c
      qk0(j,i,1,1) = q(j,1,i,1)
      qk0(j,i,2,1) = q(j,1,i,2) - 2.*wcont1*sk(j,1,i,1)
      qk0(j,i,3,1) = q(j,1,i,3) - 2.*wcont1*sk(j,1,i,2)
      qk0(j,i,4,1) = q(j,1,i,4) - 2.*wcont1*sk(j,1,i,3)
      qk0(j,i,5,1) = q(j,1,i,5)
c
      qk0(j,i,1,2) = q(j,2,i,1)
      qk0(j,i,2,2) = q(j,2,i,2) - 2.*wcont2*sk(j,1,i,1)
      qk0(j,i,3,2) = q(j,2,i,3) - 2.*wcont2*sk(j,1,i,2)
      qk0(j,i,4,2) = q(j,2,i,4) - 2.*wcont2*sk(j,1,i,3)
      qk0(j,i,5,2) = q(j,2,i,5)
c
      bck(j,i,1) = 0.0
c
   48 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = vist3d(j,2,i)
  391   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if(ivisc(3)>=70.or.ivisc(2)>=70.or.ivisc(1)>=70) then

#   ifdef CMPLX
#   else
         DO i=ista,iend1
            DO j= jsta, jend1
               ! get local orthogonal coordinate system of sym-plane
               ! i.e., the rotational matrix, rn
               CALL get_n1n2n3(sk(j,1,i,1), sk(j,1,i,2),sk(j,1,i,3),rn)
               ! rnt: the transpose of rn
               rnt(:,1) = rn(1,:)
               rnt(:,2) = rn(2,:)
               rnt(:,3) = rn(3,:)
               
               ! rotate the stress tensor to the local coordinate system
               ! tloc = RNT * T  * RN
               CALL tensor_rotate(tursav(j,1,i,1:6),rn, rnt,tloc)
               ! apply symmetry condition within the local coordinate system
               CALL tensor_sym(tloc, tghost)
               ! rotate the tensor back the current coordinate system
               CALL tensor_rotate(tghost,rnt,rn,tloc)
               
               tk0(j,i,1:6,1) = tloc(1:6)
               tk0(j,i,7,1)   = tursav(j,1,i,7) 
               
               CALL tensor_rotate(tursav(j,2,i,1:6),rn, rnt,tloc)
               CALL tensor_sym(tloc, tghost)
               CALL tensor_rotate(tghost,rnt,rn,tloc)
               
               tk0(j,i,1:6,2) = tloc(1:6)
               tk0(j,i,7,2)   = tursav(j,2,i,7) 
            ENDDO
         ENDDO
#   endif
      ELSE IF (ivisc(3).GE.4 .OR. ivisc(2).GE.4 .OR. ivisc(1).GE.4) THEN
        do l=1,nummem
        do 301 i=ista,iend1
        do 301 j=jsta,jend1
          tk0(j,i,l,1) = tursav(j,1,i,l)
          tk0(j,i,l,2) = tursav(j,2,i,l)
  301   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=kdim boundary         symmetry plane                       bctype 1001
c******************************************************************************
c
      if (nface.eq.6) then
c
c   Special for when i_specialtop_kmax1001=1
c   All sections with i_specialtop_kmax1001=1 are undocumented,
c   because they were inserted only for a particular application;
c   they could be commented out or deleted in general, if desired.
c   Adds the following Vtop (in terms of Uref):
c     = a*{fac*(xc-x)+(1-fac)}*exp[cc-({(xc-x)/sig}**2)] + 
c       vtp*{0.5*(1-erf[(x-xerf)/sigerf])}
c   with constant crossflow: wc
c   (note wc should be negative if input beta is positive)
c
      a=a_specialtop_kmax1001
      xc=xc_specialtop_kmax1001
      sig=sig_specialtop_kmax1001
      vtp=vtp_specialtop_kmax1001
      wc=wc_specialtop_kmax1001
      fac=fac_specialtop_kmax1001
      cc=cc_specialtop_kmax1001
      xerf=xerf_specialtop_kmax1001
      sigerf=sigerf_specialtop_kmax1001
      do 49 i=ista,iend1
      do 49 j=jsta,jend1
        if (i_specialtop_kmax1001 .eq. 1) then
          xxx=0.5*(x(j,kdim,i)+x(j+1,kdim,i))
          b=(xc-xxx)/sig
          dz=z(j,kdim,i)-z(j,kdim-1,i)
        end if
c
      wcont1 =  q(j,kdim-1,i,2)*sk(j,kdim,i,1) +
     .          q(j,kdim-1,i,3)*sk(j,kdim,i,2) +
     .          q(j,kdim-1,i,4)*sk(j,kdim,i,3) + sk(j,kdim,i,5) 
      wcont2 =  q(j,kdim-2,i,2)*sk(j,kdim,i,1) +
     .          q(j,kdim-2,i,3)*sk(j,kdim,i,2) +
     .          q(j,kdim-2,i,4)*sk(j,kdim,i,3) + sk(j,kdim,i,5)
c
      qk0(j,i,1,3) = q(j,kdim-1,i,1)
      if (i_specialtop_kmax1001 .eq. 1) then
        vtop=(fac*sig*a*b+(1.-fac)*a)*exp(cc-(b**2))*xmach+
     +    vtp*(0.5*(1.-erf((xxx-xerf)/sigerf)))*xmach
        dudy=fac*a*(2.0*(b**2)-1.)*exp(cc-(b**2))*xmach+
     +    (1.-fac)*a*2.0/sig*b*exp(cc-(b**2))*xmach
        qk0(j,i,2,3) = q(j,kdim-1,i,2) + dz*dudy
        if (wc .eq. 0.) then
          qk0(j,i,3,3) = q(j,kdim-1,i,3) - 2.*wcont1*sk(j,kdim,i,2)
        else
          qk0(j,i,3,3) = 2.*wc*xmach - q(j,kdim-1,i,3)
        end if
        qk0(j,i,4,3) = 2.*vtop - q(j,kdim-1,i,4)
        qk0(j,i,5,3) = q(j,kdim-1,i,5)
      else
        qk0(j,i,2,3) = q(j,kdim-1,i,2) - 2.*wcont1*sk(j,kdim,i,1)
        qk0(j,i,3,3) = q(j,kdim-1,i,3) - 2.*wcont1*sk(j,kdim,i,2)
        qk0(j,i,4,3) = q(j,kdim-1,i,4) - 2.*wcont1*sk(j,kdim,i,3)
        qk0(j,i,5,3) = q(j,kdim-1,i,5)
      end if
c
      qk0(j,i,1,4) = q(j,kdim-2,i,1)
      if (i_specialtop_kmax1001 .eq. 1) then
        qk0(j,i,2,4) = 2.*qk0(j,i,2,3) - q(j,kdim-1,i,2)
        if (wc .eq. 0.) then
          qk0(j,i,3,4) = q(j,kdim-2,i,3) - 2.*wcont2*sk(j,kdim,i,2)
        else
          qk0(j,i,3,4) = 2.*qk0(j,i,3,3) - q(j,kdim-1,i,3)
        end if
        qk0(j,i,4,4) = 2.*qk0(j,i,4,3) - q(j,kdim-1,i,4)
      else
        qk0(j,i,2,4) = q(j,kdim-2,i,2) - 2.*wcont2*sk(j,kdim,i,1)
        qk0(j,i,3,4) = q(j,kdim-2,i,3) - 2.*wcont2*sk(j,kdim,i,2)
        qk0(j,i,4,4) = q(j,kdim-2,i,4) - 2.*wcont2*sk(j,kdim,i,3)
      end if
      qk0(j,i,5,4) = q(j,kdim-2,i,5)
c
      bck(j,i,2) = 0.0
c
   49 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
          vk0(j,i,1,3) = vist3d(j,kdim-1,i)
          vk0(j,i,1,4) = vist3d(j,kdim-2,i)
  491   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if(ivisc(3)>=70.or.ivisc(2)>=70.or.ivisc(1)>=70) then

#   ifdef CMPLX
#   else
         DO i=ista,iend1
            DO j= jsta, jend1
               ! get local orthogonal coordinate system of sym-plane
               ! i.e., the rotational matrix, rn
               CALL get_n1n2n3(sk(j,kdim,i,1),
     $                         sk(j,kdim,i,2),
     $                         sk(j,kdim,i,3),rn)
               ! rnt: the transpose of rn
               rnt(:,1) = rn(1,:)
               rnt(:,2) = rn(2,:)
               rnt(:,3) = rn(3,:)
               
               ! rotate the stress tensor to the local coordinate system
               ! tloc = RNT * T  * RN
               CALL tensor_rotate(tursav(j,kdim-1,i,1:6),rn, rnt,tloc)
               ! apply symmetry condition within the local coordinate system
               CALL tensor_sym(tloc, tghost)
               ! rotate the tensor back the current coordinate system
               CALL tensor_rotate(tghost,rnt,rn,tloc)
               
               tk0(j,i,1:6,3) = tloc(1:6)
               tk0(j,i,7,3)   = tursav(j,kdim-1,i,7) 
               
               CALL tensor_rotate(tursav(j,kdim-2,i,1:6),rn, rnt,tloc)
               CALL tensor_sym(tloc, tghost)
               CALL tensor_rotate(tghost,rnt,rn,tloc)
               
               tk0(j,i,1:6,4) = tloc(1:6)
               tk0(j,i,7,4)   = tursav(j,kdim-2,i,7) 
            ENDDO
         ENDDO
#   endif
      ELSE if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 401 i=ista,iend1
        do 401 j=jsta,jend1
          tk0(j,i,l,3) = tursav(j,kdim-1,i,l)
          tk0(j,i,l,4) = tursav(j,kdim-2,i,l)
  401   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=1 boundary            symmetry plane                       bctype 1001
c******************************************************************************
c
      if (nface.eq.1) then
c
      i2 = min(2,idim1)

      do 58 k=ksta,kend1
      do 58 j=jsta,jend1
c
      ucont1 =  q(j,k,1,2)*si(j,k,1,1) +
     .          q(j,k,1,3)*si(j,k,1,2) +
     .          q(j,k,1,4)*si(j,k,1,3) + si(j,k,1,5)
      ucont2 =  q(j,k,i2,2)*si(j,k,1,1) +
     .          q(j,k,i2,3)*si(j,k,1,2) +
     .          q(j,k,i2,4)*si(j,k,1,3) + si(j,k,1,5)
c
      qi0(j,k,1,1) = q(j,k,1,1)
      qi0(j,k,2,1) = q(j,k,1,2) - 2.*ucont1*si(j,k,1,1)
      qi0(j,k,3,1) = q(j,k,1,3) - 2.*ucont1*si(j,k,1,2)
      qi0(j,k,4,1) = q(j,k,1,4) - 2.*ucont1*si(j,k,1,3)
      qi0(j,k,5,1) = q(j,k,1,5)
c
      qi0(j,k,1,2) = q(j,k,i2,1)
      qi0(j,k,2,2) = q(j,k,i2,2) - 2.*ucont2*si(j,k,1,1)
      qi0(j,k,3,2) = q(j,k,i2,3) - 2.*ucont2*si(j,k,1,2)
      qi0(j,k,4,2) = q(j,k,i2,4) - 2.*ucont2*si(j,k,1,3)
      qi0(j,k,5,2) = q(j,k,i2,5)
c
      bci(j,k,1) = 0.0
c
   58 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = vist3d(j,k,i2)
  591   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.70 .or. ivisc(2).ge.70 .or. ivisc(1).ge.70) then
#   ifdef CMPLX
#   else
         DO k=ksta,kend1
            DO j=jsta,jend1
               ! get local orthogonal coordinate system of sym-plane
               ! i.e., the rotational matrix, rn
               CALL get_n1n2n3(si(j,k,1,1), si(j,k,1,2),si(j,k,1,3),rn)
               ! rnt: the transpose of rn
               rnt(:,1) = rn(1,:)
               rnt(:,2) = rn(2,:)
               rnt(:,3) = rn(3,:)
               
               ! rotate the stress tensor to the local coordinate system
               ! tloc = RNT * T  * RN
               CALL tensor_rotate(tursav(j,k,1,1:6),rn, rnt,tloc)
               ! apply symmetry condition within the local coordinate system
               CALL tensor_sym(tloc, tghost)
               ! rotate the tensor back the current coordinate system
               CALL tensor_rotate(tghost,rnt,rn,tloc)
               
               ti0(j,k,1:6,1) = tloc(1:6)
               ti0(j,k,7,1)   = tursav(j,k,1,7) 
               
               CALL tensor_rotate(tursav(j,k,i2,1:6),rn, rnt,tloc)
               CALL tensor_sym(tloc, tghost)
               CALL tensor_rotate(tghost,rnt,rn,tloc)
               
               ti0(j,k,1:6,2) = tloc(1:6)
               ti0(j,k,7,2)   = tursav(j,k,i2,7) 
            ENDDO
         ENDDO
#   endif
      else if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 501 k=ksta,kend1
        do 501 j=jsta,jend1
          ti0(j,k,l,1) = tursav(j,k,1,l)
          ti0(j,k,l,2) = tursav(j,k,i2,l)
  501   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=idim boundary         symmetry plane                       bctype 1001
c******************************************************************************
c
      if (nface.eq.2) then
c
      i2 = max(1,idim-2)
c
      do 59 k=ksta,kend1
      do 59 j=jsta,jend1
c
      ucont1 =  q(j,k,idim-1,2)*si(j,k,idim,1) +
     .          q(j,k,idim-1,3)*si(j,k,idim,2) +
     .          q(j,k,idim-1,4)*si(j,k,idim,3) + si(j,k,idim,5)
      ucont2 =  q(j,k,i2,2)*si(j,k,idim,1) +
     .          q(j,k,i2,3)*si(j,k,idim,2) +
     .          q(j,k,i2,4)*si(j,k,idim,3) + si(j,k,idim,5)

c
      qi0(j,k,1,3) = q(j,k,idim-1,1)
      qi0(j,k,2,3) = q(j,k,idim-1,2) - 2.*ucont1*si(j,k,idim,1)
      qi0(j,k,3,3) = q(j,k,idim-1,3) - 2.*ucont1*si(j,k,idim,2)
      qi0(j,k,4,3) = q(j,k,idim-1,4) - 2.*ucont1*si(j,k,idim,3)
      qi0(j,k,5,3) = q(j,k,idim-1,5)
c
      qi0(j,k,1,4) = q(j,k,i2,1)
      qi0(j,k,2,4) = q(j,k,i2,2) - 2.*ucont2*si(j,k,idim,1)
      qi0(j,k,3,4) = q(j,k,i2,3) - 2.*ucont2*si(j,k,idim,2)
      qi0(j,k,4,4) = q(j,k,i2,4) - 2.*ucont2*si(j,k,idim,3)
      qi0(j,k,5,4) = q(j,k,i2,5)
c
      bci(j,k,2) = 0.0
c
   59 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
          vi0(j,k,1,3) = vist3d(j,k,idim-1)
          vi0(j,k,1,4) = vist3d(j,k,i2)
  691   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.70 .or. ivisc(2).ge.70 .or. ivisc(1).ge.70) then

#   ifdef CMPLX
#   else
         DO k=ksta,kend1
            DO j=jsta,jend1
               ! get local orthogonal coordinate system of sym-plane
               ! i.e., the rotational matrix, rn
               CALL get_n1n2n3(si(j,k,idim-1,1), 
     $                         si(j,k,idim-1,2),
     $                         si(j,k,idim-1,3),rn)
               ! rnt: the transpose of rn
               rnt(:,1) = rn(1,:)
               rnt(:,2) = rn(2,:)
               rnt(:,3) = rn(3,:)
               
               ! rotate the stress tensor to the local coordinate system
               ! tloc = RNT * T  * RN
               CALL tensor_rotate(tursav(j,k,idim-1,1:6),rn, rnt,tloc)
               ! apply symmetry condition within the local coordinate system
               CALL tensor_sym(tloc, tghost)
               ! rotate the tensor back the current coordinate system
               CALL tensor_rotate(tghost,rnt,rn,tloc)
               
               ti0(j,k,1:6,3) = tloc(1:6)
               ti0(j,k,7,3)   = tursav(j,k,idim-1,7) 
               
               CALL tensor_rotate(tursav(j,k,i2,1:6),rn, rnt,tloc)
               CALL tensor_sym(tloc, tghost)
               CALL tensor_rotate(tghost,rnt,rn,tloc)
               
               ti0(j,k,1:6,4) = tloc(1:6)
               ti0(j,k,7,4)   = tursav(j,k,i2,7)
            ENDDO
         ENDDO
#   endif
      else if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 601 k=ksta,kend1
        do 601 j=jsta,jend1
          ti0(j,k,l,3) = tursav(j,k,idim-1,l)
          ti0(j,k,l,4) = tursav(j,k,i2,l)
  601   continue
        enddo
      end if
      end if
      end if
c
      return
#   ifdef CMPLX
#   else
      contains

      SUBROUTINE get_n1n2n3(x1,x2,x3, rn)
      implicit none
      real, intent(in) :: x1,x2,x3
      real, intent(out) :: rn(3,3)
      
      real :: dv1(3), dv2(3), dv3(3)
      real :: rm,rmin
      integer :: i_min,i
      
      rn(1:3,1) = (/x1,x2,x3/)
      rm = sqrt(rn(1,1)**2+rn(2,1)**2+rn(3,1)**2)
      rn(1:3,1) = rn(1:3,1)/rm
      rmin = huge(rmin)
      do i=1,3
         if(rmin>abs(rn(i,1))) then
            rmin = abs(rn(i,1))
            i_min = i
         endif
      enddo
      
      rn(1:3,3) = 0
      rn(i_min,3) = 1.
      call cross(rn(1,1),rn(1,3),rn(1,2))
      rm = sqrt(rn(1,2)**2+rn(2,2)**2+rn(3,2)**2)
      rn(:,2) = rn(:,2)/rm

      call cross(rn(1,1),rn(1,2), rn(1,3))

      return
      end  subroutine

      subroutine cross(v1,v2,v3)
      implicit none
      real, intent(in) ::v1(3), v2(3)
      real, intent(out) :: v3(3)
      
      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
      
      return
      end subroutine

      subroutine tensor_rotate(t, rn, rnt, tloc)
      implicit none
      real, intent(in) :: t(6),rn(3,3),rnt(3,3)
      real, intent(out) :: tloc(6)

      real :: t33(3,3),s(3,3),v(3,3)
      
      t33(1,1)= t(1)
      t33(2,2)= t(2)
      t33(3,3)= t(3)

      t33(1,2)= t(4)
      t33(2,3)= t(5)
      t33(1,3)= t(6)

      t33(2,1)= t33(1,2)
      t33(3,2)= t33(2,3)
      t33(3,1)= t33(1,3)

      s(:,1) = t33(:,1)*rn(1,1) + t33(:,2)*rn(2,1) + t33(:,3)*rn(3,1)
      s(:,2) = t33(:,1)*rn(1,2) + t33(:,2)*rn(2,2) + t33(:,3)*rn(3,2)
      s(:,3) = t33(:,1)*rn(1,3) + t33(:,2)*rn(2,3) + t33(:,3)*rn(3,3)

      v(:,1) = rnt(:,1)*s(1,1) + rnt(:,2)*s(2,1) + rnt(:,3)*s(3,1)
      v(:,2) = rnt(:,1)*s(1,2) + rnt(:,2)*s(2,2) + rnt(:,3)*s(3,2)
      v(:,3) = rnt(:,1)*s(1,3) + rnt(:,2)*s(2,3) + rnt(:,3)*s(3,3)
      
      tloc(1) = v(1,1)
      tloc(2) = v(2,2)
      tloc(3) = v(3,3)
      tloc(4) = v(1,2)
      tloc(5) = v(2,3)
      tloc(6) = v(1,3)
      return
      end subroutine
      
      subroutine tensor_sym(tloc, tghost)
      implicit none
      real, intent(in) :: tloc(6)
      real, intent(out):: tghost(6)
      
      tghost(1) = tloc(1)
      tghost(2) = tloc(2)
      tghost(3) = tloc(3)
      tghost(4) = -tloc(4)
      tghost(5) = tloc(5)
      tghost(6) = -tloc(6)
      return

      end subroutine
#   endif
      end
            
