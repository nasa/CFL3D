      subroutine bc1011(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,isym,jsym,ksym,iuns,
     .                  nou,bou,nbuf,ibufdim)
c***********************************************************************
c     Purpose:  Set singular axis - half plane boundary conditions
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
      dimension tursav(jdim,kdim,idim,2),tj0(kdim,idim-1,2,4),
     .          tk0(jdim,idim-1,2,4),ti0(jdim,kdim,2,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
c
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c     works for half-plane symmetry only - assume checks for appropriateness
c     of this boundary condition have been made PRIOR to entering this routine
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=1011 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary            singular axis - half plane           bctype 1011
c******************************************************************************
c
      if (nface.eq.3) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1001)ista,iend,ksta,kend
      end if
c
c     symmetry in k
c
      if(ksym.gt.0) then
        do 35 i=ista,iend1
        do 35 k=ksta,kend1
c
        vcont1 =  q(1,kend-k,i,2)*sk(1,kend,i,1) +
     .            q(1,kend-k,i,3)*sk(1,kend,i,2) +
     .            q(1,kend-k,i,4)*sk(1,kend,i,3) + sk(1,kend,i,5)
        vcont2 =  q(2,kend-k,i,2)*sk(2,kend,i,1) +
     .            q(2,kend-k,i,3)*sk(2,kend,i,2) +
     .            q(2,kend-k,i,4)*sk(2,kend,i,3) + sk(2,kend,i,5)
c
        qj0(k,i,1,1) = q(1,kend-k,i,1)
        qj0(k,i,2,1) = q(1,kend-k,i,2) - 2.*vcont1*sk(1,kend,i,1)
        qj0(k,i,3,1) = q(1,kend-k,i,3) - 2.*vcont1*sk(1,kend,i,2)
        qj0(k,i,4,1) = q(1,kend-k,i,4) - 2.*vcont1*sk(1,kend,i,3)
        qj0(k,i,5,1) = q(1,kend-k,i,5)
c
        qj0(k,i,1,2) = q(2,kend-k,i,1)
        qj0(k,i,2,2) = q(2,kend-k,i,2) - 2.*vcont2*sk(2,kend,i,1)
        qj0(k,i,3,2) = q(2,kend-k,i,3) - 2.*vcont2*sk(2,kend,i,2)
        qj0(k,i,4,2) = q(2,kend-k,i,4) - 2.*vcont2*sk(2,kend,i,3)
        qj0(k,i,5,2) = q(2,kend-k,i,5)
c
        bcj(k,i,1) = 0.0
c
   35   continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 191 i=ista,iend1
          do 191 k=ksta,kend1
            vj0(k,i,1,1) = vist3d(1,kend-k,i)
            vj0(k,i,1,2) = vist3d(2,kend-k,i)
  191     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 101 i=ista,iend1
          do 101 k=ksta,kend1
            tj0(k,i,1,1) = tursav(1,kend-k,i,1)
            tj0(k,i,2,1) = tursav(1,kend-k,i,2)
            tj0(k,i,1,2) = tursav(2,kend-k,i,1)
            tj0(k,i,2,2) = tursav(2,kend-k,i,2)
  101     continue
        end if
        end if
c
      end if
c
c     symmetry in i
c
      if(isym.gt.0) then
        do 36 i=ista,iend1
        do 36 k=ksta,kend1
c
        vcont1 =  q(1,k,iend-i,2)*si(1,k,iend,1) +
     .            q(1,k,iend-i,3)*si(1,k,iend,2) +
     .            q(1,k,iend-i,4)*si(1,k,iend,3) + si(1,k,iend,5)
        vcont2 =  q(2,k,iend-i,2)*si(2,k,iend,1) +
     .            q(2,k,iend-i,3)*si(2,k,iend,2) +
     .            q(2,k,iend-i,4)*si(2,k,iend,3) + si(2,k,iend,5)
c
        qj0(k,i,1,1) = q(1,k,iend-i,1)
        qj0(k,i,2,1) = q(1,k,iend-i,2) - 2.*vcont1*si(1,k,iend,1)
        qj0(k,i,3,1) = q(1,k,iend-i,3) - 2.*vcont1*si(1,k,iend,2)
        qj0(k,i,4,1) = q(1,k,iend-i,4) - 2.*vcont1*si(1,k,iend,3)
        qj0(k,i,5,1) = q(1,k,iend-i,5)
c
        qj0(k,i,1,2) = q(2,k,iend-i,1)
        qj0(k,i,2,2) = q(2,k,iend-i,2) - 2.*vcont2*si(2,k,iend,1)
        qj0(k,i,3,2) = q(2,k,iend-i,3) - 2.*vcont2*si(2,k,iend,2)
        qj0(k,i,4,2) = q(2,k,iend-i,4) - 2.*vcont2*si(2,k,iend,3)
        qj0(k,i,5,2) = q(2,k,iend-i,5)
c
        bcj(k,i,1) = 0.0
c
   36   continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 192 i=ista,iend1
          do 192 k=ksta,kend1
            vj0(k,i,1,1) = vist3d(1,k,iend-i)
            vj0(k,i,1,2) = vist3d(2,k,iend-i)
  192     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 102 i=ista,iend1
          do 102 k=ksta,kend1
            tj0(k,i,1,1) = tursav(1,k,iend-i,1)
            tj0(k,i,2,1) = tursav(1,k,iend-i,2)
            tj0(k,i,1,2) = tursav(2,k,iend-i,1)
            tj0(k,i,2,2) = tursav(2,k,iend-i,2)
  102     continue
        end if
        end if
c
      end if
c
      end if
c
c******************************************************************************
c      j=jdim boundary         singular axis - half plane           bctype 1011
c******************************************************************************
c
      if (nface.eq.4) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1002)ista,iend,ksta,kend
      end if
c
c     symmetry in k
c
      if(ksym.gt.0) then
        do 38 i=ista,iend1
        do 38 k=ksta,kend1
c
        vcont1 =  q(jdim-1,kend-k,i,2)*sk(jdim-1,kend,i,1) +
     .            q(jdim-1,kend-k,i,3)*sk(jdim-1,kend,i,2) +
     .            q(jdim-1,kend-k,i,4)*sk(jdim-1,kend,i,3) + 
     .                                 sk(jdim-1,kend,i,5)
        vcont2 =  q(jdim-2,kend-k,i,2)*sk(jdim-2,kend,i,1) +
     .            q(jdim-2,kend-k,i,3)*sk(jdim-2,kend,i,2) +
     .            q(jdim-2,kend-k,i,4)*sk(jdim-2,kend,i,3) + 
     .                                 sk(jdim-2,kend,i,5)
c
        qj0(k,i,1,3) = q(jdim-1,kend-k,i,1)
        qj0(k,i,2,3) = q(jdim-1,kend-k,i,2) 
     .               - 2.*vcont1*sk(jdim-1,kend,i,1)
        qj0(k,i,3,3) = q(jdim-1,kend-k,i,3) 
     .               - 2.*vcont1*sk(jdim-1,kend,i,2)
        qj0(k,i,4,3) = q(jdim-1,kend-k,i,4) 
     .               - 2.*vcont1*sk(jdim-1,kend,i,3)
        qj0(k,i,5,3) = q(jdim-1,kend-k,i,5)
c
        qj0(k,i,1,4) = q(jdim-2,kend-k,i,1)
        qj0(k,i,2,4) = q(jdim-2,kend-k,i,2) 
     .               - 2.*vcont2*sk(jdim-2,kend,i,1)
        qj0(k,i,3,4) = q(jdim-2,kend-k,i,3) 
     .               - 2.*vcont2*sk(jdim-2,kend,i,2)
        qj0(k,i,4,4) = q(jdim-2,kend-k,i,4) 
     .               - 2.*vcont2*sk(jdim-2,kend,i,3)
        qj0(k,i,5,4) = q(jdim-2,kend-k,i,5)
c
        bcj(k,i,2) = 0.0
c
   38   continue
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 291 i=ista,iend1
          do 291 k=ksta,kend1
            vj0(k,i,1,3) = vist3d(jdim-1,kend-k,i)
            vj0(k,i,1,4) = vist3d(jdim-2,kend-k,i)
  291     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 201 i=ista,iend1
          do 201 k=ksta,kend1
            tj0(k,i,1,3) = tursav(jdim-1,kend-k,i,1)
            tj0(k,i,2,3) = tursav(jdim-1,kend-k,i,2)
            tj0(k,i,1,4) = tursav(jdim-2,kend-k,i,1)
            tj0(k,i,2,4) = tursav(jdim-2,kend-k,i,2)
  201     continue
        end if
        end if
c
      end if
c
c     symmetry in i
c
      if(isym.gt.0) then
        do 39 i=ista,iend1
        do 39 k=ksta,kend1
c
        vcont1 =  q(jdim-1,k,iend-i,2)*si(jdim-1,k,iend,1) +
     .            q(jdim-1,k,iend-i,3)*si(jdim-1,k,iend,2) +
     .            q(jdim-1,k,iend-i,4)*si(jdim-1,k,iend,3) + 
     .                                 si(jdim-1,k,iend,5)
        vcont2 =  q(jdim-2,k,iend-i,2)*si(jdim-2,k,iend,1) +
     .            q(jdim-2,k,iend-i,3)*si(jdim-2,k,iend,2) +
     .            q(jdim-2,k,iend-i,4)*si(jdim-2,k,iend,3) + 
     .                                 si(jdim-2,k,iend,5)
c
        qj0(k,i,1,3) = q(jdim-1,k,iend-i,1)
        qj0(k,i,2,3) = q(jdim-1,k,iend-i,2) 
     .               - 2.*vcont1*si(jdim-1,k,iend,1)
        qj0(k,i,3,3) = q(jdim-1,k,iend-i,3) 
     .               - 2.*vcont1*si(jdim-1,k,iend,2)
        qj0(k,i,4,3) = q(jdim-1,k,iend-i,4) 
     .               - 2.*vcont1*si(jdim-1,k,iend,3)
        qj0(k,i,5,3) = q(jdim-1,k,iend-i,5)
c
        qj0(k,i,1,4) = q(jdim-2,k,iend-i,1)
        qj0(k,i,2,4) = q(jdim-2,k,iend-i,2) 
     .               - 2.*vcont2*si(jdim-2,k,iend,1)
        qj0(k,i,3,4) = q(jdim-2,k,iend-i,3) 
     .               - 2.*vcont2*si(jdim-2,k,iend,2)
        qj0(k,i,4,4) = q(jdim-2,k,iend-i,4) 
     .               - 2.*vcont2*si(jdim-2,k,iend,3)
        qj0(k,i,5,4) = q(jdim-2,k,iend-i,5)
c
        bcj(k,i,2) = 0.0
c
   39   continue
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 292 i=ista,iend1
          do 292 k=ksta,kend1
            vj0(k,i,1,3) = vist3d(jdim-1,k,iend-i)
            vj0(k,i,1,4) = vist3d(jdim-2,k,iend-i)
  292     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 202 i=ista,iend1
          do 202 k=ksta,kend1
            tj0(k,i,1,3) = tursav(jdim-1,k,iend-i,1)
            tj0(k,i,2,3) = tursav(jdim-1,k,iend-i,2)
            tj0(k,i,1,4) = tursav(jdim-2,k,iend-i,1)
            tj0(k,i,2,4) = tursav(jdim-2,k,iend-i,2)
  202     continue
        end if
        end if
c
      end if
c
      end if
c
c******************************************************************************
c      k=1 boundary            singular axis - half plane           bctype 1011
c******************************************************************************
c
      if (nface.eq.5) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1003)ista,iend,jsta,jend
      end if
c
c     symmetry in j
c
      if(jsym.gt.0) then
        do 45 i=ista,iend1
        do 45 j=jsta,jend1
c
        wcont1 =  q(jend-j,1,i,2)*sj(jend,1,i,1) +
     .            q(jend-j,1,i,3)*sj(jend,1,i,2) +
     .            q(jend-j,1,i,4)*sj(jend,1,i,3) + sj(jend,1,i,5)
        wcont2 =  q(jend-j,2,i,2)*sj(jend,2,i,1) +
     .            q(jend-j,2,i,3)*sj(jend,2,i,2) +
     .            q(jend-j,2,i,4)*sj(jend,2,i,3) + sj(jend,2,i,5)
c
        qk0(j,i,1,1) = q(jend-j,1,i,1)
        qk0(j,i,2,1) = q(jend-j,1,i,2) - 2.*wcont1*sj(jend,1,i,1)
        qk0(j,i,3,1) = q(jend-j,1,i,3) - 2.*wcont1*sj(jend,1,i,2)
        qk0(j,i,4,1) = q(jend-j,1,i,4) - 2.*wcont1*sj(jend,1,i,3)
        qk0(j,i,5,1) = q(jend-j,1,i,5)
c
        qk0(j,i,1,2) = q(jend-j,2,i,1)
        qk0(j,i,2,2) = q(jend-j,2,i,2) - 2.*wcont2*sj(jend,2,i,1)
        qk0(j,i,3,2) = q(jend-j,2,i,3) - 2.*wcont2*sj(jend,2,i,2)
        qk0(j,i,4,2) = q(jend-j,2,i,4) - 2.*wcont2*sj(jend,2,i,3)
        qk0(j,i,5,2) = q(jend-j,2,i,5)
c
        bck(j,i,1) = 0.0
c
   45   continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 391 i=ista,iend1
          do 391 j=jsta,jend1
            vk0(j,i,1,1) = vist3d(jend-j,1,i)
            vk0(j,i,1,2) = vist3d(jend-j,2,i)
  391     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 301 i=ista,iend1
          do 301 j=jsta,jend1
            tk0(j,i,1,1) = tursav(jend-j,1,i,1)
            tk0(j,i,2,1) = tursav(jend-j,1,i,2)
            tk0(j,i,1,2) = tursav(jend-j,2,i,1)
            tk0(j,i,2,2) = tursav(jend-j,2,i,2)
  301     continue
        end if
        end if
c
      end if
c
c     symmetry in i
c
      if(isym.gt.0) then
        do 46 i=ista,iend1
        do 46 j=jsta,jend1
c
        wcont1 =  q(j,1,iend-i,2)*si(j,1,iend,1) +
     .            q(j,1,iend-i,3)*si(j,1,iend,2) +
     .            q(j,1,iend-i,4)*si(j,1,iend,3) + si(j,1,iend,5)
        wcont2 =  q(j,2,iend-i,2)*si(j,2,iend,1) +
     .            q(j,2,iend-i,3)*si(j,2,iend,2) +
     .            q(j,2,iend-i,4)*si(j,2,iend,3) + si(j,2,iend,5)
c
        qk0(j,i,1,1) = q(j,1,iend-i,1)
        qk0(j,i,2,1) = q(j,1,iend-i,2) - 2.*wcont1*si(j,1,iend,1)
        qk0(j,i,3,1) = q(j,1,iend-i,3) - 2.*wcont1*si(j,1,iend,2)
        qk0(j,i,4,1) = q(j,1,iend-i,4) - 2.*wcont1*si(j,1,iend,3)
        qk0(j,i,5,1) = q(j,1,iend-i,5)
c
        qk0(j,i,1,2) = q(j,2,iend-i,1)
        qk0(j,i,2,2) = q(j,2,iend-i,2) - 2.*wcont2*si(j,2,iend,1)
        qk0(j,i,3,2) = q(j,2,iend-i,3) - 2.*wcont2*si(j,2,iend,2)
        qk0(j,i,4,2) = q(j,2,iend-i,4) - 2.*wcont2*si(j,2,iend,3)
        qk0(j,i,5,2) = q(j,2,iend-i,5)
c
        bck(j,i,1) = 0.0
c
   46   continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 392 i=ista,iend1
          do 392 j=jsta,jend1
            vk0(j,i,1,1) = vist3d(j,1,iend-i)
            vk0(j,i,1,2) = vist3d(j,2,iend-i)
  392     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 302 i=ista,iend1
          do 302 j=jsta,jend1
            tk0(j,i,1,1) = tursav(j,1,iend-i,1)
            tk0(j,i,2,1) = tursav(j,1,iend-i,2)
            tk0(j,i,1,2) = tursav(j,2,iend-i,1)
            tk0(j,i,2,2) = tursav(j,2,iend-i,2)
  302     continue
        end if
        end if
c
      end if
c
      end if
c
c******************************************************************************
c      k=kdim boundary         singular axis - half plane           bctype 1011
c******************************************************************************
c
      if (nface.eq.6) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1004)ista,iend,jsta,jend
      end if
c
c     symmetry in j
c
      if(jsym.gt.0) then
        do 48 i=ista,iend1
        do 48 j=jsta,jend1
c
        wcont1 =  q(jend-j,kdim-1,i,2)*sj(jend,kdim-1,i,1) +
     .            q(jend-j,kdim-1,i,3)*sj(jend,kdim-1,i,2) +
     .            q(jend-j,kdim-1,i,4)*sj(jend,kdim-1,i,3) + 
     .                                 sj(jend,kdim-1,i,5) 
        wcont2 =  q(jend-j,kdim-2,i,2)*sj(jend,kdim-2,i,1) +
     .            q(jend-j,kdim-2,i,3)*sj(jend,kdim-2,i,2) +
     .            q(jend-j,kdim-2,i,4)*sj(jend,kdim-2,i,3) + 
     .                                 sj(jend,kdim-2,i,5)
c
        qk0(j,i,1,3) = q(jend-j,kdim-1,i,1)
        qk0(j,i,2,3) = q(jend-j,kdim-1,i,2) 
     .               - 2.*wcont1*sj(jend,kdim-1,i,1)
        qk0(j,i,3,3) = q(jend-j,kdim-1,i,3) 
     .               - 2.*wcont1*sj(jend,kdim-1,i,2)
        qk0(j,i,4,3) = q(jend-j,kdim-1,i,4) 
     .               - 2.*wcont1*sj(jend,kdim-1,i,3)
        qk0(j,i,5,3) = q(jend-j,kdim-1,i,5)
c
        qk0(j,i,1,4) = q(jend-j,kdim-2,i,1)
        qk0(j,i,2,4) = q(jend-j,kdim-2,i,2) 
     .               - 2.*wcont2*sj(jend,kdim-2,i,1)
        qk0(j,i,3,4) = q(jend-j,kdim-2,i,3) 
     .               - 2.*wcont2*sj(jend,kdim-2,i,2)
        qk0(j,i,4,4) = q(jend-j,kdim-2,i,4) 
     .               - 2.*wcont2*sj(jend,kdim-2,i,3)
        qk0(j,i,5,4) = q(jend-j,kdim-2,i,5)
c
        bck(j,i,2) = 0.0
c
   48   continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 491 i=ista,iend1
          do 491 j=jsta,jend1
            vk0(j,i,1,3) = vist3d(jend-j,kdim-1,i)
            vk0(j,i,1,4) = vist3d(jend-j,kdim-2,i)
  491     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 401 i=ista,iend1
          do 401 j=jsta,jend1
            tk0(j,i,1,3) = tursav(jend-j,kdim-1,i,1)
            tk0(j,i,2,3) = tursav(jend-j,kdim-1,i,2)
            tk0(j,i,1,4) = tursav(jend-j,kdim-2,i,1)
            tk0(j,i,2,4) = tursav(jend-j,kdim-2,i,2)
  401     continue
        end if
        end if
c
      end if
c
c     symmetry in i
c
      if(isym.gt.0) then
        do 49 i=ista,iend1
        do 49 j=jsta,jend1
c
        wcont1 =  q(j,kdim-1,iend-i,2)*si(j,kdim-1,iend,1) +
     .            q(j,kdim-1,iend-i,3)*si(j,kdim-1,iend,2) +
     .            q(j,kdim-1,iend-i,4)*si(j,kdim-1,iend,3) + 
     .                                 si(j,kdim-1,iend,5)
        wcont2 =  q(j,kdim-2,iend-i,2)*si(j,kdim-2,iend,1) +
     .            q(j,kdim-2,iend-i,3)*si(j,kdim-2,iend,2) +
     .            q(j,kdim-2,iend-i,4)*si(j,kdim-2,iend,3) + 
     .                                 si(j,kdim-2,iend,5)
c
        qk0(j,i,1,3) = q(j,kdim-1,iend-i,1)
        qk0(j,i,2,3) = q(j,kdim-1,iend-i,2) 
     .               - 2.*wcont1*si(j,kdim-1,iend,1)
        qk0(j,i,3,3) = q(j,kdim-1,iend-i,3) 
     .               - 2.*wcont1*si(j,kdim-1,iend,2)
        qk0(j,i,4,3) = q(j,kdim-1,iend-i,4) 
     .               - 2.*wcont1*si(j,kdim-1,iend,3)
        qk0(j,i,5,3) = q(j,kdim-1,iend-i,5)
c
        qk0(j,i,1,4) = q(j,kdim-2,iend-i,1)
        qk0(j,i,2,4) = q(j,kdim-2,iend-i,2) 
     .               - 2.*wcont2*si(j,kdim-2,iend,1)
        qk0(j,i,3,4) = q(j,kdim-2,iend-i,3) 
     .               - 2.*wcont2*si(j,kdim-2,iend,2)
        qk0(j,i,4,4) = q(j,kdim-2,iend-i,4) 
     .               - 2.*wcont2*si(j,kdim-2,iend,3)
        qk0(j,i,5,4) = q(j,kdim-2,iend-i,5)
c
        bck(j,i,2) = 0.0
c
   49   continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 492 i=ista,iend1
          do 492 j=jsta,jend1
            vk0(j,i,1,3) = vist3d(j,kdim-1,iend-i)
            vk0(j,i,1,4) = vist3d(j,kdim-2,iend-i)
  492     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 402 i=ista,iend1
          do 402 j=jsta,jend1
            tk0(j,i,1,3) = tursav(j,kdim-1,iend-i,1)
            tk0(j,i,2,3) = tursav(j,kdim-1,iend-i,2)
            tk0(j,i,1,4) = tursav(j,kdim-2,iend-i,1)
            tk0(j,i,2,4) = tursav(j,kdim-2,iend-i,2)
  402     continue
        end if
        end if
c
      end if
c
      end if
c
c******************************************************************************
c      i=1 boundary            singular axis - half plane           bctype 1011
c******************************************************************************
c
      if (nface.eq.1) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1005)jsta,jend,ksta,kend
      end if
c
      i2 = min(2,idim1)
c
c     symmetry in j
c
      if(jsym.gt.0) then
        do 55 k=ksta,kend1
        do 55 j=jsta,jend1
c
        ucont1 =  q(jend-j,k,1,2)*sj(jend,k,1,1) +
     .            q(jend-j,k,1,3)*sj(jend,k,1,2) +
     .            q(jend-j,k,1,4)*sj(jend,k,1,3) + sj(jend,k,1,5)
        ucont2 =  q(jend-j,k,i2,2)*sj(jend,k,i2,1) +
     .            q(jend-j,k,i2,3)*sj(jend,k,i2,2) +
     .            q(jend-j,k,i2,4)*sj(jend,k,i2,3) + sj(jend,k,i2,5)
c
        qi0(j,k,1,1) = q(jend-j,k,1,1)
        qi0(j,k,2,1) = q(jend-j,k,1,2) - 2.*ucont1*sj(jend,k,1,1)
        qi0(j,k,3,1) = q(jend-j,k,1,3) - 2.*ucont1*sj(jend,k,1,2)
        qi0(j,k,4,1) = q(jend-j,k,1,4) - 2.*ucont1*sj(jend,k,1,3)
        qi0(j,k,5,1) = q(jend-j,k,1,5)
c
        qi0(j,k,1,2) = q(jend-j,k,i2,1)
        qi0(j,k,2,2) = q(jend-j,k,i2,2) - 2.*ucont2*sj(jend,k,i2,1)
        qi0(j,k,3,2) = q(jend-j,k,i2,3) - 2.*ucont2*sj(jend,k,i2,2)
        qi0(j,k,4,2) = q(jend-j,k,i2,4) - 2.*ucont2*sj(jend,k,i2,3)
        qi0(j,k,5,2) = q(jend-j,k,i2,5)
c
        bci(j,k,1) = 0.0
c
   55   continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 591 k=ksta,kend1
          do 591 j=jsta,jend1
            vi0(j,k,1,1) = vist3d(jend-j,k,1)
            vi0(j,k,1,2) = vist3d(jend-j,k,i2)
  591     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 501 k=ksta,kend1
          do 501 j=jsta,jend1
            ti0(j,k,1,1) = tursav(jend-j,k,1,1)
            ti0(j,k,2,1) = tursav(jend-j,k,1,2)
            ti0(j,k,1,2) = tursav(jend-j,k,i2,1)
            ti0(j,k,2,2) = tursav(jend-j,k,i2,2)
  501     continue
        end if
        end if
c
      end if
c
c     symmetry in k
c
      if(ksym.gt.0) then
        do 56 k=ksta,kend1
        do 56 j=jsta,jend1
c
        ucont1 =  q(j,kend-k,1,2)*sk(j,kend,1,1) +
     .            q(j,kend-k,1,3)*sk(j,kend,1,2) +
     .            q(j,kend-k,1,4)*sk(j,kend,1,3) + sk(j,kend,1,5)
        ucont2 =  q(j,kend-k,i2,2)*sk(j,kend,i2,1) +
     .            q(j,kend-k,i2,3)*sk(j,kend,i2,2) +
     .            q(j,kend-k,i2,4)*sk(j,kend,i2,3) + sk(j,kend,i2,5)
c
        qi0(j,k,1,1) = q(j,kend-k,1,1)
        qi0(j,k,2,1) = q(j,kend-k,1,2) - 2.*ucont1*sk(j,kend,1,1)
        qi0(j,k,3,1) = q(j,kend-k,1,3) - 2.*ucont1*sk(j,kend,1,2)
        qi0(j,k,4,1) = q(j,kend-k,1,4) - 2.*ucont1*sk(j,kend,1,3)
        qi0(j,k,5,1) = q(j,kend-k,1,5)
c
        qi0(j,k,1,2) = q(j,kend-k,i2,1)
        qi0(j,k,2,2) = q(j,kend-k,i2,2) - 2.*ucont2*sk(j,kend,i2,1)
        qi0(j,k,3,2) = q(j,kend-k,i2,3) - 2.*ucont2*sk(j,kend,i2,2)
        qi0(j,k,4,2) = q(j,kend-k,i2,4) - 2.*ucont2*sk(j,kend,i2,3)
        qi0(j,k,5,2) = q(j,kend-k,i2,5)
c
        bci(j,k,1) = 0.0
c
   56   continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 592 k=ksta,kend1
          do 592 j=jsta,jend1
            vi0(j,k,1,1) = vist3d(j,kend-k,1)
            vi0(j,k,1,2) = vist3d(j,kend-k,i2)
  592     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 502 k=ksta,kend1
          do 502 j=jsta,jend1
            ti0(j,k,1,1) = tursav(j,kend-k,1,1)
            ti0(j,k,2,1) = tursav(j,kend-k,1,2)
            ti0(j,k,1,2) = tursav(j,kend-k,i2,1)
            ti0(j,k,2,2) = tursav(j,kend-k,i2,2)
  502     continue
        end if
        end if
c
      end if
c
      end if
c
c******************************************************************************
c      i=idim boundary         singular axis - half plane           bctype 1011
c******************************************************************************
c
      if (nface.eq.2) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1006)jsta,jend,ksta,kend
      end if
c
      i2 = max(1,idim-2)
c
c     symmetry in j
c
      if(jsym.gt.0) then
        do 58 k=ksta,kend1
        do 58 j=jsta,jend1
c
        ucont1 =  q(jend-j,k,idim-1,2)*sj(jend,k,idim-1,1) +
     .            q(jend-j,k,idim-1,3)*sj(jend,k,idim-1,2) +
     .            q(jend-j,k,idim-1,4)*sj(jend,k,idim-1,3) + 
     .                                 sj(jend,k,idim-1,5)
        ucont2 =  q(jend-j,k,i2,2)*sj(jend,k,i2,1) +
     .            q(jend-j,k,i2,3)*sj(jend,k,i2,2) +
     .            q(jend-j,k,i2,4)*sj(jend,k,i2,3) + 
     .                             sj(jend,k,i2,5)
c
        qi0(j,k,1,3) = q(jend-j,k,idim-1,1)
        qi0(j,k,2,3) = q(jend-j,k,idim-1,2) 
     .               - 2.*ucont1*sj(jend,k,idim-1,1)
        qi0(j,k,3,3) = q(jend-j,k,idim-1,3) 
     .               - 2.*ucont1*sj(jend,k,idim-1,2)
        qi0(j,k,4,3) = q(jend-j,k,idim-1,4) 
     .               - 2.*ucont1*sj(jend,k,idim-1,3)
        qi0(j,k,5,3) = q(jend-j,k,idim-1,5)
c
        qi0(j,k,1,4) = q(jend-j,k,i2,1)
        qi0(j,k,2,4) = q(jend-j,k,i2,2) 
     .               - 2.*ucont2*sj(jend,k,i2,1)
        qi0(j,k,3,4) = q(jend-j,k,i2,3) 
     .               - 2.*ucont2*sj(jend,k,i2,2)
        qi0(j,k,4,4) = q(jend-j,k,i2,4) 
     .               - 2.*ucont2*sj(jend,k,i2,3)
        qi0(j,k,5,4) = q(jend-j,k,i2,5)
c
        bci(j,k,2) = 0.0
c
   58   continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 691 k=ksta,kend1
          do 691 j=jsta,jend1
            vi0(j,k,1,3) = vist3d(jend-j,k,idim-1)
            vi0(j,k,1,4) = vist3d(jend-j,k,i2)
  691     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 601 k=ksta,kend1
          do 601 j=jsta,jend1
            ti0(j,k,1,3) = tursav(jend-j,k,idim-1,1)
            ti0(j,k,2,3) = tursav(jend-j,k,idim-1,2)
            ti0(j,k,1,4) = tursav(jend-j,k,i2,1)
            ti0(j,k,2,4) = tursav(jend-j,k,i2,2)
  601     continue
        end if
        end if
c
      end if
c
c     symmetry in k
c
      if(ksym.gt.0) then
        do 59 k=ksta,kend1
        do 59 j=jsta,jend1
c
        ucont1 =  q(j,kend-k,idim-1,2)*sk(j,kend,idim-1,1) +
     .            q(j,kend-k,idim-1,3)*sk(j,kend,idim-1,2) +
     .            q(j,kend-k,idim-1,4)*sk(j,kend,idim-1,3) + 
     .                                 sk(j,kend,idim-1,5)
        ucont2 =  q(j,kend-k,i2,2)*sk(j,kend,i2,1) +
     .            q(j,kend-k,i2,3)*sk(j,kend,i2,2) +
     .            q(j,kend-k,i2,4)*sk(j,kend,i2,3) + 
     .                             sk(j,kend,i2,5)
c
        qi0(j,k,1,3) = q(j,kend-k,idim-1,1)
        qi0(j,k,2,3) = q(j,kend-k,idim-1,2) 
     .               - 2.*ucont1*sk(j,kend,idim-1,1)
        qi0(j,k,3,3) = q(j,kend-k,idim-1,3) 
     .               - 2.*ucont1*sk(j,kend,idim-1,2)
        qi0(j,k,4,3) = q(j,kend-k,idim-1,4) 
     .               - 2.*ucont1*sk(j,kend,idim-1,3)
        qi0(j,k,5,3) = q(j,kend-k,idim-1,5)
c
        qi0(j,k,1,4) = q(j,kend-k,i2,1)
        qi0(j,k,2,4) = q(j,kend-k,i2,2) 
     .               - 2.*ucont2*sk(j,kend,i2,1)
        qi0(j,k,3,4) = q(j,kend-k,i2,3) 
     .               - 2.*ucont2*sk(j,kend,i2,2)
        qi0(j,k,4,4) = q(j,kend-k,i2,4) 
     .               - 2.*ucont2*sk(j,kend,i2,3)
        qi0(j,k,5,4) = q(j,kend-k,i2,5)
c
        bci(j,k,2) = 0.0
c
   59   continue
c
        if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
          do 692 k=ksta,kend1
          do 692 j=jsta,jend1
            vi0(j,k,1,3) = vist3d(j,kend-k,idim-1)
            vi0(j,k,1,4) = vist3d(j,kend-k,i2)
  692     continue
        end if
c       only need to do advanced model turbulence B.C.s on finest grid
        if (level .ge. lglobal) then
        if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
          do 602 k=ksta,kend1
          do 602 j=jsta,jend1
            ti0(j,k,1,3) = tursav(j,kend-k,idim-1,1)
            ti0(j,k,2,3) = tursav(j,kend-k,idim-1,2)
            ti0(j,k,1,4) = tursav(j,kend-k,i2,1)
            ti0(j,k,2,4) = tursav(j,kend-k,i2,2)
  602     continue
        end if
        end if
c
      end if
      end if
c
 1001 format(' ','  j=   1  singular axis - half plane     type 1011',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
 1002 format(' ','  j=jdim  singular axis - half plane     type 1011',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
 1003 format(' ','  k=   1  singular axis - half plane     type 1011',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
 1004 format(' ','  k=kdim  singular axis - half plane     type 1011',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
 1005 format(' ','  i=   1  singular axis - half plane     type 1011',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
 1006 format(' ','  i=idim  singular axis - half plane     type 1011',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
c
      return
      end
