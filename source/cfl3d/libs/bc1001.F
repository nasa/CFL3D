      subroutine bc1001(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,iuns,nou,bou,nbuf,
     .                  ibufdim)
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
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=1001 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary            symmetry plane                       bctype 1001
c******************************************************************************
c
      if (nface.eq.3) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1001)ista,iend,ksta,kend
      end if
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
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 101 i=ista,iend1
        do 101 k=ksta,kend1
          tj0(k,i,1,1) = tursav(1,k,i,1)
          tj0(k,i,2,1) = tursav(1,k,i,2)
          tj0(k,i,1,2) = tursav(2,k,i,1)
          tj0(k,i,2,2) = tursav(2,k,i,2)
  101   continue
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
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1002)ista,iend,ksta,kend
      end if
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
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 201 i=ista,iend1
        do 201 k=ksta,kend1
          tj0(k,i,1,3) = tursav(jdim-1,k,i,1)
          tj0(k,i,2,3) = tursav(jdim-1,k,i,2)
          tj0(k,i,1,4) = tursav(jdim-2,k,i,1)
          tj0(k,i,2,4) = tursav(jdim-2,k,i,2)
  201   continue
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
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1003)ista,iend,jsta,jend
      end if
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
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 301 i=ista,iend1
        do 301 j=jsta,jend1
          tk0(j,i,1,1) = tursav(j,1,i,1)
          tk0(j,i,2,1) = tursav(j,1,i,2)
          tk0(j,i,1,2) = tursav(j,2,i,1)
          tk0(j,i,2,2) = tursav(j,2,i,2)
  301   continue
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
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1004)ista,iend,jsta,jend
      end if
c
      do 49 i=ista,iend1
      do 49 j=jsta,jend1
c
      wcont1 =  q(j,kdim-1,i,2)*sk(j,kdim,i,1) +
     .          q(j,kdim-1,i,3)*sk(j,kdim,i,2) +
     .          q(j,kdim-1,i,4)*sk(j,kdim,i,3) + sk(j,kdim,i,5) 
      wcont2 =  q(j,kdim-2,i,2)*sk(j,kdim,i,1) +
     .          q(j,kdim-2,i,3)*sk(j,kdim,i,2) +
     .          q(j,kdim-2,i,4)*sk(j,kdim,i,3) + sk(j,kdim,i,5)
c
      qk0(j,i,1,3) = q(j,kdim-1,i,1)
      qk0(j,i,2,3) = q(j,kdim-1,i,2) - 2.*wcont1*sk(j,kdim,i,1)
      qk0(j,i,3,3) = q(j,kdim-1,i,3) - 2.*wcont1*sk(j,kdim,i,2)
      qk0(j,i,4,3) = q(j,kdim-1,i,4) - 2.*wcont1*sk(j,kdim,i,3)
      qk0(j,i,5,3) = q(j,kdim-1,i,5)
c
      qk0(j,i,1,4) = q(j,kdim-2,i,1)
      qk0(j,i,2,4) = q(j,kdim-2,i,2) - 2.*wcont2*sk(j,kdim,i,1)
      qk0(j,i,3,4) = q(j,kdim-2,i,3) - 2.*wcont2*sk(j,kdim,i,2)
      qk0(j,i,4,4) = q(j,kdim-2,i,4) - 2.*wcont2*sk(j,kdim,i,3)
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
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 401 i=ista,iend1
        do 401 j=jsta,jend1
          tk0(j,i,1,3) = tursav(j,kdim-1,i,1)
          tk0(j,i,2,3) = tursav(j,kdim-1,i,2)
          tk0(j,i,1,4) = tursav(j,kdim-2,i,1)
          tk0(j,i,2,4) = tursav(j,kdim-2,i,2)
  401   continue
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
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1005)jsta,jend,ksta,kend
      end if
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
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 501 k=ksta,kend1
        do 501 j=jsta,jend1
          ti0(j,k,1,1) = tursav(j,k,1,1)
          ti0(j,k,2,1) = tursav(j,k,1,2)
          ti0(j,k,1,2) = tursav(j,k,i2,1)
          ti0(j,k,2,2) = tursav(j,k,i2,2)
  501   continue
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
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1006)jsta,jend,ksta,kend
      end if
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
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 601 k=ksta,kend1
        do 601 j=jsta,jend1
          ti0(j,k,1,3) = tursav(j,k,idim-1,1)
          ti0(j,k,2,3) = tursav(j,k,idim-1,2)
          ti0(j,k,1,4) = tursav(j,k,i2,1)
          ti0(j,k,2,4) = tursav(j,k,i2,2)
  601   continue
      end if
      end if
      end if
c
 1001 format(' ','  j=   1  symmetry plane                 type 1001',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
 1002 format(' ','  j=jdim  symmetry plane                 type 1001',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
 1003 format(' ','  k=   1  symmetry plane                 type 1001',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
 1004 format(' ','  k=kdim  symmetry plane                 type 1001',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
 1005 format(' ','  i=   1  symmetry plane                 type 1001',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
 1006 format(' ','  i=idim  symmetry plane                 type 1001',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
c
      return
      end
