      subroutine bc1008(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,iuns,nou,bou,nbuf,
     .                  ibufdim)
c***********************************************************************
c     Purpose:  Set tunnel inflow boundary conditions
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
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10,tur20
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
c            * standard boundary condition bctype=1008 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary            tunnel inflow                        bctype 1008
c******************************************************************************
c
      if (nface.eq.3) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1001)ista,iend,ksta,kend
      end if
c
      ent0 = p0/rho0**gamma
      gami = 1.0/gamma
c
      do 3031 i=ista,iend1
      do 3031 k=ksta,kend1
c
      qj0(k,i,5,1) = q(1,k,i,5)
      qj0(k,i,1,1) = (qj0(k,i,5,1)/ent0)**gami
      qj0(k,i,2,1) = sqrt(2.0*(h0-qj0(k,i,5,1)/qj0(k,i,1,1)/gm1g))
      qj0(k,i,3,1) = 0.0
      qj0(k,i,4,1) = 0.0
c
      qj0(k,i,1,2) = qj0(k,i,1,1)
      qj0(k,i,2,2) = qj0(k,i,2,1)
      qj0(k,i,3,2) = qj0(k,i,3,1)
      qj0(k,i,4,2) = qj0(k,i,4,1)
      qj0(k,i,5,2) = qj0(k,i,5,1)
c
      bcj(k,i,1) = 0.0
c
 3031 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,2) = 0.0
  191   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 101 i=ista,iend1
        do 101 k=ksta,kend1
          tj0(k,i,1,1) = tur10
          tj0(k,i,2,1) = tur20
          tj0(k,i,1,2) = tur10
          tj0(k,i,2,2) = tur20
  101   continue
      end if
      end if
c
      end if
c
c******************************************************************************
c      j=jdim boundary         tunnel inflow                        bctype 1008
c******************************************************************************
c
      if (nface.eq.4) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1002)ista,iend,ksta,kend
      end if
c
      ent0 = p0/rho0**gamma
      gami = 1.0/gamma
c
      do 3131 i=ista,iend1
      do 3131 k=ksta,kend1
c
      qj0(k,i,5,3) = q(jdim1,k,i,5)
      qj0(k,i,1,3) = (qj0(k,i,5,3)/ent0)**gami
      qj0(k,i,2,3) = sqrt(2.0*(h0-qj0(k,i,5,3)/qj0(k,i,1,3)/gm1g))
      qj0(k,i,3,3) = 0.0
      qj0(k,i,4,3) = 0.0
c
      qj0(k,i,1,4) = qj0(k,i,1,3)
      qj0(k,i,2,4) = qj0(k,i,2,3)
      qj0(k,i,3,4) = qj0(k,i,3,3)
      qj0(k,i,4,4) = qj0(k,i,4,3)
      qj0(k,i,5,4) = qj0(k,i,5,3)
c
      bcj(k,i,2) = 0.0
c
 3131 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
          vj0(k,i,1,3) = vist3d(jdim-1,k,i)
          vj0(k,i,1,4) = 0.0
  291   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 201 i=ista,iend1
        do 201 k=ksta,kend1
          tj0(k,i,1,3) = tur10
          tj0(k,i,2,3) = tur20
          tj0(k,i,1,4) = tur10
          tj0(k,i,2,4) = tur20
  201   continue
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=1 boundary            tunnel inflow                        bctype 1008
c******************************************************************************
c
      if (nface.eq.5) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1003)ista,iend,jsta,jend
      end if
c
      ent0 = p0/rho0**gamma
      gami = 1.0/gamma
c
      do 3231 i=ista,iend1
      do 3231 j=jsta,jend1
c
      qk0(j,i,5,1) = q(j,1,i,5)
      qk0(j,i,1,1) = (qk0(j,i,5,1)/ent0)**gami
      qk0(j,i,2,1) = sqrt(2.0*(h0-qk0(j,i,5,1)/qk0(j,i,1,1)/gm1g))
      qk0(j,i,3,1) = 0.0
      qk0(j,i,4,1) = 0.0
c
      qk0(j,i,1,2) = qk0(j,i,1,1)
      qk0(j,i,2,2) = qk0(j,i,2,1)
      qk0(j,i,3,2) = qk0(j,i,3,1)
      qk0(j,i,4,2) = qk0(j,i,4,1)
      qk0(j,i,5,2) = qk0(j,i,5,1)
c
      bck(j,i,1) = 0.0
c
 3231 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = 0.0
  391   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 301 i=ista,iend1
        do 301 j=jsta,jend1
          tk0(j,i,1,1) = tur10
          tk0(j,i,2,1) = tur20
          tk0(j,i,1,2) = tur10
          tk0(j,i,2,2) = tur20
  301   continue
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=kdim boundary         tunnel inflow                        bctype 1008
c******************************************************************************
c
      if (nface.eq.6) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1004)ista,iend,jsta,jend
      end if
c
      ent0 = p0/rho0**gamma
      gami = 1.0/gamma
c
      do 3331 i=ista,iend1
      do 3331 j=jsta,jend1
c
      qk0(j,i,5,3) = q(j,kdim1,i,5)
      qk0(j,i,1,3) = (qk0(j,i,5,3)/ent0)**gami
      qk0(j,i,2,3) = sqrt(2.0*(h0-qk0(j,i,5,3)/qk0(j,i,1,3)/gm1g))
      qk0(j,i,3,3) = 0.0
      qk0(j,i,4,3) = 0.0
c
      qk0(j,i,1,4) = qk0(j,i,1,3)
      qk0(j,i,2,4) = qk0(j,i,2,3)
      qk0(j,i,3,4) = qk0(j,i,3,3)
      qk0(j,i,4,4) = qk0(j,i,4,3)
      qk0(j,i,5,4) = qk0(j,i,5,3)
c
      bck(j,i,2) = 0.0
c
 3331 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
          vk0(j,i,1,3) = vist3d(j,kdim-1,i)
          vk0(j,i,1,4) = 0.0
  491   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 401 i=ista,iend1
        do 401 j=jsta,jend1
          tk0(j,i,1,3) = tur10
          tk0(j,i,2,3) = tur20
          tk0(j,i,1,4) = tur10
          tk0(j,i,2,4) = tur20
  401   continue
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=1 boundary            tunnel inflow                        bctype 1008
c******************************************************************************
c
      if (nface.eq.1) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1005)jsta,jend,ksta,kend
      end if
c
      ent0 = p0/rho0**gamma
      gami = 1.0/gamma
c
      do 3431 j=jsta,jend1
      do 3431 k=ksta,kend1
c
      qi0(j,k,5,1) = q(j,k,1,5)
      qi0(j,k,1,1) = (qi0(j,k,5,1)/ent0)**gami
      qi0(j,k,2,1) = sqrt(2.0*(h0-qi0(j,k,5,1)/qi0(j,k,1,1)/gm1g))
      qi0(j,k,3,1) = 0.0
      qi0(j,k,4,1) = 0.0
c
      qi0(j,k,1,2) = qi0(j,k,1,1)
      qi0(j,k,2,2) = qi0(j,k,2,1)
      qi0(j,k,3,2) = qi0(j,k,3,1)
      qi0(j,k,4,2) = qi0(j,k,4,1)
      qi0(j,k,5,2) = qi0(j,k,5,1)
c
      bci(j,k,1) = 0.0
c
 3431 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = 0.0
  591   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 501 k=ksta,kend1
        do 501 j=jsta,jend1
          ti0(j,k,1,1) = tur10
          ti0(j,k,2,1) = tur20
          ti0(j,k,1,2) = tur10
          ti0(j,k,2,2) = tur20
  501   continue
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=idim boundary          tunnel inflow                       bctype 1008
c******************************************************************************
c
      if (nface.eq.2) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1006)jsta,jend,ksta,kend
      end if
c
      ent0 = p0/rho0**gamma
      gami = 1.0/gamma
c
      do 3531 j=jsta,jend1
      do 3531 k=ksta,kend1
c
      qi0(j,k,5,3) = q(j,k,idim1,5)
      qi0(j,k,1,3) = (qi0(j,k,5,3)/ent0)**gami
      qi0(j,k,2,3) = sqrt(2.0*(h0-qi0(j,k,5,3)/qi0(j,k,1,3)/gm1g))
      qi0(j,k,3,3) = 0.0
      qi0(j,k,4,3) = 0.0
c
      qi0(j,k,1,4) = qi0(j,k,1,1)
      qi0(j,k,2,4) = qi0(j,k,2,1)
      qi0(j,k,3,4) = qi0(j,k,3,1)
      qi0(j,k,4,4) = qi0(j,k,4,1)
      qi0(j,k,5,4) = qi0(j,k,5,1)
c
      bci(j,k,2) = 0.0
c
 3531 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
          vi0(j,k,1,3) = vist3d(j,k,idim-1)
          vi0(j,k,1,4) = 0.0
  691   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do 601 k=ksta,kend1
        do 601 j=jsta,jend1
          ti0(j,k,1,3) = tur10
          ti0(j,k,2,3) = tur20
          ti0(j,k,1,4) = tur10
          ti0(j,k,2,4) = tur20
  601   continue
      end if
      end if
      end if
c
 1001 format(' ','  j=   1  tunnel inflow                  type 1008',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
 1002 format(' ','  j=jdim  tunnel inflow                  type 1008',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
 1003 format(' ','  k=   1  tunnel inflow                  type 1008',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
 1004 format(' ','  k=kdim  tunnel inflow                  type 1008',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
 1005 format(' ','  i=   1  tunnel inflow                  type 1008',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
 1006 format(' ','  i=idim  tunnel inflow                  type 1008',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
c
      return
      end
