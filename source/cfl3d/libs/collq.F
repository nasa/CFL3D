      subroutine collq(q,qq,vol,volc,jdim,kdim,idim,jj2,kk2,ii2,
     .                 res,qr,qw,vistf,vistc,nbl,nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Restrict q (the primative variables) with a volume-
c     weighted interpolation and residuals to coarser meshes.  Also
c     restrict turbulent eddy viscosity in the case of turbulent flows
c     to coarser meshes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5),qq(jj2,kk2,ii2,5),
     .          qr(jj2,kk2,ii2-1,5),qw(jdim,kdim,idim,5)
      dimension vol(jdim,kdim,idim-1),volc(jj2,kk2,ii2-1),
     .          res(jdim,kdim,idim-1,5)
      dimension vistf(jdim,kdim,idim),vistc(jj2,kk2,ii2)
c
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
c
c      restrict q(volume-weighted) and r to coarser mesh
c      restrict viscosity to coarser mesh
c
c      jdim,kdim,idim  finer mesh
c      jj2,kk2,ii2  coarser mesh
c
      nbl1  = nbl+1
      if (isklton.gt.0)  then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),7) nbl,nbl1
      end if
    7 format(1x,45hrestricting variables and residual from finer,
     .       6h block,i4,1x,16hto coarser block,i4)
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
      jjl   = jj2-1
      kkl   = kk2-1
      iil   = ii2-1
      n     = jj2*kk2
      do 6 i=1,iil
      do 6 l=1,5
cdir$ ivdep
      do 1000 izz=1,n
      qr(izz,1,i,l) = 0.0
 1000 continue
    6 continue
      do 60 n=1,5
      nv    = jdim*kdim
      do 30 i=1,idim1
cdir$ ivdep
      do 1001 izz=1,nv
      qw(izz,1,i,n) = q(izz,1,i,n)*vol(izz,1,i)
 1001 continue
   30 continue
cdir$ ivdep
      do 1002 izz=1,nv
      qw(izz,1,idim,n) = q(izz,1,idim,n)*vol(izz,1,idim1)
 1002 continue
      if (idim.gt.2) then
      kk = 0
      do 40 k=1,kdim1,2
      kk = kk+1
      jj = 0
      do 40 j=1,jdim1,2
      jj = jj+1
      ii = 0
      do 40 i=1,idim1,2
      ii = ii+1
      qq(jj,kk,ii,n)  = (qw(j,k,i,n)      +qw(j,k,i+1,n)
     .                  +qw(j+1,k,i,n)    +qw(j+1,k,i+1,n)
     .                  +qw(j,k+1,i,n)    +qw(j,k+1,i+1,n)
     .                  +qw(j+1,k+1,i,n)  +qw(j+1,k+1,i+1,n))/
     .                   volc(jj,kk,ii)
      qr(jj,kk,ii,n)  = (res(j,k,i,n)     +res(j,k,i+1,n)
     .                  +res(j+1,k,i,n)   +res(j+1,k,i+1,n)
     .                  +res(j,k+1,i,n)   +res(j,k+1,i+1,n)
     .                  +res(j+1,k+1,i,n) +res(j+1,k+1,i+1,n))
      if ((ivisc(1).gt.1 .and. n.eq.1) .or. (ivisc(2).gt.1 .and. n.eq.1)
     .   .or. (ivisc(3).gt.1 .and. n.eq.1)) then
         vistc(jj,kk,ii) = .125e0*(vistf(j,k,i)     +vistf(j,k,i+1)
     .                            +vistf(j+1,k,i)   +vistf(j+1,k,i+1)
     .                            +vistf(j,k+1,i)   +vistf(j,k+1,i+1)
     .                            +vistf(j+1,k+1,i) +vistf(j+1,k+1,i+1))
      end if
   40 continue
c
      else
c
      ii = 1
      i  = 1
      kk = 0
      do 404 k=1,kdim1,2
      kk = kk+1
      jj = 0
      do 404 j=1,jdim1,2
      jj = jj+1
      qq(jj,kk,ii,n)  = (qw(j,k,i,n)    +qw(j+1,k,i,n)    
     .                  +qw(j,k+1,i,n)  +qw(j+1,k+1,i,n))/
     .                   volc(jj,kk,ii)
      qr(jj,kk,ii,n)  = (res(j,k,i,n)   +res(j+1,k,i,n)    
     .                  +res(j,k+1,i,n) +res(j+1,k+1,i,n))
      if ((ivisc(1).gt.1 .and. n.eq.1) .or. (ivisc(2).gt.1 .and. n.eq.1)
     .   .or. (ivisc(3).gt.1 .and. n.eq.1)) then
         vistc(jj,kk,ii) = .25e0*(vistf(j,k,i)      +vistf(j+1,k,i)
     .                           +vistf(j,k+1,i)    +vistf(j+1,k+1,i))
      end if
  404 continue
      end if
   60 continue
      call fill(jj2,kk2,ii2,qq,5)
      return
      end
