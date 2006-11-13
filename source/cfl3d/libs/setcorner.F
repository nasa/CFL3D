      subroutine setcorner(j,k,i,xnm1,ynm1,znm1,dx,dy,dz,jdim,kdim,
     .                     idim,wk,nsurf,iflag,ivert,islavept,nslave,
     .                     nou,bou,ibufdim,nbuf,myid,ibl,nmaster,iseq)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute dispacements for subgrid vertex points to use
c     for transfinite interpolation. Points off the body are given 
c     displacements that are reduced via a decay function from the 
c     displacements of the closest point on the surface. The decay
c     function is the one proposed by Peter Hartwich, with minor
c     modifications to the parameters, and to allow for more than 
c     one master point to influence the slave point (originally, only
c     the nearest master point was used; now, the nmaster nearest master
c     points can be used)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension xnm1(jdim,kdim,idim),ynm1(jdim,kdim,idim),
     .          znm1(jdim,kdim,idim)
      dimension dx(jdim,kdim,idim),dy(jdim,kdim,idim),
     .          dz(jdim,kdim,idim)
      dimension wk(9*nsurf),islavept(nslave,nmaster,5),ibl(nsurf)
c
      common /zero/ iexp
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt 
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
c
c
      return
      end
