      subroutine resnonin(nbl,jdim,kdim,idim,q,x,y,z,sj,sk,si,vol,res,
     .                    nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Compute the residual contributions to the right-hand-side
c              arising from source terms in a steadily rotating 
c              noninertial reference frame
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim*kdim*idim,5)
      dimension x(jdim*kdim*idim),y(jdim*kdim*idim),z(jdim*kdim*idim)
      dimension si(jdim*kdim*idim,5),sj(jdim*kdim*(idim-1),5),
     .          sk(jdim*kdim*(idim-1),5)
      dimension vol(jdim*kdim*(idim-1))
      dimension res(jdim*kdim*(idim-1),5)
 
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10,tur20
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /sklton/ isklton
 
      if (isklton.gt.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),101)
 101        format(34h   adding NONINERTIAL source terms)
      end if
 
c  rename xrotrate, yrotrate, zrotrate -> wx, wy, wz to save typing
                  
      wx = xrotrate
      wy = yrotrate
      wz = zrotrate
 
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
 
c*********************
c     centripetal part
c*********************
 
c compute Omega x Uinf
c  (constant with respect to solution or q) 
 
      OmegaxUx =  ( wy * qiv(4) - wz * qiv(3) )      
      OmegaxUy =  ( wz * qiv(2) - wx * qiv(4) )
      OmegaxUz =  ( wx * qiv(3) - wy * qiv(2) )
 
      do i = 1,idim1
        do k = 1,kdim1
          do j = 1,jdim1
 
c Find the indices corresponding to the 8 corners of a volume
 
          ind000 = (i-1) * jdim * kdim + (k-1) * jdim + (j-1) + 1
          ind001 = (i-1) * jdim * kdim + (k-1) * jdim + (j  ) + 1
          ind010 = (i-1) * jdim * kdim + (k  ) * jdim + (j-1) + 1
          ind011 = (i-1) * jdim * kdim + (k  ) * jdim + (j  ) + 1
          ind100 = (i  ) * jdim * kdim + (k-1) * jdim + (j-1) + 1
          ind101 = (i  ) * jdim * kdim + (k-1) * jdim + (j  ) + 1
          ind110 = (i  ) * jdim * kdim + (k  ) * jdim + (j-1) + 1
          ind111 = (i  ) * jdim * kdim + (k  ) * jdim + (j  ) + 1
 
c  Coriolis part of pseudo-acceleration
 
          corx = ( wy * q(ind000,4) - wz * q(ind000,3) )
          cory = ( wz * q(ind000,2) - wx * q(ind000,4) )
          corz = ( wx * q(ind000,3) - wy * q(ind000,2) )
 
c  increment residual (in conserved vars.)
 
          res(ind000,2) = res(ind000,2) + q(ind000,1)*corx*vol(ind000)
          res(ind000,3) = res(ind000,3) + q(ind000,1)*cory*vol(ind000)
          res(ind000,4) = res(ind000,4) + q(ind000,1)*corz*vol(ind000)
 
          enddo
        enddo
      enddo

      return 
      end
