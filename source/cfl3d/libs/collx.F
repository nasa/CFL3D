      subroutine collx(x,y,z,xx,yy,zz,jdim,kdim,idim,jj2,kk2,ii2)
c***********************************************************************
c     Purpose:  Restrict x, y, and z values to coarser meshes.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension xx(jj2,kk2,ii2),yy(jj2,kk2,ii2),zz(jj2,kk2,ii2)
c
      ii   = 0
      iinc = 2
      if (idim.eq.2) iinc = 1
      do 10 i=1,idim,iinc
      ii   = ii+1
      kk   = 0
      do 10 k=1,kdim,2
      kk   = kk+1
      jj   = 0
      do 10 j=1,jdim,2
      jj   = jj+1
      xx(jj,kk,ii) = x(j,k,i)
      yy(jj,kk,ii) = y(j,k,i)
      zz(jj,kk,ii) = z(j,k,i)
   10 continue
      return
      end
