      subroutine outgrd(nf,ib,x,y,z,ni,nj,nk,i0,i1,isk,j0,j1,jsk,
     .                  k0,k1,ksk,ip3dgrd)
c
      dimension x(ni,nj,nk), y(ni,nj,nk), z(ni,nj,nk)
c
      if (ip3dgrd.gt.0) then      
         if (ib.eq.0) then
            write(nf,*)          (((x(i,j,k),i=i0,i1,isk),
     .                                       j=j0,j1,jsk),
     .                                       k=k0,k1,ksk),
     .                           (((y(i,j,k),i=i0,i1,isk),
     .                                       j=j0,j1,jsk),
     .                                       k=k0,k1,ksk),
     .                           (((z(i,j,k),i=i0,i1,isk),
     .                                       j=j0,j1,jsk),
     .                                       k=k0,k1,ksk)
         else
            write(nf)            (((x(i,j,k),i=i0,i1,isk),
     .                                       j=j0,j1,jsk),
     .                                       k=k0,k1,ksk),
     .                           (((y(i,j,k),i=i0,i1,isk),
     .                                       j=j0,j1,jsk),
     .                                       k=k0,k1,ksk),
     .                           (((z(i,j,k),i=i0,i1,isk),
     .                                       j=j0,j1,jsk),
     .                                       k=k0,k1,ksk)
         end if
      else
         if (ib.eq.0) then
            write(nf,*)          (((x(i,j,k),j=j0,j1,jsk),
     .                                       k=k0,k1,ksk),
     .                                       i=i0,i1,isk),
     .                           (((y(i,j,k),j=j0,j1,jsk),
     .                                       k=k0,k1,ksk),
     .                                       i=i0,i1,isk),
     .                           (((z(i,j,k),j=j0,j1,jsk),
     .                                       k=k0,k1,ksk),
     .                                       i=i0,i1,isk)
         else
            write(nf)            (((x(i,j,k),j=j0,j1,jsk),
     .                                       k=k0,k1,ksk),
     .                                       i=i0,i1,isk),
     .                           (((y(i,j,k),j=j0,j1,jsk),
     .                                       k=k0,k1,ksk),
     .                                       i=i0,i1,isk),
     .                           (((z(i,j,k),j=j0,j1,jsk),
     .                                       k=k0,k1,ksk),
     .                                       i=i0,i1,isk)
         end if
      end if
c
      return
      end
