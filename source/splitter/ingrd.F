      subroutine ingrd(nf,ib,x,y,z,ni,nj,nk,ip3dgrd)
c
      dimension x(ni,nj,nk), y(ni,nj,nk), z(ni,nj,nk)
c
      if (ip3dgrd.gt.0) then      
         if (ib.eq.0) then
            read(nf,*)           (((x(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk),
     .                           (((y(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk),
     .                           (((z(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk)
         else
            read(nf)             (((x(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk),
     .                           (((y(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk),
     .                           (((z(i,j,k),i=1,ni),
     .                                       j=1,nj),
     .                                       k=1,nk)
         end if
      else
         if (ib.eq.0) then
            read(nf,*)           (((x(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni),
     .                           (((y(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni),
     .                           (((z(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni)
         else
            read(nf)             (((x(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni),
     .                           (((y(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni),
     .                           (((z(i,j,k),j=1,nj),
     .                                       k=1,nk),
     .                                       i=1,ni)
         end if
      end if
c
      return
      end
