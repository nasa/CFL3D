      subroutine abckz(i,npl,jdim,kdim,idim,ak,bk,ck,blank)
c***********************************************************************
c     Purpose:  Zero out the LHS matrix element with the help of the 
c     blank array.  For a point with blank=0, all elements of matrices
c     a and c become zero. Only diagonal elements matrix b is changed
c     to 1.0 for K-implicit ; J-sweep.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension ak(jdim-1,npl,kdim-1,5,5),bk(jdim-1,npl,kdim-1,5,5),
     .          ck(jdim-1,npl,kdim-1,5,5),blank(jdim,kdim,idim)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
c
      do 2750 m=1,5 
      do 2750 n=1,5 
      do 2750 k=1,kdim1 
      do 2750 lpl=1,npl 
      ii = i+lpl-1
cdir$ ivdep
      do 2765 j=1,jdim1 
      ak(j,lpl,k,m,n)=ak(j,lpl,k,m,n)*blank(j,k,ii)
      ck(j,lpl,k,m,n)=ck(j,lpl,k,m,n)*blank(j,k,ii)
 2765 continue
 2750 continue
c
      do 2850 m=1,5 
      do 2850 k=1,kdim1 
      do 2850 lpl=1,npl 
      ii = i+lpl-1
cdir$ ivdep
      do 2865 j=1,jdim1 
      bk(j,lpl,k,m,m) =  (bk(j,lpl,k,m,m)*blank(j,k,ii))
     .                  +(1.0-blank(j,k,ii)) 
 2865 continue
 2850 continue
      return
      end 
