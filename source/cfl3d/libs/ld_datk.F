      subroutine ld_datk(data,jdim,kdim,idim,datak,ldim,mp,np,
     .                   ista,iend,jsta,jend,ksta,kend,nou,bou,nbuf,
     .                   ibufdim,myid)
c***********************************************************************
c     Purpose:  Load data from full 3d array to smaller work array for
c               efficient message passing
c***********************************************************************
c
c     Description of variables:
c       data              = full 3D array of data array
c       jdim,kdim,idim    = dimensions of 3D data array
c       datak             = np-plane subset of data in k direction
c       mp                = array containing plane indicies to store
c       np                = number of planes of data needed (up to 4)
c       ista,iend,etc.    = starting and ending values for the
c                           boundary segment of data array to be loaded
c                           into the subset arrays
c       ldim              = number of variables in data array; e.g.
c                           ldim = 5 for q data, ldim = 2 for turb.
c                           variables, ldim = 1 for eddy viscosity,
c                           ldim = 3 for grid data
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension data(jdim,kdim,idim,ldim),
     .          datak(jdim,np,idim,ldim)
      dimension mp(4)
c
      if (np.gt.4) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)'stopping...increase dimension of mp ',
     .                'array in ld_datk and its calling routine'
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c     load data values into dataX array:
c
      do i=ista,iend
         do j=jsta,jend
           do m=1,np
              do l=1,ldim
                 datak(j,m,i,l)=data(j,mp(m),i,l)
              end do
            end do
         end do
      end do
c
      return
      end
