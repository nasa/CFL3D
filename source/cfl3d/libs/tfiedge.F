      subroutine tfiedge(idim,jdim,kdim,x,y,z,i1,i2,j1,j2,k1,k2,
     .                   arci,arcj,arck,nou,bou,nbuf,ibufdim,myid)
c***********************************************************************
c     Purpose: compute transfinite interpolation on block edges, using
c              arc-length blending functions
c
c     this subroutine expects one and only pair of (i1,i2), (j1,j2), or
c     (k1,k2) to contain distinct indicies; the distinct indicies
c     determine which edge the TFI is carried out on.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension arci(jdim,kdim,idim),arcj(jdim,kdim,idim),
     .          arck(jdim,kdim,idim)
c
      dimension psi(2)
c
c     check that one and only only one of the pairs i1,i2 j1,j2 k1,k2
c     varies
c
      nvar = 0
      if (i1.ne.i2) nvar = nvar + 1
      if (j1.ne.j2) nvar = nvar + 1
      if (k1.ne.k2) nvar = nvar + 1
      if (nvar.eq.0 .or. nvar.gt.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'('' stopping...input error to'',
     .   '' subroutine tfiedge'')') 
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      if (i1.ne.i2) then
         j = j1
         k = k1
         do i=i1,i2
            eta = (arci(j,k,i) - arci(j,k,i1))
     .          / (arci(j,k,i2)- arci(j,k,i1))
            psi(1) = eta
            psi(2) = 1.-eta
            x(j,k,i) = psi(2)*x(j,k,i1)
     .               + psi(1)*x(j,k,i2)
            y(j,k,i) = psi(2)*y(j,k,i1)
     .               + psi(1)*y(j,k,i2)
            z(j,k,i) = psi(2)*z(j,k,i1)
     .               + psi(1)*z(j,k,i2)
         end do
      end if
c
      if (j1.ne.j2) then
         i = i1
         k = k1
         do j=j1,j2
            eta = (arcj(j,k,i) - arcj(j1,k,i))
     .          / (arcj(j2,k,i)- arcj(j1,k,i))
            psi(1) = eta
            psi(2) = 1.-eta
            x(j,k,i) = psi(2)*x(j1,k,i)
     .               + psi(1)*x(j2,k,i)
            y(j,k,i) = psi(2)*y(j1,k,i)
     .               + psi(1)*y(j2,k,i)
            z(j,k,i) = psi(2)*z(j1,k,i)
     .               + psi(1)*z(j2,k,i)
         end do
      end if
c
      if (k1.ne.k2) then
         i = i1
         j = j1
         do k=k1,k2
            eta = (arck(j,k,i) - arck(j,k1,i))
     .          / (arck(j,k2,i)- arck(j,k1,i))
            psi(1) = eta
            psi(2) = 1.-eta
            x(j,k,i) = psi(2)*x(j,k1,i)
     .               + psi(1)*x(j,k2,i)
            y(j,k,i) = psi(2)*y(j,k1,i)
     .               + psi(1)*y(j,k2,i)
            z(j,k,i) = psi(2)*z(j,k1,i)
     .               + psi(1)*z(j,k2,i)
         end do
      end if
c
      return
      end
