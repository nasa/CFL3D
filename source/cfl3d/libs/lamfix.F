      subroutine lamfix(jdim,kdim,idim,smin,mdim,ndim,bcdata,
     .             ibeg,iend,jbeg,jend,kbeg,kend,nface)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Make smin negative in regions where soln is supposed
c     to be laminar
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension smin(jdim-1,kdim-1,idim-1)
      dimension bcdata(mdim,ndim,2,12)
c
c  J-dir:
      if (nface .eq. 3) then
        jbeg=1
        do i=ibeg,iend-1
          ii=i-ibeg+1
          do k=kbeg,kend-1
            kk=k-kbeg+1
            if(real(bcdata(kk,ii,1,3)) .lt. 1.) then
              jend=jdim
            else
              jend=int(bcdata(kk,ii,1,3)+.001)
            end if
            do j=jbeg,jend-1
              smin(j,k,i)=-(ccabs(smin(j,k,i)))
            enddo
          enddo
        enddo
      end if
      if (nface .eq. 4) then
        jend=jdim
        do i=ibeg,iend-1
          ii=i-ibeg+1
          do k=kbeg,kend-1
            kk=k-kbeg+1
            if(real(bcdata(kk,ii,1,3)) .lt. 1.) then
              jbeg=1
            else
              jbeg=jdim-int(bcdata(kk,ii,1,3)+.001)+1
            end if
            do j=jbeg,jend-1
              smin(j,k,i)=-(ccabs(smin(j,k,i)))
            enddo
          enddo
        enddo
      end if
c  K-dir:
      if (nface .eq. 5) then
        kbeg=1
        do i=ibeg,iend-1
          ii=i-ibeg+1
          do j=jbeg,jend-1
            jj=j-jbeg+1
            if(real(bcdata(jj,ii,1,3)) .lt. 1.) then
              kend=kdim
            else
              kend=int(bcdata(jj,ii,1,3)+.001)
            end if
            do k=kbeg,kend-1
              smin(j,k,i)=-(ccabs(smin(j,k,i)))
            enddo
          enddo
        enddo
      end if
      if (nface .eq. 6) then
        kend=kdim
        do i=ibeg,iend-1
          ii=i-ibeg+1
          do j=jbeg,jend-1
            jj=j-jbeg+1
            if(real(bcdata(jj,ii,1,3)) .lt. 1.) then
              kbeg=1
            else
              kbeg=kdim-int(bcdata(jj,ii,1,3)+.001)+1
            end if
            do k=kbeg,kend-1
              smin(j,k,i)=-(ccabs(smin(j,k,i)))
            enddo
          enddo
        enddo
      end if
c  I-dir:
      if (nface .eq. 1) then
        ibeg=1
        do k=kbeg,kend-1
          kk=k-kbeg+1
          do j=jbeg,jend-1
            jj=j-jbeg+1
            if(real(bcdata(jj,kk,1,3)) .lt. 1.) then
              iend=idim
            else
              iend=int(bcdata(jj,kk,1,3)+.001)
            end if
            do i=ibeg,iend-1
              smin(j,k,i)=-(ccabs(smin(j,k,i)))
            enddo
          enddo
        enddo
      end if
      if (nface .eq. 2) then
        iend=idim
        do k=kbeg,kend-1
          kk=k-kbeg+1
          do j=jbeg,jend-1
            jj=j-jbeg+1
            if(real(bcdata(jj,kk,1,3)) .lt. 1.) then
              ibeg=1
            else
              ibeg=idim-int(bcdata(jj,kk,1,3)+.001)+1
            end if
            do i=ibeg,iend-1
              smin(j,k,i)=-(ccabs(smin(j,k,i)))
            enddo
          enddo
        enddo
      end if
      return
      end
