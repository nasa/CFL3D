      subroutine blnkfr(nbl,iibg,kkbg,jjbg,ibpntsg,lbg,iitot,blank,
     .                  jdim,kdim,idim,maxbl,blnkval)
c***********************************************************************
c     Purpose: Undo/redo blanking at fringe points (but not hole points)
c     blnkval = 1...set blank=1 at fringe points for plotting
c               0...set blank=0 (normal setting for fringe points)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension iibg(iitot),kkbg(iitot),jjbg(iitot),ibpntsg(maxbl,4),
     .          lbg(maxbl),blank(jdim,kdim,idim)
c
      lsta = lbg(nbl)
      lend = lsta-1
      if (ibpntsg(nbl,1).gt.0) then
         lend = lsta+ibpntsg(nbl,1)-1
         do l=lsta,lend
            blank(jjbg(l),kkbg(l),iibg(l)) = blnkval
         end do
      end if
c
      return
      end
