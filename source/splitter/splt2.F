      subroutine splt2(iseg,ix1,ix2,indx,is,it,mbloc,msegt,mtot,
     .                 imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                 nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                 nseg,idno,ijk,idseg,idnext)
c
      dimension nsubbl(mbloc),idobl(mbloc),nseg(mbloc),idno(mbloc),
     .          ijk(6,mbloc),idseg(mbloc),idnext(mbloc)
      dimension imap(msegt,mtot),idbloc(mtot),
     .          ivisb(msegt,mbloc),itrb(7,mbloc),
     .          val(mtot),xdum(msegt,mtot),
     .          iold(4,mtot),nxtseg(mtot),intrfc(mtot),
     .          ipatch(mtot)
c
      call mkseg(iseg,it,ns,mbloc,msegt,mtot,
     .           imap,idbloc,ivisb,itrb,val,xdum,iold,
     .           nxtseg,intrfc,ipatch,nsubbl,idobl,
     .           nseg,idno,ijk,idseg,idnext)
c
      imap(ix2,iseg) = indx
      imap(ix1,ns) = indx
      if (imap(13,iseg).ne.0) then
         if (imap(ix1+11,iseg).ge.indx) then
            imap(13,iseg) = 0
         else if (imap(ix2+11,iseg).gt.indx) then
            imap(ix2+11,iseg) = indx
         endif
         if (imap(ix2+11,ns).le.indx) then
            imap(13,ns) = 0
         else if (imap(ix1+11,ns).lt.indx) then
            imap(ix1+11,ns) = indx
         endif
      endif
      iold(ix2-2,iseg) = iold(ix1-2,iseg) + (indx-imap(ix1,iseg))
      iold(ix1-2,ns) = iold(ix2-2,iseg)
c
      return
      end
