      subroutine renmbr(iseg,ix1,ix2,indx,mbloc,msegt,mtot,
     .                  imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                  nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                  nseg,idno,ijk,idseg,idnext)
c
c     $Id$
c
      dimension nsubbl(mbloc),idobl(mbloc),nseg(mbloc),idno(mbloc),
     .          ijk(6,mbloc),idseg(mbloc),idnext(mbloc)
      dimension imap(msegt,mtot),idbloc(mtot),
     .          ivisb(msegt,mbloc),itrb(7,mbloc),
     .          val(mtot),xdum(msegt,mtot),
     .          iold(4,mtot),nxtseg(mtot),intrfc(mtot),
     .          ipatch(mtot)
c
      imap(ix1,iseg) = imap(ix1,iseg) - (indx-1)
      imap(ix2,iseg) = imap(ix2,iseg) - (indx-1)
      if (imap(13,iseg).ne.0) then
         imap(ix1+11,iseg) = imap(ix1+11,iseg) - (indx-1)
         imap(ix2+11,iseg) = imap(ix2+11,iseg) - (indx-1)
      endif
      jseg = intrfc(iseg)
      if (jseg.gt.0) then
         jx1 = ix1 + 6
         jx2 = ix2 + 6
         if (imap(jx1,jseg).lt.imap(jx2,jseg)) then
            imap(jx1,jseg) = imap(ix1,iseg)
            imap(jx2,jseg) = imap(ix2,iseg)
         else
            imap(jx1,jseg) = imap(ix2,iseg)
            imap(jx2,jseg) = imap(ix1,iseg)
         endif
      endif
c
      return
      end
